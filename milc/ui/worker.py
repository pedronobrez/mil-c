"""
Work off the window's thread: reading files, counting gates, running an
embedding. One function, one answer; the window stays alive meanwhile.
"""

from __future__ import annotations

import traceback

from PyQt6 import QtCore


class _Signals(QtCore.QObject):
    done = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str)
    progress = QtCore.pyqtSignal(int, int, str)


class _Task(QtCore.QRunnable):
    def __init__(self, fn, wants_progress: bool):
        super().__init__()
        self.fn, self.wants_progress = fn, wants_progress
        self.signals = _Signals()

    def run(self):
        try:
            result = (self.fn(self.signals.progress.emit) if self.wants_progress
                      else self.fn())
        except Exception as exc:
            traceback.print_exc()
            self.signals.failed.emit(f"{type(exc).__name__}: {exc}")
        else:
            self.signals.done.emit(result)


#: a task's signals die with it unless something holds them until it ends
_alive: set = set()


#: numerical work goes through here one job at a time: numba's default
#: threading layer does not allow parallel regions started from two
#: threads at once, and three embeddings fighting for the same cores are
#: no faster than three in a row.
_compute = QtCore.QThreadPool()
_compute.setMaxThreadCount(1)
# a secondary thread on macOS gets 512 KB of stack, and numba compiling a
# function recurses deeper than that: the process died with SIGILL the
# first time an algorithm ran off the main thread
_compute.setStackSize(256 * 1024 * 1024)


def run(fn, on_done=None, on_failed=None, on_progress=None,
        compute: bool = False) -> None:
    """`fn()` — or `fn(progress)` when `on_progress` is given, where
    `progress(done, total, what)` may be called from the work.
    `compute=True` for anything that reaches numba: those queue."""
    task = _Task(fn, on_progress is not None)
    _alive.add(task)

    def finish(*_):
        _alive.discard(task)

    if on_done:
        task.signals.done.connect(on_done)
    if on_failed:
        task.signals.failed.connect(on_failed)
    if on_progress:
        task.signals.progress.connect(on_progress)
    task.signals.done.connect(finish)
    task.signals.failed.connect(finish)
    task.setAutoDelete(False)
    (_compute if compute else QtCore.QThreadPool.globalInstance()).start(task)


def wait_all(timeout_ms: int = 600_000) -> None:
    """For scripts and tests: block until the pool is idle, then let the
    answers be delivered."""
    QtCore.QThreadPool.globalInstance().waitForDone(timeout_ms)
    _compute.waitForDone(timeout_ms)
    QtCore.QCoreApplication.processEvents()
