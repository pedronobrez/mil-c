"""The window, started: `python -m milc`, or `python -m milc gui FOLDER`."""

from __future__ import annotations

import sys


def main(folder: str | None = None, quit_after: float | None = None) -> int:
    """`quit_after` closes the window by itself after that many seconds:
    how a build is asked whether it starts at all."""
    from PyQt6 import QtWidgets

    from .ui import style
    from .ui.theme import apply_defaults

    app = QtWidgets.QApplication(sys.argv[:1])
    app.setApplicationName("MIL-C")
    # the suite's mark for the window, the taskbar and the Dock
    from pathlib import Path

    from PyQt6 import QtGui
    icon = Path(__file__).with_name("icon.png")
    if icon.is_file():
        app.setWindowIcon(QtGui.QIcon(str(icon)))
    style.apply(app)
    apply_defaults()

    from .ui.shell import MainShell

    window = MainShell()
    window.show()
    if folder:
        window.open_folder(folder)
    if quit_after is not None:
        from PyQt6 import QtCore
        QtCore.QTimer.singleShot(int(quit_after * 1000), app.quit)
    status = app.exec()
    # take the window down before the interpreter does: PyQt's exit hook
    # otherwise walks wrappers whose Qt objects are already gone
    window.close()
    window.deleteLater()
    app.processEvents()
    return status
