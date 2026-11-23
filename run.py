"""What a packaged MIL-C starts from: the command line, which opens the
window when it is given nothing."""

import multiprocessing
import os
import sys


def _prepare() -> None:
    """A packaged application's own folder is read-only, and matplotlib
    keeps its list of the system's fonts on disk between runs. PyInstaller
    points it at a temporary folder that is deleted on exit, so the list
    was rebuilt at every start: ten seconds before the first figure,
    measured. This is set after that and wins.

    numba needs no such help: it finds a cache folder of the user's by
    itself (measured in the bundle: the first UMAP after installing takes
    16 s, compiling; every one after that, 5 s).
    """
    if not getattr(sys, "frozen", False):
        return
    if sys.platform == "darwin":
        base = os.path.expanduser("~/Library/Caches")
    elif sys.platform.startswith("win"):
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
    else:
        base = os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache")
    os.environ["MPLCONFIGDIR"] = os.path.join(base, "MIL-C", "matplotlib")
    os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)


if __name__ == "__main__":
    # scikit-learn may start worker processes; in a packaged application
    # each of them would otherwise open a window of its own
    multiprocessing.freeze_support()
    _prepare()
    from milc.cli import main
    # Finder hands a double-clicked application a process serial number
    sys.exit(main([a for a in sys.argv[1:] if not a.startswith("-psn_")]))
