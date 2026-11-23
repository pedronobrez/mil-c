"""
The one place a settings object is made. `MILC_SETTINGS`, set by the test
suite before any widget exists, points it at an INI file so a test never
writes into the person's real preferences.
"""

from __future__ import annotations

import os

from PyQt6 import QtCore

ENV_FILE = "MILC_SETTINGS"


def settings() -> QtCore.QSettings:
    path = os.environ.get(ENV_FILE)
    if path:
        return QtCore.QSettings(path, QtCore.QSettings.Format.IniFormat)
    return QtCore.QSettings("MIL-C", "MIL-C")
