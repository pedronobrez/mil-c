"""
Colours for the data, apart from the interface's own.

Categories take Paul Tol's qualitative set, which stays distinguishable
with the common colour-vision deficiencies. The ramp — for expression on
a map, the heatmap, the tree — is the user's choice: any of the named
ones or one made of their own colours, remembered between sessions, and
the same one a figure starts from.
"""

from __future__ import annotations

import json

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore

from .. import figures
from ..figures import CATEGORICAL  # noqa: F401  (one list, for screen and paper)
from .settings import settings

SETTING_RAMP = "colors/ramp"
SETTING_CUSTOM = "colors/custom"


class _Notifier(QtCore.QObject):
    changed = QtCore.pyqtSignal()


#: emits `changed` when the ramp is: every plot that uses it redraws
notifier = _Notifier()


def categorical(n: int) -> list[str]:
    return [CATEGORICAL[i % len(CATEGORICAL)] for i in range(n)]


def ramp_name() -> str:
    name = str(settings().value(SETTING_RAMP, "viridis") or "viridis")
    return name if name in figures.COLORMAPS or name == figures.CUSTOM else "viridis"


def custom_colors() -> list[str]:
    try:
        colors = json.loads(settings().value(SETTING_CUSTOM, "") or "[]")
    except (TypeError, ValueError):
        colors = []
    return [str(c) for c in colors] or list(figures.DEFAULT_CUSTOM)


def set_ramp(name: str, colors: list[str] | None = None) -> None:
    store = settings()
    store.setValue(SETTING_RAMP, name)
    if colors is not None:
        store.setValue(SETTING_CUSTOM, json.dumps(list(colors)))
    notifier.changed.emit()


def ramp(steps: int = 256) -> np.ndarray:
    """The current ramp, steps x 4 uint8."""
    return figures.lookup_table(ramp_name(), custom_colors(), steps)


def colormap():
    """The current ramp for pyqtgraph."""
    lut = ramp(64)
    return pg.ColorMap(np.linspace(0.0, 1.0, len(lut)), lut)
