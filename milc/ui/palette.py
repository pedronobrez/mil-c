"""
Colours for the data, apart from the interface's own.

Categories take Paul Tol's qualitative set, which stays distinguishable
with the common colour-vision deficiencies; expression takes viridis,
which is perceptually even and survives greyscale printing.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg

from ..figures import CATEGORICAL  # noqa: E402,F401  (one list, for screen and paper)

CONTINUOUS = "viridis"


def categorical(n: int) -> list[str]:
    return [CATEGORICAL[i % len(CATEGORICAL)] for i in range(n)]


def colormap(name: str = CONTINUOUS):
    return pg.colormap.get(name)


def ramp(name: str = CONTINUOUS, steps: int = 256) -> np.ndarray:
    """steps x 4 uint8."""
    return colormap(name).getLookupTable(0.0, 1.0, steps, alpha=True)
