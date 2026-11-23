"""
Where the events are, as curves: a smoothed histogram of one parameter and
the contour lines of two. What the overlays are drawn from, on the screen
and on paper.

Everything here works in the plot's own coordinates — after the display
scale — because a histogram of a marker is a histogram of where its events
fall on the axis that shows them.
"""

from __future__ import annotations

import numpy as np

#: how a histogram's height is expressed
TO_MODE = "mode"      # its tallest point is 100: samples of any size compare
COUNT = "count"       # events per bin

#: the share of the events inside each contour line, outermost first
LEVELS = (0.90, 0.75, 0.50, 0.25, 0.10)


def _smooth(values: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return values
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(values.astype(np.float64), sigma, mode="constant")


def auto_smooth(n: int) -> float:
    """How many bins to smooth a histogram of `n` events over: 1.5 where
    there are thousands, more where there are few — a population of 300
    cells drawn at the resolution of 50,000 is a comb, and its teeth are
    read as peaks. Never more than 6 bins of 256, about 2% of the axis."""
    return float(np.clip(1.5 * (5000.0 / max(n, 1)) ** 0.3, 1.5, 6.0))


def histogram(x: np.ndarray, lo: float, hi: float, bins: int = 256,
              smooth: float | None = None,
              normalise: str = TO_MODE) -> tuple[np.ndarray, np.ndarray]:
    """(centres, heights). Events outside [lo, hi] pile up in the first and
    last bin, as they do on an axis: they are not dropped. `smooth` is in
    bins; None chooses by the number of events."""
    x = np.clip(np.asarray(x, dtype=np.float64), lo, hi)
    counts, edges = np.histogram(x, bins=bins, range=(lo, hi))
    if smooth is None:
        smooth = auto_smooth(x.size)
    heights = _smooth(counts, smooth)
    if normalise == TO_MODE:
        top = heights.max() if heights.size else 0.0
        heights = heights / top * 100.0 if top > 0 else heights
    return (edges[:-1] + edges[1:]) / 2, heights


def contours(x: np.ndarray, y: np.ndarray, xr: tuple[float, float],
             yr: tuple[float, float], bins: int = 96, smooth: float = 1.6,
             levels=LEVELS) -> list[tuple[float, list[np.ndarray]]]:
    """For each level, the lines that enclose that share of the events —
    probability contours, so a sample of 500 cells and one of 50,000 draw
    lines that mean the same thing. Each line is an n x 2 array."""
    from contourpy import contour_generator
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if x.size < 10:
        return []
    H, ex, ey = np.histogram2d(x, y, bins=bins, range=[xr, yr])
    D = _smooth(H, smooth)
    total = D.sum()
    if total <= 0:
        return []
    ordered = np.sort(D.ravel())[::-1]
    cumulative = np.cumsum(ordered) / total
    cx, cy = (ex[:-1] + ex[1:]) / 2, (ey[:-1] + ey[1:]) / 2
    generator = contour_generator(x=cx, y=cy, z=D.T)
    out = []
    for share in levels:
        i = min(int(np.searchsorted(cumulative, share)), ordered.size - 1)
        threshold = ordered[i]
        if threshold <= 0:
            continue
        lines = [np.asarray(line) for line in generator.lines(threshold) if len(line) > 2]
        if lines:
            out.append((float(share), lines))
    return out


def outside(x: np.ndarray, y: np.ndarray, xr, yr, bins: int = 96, smooth: float = 1.6,
            share: float = LEVELS[0]) -> np.ndarray:
    """Which events fall outside the outermost contour: drawn as dots, so
    a contour plot does not hide the few that are far from everyone."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if x.size < 10:
        return np.ones(x.size, dtype=bool)
    H, ex, ey = np.histogram2d(x, y, bins=bins, range=[xr, yr])
    D = _smooth(H, smooth)
    ordered = np.sort(D.ravel())[::-1]
    cumulative = np.cumsum(ordered) / max(D.sum(), 1e-300)
    threshold = ordered[min(int(np.searchsorted(cumulative, share)), ordered.size - 1)]
    ix = np.clip(np.searchsorted(ex, x, side="right") - 1, 0, bins - 1)
    iy = np.clip(np.searchsorted(ey, y, side="right") - 1, 0, bins - 1)
    return D[ix, iy] < threshold
