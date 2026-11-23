"""
Gates the program proposes: where a histogram dips between its two
peaks, the outline of a population on two parameters, the band along a
diagonal. A proposal is an ordinary gate of `milc.gating` — a range, a
polygon, an ellipse — to accept, move or throw away like one drawn by
hand; and when the events do not show what a method looks for, it says so
and proposes nothing.

Like `milc.density`, everything is measured in the plot's own coordinates
— after the display scale — because a population is a peak where it is
looked at, not in the units the file stores.

A proposal is described by a `spec`, plain values that say how it was
made, so the same gate can be fitted again on another sample's own events:
that is what these are for when samples shift from one day to the next.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .density import _smooth
from .gating import EllipseGate, PolygonGate, RectangleGate, points_in_polygon

VALLEY, POPULATION, ELLIPSE, BAND = "valley", "population", "ellipse", "band"
ABOVE, BELOW = "above", "below"

#: what each method is called where a person reads it
NAMES = {VALLEY: "valley between two peaks", POPULATION: "density contour",
         ELLIPSE: "robust Gaussian ellipse", BAND: "band along a robust line"}


# --- one parameter: the valley -------------------------------------------------

@dataclass
class Valley:
    """Where a histogram splits in two — or `threshold` None, and why."""
    threshold: float | None
    why: str
    #: the two peaks, low then high
    modes: tuple[float, float] | None = None
    #: the height of the valley as a share of the lower of the two peaks
    depth: float | None = None
    #: the share of the events above the threshold
    above: float | None = None


def valley(x, lo: float, hi: float, bins: int = 256, max_depth: float = 0.7,
           min_share: float = 0.01) -> Valley:
    """The density minimum between the two main peaks of `x`, on an axis
    that runs from `lo` to `hi`.

    The histogram is taken over where the events are (their 0.1% and
    99.9% points, not the whole axis), smoothed, and its peaks are those
    that stand out from their surroundings by five times what chance
    alone does to a bin — not by a share of the tallest: a population of
    a few per cent is a low, wide peak beside a tall one, and a real one.
    The two that stand out most are the two populations. Nothing is
    proposed when there is one peak only, when the valley is not lower
    than `max_depth` of the lower peak, or when one side would hold less
    than `min_share` of the events: a threshold is then a matter of
    judgement, or of a control."""
    from scipy.signal import find_peaks
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x) & (x > lo) & (x < hi)]     # not the piles at the ends
    if x.size < 100:
        return Valley(None, f"too few events on the axis to look for two peaks ({x.size})")
    a, b = np.quantile(x, [0.001, 0.999])
    pad = 0.05 * (b - a)
    a, b = max(lo, a - pad), min(hi, b + pad)
    if not b > a:
        return Valley(None, "every event has the same value")
    counts, edges = np.histogram(x, bins=bins, range=(a, b))
    sigma = float(np.clip(3.0 * (5000.0 / x.size) ** 0.3, 3.0, 10.0))
    h = np.concatenate([[0.0], _smooth(counts, sigma), [0.0]])   # an end can be a peak
    peaks, found = find_peaks(h, prominence=0.0)
    # how far a bin of the smoothed histogram strays by chance alone: counts
    # are Poisson, and smoothing over `sigma` bins averages 2·sqrt(pi)·sigma
    chance = np.sqrt(h[peaks] / (2.0 * np.sqrt(np.pi) * sigma))
    real = found["prominences"] >= np.maximum(5.0 * chance, 0.001 * h.max())
    peaks, rise = peaks[real], found["prominences"][real]
    centres = (edges[:-1] + edges[1:]) / 2
    if len(peaks) < 2:
        return Valley(None, "one peak only: there is no valley to split at")
    first, second = sorted(peaks[np.argsort(rise)[-2:]])
    between = h[first:second + 1]
    low = np.flatnonzero(between <= between.min() * (1 + 1e-9)) + first
    at = int(low[len(low) // 2])
    depth = float(h[at] / min(h[first], h[second]))
    threshold = float(centres[at - 1])
    above = float((x >= threshold).mean())
    modes = (float(centres[first - 1]), float(centres[second - 1]))
    if depth > max_depth:
        return Valley(None, f"two peaks, but the dip between them is shallow: at its lowest "
                            f"the histogram is still {100 * depth:.0f}% of the lower peak",
                      modes, depth, above)
    if min(above, 1 - above) < min_share:
        return Valley(None, f"two peaks, but one side of the valley would hold only "
                            f"{100 * min(above, 1 - above):.1f}% of the events",
                      modes, depth, above)
    return Valley(threshold, f"valley at {100 * depth:.0f}% of the lower peak; "
                             f"{100 * above:.1f}% of the events above it", modes, depth, above)


# --- two parameters: where the events are ---------------------------------------

@dataclass
class Mode:
    """One population of a plot: a peak of the density and the ground
    that slopes up to it."""
    #: where the peak is
    x: float
    y: float
    #: the share of the events on its ground
    share: float
    peak: float


class Landscape:
    """The density of two parameters on a grid, cut into the grounds of
    its peaks: each cell belongs to the peak reached by walking uphill
    from it. A peak that barely stands above the pass to a neighbour (the
    pass is higher than `min_dip` of it), or whose ground holds less than
    `min_share` of the events, is a ripple of the neighbour and joins it."""

    def __init__(self, x, y, xr, yr, bins: int = 128, smooth: float | None = None,
                 min_dip: float = 0.85, min_share: float = 0.01):
        x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
        ok = (np.isfinite(x) & np.isfinite(y)
              & (x > xr[0]) & (x < xr[1]) & (y > yr[0]) & (y < yr[1]))
        self.x, self.y, self.xr, self.yr, self.bins = x[ok], y[ok], xr, yr, bins
        self.n = int(ok.sum())
        H, self.ex, self.ey = np.histogram2d(self.x, self.y, bins=bins, range=[xr, yr])
        # the events piled at the ends of the axes are not a population
        H[0, :] = H[-1, :] = H[:, 0] = H[:, -1] = 0
        self.H = H
        if smooth is None:
            smooth = float(np.clip(2.0 * (20000.0 / max(self.n, 1)) ** 0.25, 2.0, 6.0))
        self.D = _smooth(H, smooth)
        self.cx = (self.ex[:-1] + self.ex[1:]) / 2
        self.cy = (self.ey[:-1] + self.ey[1:]) / 2
        self.labels, self.modes = self._grounds(min_dip, min_share)

    def _grounds(self, min_dip: float, min_share: float):
        D, n = self.D, self.bins
        total = D.sum()
        if total <= 0:
            return np.full(D.shape, -1), []
        floor = D.max() * 1e-6
        padded = np.full((n + 2, n + 2), -np.inf)
        padded[1:-1, 1:-1] = D
        # each cell points at the highest of itself and its eight neighbours
        best, target = D.copy(), np.arange(n * n).reshape(n, n)
        index = target.copy()
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                other = padded[1 + di:n + 1 + di, 1 + dj:n + 1 + dj]
                higher = other > best
                best = np.where(higher, other, best)
                shifted = index + di * n + dj
                target = np.where(higher, shifted, target)
        target = target.ravel()
        while True:
            further = target[target]
            if np.array_equal(further, target):
                break
            target = further
        root = target.reshape(n, n)
        root = np.where(D > floor, root, -1)
        ids = [int(r) for r in np.unique(root) if r >= 0]
        peak = {r: float(D.ravel()[r]) for r in ids}
        mass = {r: float(D[root == r].sum()) for r in ids}
        # the pass between two neighbours: the highest point of their border
        passes: dict[tuple[int, int], float] = {}
        for a, b, da, db in ((root[:-1], root[1:], D[:-1], D[1:]),
                             (root[:, :-1], root[:, 1:], D[:, :-1], D[:, 1:])):
            differ = (a != b) & (a >= 0) & (b >= 0)
            lo_, hi_ = np.minimum(a[differ], b[differ]), np.maximum(a[differ], b[differ])
            height = np.minimum(da[differ], db[differ])
            for i, j, h in zip(lo_.tolist(), hi_.tolist(), height.tolist()):
                if h > passes.get((i, j), -1.0):
                    passes[(i, j)] = h
        parent = {r: r for r in ids}
        alive = set(ids)
        while len(alive) > 1:
            merged = False
            for r in sorted(alive, key=lambda r: peak[r]):
                around = [(h, j if i == r else i) for (i, j), h in passes.items()
                          if r in (i, j)]
                if not around:
                    continue
                h, into = max(around)
                if h > min_dip * peak[r] or mass[r] < min_share * total:
                    for (i, j), v in list(passes.items()):
                        if r in (i, j):
                            del passes[(i, j)]
                            k = j if i == r else i
                            if k != into:
                                key = (min(k, into), max(k, into))
                                passes[key] = max(passes.get(key, -1.0), v)
                    parent[r] = into
                    mass[into] += mass[r]
                    if peak[r] > peak[into]:
                        peak[into] = peak[r]
                    alive.discard(r)
                    merged = True
                    break
            if not merged:
                break

        def top(r):
            while parent[r] != r:
                r = parent[r]
            return r
        final = {r: top(r) for r in ids}
        kept = sorted(alive, key=lambda r: -mass[r])
        number = {r: k for k, r in enumerate(kept)}
        labels = np.full(D.shape, -1)
        for r in ids:
            labels[root == r] = number[final[r]]
        modes = []
        for k, r in enumerate(kept):
            own = np.where(labels == k, D, -1.0)
            i, j = np.unravel_index(int(np.argmax(own)), D.shape)
            modes.append(Mode(float(self.cx[i]), float(self.cy[j]), mass[r] / total,
                              float(D[i, j])))
        return labels, modes

    def cells(self, x, y) -> tuple[np.ndarray, np.ndarray]:
        i = np.clip(np.searchsorted(self.ex, x, side="right") - 1, 0, self.bins - 1)
        j = np.clip(np.searchsorted(self.ey, y, side="right") - 1, 0, self.bins - 1)
        return i, j

    def choose(self, near=None) -> int | None:
        """The population a point is on — the nearest peak when it is on
        none — or, without a point, the one with the most events."""
        if not self.modes:
            return None
        if near is None:
            return 0
        i, j = self.cells(float(near[0]), float(near[1]))
        if self.labels[i, j] >= 0:
            return int(self.labels[i, j])
        w, h = self.xr[1] - self.xr[0], self.yr[1] - self.yr[0]
        return int(np.argmin([((m.x - near[0]) / w) ** 2 + ((m.y - near[1]) / h) ** 2
                              for m in self.modes]))

    def core(self, k: int, share: float) -> tuple[np.ndarray, float]:
        """The cells of population `k` that hold `share` of its events,
        the densest first, and the density at their edge. The events are
        counted as they are, not as smoothing spreads them."""
        mine = self.labels == k
        own = np.where(mine, self.D, 0.0)
        order = np.argsort(own[mine])[::-1]
        counts = self.H[mine][order]
        cumulative = np.cumsum(counts) / max(counts.sum(), 1.0)
        at = min(int(np.searchsorted(cumulative, share)), order.size - 1)
        level = float(own[mine][order][at])
        return own >= level, level


def _simplify(line: np.ndarray, tolerance: float) -> np.ndarray:
    """Douglas and Peucker's: drop the points a straight side passes
    within `tolerance` of."""
    keep = np.zeros(len(line), dtype=bool)
    keep[[0, -1]] = True
    stack = [(0, len(line) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        seg = line[b] - line[a]
        rel = line[a + 1:b] - line[a]
        length = np.hypot(*seg)
        d = (np.abs(seg[0] * rel[:, 1] - seg[1] * rel[:, 0]) / length if length > 0
             else np.hypot(rel[:, 0], rel[:, 1]))
        k = int(np.argmax(d))
        if d[k] > tolerance:
            keep[a + 1 + k] = True
            stack += [(a, a + 1 + k), (a + 1 + k, b)]
    return line[keep]


@dataclass
class Outline:
    """A population's outline as a polygon — or `vertices` None, and why."""
    vertices: np.ndarray | None
    why: str
    mode: Mode | None = None
    #: how many populations the plot was cut into
    populations: int = 0


def population(x, y, xr, yr, near=None, share: float = 0.9, max_vertices: int = 32,
               **settings) -> Outline:
    """The contour that holds `share` of one population's events: the one
    `near` (a point of the plot) is on, or the one with the most events.

    The plot is cut into populations first (`Landscape`), so a contour
    stops at the valley to a neighbour instead of running round both;
    the line is then thinned to at most `max_vertices` corners, to be
    moved by hand like any polygon."""
    from contourpy import contour_generator
    land = Landscape(x, y, xr, yr, **settings)
    if land.n < 100:
        return Outline(None, f"too few events on the plot to find a population ({land.n})")
    k = land.choose(near)
    if k is None:
        return Outline(None, "no events on the plot")
    mode = land.modes[k]
    inside, level = land.core(k, share)
    own = np.where(land.labels == k, land.D, 0.0)
    # a ring of nothing around the grid, so a contour at its edge closes
    z = np.zeros((land.bins + 2, land.bins + 2))
    z[1:-1, 1:-1] = own
    step_x, step_y = land.cx[1] - land.cx[0], land.cy[1] - land.cy[0]
    gx = np.concatenate([[land.cx[0] - step_x], land.cx, [land.cx[-1] + step_x]])
    gy = np.concatenate([[land.cy[0] - step_y], land.cy, [land.cy[-1] + step_y]])
    lines = [np.asarray(line) for line in
             contour_generator(x=gx, y=gy, z=z.T).lines(level) if len(line) > 3]
    if not lines:
        return Outline(None, "the population is too small to draw a contour round", mode,
                       len(land.modes))
    peak = (np.array([mode.x]), np.array([mode.y]))
    holding = [line for line in lines if points_in_polygon(*peak, line)[0]]
    line = max(holding or lines, key=len)
    line = np.column_stack([np.clip(line[:, 0], *xr), np.clip(line[:, 1], *yr)])
    scale = np.array([xr[1] - xr[0], yr[1] - yr[0]], dtype=np.float64)
    unit, tolerance = line / scale, 0.002
    thin = _simplify(unit, tolerance)
    while len(thin) > max_vertices + 1:
        tolerance *= 1.5
        thin = _simplify(unit, tolerance)
    vertices = thin[:-1] * scale if np.allclose(thin[0], thin[-1]) else thin * scale
    if len(vertices) < 3:
        return Outline(None, "the population is too small to draw a contour round", mode,
                       len(land.modes))
    held = float(points_in_polygon(land.x, land.y, vertices).mean())
    others = len(land.modes) - 1
    return Outline(vertices,
                   f"the {100 * share:.0f}% contour of a population holding "
                   f"{100 * mode.share:.0f}% of the events, {100 * held:.1f}% of them "
                   f"inside it; "
                   + (f"{others} other population{'s' if others > 1 else ''} on the plot"
                      if others else "the only population on the plot"),
                   mode, len(land.modes))


@dataclass
class Gaussian:
    """A population as an ellipse — or `center` None, and why."""
    center: tuple[float, float] | None
    radii: tuple[float, float] | None
    #: counter-clockwise, in radians, of the first half-axis from the X axis
    angle: float
    why: str
    mode: Mode | None = None


def gaussian(x, y, xr, yr, near=None, share: float = 0.95, trim: float = 0.9,
             **settings) -> Gaussian:
    """One population as a two-dimensional Gaussian fitted to its core,
    and the ellipse that holds `share` of such a Gaussian.

    The fit starts from the densest half of the population `near` is on
    (or of the largest), then repeats: keep the events within the
    Mahalanobis distance that holds `trim` of a Gaussian, take their mean
    and covariance — the covariance made larger by what cutting the tails
    off a Gaussian takes away — until nothing moves. Events of other
    populations and stragglers are outside the cut and do not pull the
    ellipse towards them, which is what `flowCore::norm2Filter` (a robust
    covariance) and a one-cluster `flowClust` are used for."""
    land = Landscape(x, y, xr, yr, **settings)
    if land.n < 100:
        return Gaussian(None, None, 0.0,
                        f"too few events on the plot to find a population ({land.n})")
    k = land.choose(near)
    if k is None:
        return Gaussian(None, None, 0.0, "no events on the plot")
    mode = land.modes[k]
    # the events on the population's own ground: its neighbours are not
    # tails to trim but other populations, denser ones included
    cells = land.cells(land.x, land.y)
    own = land.labels[cells] == k
    P = np.column_stack([land.x[own], land.y[own]])
    core, _ = land.core(k, 0.5)
    start = P[core[cells][own]]
    if len(start) < 20:
        return Gaussian(None, None, 0.0, "the population is too small to fit", mode)
    cut = -2.0 * np.log(1.0 - trim)                    # chi-square, two dimensions
    # what is left of a Gaussian's covariance inside that distance
    kept = (1.0 - np.exp(-cut / 2) * (1.0 + cut / 2)) / (1.0 - np.exp(-cut / 2))
    half = -2.0 * np.log(0.5)
    kept_half = (1.0 - np.exp(-half / 2) * (1.0 + half / 2)) / (1.0 - np.exp(-half / 2))
    mean, cov = start.mean(axis=0), np.cov(start.T) / kept_half
    for _ in range(200):
        try:
            inv = np.linalg.inv(cov)
        except np.linalg.LinAlgError:
            return Gaussian(None, None, 0.0, "the population has no width to fit", mode)
        d = P - mean
        within = np.einsum("ij,jk,ik->i", d, inv, d) <= cut
        if within.sum() < 20:
            return Gaussian(None, None, 0.0, "the population is too small to fit", mode)
        new_mean, new_cov = P[within].mean(axis=0), np.cov(P[within].T) / kept
        moved = max(np.abs(new_mean - mean).max(), np.abs(new_cov - cov).max() ** 0.5)
        mean, cov = new_mean, new_cov
        if moved < 1e-6 * max(xr[1] - xr[0], yr[1] - yr[0]):
            break
    values, vectors = np.linalg.eigh(cov)
    if values.min() <= 0:
        return Gaussian(None, None, 0.0, "the population has no width to fit", mode)
    reach = np.sqrt(-2.0 * np.log(1.0 - share))
    major = vectors[:, 1]
    radii = (float(reach * np.sqrt(values[1])), float(reach * np.sqrt(values[0])))
    d = np.column_stack([land.x, land.y]) - mean
    held = float((np.einsum("ij,jk,ik->i", d, np.linalg.inv(cov), d) <= reach ** 2).mean())
    return Gaussian((float(mean[0]), float(mean[1])), radii,
                    float(np.arctan2(major[1], major[0])),
                    f"the ellipse that holds {100 * share:.0f}% of a Gaussian fitted to the "
                    f"population's core; {100 * held:.1f}% of the events inside", mode)


@dataclass
class Band:
    """The events along a line, as a polygon of four corners — or
    `vertices` None, and why."""
    vertices: np.ndarray | None
    why: str
    slope: float | None = None
    intercept: float | None = None
    #: half the band's height
    half_width: float | None = None


def band(x, y, xr, yr, width: float = 3.0, tight: bool = False) -> Band:
    """The band of events along a straight line of Y against X: singlets
    on a pulse's area against its height, where a doublet has more area
    for its height and falls off the line.

    The line is fitted with Tukey's biweight, which gives the events far
    from it no say, and the spread about it is the median distance of
    the events the fit listens to, scaled to a standard deviation. The
    band reaches `width` of those either side, from the lowest to the
    highest X of the events it holds. That is what `flowStats::
    singletGate` does (a robust linear model and a prediction band).

    `tight` is for a line with a fan on one side of it, which widens
    that spread: the line is moved onto the ridge, where the events are
    densest across it, and the spread is that of its cleaner side alone.

    Nothing is proposed when the events do not follow a line."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    ok = (np.isfinite(x) & np.isfinite(y)
          & (x > xr[0]) & (x < xr[1]) & (y > yr[0]) & (y < yr[1]))
    x, y = x[ok], y[ok]
    if x.size < 100:
        return Band(None, f"too few events on the plot to fit a line ({x.size})")
    spread_x = np.subtract(*np.quantile(x, [0.75, 0.25]))
    if spread_x <= 0:
        return Band(None, "every event has the same X")
    mx, my = np.median(x), np.median(y)
    steep = (x - mx) != 0
    slope = float(np.median((y - my)[steep] / (x - mx)[steep]))
    intercept = float(my - slope * mx)

    def sigma_of(r, s=None):
        # of the events the fit listens to: the fan of doublets above a
        # line of singlets is not part of the line's own spread
        if s:
            r = r[np.abs(r) < 4.685 * s]
        return 1.4826 * float(np.median(np.abs(r))) if r.size else 0.0
    s = sigma_of(y - (intercept + slope * x) - np.median(y - (intercept + slope * x)))
    for _ in range(100):
        r = y - (intercept + slope * x)
        s = sigma_of(r, s)
        if s <= 0:
            break
        u = r / (4.685 * s)
        w = np.where(np.abs(u) < 1, (1 - u ** 2) ** 2, 0.0)
        sw = w.sum()
        xm, ym = (w * x).sum() / sw, (w * y).sum() / sw
        sxx = (w * (x - xm) ** 2).sum()
        if sxx <= 0:
            break
        new_slope = float((w * (x - xm) * (y - ym)).sum() / sxx)
        new_intercept = float(ym - new_slope * xm)
        done = (abs(new_slope - slope) <= 1e-9 * max(abs(slope), 1e-12)
                and abs(new_intercept - intercept) <= 1e-9 * (yr[1] - yr[0]))
        slope, intercept = new_slope, new_intercept
        if done:
            break
    r = y - (intercept + slope * x)
    s = sigma_of(r, s)
    if s > 0 and tight:
        # the line is moved onto the ridge — where the events are densest
        # across it — and its spread is taken from the cleaner side: what
        # falls off a line of singlets falls off one side of it only
        near = r[np.abs(r) < 4.685 * s]
        counts, edges = np.histogram(near, bins=101)
        ridge = float(((edges[:-1] + edges[1:]) / 2)[int(np.argmax(_smooth(counts, 2.0)))])
        intercept += ridge
        r = r - ridge
        sides = [1.4826 * float(np.median(np.abs(side))) for side in (r[r < 0], r[r > 0])
                 if side.size >= 20]
        s = min(sides) if sides else s
    spread_y = 1.4826 * float(np.median(np.abs(y - my)))
    if s <= 0 or spread_y <= 0:
        return Band(None, "the events have no spread to put a band round")
    if s > 0.5 * spread_y:
        return Band(None, "the events do not follow a line: the spread about the best one "
                          f"is {100 * s / spread_y:.0f}% of the spread of Y itself")
    half = width * s
    held = np.abs(r) <= half
    x0, x1 = float(x[held].min()), float(x[held].max())
    pad = 0.005 * (xr[1] - xr[0])
    x0, x1 = x0 - pad, x1 + pad
    vertices = np.array([[x0, intercept + slope * x0 - half],
                         [x1, intercept + slope * x1 - half],
                         [x1, intercept + slope * x1 + half],
                         [x0, intercept + slope * x0 + half]])
    return Band(vertices,
                f"a band along the line Y = {slope:.3g}·X "
                f"{'+' if intercept >= 0 else '−'} {abs(intercept):.3g}, "
                f"{width:g} standard deviations {'of its cleaner side ' if tight else ''}"
                f"to each side; {100 * held.mean():.1f}% of the events inside",
                slope, intercept, half)


# --- from a description to a gate, on any sample -----------------------------------

def make_spec(method: str, dims: list[str], **settings) -> dict:
    """How a gate is fitted, as plain values: the method, the parameters
    in the order of the plot's axes, and the method's settings — `side`
    (ABOVE or BELOW) for a valley, `near` (a point of the plot, or None
    for the largest population) and `share` for a contour or an ellipse,
    `width` and `tight` for a band."""
    spec = {"method": method, "dims": list(dims)}
    for key, value in settings.items():
        if value is not None:
            spec[key] = ([float(v) for v in value] if isinstance(value, (list, tuple, np.ndarray))
                         else value)
    return spec


#: begins the `note` of a population whose gate was to be fitted on the
#: sample's own events and could not be
NOT_FITTED = "not fitted to this sample"


def note(spec: dict) -> str:
    """What to say of a population whose gate each sample has its own of."""
    return (f"fitted to this sample's own events ({NAMES[spec['method']]}): each sample "
            "has its own, and moving this one moves it here only")


@dataclass
class Proposal:
    """A gate to accept or not — or `gate` None; `why` says what was
    found, or why nothing is proposed."""
    gate: object | None
    why: str
    #: where the population found is, on the plot: the point the same
    #: population is looked for near on another sample
    near: tuple[float, float] | None = None
    #: the outline on the plot, n x 2 (two points for a threshold: its ends
    #: on X, Y left for the plot to give)
    shown: np.ndarray | None = None

    @property
    def found(self) -> bool:
        return self.gate is not None


def _stored(points: np.ndarray, scales: list) -> np.ndarray:
    out = np.array(points, dtype=np.float64)
    for k, scale in enumerate(scales):
        if scale is not None:
            out[:, k] = scale.inverse(out[:, k])
    return out


def propose(spec: dict, values: Callable[[str], np.ndarray], ranges: list[tuple[float, float]],
            scales: list) -> Proposal:
    """Fit the gate `spec` describes. `values(dim)` gives the events to
    fit it on — the parent population's — on the display scale; `ranges`
    and `scales` are the ends and the scale (None: as stored) of each of
    the spec's parameters' axes. The gate is in the file's units where a
    hand-drawn one would be."""
    method, dims = spec["method"], list(spec["dims"])
    if method not in NAMES:
        raise ValueError(f"no such way to fit a gate: {method}")
    any_scale = any(s is not None for s in scales)
    if method == VALLEY:
        found = valley(values(dims[0]), *ranges[0])
        if found.threshold is None:
            return Proposal(None, found.why)
        at = float(_stored([[found.threshold]], scales[:1])[0, 0])
        above = spec.get("side", ABOVE) == ABOVE
        gate = RectangleGate([dims[0]], [at if above else None], [None if above else at])
        ends = (found.threshold, ranges[0][1]) if above else (ranges[0][0], found.threshold)
        return Proposal(gate, found.why, shown=np.array(ends))
    x, y = values(dims[0]), values(dims[1])
    if method == POPULATION:
        found = population(x, y, ranges[0], ranges[1], spec.get("near"),
                           float(spec.get("share", 0.9)))
        if found.vertices is None:
            return Proposal(None, found.why)
        gate = PolygonGate(dims, _stored(found.vertices, scales),
                           scales=list(scales) if any_scale else None)
        return Proposal(gate, found.why, (found.mode.x, found.mode.y),
                        np.vstack([found.vertices, found.vertices[:1]]))
    if method == ELLIPSE:
        found = gaussian(x, y, ranges[0], ranges[1], spec.get("near"),
                         float(spec.get("share", 0.95)))
        if found.center is None:
            return Proposal(None, found.why)
        gate = EllipseGate(dims, list(found.center), list(found.radii), found.angle,
                           scales=list(scales) if any_scale else None)
        t = np.linspace(0, 2 * np.pi, 121)
        c, s = np.cos(found.angle), np.sin(found.angle)
        u, v = found.radii[0] * np.cos(t), found.radii[1] * np.sin(t)
        return Proposal(gate, found.why, found.center,
                        np.column_stack([found.center[0] + u * c - v * s,
                                         found.center[1] + u * s + v * c]))
    found = band(x, y, ranges[0], ranges[1], float(spec.get("width", 3.0)),
                 bool(spec.get("tight", False)))
    if found.vertices is None:
        return Proposal(None, found.why)
    gate = PolygonGate(dims, _stored(found.vertices, scales),
                       scales=list(scales) if any_scale else None)
    return Proposal(gate, found.why, shown=np.vstack([found.vertices, found.vertices[:1]]))


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    """How much two sets of events are the same one: those in both over
    those in either. 1 when both are empty."""
    either = int((a | b).sum())
    return float((a & b).sum() / either) if either else 1.0
