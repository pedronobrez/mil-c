"""
Kinetics: a parameter followed through the time of the acquisition, as
in a calcium flux — cells run for a while, a stimulus is added, and the
signal (or the ratio of two signals, as with Indo-1) rises and falls.

The acquisition is cut into slices of time of equal length. In each, one
number is taken from the events of the population: the median or the
mean of the parameter, or the percentage of events above a threshold.
With a second parameter the number is of the ratio of the two, event by
event (events whose denominator is not above zero are left out). The
series is then smoothed, by a moving average or a Gaussian, each slice
weighted alike and an empty slice skipped.

Between two pairs of marks the user sets — a baseline, before the
stimulus, and a response, after it — these are reported:

  Baseline        mean of the slices of the baseline, not smoothed, and
                  their SD: how much a slice moves by itself
  Peak            the highest value of the smoothed series in the
                  response, with its standard error
  Peak − baseline and Peak / baseline
  Time of peak    on the axis of the plot; Time to peak is counted from
                  the start of the response, where the stimulus is
  Slope           of a straight line through the slices (not smoothed)
                  from the start of the response to the peak, by least
                  squares, per unit of time, with its standard error
  AUC             area between the smoothed series and the baseline
                  over the response (trapezia): signal × time

The standard error of a slice is that of its statistic over its events:
SD/√n for a mean, 1.2533 × robust SD/√n for a median (exact for a normal
distribution, a guide otherwise), √(p(1−p)/n) for a percentage. Smoothing
is a weighted mean, and the error of a smoothed slice follows from those
of the slices in it.

Time is in seconds where the file says how long a unit of its Time
parameter is (`$TIMESTEP`), and in the file's own units where it does
not; the axis says which.

No Qt here.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import figures, order
from .qc import time_channel  # noqa: F401  (the Time parameter of a file)
from .qc_plots import timestep  # noqa: F401  (seconds in one unit of it)

MEDIAN, MEAN, ABOVE = "Median", "Mean", "% above threshold"
STATISTICS = [MEDIAN, MEAN, ABOVE]
NONE, MOVING, GAUSSIAN = "none", "moving", "gaussian"
SMOOTHINGS = {NONE: "None", MOVING: "Moving average", GAUSSIAN: "Gaussian"}
#: slices an acquisition is cut into when no length is given
SLICES = 100
#: a slice of fewer events than this has no number: one cell is not a median
MIN_EVENTS = 5


@dataclass
class Curve:
    """One sample through time."""
    label: str
    #: what its colour stands for: its own name, or the value of a column
    group: str
    #: the middle of each slice
    t: np.ndarray
    #: events in each slice
    n: np.ndarray
    #: the statistic of each slice (not a number where it has too few
    #: events), its standard error, and the two after smoothing
    y: np.ndarray
    se: np.ndarray
    smooth: np.ndarray
    smooth_se: np.ndarray


def _statistic(x: np.ndarray, statistic: str, threshold: float) -> tuple[float, float]:
    n = x.size
    if statistic == MEAN:
        return float(x.mean()), float(x.std(ddof=1) / np.sqrt(n))
    if statistic == ABOVE:
        p = float((x > threshold).mean())
        return 100.0 * p, 100.0 * float(np.sqrt(p * (1 - p) / n))
    lo, hi = np.percentile(x, [15.87, 84.13])
    return float(np.median(x)), float(1.2533 * (hi - lo) / 2 / np.sqrt(n))


def _weights(smoothing: str, window: float) -> np.ndarray:
    if smoothing == MOVING:
        size = max(int(round(window)), 1) | 1          # odd: centred on its slice
        return np.ones(size)
    if smoothing == GAUSSIAN and window > 0:
        reach = max(int(np.ceil(3 * window)), 1)
        return np.exp(-0.5 * (np.arange(-reach, reach + 1) / window) ** 2)
    return np.ones(1)


def smoothed(y: np.ndarray, se: np.ndarray, smoothing: str = MOVING,
             window: float = 5.0) -> tuple[np.ndarray, np.ndarray]:
    """The series smoothed, and the standard error of each smoothed
    value. `window` is in slices: the width of a moving average, the SD
    of a Gaussian. Slices without a number are skipped, and stay without
    one."""
    w = _weights(smoothing, window)
    have = np.isfinite(y)
    half = (w.size - 1) // 2

    def over(values, weights):      # centred on each slice, whatever the two lengths
        return np.convolve(values, weights)[half:half + y.size]

    total = over(have.astype(np.float64), w)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = over(np.where(have, y, 0.0), w) / total
        error = np.sqrt(over(np.where(have, se, 0.0) ** 2, w ** 2)) / total
    return np.where(have, mean, np.nan), np.where(have, error, np.nan)


def curve(time, values, label: str = "", group: str = "", step: float = 0.0,
          width: float | None = None, statistic: str = MEDIAN, threshold: float = 0.0,
          smoothing: str = MOVING, window: float = 5.0, denominator=None,
          span: tuple | None = None) -> Curve:
    """One sample's series. `time` and `values` are of the events of the
    population, as the file stores them; `step` the file's `$TIMESTEP`
    (0: time stays in the file's units); `width` the length of a slice in
    the units of the axis (None: `SLICES` of them over the acquisition);
    `denominator` a second parameter to divide by; `span` the (first,
    last) time to cut between, where several samples are to share their
    slices."""
    if statistic not in STATISTICS:
        raise ValueError(f"unknown statistic {statistic!r}")
    t = np.asarray(time, dtype=np.float64) * (step or 1.0)
    v = np.asarray(values, dtype=np.float64)
    if denominator is not None:
        d = np.asarray(denominator, dtype=np.float64)
        keep = d > 0
        t, v = t[keep], v[keep] / d[keep]
    keep = np.isfinite(t) & np.isfinite(v)
    t, v = t[keep], v[keep]
    if t.size == 0:
        raise ValueError("no events to follow")
    first, last = span if span is not None else (float(t.min()), float(t.max()))
    if not last > first:
        raise ValueError("every event has the same time")
    width = float(width) if width else (last - first) / SLICES
    count = max(int(np.ceil((last - first) / width - 1e-9)), 1)
    which = np.clip(((t - first) / width).astype(np.int64), 0, count - 1)
    which[(t < first) | (t > last)] = -1
    by = np.argsort(which, kind="stable")
    which, v = which[by], v[by]
    starts = np.searchsorted(which, np.arange(count + 1))
    n = np.diff(starts)
    y, se = np.full(count, np.nan), np.full(count, np.nan)
    for k in np.flatnonzero(n >= MIN_EVENTS):
        y[k], se[k] = _statistic(v[starts[k]:starts[k + 1]], statistic, threshold)
    smooth, smooth_se = smoothed(y, se, smoothing, window)
    return Curve(label, group or label, first + (np.arange(count) + 0.5) * width, n, y, se,
                 smooth, smooth_se)


#: what `measure` gives, in the order of the table
COLUMNS = ["Events", "Baseline", "Baseline SD", "Peak", "Peak SE", "Peak − baseline",
           "Peak / baseline", "Time of peak", "Time to peak", "Slope", "Slope SE", "AUC"]


def default_marks(first: float, last: float) -> tuple[tuple, tuple]:
    """(baseline, response) until the user says: the first sixth of the
    acquisition, and the rest."""
    cut = first + (last - first) / 6
    return (first, cut), (cut, last)


def measure(c: Curve, baseline: tuple, response: tuple) -> dict[str, float]:
    """The numbers of one series between the marks: (from, to) of the
    baseline and of the response, on the axis of the plot."""
    nan = float("nan")
    out = dict.fromkeys(COLUMNS, nan)
    out["Events"] = int(c.n.sum())
    have = np.isfinite(c.y)
    before = have & (c.t >= min(baseline)) & (c.t <= max(baseline))
    after = have & (c.t >= min(response)) & (c.t <= max(response))
    base = float(c.y[before].mean()) if before.any() else nan
    out["Baseline"] = base
    out["Baseline SD"] = float(c.y[before].std(ddof=1)) if before.sum() > 1 else nan
    if not after.any():
        return out
    slices = np.flatnonzero(after)
    top = slices[int(np.argmax(c.smooth[slices]))]
    peak = float(c.smooth[top])
    out.update({"Peak": peak, "Peak SE": float(c.smooth_se[top]),
                "Peak − baseline": peak - base,
                "Peak / baseline": peak / base if base else nan,
                "Time of peak": float(c.t[top]),
                "Time to peak": float(c.t[top] - min(response))})
    rise = slices[slices <= top]
    if rise.size >= 3:
        x, y = c.t[rise], c.y[rise]
        centred = x - x.mean()
        slope = float((centred * y).sum() / (centred ** 2).sum())
        left = y - y.mean() - slope * centred
        out["Slope"] = slope
        out["Slope SE"] = float(np.sqrt((left ** 2).sum() / (rise.size - 2)
                                        / (centred ** 2).sum()))
    if np.isfinite(base) and slices.size > 1:
        height, at = c.smooth[slices] - base, c.t[slices]
        out["AUC"] = float(((height[1:] + height[:-1]) / 2 * np.diff(at)).sum())
    return out


# --- the figure -----------------------------------------------------------------------

@dataclass
class KineticsContent:
    """Several samples through time, with the marks."""
    series: list[Curve] = field(default_factory=list)
    xlabel: str = "Time"
    ylabel: str = ""
    baseline: tuple | None = None
    response: tuple | None = None
    #: the slices before smoothing, as dots under each line
    points: bool = True
    uses_colormap = False

    def groups(self) -> list[str]:
        return order.ordered(s.group for s in self.series)

    def category_names(self) -> list[tuple[str, str]]:
        names = self.groups()
        return list(zip(names, figures.categorical(len(names))))


def _draw_kinetics(fig, c: KineticsContent, style: figures.FigureStyle) -> None:
    from matplotlib.lines import Line2D
    ax = fig.subplots()
    colour = {g: style.color(g, default) for g, default in c.category_names()}
    if c.baseline is not None:
        ax.axvspan(min(c.baseline), max(c.baseline), color="#888888", alpha=0.13, linewidth=0)
    if c.response is not None:
        for at in c.response:
            ax.axvline(at, color="#888888", linewidth=0.6, linestyle="--")
    for s in c.series:
        if c.points:
            ax.plot(s.t, s.y, "o", color=colour[s.group], markersize=style.point_size * 0.7,
                    alpha=0.35, markeredgewidth=0, rasterized=style.rasterize)
        ax.plot(s.t, s.smooth, color=colour[s.group], linewidth=1.1)
    ax.set_xlabel(c.xlabel if style.xlabel is None else style.xlabel)
    ax.set_ylabel(c.ylabel if style.ylabel is None else style.ylabel)
    ax.margins(x=0.01)
    figures._frame(ax, style)
    handles = [Line2D([], [], color=colour[g], linewidth=1.6, label=style.label(g))
               for g in colour]
    if len(handles) > 1:
        figures._legend(fig, [ax], handles, style, "")


figures._DRAW[KineticsContent] = _draw_kinetics
