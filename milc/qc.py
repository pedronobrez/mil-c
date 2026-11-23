"""
Quality control of an acquisition: which stretches of time went wrong.

A clog, an air bubble, a pressure change or a drifting laser spoils the
events acquired while it lasts. They are looked for the way flowAI and
PeacoQC look for them — in the order of acquisition, a fixed number of
events at a time — and what is found is a list of time intervals, so that
it can be a gate like any other (`gating.TimeGate`).

The method. The events are put in the order of Time and cut into bins of
the same number of events. A bin is suspect when

- its flow rate is far from the acquisition's usual one;
- the running median, over five bins, of the bin medians of any population
  of any parameter lies more than `mad` MADs from that population's median
  over the whole acquisition (the signal: steps, bubbles, anything abrupt);
- or the running median of the same bin medians over a long window lies
  far from the level the population keeps for longest (the drift: a
  moderate shift that lasts).

Suspect bins are removed, and so is a stretch of fewer than `consecutive`
good bins between two removed ones.

What follows from that:

- A parameter is followed one population at a time: it is cut at the
  valleys of its histogram over the whole file, and each bin gives the
  median of its events between two cuts. The median of all its events will
  not do for a marker that half of the cells carry: it jumps from one
  population to the other with the bin's share of each, and a quarter of
  a clean file was removed that way (two populations of 48 and 52%).
  How many cells a population has is not judged, only where it lies.
- The reference of the signal is flat: what most of the acquisition looked
  like. The running median (five bins) is what keeps one or two noisy bins
  among good ones: they do not move it. A defect has to last three bins.
- The drift is what five bins cannot show. A shift of a third of a
  population's spread is within the noise of a bin of 150 to 500 events,
  and plain over 10,000: the window of the slow pass (never more than a
  quarter of the bins). Its reference is the level of the longest stable
  stretch — the middle of the narrowest range of values that holds half of
  the bins — and not the median, which a shift over a third of the file
  pulls towards itself (by the median, a real file whose shift lasts a
  third of it lost 18%). A stretch is a drift when the slow level goes
  further from that level than both `mad` MADs of the slow level itself
  and `drift` times the spread of the parameter's events. The first limit
  is the file's own: an acquisition that wanders all along has a wide one,
  and loses nothing. The second keeps a file whose level hardly moves from
  losing the bins where it moved most: without it, a quarter of two real
  files and a sixth of a third was removed, as PeacoQC removes it, where
  no parameter's median moves by 0.05 of its spread. What is removed of a
  drift is all of the departure, from where the slow level is seven tenths
  of the way to the limit: the part beyond the limit alone was a fifth of
  that real file in which the shift lasts a third.
- Either way a defect is found only while less than half of the
  acquisition is affected. With more, the bad part is the reference. And
  a slow ramp is found from where it nears the limit, not from where it
  began.
- The flow rate is not smoothed. A clog holds few events and therefore few
  bins, however long it lasts, and its bins are removed on their own.
- A population too small to give a median in most bins (fewer than 20
  events in half of them) is not followed, and a defect that touches only
  it is not seen. One that is followed and then missing from a bin is
  suspect there.

How this differs from the programs it takes after (it was run beside both:
`docs/comparacao.md`, section 3). PeacoQC finds the density peaks of each
parameter in every bin and groups them into trajectories, where this cuts
the parameter once for the whole file; it judges one trajectory, smoothed
by a spline, with the median and the MAD of the smoothed curve, where this
has two passes and a floor under the second; it adds an isolation tree, its
bins overlap by half, and it does not look at the flow rate. flowAI judges
the flow rate in bins of time after removing its trend, finds the signal's
changes by changepoint detection on the bin medians and keeps the longest
stable stretch, and has a third check, of the dynamic range.

The defaults: `mad` 6 and `consecutive` 5 are PeacoQC's. The bins are not:
here they hold 1/500 of the file and at least 150 events, so 150 to 500
for files of up to 250,000 events; PeacoQC's hold twice 1/500 of the file
rounded up to the next 500 (150 is only its floor, never reached), so 500
or 1,000 for the same files. `drift` 0.25 and the window were chosen on
the 32 real files and the simulated defects of that comparison, not from
a principle.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import find_peaks

from .gating import TimeGate

#: a bin's median has to mean something, and 500 bins draw a file of any size
MIN_BIN, MAX_BINS = 150, 500
#: fewer bins than this and there is no "usual" to compare a bin with
MIN_BINS = 20
#: bins in the running median: two noisy ones in a row do not move it
SMOOTH = 5
#: a flow rate within this factor of the usual one is never a defect,
#: however steady the rest of the acquisition was
RATE_FACTOR = 1.5
#: Time may step back by this share of its span (events written slightly
#: out of order) before the file is taken for several acquisitions in one
BACKWARDS = 0.01

#: a population is followed when half of the bins hold this many of its
#: events, and is taken to be missing from a bin that holds fewer than FEW
FOLLOWED, FEW = 20, 5
#: the histogram a parameter's populations are told apart on: its bins,
#: the Gaussian it is smoothed with (in bins), and how far a peak has to
#: rise above the valley beside it, as a share of the tallest
GRID, BLUR, PROMINENCE = 64, 1.5, 0.15

#: the slow pass: the events its running median is taken over, and the
#: share of the bins that may be at most
SLOW_EVENTS, SLOW_SHARE = 10000, 0.25
#: a drift, once found, is taken from where the slow level is this share
#: of the way to the limit
CARRY = 0.7

RATE, SIGNAL, DRIFT, ISLANDS = "flow rate", "signal", "drift", "islands"
#: the population quality control makes, above every other
POPULATION = "Clean"
#: the default of `check`'s `drift`
DRIFT_LIMIT = 0.25


@dataclass
class Bins:
    """The acquisition a bin at a time: what the plots are drawn from.
    Times are in the file's units, rates in events per one of them."""
    start: np.ndarray
    end: np.ndarray
    events: np.ndarray
    rate: np.ndarray
    #: lowest and highest flow rate that is not suspect
    rate_limits: tuple[float, float]
    #: bins x populations, each on its parameter's scale; infinite where
    #: the population is missing from the bin
    medians: np.ndarray
    #: for each column of `medians`, which of the parameters it is one of
    channel: np.ndarray
    #: populations x 2: the lowest and highest running median not suspect
    limits: np.ndarray
    bad_rate: np.ndarray
    #: bins x populations
    bad_signal: np.ndarray
    #: bins x populations: the running median over `window` bins, which
    #: the drift is judged on
    slow: np.ndarray
    window: int
    #: populations x 2: how far `slow` may go from the stable level
    slow_limits: np.ndarray
    #: bins x populations
    bad_drift: np.ndarray
    #: what is removed in the end, the islands included
    bad: np.ndarray


@dataclass
class Result:
    #: False where the file could not be checked; `reason` says why
    checked: bool
    #: per event, in the file's order: True for the ones to analyse
    keep: np.ndarray = field(repr=False)
    #: (first, last) Time of every stretch to remove, in the file's units
    intervals: list[tuple[float, float]] = field(default_factory=list)
    #: events in the bins each check found: RATE, SIGNAL, DRIFT, ISLANDS.
    #: A bin found by two checks is counted by both.
    counts: dict[str, int] = field(default_factory=dict)
    #: share of the events removed
    removed: float = 0.0
    summary: str = ""
    reason: str = ""
    channels: list[str] = field(default_factory=list)
    #: a name for each population followed: the parameter's, and where the
    #: population lies when the parameter has several
    populations: list[str] = field(default_factory=list)
    #: the scale each was judged on, None for the values as stored
    scales: list = field(default_factory=list, repr=False)
    #: the parameters with a suspect population anywhere, the worst first
    blamed: list[str] = field(default_factory=list)
    #: what found each stretch of `intervals`: RATE, SIGNAL and DRIFT, as
    #: many as found a bin of it
    found: list[list[str]] = field(default_factory=list)
    bins: Bins | None = field(default=None, repr=False)

    def gate(self, dim: str) -> TimeGate:
        """The gate that keeps what this check keeps."""
        return TimeGate([dim], [list(i) for i in self.intervals], inside=False)


def time_channel(channels: list[str]) -> str | None:
    """The parameter that holds the time of acquisition, if there is one."""
    return next((c for c in channels if c.strip().lower() == "time"), None)


def unchecked(reason: str, n_events: int) -> Result:
    return Result(False, np.ones(n_events, dtype=bool), reason=reason,
                  summary=f"Cannot be checked: {reason}.")


def _unusable(time, n_events: int, size: int) -> str:
    """Why this Time cannot order the events, or "" where it can."""
    if time is None:
        return "the file has no Time parameter"
    if n_events < MIN_BINS * size:
        return f"{n_events:,d} events are too few (it takes {MIN_BINS * size:,d})"
    if not np.isfinite(time).all():
        return "Time is not a number for some events"
    span = float(time.max() - time.min())
    if span <= 0:
        return "Time is the same for every event"
    if float(-np.diff(time).min()) > BACKWARDS * span:
        return "Time runs backwards: several acquisitions in one file, or a clock that was reset"
    return ""


def _spread(x: np.ndarray) -> float:
    """The MAD, scaled to a standard deviation; the mean deviation where
    more than half of the values are the same and the MAD is zero."""
    x = x[np.isfinite(x)]
    away = np.abs(x - np.median(x))
    return float(1.4826 * np.median(away) or 1.2533 * away.mean())


def _spread_about(x: np.ndarray, centre: float) -> float:
    """The MAD about a centre that is given, scaled to a standard deviation."""
    return float(1.4826 * np.median(np.abs(x[np.isfinite(x)] - centre)))


def _valleys(v: np.ndarray) -> np.ndarray:
    """Where a parameter's populations part: the lowest point between
    every two peaks of its histogram that stand clear of each other.
    Empty for a parameter with one population."""
    lo, hi = np.percentile(v, [0.5, 99.5])
    if not hi > lo:
        return np.array([])
    heights, edges = np.histogram(v, bins=GRID, range=(lo, hi))
    # a row of zeros at each end, so a population at the edge is a peak too
    smooth = np.pad(gaussian_filter1d(heights.astype(float), BLUR, mode="constant"), 1)
    peaks, _ = find_peaks(smooth, prominence=PROMINENCE * smooth.max())
    lowest = [a + int(np.argmin(smooth[a:b + 1])) for a, b in zip(peaks, peaks[1:])]
    centres = (edges[:-1] + edges[1:]) / 2
    return centres[np.array(lowest, dtype=int) - 1]


def _follow(v: np.ndarray, edges: np.ndarray) -> list[np.ndarray]:
    """The median of each population of one parameter in every bin: one
    array for each population that at least half of the bins hold enough
    of. One seen in fewer is not a population of the sample — it may be
    what a defect made of another — and is not followed."""
    which = np.searchsorted(_valleys(v), v)
    out = []
    for p in range(int(which.max()) + 1):
        inside = which == p
        held = np.add.reduceat(inside, edges[:-1])
        if np.median(held) < FOLLOWED:
            continue
        out.append(np.array([np.median(v[a:b][inside[a:b]]) if n >= FEW else np.inf
                             for a, b, n in zip(edges[:-1], edges[1:], held)]))
    return out


def _level(slow: np.ndarray) -> float:
    """Where a population stays for longest: the middle of the narrowest
    range of values that holds half of the bins. Not the median of them
    all, which a drift over a third of the acquisition pulls a third of
    the way towards itself."""
    v = np.sort(slow[np.isfinite(slow)])
    half = (len(v) + 1) // 2
    first = int(np.argmin(v[half - 1:] - v[:len(v) - half + 1]))
    return float(np.median(v[first:first + half]))


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """(first, one past the last) of every run of True."""
    edges = np.flatnonzero(np.diff(np.concatenate(([0], mask.astype(np.int8), [0]))))
    return [(int(a), int(b)) for a, b in zip(edges[::2], edges[1::2])]


def _with_islands(flagged: np.ndarray, consecutive: int) -> np.ndarray:
    bad = flagged.copy()
    runs = _runs(flagged)
    for (_, b), (a, _) in zip(runs, runs[1:]):
        if a - b < consecutive:
            bad[b:a] = True
    return bad


def check(events: np.ndarray, time: np.ndarray | None, channels: list[str],
          scales: list | None = None, mad: float = 6.0, consecutive: int = 5,
          events_per_bin: int | None = None, *, drift: float = DRIFT_LIMIT) -> Result:
    """Check one acquisition.

    `events` holds one column for each of `channels`, the parameters to
    judge, and `scales` the scale each is drawn on (None: as stored) —
    fluorescence is judged where it is looked at, not in raw units where
    the brightest events would decide. `time` is the Time parameter, per
    event, in the file's order.

    `mad`: how many MADs from the acquisition's median a bin has to lie to
    be suspect; larger removes less. `consecutive`: good bins between two
    removed stretches are kept only when there are at least this many.
    `events_per_bin`: by default 1/500 of the file and at least 150.
    `drift`: the smallest move of a population's slow level away from
    where it stays for longest that is taken for a drift, as a share of
    the spread of the parameter's events; larger removes less, and 0
    leaves the drift unjudged.
    """
    n = len(events)
    size = int(events_per_bin or max(MIN_BIN, -(-n // MAX_BINS)))
    why = _unusable(time, n, size)
    if why:
        return unchecked(why, n)
    time = np.asarray(time, dtype=np.float64)
    order = np.argsort(time, kind="stable")
    t = time[order]
    edges = np.append(np.arange(n // size) * size, n)   # the last bin takes the remainder
    start, end, count = t[edges[:-1]], t[edges[1:] - 1], np.diff(edges)

    # a bin lasts until the next one begins, so a pause belongs to a bin;
    # and never less than the clock's own step, where it is coarse
    steps = np.diff(t)
    lasted = np.maximum(np.append(np.diff(start), end[-1] - start[-1]), steps[steps > 0].min())
    rate = count / lasted
    log_rate = np.log(rate)
    usual = float(np.median(log_rate))
    reach = max(mad * _spread(log_rate), float(np.log(RATE_FACTOR)))
    bad_rate = np.abs(log_rate - usual) > reach

    scales = list(scales) if scales is not None else [None] * len(channels)
    tracks, channel, names, width = [], [], [], []
    for j, scale in enumerate(scales):
        v = np.asarray(events[:, j], dtype=np.float64)[order]
        v = scale(v) if scale is not None else v
        followed = _follow(v, edges)
        width += [_spread(v)] * len(followed)
        for track in followed:
            where = float(np.median(track[np.isfinite(track)]))
            where = float(scale.inverse(where)) if scale is not None else where
            names.append(channels[j] if len(followed) == 1
                         else f"{channels[j]}, around {float(f'{where:.2g}'):,.6g}")
        tracks += followed
        channel += [j] * len(followed)
    medians = np.column_stack(tracks) if tracks else np.empty((len(count), 0))
    channel = np.array(channel, dtype=int)
    smooth = median_filter(medians, size=(SMOOTH, 1), mode="mirror")
    centre = np.array([np.median(c[np.isfinite(c)]) for c in smooth.T])
    spread = np.array([_spread(c) for c in smooth.T])
    # a population missing from three bins of five is as wrong as one that moved
    bad_signal = ~np.isfinite(smooth) | ((np.abs(smooth - centre) > mad * spread) & (spread > 0))

    # the slow pass: the same medians over a window long enough for a
    # shift of a quarter of a spread to stand clear of the bins' noise
    window = max(SMOOTH, min(round(SLOW_EVENTS / size), int(SLOW_SHARE * len(count)))) | 1
    slow = median_filter(medians, size=(window, 1), mode="mirror")
    level = np.array([_level(c) for c in slow.T])
    wander = np.array([_spread_about(c, at) for c, at in zip(slow.T, level)])
    reach_slow = np.maximum(mad * wander, drift * np.array(width))
    away = np.where(np.isfinite(slow) & (reach_slow > 0) & (drift > 0),
                    np.abs(slow - level), 0.0)
    bad_drift = away > reach_slow
    # the whole of a departure that went beyond the limit, not only its peak
    for j, near in enumerate((away > CARRY * reach_slow).T):
        for a, b in _runs(near):
            bad_drift[a:b, j] = bad_drift[a:b, j].any()

    flagged = bad_rate | bad_signal.any(axis=1) | bad_drift.any(axis=1)
    bad = _with_islands(flagged, consecutive)
    intervals = [(float(start[a]), float(end[b - 1])) for a, b in _runs(bad)]
    result = Result(True, np.ones(n, dtype=bool), intervals, channels=list(channels),
                    populations=names, scales=scales)
    # by the gate itself, so that what is reported is what the gate removes:
    # events of a neighbouring bin on the very tick an interval ends at go too
    result.keep = result.gate("Time").contains(lambda _: time)
    result.removed = float(1.0 - result.keep.mean())
    by = {RATE: bad_rate, SIGNAL: bad_signal.any(axis=1), DRIFT: bad_drift.any(axis=1)}
    result.counts = {**{what: int(count[where].sum()) for what, where in by.items()},
                     ISLANDS: int(count[bad & ~flagged].sum())}
    result.found = [[what for what, where in by.items() if where[a:b].any()]
                    for a, b in _runs(bad)]
    result.blamed = _blamed(bad_signal | bad_drift, channel, channels)
    result.bins = Bins(start, end, count, rate,
                       (float(np.exp(usual - reach)), float(np.exp(usual + reach))),
                       medians, channel,
                       np.column_stack([centre - mad * spread, centre + mad * spread]),
                       bad_rate, bad_signal, slow, window,
                       np.column_stack([level - reach_slow, level + reach_slow]),
                       bad_drift, bad)
    result.summary = _summary(result, n, size)
    return result


def _blamed(bad: np.ndarray, channel: np.ndarray, channels: list[str]) -> list[str]:
    """The parameters with a population among `bad` (bins x populations),
    the one with most bins first."""
    worst = np.array([bad[:, channel == j].any(axis=1).sum() for j in range(len(channels))])
    return [channels[j] for j in np.argsort(-worst, kind="stable") if worst[j]]


def _named(names: list[str], most: int = 3) -> str:
    return ", ".join(names[:most]) + ("…" if len(names) > most else "")


def _summary(r: Result, n: int, size: int) -> str:
    if not r.intervals:
        return (f"Nothing to remove: the flow rate and the signal of {len(r.channels)} "
                f"parameter(s) are steady over {len(r.bins.bad)} bins of {size:,d} events.")
    found = []
    if r.counts[RATE]:
        found.append(f"flow rate {r.counts[RATE]:,d} events")
    for what, bad in ((SIGNAL, r.bins.bad_signal), (DRIFT, r.bins.bad_drift)):
        if r.counts[what]:
            found.append(f"{what} {r.counts[what]:,d} "
                         f"({_named(_blamed(bad, r.bins.channel, r.channels))})")
    if r.counts[ISLANDS]:
        found.append(f"short stretches between them {r.counts[ISLANDS]:,d}")
    text = (f"To remove: {100 * r.removed:.1f}% of {n:,d} events, in {len(r.intervals)} "
            f"stretch(es) of time — {'; '.join(found)}.")
    if r.removed > 0.2:
        text += " That is more than a fifth of the file: look at the plots before applying."
    return text
