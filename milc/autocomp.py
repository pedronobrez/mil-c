"""
A spillover matrix computed from single-stain controls: one FCS file per
fluorochrome, each stained with that fluorochrome alone, and optionally
an unstained one.

Row i of a spillover matrix is how much fluorochrome i shows in every
detector (see `compensation.py`), so a control stained with fluorochrome
i alone measures row i: in its events, every detector is a straight line
against the control's own — its primary — detector, and the slope of
that line is the spillover.

Two ways of measuring it:

`classic` — what cytometer software has always done: a positive and a
negative population of the control, and per detector the difference of
their medians over the same difference in the primary detector. The
positive population is what lies above the valley of the primary
detector's histogram; the negative, the unstained control's events where
one is given, else the control's own events below the valley.

`autospill` — the method of Roca et al. (2021, Nat Commun 12: 2890,
"AutoSpill is a principled framework that simplifies the analysis of
multichromatic flow cytometry data"), written from the authors' R package
(`carlosproca/autospill` 0.2.0): a robust linear regression (Huber's M
estimator, as `MASS::rlm`) of each detector on the primary one over all
the control's events gives a first matrix; then the controls are
compensated with it, the same regression on the compensated events gives
what spill is left, and the matrix is corrected by it, again and again,
first on the linear scale until no slope left is above 1e-2, then on the
biexponential one until none is above 1e-4 or nothing improves any more.
Their thresholds, step sizes and stopping rules are kept as they are.

Nothing here knows about the window: `compute()` is called by the
compensation dialog through the worker and by `milc autocomp`.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .compensation import Compensation
from .fcs import read_fcs
from .transforms import Biex

CLASSIC, AUTOSPILL = "classic", "autospill"
METHODS = {CLASSIC: "classic: positive against negative population",
           AUTOSPILL: "AutoSpill: robust regression, refined until no spill is left"}

#: words in a file name that say a control is the unstained one
_UNSTAINED = re.compile(r"(unstain|un-stain|no[ _-]?stain|blank|negative[ _-]?control|"
                        r"(^|[^a-z])us([^a-z]|$)|(^|[^a-z])neg([^a-z]|$))", re.I)


# --- the controls --------------------------------------------------------------

@dataclass
class Control:
    path: str
    #: events x detectors, as acquired (not compensated), after the scatter gate
    events: np.ndarray
    #: how many events the file has, and how many are used
    total: int
    #: the detector this control's fluorochrome is read in; None: unstained
    primary: str | None = None
    unstained: bool = False
    #: share of the file's events inside the scatter gate (None: no gate)
    gated: float | None = None
    #: why the gate was not applied, where it was asked for and could not be
    gate_note: str = ""

    @property
    def name(self) -> str:
        return os.path.splitext(os.path.basename(self.path))[0]


def fluorescence(channels: list[str], keywords: dict[str, str] | None = None) -> list[str]:
    """The detectors a matrix is for: those of the file's own `$SPILL`,
    else every channel but scatter and time (the areas, where a file has
    areas and heights)."""
    from .compensation import from_keywords
    own = from_keywords(keywords or {})
    if own is not None:
        return list(own.channels)
    skip = re.compile(r"^(FSC|SSC|Time|TIME|Event|EVENT)", re.I)
    out = [c for c in channels if not skip.match(c)]
    areas = [c for c in out if c.endswith("-A")]
    return areas or out


def looks_unstained(path: str) -> bool:
    return bool(_UNSTAINED.search(os.path.splitext(os.path.basename(path))[0]))


def _scatter(channels: list[str]) -> tuple[int, int] | None:
    def first(prefix):
        for suffix in ("-A", "-H", ""):
            for i, c in enumerate(channels):
                if c.upper().startswith(prefix) and c.upper().endswith(suffix):
                    return i
        return None
    f, s = first("FSC"), first("SSC")
    return (f, s) if f is not None and s is not None else None


def scatter_gate(events: np.ndarray, channels: list[str]) -> tuple[np.ndarray | None, str]:
    """The main population on forward against side scatter, as the robust
    Gaussian ellipse of `autogate.gaussian` that holds 95% of it: debris
    and clumps out. (None, why) where there is no scatter or no population."""
    from . import autogate
    which = _scatter(channels)
    if which is None:
        return None, "no forward and side scatter"
    x = events[:, which[0]].astype(np.float64)
    y = events[:, which[1]].astype(np.float64)
    xr = tuple(np.percentile(x, [0.1, 99.9]))
    yr = tuple(np.percentile(y, [0.1, 99.9]))
    if xr[1] <= xr[0] or yr[1] <= yr[0]:
        return None, "the scatter has no spread"
    g = autogate.gaussian(x, y, xr, yr, share=0.95)
    if g.center is None:
        return None, g.why
    dx, dy = x - g.center[0], y - g.center[1]
    c, s = np.cos(g.angle), np.sin(g.angle)
    u, v = dx * c + dy * s, -dx * s + dy * c
    return (u / g.radii[0]) ** 2 + (v / g.radii[1]) ** 2 <= 1.0, ""


def load_controls(paths: list[str], detectors: list[str], gate: bool = False,
                  max_events: int | None = 100_000, seed: int = 0) -> list[Control]:
    """Read the controls: the detectors' columns as acquired, inside the
    scatter gate where `gate`, at most `max_events` of them (a random
    choice, the same every time). A file that lacks a detector is an
    error that names it."""
    out = []
    for path in paths:
        f = read_fcs(path)
        missing = [d for d in detectors if d not in f.channels]
        if missing:
            raise ValueError(f"{os.path.basename(path)} has no {', '.join(missing)}")
        events = f.events
        kept, note = None, ""
        if gate:
            mask, note = scatter_gate(events, f.channels)
            if mask is not None and mask.sum() >= 100:
                kept = float(mask.mean())
                events = events[mask]
            elif mask is not None:
                note = "fewer than 100 events in the scatter gate"
        cols = [f.index(d) for d in detectors]
        x = np.asarray(events[:, cols], dtype=np.float64)
        if max_events and len(x) > max_events:
            pick = np.random.default_rng(seed).choice(len(x), max_events, replace=False)
            x = x[np.sort(pick)]
        out.append(Control(path, x, len(f), unstained=looks_unstained(path), gated=kept,
                           gate_note=note))
    return out


# --- which detector is each control's --------------------------------------------

def brightness(controls: list[Control], detectors: list[str]) -> np.ndarray:
    """controls x detectors: how far the brightest events of a control
    (its 99.5th percentile) stand above the negative events (the
    unstained control's median where there is one, else the control's own
    20th percentile), in each detector's units."""
    unstained = [c for c in controls if c.unstained]
    base = (np.median(unstained[0].events, axis=0) if unstained else None)
    out = np.zeros((len(controls), len(detectors)))
    for i, c in enumerate(controls):
        low = base if base is not None else np.percentile(c.events, 20, axis=0)
        out[i] = np.percentile(c.events, 99.5, axis=0) - low
    return out


def assign(controls: list[Control], detectors: list[str]) -> list[str | None]:
    """Each stained control's primary detector: the one it is brightest
    in, no two controls on the same detector — where two would be, the
    choice that keeps the most of each one's own brightness (an
    assignment problem over brightness relative to each control's
    brightest). Sets `primary` and returns them; None for the unstained."""
    from scipy.optimize import linear_sum_assignment
    stained = [i for i, c in enumerate(controls) if not c.unstained]
    if len(stained) > len(detectors):
        raise ValueError(f"{len(stained)} stained controls for {len(detectors)} detectors")
    b = brightness(controls, detectors)
    rel = b[stained] / np.maximum(b[stained].max(axis=1, keepdims=True), 1e-12)
    rows, cols = linear_sum_assignment(-rel)
    for c in controls:
        c.primary = None
    for r, j in zip(rows, cols):
        controls[stained[r]].primary = detectors[j]
    return [c.primary for c in controls]


# --- robust regression --------------------------------------------------------------

def huber(x: np.ndarray, y: np.ndarray, k: float = 1.345, max_iter: int = 100,
          acc: float = 1e-4) -> tuple[float, float, bool]:
    """(intercept, slope, converged) of y on x by Huber's M estimator, as
    `MASS::rlm(y ~ x)` computes it by default: least squares to start,
    then weights min(1, k / |r / s|) with s the residuals' MAD / 0.6745,
    until the residuals move by less than `acc` of their size."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if len(x) < 3:
        return 0.0, 0.0, False

    def fit(w):
        sw = w.sum()
        mx, my = (w * x).sum() / sw, (w * y).sum() / sw
        dx = x - mx
        sxx = (w * dx * dx).sum()
        b = (w * dx * (y - my)).sum() / sxx if sxx > 0 else 0.0
        return my - b * mx, b

    a, b = fit(np.ones_like(x))
    r = y - a - b * x
    for _ in range(max_iter):
        s = np.median(np.abs(r)) / 0.6745
        if s == 0:
            return a, b, True
        u = np.abs(r) / s
        w = np.minimum(1.0, k / np.maximum(u, 1e-300))
        a, b = fit(w)
        new = y - a - b * x
        change = np.sqrt(((r - new) ** 2).sum() / max(1e-20, (r ** 2).sum()))
        r = new
        if change <= acc:
            return a, b, True
    return a, b, False


def _trim_bounds(v: np.ndarray, share: float) -> tuple[float, float]:
    """The values `round(n * share)` from each end, as AutoSpill takes them."""
    n = max(int(round(len(v) * share)), 1)
    s = np.sort(v)
    return s[n - 1], s[len(s) - n]


def _slopes(x: np.ndarray, primary: int, trim: float = 0.01,
            transform=None) -> tuple[np.ndarray, int]:
    """The slope of every detector on the primary one over one control's
    events, on the linear scale — or, with `transform`, fitted on that
    scale and read back as the linear slope between the ends of the
    fitted line (AutoSpill's `get.compensation.error`). Events in the top
    and bottom `trim` of either detector are left out of each fit.
    Returns the slopes and how many fits did not converge."""
    d = x.shape[1]
    data = transform(x) if transform is not None else x
    p = data[:, primary]
    plo, phi = _trim_bounds(p, trim)
    pin = (p > plo) & (p < phi)
    out = np.zeros(d)
    failed = 0
    for j in range(d):
        if j == primary:
            out[j] = 1.0
            continue
        v = data[:, j]
        lo, hi = _trim_bounds(v, trim)
        keep = pin & (v > lo) & (v < hi)
        if keep.sum() < 10:
            failed += 1
            continue
        a, b, ok = huber(p[keep], v[keep])
        if not ok:
            # as the package does: least squares where rlm does not converge
            failed += 1
            xx, yy = p[keep], v[keep]
            b = np.polyfit(xx, yy, 1)[0]
            a = yy.mean() - b * xx.mean()
        if transform is None:
            out[j] = b
        else:
            y1, y2 = p[keep].min(), p[keep].max()
            x1, x2 = a + b * y1, a + b * y2
            if y1 == y2 or x1 == x2:
                out[j] = 0.0
            else:
                inv = transform.inverse
                out[j] = (inv(x2) - inv(x1)) / (inv(y2) - inv(y1))
    return out, failed


# --- the result ----------------------------------------------------------------------

@dataclass
class Step:
    iteration: int
    scale: str          # "linear" or "biexponential"
    step: float         # the share of the correction taken (AutoSpill's lambda)
    #: the spread (standard deviation) and the largest of the slopes left
    spread: float
    largest: float


@dataclass
class Result:
    method: str
    #: the detectors, the matrix's rows and columns
    channels: list[str]
    #: square over `channels`; rows of detectors without a control are the base's
    matrix: np.ndarray
    #: detector → the control file it was computed from
    controls: dict[str, str]
    #: rows (fluorochromes with a control) x channels: the slopes still
    #: left in the controls once compensated with `matrix`, on the linear scale
    residual: np.ndarray | None = None
    #: the same, compensated with the matrix the refinement started from
    residual_start: np.ndarray | None = None
    history: list[Step] = field(default_factory=list)
    converged: bool = True
    notes: list[str] = field(default_factory=list)

    def compensation(self) -> Compensation:
        return Compensation(list(self.channels), self.matrix.copy(), name="controls")

    def largest_residual(self) -> float | None:
        return None if self.residual is None else float(np.abs(self.residual).max())


def _rows(controls: list[Control], detectors: list[str]) -> list[tuple[int, Control]]:
    rows = []
    seen = set()
    for c in controls:
        if c.unstained or c.primary is None:
            continue
        if c.primary not in detectors:
            raise ValueError(f"{c.name}: {c.primary} is not one of the detectors")
        if c.primary in seen:
            raise ValueError(f"two controls are given {c.primary}")
        seen.add(c.primary)
        rows.append((detectors.index(c.primary), c))
    if not rows:
        raise ValueError("no stained control has a detector")
    return rows


def residual_spill(controls: list[Control], detectors: list[str], matrix: np.ndarray,
                   trim: float = 0.01) -> tuple[np.ndarray, list[str]]:
    """What spill each control still shows once compensated with `matrix`:
    per control, the robust slope of every detector on its primary one,
    the primary's own 1 taken away. Rows in the order of the controls with
    a detector; also returns those detectors."""
    rows = _rows(controls, detectors)
    inv = np.linalg.inv(matrix)
    out = np.zeros((len(rows), len(detectors)))
    for k, (i, c) in enumerate(rows):
        s, _ = _slopes(c.events @ inv, i, trim)
        s[i] -= 1.0
        out[k] = s
    return out, [detectors[i] for i, _ in rows]


def _base(detectors: list[str], base: Compensation | None) -> np.ndarray:
    S = np.eye(len(detectors))
    if base is not None:
        for i, a in enumerate(detectors):
            for j, b in enumerate(detectors):
                if a in base.channels and b in base.channels:
                    S[i, j] = base.matrix[base.channels.index(a), base.channels.index(b)]
    return S


# --- classic ------------------------------------------------------------------------

def _split(p: np.ndarray, transform) -> float | None:
    """The valley between the negative and positive events of a primary
    detector, on the biexponential scale: two means first, then the lowest
    point of the smoothed histogram between them (AutoSpill's
    `get.marker.spillover.posnegpop`). None where they are not two."""
    t = transform(p)
    lo, hi = np.percentile(t, [10, 90])
    if hi - lo < 1e-6:
        return None
    for _ in range(100):
        cut = (lo + hi) / 2
        a, b = t[t <= cut], t[t > cut]
        if len(a) == 0 or len(b) == 0:
            return None
        nlo, nhi = a.mean(), b.mean()
        if abs(nlo - lo) + abs(nhi - hi) < 1e-9:
            break
        lo, hi = nlo, nhi
    edges = np.linspace(np.percentile(t, 0.5), np.percentile(t, 99.5), 257)
    h, _ = np.histogram(t, edges)
    centres = (edges[:-1] + edges[1:]) / 2
    sd = t.std()
    bw = 1.06 * sd * len(t) ** -0.2 / (edges[1] - edges[0])   # Silverman, in bins
    k = np.arange(-int(4 * bw) - 1, int(4 * bw) + 2)
    kernel = np.exp(-0.5 * (k / max(bw, 0.5)) ** 2)
    smooth = np.convolve(h, kernel / kernel.sum(), mode="same")
    between = (centres > lo) & (centres < hi)
    if not between.any():
        return None
    valley = centres[between][np.argmin(smooth[between])]
    # a valley as high as the lower peak is no valley
    peak_lo = smooth[np.argmin(np.abs(centres - lo))]
    peak_hi = smooth[np.argmin(np.abs(centres - hi))]
    if smooth[between].min() > 0.8 * min(peak_lo, peak_hi):
        return None
    return float(transform.inverse(valley))


def classic(controls: list[Control], detectors: list[str], base: Compensation | None = None,
            min_events: int = 50, universal_negative: bool = False) -> Result:
    """Each control's row: (median of the positive events − median of the
    negative ones) in every detector, over the same in its primary.

    The negative events are the control's own, below the valley — the
    same cells, the same tube, the same autofluorescence — unless
    `universal_negative`, or unless the control has no negative events
    (beads, often): then the unstained control's, where there is one."""
    rows = _rows(controls, detectors)
    unstained = [c for c in controls if c.unstained]
    neg_un = unstained[0].events if unstained else None
    t = Biex()
    S = _base(detectors, base)
    notes = []
    for i, c in rows:
        p = c.events[:, i]
        cut = _split(p, t)
        if cut is None and neg_un is not None:
            # every event stained (beads, often): above the unstained's own
            cut = float(np.percentile(neg_un[:, i], 99.5))
            notes.append(f"{c.name}: one population only; positive above the unstained "
                         "control's 99.5th percentile")
        if cut is None:
            raise ValueError(f"{c.name}: no negative and positive populations in {detectors[i]}, "
                             "and no unstained control to compare with")
        plo, phi = _trim_bounds(p, 0.01)
        pos = c.events[(p > cut) & (p < phi)]
        own = c.events[(p < cut) & (p > plo)]
        if neg_un is not None and (universal_negative or len(own) < min_events):
            neg = neg_un
        else:
            neg = own
        if len(pos) < min_events or len(neg) < min_events:
            raise ValueError(f"{c.name}: {len(pos)} positive and {len(neg)} negative events "
                             f"in {detectors[i]}; at least {min_events} of each are needed")
        diff = np.median(pos, axis=0) - np.median(neg, axis=0)
        if diff[i] <= 0:
            raise ValueError(f"{c.name}: the positive events are not brighter in {detectors[i]}")
        S[i] = diff / diff[i]
    res, _ = residual_spill(controls, detectors, S)
    return Result(CLASSIC, list(detectors), S, {detectors[i]: c.path for i, c in rows},
                  residual=res, notes=notes)


# --- AutoSpill ----------------------------------------------------------------------

def autospill(controls: list[Control], detectors: list[str], base: Compensation | None = None,
              max_iter: int = 100, tolerance_linear: float = 1e-2,
              tolerance_biex: float = 1e-4, tolerance_change: float = 1e-6,
              coarse: float = 1.0, fine: float = 0.1, history: int = 10, trim: float = 0.01,
              progress: Callable[[int, int, str], None] | None = None) -> Result:
    """AutoSpill's two steps: `get.marker.spillover` (robust slopes on the
    events as acquired) and `refine.spillover` (correct by what is left
    once compensated, until it is below the tolerances). The defaults are
    the package's `get.autospill.param("minimal")`."""
    rows = _rows(controls, detectors)
    idx = [i for i, _ in rows]
    n = len(detectors)
    S0 = _base(detectors, base)
    t = Biex()
    # the first matrix: robust slopes on the events as acquired
    first = S0.copy()
    failed = 0
    for i, c in rows:
        s, f = _slopes(c.events, i, trim)
        first[i] = s
        failed += f
    notes = []
    if failed:
        notes.append(f"{failed} robust fit(s) did not converge; least squares used there")

    def error(S, transform):
        inv = np.linalg.inv(S)
        out = np.zeros((len(rows), n))
        for k, (i, c) in enumerate(rows):
            s, _ = _slopes(c.events @ inv, i, trim, transform)
            s[i] -= 1.0
            out[k] = s
        return out

    linear = True
    step = coarse
    delta = -1.0
    threshold = tolerance_linear
    hist = [-1.0] * history
    # as the package starts: the identity plus (first - identity); the rows
    # of detectors without a control are the base's and stay so throughout
    S = S0.copy()
    S[idx] = np.eye(n)[idx]
    update = np.zeros((n, n))
    update[idx] = first[idx] - S[idx]
    it, last, done, converged = 0, False, False, False
    steps: list[Step] = []
    start_residual = None
    while not done:
        S = S + update
        S = S / np.diag(S)[:, None]
        err = error(S, None if linear else t)
        if it == 0:
            start_residual = err if linear else None
        previous = delta
        delta = float(np.std(err, ddof=1)) if err.size > 1 else float(np.abs(err).max())
        largest = float(np.abs(err).max())
        hist[it % history] = delta - previous if previous >= 0 else -1.0
        change = float(np.mean(hist))
        steps.append(Step(it, "linear" if linear else "biexponential", step, delta, largest))
        if progress is not None:
            progress(min(it + 1, max_iter), max_iter,
                     f"iteration {it}, {'linear' if linear else 'biexponential'} scale, "
                     f"largest spill left {largest:.2g}")
        if linear and largest < threshold:
            linear, threshold, step, delta = False, tolerance_biex, coarse, -1.0
            hist, change = [-1.0] * history, -1.0
        if change > -tolerance_change and step == coarse:
            step, delta = fine, -1.0
            hist, change = [-1.0] * history, -1.0
        converged = (not linear) and (largest < threshold or change > -tolerance_change)
        done = ((converged and last) or (change > -tolerance_change and linear)
                or (not converged and it == max_iter) or it > max_iter)
        last = converged
        it += 1
        full = np.zeros((n, n))
        full[idx] = err
        update = step * (full @ S)
        update[[j for j in range(n) if j not in idx]] = 0.0
    if not converged:
        notes.append(f"the refinement did not converge in {max_iter} iterations; "
                     "the last matrix is given")
    res, _ = residual_spill(controls, detectors, S, trim)
    return Result(AUTOSPILL, list(detectors), S, {detectors[i]: c.path for i, c in rows},
                  residual=res, residual_start=start_residual, history=steps,
                  converged=converged, notes=notes)


def compute(paths: list[str], detectors: list[str], method: str = AUTOSPILL,
            primaries: dict[str, str | None] | None = None, gate: bool = False,
            base: Compensation | None = None, max_events: int | None = 100_000,
            progress: Callable[[int, int, str], None] | None = None,
            controls: list[Control] | None = None,
            universal_negative: bool = False) -> Result:
    """Everything from the files: read them, give each its detector — the
    brightest, except where `primaries` (file path → detector, or None for
    the unstained one) says otherwise — and compute the matrix."""
    if controls is None:
        if progress is not None:
            progress(0, 1, "reading the controls")
        controls = load_controls(paths, detectors, gate, max_events)
    assign(controls, detectors)
    for c in controls:
        if primaries and c.path in primaries:
            c.primary = primaries[c.path]
            c.unstained = c.primary is None
    if method == CLASSIC:
        result = classic(controls, detectors, base, universal_negative=universal_negative)
    elif method == AUTOSPILL:
        result = autospill(controls, detectors, base, progress=progress)
    else:
        raise ValueError(f"no method {method!r}")
    for c in controls:
        if c.gate_note:
            result.notes.append(f"{c.name}: not gated on scatter ({c.gate_note})")
    missing = [d for d in detectors if d not in result.controls]
    if missing:
        result.notes.append("no control for " + ", ".join(missing)
                            + (": those rows are the acquisition matrix's" if base is not None
                               else ": those rows are left uncompensated"))
    return result


def simulate(spill: np.ndarray, detectors: list[str], n: int = 20_000, seed: int = 0,
             positive: float = 0.6, unstained: bool = True, debris: float = 0.0,
             outliers: float = 0.0) -> list[tuple[str, np.ndarray, list[str]]]:
    """Invented single-stain controls with a known spillover matrix, for
    the tests and the comparisons: (name, events, channels) per control,
    FSC-A and SSC-A first, then the detectors. Each control has a negative
    and a positive population (`positive` of the cells, lognormal, a
    different brightness per fluorochrome); every cell has
    autofluorescence; the detectors add photon (Poisson-like) and
    electronic noise and clip at 262143. `debris`: a share of small
    events, low scatter, dim and unrelated to the stain; `outliers`: a
    share of events with random bright values in one detector."""
    rng = np.random.default_rng(seed)
    d = len(detectors)
    autof = 120.0 * 0.6 ** np.arange(d)                 # brightest in the first detectors
    names = [(f"{det.rsplit('-', 1)[0]} stained control", k) for k, det in enumerate(detectors)]
    if unstained:
        names.append(("Unstained control", None))
    out = []
    for name, k in names:
        cells = rng.lognormal(0.0, 0.3, n)
        true = np.zeros((n, d))
        if k is not None:
            pos = rng.random(n) < positive
            level = 4000.0 * 2.0 ** (k % 4)
            true[pos, k] = level * rng.lognormal(0.0, 0.6, pos.sum())
        signal = true @ spill + cells[:, None] * autof
        noisy = signal + rng.normal(0.0, 1.0, signal.shape) * np.sqrt(
            2.0 * np.maximum(signal, 0.0) + 25.0 ** 2)
        fsc = 60000.0 * cells + rng.normal(0, 4000, n)
        ssc = 25000.0 * cells * rng.lognormal(0, 0.15, n)
        if debris:
            m = rng.random(n) < debris
            fsc[m] = rng.uniform(1000, 12000, m.sum())
            ssc[m] = rng.uniform(500, 8000, m.sum())
            noisy[m] = rng.lognormal(5.0, 1.0, (m.sum(), d))
        if outliers:
            m = rng.random(n) < outliers
            where = rng.integers(0, d, m.sum())
            noisy[np.flatnonzero(m), where] = rng.uniform(5e4, 2.5e5, m.sum())
        events = np.column_stack([fsc, ssc, np.minimum(noisy, 262143.0)])
        out.append((name, events, ["FSC-A", "SSC-A"] + list(detectors)))
    return out


def write_csv(path: str, result: Result) -> None:
    """The matrix as a table, detectors in the first row and column, in
    the proportions a `$SPILL` keyword holds (not percentages)."""
    import csv
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow([""] + result.channels)
        for c, row in zip(result.channels, result.matrix):
            w.writerow([c] + [f"{v:.6g}" for v in row])
