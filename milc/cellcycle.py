"""
The cell cycle from the DNA content of each cell (propidium iodide, DAPI,
Hoechst…), on a linear scale: a peak of cells with one complement of DNA
(G0/G1), a peak with two (G2/M) and, between them, the cells that are
copying it (S).

G0/G1 and G2/M are Gaussians. The mean of G2/M is held between two
multiples of that of G0/G1 (1.9 to 2.1 unless told otherwise: a dye is
seldom exactly proportional to DNA), and its CV is its own or, where
asked, that of G0/G1. S is found in one of two ways:

  Watson pragmatic   (after Watson, Chambers & Smith, Cytometry 1987)
      assumes nothing of the shape of S. G0/G1 is fitted on the side of
      its peak that S does not reach, the left, and G2/M on its right
      side; S is what is left of the histogram between the two means
      once both Gaussians are taken away. The S that lies under the
      sides fitted is allowed for: the level of S just clear of each
      peak, blurred as the peak is, is taken from the side before it is
      fitted again. Its RMS is that of the two sides fitted.
  Dean-Jett-Fox      (Dean & Jett 1974; Fox 1980) gives S a shape: a
      second-degree polynomial of the DNA content between the two means,
      never negative, broadened by the same measurement error as the
      peaks, the CV running from that of G0/G1 to that of G2/M. All three
      are fitted at once to the whole histogram. (Fox's extra Gaussian
      for a synchronised population in S is not fitted.)

Both are fitted between 4 SD below G0/G1 and 4 SD above G2/M. What is
below is reported as sub-G1 and what is above as above G2 (">4N"), each
as a percentage of all events; %G1, %S and %G2 are of the cells in the
three phases. Unless aggregates are modelled (below), gate single cells
first (area against width or height of the DNA signal), or two G1 cells
stuck together are counted in G2/M.

With Dean-Jett-Fox, three more components can be switched on (`Extras`),
each fitted with the rest to the whole histogram:

  a second population   (aneuploid) its own G1, S and G2, its G1 at the
      DNA index times the first's (searched within 10% of where it is
      found, or of where the user says), G2 at the same multiple of its
      G1 as the first, its own CV. The lower G1 is the diploid one unless
      `g1` says otherwise; a G1 near 2 cannot be told from the diploid
      G2 and is only taken where its DNA index is given.
  debris                exponential (a falling curve, its rate fitted) or
      nuclei cut once (after Bagwell et al., Cytometry 1991), worked out
      from the intact particles fitted. The fit then reaches down to a
      fifth of G1.
  aggregates            two and three cells stuck together, the
      distribution of single cells convolved with itself (as ModFit's
      aggregate model does; the chance of a clump is not made to grow
      with the size of its cells), each fitted for its number.
Debris and aggregates take their shape from the cells found, so the fit
is repeated four times, each with the shapes the one before gave.
% debris and % aggregates are of all events modelled (an aggregate is
one event); % aneuploid is of the cells of the two populations, and its
phases of its own cells.

  On 30,000 events simulated (`simulate_sample`) at a CV of 3 to 5%:
      debris (25%, either shape) within 0.7 of a point, aggregates (15%)
      within 0.5, the diploid phases then within 1.5; a population at a
      DNA index of 1.5 holding 40% of the cells within about 3 points, its
      phases within 4.5 and its index within 0.005. Where the second
      population lies on the first's S (index 1.25, or 0.8) its phases
      move by 9 to 13 points and the first's by up to 10: the histogram
      cannot say which S a cell is in.
  On 223 files of four ovarian cancer lines stained for BrdU and PI
  (Yang & Hill, Dryad doi:10.5061/dryad.sn02v6x89), PI alone against the
  phases BrdU gives in the same file: see docs/comparacao.md. Fitted
  without gating single cells, the aggregates found (24.8% on average)
  match the events a PI-A/PI-H gate drops (25.6%; r = 0.90), and the
  phases match those of the singlets.

How far to trust them, from histograms of 20,000 cells simulated with
known phases (`simulate`; the tests hold them to it), in percentage
points of a phase:

  S spread evenly, or leaning to one side (what a polynomial can follow):
      Dean-Jett-Fox is unbiased (under 0.3) at a CV of 3 to 5% and varies
      by 0.5 to 1.5 from one acquisition to the next. At a CV of 8% the
      peaks run into S: the bias stays about 1 but S and G2 trade 2 to 5
      points between them from one acquisition to the next.
      Watson is within 1 of the truth on average at any of those CVs
      where S is even (1.5 where most cells are in G2), and noisier: 1.5
      to 3 points, since it uses half of each peak.
  S piled up against one of the peaks:
      Watson is off by up to 2 points at a CV of 5% and up to 5 at 8%,
      between S and G2: it takes S to be level under the side of a peak.
      Dean-Jett-Fox stays within 0.5 at 5% and within 3 at 8%.
  S with a shape no parabola has (a fifth of the cells held at one point
  of S):
      held in the middle of S, both stay within 1 point at a CV of 3% and
      within 3 at 5%. Held a quarter of the way from G1, Dean-Jett-Fox
      takes 6 points from G1 at a CV of 3% and 12 at 5%; Watson 1 and 5.
      The chi-square of Dean-Jett-Fox is then 9 to 30 instead of 1, which
      is how it is known; Watson's stays at 1 whatever S is.

Each number comes with a standard error; see `modelfit` for what it does
and does not cover. Those of Dean-Jett-Fox match the scatter of repeated
simulations. Those of Watson are approximate — the two peaks are taken
to be independent and S to vary as a count — and that of G2 is about
half of the scatter seen.

No Qt here.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from scipy.special import ndtr

from . import figures, modelfit
from .modelfit import Component, FitContent

WATSON, DJF = "watson", "djf"
MODELS = {WATSON: "Watson pragmatic", DJF: "Dean-Jett-Fox"}
#: the multiples of the mean of G0/G1 the mean of G2/M is held between
RATIO = (1.9, 2.1)
#: how far beyond the peaks the fit reaches, in their SD
REACH = 4.0
#: slices S is cut into before each is broadened
_SLICES = 48
#: the parts of the side of a peak that Watson fits, in SD from its mean:
#: all of the side S does not reach, and the top
_SIDE = 0.35
#: times Watson's peaks are fitted again once S under their sides is known
_ROUNDS = 4
PHASES = ["G1", "S", "G2"]
NO_DEBRIS, EXPONENTIAL, SINGLE_CUT = "", "exponential", "single cut"
DEBRIS = {NO_DEBRIS: "None", EXPONENTIAL: "Exponential", SINGLE_CUT: "Single cut"}
ANEUPLOID = ["Aneuploid G1", "Aneuploid S", "Aneuploid G2"]
#: how far down the fit reaches where debris is modelled, as a multiple
#: of the mean of G0/G1
DEBRIS_FROM = 0.2
#: a second G1 is looked for no nearer than this to the first or its G2
_APART = 0.08


@dataclass(frozen=True)
class Extras:
    """What is fitted besides G0/G1, S and G2/M of one population: a
    second (aneuploid) population with its own G1, S and G2, near the
    DNA index given or wherever its G1 is found; debris of one of the two
    shapes; and aggregates of two and three cells."""
    aneuploid: bool = False
    dna_index: float | None = None
    debris: str = NO_DEBRIS
    aggregates: bool = False

    def __bool__(self) -> bool:
        return self.aneuploid or bool(self.debris) or self.aggregates


@dataclass
class Fit:
    x: np.ndarray
    counts: np.ndarray
    model: np.ndarray
    #: the curves of G1, S and G2 over the bins
    components: np.ndarray
    kind: str
    #: events in G1, S and G2 and the standard error of each
    cells: np.ndarray
    cells_se: np.ndarray
    #: percent of the three phases together, and standard errors
    percent: np.ndarray
    percent_se: np.ndarray
    g1_mean: float
    g1_cv: float
    g2_mean: float
    g2_cv: float
    ratio: float
    ratio_se: float
    #: the DNA content the fit reaches from and to
    fitted: tuple
    #: percent of all events below and above that
    sub_g1: float
    above_g2: float
    events: int
    rms: float
    rms_percent: float
    chi2: float
    #: what else was modelled, and the name of each curve of `components`
    extras: Extras = field(default_factory=lambda: Extras())
    names: list = field(default_factory=lambda: list(PHASES))
    #: the share each curve is labelled with: a phase as a percentage of
    #: its population, debris and aggregates as one of all events modelled
    shares: np.ndarray | None = None
    #: the second population: its phases (percent of it), its DNA index
    #: (mean of its G1 over that of the first), the CV of its G1, and its
    #: share of the cells of the two populations
    aneuploid_percent: np.ndarray = field(default_factory=lambda: np.full(3, np.nan))
    aneuploid_percent_se: np.ndarray = field(default_factory=lambda: np.full(3, np.nan))
    dna_index: float = float("nan")
    dna_index_se: float = float("nan")
    aneuploid_cv: float = float("nan")
    aneuploid: float = float("nan")
    aneuploid_se: float = float("nan")
    #: percent of all events modelled that are debris, and aggregates
    debris: float = float("nan")
    debris_se: float = float("nan")
    aggregates: float = float("nan")
    aggregates_se: float = float("nan")


def _peaks(values: np.ndarray, bins: int) -> tuple:
    """(centres, the histogram smoothed, the bins of its peaks worth the
    name), up to the 99.5th percentile."""
    top = float(np.percentile(values, 99.5))
    x, counts, _ = modelfit.histogram(values, 0.0, top, bins)
    kernel = np.exp(-0.5 * (np.arange(-6, 7) / 2.0) ** 2)
    smooth = np.convolve(counts, kernel / kernel.sum(), mode="same")
    peaks, _ = find_peaks(smooth, prominence=0.05 * smooth.max())
    if peaks.size == 0:
        peaks = np.array([int(smooth.argmax())])
    return x, smooth, peaks


def _find_g1(values: np.ndarray, bins: int) -> tuple[float, float]:
    """Roughly where G0/G1 is and how wide: the tallest peak, or the one
    at half its DNA content if there is one worth the name (cells held in
    G2/M may outnumber those in G1)."""
    x, smooth, peaks = _peaks(values, bins)
    tallest = peaks[smooth[peaks].argmax()]
    half = [p for p in peaks if abs(x[p] / x[tallest] - 0.5) < 0.07
            and smooth[p] >= 0.1 * smooth[tallest]]
    g1 = max(half, key=lambda p: smooth[p]) if half else tallest
    below = np.flatnonzero(smooth[:g1] < smooth[g1] / 2)       # half its height, on the left
    hwhm = x[g1] - x[below[-1]] if below.size else 0.05 * x[g1]
    return float(x[g1]), float(np.clip(hwhm / 1.177 / x[g1], 0.01, 0.15))


def _s_basis(x, width, m1, m2, cv1, cv2, degree: int = 2) -> np.ndarray:
    """A column for each term of the polynomial of S, each of area 1:
    Bernstein's, which are never negative, so that heights that are not
    negative give a polynomial that is not."""
    t = (np.arange(_SLICES) + 0.5) / _SLICES
    at = m1 + t * (m2 - m1)
    sd = (cv1 + t * (cv2 - cv1)) * at
    spread = modelfit.gaussian(x[:, None], at[None, :], sd[None, :], width)
    terms = {0: [np.ones_like(t)],
             2: [3 * (1 - t) ** 2, 6 * t * (1 - t), 3 * t ** 2]}[degree]
    return np.column_stack([spread @ (term / _SLICES) for term in terms])


def _djf_basis(x, width, p, same_cv: bool, degree: int = 2) -> np.ndarray:
    m1, cv1, ratio = p[:3]
    cv2 = cv1 if same_cv else p[3]
    m2 = ratio * m1
    return np.column_stack([modelfit.gaussian(x, m1, cv1 * m1, width),
                            _s_basis(x, width, m1, m2, cv1, cv2, degree),
                            modelfit.gaussian(x, m2, cv2 * m2, width)])


def _djf(x, counts, width, start, bounds, same_cv: bool, degree: int = 2):
    def residual(p):
        b = _djf_basis(x, width, p, same_cv, degree)
        return counts - b @ modelfit.heights(b, counts)
    found = least_squares(residual, start, bounds=bounds, x_scale=[start[0] * 0.02, 0.01,
                                                                   0.02, 0.01][:len(start)],
                          diff_step=1e-3)
    return found.x


def _window(m1, cv1, m2, cv2) -> tuple[float, float]:
    return m1 * (1 - REACH * cv1), m2 * (1 + REACH * cv2)


def _ceiling(v: np.ndarray) -> float:
    """The top of the scale where events pile up on it (the detector
    saturated): they are left out of every bin. Infinite where they do
    not."""
    top = float(v.max())
    return top if np.count_nonzero(v >= top) >= max(0.002 * v.size, 3) else np.inf


def fit(values, kind: str = WATSON, ratio: tuple = RATIO, same_cv: bool = False,
        g1: float | None = None, bins: int = modelfit.BINS,
        extras: Extras | None = None) -> Fit:
    """The phases of one sample. `values` are the DNA signal of each event
    as the file stores it, linear; `g1` is about where G0/G1 is, when the
    tallest peak (or the one at half its DNA) is not it. `extras`: a
    second population, debris, aggregates (Dean-Jett-Fox only)."""
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size < 100:
        raise ValueError("too few events")
    if kind not in MODELS:
        raise ValueError(f"unknown model {kind!r}")
    extras = extras or Extras()
    if extras.debris not in DEBRIS:
        raise ValueError(f"unknown debris {extras.debris!r}")
    if extras and kind != DJF:
        raise ValueError("a second population, debris and aggregates are fitted with "
                         "Dean-Jett-Fox only: Watson's S has no shape to add them to")
    lo_ratio, hi_ratio = sorted(float(r) for r in ratio)
    hi_ratio = max(hi_ratio, lo_ratio + 1e-6)
    m1, cv1 = _find_g1(v[v > 0], bins)
    if g1:
        near = v[(v > 0.8 * g1) & (v < 1.2 * g1)]
        m1 = float(np.median(near)) if near.size > 20 else float(g1)
    if extras:
        return _extended(v, m1, cv1, (lo_ratio, hi_ratio), same_cv, bool(g1), bins, extras)
    # the histogram: G0/G1 at about a third of the way, G2/M and its tail on
    # it; not beyond where a saturated detector piles events up
    top = m1 * hi_ratio * 1.35
    ceiling = _ceiling(v)
    if ceiling < top:
        bins = max(round(bins * ceiling / top), 64)
        top = ceiling * (1 - 1e-9)
    x, counts, width = modelfit.histogram(v, 0.0, top, bins)
    start = [m1, cv1, min(max(2.0, lo_ratio), hi_ratio)] + ([] if same_cv else [cv1])
    lower = [m1 * 0.85, 0.005, lo_ratio] + ([] if same_cv else [0.005])
    upper = [m1 * 1.15, 0.2, hi_ratio] + ([] if same_cv else [0.2])
    inside = np.ones(x.size, dtype=bool)
    p = np.array(start)
    for _ in range(2):   # the second time, on the part of the histogram the first says to fit
        p = _djf(x[inside], counts[inside], width, p, (lower, upper), same_cv)
        cv2 = p[1] if same_cv else p[3]
        lo, hi = _window(p[0], p[1], p[0] * p[2], cv2)
        inside = (x >= lo) & (x <= hi)
    if kind == DJF:
        found = _finish_djf(x, counts, width, inside, p, same_cv)
    else:
        found = _watson(x, counts, width, p, same_cv, (lo_ratio, hi_ratio))
    return _report(v, x, counts, kind, *found)


def _finish_djf(x, counts, width, inside, p, same_cv):
    xi, ci = x[inside], counts[inside]
    b = _djf_basis(xi, width, p, same_cv)
    h = modelfit.heights(b, ci)
    active = np.flatnonzero(h > 0)
    theta = np.concatenate([p, h[active]])

    def curve(t):
        return _djf_basis(xi, width, t[:p.size], same_cv)[:, active] @ t[p.size:]

    def numbers(t):
        full = np.zeros(h.size)
        full[active] = t[p.size:]
        cells = np.array([full[0], full[1:-1].sum(), full[-1]])
        return np.concatenate([cells, 100.0 * cells / max(cells.sum(), 1e-12), [t[2]]])

    se = modelfit.propagate(numbers, theta, modelfit.covariance(curve, theta))
    full = _djf_basis(x, width, p, same_cv) * h
    parts = np.array([full[:, 0], full[:, 1:-1].sum(axis=1), full[:, -1]])
    cells = np.array([h[0], h[1:-1].sum(), h[-1]])
    m1, cv1, ratio = p[:3]
    cv2 = cv1 if same_cv else p[3]
    return (parts, parts.sum(axis=0), cells, se[:3], se[3:6], (m1, cv1, m1 * ratio, cv2),
            float(se[6]), inside, theta.size)


def _side(x, counts, width, mean, sd, left: bool, bounds=None, cv: float | None = None,
          again: int = 3):
    """A Gaussian fitted to one side of a peak and its top. Returns
    (events, mean, sd, (its curve, what was fitted), the bins it was
    fitted to)."""
    found = None
    for _ in range(again):
        keep = ((x <= mean + _SIDE * sd) & (x >= mean - REACH * sd) if left
                else (x >= mean - _SIDE * sd) & (x <= mean + REACH * sd))
        if keep.sum() < 4:
            break
        xi, ci = x[keep], counts[keep]
        area = max(float(ci.sum()) * 2 / (1 + _SIDE), 1.0)

        def curve(t, xi=xi):
            s = cv * t[1] if cv is not None else t[2]
            return t[0] * modelfit.gaussian(xi, t[1], s, width)

        lo_mean, hi_mean = bounds or (mean - 2 * sd, mean + 2 * sd)
        start = [area, float(np.clip(mean, lo_mean, hi_mean))] + ([] if cv is not None else [sd])
        lower = [0.0, lo_mean] + ([] if cv is not None else [sd / 4])
        upper = [np.inf, hi_mean] + ([] if cv is not None else [sd * 4])
        found = least_squares(lambda t: ci - curve(t), start, bounds=(lower, upper),
                              x_scale=[area, sd, sd][:len(start)])
        mean = float(found.x[1])
        sd = cv * mean if cv is not None else float(found.x[2])
    if found is None:
        raise ValueError("no peak to fit: too few bins on its side")
    return float(found.x[0]), mean, sd, (curve, found.x), keep


def _watson(x, counts, width, p, same_cv, ratio):
    m1, cv1, r = p[:3]
    cv2 = cv1 if same_cv else p[3]
    sd1, m2, sd2 = cv1 * m1, r * m1, cv2 * r * m1
    under = np.zeros(x.size)
    for k in range(_ROUNDS):
        rest = counts - under
        again = 3 if k == 0 else 1
        n1, m1, sd1, curve1, kept1 = _side(x, rest, width, m1, sd1, left=True, again=again)
        n2, m2, sd2, curve2, kept2 = _side(x, rest, width, m2, sd2, left=False,
                                           bounds=(ratio[0] * m1, ratio[1] * m1),
                                           cv=sd1 / m1 if same_cv else None, again=again)
        g1 = n1 * modelfit.gaussian(x, m1, sd1, width)
        g2 = n2 * modelfit.gaussian(x, m2, sd2, width)
        left = counts - g1 - g2
        # S just above G1 and just below G2, where the peaks have ended:
        # that much of it, blurred as the peaks are, lies under their sides
        middle = (m1 + m2) / 2
        first = (x >= m1 + 2 * sd1) & (x <= min(m1 + REACH * sd1, middle))
        last = (x <= m2 - 2 * sd2) & (x >= max(m2 - REACH * sd2, middle))
        s1 = max(float(left[first].mean()), 0.0) if first.any() else 0.0
        s2 = max(float(left[last].mean()), 0.0) if last.any() else 0.0
        under = np.where(x <= middle, s1 * ndtr((x - m1) / sd1), s2 * ndtr((m2 - x) / sd2))
    early, late = under * (x <= m1), under * (x >= m2)
    cov1, cov2 = modelfit.covariance(*curve1), modelfit.covariance(*curve2)
    cv1 = sd1 / m1
    between = (x > m1) & (x < m2)
    s = np.where(between, np.clip(left, 0.0, None), 0.0) + early + late
    # counted without the clipping: noise below the curve of a peak takes
    # from S what noise above it gives
    in_s = max(float(left[between].sum() + early.sum() + late.sum()), 0.0)
    cells = np.array([n1, in_s, n2])
    # S is what is left of the events between the means: it has their
    # counting error, and that of the two halves taken from it
    var = np.array([cov1[0, 0], 0.0, cov2[0, 0]])
    var[1] = float(counts[between].sum()) + (var[0] + var[2]) / 4
    cells_se = np.sqrt(np.maximum(var, 0.0))
    total = cells.sum()
    grad = (np.eye(3) * total - cells[:, None]) * 100.0 / max(total, 1e-12) ** 2
    cov = np.diag(var)
    cov[0, 1] = cov[1, 0] = -var[0] / 2
    cov[2, 1] = cov[1, 2] = -var[2] / 2
    percent_se = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", grad, cov, grad), 0.0))
    ratio_se = (m2 / m1) * np.sqrt(cov1[1, 1] / m1 ** 2 + cov2[1, 1] / m2 ** 2)
    lo, hi = _window(m1, cv1, m2, sd2 / m2)
    inside = (x >= lo) & (x <= hi)
    parts = np.array([g1, s, g2])
    cv2 = sd2 / m2
    return (parts, parts.sum(axis=0), cells, cells_se, percent_se, (m1, cv1, m2, cv2),
            float(ratio_se), inside & (kept1 | kept2), 5 if same_cv else 6)


def _report(v, x, counts, kind, parts, model, cells, cells_se, percent_se, peaks, ratio_se,
            judged, parameters) -> Fit:
    m1, cv1, m2, cv2 = (float(p) for p in peaks)
    lo, hi = _window(m1, cv1, m2, cv2)
    total = cells.sum()
    return Fit(x, counts, model, parts, kind, cells, cells_se,
               100.0 * cells / total if total > 0 else np.full(3, np.nan), percent_se,
               m1, 100.0 * cv1, m2, 100.0 * cv2, m2 / m1, ratio_se, (float(lo), float(hi)),
               100.0 * float((v < lo).sum()) / v.size, 100.0 * float((v > hi).sum()) / v.size,
               int(v.size), modelfit.rms(counts[judged], model[judged]),
               100.0 * modelfit.rms(counts[judged], model[judged]) / max(float(counts.max()), 1.0),
               modelfit.chi2(counts[judged], model[judged], parameters))


def simulate(n: int, g1: float, s: float, g2: float, mean: float = 50000.0, cv: float = 0.04,
             rng=None, ratio: float = 2.0, early: float = 0.0, sub_g1: float = 0.0) -> np.ndarray:
    """The DNA signal of `n` cells of which the shares `g1`, `s`, `g2`
    are in each phase, measured with the coefficient of variation `cv`.
    The cells in S have a DNA content spread evenly between the two
    peaks — or, with `early` (−1 to 1), more of them towards G1 (above 0)
    or G2 (below). `sub_g1` adds that share of all events as debris."""
    rng = rng or np.random.default_rng(0)
    shares = np.array([g1, s, g2], dtype=np.float64)
    phase = rng.choice(3, size=n, p=shares / shares.sum())
    u = rng.random(n)
    if early:   # a density 1 + early·(1 − 2t) over the way through S
        u = ((1 + early) - np.sqrt((1 + early) ** 2 - 4 * early * u)) / (2 * early)
    dna = np.where(phase == 0, 1.0, np.where(phase == 2, ratio, 1.0 + u * (ratio - 1.0)))
    out = mean * dna * (1.0 + cv * rng.standard_normal(n))
    if sub_g1 > 0:
        debris = rng.random(n) < sub_g1
        out = np.where(debris, mean * rng.uniform(0.1, 0.7, n), out)
    return out


# --- a second population, debris and aggregates ---------------------------------------

#: times the shapes of debris and aggregates are worked out again from the
#: cells the fit before found
_EXTENDED_ROUNDS = 4


def _second_g1(v: np.ndarray, m1: float, ratio: tuple, bins: int) -> float | None:
    """Where the G0/G1 of a second population is: the tallest peak that is
    neither the G0/G1 at `m1`, nor its G2/M, nor debris. None: no such peak."""
    x, smooth, peaks = _peaks(v[v > 0], bins)
    tallest = float(smooth[peaks].max())
    keep = [p for p in peaks if x[p] > 0.3 * m1 and smooth[p] >= 0.05 * tallest
            and abs(x[p] / m1 - 1) > _APART
            and not ratio[0] * (1 - _APART) < x[p] / m1 < ratio[1] * (1 + _APART)]
    return float(x[max(keep, key=lambda p: smooth[p])]) if keep else None


class _Layout:
    """Which number of the vector searched for is which, and which columns
    of the basis are each component."""

    def __init__(self, same_cv: bool, extras: Extras):
        self.same_cv, self.extras = same_cv, extras
        names = ["m1", "cv1", "ratio"] + ([] if same_cv else ["cv2"])
        names += ["index", "cva"] if extras.aneuploid else []
        names += ["decay"] if extras.debris == EXPONENTIAL else []
        self.at = {n: i for i, n in enumerate(names)}
        self.size = len(names)
        self.groups: dict[str, list[int]] = {"G1": [0], "S": [1, 2, 3], "G2": [4]}
        k = 5
        if extras.aneuploid:
            self.groups.update({"Aneuploid G1": [5], "Aneuploid S": [6, 7, 8],
                                "Aneuploid G2": [9]})
            k = 10
        if extras.debris:
            self.groups["Debris"] = [k]
            k += 1
        if extras.aggregates:
            self.groups["Aggregates"] = [k, k + 1]
            k += 2
        self.columns = k

    def peaks(self, p) -> dict:
        """The means and CVs of the peaks."""
        m1, cv1, r = (float(t) for t in p[:3])
        cv2 = cv1 if self.same_cv else float(p[self.at["cv2"]])
        out = {"m1": m1, "cv1": cv1, "m2": r * m1, "cv2": cv2}
        if self.extras.aneuploid:
            ma, cva = float(p[self.at["index"]]) * m1, float(p[self.at["cva"]])
            out.update(ma1=ma, cva1=cva, ma2=r * ma, cva2=cva * cv2 / cv1)
        return out

    def window(self, p, top: float) -> tuple[float, float]:
        """The DNA content the fit reaches from and to."""
        k = self.peaks(p)
        means = [(k["m1"], k["cv1"], k["m2"], k["cv2"])]
        if self.extras.aneuploid:
            means.append((k["ma1"], k["cva1"], k["ma2"], k["cva2"]))
        lo = min(a * (1 - REACH * c) for a, c, _, _ in means)
        hi = max(b * (1 + REACH * c) for _, _, b, c in means)
        if self.extras.debris:
            lo = DEBRIS_FROM * k["m1"]
        if self.extras.aggregates:   # two cells in G2/M stuck together
            hi = max(hi, max(2 * b * (1 + REACH * c / np.sqrt(2)) for _, _, b, c in means))
        return lo, min(hi, top)

    def basis(self, x, width, p, shapes: dict) -> np.ndarray:
        k = self.peaks(p)
        cols = [modelfit.gaussian(x, k["m1"], k["cv1"] * k["m1"], width),
                _s_basis(x, width, k["m1"], k["m2"], k["cv1"], k["cv2"]),
                modelfit.gaussian(x, k["m2"], k["cv2"] * k["m2"], width)]
        if self.extras.aneuploid:
            cols += [modelfit.gaussian(x, k["ma1"], k["cva1"] * k["ma1"], width),
                     _s_basis(x, width, k["ma1"], k["ma2"], k["cva1"], k["cva2"]),
                     modelfit.gaussian(x, k["ma2"], k["cva2"] * k["ma2"], width)]
        if self.extras.debris == EXPONENTIAL:
            decay = float(p[self.at["decay"]]) * k["m1"]
            cols.append(width * np.exp(-x / decay) / decay)
        elif self.extras.debris == SINGLE_CUT:
            cols.append(shapes["cut"])
        if self.extras.aggregates:
            cols += [shapes["doublets"], shapes["triplets"]]
        return np.column_stack(cols)


def _aggregates(x: np.ndarray, width: float, singlets: np.ndarray) -> tuple:
    """(doublets, triplets), each of area 1 over all DNA contents (what
    lies beyond the histogram is in no bin): the DNA of two or three cells
    drawn at random from the single cells, each as measured — the
    distribution of single cells convolved with itself."""
    s = np.clip(singlets, 0.0, None)
    s = s / max(float(s.sum()), 1e-12)
    two = np.convolve(s, s)        # bins i and j: a content of (i + j + 1) widths
    three = np.convolve(two, s)    # (i + j + k + 1.5) widths
    k = np.arange(two.size)
    doublets = np.interp(x, (k + 1.0) * width, two, left=0.0, right=0.0)
    k = np.arange(three.size)
    triplets = np.interp(x, (k + 1.5) * width, three, left=0.0, right=0.0)
    return doublets, triplets


def _single_cut(x: np.ndarray, intact: np.ndarray) -> np.ndarray:
    """Debris of nuclei cut once at a random point (after Bagwell et al.,
    Cytometry 1991), of area 1: a nucleus of DNA y gives two pieces each
    anywhere from 0 to y, and is cut with a chance taken to grow with its
    surface, y^(2/3). The density at x is the sum over the nuclei above
    it of their number times y^(2/3) / y."""
    y = np.clip(intact, 0.0, None) * np.maximum(x, 1e-12) ** (-1.0 / 3.0)
    above = np.cumsum(y[::-1])[::-1] - y / 2       # all bins above, and half of its own
    return above / max(float(above.sum()), 1e-12)


def _extended(v, m1, cv1, ratio, same_cv: bool, told: bool, bins: int, extras: Extras) -> Fit:
    """Dean-Jett-Fox with a second population, debris or aggregates."""
    lo_ratio, hi_ratio = ratio
    index = None
    if extras.aneuploid:
        if extras.dna_index:
            index = float(extras.dna_index)
        else:
            ma = _second_g1(v, m1, ratio, bins)
            if ma is None:
                raise ValueError("no second G1 peak found: give the DNA index of the "
                                 "aneuploid population")
            if ma < m1 and not told:     # the lower of the two is taken as diploid
                m1, ma = ma, m1
            index = ma / m1
    # the histogram: as far as the last G2/M, or two of them stuck together,
    # with bins as wide as without them
    furthest = m1 * hi_ratio * max(1.0, index or 1.0)
    top = furthest * (2.3 if extras.aggregates else 1.35)
    plain = m1 * hi_ratio * 1.35
    ceiling = _ceiling(v)
    if ceiling < top:
        top = ceiling * (1 - 1e-9)
    bins = int(np.clip(round(bins * top / plain), 64, 1024))
    x, counts, width = modelfit.histogram(v, 0.0, top, bins)

    layout = _Layout(same_cv, extras)
    start = [m1, cv1, min(max(2.0, lo_ratio), hi_ratio)] + ([] if same_cv else [cv1])
    lower = [m1 * 0.85, 0.005, lo_ratio] + ([] if same_cv else [0.005])
    upper = [m1 * 1.15, 0.2, hi_ratio] + ([] if same_cv else [0.2])
    scale = [m1 * 0.02, 0.01, 0.02] + ([] if same_cv else [0.01])
    if extras.aneuploid:
        start += [index, cv1]
        lower += [index * 0.9, 0.005]
        upper += [index * 1.1, 0.2]
        scale += [0.02, 0.01]
    if extras.debris == EXPONENTIAL:
        start += [0.3]
        lower += [0.02]
        upper += [3.0]
        scale += [0.05]
    p = np.array(start, dtype=np.float64)

    # to start with, the events above the debris are taken for cells
    first = min(m1, m1 * (index or 1.0))
    cells = np.where(x > 0.75 * first, counts, 0.0)
    shapes: dict = {}

    def reshape(singlets, intact):
        if extras.aggregates:
            shapes["doublets"], shapes["triplets"] = _aggregates(x, width, singlets)
        if extras.debris == SINGLE_CUT:
            shapes["cut"] = _single_cut(x, intact)

    reshape(cells, cells)
    lo, hi = layout.window(p, top)
    inside = (x >= lo) & (x <= hi)
    rounds = _EXTENDED_ROUNDS if (extras.aggregates or extras.debris == SINGLE_CUT) else 2
    for k in range(rounds):
        ci = counts[inside]

        def residual(t, inside=inside, ci=ci):
            b = layout.basis(x, width, t, shapes)[inside]
            return ci - b @ modelfit.heights(b, ci)

        p = least_squares(residual, p, bounds=(lower, upper), x_scale=scale,
                          diff_step=1e-3).x
        if k == rounds - 1:
            break
        b = layout.basis(x, width, p, shapes)
        parts = b * modelfit.heights(b[inside], ci)
        g = layout.groups
        single = [c for name in g if name in PHASES or name in ANEUPLOID for c in g[name]]
        singlets = parts[:, single].sum(axis=1)
        intact = singlets + (parts[:, g["Aggregates"]].sum(axis=1) if extras.aggregates else 0)
        reshape(singlets, intact)
        lo, hi = layout.window(p, top)
        inside = (x >= lo) & (x <= hi)
    return _finish_extended(v, x, counts, width, inside, p, layout, shapes, (lo, hi))


def _finish_extended(v, x, counts, width, inside, p, layout: _Layout, shapes, window) -> Fit:
    ci = counts[inside]
    b = layout.basis(x, width, p, shapes)
    h = modelfit.heights(b[inside], ci)
    active = np.flatnonzero(h > 0)
    theta = np.concatenate([p, h[active]])
    g = layout.groups
    extras = layout.extras

    def curve(t):
        return layout.basis(x, width, t[:p.size], shapes)[inside][:, active] @ t[p.size:]

    def numbers(t):
        full = np.zeros(h.size)
        full[active] = t[p.size:]
        total = max(full.sum(), 1e-12)
        dip = np.array([full[g[n]].sum() for n in PHASES])
        out = [dip, 100.0 * dip / max(dip.sum(), 1e-12), [t[2]]]
        if extras.aneuploid:
            ane = np.array([full[g[n]].sum() for n in ANEUPLOID])
            out += [100.0 * ane / max(ane.sum(), 1e-12), [t[layout.at["index"]]],
                    [100.0 * ane.sum() / max(ane.sum() + dip.sum(), 1e-12)]]
        if extras.debris:
            out.append([100.0 * full[g["Debris"]].sum() / total])
        if extras.aggregates:
            out.append([100.0 * full[g["Aggregates"]].sum() / total])
        return np.concatenate(out)

    value = numbers(theta)
    se = modelfit.propagate(numbers, theta, modelfit.covariance(curve, theta))
    full = b * h
    names = list(g)
    parts = np.array([full[:, g[n]].sum(axis=1) for n in names])
    k = layout.peaks(p)
    f = _report(v, x, counts, DJF, parts, parts.sum(axis=0), value[:3], se[:3], se[3:6],
                (k["m1"], k["cv1"], k["m2"], k["cv2"]), float(se[6]), inside, theta.size)
    lo, hi = window
    f.fitted = (float(lo), float(hi))
    f.sub_g1 = 100.0 * float((v < lo).sum()) / v.size
    f.above_g2 = 100.0 * float((v > hi).sum()) / v.size
    f.extras, f.names = extras, names
    shares = dict(zip(PHASES, f.percent))
    at = 7
    if extras.aneuploid:
        f.aneuploid_percent, f.aneuploid_percent_se = value[at:at + 3], se[at:at + 3]
        f.dna_index, f.dna_index_se = float(value[at + 3]), float(se[at + 3])
        f.aneuploid, f.aneuploid_se = float(value[at + 4]), float(se[at + 4])
        f.aneuploid_cv = 100.0 * k["cva1"]
        shares.update(zip(ANEUPLOID, f.aneuploid_percent))
        at += 5
    if extras.debris:
        f.debris, f.debris_se = float(value[at]), float(se[at])
        shares["Debris"] = f.debris
        at += 1
    if extras.aggregates:
        f.aggregates, f.aggregates_se = float(value[at]), float(se[at])
        shares["Aggregates"] = f.aggregates
    f.shares = np.array([shares[n] for n in names])
    return f


def simulate_sample(n: int, phases=(60, 25, 15), mean: float = 50000.0, cv: float = 0.04,
                    rng=None, ratio: float = 2.0, aneuploid: tuple | None = None,
                    debris: float = 0.0, debris_kind: str = SINGLE_CUT, decay: float = 0.3,
                    aggregates: float = 0.0, triplets: float = 0.2) -> np.ndarray:
    """The DNA signal of `n` events of which the share `debris` is debris
    and `aggregates` two or three cells stuck together (`triplets` of them
    three), the rest single cells. `aneuploid` = (share of the single
    cells, DNA index, its phases, its CV) adds a second population.
    Debris is of nuclei cut once (`SINGLE_CUT`, as `_single_cut` models
    it) or exponential, falling by e every `decay` × `mean`. An aggregate
    is the sum of the signals of cells drawn at random, each as measured."""
    rng = rng or np.random.default_rng(0)

    def cells(k: int) -> np.ndarray:
        if not aneuploid:
            return simulate(k, *phases, mean=mean, cv=cv, rng=rng, ratio=ratio)
        share, index, its, its_cv = aneuploid
        second = rng.random(k) < share
        out = simulate(k, *phases, mean=mean, cv=cv, rng=rng, ratio=ratio)
        other = simulate(k, *its, mean=mean * index, cv=its_cv, rng=rng, ratio=ratio)
        return np.where(second, other, out)

    n_debris = round(n * debris)
    n_aggregates = round(n * aggregates)
    n_three = round(n_aggregates * triplets)
    n_two = n_aggregates - n_three
    singles = cells(n - n_debris - n_aggregates)
    clumps = np.concatenate([cells(n_two) + cells(n_two), cells(n_three) + cells(n_three)
                             + cells(n_three)])
    if debris_kind == EXPONENTIAL:
        pieces = rng.exponential(decay * mean, n_debris)
    else:   # a nucleus picked with a chance of y^(2/3), cut anywhere
        intact = np.concatenate([singles, clumps])
        weight = np.clip(intact, 0, None) ** (2 / 3)
        picked = rng.choice(intact, size=n_debris, p=weight / weight.sum())
        pieces = picked * rng.random(n_debris)
    return rng.permutation(np.concatenate([singles, clumps, pieces]))


# --- the table and the figure ---------------------------------------------------------

COLUMNS = ["Events", "% G1", "% G1 SE", "% S", "% S SE", "% G2", "% G2 SE", "G1 mean", "G1 CV",
           "G2 mean", "G2 CV", "G2/G1", "G2/G1 SE", "% sub-G1", "% above G2", "RMS", "RMS %",
           "Chi²"]


ANEUPLOID_COLUMNS = ["% aneuploid", "% aneuploid SE", "DNA index", "DNA index SE",
                     "Aneuploid % G1", "Aneuploid % G1 SE", "Aneuploid % S", "Aneuploid % S SE",
                     "Aneuploid % G2", "Aneuploid % G2 SE", "Aneuploid G1 CV"]
DEBRIS_COLUMNS = ["% debris", "% debris SE"]
AGGREGATE_COLUMNS = ["% aggregates", "% aggregates SE"]


def columns(extras: Extras | None = None) -> list[str]:
    """The columns of the table, for what is modelled."""
    extras = extras or Extras()
    return (list(COLUMNS) + (ANEUPLOID_COLUMNS if extras.aneuploid else [])
            + (DEBRIS_COLUMNS if extras.debris else [])
            + (AGGREGATE_COLUMNS if extras.aggregates else []))


def row(f: Fit) -> list:
    paired = [v for k in range(3) for v in (float(f.percent[k]), float(f.percent_se[k]))]
    out = [f.events, *paired, f.g1_mean, f.g1_cv, f.g2_mean, f.g2_cv, f.ratio, f.ratio_se,
           f.sub_g1, f.above_g2, f.rms, f.rms_percent, f.chi2]
    if f.extras.aneuploid:
        out += [f.aneuploid, f.aneuploid_se, f.dna_index, f.dna_index_se]
        out += [float(v) for k in range(3)
                for v in (f.aneuploid_percent[k], f.aneuploid_percent_se[k])]
        out.append(f.aneuploid_cv)
    if f.extras.debris:
        out += [f.debris, f.debris_se]
    if f.extras.aggregates:
        out += [f.aggregates, f.aggregates_se]
    return out


def summary(f: Fit) -> str:
    g1, s, g2 = f.percent
    out = (f"{MODELS[f.kind]} · G1 {g1:.1f}% · S {s:.1f}% · G2 {g2:.1f}% · "
           f"G2/G1 {f.ratio:.2f} · CV {f.g1_cv:.1f}%")
    if f.extras.aneuploid:
        a1, a_s, a2 = f.aneuploid_percent
        out += (f" · aneuploid {f.aneuploid:.1f}%, DNA index {f.dna_index:.2f}, "
                f"G1 {a1:.1f}% S {a_s:.1f}% G2 {a2:.1f}%")
    if f.extras.debris:
        out += f" · debris {f.debris:.1f}%"
    if f.extras.aggregates:
        out += f" · aggregates {f.aggregates:.1f}%"
    return out + f" · RMS {f.rms:.1f} · χ² {f.chi2:.2f}"


@dataclass
class CellCycleContent(FitContent):
    """The DNA content on a linear axis, with G1, S and G2."""


def content(f: Fit, xlabel: str = "", title: str = "") -> CellCycleContent:
    shares = f.percent if f.shares is None else f.shares
    parts = [Component(name, f.components[k], float(shares[k]))
             for k, name in enumerate(f.names)]
    lines = [title, MODELS[f.kind], f"G2/G1 {f.ratio:.2f}", f"CV of G1 {f.g1_cv:.1f}%"]
    if f.extras.aneuploid:
        lines += [f"Aneuploid {f.aneuploid:.1f}% of cells", f"DNA index {f.dna_index:.2f}",
                  f"CV of aneuploid G1 {f.aneuploid_cv:.1f}%"]
    if f.extras:
        lines.append(f"Reduced χ² {f.chi2:.2f}")
    return CellCycleContent(f.x, f.counts, f.model, parts, xlabel, "Events", None,
                            f.fitted, "\n".join(t for t in lines if t))


figures._DRAW[CellCycleContent] = modelfit._draw_fit
