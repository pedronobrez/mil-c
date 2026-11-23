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
three phases. Doublets are not modelled: gate single cells first (area
against width or height of the DNA signal), or two G1 cells stuck
together are counted in G2/M.

How far to trust them, from histograms of 20,000 cells simulated with
known phases (`simulate`; the tests hold them to it), in percentage
points of a phase:

  S spread evenly, or leaning to one side (what a polynomial can follow):
      Dean-Jett-Fox is unbiased (under 0.3) at a CV of 3 to 5% and varies
      by 0.5 to 1.5 from one acquisition to the next. At a CV of 8% the
      peaks run into S: the bias stays about 1 but S and G2 trade 3 to 5
      points between them from one acquisition to the next.
      Watson is within 1 of the truth on average at any of those CVs
      where S is even (1.5 where most cells are in G2), and noisier: 1.5
      to 3 points, since it uses half of each peak.
  S piled up against one of the peaks:
      Watson moves 1 to 2 points from S to that peak at a CV of 5%, and
      4 to 5 at 8%: it takes S to be level under the side of the peak.
      Dean-Jett-Fox stays within 0.5 at 5% and within 3 at 8%.
  S with a shape no parabola has (a fifth of the cells held at one point
  of S):
      held in the middle of S, both stay within 1 point at a CV of 3% and
      2.5 at 5%. Held a quarter of the way from G1, Dean-Jett-Fox takes
      6 points from G1 at a CV of 3% and 12 at 5%; Watson 1 and 4. The
      chi-square of Dean-Jett-Fox is then 10 to 30 instead of 1, which
      is how it is known; Watson's stays at 1 whatever S is.

Each number comes with a standard error; see `modelfit` for what it does
and does not cover. Those of Dean-Jett-Fox match the scatter of repeated
simulations. Those of Watson are approximate — the two peaks are taken
to be independent and S to vary as a count — and that of G2 is about
half of the scatter seen.

No Qt here.
"""

from __future__ import annotations

from dataclasses import dataclass

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


def _find_g1(values: np.ndarray, bins: int) -> tuple[float, float]:
    """Roughly where G0/G1 is and how wide: the tallest peak, or the one
    at half its DNA content if there is one worth the name (cells held in
    G2/M may outnumber those in G1)."""
    top = float(np.percentile(values, 99.5))
    x, counts, _ = modelfit.histogram(values, 0.0, top, bins)
    kernel = np.exp(-0.5 * (np.arange(-6, 7) / 2.0) ** 2)
    smooth = np.convolve(counts, kernel / kernel.sum(), mode="same")
    peaks, _ = find_peaks(smooth, prominence=0.05 * smooth.max())
    if peaks.size == 0:
        peaks = np.array([int(smooth.argmax())])
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


def fit(values, kind: str = WATSON, ratio: tuple = RATIO, same_cv: bool = False,
        g1: float | None = None, bins: int = modelfit.BINS) -> Fit:
    """The phases of one sample. `values` are the DNA signal of each event
    as the file stores it, linear; `g1` is about where G0/G1 is, when the
    tallest peak (or the one at half its DNA) is not it."""
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size < 100:
        raise ValueError("too few events")
    if kind not in MODELS:
        raise ValueError(f"unknown model {kind!r}")
    lo_ratio, hi_ratio = sorted(float(r) for r in ratio)
    hi_ratio = max(hi_ratio, lo_ratio + 1e-6)
    m1, cv1 = _find_g1(v[v > 0], bins)
    if g1:
        near = v[(v > 0.8 * g1) & (v < 1.2 * g1)]
        m1 = float(np.median(near)) if near.size > 20 else float(g1)
    # the histogram: G0/G1 at about a third of the way, G2/M and its tail on it
    x, counts, width = modelfit.histogram(v, 0.0, m1 * hi_ratio * 1.35, bins)
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


# --- the table and the figure ---------------------------------------------------------

COLUMNS = ["Events", "% G1", "% G1 SE", "% S", "% S SE", "% G2", "% G2 SE", "G1 mean", "G1 CV",
           "G2 mean", "G2 CV", "G2/G1", "G2/G1 SE", "% sub-G1", "% above G2", "RMS", "RMS %",
           "Chi²"]


def row(f: Fit) -> list:
    paired = [v for k in range(3) for v in (float(f.percent[k]), float(f.percent_se[k]))]
    return [f.events, *paired, f.g1_mean, f.g1_cv, f.g2_mean, f.g2_cv, f.ratio, f.ratio_se,
            f.sub_g1, f.above_g2, f.rms, f.rms_percent, f.chi2]


def summary(f: Fit) -> str:
    g1, s, g2 = f.percent
    return (f"{MODELS[f.kind]} · G1 {g1:.1f}% · S {s:.1f}% · G2 {g2:.1f}% · "
            f"G2/G1 {f.ratio:.2f} · CV {f.g1_cv:.1f}% · RMS {f.rms:.1f}")


@dataclass
class CellCycleContent(FitContent):
    """The DNA content on a linear axis, with G1, S and G2."""


def content(f: Fit, xlabel: str = "", title: str = "") -> CellCycleContent:
    parts = [Component(name, f.components[k], float(f.percent[k]))
             for k, name in enumerate(PHASES)]
    said = "\n".join(t for t in (title, MODELS[f.kind], f"G2/G1 {f.ratio:.2f}",
                                 f"CV of G1 {f.g1_cv:.1f}%") if t)
    return CellCycleContent(f.x, f.counts, f.model, parts, xlabel, "Events", None,
                            f.fitted, said)


figures._DRAW[CellCycleContent] = modelfit._draw_fit
