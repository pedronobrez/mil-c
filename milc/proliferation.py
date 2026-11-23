"""
Proliferation from the dilution of a dye (CFSE, CellTrace Violet, …): a
cell that divides gives each daughter half its dye, so the generations
are peaks a factor of two apart — evenly spaced on a logarithmic axis.

The model. The histogram of log10 of the dye is fitted as a sum of
Gaussians, one for each generation g = 0, 1, 2, …, at

    position(g) = position(0) − g · log10(2)

all of one width (or, where asked, a width that grows by the same
fraction with each generation). Generation 0, the undivided peak, is
looked for within `ANCHOR` decades — a third of a division — of the peak
of a control that was not stimulated, if one is given: cells that were
stimulated and did not divide are seldom exactly as bright as cells that
were left alone, but they are never a whole division away, so the
control says which peak is the undivided one and the sample says where
it is. Without a control the brightest peak that holds cells is taken to
be the undivided one, which is wrong where every cell has divided. The
width is always the sample's own. With a
`background` (the autofluorescence of unlabelled cells, in the units of
the file) the dye above it is what halves:
position(g) = log10(background + (10^position(0) − background) / 2^g).

What is reported, with FlowJo's definitions. N(g) is the number of events
in generation g; the cells that were there at the start and whose
descendants those are number N(g) / 2^g.

  Generations        the highest g holding at least `MIN_SHARE` percent
                     of the modelled events
  % Divided          the share of the cells at the start that divided at
                     least once (the precursor frequency):
                     100 · Σ(g≥1) N(g)/2^g  /  Σ(g≥0) N(g)/2^g
  Division Index     divisions per cell at the start, those that never
                     divided included:
                     Σ(g≥1) g·N(g)/2^g  /  Σ(g≥0) N(g)/2^g
  Proliferation      divisions per cell that divided:
  Index              Σ(g≥1) g·N(g)/2^g  /  Σ(g≥1) N(g)/2^g
  Expansion Index    how many times the culture grew:
                     Σ(g≥0) N(g)  /  Σ(g≥0) N(g)/2^g
  Replication Index  how many times the cells that divided grew:
                     Σ(g≥1) N(g)  /  Σ(g≥1) N(g)/2^g

(R's flowFit calls the Expansion Index its "proliferation index".)
Cells that died are not seen, so every index is of the cells that are
there to be counted. RMS is that of the fit, in events per bin and as a
percentage of the tallest bin. Each number comes with a standard error:
see `modelfit` for what it does and does not cover.

No Qt here.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from . import figures, modelfit
from .modelfit import Component, FitContent

HALVING = float(np.log10(2.0))
#: a generation is counted among those found when it holds this share of
#: the modelled events, in percent
MIN_SHARE = 0.5
MAX_GENERATIONS = 8
#: how far from a control's peak the undivided peak of a sample may be,
#: in decades: a third of the distance between two generations
ANCHOR = HALVING / 3
#: the indices, in the order of the table
INDICES = ["% Divided", "Division Index", "Proliferation Index", "Expansion Index",
           "Replication Index"]
#: widths tried for a peak, in decades, before the search starts
_WIDTHS = (0.04, 0.06, 0.08, 0.11, 0.15)
_SUPER = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def indices(cells) -> dict[str, float]:
    """FlowJo's indices from the events of each generation, 0 first."""
    n = np.asarray(cells, dtype=np.float64)
    g = np.arange(n.size)
    start = n / 2.0 ** g                   # the cells at the start, by where they ended
    def over(a, b):
        return float(a / b) if b > 0 else float("nan")
    return {"% Divided": 100.0 * over(start[1:].sum(), start.sum()),
            "Division Index": over((g * start).sum(), start.sum()),
            "Proliferation Index": over((g * start)[1:].sum(), start[1:].sum()),
            "Expansion Index": over(n.sum(), start.sum()),
            "Replication Index": over(n[1:].sum(), start[1:].sum())}


@dataclass
class Peak:
    """The undivided peak of a control: where it is and how wide, in
    log10 of the dye."""
    position: float
    sd: float
    events: int = 0
    rms: float = float("nan")


@dataclass
class Fit:
    #: the histogram, on log10 of the dye, and the model on the same bins
    x: np.ndarray
    counts: np.ndarray
    model: np.ndarray
    #: one curve for each generation, 0 first
    components: np.ndarray
    #: events in each generation, and the standard error of each
    cells: np.ndarray
    cells_se: np.ndarray
    #: log10 of the dye at the undivided peak, the width of a peak in
    #: decades and by how much of it the width grows each generation
    position: float
    sd: float
    growth: float
    position_se: float
    sd_se: float
    #: a control said which peak is the undivided one
    anchored: bool
    generations: int
    #: name → (value, standard error), for `INDICES`
    indices: dict = field(default_factory=dict)
    #: events with a dye above zero, which are all that a logarithm can hold
    events: int = 0
    rms: float = float("nan")
    rms_percent: float = float("nan")
    chi2: float = float("nan")

    @property
    def shares(self) -> np.ndarray:
        """Percent of the modelled events in each generation."""
        total = self.cells.sum()
        return 100.0 * self.cells / total if total > 0 else np.full(self.cells.size, np.nan)

    @property
    def modelled(self) -> float:
        """Percent of the events the generations account for."""
        return 100.0 * float(self.cells.sum()) / self.events if self.events else float("nan")


def _log(values) -> np.ndarray:
    v = np.asarray(values, dtype=np.float64)
    return np.log10(v[np.isfinite(v) & (v > 0)])


def _positions(position: float, generations: int, background: float) -> np.ndarray:
    g = np.arange(generations + 1)
    if background > 0 and 10.0 ** position > background:
        return np.log10(background + (10.0 ** position - background) / 2.0 ** g)
    return position - g * HALVING


def _basis(x, width, position, sd, growth, generations, background,
           floor: float = -np.inf) -> np.ndarray:
    """A column for each generation. One whose peak would be below `floor`,
    under all the events there are, is left empty: only its tail would
    reach them, and a tail fits a handful of events with any area."""
    at = _positions(position, generations, background)
    return np.column_stack([modelfit.gaussian(x, m, sd * (1.0 + growth * g), width)
                            if m >= floor else np.zeros(x.size) for g, m in enumerate(at)])


def undivided(values, bins: int = modelfit.BINS) -> Peak:
    """The peak of a control that did not divide: one Gaussian fitted to
    the histogram of log10 of its dye, around its mode."""
    u = _log(values)
    if u.size < 20:
        raise ValueError("the control has too few events with a dye above zero")
    median = float(np.median(u))
    spread = float(np.diff(np.percentile(u, [15.87, 84.13]))[0] / 2) or 0.05
    x, counts, width = modelfit.histogram(u, median - 6 * spread, median + 6 * spread, bins)

    def residual(p):
        return counts - p[0] * modelfit.gaussian(x, p[1], p[2], width)

    found = least_squares(residual, [float(counts.sum()), median, spread],
                          bounds=([0, x[0], spread / 20], [np.inf, x[-1], spread * 5]))
    return Peak(float(found.x[1]), float(found.x[2]), int(u.size),
                modelfit.rms(counts, counts - found.fun))


def fit(values, control: Peak | None = None, max_generations: int = MAX_GENERATIONS,
        grow: bool = False, background: float = 0.0, bins: int = modelfit.BINS) -> Fit:
    """The generations of one sample. `values` are the dye of each event
    as the file stores it (compensated, not on a display scale);
    `control` the undivided peak of a sample that was not stimulated."""
    u = _log(values)
    if u.size < 50:
        raise ValueError("too few events with a dye above zero")
    low, high = np.percentile(u, [0.5, 99.5])
    if control is not None:
        high = max(high, control.position + 2 * control.sd)
    x, counts, width = modelfit.histogram(u, low - 0.25, high + 0.25, bins)
    generations = max(int(max_generations), 0)
    growth_limit = 0.3 if grow else 0.0

    def shape(p):
        return p[0], p[1], (p[2] if grow else 0.0)

    def basis(p):
        return _basis(x, width, *shape(p), generations, background, low - 0.1)

    def residual(p):
        b = basis(p)
        return counts - b @ modelfit.heights(b, counts)

    def search(start, lower, upper):
        found = least_squares(residual, start, bounds=(lower, upper), x_scale=0.05,
                              diff_step=1e-3)
        return found.x, float(np.sum(found.fun ** 2))

    tail = ([0.0], [0.0], [growth_limit]) if grow else ([], [], [])
    if control is not None:
        at = control.position
        tried = [search([at + shift, sd, *tail[0]], [at - ANCHOR, 0.015, *tail[1]],
                        [at + ANCHOR, 0.3, *tail[2]])
                 for shift in (-ANCHOR / 2, 0.0, ANCHOR / 2)
                 for sd in (control.sd, control.sd * 1.6)]
        best = min(tried, key=lambda t: t[1])[0]
    else:
        best = _unanchored(u, residual, search, tail)
    position, sd, growth = shape(best)
    b = basis(best)
    cells = modelfit.heights(b, counts)
    model = b @ cells

    # the uncertainty: of the shape and of the generations that hold cells
    active = np.flatnonzero(cells > 0)
    free = np.asarray(best, dtype=np.float64)
    theta = np.concatenate([free, cells[active]])

    def curve(t):
        return basis(t[:free.size])[:, active] @ t[free.size:]

    def numbers(t):
        n = np.zeros(cells.size)
        n[active] = t[free.size:]
        return np.concatenate([n, list(indices(n).values())])

    cov = modelfit.covariance(curve, theta)
    se = modelfit.propagate(numbers, theta, cov)
    shape_se = np.sqrt(np.maximum(np.diag(cov)[:free.size], 0.0))
    found = indices(cells)
    shares = 100.0 * cells / cells.sum() if cells.sum() > 0 else np.zeros(cells.size)
    holding = np.flatnonzero(shares >= MIN_SHARE)
    return Fit(x, counts, model, (b * cells).T, cells, se[:cells.size],
               float(position), float(sd), float(growth),
               float(shape_se[0]), float(shape_se[1]),
               control is not None, int(holding.max()) if holding.size else 0,
               {name: (found[name], float(s) if np.isfinite(found[name]) else float("nan"))
                for name, s in zip(found, se[cells.size:])},
               int(u.size), modelfit.rms(counts, model),
               100.0 * modelfit.rms(counts, model) / max(float(counts.max()), 1.0),
               modelfit.chi2(counts, model, theta.size))


def _unanchored(u, residual, search, tail) -> np.ndarray:
    """Where the undivided peak is, when no control says: every place it
    could be is tried, and of those that fit as well as the best the
    dimmest is kept — a peak further up that holds no cells explains
    nothing. So the brightest peak with cells in it is generation 0."""
    low, high = np.percentile(u, [25, 99.9])
    grid = np.arange(high + 0.1, low - 0.01, -0.02)
    errors = np.array([min(float(np.sum(residual([m, s, *tail[0]]) ** 2)) for s in _WIDTHS)
                       for m in grid])
    last = grid.size - 1
    dips = [k for k in range(grid.size)
            if errors[k] <= errors[max(k - 1, 0)] and errors[k] <= errors[min(k + 1, last)]]
    dips = sorted(dips, key=lambda k: errors[k])[:6]
    tried = []
    for k in dips:
        widths = [float(np.sum(residual([grid[k], s, *tail[0]]) ** 2)) for s in _WIDTHS]
        start = [grid[k], _WIDTHS[int(np.argmin(widths))], *tail[0]]
        tried.append(search(start, [grid[k] - HALVING / 2, 0.015, *tail[1]],
                            [grid[k] + HALVING / 2, 0.3, *tail[2]]))
    least = min(t[1] for t in tried)
    good = [t for t in tried if t[1] <= 1.05 * least]
    return min(good, key=lambda t: t[0][0])[0]


def simulate(n: int, shares, position: float = 4.5, sd: float = 0.07, rng=None,
             growth: float = 0.0, background: float = 0.0) -> np.ndarray:
    """The dye of `n` events of which `shares` (generation 0 first, any
    total) are in each generation: what `fit` is checked against."""
    rng = rng or np.random.default_rng(0)
    shares = np.asarray(shares, dtype=np.float64)
    g = rng.choice(shares.size, size=n, p=shares / shares.sum())
    at = _positions(position, shares.size - 1, background)
    return 10.0 ** rng.normal(at[g], sd * (1.0 + growth * g))


# --- the table and the figure ---------------------------------------------------------

def columns(max_generations: int = MAX_GENERATIONS) -> list[str]:
    """The names of what `row` gives for one sample."""
    paired = [c for name in INDICES for c in (name, name + " SE")]
    return (["Events", "Generations", *paired, "Undivided peak", "Peak SD (decades)",
             "RMS", "RMS %", "% modelled"]
            + [f"% in generation {g}" for g in range(max_generations + 1)])


def row(f: Fit, max_generations: int = MAX_GENERATIONS) -> list:
    shares = list(f.shares[:max_generations + 1])
    shares += [float("nan")] * (max_generations + 1 - len(shares))
    paired = [v for name in INDICES for v in f.indices[name]]
    return [f.events, f.generations, *paired, 10.0 ** f.position, f.sd, f.rms,
            f.rms_percent, f.modelled, *shares]


def summary(f: Fit) -> str:
    div, pro = f.indices["Division Index"], f.indices["Proliferation Index"]
    return (f"{f.generations} generation(s) · {f.indices['% Divided'][0]:.1f}% divided · "
            f"Division Index {div[0]:.2f} · Proliferation Index {pro[0]:.2f} · "
            f"RMS {f.rms:.1f}")


def decades(lo: float, hi: float) -> list[tuple[float, str]]:
    """Ticks of an axis in log10: one at each power of ten within it."""
    return [(float(k), "10" + str(k).translate(_SUPER))
            for k in range(int(np.ceil(lo)), int(np.floor(hi)) + 1)]


@dataclass
class ProliferationContent(FitContent):
    """The dye on a logarithmic axis, with each generation's peak."""


def content(f: Fit, xlabel: str = "", title: str = "") -> ProliferationContent:
    used = max(f.generations, int(np.flatnonzero(f.cells > 0).max(initial=0)))
    parts = [Component(f"Generation {g}", f.components[g], float(f.shares[g]))
             for g in range(used + 1)]
    said = "\n".join(t for t in (
        title, f"Division Index {f.indices['Division Index'][0]:.2f}",
        f"Proliferation Index {f.indices['Proliferation Index'][0]:.2f}",
        f"{f.indices['% Divided'][0]:.1f}% divided") if t)
    return ProliferationContent(f.x, f.counts, f.model, parts, xlabel, "Events",
                                decades(f.x[0], f.x[-1]), None, said)


figures._DRAW[ProliferationContent] = modelfit._draw_fit
