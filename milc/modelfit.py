"""
What the models of one histogram — the generations of a dye dilution, the
phases of the cell cycle — have in common: the histogram that is fitted,
the fit of the heights that enter the model linearly, the uncertainty of
what is fitted, and the figure of data and components.

How a model is fitted. A model is a sum of components, each a shape of
known area times how many events it holds. For given shapes the events of
each are found by non-negative least squares, which has one answer; only
the few numbers that set the shapes (a position, a width) are searched
for. The fit is by least squares on the counts, unweighted, as FlowJo's
and ModFit's are, and its RMS is in events per bin.

The uncertainty that is reported is that of counting: a bin of n events
varies by about the square root of n, and the standard error of every
number fitted, and of everything computed from them, follows from that
(the sandwich estimator of a least-squares fit whose errors are not all
alike, then the delta method). It says how far the same cells, acquired
again, would move the answer. It does not say whether the model is the
right one: a peak that is not Gaussian biases the answer by more than
this, and the RMS is where that shows.

No Qt here, and matplotlib is imported only when a figure is drawn.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import nnls

from . import figures

BINS = 256


def histogram(x: np.ndarray, lo: float, hi: float, bins: int = BINS) -> tuple:
    """(centres, counts as numbers that can be fitted, width of a bin).
    Events outside [lo, hi] are in no bin."""
    counts, edges = np.histogram(x, bins=bins, range=(lo, hi))
    return (edges[:-1] + edges[1:]) / 2, counts.astype(np.float64), float(edges[1] - edges[0])


def gaussian(x: np.ndarray, mean: float, sd: float, width: float = 1.0) -> np.ndarray:
    """A Gaussian of area 1 in events: what one event of the component
    adds to each bin of this `width`."""
    return width * np.exp(-0.5 * ((x - mean) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


def heights(basis: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """The events of each component (a column of `basis`), none negative."""
    return nnls(basis, counts)[0]


def rms(counts: np.ndarray, model: np.ndarray) -> float:
    return float(np.sqrt(np.mean((counts - model) ** 2))) if counts.size else float("nan")


def chi2(counts: np.ndarray, model: np.ndarray, parameters: int) -> float:
    """Reduced chi-square with the variance of a count taken as the model's
    (at least 1): near 1 where the model leaves nothing but counting."""
    freedom = counts.size - parameters
    if freedom <= 0:
        return float("nan")
    return float(np.sum((counts - model) ** 2 / np.maximum(model, 1.0)) / freedom)


def jacobian(f, theta: np.ndarray) -> np.ndarray:
    """d f / d theta by central differences, a column each."""
    theta = np.asarray(theta, dtype=np.float64)
    columns = []
    for k in range(theta.size):
        h = 1e-4 * max(abs(theta[k]), 1e-3)
        up, down = theta.copy(), theta.copy()
        up[k] += h
        down[k] -= h
        columns.append((np.asarray(f(up), dtype=np.float64)
                        - np.asarray(f(down), dtype=np.float64)) / (2 * h))
    return np.column_stack(columns) if columns else np.zeros((0, 0))


def covariance(curve, theta: np.ndarray) -> np.ndarray:
    """The covariance of the fitted numbers `theta` of `curve(theta)`, a
    least-squares fit to counts: (JᵀJ)⁻¹ Jᵀ V J (JᵀJ)⁻¹, V the variance of
    each count (the model's own value, at least 1)."""
    theta = np.asarray(theta, dtype=np.float64)
    if theta.size == 0:
        return np.zeros((0, 0))
    j = jacobian(curve, theta)
    v = np.maximum(np.asarray(curve(theta), dtype=np.float64), 1.0)
    bread = np.linalg.pinv(j.T @ j)
    return bread @ (j.T * v) @ j @ bread


def propagate(f, theta: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """The standard error of each number of `f(theta)`."""
    theta = np.asarray(theta, dtype=np.float64)
    if theta.size == 0:
        return np.zeros(np.size(f(theta)))
    g = jacobian(f, theta)
    return np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", g, cov, g), 0.0))


# --- the figure ---------------------------------------------------------------------

@dataclass
class Component:
    """One part of the model: what it is called, its curve over the
    histogram's bins and its share of the modelled events, in percent."""
    name: str
    y: np.ndarray
    share: float = float("nan")


@dataclass
class FitContent:
    """One histogram with the model fitted to it."""
    x: np.ndarray
    counts: np.ndarray
    model: np.ndarray
    components: list[Component] = field(default_factory=list)
    xlabel: str = ""
    ylabel: str = "Events"
    #: [(position, label)] where the axis is not a plain number
    xticks: list | None = None
    #: the bins the model was fitted to, as (from, to); None: all of them
    fitted: tuple | None = None
    #: written in the corner of the plot, a line each: the sample and what
    #: the legend does not say
    summary: str = ""
    uses_colormap = False

    def category_names(self) -> list[tuple[str, str]]:
        names = [c.name for c in self.components]
        return list(zip(names, figures.categorical(len(names))))


def _draw_fit(fig, c: FitContent, style: figures.FigureStyle) -> None:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax = fig.subplots()
    width = float(c.x[1] - c.x[0]) if c.x.size > 1 else 1.0
    edges = np.concatenate([c.x - width / 2, c.x[-1:] + width / 2])
    ax.stairs(c.counts, edges, fill=True, color="#c9ccd3", linewidth=0, zorder=1)
    handles = [Patch(facecolor="#c9ccd3", label="Events")]
    for comp, (name, default) in zip(c.components, c.category_names()):
        colour = style.color(name, default)
        ax.fill_between(c.x, 0, comp.y, color=colour, alpha=0.30, linewidth=0, zorder=2)
        ax.plot(c.x, comp.y, color=colour, linewidth=0.8, zorder=3)
        share = f" · {comp.share:.1f}%" if np.isfinite(comp.share) else ""
        handles.append(Line2D([], [], color=colour, linewidth=1.6,
                              label=style.label(name) + share))
    ax.plot(c.x, c.model, color="black", linewidth=0.9, zorder=4)
    handles.append(Line2D([], [], color="black", linewidth=1.2, label="Model"))
    if c.fitted is not None:
        for at in c.fitted:
            ax.axvline(at, color="#888888", linewidth=0.6, linestyle=":", zorder=1)
    if c.xticks:
        ax.set_xticks([p for p, _ in c.xticks], [figures._tick_text(t) for _, t in c.xticks])
    ax.set_xlim(edges[0], edges[-1])
    top = max(float(c.counts.max(initial=1.0)), float(c.model.max(initial=1.0))) * 1.06
    ax.set_xlabel(c.xlabel if style.xlabel is None else style.xlabel)
    ax.set_ylabel(c.ylabel if style.ylabel is None else style.ylabel)
    if c.summary:
        # in the upper corner the data leave emptier, the axis raised until
        # the lines clear the data under them
        third = max(c.x.size * 2 // 5, 1)
        tallest = np.maximum(c.counts, c.model)
        left, right = float(tallest[:third].max()), float(tallest[-third:].max())
        corner = 0.02 if left <= right else 0.98
        lines = c.summary.count("\n") + 1
        room = 1.0 - 0.02 - lines * 0.062 * (style.font_size - 1) / 9
        top = max(top, min(left, right) * 1.04 / max(room, 0.3))
        ax.text(corner, 0.98, c.summary, transform=ax.transAxes, va="top",
                ha="left" if corner < 0.5 else "right", fontsize=style.font_size - 1,
                linespacing=1.3)
    ax.set_ylim(0, top)
    figures._frame(ax, style)
    figures._legend(fig, [ax], handles, style, "")
