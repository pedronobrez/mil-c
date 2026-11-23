"""
A differential-expression analysis as figures: the volcano, every cluster
and marker a point, and the heatmap — clusters down, markers across, each
cell the difference in median between the groups, starred where the
adjusted p says it is one.

Contents like those of `figures`, drawn by its machinery: importing this
module teaches `figures.save` and `figures.preview_png` to draw them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import compare, figures
from .expression import Result
from .volcano import DOWN, REST, UP, VolcanoContent, _draw_volcano

#: the ramp of the differences: white at none
DIVERGING = "RdBu_r"
#: a cell that was not tested
UNTESTED = "#D9D9D9"


def short(marker: str) -> str:
    """`CD8 FITC-A` → `CD8`."""
    return marker.split()[0] if marker else marker


def stars(p: float, alpha: float = 0.05) -> str:
    """What is written in a cell: nothing unless significant."""
    return compare.stars(p) if np.isfinite(p) and p < alpha else ""


@dataclass
class ExpressionVolcano(VolcanoContent):
    """The difference in median across, significance up; a point is a
    cluster and a marker, and is called both."""

    #: one guide line, at no difference: the scale has no "twofold"
    fold: float = 0.0

    def limits(self) -> tuple[float, float]:
        x, y = super().limits()
        return x or 1.0, y

    def points(self, panel: str | None = None) -> list[tuple[float, float, str, str]]:
        out = []
        for r in self.results:
            if r.panel == panel and np.isfinite(r.p_adjusted) and np.isfinite(r.difference):
                kind = REST if r.p_adjusted >= self.alpha else UP if r.difference > 0 else DOWN
                out.append((r.difference, -np.log10(max(r.p_adjusted, 1e-300)),
                            f"{r.cluster} · {short(r.marker)}", kind))
        return out

    def names(self) -> dict[str, str]:
        a, b = self.results[0].groups if self.results else ("", "")
        return {UP: f"Higher in {a or '—'}", DOWN: f"Higher in {b or '—'}",
                REST: "Not significant"}

    @property
    def xlabel(self) -> str:
        a, b = self.results[0].groups if self.results else ("", "")
        return f"Difference in median ({a or '—'} − {b or '—'})"


class DifferenceHeatmap(figures.HeatmapContent):
    """Cluster x marker: the difference in median, a panel for each value
    of "within each". A heatmap to the figure dialog, with a ramp of its
    own: a difference has a zero, and the ramp its white there."""

    uses_colormap = False
    per_column = False

    def __init__(self, results: list[Result], alpha: float = 0.05):
        self.results, self.alpha = results, alpha

    # the parent's are of a table this does not have
    __eq__ = object.__eq__
    __hash__ = object.__hash__
    __repr__ = object.__repr__

    def panels(self) -> list[str | None]:
        return compare.ordered(r.panel for r in self.results) or [None]

    def clusters(self) -> list[str]:
        return list(dict.fromkeys(r.cluster for r in self.results))

    def markers(self) -> list[str]:
        return list(dict.fromkeys(r.marker for r in self.results))

    def cells(self, panel: str | None = None) -> tuple[np.ndarray, np.ndarray]:
        """(difference, adjusted p), clusters x markers. A cell that was
        not tested has no p, and no difference either: it is drawn grey."""
        row = {c: i for i, c in enumerate(self.clusters())}
        col = {m: i for i, m in enumerate(self.markers())}
        difference = np.full((len(row), len(col)), np.nan)
        p = np.full(difference.shape, np.nan)
        for r in self.results:
            if r.panel == panel and np.isfinite(r.p_adjusted):
                difference[row[r.cluster], col[r.marker]] = r.difference
                p[row[r.cluster], col[r.marker]] = r.p_adjusted
        return difference, p

    def limit(self) -> float:
        """The largest difference either way: the ends of the ramp, the
        same in every panel."""
        found = [abs(r.difference) for r in self.results
                 if np.isfinite(r.p_adjusted) and np.isfinite(r.difference)]
        return float(max(found)) if found and max(found) > 0 else 1.0

    @property
    def label(self) -> str:
        a, b = self.results[0].groups if self.results else ("", "")
        return f"Difference in median ({a or '—'} − {b or '—'})"


def _draw_differences(fig, c: DifferenceHeatmap, style: figures.FigureStyle) -> None:
    panels, clusters, markers = c.panels(), c.clusters(), c.markers()
    axes = fig.subplots(1, len(panels), sharey=True, squeeze=False)[0]
    cmap = figures.colormap(DIVERGING).with_extremes(bad=UNTESTED)
    limit = c.limit()
    turned = len(markers) * len(panels) > 6
    for k, (ax, panel) in enumerate(zip(axes, panels)):
        difference, p = c.cells(panel)
        im = ax.imshow(np.ma.masked_invalid(difference), cmap=cmap, vmin=-limit, vmax=limit,
                       aspect="auto")
        ax.set_xticks(range(len(markers)), [style.label(short(m)) for m in markers],
                      rotation=45 if turned else 0, ha="right" if turned else "center")
        ax.set_yticks(range(len(clusters)), [style.label(n) for n in clusters])
        if k == 0:
            ax.set_ylabel("Cluster" if style.ylabel is None else style.ylabel)
        ax.set_xlabel("" if style.xlabel is None else style.xlabel)
        if panel is not None:
            ax.set_title(style.label(panel))
        ax.tick_params(length=0)
        ax.spines[:].set_visible(style.frame == "box")
        if style.annotate:
            for i, j in zip(*np.nonzero(np.isfinite(p))):
                mark = stars(p[i, j], c.alpha)
                if mark:
                    dark = abs(difference[i, j]) > 0.6 * limit
                    # lowered: an asterisk sits high in its box
                    ax.text(j, i + 0.12, mark, ha="center", va="center",
                            fontsize=style.font_size, color="white" if dark else "black")
    if style.legend != "none":
        bar = fig.colorbar(im, ax=list(axes), shrink=0.7, aspect=22, pad=0.02,
                           location="bottom" if style.legend == "bottom" else "right")
        bar.set_label(style.legend_title if style.legend_title is not None else c.label)
        bar.outline.set_linewidth(0.5)


figures._DRAW[ExpressionVolcano] = _draw_volcano
figures._DRAW[DifferenceHeatmap] = _draw_differences
