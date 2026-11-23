"""
The volcano of a differential-abundance analysis: every cluster a point,
how much it changed across (log2 fold change) against how sure that is up
(-log10 of the adjusted p), so that what matters ends in the top corners.

A content like those of `figures`, and drawn by its machinery: importing
this module teaches `figures.save` and `figures.preview_png` to draw it.
What is a point, which side it is on and whether it is named is decided
here once, for the paper and for the screen.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import compare, figures
from .abundance import Result

UP, DOWN, REST = "up", "down", "rest"
COLORS = {UP: figures.CATEGORICAL[1], DOWN: figures.CATEGORICAL[0], REST: "#BBBBBB"}


@dataclass
class VolcanoContent:
    results: list[Result]
    #: an adjusted p below this is significant: coloured and named
    alpha: float = 0.05
    #: guide lines at this log2 fold change either side of 0; 1 is twofold
    fold: float = 1.0
    #: how the p values were adjusted: what the y axis is called
    correction: str = compare.BH

    def panels(self) -> list[str | None]:
        return compare.ordered(r.panel for r in self.results) or [None]

    def points(self, panel: str | None = None) -> list[tuple[float, float, str, str]]:
        """(x, y, cluster, UP / DOWN / REST) for each cluster of the panel
        that was tested. A p of zero is drawn at 1e-300."""
        out = []
        for r in self.results:
            if r.panel == panel and np.isfinite(r.p_adjusted) and np.isfinite(r.log2_fold):
                kind = REST if r.p_adjusted >= self.alpha else UP if r.log2_fold > 0 else DOWN
                out.append((r.log2_fold, -np.log10(max(r.p_adjusted, 1e-300)), r.cluster, kind))
        return out

    def names(self) -> dict[str, str]:
        """What each kind of point is called in the legend."""
        a, b = self.results[0].groups if self.results else ("", "")
        return {UP: f"More in {a or '—'}", DOWN: f"More in {b or '—'}",
                REST: "Not significant"}

    @property
    def xlabel(self) -> str:
        a, b = self.results[0].groups if self.results else ("", "")
        return f"log2 fold change ({a or '—'} / {b or '—'})"

    @property
    def ylabel(self) -> str:
        return {compare.NONE: "−log10 p", compare.BH: "−log10 FDR"}.get(
            self.correction, f"−log10 p ({self.correction})")

    def limits(self) -> tuple[float, float]:
        """(half-width of x, top of y): the same for every panel, with the
        guide lines always inside."""
        pts = [p for panel in self.panels() for p in self.points(panel)]
        x = max([abs(p[0]) for p in pts] + [self.fold]) * 1.15
        y = max([p[1] for p in pts] + [-np.log10(self.alpha)]) * 1.15
        return float(x), float(y)


def _draw_volcano(fig, c: VolcanoContent, style: figures.FigureStyle) -> None:
    panels = c.panels()
    axes = fig.subplots(1, len(panels), sharex=True, sharey=True, squeeze=False)[0]
    half, top = c.limits()
    names = c.names()
    colour = {k: style.color(names[k], COLORS[k]) for k in COLORS}
    for k, (ax, panel) in enumerate(zip(axes, panels)):
        ax.axhline(-np.log10(c.alpha), color="#888888", linewidth=0.6, linestyle="--", zorder=1)
        for x in (-c.fold, c.fold):
            ax.axvline(x, color="#888888", linewidth=0.6, linestyle=":", zorder=1)
        pts = c.points(panel)
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=(style.point_size * 1.6) ** 2,
                   c=[colour[p[3]] for p in pts], linewidths=0, zorder=2)
        for x, y, cluster, kind in pts:
            if kind != REST:
                ax.annotate(style.label(cluster), (x, y), xytext=(3, 3),
                            textcoords="offset points", fontsize=style.font_size - 1)
        ax.set_xlim(-half, half)
        ax.set_ylim(0, top)
        if panel is not None:
            ax.set_title(style.label(panel))
        if len(panels) == 1:
            ax.set_xlabel(c.xlabel if style.xlabel is None else style.xlabel)
        if k == 0:
            ax.set_ylabel(c.ylabel if style.ylabel is None else style.ylabel)
        figures._frame(ax, style)
    if len(panels) > 1:      # one label under them all: side by side they run into each other
        fig.supxlabel(c.xlabel if style.xlabel is None else style.xlabel,
                      fontsize=axes[0].xaxis.label.get_size())
    handles = [figures._dot(colour[k], style.label(names[k])) for k in COLORS]
    figures._legend(fig, list(axes), handles, style, "")


figures._DRAW[VolcanoContent] = _draw_volcano
