"""
The 2D plot: events as points, coloured by a category or by a value, in
one panel or split into one panel per category with the rest of the events
drawn faint behind — so where a group sits is read against where everyone
sits.

Points are drawn in a fixed shuffled order. In file order the last sample
would be painted over all the others and a map coloured by group would
show mostly the group that happened to be concatenated last.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui

from ..figures import Categories, MapContent, Values, categories_of  # noqa: F401
from ..gating import points_in_polygon
from . import palette, theme


class _LassoBox(pg.ViewBox):
    """A view that pans and zooms as usual, and with Shift held draws a
    freehand outline instead."""

    def __init__(self, owner: "MapView"):
        super().__init__()
        self.owner = owner
        self._points: list[tuple[float, float]] = []
        self._curve = None

    def mouseDragEvent(self, ev, axis=None):
        shift = bool(ev.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier)
        if not (shift or self._points) or ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return super().mouseDragEvent(ev, axis)
        ev.accept()
        p = self.mapSceneToView(ev.scenePos())
        if ev.isStart():
            self._points = []
            pen = pg.mkPen(theme.foreground(), width=1.5, style=QtCore.Qt.PenStyle.DashLine)
            self._curve = pg.PlotCurveItem(pen=pen)
            self.addItem(self._curve, ignoreBounds=True)
        self._points.append((p.x(), p.y()))
        xs, ys = zip(*(self._points + self._points[:1]))
        self._curve.setData(np.array(xs), np.array(ys))
        if ev.isFinish():
            poly = np.array(self._points)
            self.removeItem(self._curve)
            self._points, self._curve = [], None
            if len(poly) >= 3:
                self.owner.lasso(poly)


class MapView(pg.GraphicsLayoutWidget):
    #: a boolean mask over the events: what a lasso enclosed (Shift+drag)
    sigSelected = QtCore.pyqtSignal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        #: events drawn ringed, on top: the current selection
        self.selected: np.ndarray | None = None
        self.panels: list[pg.PlotItem] = []
        self.point_size = 4.0
        self.opacity = 0.8
        self._args = None
        self._ticks = (None, None)
        self._columns = None
        self.setBackground(theme.background())

    # --- what to draw ------------------------------------------------------

    def show_points(self, x, y, xlabel: str, ylabel: str, color=None,
                    split: Categories | None = None, xticks=None, yticks=None,
                    columns: int | None = None) -> None:
        self._args = (np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64),
                      xlabel, ylabel, color, split)
        self._ticks = (xticks, yticks)
        self._columns = columns
        self.redraw()

    def content(self) -> MapContent | None:
        """What is on the plot, for a figure."""
        if self._args is None:
            return None
        x, y, xlabel, ylabel, color, split = self._args
        return MapContent(x, y, xlabel, ylabel, color, split, *self._ticks,
                          columns=self._columns)

    def lasso(self, polygon: np.ndarray) -> None:
        if self._args is None:
            return
        x, y = self._args[0], self._args[1]
        self.sigSelected.emit(points_in_polygon(x, y, polygon))

    def set_selected(self, mask: np.ndarray | None) -> None:
        self.selected = mask
        self.redraw()

    def clear_points(self) -> None:
        self._args = None
        self.selected = None
        self.clear()
        self.panels = []

    def set_point_size(self, size: float) -> None:
        self.point_size = float(size)
        self.redraw()

    def set_opacity(self, opacity: float) -> None:
        self.opacity = float(opacity)
        self.redraw()

    def retheme(self) -> None:
        self.setBackground(theme.background())
        self.redraw()

    # --- drawing -----------------------------------------------------------

    def _brushes(self, color, n: int) -> tuple[list, np.ndarray]:
        """A small table of brushes and, per event, which one."""
        alpha = int(round(255 * self.opacity))
        if isinstance(color, Categories):
            table = []
            for c in color.colors:
                q = QtGui.QColor(c)
                q.setAlpha(alpha)
                table.append(pg.mkBrush(q))
            return table, color.codes
        if isinstance(color, Values):
            lo, hi = color.limits()
            lut = palette.ramp()
            table = [pg.mkBrush(int(r), int(g), int(b), alpha) for r, g, b, _ in lut]
            v = np.nan_to_num(color.values, nan=lo)
            idx = np.clip((v - lo) / (hi - lo) * (len(lut) - 1), 0, len(lut) - 1)
            return table, idx.astype(np.int64)
        q = QtGui.QColor(theme.accent())
        q.setAlpha(alpha)
        return [pg.mkBrush(q)], np.zeros(n, dtype=np.int64)

    def redraw(self) -> None:
        self.clear()
        self.panels = []
        if self._args is None:
            return
        x, y, xlabel, ylabel, color, split = self._args
        n = x.size
        if n == 0:
            return
        table, which = self._brushes(color, n)
        order = np.random.default_rng(0).permutation(n)
        faint = QtGui.QColor(theme.faint_axis())
        faint.setAlpha(110)

        groups = ([(None, np.ones(n, dtype=bool))] if split is None else
                  [(name, split.codes == i) for i, name in enumerate(split.names)])
        cols = self._columns or (int(np.ceil(np.sqrt(len(groups)))) if len(groups) > 1 else 1)
        cols = max(1, min(cols, len(groups)))
        for k, (name, member) in enumerate(groups):
            plot = self.addPlot(row=k // cols, col=k % cols, viewBox=_LassoBox(self))
            plot.setMenuEnabled(False)
            plot.hideButtons()
            theme.style_axes(plot)
            plot.setLabel("bottom", xlabel)
            plot.setLabel("left", ylabel)
            for side, ticks in zip(("bottom", "left"), self._ticks):
                if ticks:
                    plot.getAxis(side).setTicks([ticks])
            if name is not None:
                plot.setTitle(f"{name}  <span style='color:{theme.muted()}'>"
                              f"{int(member.sum()):,d}</span>",
                              color=theme.foreground(), size="10pt")
                rest = order[~member[order]]
                plot.addItem(pg.ScatterPlotItem(
                    x[rest], y[rest], size=self.point_size, pen=None,
                    brush=pg.mkBrush(faint), pxMode=True))
            rows = order[member[order]]
            plot.addItem(pg.ScatterPlotItem(
                x[rows], y[rows], size=self.point_size, pen=None, pxMode=True,
                brush=[table[i] for i in which[rows]]))
            if self.selected is not None and self.selected.size == n:
                rows = np.flatnonzero(self.selected & member)
                plot.addItem(pg.ScatterPlotItem(
                    x[rows], y[rows], size=self.point_size + 3, pxMode=True,
                    pen=pg.mkPen(theme.foreground(), width=1.1), brush=None))
            if self.panels:
                plot.setXLink(self.panels[0])
                plot.setYLink(self.panels[0])
            self.panels.append(plot)

        if isinstance(color, Values):
            lo, hi = color.limits()
            bar = pg.ColorBarItem(values=(lo, hi), colorMap=palette.colormap(),
                                  interactive=False, width=12, label=color.title)
            for ax in (bar.axis,):
                ax.setPen(pg.mkPen(theme.axis()))
                ax.setTextPen(pg.mkPen(theme.foreground()))
            self.addItem(bar, row=0, col=cols, rowspan=max(1, -(-len(groups) // cols)))
        if self.panels:
            self.panels[0].enableAutoRange()

    def snapshot(self) -> QtGui.QImage:
        return self.grab().toImage()
