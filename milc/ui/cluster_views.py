"""
The three pictures of a clustering beside the map: what each cluster
expresses (heatmap), where it is found (frequencies by group) and how the
map's nodes relate (the FlowSOM tree). Each says which cluster was
clicked; none decides what a click means.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui

from .. import cluster_stats
from ..figures import FrequencyContent, HeatmapContent, TreeContent
from . import palette, theme


def _short(name: str) -> str:
    """`CD8 FITC-A` → `CD8`: an axis has no room for the fluorochrome."""
    return name.split()[0] if name else name


class _Clickable(pg.PlotWidget):
    #: the cluster (its name in the table, an int) under a click
    sigCluster = QtCore.pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setBackground(theme.background())
        self.setMenuEnabled(False)
        self.hideButtons()
        self.setMouseEnabled(False, False)
        self.scene().sigMouseClicked.connect(self._clicked)

    def _clicked(self, ev) -> None:
        if ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return
        p = self.plotItem.vb.mapSceneToView(ev.scenePos())
        cluster = self.cluster_at(p.x(), p.y())
        if cluster is not None:
            self.sigCluster.emit(int(cluster))

    def cluster_at(self, x: float, y: float):
        return None


class HeatmapView(_Clickable):
    """Clusters down, markers across."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.table = None

    def show_table(self, table, per_column: bool = True, selected: int | None = None) -> None:
        self.table = table
        self._per_column = per_column
        self.clear()
        self.setBackground(theme.background())
        if table is None:
            return
        v = table.values
        rows, cols = v.shape
        shade = cluster_stats.scaled(table, per_column)
        lut = palette.ramp()
        img = pg.ImageItem(shade.T[:, ::-1], levels=(0.0, 1.0), lut=lut)
        img.setRect(QtCore.QRectF(0, 0, cols, rows))
        self.addItem(img)
        if rows * cols <= 400:
            font = QtGui.QFont()
            font.setPointSize(9)
            for r in range(rows):
                for c in range(cols):
                    rgb = lut[int(np.clip(shade[r, c], 0, 1) * (len(lut) - 1))]
                    dark = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2] < 140
                    t = pg.TextItem(f"{v[r, c]:.0f}", anchor=(0.5, 0.5),
                                    color="#ffffff" if dark else "#111111")
                    t.setFont(font)
                    t.setPos(c + 0.5, rows - r - 0.5)
                    self.addItem(t)
        theme.style_axes(self.plotItem)
        self.getAxis("bottom").setTicks([[(c + 0.5, _short(n))
                                          for c, n in enumerate(table.columns)]])
        self.getAxis("left").setTicks([[(rows - r - 0.5, n)
                                        for r, n in enumerate(table.rows)]])
        self.setLabel("left", "Cluster")
        self.setLabel("bottom", "")
        if selected is not None and str(selected) in table.rows:
            r = table.rows.index(str(selected))
            box = pg.QtWidgets.QGraphicsRectItem(0, rows - r - 1, cols, 1)
            pen = pg.mkPen(theme.foreground(), width=2.5)
            pen.setCosmetic(True)
            box.setPen(pen)
            self.addItem(box)
        self.setXRange(0, cols, padding=0)
        self.setYRange(0, rows, padding=0)

    def content(self):
        return (HeatmapContent(self.table, self._per_column)
                if self.table is not None else None)

    def cluster_at(self, x, y):
        if self.table is None:
            return None
        rows = len(self.table.rows)
        r = rows - 1 - int(np.floor(y))
        if 0 <= r < rows and 0 <= x <= len(self.table.columns):
            return int(self.table.rows[r])
        return None


class FrequencyView(pg.GraphicsLayoutWidget):
    """For each cluster, the share of each sample's events in it: a bar
    per group at the mean, a whisker at one standard deviation, and every
    sample as a dot — five animals are five dots, not just a bar.

    With a second subdivision there is one panel for each of its values,
    on one scale: the groups compared within each day, the days side by
    side."""

    sigCluster = QtCore.pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.freq = None
        self._groups: list[str] = []
        self._facets: list[str] | None = None
        self.panels: list[pg.PlotItem] = []
        self.setBackground(theme.background())
        self.scene().sigMouseClicked.connect(self._clicked)

    def content(self):
        return (FrequencyContent(self.freq, self._groups, self._facets)
                if self.freq is not None else None)

    def show_frequencies(self, freq, groups: list[str], selected: int | None = None,
                         facets: list[str] | None = None) -> None:
        self.freq = freq
        self._groups = list(groups)
        self._facets = list(facets) if facets else None
        self.clear()
        self.panels = []
        self.setBackground(theme.background())
        if freq is None:
            return
        names = cluster_stats.ordered(groups)
        colors = palette.categorical(len(names))
        n_clusters = len(freq.rows)
        width = 0.8 / max(len(names), 1)
        top = np.nanmax(freq.values) if np.isfinite(freq.values).any() else 1.0
        values = cluster_stats.ordered(self._facets) if self._facets else [None]
        for k, value in enumerate(values):
            plot = self.addPlot(row=0, col=k)
            plot.setMenuEnabled(False)
            plot.hideButtons()
            plot.setMouseEnabled(False, False)
            theme.style_axes(plot)
            if value is not None:
                plot.setTitle(str(value) or "—", color=theme.foreground(), size="10pt")
            if selected is not None and str(selected) in freq.rows:
                i = freq.rows.index(str(selected))
                band = pg.LinearRegionItem((i - 0.48, i + 0.48), movable=False,
                                           brush=pg.mkBrush(QtGui.QColor(theme.line())),
                                           pen=pg.mkPen(None))
                band.setZValue(-10)
                plot.addItem(band)
            legend = None
            if k == len(values) - 1:
                legend = plot.addLegend(offset=(-8, 8))
                legend.setBrush(theme.legend_brush())
                legend.setPen(theme.legend_pen())
                legend.setLabelTextColor(theme.foreground())
            rng = np.random.default_rng(0)
            for j, (g, color) in enumerate(zip(names, colors)):
                cols = [c for c, x in enumerate(groups)
                        if x == g and (value is None or self._facets[c] == value)]
                vals = freq.values[:, cols]
                vals = vals[:, ~np.isnan(vals).all(axis=0)] if vals.size else vals
                fill = QtGui.QColor(color)
                fill.setAlpha(150)
                if vals.shape[1] == 0:
                    if legend is not None:
                        ghost = pg.BarGraphItem(x=[0], height=[0], width=0,
                                                brush=pg.mkBrush(fill), pen=pg.mkPen(color))
                        legend.addItem(ghost, f"{g or '—'} (n=0)")
                    continue
                x = np.arange(n_clusters) - 0.4 + width * (j + 0.5)
                mean = vals.mean(axis=1)
                bars = pg.BarGraphItem(x=x, height=mean, width=width * 0.9,
                                       brush=pg.mkBrush(fill), pen=pg.mkPen(color))
                plot.addItem(bars)
                if legend is not None:
                    legend.addItem(bars, f"{g or '—'} (n={vals.shape[1]})")
                if vals.shape[1] > 1:
                    sd = vals.std(axis=1, ddof=1)
                    plot.addItem(pg.ErrorBarItem(x=x, y=mean, top=sd, bottom=sd,
                                                 beam=width * 0.3,
                                                 pen=pg.mkPen(theme.foreground(), width=1)))
                jitter = (rng.random(vals.shape) - 0.5) * width * 0.5
                plot.addItem(pg.ScatterPlotItem(
                    (x[:, None] + jitter).ravel(), vals.ravel(), size=5,
                    pen=pg.mkPen(theme.foreground(), width=0.8),
                    brush=pg.mkBrush(theme.background())))
            plot.getAxis("bottom").setTicks([[(i, n) for i, n in enumerate(freq.rows)]])
            plot.setLabel("bottom", "Cluster")
            if k == 0:
                plot.setLabel("left", "% of the sample's events")
            plot.setXRange(-0.6, n_clusters - 0.4, padding=0)
            plot.setYRange(0, top * 1.12, padding=0)
            if self.panels:
                plot.setYLink(self.panels[0])
            self.panels.append(plot)

    def _clicked(self, ev) -> None:
        if ev.button() != QtCore.Qt.MouseButton.LeftButton or self.freq is None:
            return
        for plot in self.panels:
            if plot.vb.sceneBoundingRect().contains(ev.scenePos()):
                cluster = self.cluster_at(plot.vb.mapSceneToView(ev.scenePos()).x())
                if cluster is not None:
                    self.sigCluster.emit(int(cluster))
                return

    def cluster_at(self, x, y=0.0):
        if self.freq is None:
            return None
        i = int(round(x))
        return int(self.freq.rows[i]) if 0 <= i < len(self.freq.rows) else None


class TreeView(_Clickable):
    """The FlowSOM tree: a node per map unit that holds events, as large
    as it is populated, coloured by its metacluster or by a marker."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tree = None
        self.setAspectLocked(True)
        self.hideAxis("bottom")
        self.hideAxis("left")

    def show_tree(self, tree, color_by: str | None = None,
                  selected: int | None = None) -> None:
        """`color_by` is one of `tree["markers"]`, or None for the
        metacluster."""
        self.tree = tree
        self._color_by = color_by if tree and color_by in tree["markers"] else None
        self.clear()
        self.setBackground(theme.background())
        if tree is None:
            return
        pos, sizes, cluster = tree["pos"], tree["sizes"], tree["cluster"]
        used = np.flatnonzero(sizes > 0)
        self._used = used
        if used.size == 0:
            return
        if len(tree["edges"]):
            pairs = pos[tree["edges"].ravel()]
            self.addItem(pg.PlotCurveItem(pairs[:, 0], pairs[:, 1], connect="pairs",
                                          pen=pg.mkPen(theme.faint_axis(), width=1.2)))
        if color_by in tree["markers"]:
            m = tree["medians"][used, tree["markers"].index(color_by)]
            lo, hi = np.nanmin(m), np.nanmax(m)
            lut = palette.ramp()
            idx = ((m - lo) / (hi - lo if hi > lo else 1.0) * (len(lut) - 1)).astype(int)
            brushes = [pg.mkBrush(*[int(c) for c in lut[i][:3]]) for i in idx]
        else:
            ids = sorted(set(int(c) for c in cluster[used]))
            colors = dict(zip(ids, palette.categorical(len(ids))))
            brushes = [pg.mkBrush(colors[int(c)]) for c in cluster[used]]
        radius = 7 + 22 * np.sqrt(sizes[used] / sizes[used].max())
        pens = [pg.mkPen(theme.foreground(), width=2.4) if selected is not None
                and int(c) == selected else pg.mkPen(theme.background(), width=1)
                for c in cluster[used]]
        self.addItem(pg.ScatterPlotItem(pos[used, 0], pos[used, 1], size=radius,
                                        brush=brushes, pen=pens, pxMode=True))
        self.autoRange(padding=0.08)

    def content(self):
        return TreeContent(self.tree, self._color_by) if self.tree is not None else None

    def cluster_at(self, x, y):
        if self.tree is None or self._used.size == 0:
            return None
        p = self.tree["pos"][self._used]
        d = np.hypot(p[:, 0] - x, p[:, 1] - y)
        i = int(np.argmin(d))
        span = np.nanmax(p, axis=0) - np.nanmin(p, axis=0)
        # a click in the empty part of the picture selects nothing
        if d[i] > 0.06 * max(float(span.max()), 1e-9) + 1e-9:
            return None
        return int(self.tree["cluster"][self._used[i]])
