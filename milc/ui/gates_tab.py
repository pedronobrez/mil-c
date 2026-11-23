"""
The gates: open one sample, look at any two of its parameters or at the
histogram of one, follow the gating tree down, and draw, move, rename or
delete a gate.

A gate drawn here belongs to the population being looked at, as in
FlowJo, and — unless the box is unticked — to every sample at once, as a
gate put on a group does. What is drawn is kept beside the files
(`MIL-C gates.json`) and takes the place of the workspace's from then on.

A marker is shown on FlowJo's biexponential scale and a polygon drawn
there is tested there: a side that is straight on the screen is the side
the events are compared with.
"""

from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets

from .. import density, figures, order, views
from ..gating import (AND, NOT, OR, QUADRANTS, BooleanGate, EllipseGate, GateNode,
                      PolygonGate, QuadrantGate, RectangleGate, depends_on)
from ..gating import TimeGate, time_node
from ..session import Session
from ..transforms import Arcsinh, Biex, Log, Logicle, decade_ticks, display_range, same
from . import palette, style, theme
from .flow import captioned, flow, together
from .view_dialogs import BackgateDialog, MatrixDialog, OverlayChoice, draw_distributions

HISTOGRAM = "#histogram"
PSEUDOCOLOUR, CONTOUR = "pseudocolour", "contour"
BY_SAMPLE = "#sample"
MAX_POINTS = 40_000


class _DrawBox(pg.ViewBox):
    """Pans and zooms as usual; in a drawing mode the mouse draws."""

    sigRectangle = QtCore.pyqtSignal(float, float, float, float)
    sigEllipse = QtCore.pyqtSignal(float, float, float, float)
    sigPolygon = QtCore.pyqtSignal(object)
    sigPoint = QtCore.pyqtSignal(float, float)

    def __init__(self):
        super().__init__()
        self.mode: str | None = None
        self._points: list[tuple[float, float]] = []
        self._shape = None

    def set_mode(self, mode: str | None) -> None:
        self.cancel()
        self.mode = mode

    def cancel(self) -> None:
        if self._shape is not None:
            self.removeItem(self._shape)
        self._shape, self._points = None, []

    def _outline(self, xs, ys) -> None:
        if self._shape is None:
            pen = pg.mkPen(theme.foreground(), width=1.6, style=QtCore.Qt.PenStyle.DashLine)
            self._shape = pg.PlotCurveItem(pen=pen)
            self.addItem(self._shape, ignoreBounds=True)
        self._shape.setData(np.asarray(xs, dtype=float), np.asarray(ys, dtype=float))

    def mouseDragEvent(self, ev, axis=None):
        if self.mode not in ("rect", "ellipse") or ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return super().mouseDragEvent(ev, axis)
        ev.accept()
        a = self.mapSceneToView(ev.buttonDownScenePos())
        b = self.mapSceneToView(ev.scenePos())
        if self.mode == "ellipse":
            t = np.linspace(0, 2 * np.pi, 73)
            self._outline((a.x() + b.x()) / 2 + abs(b.x() - a.x()) / 2 * np.cos(t),
                          (a.y() + b.y()) / 2 + abs(b.y() - a.y()) / 2 * np.sin(t))
        else:
            self._outline([a.x(), b.x(), b.x(), a.x(), a.x()],
                          [a.y(), a.y(), b.y(), b.y(), a.y()])
        if ev.isFinish():
            signal = self.sigEllipse if self.mode == "ellipse" else self.sigRectangle
            self.cancel()
            if a.x() != b.x():
                signal.emit(min(a.x(), b.x()), min(a.y(), b.y()),
                            max(a.x(), b.x()), max(a.y(), b.y()))

    def mouseClickEvent(self, ev):
        if self.mode not in ("poly", "quad") or ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return super().mouseClickEvent(ev)
        ev.accept()
        p = self.mapSceneToView(ev.scenePos())
        if self.mode == "quad":
            self.sigPoint.emit(p.x(), p.y())
            return
        if ev.double() and len(self._points) >= 3:
            points = np.array(self._points)
            self.cancel()
            self.sigPolygon.emit(points)
            return
        if not ev.double():
            self._points.append((p.x(), p.y()))
            xs, ys = zip(*self._points)
            self._outline(xs, ys)


class GatesTab(QtWidgets.QWidget):
    sigStatus = QtCore.pyqtSignal(str)
    #: the gates of one sample or of all changed: counts elsewhere are stale
    sigGatesChanged = QtCore.pyqtSignal()
    #: an axis was given another scale: other plots of that parameter are stale
    sigScalesChanged = QtCore.pyqtSignal()

    def __init__(self, session: Session, parent=None):
        super().__init__(parent)
        self.session = session
        self.entry = None
        self.sample = None
        self.node: GateNode | None = None
        self._filling = False
        self._keep_axes = False
        self._overlays: list = []

        # --- left: the samples and the tree
        self.samples = QtWidgets.QListWidget()
        self.samples.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.samples.setToolTip(
            "Click a sample to look at it. Select several — Shift or Cmd/Ctrl-click — "
            "to overlay them: histograms on a histogram, contours on two parameters.")
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderLabels(["Population", "Events", "% parent"])
        self.tree.setRootIsDecorated(True)
        # six gates deep is ordinary, and at the default indentation the
        # names at the bottom of the tree had no room left at all
        self.tree.setIndentation(11)
        self.tree.header().setStretchLastSection(False)
        self.tree.header().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.rect_button = QtWidgets.QPushButton("Rectangle")
        self.poly_button = QtWidgets.QPushButton("Polygon")
        for b, tip in ((self.rect_button, "Drag on the plot. On a histogram, a range."),
                       (self.poly_button, "Click each corner; double-click to close.")):
            b.setCheckable(True)
            b.setToolTip(tip)
        self.ellipse_button = QtWidgets.QPushButton("Ellipse")
        self.ellipse_button.setCheckable(True)
        self.ellipse_button.setToolTip("Drag the box it fits in; then turn and stretch it by its handles.")
        self.quad_button = QtWidgets.QPushButton("Quadrant")
        self.quad_button.setCheckable(True)
        self.quad_button.setToolTip("Click where the two lines cross: four populations, "
                                    "named by which side of each line they are on.")
        self.combine_button = QtWidgets.QPushButton("Combine…")
        self.combine_button.setToolTip("A population made of others of this one: the events "
                                       "in all of them, in any of them, or in none of them.")
        self.rename_button = QtWidgets.QPushButton("Rename")
        self.delete_button = QtWidgets.QPushButton("Delete")
        self.reset_button = QtWidgets.QPushButton("Back to the workspace's gates")
        self.reset_button.setToolTip("Forget every gate drawn or moved here")
        self.everyone = QtWidgets.QCheckBox("Apply to every sample")
        self.everyone.setChecked(True)
        self.everyone.setToolTip(
            "Ticked, a gate drawn, moved or deleted on this sample is drawn, moved or "
            "deleted on all of them, as a gate on a FlowJo group is. Unticked, only here.")
        tools = QtWidgets.QGridLayout()
        tools.addWidget(self.rect_button, 0, 0)
        tools.addWidget(self.poly_button, 0, 1)
        tools.addWidget(self.ellipse_button, 0, 2)
        tools.addWidget(self.quad_button, 1, 0)
        tools.addWidget(self.combine_button, 1, 1, 1, 2)
        tools.addWidget(self.rename_button, 2, 0)
        tools.addWidget(self.delete_button, 2, 1, 1, 2)
        left = QtWidgets.QVBoxLayout()
        left.addWidget(self._eyebrow("Samples"))
        left.addWidget(self.samples, 2)
        left.addWidget(self._eyebrow("Populations"))
        left.addWidget(self.tree, 3)
        left.addWidget(self._eyebrow("New gate on the population shown"))
        left.addLayout(tools)
        left.addWidget(self.everyone)
        left.addWidget(self.reset_button)
        left_panel = QtWidgets.QWidget()
        left_panel.setLayout(left)
        left_panel.setFixedWidth(370)

        # --- right: the plot
        self.x_axis = QtWidgets.QComboBox()
        self.y_axis = QtWidgets.QComboBox()
        self.caption = QtWidgets.QLabel("")
        self.caption.setProperty("role", "caption")
        self.caption.setWordWrap(True)
        self.kind = QtWidgets.QComboBox()
        self.kind.addItem("Pseudocolour", PSEUDOCOLOUR)
        self.kind.addItem("Contour", CONTOUR)
        self.kind.setToolTip("Contour lines enclose 90, 75, 50, 25 and 10% of the events; "
                             "the events outside the outermost are drawn as dots")
        self.colour_by = QtWidgets.QComboBox()
        self.colour_by.setToolTip("With several samples overlaid: a colour for each "
                                  "sample, or for each value of one of your columns")
        self.overlay = OverlayChoice()
        #: histograms one above the other instead of on one baseline
        self.offset = self.overlay.staggered
        self.matrix_button = QtWidgets.QPushButton("Matrix…")
        self.matrix_button.setToolTip("Every pair of parameters of the population shown, "
                                      "in one grid: where spillover shows")
        self.backgate_button = QtWidgets.QPushButton("Backgate…")
        self.backgate_button.setToolTip("Where the events of the population shown fall on "
                                        "each gate above it")
        self.figure_button = QtWidgets.QPushButton("Figure…")
        self.x_scale = QtWidgets.QToolButton()
        self.y_scale = QtWidgets.QToolButton()
        for b, which in ((self.x_scale, "x"), (self.y_scale, "y")):
            b.setText("Scale…")
            b.setToolTip("How this parameter's axis is drawn: linear, biexponential, "
                         "Logicle, arcsinh or logarithmic, and with what settings")
            b.clicked.connect(lambda _=False, w=which: self._ask_scale(w))
        self.compensation_button = QtWidgets.QPushButton("Compensation…")
        self.compensation_button.setToolTip(
            "The spillover matrix each sample is compensated with — its file's, the "
            "workspace's, or one edited here — and how to change it for one sample, "
            "some, or all")
        self.quality_button = QtWidgets.QPushButton("Quality…")
        self.quality_button.setToolTip("Look for the stretches of the acquisition that went "
                                       "wrong — a clog, a bubble, a drift — and remove them")
        self.export_button = QtWidgets.QPushButton("Export…")
        self.export_button.setToolTip("The events of the population shown, of the samples "
                                      "selected or of all, as FCS or CSV files")
        bar = flow(captioned("X", self.x_axis, self.x_scale),
                   captioned("Y", self.y_axis, self.y_scale),
                   captioned("Plot", self.kind), captioned("Colour by", self.colour_by),
                   self.overlay,
                   together(self.quality_button, self.compensation_button,
                            self.export_button, self.figure_button))
        self.box = _DrawBox()
        self.plot = pg.PlotWidget(viewBox=self.box)
        self.plot.setBackground(theme.background())
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.empty = QtWidgets.QLabel("Open a folder in the Samples tab, then choose a sample here.")
        self.empty.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.empty.setWordWrap(True)      # or the sentence sets the window's width
        self.empty.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                                 QtWidgets.QSizePolicy.Policy.Preferred)
        self.empty.setProperty("role", "hint")
        self.stack = QtWidgets.QStackedWidget()
        self.stack.addWidget(self.empty)
        self.stack.addWidget(self.plot)
        right = QtWidgets.QVBoxLayout()
        self.review = QtWidgets.QLabel("")
        self.review.setProperty("role", "warning")
        self.review.setWordWrap(True)
        self.reviewed_button = QtWidgets.QPushButton("Reviewed")
        self.reviewed_button.setToolTip("The gates were looked at and are where they "
                                        "should be: stop saying so")
        self.reviewed_button.clicked.connect(self.reviewed)
        self.review_bar = QtWidgets.QWidget()
        said = QtWidgets.QHBoxLayout(self.review_bar)
        said.setContentsMargins(0, 0, 0, 0)
        said.addWidget(self.review, 1)
        said.addWidget(self.reviewed_button)
        self.review_bar.hide()
        right.addWidget(self.review_bar)
        right.addWidget(bar)
        looks = QtWidgets.QHBoxLayout()
        looks.addWidget(self.caption, 1)
        looks.addWidget(self.matrix_button)
        looks.addWidget(self.backgate_button)
        self.model_button = QtWidgets.QPushButton("Model…")
        self.model_button.setToolTip(
            "Fit a model to the population shown, in every sample selected: the "
            "generations of a dye's dilution, the phases of the cell cycle, or a "
            "parameter through the time of the acquisition")
        models = QtWidgets.QMenu(self.model_button)
        for kind in ("Proliferation", "Cell cycle", "Kinetics"):
            models.addAction(kind + "…").triggered.connect(
                lambda _checked=False, k=kind: self.open_model(k))
        self.model_button.setMenu(models)
        looks.addWidget(self.model_button)
        right.addLayout(looks)
        right.addWidget(self.stack, 1)

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.addWidget(left_panel)
        layout.addLayout(right, 1)

        self.samples.currentRowChanged.connect(self._sample_chosen)
        self.samples.itemSelectionChanged.connect(self._selection_changed)
        self.kind.currentIndexChanged.connect(self._axes_changed)
        self.colour_by.currentIndexChanged.connect(self._axes_changed)
        self.overlay.sigChanged.connect(self._axes_changed)
        self.matrix_button.clicked.connect(self.open_matrix)
        self.backgate_button.clicked.connect(self.open_backgate)
        self.figure_button.clicked.connect(self.make_figure)
        self.compensation_button.clicked.connect(self._ask_compensation)
        self.export_button.clicked.connect(self._ask_export)
        self.quality_button.clicked.connect(self._ask_quality)
        self.tree.currentItemChanged.connect(self._node_chosen)
        self.x_axis.currentIndexChanged.connect(self._axes_changed)
        self.y_axis.currentIndexChanged.connect(self._axes_changed)
        self._tools = {"rect": self.rect_button, "poly": self.poly_button,
                       "ellipse": self.ellipse_button, "quad": self.quad_button}
        for mode, button in self._tools.items():
            button.toggled.connect(lambda on, m=mode: self._tool(m, on))
        self.box.sigRectangle.connect(self._rectangle_drawn)
        self.box.sigEllipse.connect(self._ellipse_drawn)
        self.box.sigPolygon.connect(self._polygon_drawn)
        self.box.sigPoint.connect(self._quadrant_drawn)
        self.combine_button.clicked.connect(self._ask_combination)
        self.rename_button.clicked.connect(self.rename_selected)
        self.delete_button.clicked.connect(self.delete_selected)
        self.reset_button.clicked.connect(self.reset)
        cancel = QtGui.QShortcut(QtGui.QKeySequence("Escape"), self)
        cancel.activated.connect(self._stop_drawing)
        self.refresh()

    @staticmethod
    def _eyebrow(text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setFont(style.eyebrow_font())
        label.setProperty("role", "caption")
        label.setContentsMargins(0, 8, 0, 2)
        return label

    # --- the folder -----------------------------------------------------------

    def refresh(self) -> None:
        """After a folder or a workspace was opened."""
        self.show_review()
        self._filling = True
        self.samples.clear()
        self.entries = self.session.usable()
        for e in self.entries:
            self.samples.addItem(f"{e.name}    {e.events:,d}")
        keep = self.colour_by.currentData()
        self.colour_by.clear()
        self.colour_by.addItem("Sample", BY_SAMPLE)
        for column in self.session.columns:
            self.colour_by.addItem(column, column)
        self.colour_by.setCurrentIndex(max(self.colour_by.findData(keep), 0))
        self._filling = False
        self._content = None
        self.entry = self.sample = self.node = None
        self.tree.clear()
        self.stack.setCurrentIndex(0)
        self.reset_button.setEnabled(any(e.edited for e in self.entries))

    def showEvent(self, ev) -> None:
        super().showEvent(ev)
        # a file is read when someone wants to see it, not when the folder opens
        if self.entry is None and self.samples.count():
            self.samples.setCurrentRow(0)

    def _sample_chosen(self, row: int) -> None:
        if self._filling or row < 0 or row >= len(self.entries):
            return
        keep = self.node.path if self.node is not None else ""
        self.entry = self.entries[row]
        self.sigStatus.emit(f"Reading {self.entry.name}…")
        QtWidgets.QApplication.processEvents()
        self.sample = self.session.load(self.entry)
        self.sigStatus.emit("")
        self._fill_axes()
        self._fill_tree(select=keep)

    # --- axes -----------------------------------------------------------------

    def _dim(self, channel: str) -> str:
        """The name a gate uses for a channel: `Comp-…` where it is
        compensated, as FlowJo writes it."""
        comp = self.sample.compensation
        return "Comp-" + channel if comp is not None and channel in comp.channels else channel

    def _scale(self, dim: str):
        return self.session.scale_of(dim)

    def _display(self, dim: str, sample=None) -> np.ndarray:
        x = (sample or self.sample).column(dim).astype(np.float64)
        scale = self._scale(dim)
        return scale(x) if scale is not None else x

    # --- several samples at once ------------------------------------------------

    def overlaid(self) -> list:
        """The samples selected together, in the list's order; just the
        one being looked at when nothing else is selected."""
        rows = sorted({i.row() for i in self.samples.selectedIndexes()})
        chosen = [self.entries[r] for r in rows if r < len(self.entries)]
        return chosen if len(chosen) > 1 else ([self.entry] if self.entry else [])

    def select_samples(self, names: list[str]) -> None:
        """Overlay these; the first becomes the one whose gates are shown."""
        rows = [i for i, e in enumerate(self.entries) if e.name in names]
        if not rows:
            return
        self._filling = True
        self.samples.clearSelection()
        for r in rows:
            self.samples.item(r).setSelected(True)
        self._filling = False
        if self.entry is None or self.entry.name not in names:
            self.samples.setCurrentRow(rows[0], QtCore.QItemSelectionModel.SelectionFlag.NoUpdate)
        self._draw()

    def _selection_changed(self) -> None:
        if not self._filling and self.sample is not None:
            self._draw()

    def _series(self, entries, xdim, ydim, xr, yr) -> list:
        """Each sample's histogram or contours of the population shown."""
        by = self.colour_by.currentData()
        path = self.node.path
        out = []
        for e in entries:
            sample = self.session.load(e)
            masks = sample.masks()
            if path not in masks:
                continue   # gated its own way: this population is not in it
            rows = masks[path]
            group = e.name if by == BY_SAMPLE else (e.meta.get(by, "") or "—")
            x = self._display(xdim, sample)[rows]
            if ydim == HISTOGRAM:
                cx, h = density.histogram(x, *xr)
                out.append(figures.Series(e.name, group, cx, h))
            else:
                y = self._display(ydim, sample)[rows]
                far = density.outside(x, y, xr, yr)
                pts = np.column_stack([x[far], y[far]])
                if len(pts) > 4000:
                    pts = pts[np.random.default_rng(0).choice(len(pts), 4000, replace=False)]
                out.append(figures.Series(e.name, group, outliers=pts,
                                          levels=density.contours(x, y, xr, yr)))
        return out

    def _draw_series(self, series, ydim, yr) -> tuple[float, float]:
        """Draw them, and say what range of Y they need."""
        names = order.ordered(s.group for s in series)
        colour = dict(zip(names, palette.categorical(len(names))))
        legend = self.plot.plotItem.legend or self.plot.plotItem.addLegend(offset=(-8, 8))
        legend.clear()
        legend.setBrush(theme.legend_brush())
        legend.setPen(theme.legend_pen())
        legend.setLabelTextColor(theme.foreground())
        legend.setVisible(1 < len(names) <= 14)
        shown = set()
        if ydim == HISTOGRAM:
            top = max((float(s.y.max()) for s in series), default=1.0) or 1.0
            step = figures.STAGGER * top if self.offset.isChecked() else 0.0
            for k, s in enumerate(series):
                fill = QtGui.QColor(colour[s.group])
                fill.setAlpha(95 if step else 45)
                curve = pg.PlotCurveItem(s.x, s.y + k * step, fillLevel=k * step,
                                         brush=pg.mkBrush(fill),
                                         pen=pg.mkPen(colour[s.group], width=1.4))
                curve.setZValue(len(series) - k)
                self.plot.addItem(curve)
                if s.group not in shown:
                    shown.add(s.group)
                    legend.addItem(curve, s.group)
            if step:
                self.plot.getAxis("left").setTicks(
                    [[(k * step, s.label) for k, s in enumerate(series)]])
                self.plot.setLabel("left", "")
            else:
                self.plot.getAxis("left").setTicks(None)
                self.plot.setLabel("left", "% of the tallest bin")
            return (0.0, (len(series) - 1) * step + top * 1.06)
        widths = {0.90: 0.8, 0.75: 1.0, 0.50: 1.3, 0.25: 1.6, 0.10: 1.9}
        for s in series:
            first = None
            for share, lines in s.levels or []:
                pen = pg.mkPen(colour[s.group], width=widths.get(round(share, 2), 1.0))
                for line in lines:
                    item = pg.PlotCurveItem(line[:, 0], line[:, 1], pen=pen)
                    self.plot.addItem(item)
                    first = first or item
            if s.outliers is not None and len(s.outliers):
                dots = QtGui.QColor(colour[s.group])
                dots.setAlpha(150)
                self.plot.addItem(pg.ScatterPlotItem(s.outliers[:, 0], s.outliers[:, 1],
                                                     size=2.2, pen=None, pxMode=True,
                                                     brush=pg.mkBrush(dots)))
            if first is not None and s.group not in shown:
                shown.add(s.group)
                legend.addItem(first, s.group)
        return yr

    # --- violins and boxes, every pair, the gates above ---------------------------

    def set_overlay(self, mode: str, pooled: bool = False) -> None:
        """How several samples are shown on a histogram: "overlaid",
        "staggered", "violin" or "box" — the last two one for each sample
        or, `pooled`, one for each value of the "Colour by" column."""
        self.overlay.select(mode, pooled)

    def _distributions(self, entries, xdim, xr, kind) -> views.DistributionContent:
        """The parameter on X in each sample's population shown — or in
        each group's, its samples' events taken together."""
        by = self.colour_by.currentData()
        path = self.node.path
        parts: dict[str, tuple[str, list]] = {}
        for e in entries:
            sample = self.session.load(e)
            masks = sample.masks()
            if path not in masks:
                continue   # gated its own way: this population is not in it
            group = e.name if by == BY_SAMPLE else (e.meta.get(by, "") or "—")
            label = group if self.overlay.pooled() else e.name
            parts.setdefault(label, (group, []))[1].append(
                self._display(xdim, sample)[masks[path]])
        series = [views.distribution(np.concatenate(values), *xr, label, group)
                  for label, (group, values) in sorted(
                      parts.items(), key=lambda part: order.key(part[1][0]))]
        return views.DistributionContent(
            series, kind, self.x_axis.currentText(),
            decade_ticks(self._scale(xdim)) if self._scale(xdim) else None, xr)

    def _label(self, dim: str) -> str:
        i = self.x_axis.findData(dim)
        return self.x_axis.itemText(i) if i >= 0 else dim

    def matrix_default(self) -> list[str]:
        """The stained fluorescence parameters — those the files name a
        marker for; every fluorescence one where none is named, and every
        parameter where nothing is compensated."""
        f, fluo = self.sample.fcs, self.session.fluorescence()
        stained = [self._dim(c) for c, m in zip(f.channels, f.markers) if c in fluo and m]
        return (stained or [self._dim(c) for c in f.channels if c in fluo]
                or list(f.channels))[:views.MAX_PARAMETERS]

    def matrix_content(self, dims: list[str] | None = None) -> views.MatrixContent:
        """Every pair of `dims` for the population shown of this sample."""
        dims = (dims or self.matrix_default())[:views.MAX_PARAMETERS]
        rows = self.sample.mask(self.node.path)
        ranges = [self._range(d) for d in dims]
        ticks = [decade_ticks(self._scale(d)) if self._scale(d) else views.plain_ticks(*r)
                 for d, r in zip(dims, ranges)]
        return views.MatrixContent([self._label(d) for d in dims],
                                   views.pairs([self._display(d)[rows] for d in dims], ranges),
                                   ranges, ticks)

    def backgate_content(self) -> views.BackgateContent:
        """The population shown on each gate above it, and on its own."""
        return views.backgate(self.entry.gates, self.node.path, self.sample.masks(),
                              self._display, self._scale, self._range, self._label)

    def _shown_name(self) -> str:
        name = self.node.name if self.node.gate is not None else "All events"
        return f"{self.entry.name} · {name}"

    def matrix_dialog(self) -> MatrixDialog | None:
        if self.sample is None or self.node is None:
            return None
        parameters = [(self.x_axis.itemText(i), self.x_axis.itemData(i))
                      for i in range(self.x_axis.count())]
        return MatrixDialog(f"Every pair — {self._shown_name()}", parameters,
                            self.matrix_default(), self.matrix_content, self)

    def backgate_dialog(self) -> BackgateDialog | None:
        if self.sample is None or self.node is None:
            return None
        if self.node.gate is None:
            self.sigStatus.emit("Choose a gated population: all events have no gate above them")
            return None
        return BackgateDialog(f"Backgate — {self._shown_name()}", self.backgate_content(), self)

    def open_matrix(self) -> None:
        dialog = self.matrix_dialog()
        if dialog is not None:
            dialog.exec()

    def open_backgate(self) -> None:
        dialog = self.backgate_dialog()
        if dialog is not None:
            dialog.exec()

    def open_model(self, kind: str, show: bool = True):
        """Proliferation, Cell cycle or Kinetics of the population shown,
        on the parameter across, for the samples selected."""
        from .model_dialogs import dialog_for
        dialog = dialog_for(self, kind)
        if dialog is not None and show:
            dialog.exec()
        return dialog

    def content(self):
        """What is on the plot, for a figure."""
        return self._content

    # --- scales, compensation, export ----------------------------------------------

    def set_scale(self, dim: str, scale, every_marker: bool = False) -> None:
        """Draw this parameter on this scale (None: linear, as stored) —
        or, with `every_marker`, every fluorescence parameter."""
        s = self.session
        s.set_scale(sorted(s.fluorescence()) if every_marker else [dim], scale)
        self._scales_changed()

    def reset_scale(self, dim: str, every_marker: bool = False) -> None:
        s = self.session
        s.reset_scale(sorted(s.fluorescence()) if every_marker else [dim])
        self._scales_changed()

    def _scales_changed(self) -> None:
        if self.sample is not None:
            self._keep_axes = True
            try:
                self._fill_tree(select=self.node.path if self.node is not None else "")
            finally:
                self._keep_axes = False
        self.sigScalesChanged.emit()
        self.sigGatesChanged.emit()   # a polygon follows its axis: counts may have moved

    def _ask_scale(self, which: str) -> None:
        combo = self.x_axis if which == "x" else self.y_axis
        dim = combo.currentData()
        if self.sample is None or dim in (None, HISTOGRAM):
            return
        dialog = ScaleDialog(combo.currentText(), self._scale(dim), self)
        if dialog.exec():
            if dialog.reset_asked:
                self.reset_scale(dim, dialog.every.isChecked())
            else:
                self.set_scale(dim, dialog.scale(), dialog.every.isChecked())

    def set_compensation(self, matrix, entries: list | None = None) -> None:
        """Compensate every sample with this matrix; None goes back to
        the workspace's or the files' own. With `entries`, only those
        samples."""
        row = self.samples.currentRow()
        path = self.node.path if self.node is not None else ""
        self.session.set_matrix(matrix, entries)
        if self.entry is not None:
            self.sample = self.session.load(self.entry)
            self._keep_axes = True
            try:
                self._fill_tree(select=path)
            finally:
                self._keep_axes = False
        self.samples.setCurrentRow(row, QtCore.QItemSelectionModel.SelectionFlag.NoUpdate)
        self.show_review()
        self.sigGatesChanged.emit()

    def show_review(self) -> None:
        """Say that the gates were drawn on other values than those shown
        now, for as long as the session has a reason to."""
        why = self.session.review
        self.review.setText(
            f"These gates were drawn before {why}: each cuts the values as they are now, "
            "and may no longer sit on its population. Look at them, then press Reviewed."
            if why else "")
        self.review_bar.setVisible(bool(why))

    def reviewed(self) -> None:
        self.session.clear_review()
        self.show_review()

    def matrix_source(self) -> str:
        """Which matrix the sample being looked at is compensated with:
        `file`, `workspace`, `edited`, `shared`, or nothing."""
        return self.session.matrix_in_use(self.entry)[1] if self.entry is not None else ""

    def _ask_compensation(self) -> None:
        if self.entry is None or not any(self.session.matrix_in_use(e)[0] is not None
                                         for e in self.entries):
            self.sigStatus.emit("These files carry no compensation matrix")
            return
        labels = dict((c, m) for c, m in zip(self.sample.fcs.channels, self.sample.fcs.markers)) \
            if self.sample is not None else {}
        CompensationDialog(self.session, self.entries, self.entry, self.overlaid(), labels,
                           self.set_compensation, self).exec()

    def export_population(self, folder: str, fmt: str = "fcs", compensated: bool = True,
                          every_sample: bool = False) -> list[str]:
        """The population shown, of the samples selected here — or of all
        that are ticked in Samples."""
        entries = ([e for e in self.session.usable() if e.include] if every_sample
                   else self.overlaid())
        written = self.session.export_population(entries, self.node.path, folder, fmt,
                                                 compensated)
        self.sigStatus.emit(f"Wrote {len(written)} file(s) to {folder}")
        return written

    def _ask_export(self) -> None:
        if self.node is None:
            return
        name = self.node.name if self.node.gate is not None else "All events"
        dialog = ExportDialog(name, len(self.overlaid()), self.session.folder, self)
        if dialog.exec() and dialog.folder.text():
            self.export_population(dialog.folder.text(), dialog.format.currentData(),
                                   dialog.compensated.isChecked(), dialog.all_samples.isChecked())

    # --- quality of the acquisition ------------------------------------------------

    def check_quality(self):
        """What `qc.check` finds in the sample being looked at."""
        return self.session.quality(self.entry)

    def apply_quality(self, found: list) -> None:
        """Remove from each sample its own bad stretches of time —
        (entry, intervals) pairs — as a first population every gate of
        every sample is then under."""
        self.session.apply_quality(found)
        self._quality_changed()

    def remove_quality(self) -> None:
        """Take that population out of every sample's tree again."""
        self.session.remove_quality()
        self._quality_changed()

    def _quality_changed(self) -> None:
        self._keep_axes = True
        try:
            self._fill_tree(select=self.node.path if self.node is not None else "")
        finally:
            self._keep_axes = False
        self.sigGatesChanged.emit()

    def open_quality(self):
        """The dialog, not yet shown; None before a sample is chosen."""
        if self.sample is None:
            return None
        from .quality_dialog import QualityDialog
        return QualityDialog(self, self)

    def _ask_quality(self) -> None:
        dialog = self.open_quality()
        if dialog is not None:
            dialog.exec()

    def make_figure(self) -> None:
        from .figure_dialog import open_for
        if open_for(self, self, "gates") is None:
            self.sigStatus.emit("Nothing on the plot yet")

    def _fill_axes(self) -> None:
        keep = (self.x_axis.currentData(), self.y_axis.currentData())
        self._filling = True
        f = self.sample.fcs
        for combo in (self.x_axis, self.y_axis):
            combo.clear()
            for c, m in zip(f.channels, f.markers):
                combo.addItem(f"{m}" if m else c, self._dim(c))
        self.y_axis.insertSeparator(self.y_axis.count())
        self.y_axis.addItem("Histogram", HISTOGRAM)
        for combo, value, default in ((self.x_axis, keep[0], 0), (self.y_axis, keep[1], 2)):
            i = combo.findData(value) if value else -1
            combo.setCurrentIndex(i if i >= 0 else min(default, combo.count() - 1))
        self._filling = False

    def set_axes(self, x: str, y: str) -> None:
        self._filling = True
        for combo, value in ((self.x_axis, x), (self.y_axis, y)):
            i = combo.findData(value)
            if i >= 0:
                combo.setCurrentIndex(i)
        self._filling = False
        self._draw()

    def _axes_changed(self) -> None:
        if not self._filling:
            self._draw()

    # --- the tree ---------------------------------------------------------------

    def _fill_tree(self, select: str = "") -> None:
        self._filling = True
        self.tree.clear()
        masks = self.sample.masks()
        total = len(self.sample)
        found = None

        def add(node: GateNode, parent_item, parent_count: int):
            nonlocal found
            n = int(masks[node.path].sum())
            pct = f"{100 * n / parent_count:.1f}" if node.gate is not None and parent_count else ""
            item = QtWidgets.QTreeWidgetItem(
                [node.name if node.gate is not None else "All events", f"{n:,d}", pct])
            item.setData(0, QtCore.Qt.ItemDataRole.UserRole, node.path)
            for col in (1, 2):
                item.setTextAlignment(col, QtCore.Qt.AlignmentFlag.AlignRight
                                      | QtCore.Qt.AlignmentFlag.AlignVCenter)
            if node.note:
                # read from a workspace by a rule not checked against FlowJo:
                # set apart, with FlowJo's own count to compare
                font = item.font(0)
                font.setItalic(True)
                item.setFont(0, font)
                theirs = (f" FlowJo counted {node.reference_count:,d}; here {n:,d}."
                          if node.reference_count is not None else "")
                for col in range(3):
                    item.setToolTip(col, node.note + theirs)
            (parent_item.addChild if parent_item is not None
             else self.tree.addTopLevelItem)(item)
            if node.path == select:
                found = item
            for child in node.children:
                add(child, item, n)

        add(self.entry.gates, None, total)
        self.tree.expandAll()
        for col in (1, 2):
            self.tree.resizeColumnToContents(col)
        self._filling = False
        self.tree.setCurrentItem(found or self.tree.topLevelItem(0))
        self.reset_button.setEnabled(any(e.edited for e in self.entries))

    def _node_chosen(self, item, _previous=None) -> None:
        if self._filling or item is None or self.entry is None:
            return
        path = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
        self.node = self.entry.gates.find(path) if path else self.entry.gates
        # arriving at a population means looking at it the way its gates were
        # drawn; staying on it — after drawing one more — means staying put
        for child in ([] if self._keep_axes else self.node.children):
            dims = child.gate.dims
            if dims and self.x_axis.findData(dims[0]) >= 0:
                self.set_axes(dims[0], dims[1] if len(dims) > 1 else HISTOGRAM)
                break
        else:
            self._draw()
        self.rename_button.setEnabled(self.node.gate is not None)
        self.delete_button.setEnabled(self.node.gate is not None)

    def select(self, path: str) -> None:
        it = QtWidgets.QTreeWidgetItemIterator(self.tree)
        while it.value() is not None:
            if it.value().data(0, QtCore.Qt.ItemDataRole.UserRole) == path:
                self.tree.setCurrentItem(it.value())
                return
            it += 1

    # --- drawing the events ------------------------------------------------------

    def _draw(self) -> None:
        if self.sample is None or self.node is None:
            return
        self.stack.setCurrentIndex(1)
        self.plot.clear()
        self.plot.setBackground(theme.background())
        self._overlays = []
        theme.style_axes(self.plot.plotItem)
        xdim, ydim = self.x_axis.currentData(), self.y_axis.currentData()
        mask = self.sample.mask(self.node.path)
        rows = np.flatnonzero(mask)
        shown = rows
        if rows.size > MAX_POINTS:
            shown = np.sort(np.random.default_rng(0).choice(rows, MAX_POINTS, replace=False))
        x = self._display(xdim)
        xr = self._range(xdim)
        self.plot.setLabel("bottom", self.x_axis.currentText())
        self._ticks("bottom", xdim)
        if self.plot.plotItem.legend is not None:
            self.plot.plotItem.legend.setVisible(False)
        overlay = self.overlaid()
        as_lines = len(overlay) > 1 or (ydim != HISTOGRAM
                                        and self.kind.currentData() == CONTOUR)
        self.kind.setEnabled(ydim != HISTOGRAM and len(overlay) <= 1)
        self.colour_by.setEnabled(len(overlay) > 1)
        self.overlay.enable(len(overlay) > 1 and ydim == HISTOGRAM)
        self._content = None
        shape = self.overlay.mode() if len(overlay) > 1 and ydim == HISTOGRAM else None
        if shape in (views.VIOLIN, views.BOX):
            self._content = self._distributions(overlay, xdim, xr, shape)
            yr = draw_distributions(self.plot, self._content)
        elif as_lines:
            yr = (0.0, 1.0) if ydim == HISTOGRAM else self._range(ydim)
            if ydim != HISTOGRAM:
                self.plot.setLabel("left", self.y_axis.currentText())
                self._ticks("left", ydim)
            series = self._series(overlay, xdim, ydim, xr, yr)
            yr = self._draw_series(series, ydim, yr)
            hist = ydim == HISTOGRAM
            self._content = figures.OverlayContent(
                series, self.x_axis.currentText(),
                "% of the tallest bin" if hist else self.y_axis.currentText(),
                decade_ticks(self._scale(xdim)) if self._scale(xdim) else None,
                None if hist or not self._scale(ydim) else decade_ticks(self._scale(ydim)),
                xr, None if hist else yr, offset=hist and self.offset.isChecked())
        elif ydim == HISTOGRAM:
            counts, edges = np.histogram(x[rows], bins=256, range=xr)
            fill = QtGui.QColor(theme.accent())
            fill.setAlpha(90)
            self.plot.addItem(pg.PlotCurveItem(edges, counts, stepMode="center", fillLevel=0,
                                               brush=pg.mkBrush(fill),
                                               pen=pg.mkPen(theme.accent(), width=1.2)))
            self.plot.setLabel("left", "Events")
            self.plot.getAxis("left").setTicks(None)
            yr = (0.0, float(counts.max()) * 1.08 + 1)
            self._content = figures.OverlayContent(
                [figures.Series(self.entry.name, self.entry.name,
                                (edges[:-1] + edges[1:]) / 2, counts.astype(float))],
                self.x_axis.currentText(), "Events",
                decade_ticks(self._scale(xdim)) if self._scale(xdim) else None, xrange=xr)
        else:
            y = self._display(ydim)
            yr = self._range(ydim)
            self.plot.setLabel("left", self.y_axis.currentText())
            self._ticks("left", ydim)
            if shown.size:
                # pseudocolour: each event takes the colour of how crowded its
                # neighbourhood is, on a log scale
                H, ex, ey = np.histogram2d(x[shown], y[shown], bins=160, range=[xr, yr])
                ix = np.clip(np.searchsorted(ex, x[shown], side="right") - 1, 0, 159)
                iy = np.clip(np.searchsorted(ey, y[shown], side="right") - 1, 0, 159)
                d = np.log1p(H[ix, iy])
                level = (d / d.max() * 63).astype(int) if d.max() > 0 else np.zeros(d.size, int)
                lut = figures.lookup_table("turbo", steps=64)
                brushes = [pg.mkBrush(int(r), int(g), int(b)) for r, g, b, _ in lut]
                order = np.argsort(level, kind="stable")   # the dense on top
                self.plot.addItem(pg.ScatterPlotItem(
                    x[shown][order], y[shown][order], size=2.5, pen=None, pxMode=True,
                    brush=[brushes[i] for i in level[order]]))
                self._content = figures.MapContent(
                    x[shown][order], y[shown][order], self.x_axis.currentText(),
                    self.y_axis.currentText(),
                    color=figures.Values(d[order], "Density (log)", float(d.min()), float(d.max())),
                    xticks=decade_ticks(self._scale(xdim)) if self._scale(xdim) else None,
                    yticks=decade_ticks(self._scale(ydim)) if self._scale(ydim) else None)
        self.plot.setXRange(*xr, padding=0.01)
        self.plot.setYRange(*yr, padding=0.01)
        self._crossings = []
        self._draw_gates(xdim, ydim, xr, yr)
        self._draw_quadrants(xdim, ydim, xr, yr, self.sample.masks(),
                             max(int(self.sample.mask(self.node.path).sum()), 1))
        parent = self.node.parent
        of = (f", {100 * rows.size / max(int(self.sample.mask(parent.path).sum()), 1):.1f}% "
              f"of {parent.name if parent.gate is not None else 'all events'}"
              if parent is not None else "")
        more = f" — {shown.size:,d} drawn" if shown.size < rows.size and not as_lines else ""
        name = self.node.name if self.node.gate is not None else "All events"
        together = (f" · {len(overlay)} samples overlaid; gates are {self.entry.name}'s"
                    if len(overlay) > 1 else "")
        if isinstance(self._content, views.DistributionContent):
            together += " · " + self._content.note()
        if getattr(self.sample, "normalised", False):
            together += " · batch-normalised values"
        matrix = {"file": "its file's matrix", "workspace": "the workspace's matrix",
                  "edited": "a matrix edited for this sample",
                  "shared": "the edited matrix every sample shares"}.get(self.matrix_source())
        if matrix:
            together += f" · compensated with {matrix}"
        self.caption.setText(
            f"{self.entry.name} · {name}: {rows.size:,d} events{of}{more}{together}")

    def _range(self, dim: str) -> tuple[float, float]:
        scale = self._scale(dim)
        if scale is not None:
            return display_range(scale)
        base = dim
        f = self.sample.fcs
        top = f.ranges[f.index(base)] if base in f.channels else float(self._display(dim).max())
        return (0.0, float(top))

    def _ticks(self, side: str, dim: str) -> None:
        scale = self._scale(dim)
        self.plot.getAxis(side).setTicks([decade_ticks(scale)] if scale is not None else None)

    # --- the gates on the plot -----------------------------------------------------

    def _draw_gates(self, xdim, ydim, xr, yr) -> None:
        masks = self.sample.masks()
        here = max(int(masks[self.node.path].sum()), 1)
        pen = pg.mkPen(theme.foreground(), width=1.6)
        hover = pg.mkPen(theme.accent(), width=2.2)
        for child in self.node.children:
            g = child.gate
            item = None
            anchor = None
            if isinstance(g, PolygonGate) and ydim != HISTOGRAM and set(g.dims) == {xdim, ydim}:
                pts = self._to_display(g.vertices, g.dims)
                if g.dims[0] != xdim:
                    pts = pts[:, ::-1]
                item = pg.PolyLineROI(pts.tolist(), closed=True, pen=pen, hoverPen=hover)
                anchor = (pts[:, 0].min(), pts[:, 1].max())
            elif (isinstance(g, EllipseGate) and ydim != HISTOGRAM and g.dims == [xdim, ydim]
                    and not all(same(a, self._scale(d))
                                for a, d in zip(g.scales or [None, None], g.dims))):
                # drawn on another scale than the axes have now: shown where it
                # is, as the curve it has become, but not for dragging
                t = np.linspace(0, 2 * np.pi, 121)
                c, sn = np.cos(g.angle), np.sin(g.angle)
                u, v = g.radii[0] * np.cos(t), g.radii[1] * np.sin(t)
                own = np.column_stack([g.center[0] + u * c - v * sn,
                                       g.center[1] + u * sn + v * c])
                for k, scale in enumerate(g.scales or [None, None]):
                    if scale is not None:
                        own[:, k] = scale.inverse(own[:, k])
                shown = self._to_display(own, g.dims)
                outline = pg.PlotCurveItem(shown[:, 0], shown[:, 1], pen=pen)
                self.plot.addItem(outline)
                n = int(masks[child.path].sum())
                label = pg.TextItem(f"{child.name}  {100 * n / here:.1f}%", anchor=(0, 1),
                                    color=theme.foreground(),
                                    fill=pg.mkBrush(theme.legend_brush().color()))
                label.setPos(float(shown[:, 0].min()), float(shown[:, 1].max()))
                self.plot.addItem(label)
                self._overlays.append((child, outline, label))
                continue
            elif (isinstance(g, EllipseGate) and ydim != HISTOGRAM
                    and g.dims == [xdim, ydim]):
                # pyqtgraph turns an ellipse about its corner, not its centre
                w, h = 2 * g.radii[0], 2 * g.radii[1]
                c, sn = np.cos(g.angle), np.sin(g.angle)
                corner = [g.center[0] - (w / 2 * c - h / 2 * sn),
                          g.center[1] - (w / 2 * sn + h / 2 * c)]
                item = pg.EllipseROI(corner, [w, h], angle=float(np.degrees(g.angle)),
                                     pen=pen, hoverPen=hover)
                anchor = (g.center[0] - g.radii[0], g.center[1] + g.radii[1])
            elif isinstance(g, QuadrantGate):
                continue   # the four corners of a crossing are drawn once, below
            elif isinstance(g, TimeGate):
                # what quality control removed, to see and not to drag
                shade = QtGui.QColor(theme.danger())
                shade.setAlpha(46)
                for lo, hi in (g.intervals if xdim == g.dims[0] else []):
                    self.plot.addItem(pg.LinearRegionItem(
                        (lo, hi), movable=False, brush=pg.mkBrush(shade), pen=pg.mkPen(None)))
                continue
            elif isinstance(g, RectangleGate):
                bounds = dict(zip(g.dims, zip(g.mins, g.maxs)))
                if xdim in bounds and (ydim in bounds or len(g.dims) == 1):
                    x0, x1 = self._bounds_display(xdim, bounds[xdim], xr)
                    if ydim in bounds:
                        y0, y1 = self._bounds_display(ydim, bounds[ydim], yr)
                        item = pg.RectROI([x0, y0], [x1 - x0, y1 - y0], pen=pen, hoverPen=hover)
                        item.addScaleHandle([0, 0], [1, 1])
                        anchor = (x0, y1)
                    else:
                        item = pg.LinearRegionItem((x0, x1), pen=pen, hoverPen=hover,
                                                   brush=pg.mkBrush(None))
                        anchor = (x0, yr[1])
            if item is None:
                continue
            self.plot.addItem(item)
            n = int(masks[child.path].sum())
            label = pg.TextItem(f"{child.name}  {100 * n / here:.1f}%", anchor=(0, 1),
                                color=theme.foreground(),
                                fill=pg.mkBrush(theme.legend_brush().color()))
            label.setPos(*anchor)
            self.plot.addItem(label)
            item.sigRegionChangeFinished.connect(
                lambda _=None, c=child, it=item: self._gate_moved(c, it))
            self._overlays.append((child, item, label))

    def _draw_quadrants(self, xdim, ydim, xr, yr, masks, here) -> None:
        """One crossing for each set of four made together, with each
        corner's share written in it."""
        if ydim == HISTOGRAM:
            return
        groups: dict[str, list[GateNode]] = {}
        for child in self.node.children:
            g = child.gate
            if isinstance(g, QuadrantGate) and set(g.dims) == {xdim, ydim}:
                groups.setdefault(g.group or child.name, []).append(child)
        pen = pg.mkPen(theme.foreground(), width=1.6)
        hover = pg.mkPen(theme.accent(), width=2.2)
        for members in groups.values():
            g = members[0].gate
            flipped = g.dims[0] != xdim
            px, py = self._to_display([[g.x, g.y] if not flipped else [g.y, g.x]],
                                      [xdim, ydim])[0]
            vertical = pg.InfiniteLine(px, angle=90, movable=True, pen=pen, hoverPen=hover)
            horizontal = pg.InfiniteLine(py, angle=0, movable=True, pen=pen, hoverPen=hover)
            labels = []
            for child in members:
                high_x, high_y = QUADRANTS[child.gate.quadrant]
                if flipped:
                    high_x, high_y = high_y, high_x
                n = int(masks[child.path].sum())
                label = pg.TextItem(f"{child.name}  {100 * n / here:.1f}%",
                                    anchor=(1 if high_x else 0, 0 if high_y else 1),
                                    color=theme.foreground(),
                                    fill=pg.mkBrush(theme.legend_brush().color()))
                label.setPos(xr[1] if high_x else xr[0], yr[1] if high_y else yr[0])
                self.plot.addItem(label)
                labels.append((child, label))
            for line in (vertical, horizontal):
                self.plot.addItem(line)
                line.sigPositionChangeFinished.connect(
                    lambda _=None, m=members, v=vertical, h=horizontal:
                    self._crossing_moved(m, v.value(), h.value()))
            self._crossings.append((members, vertical, horizontal, labels))

    def _crossing_moved(self, members, px: float, py: float) -> None:
        xdim, ydim = self.x_axis.currentData(), self.y_axis.currentData()
        sx, sy = self._to_stored([[px, py]], [xdim, ydim])[0]
        for child in members:
            g = child.gate
            g.x, g.y = (float(sx), float(sy)) if g.dims[0] == xdim else (float(sy), float(sx))
        self._commit(redraw=False)

    def _to_display(self, vertices, dims) -> np.ndarray:
        out = np.array(vertices, dtype=np.float64)
        for k, dim in enumerate(dims):
            scale = self._scale(dim)
            if scale is not None:
                out[:, k] = scale(out[:, k])
        return out

    def _to_stored(self, points, dims) -> np.ndarray:
        out = np.array(points, dtype=np.float64)
        for k, dim in enumerate(dims):
            scale = self._scale(dim)
            if scale is not None:
                out[:, k] = scale.inverse(out[:, k])
        return out

    def _bounds_display(self, dim, bounds, view) -> tuple[float, float]:
        scale = self._scale(dim)
        lo, hi = bounds
        f = (lambda v: float(scale(v))) if scale is not None else float
        return (view[0] if lo is None else f(lo), view[1] if hi is None else f(hi))

    def _gate_moved(self, child: GateNode, item) -> None:
        """A handle was let go: the gate is where the outline now is."""
        g = child.gate
        xdim, ydim = self.x_axis.currentData(), self.y_axis.currentData()
        if isinstance(item, pg.EllipseROI):
            w, h = item.size().x(), item.size().y()
            centre = item.mapToParent(QtCore.QPointF(w / 2, h / 2))
            g.center = [float(centre.x()), float(centre.y())]
            g.radii = [abs(float(w)) / 2, abs(float(h)) / 2]
            g.angle = float(np.radians(item.angle()))
        elif isinstance(item, pg.PolyLineROI):
            pts = np.array([[p.x(), p.y()] for _, p in item.getLocalHandlePositions()])
            pts += np.array([item.pos().x(), item.pos().y()])
            dims = [xdim, ydim]
            if g.dims[0] != xdim:
                pts, dims = pts[:, ::-1], [ydim, xdim]
            g.vertices = self._to_stored(pts, dims)
        elif isinstance(item, pg.RectROI):
            x0, y0 = item.pos().x(), item.pos().y()
            x1, y1 = x0 + item.size().x(), y0 + item.size().y()
            lo = self._to_stored([[x0, y0]], [xdim, ydim])[0]
            hi = self._to_stored([[x1, y1]], [xdim, ydim])[0]
            for k, dim in enumerate((xdim, ydim)):
                i = g.dims.index(dim)
                g.mins[i], g.maxs[i] = float(lo[k]), float(hi[k])
        else:
            x0, x1 = item.getRegion()
            lo, hi = self._to_stored([[x0], [x1]], [xdim])[:, 0]
            g.mins[0], g.maxs[0] = float(lo), float(hi)
        self._commit(redraw=False)

    # --- making, renaming, deleting ---------------------------------------------------

    def _tool(self, mode: str, on: bool) -> None:
        hints = {"rect": "Drag across the events to gate",
                 "poly": "Click each corner; double-click to close. Esc cancels.",
                 "ellipse": "Drag the box the ellipse fits in",
                 "quad": "Click where the two lines should cross"}
        if on:
            for other, button in self._tools.items():
                if other != mode:
                    button.setChecked(False)
            if mode != "rect" and self.y_axis.currentData() == HISTOGRAM:
                self.sigStatus.emit("That gate needs two parameters: choose a Y axis")
                self._tools[mode].setChecked(False)
                return
            self.box.set_mode(mode)
            self.sigStatus.emit(hints[mode])
        elif self.box.mode == mode:
            self.box.set_mode(None)
            self.sigStatus.emit("")

    def _stop_drawing(self) -> None:
        for button in self._tools.values():
            button.setChecked(False)
        self.box.set_mode(None)

    def _ask_name(self) -> str | None:
        default = self._free_name("Gate")
        name, ok = QtWidgets.QInputDialog.getText(self, "New gate", "Name", text=default)
        return name.strip() if ok and name.strip() else None

    def _free_name(self, wanted: str) -> str:
        taken = {c.name for c in self.node.children}
        name, k = wanted, 2
        while name in taken:
            name, k = f"{wanted} {k}", k + 1
        return name

    def _rectangle_drawn(self, x0, y0, x1, y1) -> None:
        self._stop_drawing()
        name = self._ask_name()
        if name:
            self.add_rectangle(x0, y0, x1, y1, name)

    def _ellipse_drawn(self, x0, y0, x1, y1) -> None:
        self._stop_drawing()
        name = self._ask_name()
        if name:
            self.add_ellipse((x0 + x1) / 2, (y0 + y1) / 2, (x1 - x0) / 2, (y1 - y0) / 2, name)

    def _quadrant_drawn(self, x, y) -> None:
        self._stop_drawing()
        self.add_quadrants(x, y)

    def _polygon_drawn(self, points) -> None:
        self._stop_drawing()
        name = self._ask_name()
        if name:
            self.add_polygon(points, name)

    def add_rectangle(self, x0, y0, x1, y1, name: str) -> GateNode:
        """In the plot's own coordinates. On a histogram, a range of X."""
        xdim, ydim = self.x_axis.currentData(), self.y_axis.currentData()
        if ydim == HISTOGRAM:
            lo, hi = self._to_stored([[x0], [x1]], [xdim])[:, 0]
            gate = RectangleGate([xdim], [float(lo)], [float(hi)])
        else:
            lo = self._to_stored([[x0, y0]], [xdim, ydim])[0]
            hi = self._to_stored([[x1, y1]], [xdim, ydim])[0]
            gate = RectangleGate([xdim, ydim], [float(lo[0]), float(lo[1])],
                                 [float(hi[0]), float(hi[1])])
        return self._add(gate, name)

    def add_polygon(self, points, name: str) -> GateNode:
        xdim, ydim = self.x_axis.currentData(), self.y_axis.currentData()
        scales = [self._scale(xdim), self._scale(ydim)]
        gate = PolygonGate([xdim, ydim], self._to_stored(points, [xdim, ydim]),
                           scales=scales if any(s is not None for s in scales) else None)
        return self._add(gate, name)

    def add_ellipse(self, cx, cy, rx, ry, name: str, angle: float = 0.0) -> GateNode:
        """Centre and half-axes in the plot's own coordinates; `angle` in
        radians, counter-clockwise."""
        xdim, ydim = self.x_axis.currentData(), self.y_axis.currentData()
        scales = [self._scale(xdim), self._scale(ydim)]
        gate = EllipseGate([xdim, ydim], [float(cx), float(cy)],
                           [abs(float(rx)), abs(float(ry))], float(angle),
                           scales=scales if any(s is not None for s in scales) else None)
        return self._add(gate, name)

    def add_quadrants(self, x, y) -> list[GateNode]:
        """Four populations around the point, named by the side of each
        line they are on: `CD4+ CD8-`."""
        xdim, ydim = self.x_axis.currentData(), self.y_axis.currentData()
        sx, sy = self._to_stored([[x, y]], [xdim, ydim])[0]
        short = [self.x_axis.currentText().split()[0], self.y_axis.currentText().split()[0]]
        group = f"{xdim}|{ydim}|{len(self.node.children)}"
        made = []
        for quadrant, (high_x, high_y) in QUADRANTS.items():
            name = f"{short[0]}{'+' if high_x else '-'} {short[1]}{'+' if high_y else '-'}"
            made.append(self.node.add(GateNode(
                self._free_name(name),
                QuadrantGate([xdim, ydim], float(sx), float(sy), quadrant, group))))
        self._commit()
        return made

    def add_combination(self, op: str, names: list[str], name: str) -> GateNode:
        """`op` is AND, OR or NOT; `names` are populations of the one shown."""
        return self._add(BooleanGate(op, list(names)), name)

    def _ask_combination(self) -> None:
        if self.node is None:
            return
        sisters = [c.name for c in self.node.children]
        if not sisters:
            self.sigStatus.emit("There is nothing to combine yet: gate this population first")
            return
        dialog = CombineDialog(sisters, self)
        if dialog.exec() and dialog.chosen():
            self.add_combination(dialog.op.currentData(), dialog.chosen(),
                                 dialog.name.text().strip() or dialog.suggestion())

    def _add(self, gate, name: str) -> GateNode:
        node = self.node.add(GateNode(self._free_name(name), gate))
        self._commit()
        return node

    def rename_selected(self) -> None:
        if self.node is None or self.node.gate is None:
            return
        name, ok = QtWidgets.QInputDialog.getText(self, "Rename", "Name", text=self.node.name)
        if ok and name.strip():
            self.rename(name.strip())

    def rename(self, name: str) -> None:
        old, parent = self.node.name, self.node.parent
        self.node.name = name
        # a combination names what it is made of: it follows the new name
        for sister in (parent.children if parent is not None else []):
            if depends_on(sister, old):
                sister.gate.operands = [name if o == old else o for o in sister.gate.operands]
        self._commit()

    def delete_selected(self) -> None:
        """The gate being looked at, and everything gated inside it."""
        node = self.node
        if node is None or node.gate is None:
            return
        if node is time_node(self.entry.gates):
            # everything is under it: deleting it would delete every gate
            self.remove_quality()
            return
        parent = node.parent
        parent.children.remove(node)
        # and what was made of it, which would otherwise be left empty
        for sister in list(parent.children):
            if depends_on(sister, node.name):
                parent.children.remove(sister)
        self.node = parent
        self._commit()

    def reset(self) -> None:
        self.session.discard_gates()
        row = self.samples.currentRow()
        self.refresh()
        self.samples.setCurrentRow(max(row, 0))
        self.sigGatesChanged.emit()

    def _commit(self, redraw: bool = True) -> None:
        """The tree of this sample changed: tell the session, count again."""
        self.session.gates_changed(self.entry, everyone=self.everyone.isChecked())
        self.sample = self.session.load(self.entry)
        path = self.node.path if self.node is not None else ""
        if redraw:
            self._keep_axes = True
            try:
                self._fill_tree(select=path)
            finally:
                self._keep_axes = False
        else:
            # a gate was dragged: its outline is already where it should be,
            # and rebuilding the plot under the mouse would drop the handle
            self._refresh_counts()
        self.sigGatesChanged.emit()

    def _refresh_counts(self) -> None:
        masks = self.sample.masks()
        it = QtWidgets.QTreeWidgetItemIterator(self.tree)
        while it.value() is not None:
            item = it.value()
            path = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
            n = int(masks[path].sum())
            item.setText(1, f"{n:,d}")
            parent = item.parent()
            if parent is not None:
                pn = int(masks[parent.data(0, QtCore.Qt.ItemDataRole.UserRole)].sum())
                item.setText(2, f"{100 * n / pn:.1f}" if pn else "")
            it += 1
        here = max(int(masks[self.node.path].sum()), 1)
        for child, _item, label in self._overlays:
            label.setText(f"{child.name}  {100 * int(masks[child.path].sum()) / here:.1f}%")
        for _members, _v, _h, labels in getattr(self, "_crossings", []):
            for child, label in labels:
                label.setText(f"{child.name}  {100 * int(masks[child.path].sum()) / here:.1f}%")
        self.reset_button.setEnabled(True)


class CombineDialog(QtWidgets.QDialog):
    """A population made of others: which, and how."""

    WORDS = {AND: "and", OR: "or", NOT: "not"}

    def __init__(self, sisters: list[str], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Combine populations")
        self.op = QtWidgets.QComboBox()
        self.op.addItem("In all of them (AND)", AND)
        self.op.addItem("In any of them (OR)", OR)
        self.op.addItem("In none of them (NOT)", NOT)
        self.list = QtWidgets.QListWidget()
        for name in sisters:
            item = QtWidgets.QListWidgetItem(name)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.list.addItem(item)
        self.name = QtWidgets.QLineEdit()
        self.box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                              | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.box.accepted.connect(self.accept)
        self.box.rejected.connect(self.reject)
        form = QtWidgets.QFormLayout(self)
        form.addRow("Events", self.op)
        form.addRow(self.list)
        form.addRow("Name", self.name)
        form.addRow(self.box)
        self.op.currentIndexChanged.connect(self._show)
        self.list.itemChanged.connect(self._show)
        self._show()

    def chosen(self) -> list[str]:
        return [self.list.item(i).text() for i in range(self.list.count())
                if self.list.item(i).checkState() == QtCore.Qt.CheckState.Checked]

    def suggestion(self) -> str:
        names, op = self.chosen(), self.op.currentData()
        if op == NOT:
            return "not " + " or ".join(names) if names else ""
        return f" {self.WORDS[op]} ".join(names)

    def _show(self, *_) -> None:
        self.name.setPlaceholderText(self.suggestion())
        self.box.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setEnabled(
            bool(self.chosen()))


class ScaleDialog(QtWidgets.QDialog):
    """The scale of one parameter's axis."""

    def __init__(self, label: str, current, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Scale of {label}")
        self.reset_asked = False
        self.kind = QtWidgets.QComboBox()
        for text, key in (("Biexponential", "biex"), ("Logicle", "logicle"),
                          ("Arcsinh", "arcsinh"), ("Logarithmic", "log"),
                          ("Linear", "linear")):
            self.kind.addItem(text, key)
        b = current if isinstance(current, Biex) else Biex()
        self.width = QtWidgets.QDoubleSpinBox()
        self.width.setRange(1, 5000)
        self.width.setDecimals(0)
        self.width.setValue(abs(b.width))
        self.width.setToolTip("How wide the nearly linear part around zero is. Larger "
                              "pulls a smeared negative population together.")
        self.neg = QtWidgets.QDoubleSpinBox()
        self.neg.setRange(0, 3)
        self.neg.setSingleStep(0.25)
        self.neg.setValue(b.neg)
        self.neg.setToolTip("More room below zero, in decades")
        self.pos = QtWidgets.QDoubleSpinBox()
        self.pos.setRange(2, 7)
        self.pos.setSingleStep(0.25)
        self.pos.setDecimals(4)
        self.pos.setValue(b.pos)
        # FlowJo's own number of decades has nine decimals; shown with four,
        # it must not come back rounded just because the dialog was opened
        self._pos = b.pos
        self.cofactor = QtWidgets.QDoubleSpinBox()
        self.cofactor.setRange(1, 10000)
        self.cofactor.setDecimals(0)
        self.cofactor.setValue(current.cofactor if isinstance(current, Arcsinh) else 150)
        self.cofactor.setToolTip("150 is usual for fluorescence, 5 for mass cytometry")
        self.lowest = QtWidgets.QDoubleSpinBox()
        self.lowest.setRange(0.01, 10000)
        self.lowest.setDecimals(2)
        self.lowest.setValue(current.min if isinstance(current, Log) else 10)
        self.lowest.setToolTip("Zero and everything below it is drawn at this edge")
        g = current if isinstance(current, Logicle) else Logicle()
        self.logicle_w = QtWidgets.QDoubleSpinBox()
        self.logicle_w.setSingleStep(0.1)
        self.logicle_w.setToolTip("W: the decades given to the nearly linear part around "
                                  "zero. Larger pulls a smeared negative population together.")
        self.logicle_m = QtWidgets.QDoubleSpinBox()
        self.logicle_m.setRange(1, 8)
        self.logicle_m.setSingleStep(0.25)
        self.logicle_m.setToolTip("M: the decades the whole scale spans; 4.5 is the usual")
        self.logicle_a = QtWidgets.QDoubleSpinBox()
        self.logicle_a.setSingleStep(0.25)
        self.logicle_a.setToolTip("A: more room below zero, in decades")
        self.logicle_t = QtWidgets.QDoubleSpinBox()
        self.logicle_t.setRange(1, 1e9)
        self.logicle_t.setDecimals(0)
        self.logicle_t.setGroupSeparatorShown(True)
        self.logicle_t.setToolTip("T: the value at the top of the scale; 262,144 for an "
                                  "18-bit cytometer")
        self.logicle_m.setValue(g.M)
        self._logicle_limits()
        self.logicle_w.setValue(g.W)
        self._logicle_limits()
        self.logicle_a.setValue(g.A)
        self.logicle_t.setValue(g.T)
        self.logicle_m.valueChanged.connect(self._logicle_limits)
        self.logicle_w.valueChanged.connect(self._logicle_limits)
        self.every = QtWidgets.QCheckBox("Every fluorescence parameter, not just this one")
        self.note = QtWidgets.QLabel(
            "Polygons drawn on this parameter are tested on the new scale. "
            "The algorithms in Explore are fed each marker on its own scale.")
        self.note.setWordWrap(True)
        self.note.setProperty("role", "caption")
        self.form = QtWidgets.QFormLayout(self)
        self.form.addRow("Scale", self.kind)
        self.rows = {"biex": [("Linear width", self.width), ("Negative decades", self.neg),
                              ("Positive decades", self.pos)],
                     "logicle": [("Linear width (W), decades", self.logicle_w),
                                 ("Decades (M)", self.logicle_m),
                                 ("Negative decades (A)", self.logicle_a),
                                 ("Top of scale (T)", self.logicle_t)],
                     "arcsinh": [("Cofactor", self.cofactor)],
                     "log": [("Lowest value", self.lowest)], "linear": []}
        self._labels = {}
        for rows in self.rows.values():
            for text, widget in rows:
                label = QtWidgets.QLabel(text)
                self._labels[widget] = label
                self.form.addRow(label, widget)
        self.form.addRow(self.every)
        self.form.addRow(self.note)
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                         | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        usual = box.addButton("The usual one", QtWidgets.QDialogButtonBox.ButtonRole.ResetRole)
        usual.clicked.connect(self._reset)
        box.accepted.connect(self.accept)
        box.rejected.connect(self.reject)
        self.form.addRow(box)
        key = ("biex" if isinstance(current, Biex) else "logicle" if isinstance(current, Logicle)
               else "arcsinh" if isinstance(current, Arcsinh) else "log" if isinstance(current, Log) else "linear")
        self.kind.setCurrentIndex(self.kind.findData(key))
        self.kind.currentIndexChanged.connect(self._show)
        self._show()

    def _show(self, *_) -> None:
        for key, rows in self.rows.items():
            for _text, widget in rows:
                on = key == self.kind.currentData()
                widget.setVisible(on)
                self._labels[widget].setVisible(on)

    def _logicle_limits(self, *_) -> None:
        """What a Logicle can be: W up to M / 2, A from −W to M − 2W."""
        m, w = self.logicle_m.value(), self.logicle_w.value()
        self.logicle_w.setRange(0, m / 2)
        w = min(w, m / 2)
        self.logicle_a.setRange(-w, m - 2 * w)

    def _reset(self) -> None:
        self.reset_asked = True
        self.accept()

    def scale(self):
        key = self.kind.currentData()
        if key == "logicle":
            return Logicle(self.logicle_t.value(), self.logicle_w.value(),
                           self.logicle_m.value(), self.logicle_a.value())
        if key == "biex":
            pos = self._pos if abs(self.pos.value() - self._pos) < 5e-5 else self.pos.value()
            return Biex(256, 262144.0, self.neg.value(), -self.width.value(), pos)
        if key == "arcsinh":
            return Arcsinh(self.cofactor.value())
        if key == "log":
            return Log(self.lowest.value())
        return None


class CompensationDialog(QtWidgets.QDialog):
    """Which spillover matrix each sample is compensated with, and the
    matrix of one of them as percentages: row, the fluorochrome; column,
    the detector it spills into. An edited matrix is given to the samples
    selected in the list, or to every sample; `apply(matrix, entries)` is
    what does it — `GatesTab.set_compensation`."""

    SOURCES = {"file": "the matrix in its own FCS file",
               "workspace": "the matrix the FlowJo workspace gives it",
               "edited": "a matrix edited here for this sample",
               "shared": "a matrix edited here, shared by every sample",
               "": "no matrix: the file carries none"}

    def __init__(self, session, entries: list, entry, selected: list, labels: dict,
                 apply, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Compensation")
        self.session, self.entries, self.labels, self._apply = session, list(entries), labels, apply
        self.channels: list[str] = []
        self._start = np.zeros((0, 0))
        self._filling = False
        #: the sample whose matrix is in the table
        self.shown = None
        head = QtWidgets.QLabel(
            "Each sample is compensated with its own matrix: the one in its file, or the "
            "FlowJo workspace's for it. Each cell is the percentage of the row's "
            "fluorochrome seen in the column's detector. Type over a number, then give "
            "the result to the samples selected in the list or to every sample.")
        head.setWordWrap(True)
        head.setProperty("role", "caption")
        self.samples = QtWidgets.QTreeWidget()
        self.samples.setHeaderLabels(["Sample", "Matrix"])
        self.samples.setRootIsDecorated(False)
        self.samples.setUniformRowHeights(True)
        self.samples.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.samples.setToolTip(
            "Where each sample's matrix comes from. file: its own FCS file. workspace: "
            "the FlowJo workspace. edited: changed here for that sample. shared: changed "
            "here and given to every sample.")
        for e in self.entries:
            self.samples.addTopLevelItem(QtWidgets.QTreeWidgetItem([e.name, ""]))
        self.samples.header().setStretchLastSection(False)
        self.samples.header().setSectionResizeMode(
            0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.samples.header().setSectionResizeMode(
            1, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)
        self.summary.setProperty("role", "caption")
        self.using = QtWidgets.QLabel("")
        self.using.setWordWrap(True)
        self.using.setProperty("role", "strong")
        self.table = QtWidgets.QTableWidget(0, 0)
        self.state = QtWidgets.QLabel("")
        self.state.setWordWrap(True)
        self.state.setProperty("role", "warning")
        self.apply_selected = QtWidgets.QPushButton("")
        self.apply_selected.setToolTip("Compensate the samples selected in the list with the "
                                       "matrix shown; the others keep theirs")
        self.apply_all = QtWidgets.QPushButton("Apply to every sample")
        self.apply_all.setToolTip("One matrix for every sample, in place of their own")
        self.back = QtWidgets.QPushButton("Back to the original")
        self.back.setToolTip("The samples selected in the list go back to the matrix of "
                             "their file, or of the workspace")
        self.discard = QtWidgets.QPushButton("Discard changes")
        self.discard.setToolTip("Forget what was typed and show the sample's matrix again")
        close = QtWidgets.QPushButton("Close")
        close.setAutoDefault(False)
        for b in (self.apply_selected, self.apply_all, self.back, self.discard):
            b.setAutoDefault(False)
        left = QtWidgets.QVBoxLayout()
        left.addWidget(self.samples, 1)
        left.addWidget(self.summary)
        right = QtWidgets.QVBoxLayout()
        right.addWidget(self.using)
        right.addWidget(self.table, 1)
        right.addWidget(self.state)
        middle = QtWidgets.QHBoxLayout()
        middle.addLayout(left, 0)
        middle.addLayout(right, 1)
        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(self.apply_selected)
        buttons.addWidget(self.apply_all)
        buttons.addWidget(self.back)
        buttons.addWidget(self.discard)
        buttons.addStretch(1)
        buttons.addWidget(close)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(head)
        layout.addLayout(middle, 1)
        layout.addLayout(buttons)
        self.samples.setMinimumWidth(230)
        self.samples.setMaximumWidth(300)

        self._filling = True
        for e in selected:
            if e in self.entries:
                self.samples.topLevelItem(self.entries.index(e)).setSelected(True)
        start = self.entries.index(entry) if entry in self.entries else 0
        self.samples.setCurrentItem(
            self.samples.topLevelItem(start),
            0, QtCore.QItemSelectionModel.SelectionFlag.Select
            | QtCore.QItemSelectionModel.SelectionFlag.Rows)
        self._filling = False
        self._sources()
        self._load(self.entries[start])
        n = len(self.channels)
        self.resize(min(560 + 78 * n, 1280), min(300 + 30 * max(n, 6), 760))

        self.samples.currentItemChanged.connect(self._current_changed)
        self.samples.itemSelectionChanged.connect(self._show_state)
        self.table.itemChanged.connect(self._show_state)
        self.apply_selected.clicked.connect(lambda: self.apply_to(self.selected()))
        self.apply_all.clicked.connect(lambda: self.apply_to(None))
        self.back.clicked.connect(lambda: self.reset(self.selected()))
        self.discard.clicked.connect(lambda: self._load(self.current()))
        close.clicked.connect(self.accept)

    # --- the samples ----------------------------------------------------------------

    def current(self):
        item = self.samples.currentItem()
        return self.entries[self.samples.indexOfTopLevelItem(item)] if item else self.entries[0]

    def selected(self) -> list:
        rows = sorted(self.samples.indexOfTopLevelItem(i) for i in self.samples.selectedItems())
        return [self.entries[r] for r in rows] or [self.current()]

    def select(self, names: list[str]) -> None:
        """Select these samples; the first becomes the one shown."""
        rows = [i for i, e in enumerate(self.entries) if e.name in names]
        self.samples.clearSelection()
        if rows:
            self.samples.setCurrentItem(self.samples.topLevelItem(rows[0]))
        for r in rows:
            self.samples.topLevelItem(r).setSelected(True)

    def source(self, name: str) -> str:
        """What the list says of this sample."""
        return self.samples.topLevelItem(
            [e.name for e in self.entries].index(name)).text(1)

    def _sources(self) -> None:
        for i, e in enumerate(self.entries):
            source = self.session.matrix_in_use(e)[1]
            item = self.samples.topLevelItem(i)
            item.setText(1, source or "none")
            font = item.font(1)
            font.setBold(source in ("edited", "shared"))
            item.setFont(1, font)
            item.setToolTip(1, f"{e.name}: {self.SOURCES[source]}")
        n, different = len(self.entries), self.session.distinct_matrices()
        self.summary.setText(
            "No sample has a matrix." if different == 0
            else f"The {n} samples all have the same matrix." if different == 1 and n > 1
            else "One sample, one matrix." if different == 1
            else f"{different} different matrices among the {n} samples.")

    def _current_changed(self, *_) -> None:
        if self._filling:
            return
        if self.changed():
            self._show_state()   # what was typed is kept, to give to other samples
        else:
            self._load(self.current())

    # --- the matrix -------------------------------------------------------------------

    def _load(self, entry) -> None:
        matrix = self.session.matrix_in_use(entry)[0]
        self.shown = entry
        self._filling = True
        self.channels = list(matrix.channels) if matrix is not None else []
        n = len(self.channels)
        self._start = (np.array(matrix.matrix, dtype=float) if matrix is not None
                       else np.zeros((0, 0)))
        self.table.clear()
        self.table.setRowCount(n)
        self.table.setColumnCount(n)
        labels = self.labels
        names = [(labels.get(c) or c).split()[0] if labels.get(c) else c for c in self.channels]
        self.table.setHorizontalHeaderLabels(self.channels)
        self.table.setVerticalHeaderLabels([f"{c}  ({m})" if m != c else c
                                            for c, m in zip(self.channels, names)])
        for i in range(n):
            for j in range(n):
                item = QtWidgets.QTableWidgetItem(self._shown(i, j))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight
                                      | QtCore.Qt.AlignmentFlag.AlignVCenter)
                if i == j:   # a fluorochrome is all of itself in its own detector
                    item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()
        self._filling = False
        self.state.setText("")
        self._show_state()

    def _shown(self, i: int, j: int) -> str:
        return f"{100 * self._start[i, j]:.5g}"

    def changed(self) -> bool:
        """Was anything typed that is not applied yet?"""
        n = len(self.channels)
        return any(self.table.item(i, j) is not None
                   and self.table.item(i, j).text() != self._shown(i, j)
                   for i in range(n) for j in range(n))

    def set_value(self, row: str, column: str, percent: float) -> None:
        self.table.item(self.channels.index(row), self.channels.index(column)).setText(f"{percent:g}")

    def matrix(self):
        from ..compensation import Compensation
        n = len(self.channels)
        values = np.empty((n, n))
        for i in range(n):
            for j in range(n):
                text = self.table.item(i, j).text()
                # a cell nobody typed in keeps its value to the last digit,
                # not the five it is shown with
                values[i, j] = (self._start[i, j] if text == self._shown(i, j)
                                else float(text.replace(",", ".")) / 100.0)
        return Compensation(self.channels, values, "edited")

    def _show_state(self, *_) -> None:
        if self._filling or self.shown is None:
            return
        changed, chosen = self.changed(), self.selected()
        source = self.session.matrix_in_use(self.shown)[1]
        if changed:
            self.using.setText(f"{self.shown.name}'s matrix with your changes — not applied "
                               "yet. Select the samples it is for, then apply it.")
        elif source == "workspace" and getattr(self.shown.wsp, "compensation", None) is None:
            self.using.setText(f"{self.shown.name} is compensated with the matrix the FlowJo "
                               "workspace's samples share: the workspace does not list it.")
        else:
            self.using.setText(f"{self.shown.name} is compensated with {self.SOURCES[source]}.")
        self.apply_selected.setText(f"Apply to {chosen[0].name}" if len(chosen) == 1
                                    else f"Apply to the {len(chosen)} selected samples")
        has = bool(self.channels)
        self.apply_selected.setEnabled(has)
        self.apply_all.setEnabled(has)
        self.discard.setEnabled(changed)
        self.back.setEnabled(any(self.session.matrix_in_use(e)[1] in ("edited", "shared")
                                 for e in chosen))

    # --- what the buttons do ------------------------------------------------------------

    def _checked(self):
        try:
            m = self.matrix()
            np.linalg.inv(m.matrix)
        except ValueError:
            self.state.setText("Every cell has to be a number")
            return None
        except np.linalg.LinAlgError:
            self.state.setText("This matrix cannot be inverted: two fluorochromes "
                               "would be indistinguishable")
            return None
        return m

    def apply_to(self, entries: list | None) -> bool:
        """Give the matrix in the table to these samples — None: to every
        sample, as the one they share."""
        m = self._checked()
        if m is None:
            return False
        try:
            self._apply(m, entries)
        except ValueError as error:
            self.state.setText(str(error)[:1].upper() + str(error)[1:])
            return False
        self._done()
        return True

    def reset(self, entries: list | None) -> None:
        """These samples back on their file's or the workspace's matrix —
        None: every sample."""
        if entries is not None and len(entries) == len(self.entries):
            entries = None
        self._apply(None, entries)
        self._done()

    def _done(self) -> None:
        self._sources()
        self._load(self.current())


class ExportDialog(QtWidgets.QDialog):
    def __init__(self, population: str, selected: int, start: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Export {population}")
        self.format = QtWidgets.QComboBox()
        self.format.addItem("FCS 3.1", "fcs")
        self.format.addItem("CSV", "csv")
        self.all_samples = QtWidgets.QCheckBox("Every sample ticked in Samples")
        self.all_samples.setToolTip(f"Unticked: the {selected} selected here")
        self.compensated = QtWidgets.QCheckBox("Compensated values")
        self.compensated.setChecked(True)
        self.compensated.setToolTip(
            "Ticked: fluorescence as this program works with it, in columns named Comp-…. "
            "Unticked: as acquired; an FCS file then carries the matrix in use.")
        self.folder = QtWidgets.QLineEdit(os.path.join(start, "exported") if start else "")
        browse = QtWidgets.QPushButton("Choose…")
        browse.clicked.connect(self._browse)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.folder, 1)
        row.addWidget(browse)
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                         | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        box.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setText("Export")
        box.accepted.connect(self.accept)
        box.rejected.connect(self.reject)
        form = QtWidgets.QFormLayout(self)
        form.addRow("Format", self.format)
        form.addRow(self.all_samples)
        form.addRow(self.compensated)
        form.addRow("Into", row)
        form.addRow(box)
        self.resize(560, 200)

    def _browse(self) -> None:
        chosen = QtWidgets.QFileDialog.getExistingDirectory(self, "Export into", self.folder.text())
        if chosen:
            self.folder.setText(chosen)
