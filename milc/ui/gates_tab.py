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

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets

from .. import figures
from ..compensation import PREFIXES
from ..gating import GateNode, PolygonGate, RectangleGate
from ..session import Session
from ..transforms import Biex, decade_ticks
from . import style, theme

HISTOGRAM = "#histogram"
MAX_POINTS = 40_000


class _DrawBox(pg.ViewBox):
    """Pans and zooms as usual; in a drawing mode the mouse draws."""

    sigRectangle = QtCore.pyqtSignal(float, float, float, float)
    sigPolygon = QtCore.pyqtSignal(object)

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
        if self.mode != "rect" or ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return super().mouseDragEvent(ev, axis)
        ev.accept()
        a = self.mapSceneToView(ev.buttonDownScenePos())
        b = self.mapSceneToView(ev.scenePos())
        self._outline([a.x(), b.x(), b.x(), a.x(), a.x()], [a.y(), a.y(), b.y(), b.y(), a.y()])
        if ev.isFinish():
            self.cancel()
            if a.x() != b.x():
                self.sigRectangle.emit(min(a.x(), b.x()), min(a.y(), b.y()),
                                       max(a.x(), b.x()), max(a.y(), b.y()))

    def mouseClickEvent(self, ev):
        if self.mode != "poly" or ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return super().mouseClickEvent(ev)
        ev.accept()
        p = self.mapSceneToView(ev.scenePos())
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
        tools.addWidget(self.rename_button, 1, 0)
        tools.addWidget(self.delete_button, 1, 1)
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
        bar = QtWidgets.QHBoxLayout()
        for text, w in (("X", self.x_axis), ("Y", self.y_axis)):
            cap = QtWidgets.QLabel(text)
            cap.setProperty("role", "caption")
            bar.addWidget(cap)
            bar.addWidget(w)
            bar.addSpacing(10)
        bar.addWidget(self.caption, 1)
        self.box = _DrawBox()
        self.plot = pg.PlotWidget(viewBox=self.box)
        self.plot.setBackground(theme.background())
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.empty = QtWidgets.QLabel("Open a folder in the Samples tab, then choose a sample here.")
        self.empty.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.empty.setProperty("role", "hint")
        self.stack = QtWidgets.QStackedWidget()
        self.stack.addWidget(self.empty)
        self.stack.addWidget(self.plot)
        right = QtWidgets.QVBoxLayout()
        right.addLayout(bar)
        right.addWidget(self.stack, 1)

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.addWidget(left_panel)
        layout.addLayout(right, 1)

        self.samples.currentRowChanged.connect(self._sample_chosen)
        self.tree.currentItemChanged.connect(self._node_chosen)
        self.x_axis.currentIndexChanged.connect(self._axes_changed)
        self.y_axis.currentIndexChanged.connect(self._axes_changed)
        self.rect_button.toggled.connect(lambda on: self._tool("rect", on))
        self.poly_button.toggled.connect(lambda on: self._tool("poly", on))
        self.box.sigRectangle.connect(self._rectangle_drawn)
        self.box.sigPolygon.connect(self._polygon_drawn)
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
        self._filling = True
        self.samples.clear()
        self.entries = self.session.usable()
        for e in self.entries:
            self.samples.addItem(f"{e.name}    {e.events:,d}")
        self._filling = False
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

    @staticmethod
    def _scale(dim: str):
        return Biex() if dim.startswith(PREFIXES) else None

    def _display(self, dim: str) -> np.ndarray:
        x = self.sample.column(dim).astype(np.float64)
        scale = self._scale(dim)
        return scale(x) if scale is not None else x

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
            if self.x_axis.findData(dims[0]) >= 0:
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
        if ydim == HISTOGRAM:
            counts, edges = np.histogram(x[rows], bins=256, range=xr)
            fill = QtGui.QColor(theme.accent())
            fill.setAlpha(90)
            self.plot.addItem(pg.PlotCurveItem(edges, counts, stepMode="center", fillLevel=0,
                                               brush=pg.mkBrush(fill),
                                               pen=pg.mkPen(theme.accent(), width=1.2)))
            self.plot.setLabel("left", "Events")
            self.plot.getAxis("left").setTicks(None)
            yr = (0.0, float(counts.max()) * 1.08 + 1)
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
        self.plot.setXRange(*xr, padding=0.01)
        self.plot.setYRange(*yr, padding=0.01)
        self._draw_gates(xdim, ydim, xr, yr)
        parent = self.node.parent
        of = (f", {100 * rows.size / max(int(self.sample.mask(parent.path).sum()), 1):.1f}% "
              f"of {parent.name if parent.gate is not None else 'all events'}"
              if parent is not None else "")
        more = f" — {shown.size:,d} drawn" if shown.size < rows.size else ""
        name = self.node.name if self.node.gate is not None else "All events"
        self.caption.setText(f"{self.entry.name} · {name}: {rows.size:,d} events{of}{more}")

    def _range(self, dim: str) -> tuple[float, float]:
        if self._scale(dim) is not None:
            return (0.0, 1.0)
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
        if isinstance(item, pg.PolyLineROI):
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
        other = self.poly_button if mode == "rect" else self.rect_button
        if on:
            other.setChecked(False)
            if mode == "poly" and self.y_axis.currentData() == HISTOGRAM:
                self.sigStatus.emit("A polygon needs two parameters: choose a Y axis")
                self.poly_button.setChecked(False)
                return
            self.box.set_mode(mode)
            self.sigStatus.emit("Drag across the events to gate" if mode == "rect" else
                                "Click each corner; double-click to close. Esc cancels.")
        elif self.box.mode == mode:
            self.box.set_mode(None)
            self.sigStatus.emit("")

    def _stop_drawing(self) -> None:
        self.rect_button.setChecked(False)
        self.poly_button.setChecked(False)
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
        self.node.name = name
        self._commit()

    def delete_selected(self) -> None:
        """The gate being looked at, and everything gated inside it."""
        node = self.node
        if node is None or node.gate is None:
            return
        parent = node.parent
        parent.children.remove(node)
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
        self.reset_button.setEnabled(True)
