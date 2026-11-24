"""
The window's side of `milc.views`: the choice among the ways several
samples are shown on one parameter, violins and boxes on the Gates plot,
and the two dialogs — the matrix of every pair and the backgate.

The dialogs show what `figures.py` draws, as the Figure dialog's preview
does: a grid of a hundred cells, or a row of small plots, is a figure
already, and drawing it a second way for the screen would be a second
thing to keep right.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets

from .help_window import describe
from .. import figures, views
from . import palette, theme

OVERLAID, STAGGERED = "overlaid", "staggered"
MODES = [(OVERLAID, "Overlaid"), (STAGGERED, "Staggered"),
         (views.VIOLIN, "Violins"), (views.BOX, "Boxes")]


class OverlayChoice(QtWidgets.QToolButton):
    """One button for how several samples' histograms are shown: it reads
    as the choice made and opens the others. A row of four more widgets
    would not fit the bar."""

    sigChanged = QtCore.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        self.setToolTip(
            "Several samples on one parameter: histograms on one baseline or one above "
            "the other, or each sample's distribution as a violin or a box")
        menu = QtWidgets.QMenu(self)
        menu.setToolTipsVisible(True)
        self._group = QtGui.QActionGroup(self)
        self._modes = {}
        for mode, text in MODES:
            action = menu.addAction(text)
            action.setCheckable(True)
            action.setData(mode)
            self._group.addAction(action)
            self._modes[mode] = action
            action.toggled.connect(self._chosen)
        self.staggered = self._modes[STAGGERED]
        menu.addSeparator()
        self.pool = menu.addAction("One for each group")
        self.pool.setCheckable(True)
        self.pool.setToolTip("Violins and boxes: the samples of each value of the "
                             "“Colour by” column pooled into one")
        self.pool.toggled.connect(lambda _on: self.sigChanged.emit())
        self.setMenu(menu)
        self._modes[OVERLAID].setChecked(True)

    def _chosen(self, on: bool) -> None:
        if on:   # the one let go says nothing: one choice, one redraw
            self.setText(self._group.checkedAction().text())
            self.sigChanged.emit()

    def mode(self) -> str:
        return self._group.checkedAction().data()

    def pooled(self) -> bool:
        return self.pool.isChecked()

    def select(self, mode: str, pooled: bool = False) -> None:
        self.pool.setChecked(pooled)
        self._modes[mode].setChecked(True)

    def enable(self, on: bool) -> None:
        self.setEnabled(on)
        self._group.setEnabled(on)
        self.pool.setEnabled(on)


def draw_distributions(plot, c: views.DistributionContent) -> tuple[float, float]:
    """Violins or boxes on a pyqtgraph plot, a row each, the parameter
    across as on a histogram — so a range gate stays where it is. Returns
    the range of Y they need."""
    names = c.groups()
    colour = dict(zip(names, palette.categorical(len(names))))
    legend = plot.plotItem.legend or plot.plotItem.addLegend(offset=(-8, 8))
    legend.clear()
    legend.setBrush(theme.legend_brush())
    legend.setPen(theme.legend_pen())
    legend.setLabelTextColor(theme.foreground())
    legend.setVisible(1 < len(names) <= 14 and len(names) != len(c.series))
    ink = pg.mkPen(theme.foreground(), width=1.4)
    listed = set()
    for k, s in enumerate(c.series):
        m = s.summary
        if not m.n:
            continue
        pen = pg.mkPen(colour[s.group], width=1.4)
        fill = QtGui.QColor(colour[s.group])
        fill.setAlpha(90)
        if c.kind == views.VIOLIN:
            x, half = views.violin(s)
            shape = pg.PlotCurveItem(np.concatenate([x, x[::-1], x[:1]]),
                                     np.concatenate([k + half, (k - half)[::-1], k + half[:1]]),
                                     pen=pen, brush=pg.mkBrush(fill), fillLevel=float(k))
            plot.addItem(shape)
            plot.addItem(pg.PlotCurveItem([m.q1, m.q3], [k, k],
                                          pen=pg.mkPen(theme.foreground(), width=4)))
            plot.addItem(pg.ScatterPlotItem([m.median], [k], size=7, pen=ink,
                                            brush=pg.mkBrush(theme.background())))
        else:
            h, cap = views._BOX_HALF, views._BOX_HALF / 2
            shape = pg.PlotCurveItem([m.q1, m.q3, m.q3, m.q1, m.q1],
                                     [k - h, k - h, k + h, k + h, k - h],
                                     pen=pen, brush=pg.mkBrush(fill), fillLevel=float(k))
            plot.addItem(shape)
            plot.addItem(pg.PlotCurveItem(
                [m.low, m.q1, m.q3, m.high, m.low, m.low, m.high, m.high],
                [k, k, k, k, k - cap, k + cap, k - cap, k + cap], pen=pen, connect="pairs"))
            plot.addItem(pg.PlotCurveItem([m.median, m.median], [k - h, k + h], pen=ink))
            dots = QtGui.QColor(colour[s.group])
            dots.setAlpha(150)
            plot.addItem(pg.ScatterPlotItem(s.outliers, np.full(s.outliers.size, float(k)),
                                            size=3, pen=None, brush=pg.mkBrush(dots)))
        if s.group not in listed:
            listed.add(s.group)
            legend.addItem(shape, s.group)
    plot.getAxis("left").setTicks([[(k, s.label) for k, s in enumerate(c.series)]])
    plot.setLabel("left", "")
    return (-0.6, len(c.series) - 0.4)


class _FigureDialog(QtWidgets.QDialog):
    """A content shown as `figures.py` draws it, with the way to Figure…"""

    #: the size the content is drawn at, before it is fitted to the window
    SIZE_MM = (170.0, 160.0)

    def __init__(self, title: str, name: str, parent=None):
        super().__init__(parent)
        describe(self, "gates-tab")
        self.setWindowTitle(title)
        self.name = name
        self.content = None
        self.preview = QtWidgets.QLabel()
        self.preview.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.preview.setMinimumSize(240, 200)
        self.note = QtWidgets.QLabel("")
        self.note.setProperty("role", "caption")
        self.note.setWordWrap(True)
        self.figure_button = QtWidgets.QPushButton("Figure…")
        self.close_button = QtWidgets.QPushButton("Close")
        self.figure_button.clicked.connect(self.make_figure)
        self.close_button.clicked.connect(self.accept)
        self.buttons = QtWidgets.QHBoxLayout()
        self.buttons.addWidget(self.note, 1)
        self.buttons.addWidget(self.figure_button)
        self.buttons.addWidget(self.close_button)

    def size_mm(self) -> tuple[float, float]:
        return self.SIZE_MM

    def show_content(self, content) -> None:
        self.content = content
        self.figure_button.setEnabled(content is not None)
        if content is None:
            self.preview.clear()
            return
        width, height = self.size_mm()
        # the ramp on the screen, which is also the one Figure… starts from
        look = figures.FigureStyle(width_mm=width, height_mm=height,
                                   colormap=palette.ramp_name(),
                                   colormap_colors=palette.custom_colors())
        pix = QtGui.QPixmap()
        pix.loadFromData(figures.preview_png(content, look), "PNG")
        room = self.preview.size()
        if self.isVisible() and (pix.width() > room.width() or pix.height() > room.height()):
            pix = pix.scaled(room, QtCore.Qt.AspectRatioMode.KeepAspectRatio,
                             QtCore.Qt.TransformationMode.SmoothTransformation)
        self.preview.setPixmap(pix)

    def make_figure(self):
        from .figure_dialog import FigureDialog
        if self.content is None:
            return None
        dialog = FigureDialog(self.content, self, self.name)
        dialog.exec()
        return dialog


class MatrixDialog(_FigureDialog):
    """Every pair of the parameters ticked, for one population of one
    sample. `build(dims)` gives the content for a list of them."""

    def __init__(self, title: str, parameters: list[tuple[str, str]], chosen: list[str],
                 build, parent=None):
        super().__init__(title, "matrix", parent)
        self._build = build
        self._timer = QtCore.QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(250)
        self._timer.timeout.connect(self.render)
        self.list = QtWidgets.QListWidget()
        self.list.setToolTip(f"At most {views.MAX_PARAMETERS}: the first ticked are drawn")
        for label, dim in parameters:
            item = QtWidgets.QListWidgetItem(label)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, dim)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.CheckState.Checked if dim in chosen
                               else QtCore.Qt.CheckState.Unchecked)
            self.list.addItem(item)
        self.list.setFixedWidth(210)
        self.list.itemChanged.connect(lambda _item: self._timer.start())
        body = QtWidgets.QHBoxLayout()
        body.addWidget(self.list)
        body.addWidget(self.preview, 1)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(body, 1)
        layout.addLayout(self.buttons)
        self.resize(1010, 800)
        self.render()

    def ticked(self) -> list[str]:
        return [self.list.item(i).data(QtCore.Qt.ItemDataRole.UserRole)
                for i in range(self.list.count())
                if self.list.item(i).checkState() == QtCore.Qt.CheckState.Checked]

    def set_parameters(self, dims: list[str]) -> None:
        self.list.blockSignals(True)
        for i in range(self.list.count()):
            item = self.list.item(i)
            on = item.data(QtCore.Qt.ItemDataRole.UserRole) in dims
            item.setCheckState(QtCore.Qt.CheckState.Checked if on
                               else QtCore.Qt.CheckState.Unchecked)
        self.list.blockSignals(False)
        self.render()

    def render(self) -> None:
        ticked = self.ticked()
        dims = ticked[:views.MAX_PARAMETERS]
        if len(dims) < 2:
            self.show_content(None)
            self.note.setText("Tick two parameters or more")
            return
        self.show_content(self._build(dims))
        left = (f" Only the first {len(dims)} of the {len(ticked)} ticked are drawn."
                if len(dims) < len(ticked) else "")
        self.note.setText("Each cell: the parameter of its column across, of its row up; "
                          "colour is the density of events, on a logarithmic ramp of the "
                          "cell's own." + left)


class BackgateDialog(_FigureDialog):
    """A population on each gate above it."""

    def __init__(self, title: str, content: views.BackgateContent, parent=None):
        super().__init__(title, "backgate", parent)
        self._panels = len(content.panels)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.preview, 1)
        layout.addLayout(self.buttons)
        cols, rows = self._grid()
        self.resize(min(110 + 250 * cols, 1160), 130 + 240 * rows)
        self.show_content(content)
        few = (f" {content.drawn:,d} of its {content.events:,d} events are drawn."
               if content.drawn < content.events else "")
        self.note.setText(
            f"Grey: every event the gate was applied to. Colour: the events of "
            f"{content.population}, wherever they fall.{few}")

    def _grid(self) -> tuple[int, int]:
        cols = min(self._panels, 4)
        return cols, -(-self._panels // cols)

    def size_mm(self) -> tuple[float, float]:
        cols, rows = self._grid()
        return (58.0 * cols + 30.0, 54.0 * rows)
