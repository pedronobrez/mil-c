"""
What a batch normalisation did, to look at: for one marker, the reference
samples of every batch drawn over one another, before and after.

The events shown were held out of the fit, so "after" is not the
agreement the curves were built to have but the one they achieve on
events they have not seen.
"""

from __future__ import annotations

import pyqtgraph as pg
from PyQt6 import QtGui, QtWidgets

from .help_window import describe
from .. import density, figures
from ..cluster_stats import ordered
from ..transforms import decade_ticks, display_range
from . import palette, style, theme


class NormalizationDialog(QtWidgets.QDialog):
    def __init__(self, check: dict, parent=None):
        super().__init__(parent)
        describe(self, "normalisation")
        self.setWindowTitle("Batch normalisation: before and after")
        self.check = check
        self.marker = QtWidgets.QComboBox()
        self.marker.addItems(check["markers"])
        note = QtWidgets.QLabel(
            "The reference samples of each batch, on events held out of the fit. Where the "
            "curves lie on one another the batches agree. Each curve is drawn to its own peak.")
        note.setWordWrap(True)
        note.setProperty("role", "caption")
        self.plots = {}
        layout = QtWidgets.QVBoxLayout(self)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("Marker"))
        top.addWidget(self.marker)
        top.addStretch(1)
        self.figure_before = QtWidgets.QPushButton("Figure: before…")
        self.figure_after = QtWidgets.QPushButton("Figure: after…")
        top.addWidget(self.figure_before)
        top.addWidget(self.figure_after)
        layout.addLayout(top)
        layout.addWidget(note)
        for key, title in (("before", "Before"), ("after", "After")):
            label = QtWidgets.QLabel(title)
            label.setFont(style.eyebrow_font())
            label.setProperty("role", "caption")
            plot = pg.PlotWidget()
            plot.setBackground(theme.background())
            plot.setMenuEnabled(False)
            plot.hideButtons()
            plot.setMouseEnabled(False, False)
            layout.addWidget(label)
            layout.addWidget(plot, 1)
            self.plots[key] = plot
        close = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        close.rejected.connect(self.reject)
        layout.addWidget(close)
        self.resize(820, 680)
        self.marker.currentIndexChanged.connect(self.draw)
        self.figure_before.clicked.connect(lambda: self._figure("before"))
        self.figure_after.clicked.connect(lambda: self._figure("after"))
        self.draw()

    def content(self, key: str) -> figures.OverlayContent:
        """One of the two panels, as something a figure can be made of."""
        k = self.marker.currentIndex()
        scale = self.check["scales"][k]
        events = self.check[key]
        lo, hi = display_range(scale) if scale is not None else (
            min(float(v[:, k].min()) for v in events.values()),
            max(float(v[:, k].max()) for v in events.values()))
        series = [figures.Series(b, b, *density.histogram(events[b][:, k], lo, hi))
                  for b in ordered(events)]
        return figures.OverlayContent(
            series, self.check["markers"][k], "% of the tallest bin",
            decade_ticks(scale) if scale is not None else None, xrange=(lo, hi))

    def draw(self) -> None:
        for key, plot in self.plots.items():
            c = self.content(key)
            plot.clear()
            theme.style_axes(plot.plotItem)
            legend = plot.plotItem.legend or plot.plotItem.addLegend(offset=(-8, 8))
            legend.clear()
            legend.setBrush(theme.legend_brush())
            legend.setPen(theme.legend_pen())
            legend.setLabelTextColor(theme.foreground())
            for s, colour in zip(c.series, palette.categorical(len(c.series))):
                fill = QtGui.QColor(colour)
                fill.setAlpha(40)
                curve = pg.PlotCurveItem(s.x, s.y, fillLevel=0, brush=pg.mkBrush(fill),
                                         pen=pg.mkPen(colour, width=1.6))
                plot.addItem(curve)
                legend.addItem(curve, s.label or "—")
            plot.getAxis("bottom").setTicks([c.xticks] if c.xticks else None)
            plot.setLabel("bottom", c.xlabel)
            plot.setXRange(*c.xrange, padding=0.01)
            plot.setYRange(0, 106, padding=0)

    def _figure(self, key: str) -> None:
        from .figure_dialog import FigureDialog
        FigureDialog(self.content(key), self, f"normalisation_{key}").exec()
