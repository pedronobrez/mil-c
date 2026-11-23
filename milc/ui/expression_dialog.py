"""
Differential expression within the clusters, asked from the Clusters tab:
which two groups, which markers, adjusted for what, and the answer as a
table to sort, copy and save, a volcano, a heatmap of the differences,
and a paragraph that says what was done and what it does not show.
"""

from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets

from .. import compare, expression, figures
from ..expression_figures import (DIVERGING, UNTESTED, DifferenceHeatmap, ExpressionVolcano,
                                  short, stars)
from . import theme
from .abundance_dialog import VolcanoView
from .results_dialog import SETTING_FOLDER, _text
from .settings import settings
from .statistics_tab import _Cell

METHOD = (
    "Each sample gives one number for each cluster and marker: the median of the marker "
    "among its cells of the cluster, {scale}. Those medians are compared between {first} "
    "and {second} by a linear model{adjusting}{weighing} and limma's moderated t, as "
    "diffcyt's DS-limma does{within}. A sample with fewer than {cells} cells of a cluster "
    "gives no median for it, and a cluster with fewer than {samples} such samples in "
    "either group is not tested. Each sample is one observation; {tests} test(s), "
    "{correction}. "
    "The limits: the cells are those that were concatenated; a marker the clusters were "
    "built on is held in place within each cluster by construction, so its differences "
    "show mostly as abundance; a median does not move when a minority of the cluster "
    "switches a marker on; and the moderated t, which agrees with limma's to 8 digits, "
    "gave p < 0.05 about as often as it should in simulations without weights, and "
    "somewhat too often (about 7%) with them.")
SCALES = {
    "biex": "on the biexponential scale of the plots as numbers from 0 to 1023 (where the "
            "scale is logarithmic, 100 is about half a decade)",
    "arcsinh": "as arcsinh(x / 150) (where the marker is bright, 0.69 is a doubling)"}
SCALE_NAMES = {"biex": "Biexponential", "arcsinh": "arcsinh (cofactor 150)"}


class _TurnedAxis(pg.AxisItem):
    """A bottom axis whose names read upwards: a dozen markers in three
    panels have no room side by side."""

    def fit(self, names: list[str]) -> None:
        # the first and last names would overlap the ends if they lay flat
        self.setStyle(hideOverlappingLabels=False)
        metrics = QtGui.QFontMetricsF(self.style["tickFont"] or QtGui.QFont())
        self.setHeight(max((metrics.horizontalAdvance(n) for n in names), default=0.0) + 14)

    def drawPicture(self, p, axisSpec, tickSpecs, textSpecs):
        super().drawPicture(p, axisSpec, tickSpecs, [])
        for rect, _, text in textSpecs:
            p.save()
            p.translate(rect.center().x(), rect.top())
            p.rotate(-90)
            p.drawText(QtCore.QRectF(-1000, -rect.height() / 2, 1000, rect.height()),
                       int(QtCore.Qt.AlignmentFlag.AlignRight
                           | QtCore.Qt.AlignmentFlag.AlignVCenter), text)
            p.restore()


class DifferenceView(pg.GraphicsLayoutWidget):
    """Clusters down, markers across, the difference in median as colour
    — red higher in the first group, blue in the second, grey not tested
    — and stars where the adjusted p is below alpha. A panel for each
    value of "within each"."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._content: DifferenceHeatmap | None = None
        self.panels: list[pg.PlotItem] = []
        self.setBackground(theme.background())

    def content(self):
        return self._content if self._content and self._content.results else None

    def show_content(self, content: DifferenceHeatmap | None) -> None:
        self._content = content
        self.clear()
        self.panels = []
        self.setBackground(theme.background())
        if self.content() is None:
            return
        clusters, markers = content.clusters(), content.markers()
        rows, cols = len(clusters), len(markers)
        lut = figures.lookup_table(DIVERGING)
        grey = np.array([*QtGui.QColor(UNTESTED).getRgb()[:3], 255], dtype=np.uint8)
        limit = content.limit()
        font = QtGui.QFont()
        font.setPointSize(11)
        panels = content.panels()
        turned = cols * len(panels) > 10
        for k, panel in enumerate(panels):
            plot = self.addPlot(row=0, col=k, axisItems=(
                {"bottom": _TurnedAxis("bottom")} if turned else {}))
            plot.setMenuEnabled(False)
            plot.hideButtons()
            plot.setMouseEnabled(False, False)
            theme.style_axes(plot)
            if turned:
                plot.getAxis("bottom").fit([short(m) for m in markers])
            if panel is not None:
                plot.setTitle(str(panel) or "—", color=theme.foreground(), size="10pt")
            difference, p = content.cells(panel)
            tested = np.isfinite(difference)
            shade = np.clip((np.where(tested, difference, 0.0) / limit + 1) / 2, 0, 1)
            rgba = lut[(shade * (len(lut) - 1)).astype(int)]
            rgba[~tested] = grey
            # the colours as they are: left to itself the item stretches them
            img = pg.ImageItem(np.ascontiguousarray(rgba.transpose(1, 0, 2)[:, ::-1]),
                               levels=(0, 255))
            img.setRect(QtCore.QRectF(0, 0, cols, rows))
            plot.addItem(img)
            for r, c in zip(*np.nonzero(tested)):
                mark = stars(p[r, c], content.alpha)
                if mark:
                    dark = abs(difference[r, c]) > 0.6 * limit
                    t = pg.TextItem(mark, anchor=(0.5, 0.4),
                                    color="#ffffff" if dark else "#111111")
                    t.setFont(font)
                    t.setPos(c + 0.5, rows - r - 0.5)
                    plot.addItem(t)
            plot.getAxis("bottom").setTicks([[(c + 0.5, short(m))
                                              for c, m in enumerate(markers)]])
            plot.getAxis("left").setTicks([[(rows - r - 0.5, n) if k == 0 else
                                            (rows - r - 0.5, "")
                                            for r, n in enumerate(clusters)]])
            if k == 0:
                plot.setLabel("left", "Cluster")
            plot.setXRange(0, cols, padding=0)
            plot.setYRange(0, rows, padding=0)
            self.panels.append(plot)


class ExpressionDialog(QtWidgets.QDialog):
    """`medians` is what each sample's cells of each cluster express;
    `groups` names each sample's group (`group_name` is the column they
    come from); `columns` are the other columns of the samples that may
    be adjusted for; `facets`, where given, are the values of the tab's
    "within each"; `chosen` the markers ticked to begin with (all of
    them where none is named); `reload`, where given, makes the medians
    again on another scale."""

    def __init__(self, medians: expression.Medians, groups: list[str],
                 group_name: str = "Group", columns: dict[str, list[str]] | None = None,
                 facets: list[str] | None = None, facet_name: str = "",
                 chosen: list[str] | None = None, reload=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Differential expression")
        self.medians, self.groups, self.columns = medians, list(groups), dict(columns or {})
        self.facets, self.facet_name, self._reload = facets, facet_name, reload
        self.results: list[expression.Result] = []
        self.result = expression.to_table([])

        self.first, self.second = QtWidgets.QComboBox(), QtWidgets.QComboBox()
        levels = compare.ordered(groups)
        for combo in (self.first, self.second):
            for level in levels:
                combo.addItem(level or "—", level)
        self.second.setCurrentIndex(min(1, len(levels) - 1))
        self.second.setToolTip("The reference: a positive difference is higher in the first")
        self.adjust = QtWidgets.QComboBox()
        self.adjust.addItem("—", "")
        for name in self.columns:
            self.adjust.addItem(name, name)
        self.adjust.setToolTip(
            "A batch or a day: its effect is estimated beside the group's and taken out. "
            "The column that says which animal or donor a sample is of makes the "
            "comparison a paired one.")
        self.markers_button = QtWidgets.QPushButton()
        self.markers_button.setToolTip(
            "The markers compared. Those the clusters were built on differ little within "
            "a cluster by construction; the others are where a cluster's state shows.")
        menu = QtWidgets.QMenu(self.markers_button)
        box = QtWidgets.QWidget()
        ticks = QtWidgets.QVBoxLayout(box)
        ticks.setContentsMargins(10, 6, 10, 6)
        ticks.setSpacing(4)
        self.marker_boxes: list[QtWidgets.QCheckBox] = []
        for marker in medians.markers:
            tick = QtWidgets.QCheckBox(marker)
            tick.setChecked(chosen is None or marker in chosen)
            ticks.addWidget(tick)
            self.marker_boxes.append(tick)
        action = QtWidgets.QWidgetAction(menu)
        action.setDefaultWidget(box)
        menu.addAction(action)
        self.markers_button.setMenu(menu)
        self.scale = QtWidgets.QComboBox()
        for key, text in SCALE_NAMES.items():
            self.scale.addItem(text, key)
        self.scale.setCurrentIndex(max(self.scale.findData(medians.scale), 0))
        self.scale.setToolTip("The scale the medians are taken on, and the difference is in")
        self.min_cells = QtWidgets.QSpinBox()
        self.min_cells.setRange(1, 100000)
        self.min_cells.setValue(expression.MIN_CELLS)
        self.min_cells.setToolTip(
            "A sample with fewer cells of a cluster than this gives no median for it")
        self.min_samples = QtWidgets.QSpinBox()
        self.min_samples.setRange(2, 1000)
        self.min_samples.setValue(expression.MIN_SAMPLES)
        self.min_samples.setToolTip(
            "A cluster with fewer samples than this in either group, after the line "
            "above, is not tested")
        self.weights = QtWidgets.QCheckBox("Weigh by cells")
        self.weights.setToolTip(
            "A median of 5,000 cells counts for more than one of 50, as in diffcyt. Right "
            "where the medians differ because of how few cells they are of; where the "
            "animals themselves differ, as they usually do, it makes the test somewhat "
            "too generous.")
        self.correction = QtWidgets.QComboBox()
        for name in compare.CORRECTIONS:
            self.correction.addItem(name, name)
        self.correction.setToolTip("Over every cluster and marker tested")

        def row(*items) -> QtWidgets.QHBoxLayout:
            lay = QtWidgets.QHBoxLayout()
            for label, widget in items:
                if label:
                    cap = QtWidgets.QLabel(label)
                    cap.setProperty("role", "caption")
                    lay.addWidget(cap)
                lay.addWidget(widget)
                lay.addSpacing(10)
            lay.addStretch(1)
            return lay

        bar = row((f"Compare {group_name}", self.first), ("against", self.second),
                  ("Adjust for", self.adjust), ("Markers", self.markers_button),
                  *([("Scale", self.scale)] if reload is not None else []))
        options = row(("Cells in a sample, at least", self.min_cells),
                      ("Samples in a group, at least", self.min_samples),
                      (None, self.weights), ("Correction", self.correction))

        self.table = QtWidgets.QTableWidget()
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.volcano = VolcanoView()
        self.heatmap = DifferenceView()
        self.legend = QtWidgets.QLabel()
        self.legend.setProperty("role", "caption")
        self.legend.setWordWrap(True)
        page = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(page)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.heatmap, 1)
        lay.addWidget(self.legend)
        self.views = QtWidgets.QTabWidget()
        self.views.addTab(page, "Heatmap")
        self.views.addTab(self.volcano, "Volcano")
        figure = QtWidgets.QToolButton()
        figure.setText("Figure…")
        figure.setAutoRaise(True)
        self.views.setCornerWidget(figure)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        split.addWidget(self.table)
        split.addWidget(self.views)
        split.setChildrenCollapsible(False)
        split.setSizes([680, 560])

        self.statement = QtWidgets.QLabel()
        self.statement.setWordWrap(True)
        self.statement.setProperty("role", "caption")
        export = QtWidgets.QPushButton("Export CSV…")
        copy = QtWidgets.QPushButton("Copy")
        close = QtWidgets.QPushButton("Close")
        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(export)
        buttons.addWidget(copy)
        buttons.addStretch(1)
        buttons.addWidget(close)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(bar)
        layout.addLayout(options)
        layout.addWidget(split, 1)
        layout.addWidget(self.statement)
        layout.addLayout(buttons)
        self.resize(1280, 760)

        for combo in (self.first, self.second, self.adjust, self.correction):
            combo.currentIndexChanged.connect(self.run)
        self.scale.currentIndexChanged.connect(self._rescale)
        for spin in (self.min_cells, self.min_samples):
            spin.valueChanged.connect(self.run)
        self.weights.toggled.connect(self.run)
        for tick in self.marker_boxes:
            tick.toggled.connect(self.run)
        figure.clicked.connect(self.make_figure)
        export.clicked.connect(self._ask_export)
        copy.clicked.connect(self.copy)
        close.clicked.connect(self.accept)
        self.run()

    @classmethod
    def for_tab(cls, tab) -> "ExpressionDialog":
        """What the Clusters tab is comparing now, in the clusters as it
        shows them. Every stained marker is offered; the ones the
        clusters were not built on, where there are any, are the ones
        ticked."""
        session = tab.session
        ds = session.dataset
        used = ds.runs.get(tab.key, {}).get("channels") or session.default_channels()
        others = tab._other_markers()
        channels = list(used) + others

        def load(scale: str | None = None) -> expression.Medians:
            return session.cluster_medians(channels, tab.key, scale)

        facets = tab._facets()
        taken = {tab.group_by.currentData(), tab.facet_by.currentData() if facets else ""}
        columns = {k: [m.get(k, "") for m in ds.meta]
                   for k in dict.fromkeys(k for m in ds.meta for k in m)
                   if k not in taken and len({m.get(k, "") for m in ds.meta}) > 1}
        return cls(load(), tab._groups(), tab.group_by.currentText(), columns, facets,
                   tab.facet_by.currentText() if facets else "",
                   [ds.label(c) for c in others] or None, load, tab)

    # --- what is asked ---------------------------------------------------------------

    def chosen_markers(self) -> list[str]:
        return [t.text() for t in self.marker_boxes if t.isChecked()]

    def set_markers(self, markers: list[str]) -> None:
        for tick in self.marker_boxes:
            tick.blockSignals(True)
            tick.setChecked(tick.text() in markers)
            tick.blockSignals(False)
        self.run()

    def _rescale(self) -> None:
        if self._reload is not None:
            self.medians = self._reload(self.scale.currentData())
        self.run()

    def run(self, *_) -> None:
        """Compute with what the boxes say now, and show it."""
        m = self.medians
        chosen = self.chosen_markers()
        self.markers_button.setText(f"{len(chosen)} of {len(m.markers)}")
        contrast = (self.first.currentData(), self.second.currentData())
        covariate = self.columns.get(self.adjust.currentData())
        keep = [i for i, name in enumerate(m.markers) if name in chosen]
        args = (m.values[:, :, keep], m.counts, m.clusters, chosen, self.groups, contrast)
        options = (covariate, self.weights.isChecked(), self.min_cells.value(),
                   self.min_samples.value(), self.correction.currentData())
        if self.first.count() < 2 or contrast[0] == contrast[1] or not chosen:
            self.results = []
        elif self.facets:
            self.results = expression.within(*args, self.facets, *options)
        else:
            self.results = expression.differential(*args, *options)
        self.result = expression.to_table(self.results, self.correction.currentData(),
                                          self.facet_name or None)
        self._fill()
        self.volcano.show_content(ExpressionVolcano(
            self.results, correction=self.correction.currentData()))
        self.heatmap.show_content(DifferenceHeatmap(self.results))
        self.legend.setText(self._legend(contrast))
        self.statement.setText(self._statement(contrast))

    def _legend(self, contrast: tuple[str, str]) -> str:
        content = self.heatmap.content()
        if content is None:
            return ""
        what = {compare.NONE: "p", compare.BH: "FDR"}.get(
            self.correction.currentData(), "adjusted p")
        return (f"Red: higher in {contrast[0] or '—'}; blue: higher in {contrast[1] or '—'}; "
                f"the deepest colour is a difference of {content.limit():.3g}. Grey: not "
                f"tested. * {what} < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001.")

    def _statement(self, contrast: tuple[str, str]) -> str:
        if not self.results:
            return ("Choose two different groups to compare." if self.chosen_markers()
                    else "Choose at least one marker.")
        correction = self.correction.currentData()
        return METHOD.format(
            scale=SCALES.get(self.medians.scale, "on the scale of the analysis"),
            first=contrast[0] or "—", second=contrast[1] or "—",
            adjusting=(f" with {self.adjust.currentData()} as a second factor"
                       if self.adjust.currentData() else ""),
            weighing=(", each median weighing as the cells it is of,"
                      if self.weights.isChecked() else ""),
            within=(f", separately within each {self.facet_name}" if self.facets else ""),
            cells=self.min_cells.value(), samples=self.min_samples.value(),
            tests=sum(1 for r in self.results if r.p == r.p),
            correction=("p not adjusted for their number" if correction == compare.NONE
                        else f"p adjusted by {correction} over all of them"))

    def _fill(self) -> None:
        table = self.result
        self.table.setSortingEnabled(False)
        self.table.clear()
        self.table.setRowCount(len(table.rows))
        self.table.setColumnCount(len(table.columns))
        self.table.setHorizontalHeaderLabels(table.columns)
        for r, row in enumerate(table.rows):
            for k, value in enumerate(row):
                self.table.setItem(r, k, _Cell(value, _text(table.columns[k], value)))
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()

    # --- taking it elsewhere ---------------------------------------------------------

    def make_figure(self) -> None:
        from .figure_dialog import open_for
        if self.views.currentIndex() == 0:
            open_for(self.heatmap, self, "expression_differences")
        else:
            open_for(self.volcano, self, "expression_volcano")

    def save_csv(self, path: str) -> None:
        self.result.to_csv(path)

    def _ask_export(self) -> None:
        folder = settings().value(SETTING_FOLDER, "") or os.path.expanduser("~")
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export", os.path.join(folder, "differential_expression.csv"),
            "Comma-separated values (*.csv)")
        if path:
            self.save_csv(path)
            settings().setValue(SETTING_FOLDER, os.path.dirname(path))

    def copy(self) -> None:
        QtWidgets.QApplication.clipboard().setText(self.result.to_text())
