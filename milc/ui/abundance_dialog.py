"""
Differential abundance of the clusters, asked from the Clusters tab: which
two groups, adjusted for what, and the answer three ways — a table to
sort, copy and save, a volcano, and a paragraph that says what was done
and what it does not show. The paragraph is part of the result: a p value
from a model is only as good as the reader's idea of the model.
"""

from __future__ import annotations

import os

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets

from .help_window import describe
from .. import abundance, cluster_stats, clusters, compare
from ..volcano import COLORS, REST, VolcanoContent
from . import style, theme
from .results_dialog import SETTING_FOLDER, _text
from .settings import settings
from .statistics_tab import _Cell

METHOD = (
    "Each cluster's number of cells in each sample is modelled as negative-binomial with the "
    "sample's total as offset, as diffcyt's DA-edgeR does, and the difference between "
    "{first} and {second} is tested by a {test}{adjusting}{within}. Each sample is one "
    "observation; {tests} test(s), p adjusted by Benjamini-Hochberg (FDR). "
    "The limits: the cells counted are those that were concatenated, after the sampling, not "
    "all that were acquired; what is tested is a cluster's share of the sample, so one "
    "cluster that grows makes every other a little smaller; and the method is written in "
    "the manner of edgeR, not checked against edgeR itself. In simulations with 5 samples a "
    "group the quasi-likelihood test gave p < 0.05 about as often as it should under no "
    "difference; with 2 or 3 a group, or with the likelihood-ratio test, somewhat too "
    "often, and very small p values more so.")
#: the tests as the sentence above names them
SAID = {abundance.QL: "quasi-likelihood F test",
        abundance.LRT: "likelihood-ratio test with cluster-wise dispersions"}


class VolcanoView(pg.GraphicsLayoutWidget):
    """Fold change across, significance up, a panel for each value of
    "within each"; the significant clusters coloured and named."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._content: VolcanoContent | None = None
        self.panels: list[pg.PlotItem] = []
        self.setBackground(theme.background())

    def content(self):
        return self._content if self._content and self._content.results else None

    def show_content(self, content: VolcanoContent | None) -> None:
        self._content = content
        self.clear()
        self.panels = []
        self.setBackground(theme.background())
        if self.content() is None:
            return
        half, top = content.limits()
        guide = pg.mkPen(theme.axis(), width=1, style=QtCore.Qt.PenStyle.DashLine)
        font = QtGui.QFont()
        font.setPointSize(9)
        for k, panel in enumerate(content.panels()):
            plot = self.addPlot(row=0, col=k)
            plot.setMenuEnabled(False)
            plot.hideButtons()
            plot.setMouseEnabled(False, False)
            theme.style_axes(plot)
            if panel is not None:
                plot.setTitle(str(panel) or "—", color=theme.foreground(), size="10pt")
            plot.addItem(pg.InfiniteLine(-np.log10(content.alpha), angle=0, pen=guide))
            for x in (-content.fold, content.fold):
                plot.addItem(pg.InfiniteLine(x, angle=90, pen=guide))
            pts = content.points(panel)
            plot.addItem(pg.ScatterPlotItem(
                [p[0] for p in pts], [p[1] for p in pts], size=9,
                brush=[pg.mkBrush(COLORS[p[3]]) for p in pts],
                pen=pg.mkPen(theme.background(), width=1)))
            for x, y, cluster, kind in pts:
                if kind != REST:
                    name = pg.TextItem(cluster, anchor=(0.0, 1.0), color=theme.foreground())
                    name.setFont(font)
                    name.setPos(x, y)
                    plot.addItem(name)
            plot.setLabel("bottom", content.xlabel)
            if k == 0:
                plot.setLabel("left", content.ylabel)
            plot.setXRange(-half, half, padding=0)
            plot.setYRange(0, top, padding=0)
            self.panels.append(plot)


class AbundanceDialog(QtWidgets.QDialog):
    """`counts` is clusters x samples; `groups` names each sample's group
    (`group_name` is the column they come from); `columns` are the other
    columns of the samples that may be adjusted for; `facets`, where
    given, are the values of the tab's "within each"."""

    def __init__(self, counts: cluster_stats.Table, groups: list[str], group_name: str = "Group",
                 columns: dict[str, list[str]] | None = None,
                 facets: list[str] | None = None, facet_name: str = "", parent=None):
        super().__init__(parent)
        describe(self, "abundance")
        self.setWindowTitle("Differential abundance")
        self.counts, self.groups, self.columns = counts, list(groups), dict(columns or {})
        self.facets, self.facet_name = facets, facet_name
        self.results: list[abundance.Result] = []
        self.result = abundance.to_table([])

        self.first, self.second = QtWidgets.QComboBox(), QtWidgets.QComboBox()
        levels = compare.ordered(groups)
        for combo in (self.first, self.second):
            for level in levels:
                combo.addItem(level or "—", level)
        self.second.setCurrentIndex(min(1, len(levels) - 1))
        self.second.setToolTip("The reference: a positive fold change is more in the first")
        self.adjust = QtWidgets.QComboBox()
        self.adjust.addItem("—", "")
        for name in self.columns:
            self.adjust.addItem(name, name)
        self.adjust.setToolTip(
            "A batch, a day, an animal's sex: its effect is estimated beside the group's "
            "and taken out, as an additive term of the model")
        self.test = QtWidgets.QComboBox()
        for key, text in abundance.TESTS.items():
            self.test.addItem(text, key)
        self.test.setToolTip(
            "The quasi-likelihood F test allows for the dispersion being an estimate from "
            "few samples, and is the better calibrated. The likelihood-ratio test with "
            "cluster-wise dispersions is what diffcyt calls; it finds more and is a little "
            "too generous with few samples.")
        bar = QtWidgets.QHBoxLayout()
        for label, widget in ((f"Compare {group_name}", self.first), ("against", self.second),
                              ("Adjust for", self.adjust), ("Test", self.test)):
            cap = QtWidgets.QLabel(label)
            cap.setProperty("role", "caption")
            bar.addWidget(cap)
            bar.addWidget(widget)
            bar.addSpacing(10)
        bar.addStretch(1)

        self.table = QtWidgets.QTableWidget()
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.volcano = VolcanoView()
        figure = QtWidgets.QToolButton()
        figure.setText("Figure…")
        figure.setAutoRaise(True)
        title = QtWidgets.QLabel("Volcano")
        title.setFont(style.eyebrow_font())
        title.setProperty("role", "caption")
        head = QtWidgets.QHBoxLayout()
        head.addWidget(title)
        head.addStretch(1)
        head.addWidget(figure)
        side = QtWidgets.QWidget()
        right = QtWidgets.QVBoxLayout(side)
        right.setContentsMargins(0, 0, 0, 0)
        right.addLayout(head)
        right.addWidget(self.volcano, 1)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        split.addWidget(self.table)
        split.addWidget(side)
        split.setChildrenCollapsible(False)
        split.setSizes([700, 500])

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
        layout.addWidget(split, 1)
        layout.addWidget(self.statement)
        layout.addLayout(buttons)
        self.resize(1240, 700)

        for combo in (self.first, self.second, self.adjust, self.test):
            combo.currentIndexChanged.connect(self.run)
        figure.clicked.connect(self.make_figure)
        export.clicked.connect(self._ask_export)
        copy.clicked.connect(self.copy)
        close.clicked.connect(self.accept)
        self.run()

    @classmethod
    def for_tab(cls, tab) -> "AbundanceDialog":
        """What the Clusters tab is comparing now, as counts."""
        ds = tab.session.dataset
        # the clusters as the tab shows them: merged, and by their names
        now = tab.session.clusters(tab.key)
        counts = clusters.named(
            cluster_stats.counts(now.labels, ds.sample_index, ds.samples), now)
        facets = tab._facets()
        used = {tab.group_by.currentData(), tab.facet_by.currentData() if facets else ""}
        columns = {k: [m.get(k, "") for m in ds.meta]
                   for k in dict.fromkeys(k for m in ds.meta for k in m)
                   if k not in used and len({m.get(k, "") for m in ds.meta}) > 1}
        return cls(counts, tab._groups(), tab.group_by.currentText(), columns, facets,
                   tab.facet_by.currentText() if facets else "", tab)

    def run(self, *_) -> None:
        """Compute with what the boxes say now, and show it."""
        contrast = (self.first.currentData(), self.second.currentData())
        covariate = self.columns.get(self.adjust.currentData())
        args = (self.counts.values, self.counts.rows, self.groups, contrast)
        if self.first.count() < 2 or contrast[0] == contrast[1]:
            self.results = []
        elif self.facets:
            self.results = abundance.within(*args, self.facets, covariate,
                                            self.test.currentData())
        else:
            self.results = abundance.differential(*args, covariate, self.test.currentData())
        self.result = abundance.to_table(self.results, panel=self.facet_name or None)
        self._fill()
        self.volcano.show_content(VolcanoContent(self.results))
        self.statement.setText(self._statement(contrast))

    def _statement(self, contrast: tuple[str, str]) -> str:
        if not self.results:
            return "Choose two different groups to compare."
        return METHOD.format(
            first=contrast[0] or "—", second=contrast[1] or "—",
            test=SAID[self.test.currentData()],
            adjusting=(f", adjusted for {self.adjust.currentData()}"
                       if self.adjust.currentData() else ""),
            within=(f", separately within each {self.facet_name} (model, dispersion and all)"
                    if self.facets else ""),
            tests=sum(1 for r in self.results if r.p == r.p))

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

    def make_figure(self) -> None:
        from .figure_dialog import open_for
        open_for(self.volcano, self, "volcano")

    def save_csv(self, path: str) -> None:
        self.result.to_csv(path)

    def _ask_export(self) -> None:
        folder = settings().value(SETTING_FOLDER, "") or os.path.expanduser("~")
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export", os.path.join(folder, "differential_abundance.csv"),
            "Comma-separated values (*.csv)")
        if path:
            self.save_csv(path)
            settings().setValue(SETTING_FOLDER, os.path.dirname(path))

    def copy(self) -> None:
        QtWidgets.QApplication.clipboard().setText(self.result.to_text())
