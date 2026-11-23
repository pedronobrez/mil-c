"""
The statistics table: choose populations, parameters and statistics, and
get a row for each population of each sample — with the samples' own
columns beside the name, so the table can be sorted here and compared by
group anywhere else.
"""

from __future__ import annotations

import os

import numpy as np
from PyQt6 import QtCore, QtWidgets

from .. import stats
from ..session import Session
from . import style, worker
from .settings import settings

SETTING_FOLDER = "statistics/folder"
SETTING_STATS = "statistics/chosen"


class _Cell(QtWidgets.QTableWidgetItem):
    """Shows a number as it should be read and sorts it as a number."""

    def __init__(self, value, text: str):
        super().__init__(text)
        self.value = value
        self.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
        if isinstance(value, (int, float, np.integer, np.floating)):
            self.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight
                                  | QtCore.Qt.AlignmentFlag.AlignVCenter)

    def __lt__(self, other):
        a, b = self.value, getattr(other, "value", None)
        numeric = (int, float, np.integer, np.floating)
        if isinstance(a, numeric) and isinstance(b, numeric):
            # no value sorts last whichever way the column is turned
            a = float("-inf") if a != a else a
            b = float("-inf") if b != b else b
            return a < b
        return str(a) < str(b)


def _format(column: str, value) -> str:
    if value is None or (isinstance(value, float) and value != value):
        return ""
    if column == "Count":
        return f"{int(value):,d}"
    if column.startswith("% of") or column.endswith("CV"):
        return f"{value:.2f}"
    if isinstance(value, (float, np.floating)):
        return f"{value:,.1f}" if abs(value) >= 100 else f"{value:.3g}"
    return str(value)


class StatisticsTab(QtWidgets.QWidget):
    sigStatus = QtCore.pyqtSignal(str)

    def __init__(self, session: Session, parent=None):
        super().__init__(parent)
        self.session = session
        self.result: stats.StatsTable | None = None
        self._busy = False

        self.populations = QtWidgets.QTreeWidget()
        self.populations.setHeaderHidden(True)
        self.populations.setIndentation(11)
        self.parameters = QtWidgets.QListWidget()
        self.parameters.setToolTip(
            "The parameters to report statistics of. Fluorescence is taken "
            "compensated and in the file's own units, as FlowJo reports it.")
        self.statistics = QtWidgets.QListWidget()
        self.statistics.setToolTip(
            "Robust SD is half the distance between the 15.87th and 84.13th "
            "percentiles; Robust CV is that over the median. The geometric mean "
            "leaves out events at or below zero, which have no logarithm.")
        remembered = settings().value(SETTING_STATS, "Median") or "Median"
        for name in stats.STATISTICS:
            item = QtWidgets.QListWidgetItem(name)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.CheckState.Checked if name in str(remembered).split("|")
                               else QtCore.Qt.CheckState.Unchecked)
            self.statistics.addItem(item)
        self.only_ticked = QtWidgets.QCheckBox("Only the samples ticked in Samples")
        self.only_ticked.setChecked(True)
        self.compute_button = QtWidgets.QPushButton("Compute")
        self.compute_button.setProperty("primary", True)
        self.export_button = QtWidgets.QPushButton("Export CSV…")
        self.copy_button = QtWidgets.QPushButton("Copy")
        self.copy_button.setToolTip("The whole table, to paste into Excel or Prism")

        left = QtWidgets.QVBoxLayout()
        left.addWidget(self._eyebrow("Populations"))
        left.addWidget(self.populations, 3)
        left.addWidget(self._eyebrow("Parameters"))
        left.addWidget(self.parameters, 3)
        left.addWidget(self._eyebrow("Statistics"))
        left.addWidget(self.statistics, 2)
        left.addWidget(self.only_ticked)
        left.addWidget(self.compute_button)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.export_button)
        row.addWidget(self.copy_button)
        left.addLayout(row)
        panel = QtWidgets.QWidget()
        panel.setLayout(left)
        panel.setFixedWidth(330)

        self.caption = QtWidgets.QLabel("")
        self.caption.setProperty("role", "caption")
        self.stale = QtWidgets.QLabel("")
        self.stale.setProperty("role", "warning")
        self.table = QtWidgets.QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.verticalHeader().setVisible(False)
        self.empty = QtWidgets.QLabel(
            "Tick populations, parameters and statistics, then Compute.\n"
            "Counts and percentages are always there.")
        self.empty.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.empty.setProperty("role", "hint")
        self.stack = QtWidgets.QStackedWidget()
        self.stack.addWidget(self.empty)
        self.stack.addWidget(self.table)
        head = QtWidgets.QHBoxLayout()
        head.addWidget(self.caption, 1)
        head.addWidget(self.stale)
        right = QtWidgets.QVBoxLayout()
        right.addLayout(head)
        right.addWidget(self.stack, 1)

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.addWidget(panel)
        layout.addLayout(right, 1)

        self.compute_button.clicked.connect(self.compute)
        self.export_button.clicked.connect(self.export)
        self.copy_button.clicked.connect(self.copy)
        self.refresh()

    @staticmethod
    def _eyebrow(text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setFont(style.eyebrow_font())
        label.setProperty("role", "caption")
        label.setContentsMargins(0, 8, 0, 2)
        return label

    # --- what can be asked for -------------------------------------------------

    def refresh(self) -> None:
        """After a folder was opened or the gates changed: the lists are
        rebuilt, keeping what was ticked where it still exists."""
        s = self.session
        ticked_pops = set(self.chosen_populations())
        known = set(getattr(self, "_known", ()))
        self._known = set()
        ticked_params = {c for _, c in self.chosen_parameters()}
        self.populations.clear()
        self.parameters.clear()
        entries = s.usable()
        has = bool(entries)
        for w in (self.compute_button, self.export_button, self.copy_button):
            w.setEnabled(has and (w is self.compute_button or self.result is not None))
        if not has:
            self.result = None
            self.stack.setCurrentIndex(0)
            self.caption.setText("")
            self.stale.setText("")
            return
        common = set(s.populations())
        checked, unchecked = QtCore.Qt.CheckState.Checked, QtCore.Qt.CheckState.Unchecked

        def add(node, parent_item):
            item = QtWidgets.QTreeWidgetItem(
                [node.name if node.gate is not None else "All events"])
            item.setData(0, QtCore.Qt.ItemDataRole.UserRole, node.path)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            # what was unticked stays unticked; a population that was not
            # there before — a gate just drawn — comes in ticked
            on = node.path in ticked_pops or node.path not in known
            self._known.add(node.path)
            item.setCheckState(0, checked if on else unchecked)
            (parent_item.addChild if parent_item is not None
             else self.populations.addTopLevelItem)(item)
            for child in node.children:
                # only what every gated sample has: a row per sample needs it
                if child.path in common:
                    add(child, item)

        gated = next((e for e in entries if e.gated), entries[0])
        add(gated.gates, None)
        self.populations.expandAll()
        for label, column in s.stat_channels():
            item = QtWidgets.QListWidgetItem(label)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, column)
            item.setToolTip(column)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(checked if column in ticked_params else unchecked)
            self.parameters.addItem(item)
        if self.result is not None:
            self.mark_stale()

    def mark_stale(self) -> None:
        """The gates, the samples or their columns changed since the table
        was computed: it still shows what it was computed from."""
        if self.result is not None:
            self.stale.setText("Computed before the last change — Compute again")

    def chosen_populations(self) -> list[str]:
        out = []
        it = QtWidgets.QTreeWidgetItemIterator(self.populations)
        while it.value() is not None:
            if it.value().checkState(0) == QtCore.Qt.CheckState.Checked:
                out.append(it.value().data(0, QtCore.Qt.ItemDataRole.UserRole))
            it += 1
        return out

    def chosen_parameters(self) -> list[tuple[str, str]]:
        return [(self.parameters.item(i).text(),
                 self.parameters.item(i).data(QtCore.Qt.ItemDataRole.UserRole))
                for i in range(self.parameters.count())
                if self.parameters.item(i).checkState() == QtCore.Qt.CheckState.Checked]

    def chosen_statistics(self) -> list[str]:
        return [self.statistics.item(i).text() for i in range(self.statistics.count())
                if self.statistics.item(i).checkState() == QtCore.Qt.CheckState.Checked]

    def tick(self, populations=None, parameters=None, statistics=None) -> None:
        """Tick exactly these (parameters by label or column name); None
        leaves a list as it is."""
        on, off = QtCore.Qt.CheckState.Checked, QtCore.Qt.CheckState.Unchecked
        if populations is not None:
            it = QtWidgets.QTreeWidgetItemIterator(self.populations)
            while it.value() is not None:
                path = it.value().data(0, QtCore.Qt.ItemDataRole.UserRole)
                it.value().setCheckState(0, on if path in populations else off)
                it += 1
        if parameters is not None:
            for i in range(self.parameters.count()):
                item = self.parameters.item(i)
                hit = item.text() in parameters or item.data(
                    QtCore.Qt.ItemDataRole.UserRole) in parameters
                item.setCheckState(on if hit else off)
        if statistics is not None:
            for i in range(self.statistics.count()):
                item = self.statistics.item(i)
                item.setCheckState(on if item.text() in statistics else off)

    # --- computing -----------------------------------------------------------------

    def compute(self) -> None:
        if self._busy or not self.session.usable():
            return
        populations = self.chosen_populations()
        if not populations:
            self.sigStatus.emit("Tick at least one population")
            return
        parameters, chosen = self.chosen_parameters(), self.chosen_statistics()
        settings().setValue(SETTING_STATS, "|".join(chosen))
        only = self.only_ticked.isChecked()
        self._busy = True
        self.compute_button.setEnabled(False)
        self.sigStatus.emit("Computing statistics…")

        def work(progress):
            return self.session.statistics(populations, parameters, chosen, progress,
                                           only_included=only)

        def done(table):
            self._busy = False
            self.compute_button.setEnabled(True)
            self.show_table(table)
            self.sigStatus.emit(f"{len(table.rows):,d} rows")

        def failed(message):
            self._busy = False
            self.compute_button.setEnabled(True)
            self.sigStatus.emit(message)

        worker.run(work, on_done=done, on_failed=failed,
                   on_progress=lambda d, t, n: self.sigStatus.emit(f"Reading {n} ({d}/{t})"))

    def show_table(self, table: stats.StatsTable) -> None:
        self.result = table
        self.stale.setText("")
        self.table.setSortingEnabled(False)
        # the path only matters where two populations share a name
        names = table.column("Population")
        ambiguous = len(set(names)) < len(set(zip(names, table.column("Path"))))
        shown = [i for i, c in enumerate(table.columns) if c != "Path" or ambiguous]
        self.table.setColumnCount(len(shown))
        self.table.setHorizontalHeaderLabels([table.columns[i] for i in shown])
        self.table.setRowCount(len(table.rows))
        for r, row in enumerate(table.rows):
            for k, i in enumerate(shown):
                self.table.setItem(r, k, _Cell(row[i], _format(table.columns[i], row[i])))
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        self.stack.setCurrentIndex(1)
        samples = len(set(table.column("Sample")))
        self.caption.setText(
            f"{len(table.rows):,d} rows: {samples} samples × "
            f"{len(table.rows) // max(samples, 1)} populations. Statistics are of "
            "compensated values in the file's own units. Click a header to sort.")
        self.export_button.setEnabled(True)
        self.copy_button.setEnabled(True)

    # --- taking it elsewhere -----------------------------------------------------------

    def export(self) -> None:
        if self.result is None:
            return
        folder = settings().value(SETTING_FOLDER, "") or self.session.folder
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export the statistics", os.path.join(folder, "statistics.csv"),
            "Comma-separated values (*.csv)")
        if path:
            self.result.to_csv(path)
            settings().setValue(SETTING_FOLDER, os.path.dirname(path))
            self.sigStatus.emit(f"Wrote {path}")

    def copy(self) -> None:
        if self.result is not None:
            QtWidgets.QApplication.clipboard().setText(self.result.to_text())
            self.sigStatus.emit(f"{len(self.result.rows):,d} rows copied")
