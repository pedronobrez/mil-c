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

from .. import compare, stats
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
        self.compare_button = QtWidgets.QPushButton("Compare groups…")
        self.compare_button.setToolTip(
            "Test any column of the table between the groups of samples: a t test or "
            "Mann-Whitney for two, ANOVA or Kruskal-Wallis for more, adjusted for the "
            "number of populations tested")

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
        left.addWidget(self.compare_button)
        self.prism_button = QtWidgets.QPushButton("Prism…")
        self.prism_button.setToolTip(
            "One column of the table as a GraphPad Prism file: populations down, groups "
            "across, a subcolumn for each sample")
        left.addWidget(self.prism_button)
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
        self.compare_button.clicked.connect(self._ask_comparison)
        self.prism_button.clicked.connect(self._ask_prism)
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
        for w in (self.compute_button, self.export_button, self.copy_button,
                  self.compare_button, self.prism_button):
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
        self.compare_button.setEnabled(True)
        self.prism_button.setEnabled(True)

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

    # --- between groups -----------------------------------------------------------------

    def value_columns(self) -> list[str]:
        """The columns of the table that hold a number per sample."""
        cols = self.result.columns
        return cols[cols.index("Count"):] if self.result is not None else []

    def comparison(self, value: str, group: list[str], within: str | None = None,
                   kind: str = compare.PARAMETRIC, correction: str = compare.BH,
                   pairwise: bool = False, paired_by: str | None = None):
        """(table of results, the comparisons): `value` between the groups
        `group` makes, for every population of the table. With `paired_by`,
        a column that says which animal or donor each sample is of, the
        groups are compared within each of them — and a ValueError says
        why where the column does not pair the samples."""
        if paired_by:
            return self._paired_comparison(value, group, within, kind, correction, pairwise,
                                           paired_by)
        values, names = compare.from_table(self.result, value, group, within)
        found = compare.compare(values, kind, correction, pairwise)
        return compare.to_table(found, names, correction), found

    def _paired_comparison(self, value, group, within, kind, correction, pairwise, paired_by):
        values, names, subjects = compare.from_table(self.result, value, group, within, paired_by)
        problem = compare.pairing_problem(values, subjects)
        if problem:
            raise ValueError(f"The samples cannot be paired by {paired_by}: {problem}.")
        found = compare.compare(values, kind, correction, pairwise, subjects)
        return compare.to_table(found, names, correction, paired_by), found

    def save_prism(self, path: str, value: str, group: list[str],
                   within: str | None = None) -> int:
        from .. import prism
        tables = prism.from_stats_table(self.result, value, group, within)
        prism.write_pzfx(path, tables)
        self.sigStatus.emit(f"Wrote {len(tables)} Prism table(s) of {value} to {path}")
        return len(tables)

    def _ask_prism(self) -> None:
        if self.result is None:
            return
        dialog = CompareDialog(self.value_columns(), self.session.columns,
                               self.session.group_by, self, for_prism=True)
        if not dialog.exec() or not dialog.groups():
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save for GraphPad Prism",
            os.path.join(self.session.folder, "statistics.pzfx"), "GraphPad Prism (*.pzfx)")
        if path:
            self.save_prism(path if path.endswith(".pzfx") else path + ".pzfx",
                            dialog.value.currentData(), dialog.groups(),
                            dialog.within.currentData() or None)

    def _ask_comparison(self) -> None:
        if self.result is None:
            return
        from .results_dialog import NOTE, ResultsDialog
        dialog = CompareDialog(self.value_columns(), self.session.columns,
                               self.session.group_by, self)
        if not dialog.exec() or not dialog.groups():
            return
        try:
            table, found = self.comparison(
                dialog.value.currentData(), dialog.groups(), dialog.within.currentData() or None,
                dialog.test.currentData(), dialog.correction.currentData(),
                dialog.pairwise.isChecked(), dialog.paired_by.currentData() or None)
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, "Compare groups", str(exc))
            return
        note = NOTE.format(tests=sum(1 for c in found if c.p == c.p),
                           correction=dialog.correction.currentData())
        ResultsDialog(table, f"{dialog.value.currentData()} between groups",
                      f"{dialog.value.currentData()} of each population, compared between "
                      f"{' + '.join(dialog.groups())}. " + note, self,
                      "group_comparison").exec()

    def copy(self) -> None:
        if self.result is not None:
            QtWidgets.QApplication.clipboard().setText(self.result.to_text())
            self.sigStatus.emit(f"{len(self.result.rows):,d} rows copied")


class CompareDialog(QtWidgets.QDialog):
    """What to compare, between which groups, and how."""

    def __init__(self, values: list[str], columns: list[str], chosen: list[str], parent=None,
                 for_prism: bool = False):
        super().__init__(parent)
        self.setWindowTitle("A table for Prism" if for_prism else "Compare groups")
        self.value = QtWidgets.QComboBox()
        for v in values:
            self.value.addItem(v, v)
        self.value.setCurrentIndex(max(self.value.findData("% of parent"), 0))
        self.group_list = QtWidgets.QListWidget()
        self.group_list.setMaximumHeight(110)
        for c in columns:
            item = QtWidgets.QListWidgetItem(c)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.CheckState.Checked if c in chosen[:1]
                               else QtCore.Qt.CheckState.Unchecked)
            self.group_list.addItem(item)
        self.within = QtWidgets.QComboBox()
        self.within.addItem("—", "")
        for c in columns:
            self.within.addItem(c, c)
        if len(chosen) > 1:
            self.within.setCurrentIndex(max(self.within.findData(chosen[1]), 0))
        self.within.setToolTip("A separate comparison for each value of this column: "
                               "the groups within each day, say")
        self.test = QtWidgets.QComboBox()
        for key, text in compare.KINDS.items():
            self.test.addItem(text, key)
        self.paired_by = QtWidgets.QComboBox()
        self.paired_by.addItem("—", "")
        for c in columns:
            self.paired_by.addItem(c, c)
        self.paired_by.setToolTip(
            "The column that says which animal or donor a sample is of, where the same ones "
            "were measured in every group. The test is then of each one's change: paired t "
            "or Wilcoxon for two groups, repeated-measures ANOVA or Friedman for more. "
            "Those not found in every group are left out.")
        self.correction = QtWidgets.QComboBox()
        for name in compare.CORRECTIONS:
            self.correction.addItem(name, name)
        self.pairwise = QtWidgets.QCheckBox("With more than two groups, every pair as well")
        self.pairwise.setChecked(True)
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                         | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        box.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setText("Compare")
        box.accepted.connect(self.accept)
        box.rejected.connect(self.reject)
        form = QtWidgets.QFormLayout(self)
        form.addRow("Compare", self.value)
        form.addRow("Between", self.group_list)
        form.addRow("Within each", self.within)
        if for_prism:   # Prism does the testing itself
            box.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setText("Save…")
            form.labelForField(self.value).setText("Values")
            form.labelForField(self.group_list).setText("Columns")
            form.labelForField(self.within).setText("A table for each")
        else:
            form.addRow("Paired by", self.paired_by)
            form.addRow("Test", self.test)
            form.addRow("Adjust p by", self.correction)
            form.addRow(self.pairwise)
        form.addRow(box)

    def groups(self) -> list[str]:
        return [self.group_list.item(i).text() for i in range(self.group_list.count())
                if self.group_list.item(i).checkState() == QtCore.Qt.CheckState.Checked
                and self.group_list.item(i).text() != self.within.currentData()]
