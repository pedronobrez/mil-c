"""
The first tab: which files, how they are grouped, how many events each
gives. The plan is on the table before anything is concatenated — what
each sample has, what it will give, and which one set the limit.
"""

from __future__ import annotations

import os
import re

from PyQt6 import QtCore, QtGui, QtWidgets

from .. import sampling
from ..session import ALL_EVENTS, Session
from . import style, worker

MODES = [
    (sampling.EQUAL_PER_SAMPLE, "Same number from every sample"),
    (sampling.EQUAL_PER_GROUP, "Same total from every group"),
    (sampling.CEILING, "At most N from each sample"),
]

FIXED = ["", "Sample"]            # before the annotation columns
TAIL = ["Events", "In gate", "Takes", "Note"]


class SamplesTab(QtWidgets.QWidget):
    #: a dataset was built
    sigConcatenated = QtCore.pyqtSignal()
    sigStatus = QtCore.pyqtSignal(str)

    def __init__(self, session: Session, parent=None):
        super().__init__(parent)
        self.session = session
        self._filling = False
        self._busy = False

        self.folder_label = QtWidgets.QLabel("No folder open")
        self.folder_label.setProperty("role", "strong")
        self.wsp_label = QtWidgets.QLabel("")
        self.wsp_label.setProperty("role", "caption")
        self.open_button = QtWidgets.QPushButton("Open folder…")
        self.wsp_button = QtWidgets.QPushButton("Workspace…")
        self.wsp_button.setToolTip("A FlowJo .wsp: its gates and its "
                                   "compensation are applied to the files it names")
        top = QtWidgets.QHBoxLayout()
        names = QtWidgets.QVBoxLayout()
        names.setSpacing(2)
        names.addWidget(self.folder_label)
        names.addWidget(self.wsp_label)
        top.addLayout(names, 1)
        top.addWidget(self.wsp_button)
        top.addWidget(self.open_button)

        self.table = QtWidgets.QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemChanged.connect(self._ticked)

        # --- the plan's settings
        self.population = QtWidgets.QComboBox()
        self.population.setSizeAdjustPolicy(
            QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.population.setMinimumContentsLength(24)
        self.extend = QtWidgets.QCheckBox("Gate the files the workspace does not list")
        self.extend.setToolTip(
            "Every sample of the workspace has the same gates and the same "
            "compensation matrix. Ticked, the files it does not list are given "
            "those too — what FlowJo does when a file is dropped into the group. "
            "Nobody has looked at those gates on those files: check them.")
        self.pattern = QtWidgets.QLineEdit(session.pattern)
        self.pattern.setToolTip(
            "A regular expression read from each file name. Every named "
            "group, (?P<name>…), becomes a column the samples can be "
            "grouped and coloured by.")
        self.group_by = QtWidgets.QLineEdit(",".join(session.group_by))
        self.group_by.setToolTip("Annotation(s) that make a group, "
                                 "comma-separated: group, or day,group")
        self.mode = QtWidgets.QComboBox()
        for key, text in MODES:
            self.mode.addItem(text, key)
        self.n = QtWidgets.QSpinBox()
        self.n.setRange(0, 10_000_000)
        self.n.setSpecialValueText("as many as the data allows")
        self.n.setGroupSeparatorShown(True)
        self.min_events = QtWidgets.QSpinBox()
        self.min_events.setRange(1, 10_000_000)
        self.min_events.setGroupSeparatorShown(True)
        self.min_events.setToolTip(
            "A sample with fewer events than this is left out instead of "
            "pulling every other sample down to its size")
        self.seed = QtWidgets.QSpinBox()
        self.seed.setRange(0, 999_999)
        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)
        self.notes = QtWidgets.QLabel("")
        self.notes.setWordWrap(True)
        self.notes.setProperty("role", "warning")
        self.build_button = QtWidgets.QPushButton("Concatenate")
        self.build_button.setProperty("primary", True)

        form = QtWidgets.QFormLayout()
        form.setLabelAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        form.addRow(self._eyebrow("Samples"))
        form.addRow("Population", self.population)
        form.addRow(self.extend)
        form.addRow("Name pattern", self.pattern)
        form.addRow("Group by", self.group_by)
        form.addRow(self._eyebrow("Downsampling"))
        form.addRow("Balance", self.mode)
        form.addRow("N", self.n)
        form.addRow("Minimum", self.min_events)
        form.addRow("Seed", self.seed)
        side = QtWidgets.QVBoxLayout()
        side.addLayout(form)
        side.addSpacing(8)
        side.addWidget(self.summary)
        side.addWidget(self.notes)
        side.addStretch(1)
        side.addWidget(self.build_button)
        panel = QtWidgets.QWidget()
        panel.setLayout(side)
        panel.setFixedWidth(360)

        body = QtWidgets.QHBoxLayout()
        body.addWidget(self.table, 1)
        body.addWidget(panel)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.addLayout(top)
        layout.addLayout(body, 1)

        self.population.currentIndexChanged.connect(self._population_changed)
        self.extend.toggled.connect(self._extend_changed)
        self.pattern.editingFinished.connect(self._pattern_changed)
        self.group_by.editingFinished.connect(self._settings_changed)
        self.mode.currentIndexChanged.connect(self._settings_changed)
        self.n.valueChanged.connect(self._settings_changed)
        self.min_events.valueChanged.connect(self._settings_changed)
        self.seed.valueChanged.connect(self._settings_changed)
        self.build_button.clicked.connect(self.concatenate)
        self.refresh()

    @staticmethod
    def _eyebrow(text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setFont(style.eyebrow_font())
        label.setProperty("role", "caption")
        label.setContentsMargins(0, 10, 0, 2)
        return label

    # --- from the session ----------------------------------------------------

    def refresh(self) -> None:
        """Everything, after the folder or the workspace changed."""
        s = self.session
        self._filling = True
        self.folder_label.setText(s.folder or "No folder open")
        self.wsp_label.setText(
            f"Gates and compensation from {os.path.basename(s.workspace.path)}"
            if s.workspace else
            ("No workspace: every event of every file, compensated with the "
             "file's own matrix" if s.folder else
             "Open a folder of .fcs files; a FlowJo workspace beside them is picked up"))
        self.population.clear()
        self.population.addItem("All events", ALL_EVENTS)
        for path in s.populations():
            depth = path.count("/")
            self.population.addItem("    " * depth + path.split("/")[-1], path)
            self.population.setItemData(self.population.count() - 1, path,
                                        QtCore.Qt.ItemDataRole.ToolTipRole)
        i = self.population.findData(s.population)
        self.population.setCurrentIndex(max(i, 0))
        self.extend.setVisible(s.can_extend()
                               and any(e.wsp is None and not e.problem for e in s.entries))
        self.extend.setChecked(s.extend_gates)
        self.pattern.setText(s.pattern)
        self._filling = False
        self._fill_table()
        self._count()

    def _columns(self) -> list[str]:
        return FIXED + self.session.meta_keys() + TAIL

    def _fill_table(self) -> None:
        s = self.session
        self._filling = True
        cols = self._columns()
        self.table.setColumnCount(len(cols))
        self.table.setHorizontalHeaderLabels([c.capitalize() if c in s.meta_keys() else c
                                              for c in cols])
        self.table.setRowCount(len(s.entries))
        for r, e in enumerate(s.entries):
            tick = QtWidgets.QTableWidgetItem()
            flags = QtCore.Qt.ItemFlag.ItemIsUserCheckable | QtCore.Qt.ItemFlag.ItemIsEnabled
            tick.setFlags(flags if not e.problem else QtCore.Qt.ItemFlag.NoItemFlags)
            tick.setCheckState(QtCore.Qt.CheckState.Checked if e.include and not e.problem
                               else QtCore.Qt.CheckState.Unchecked)
            self.table.setItem(r, 0, tick)
            self.table.setItem(r, 1, QtWidgets.QTableWidgetItem(e.name))
            for k in s.meta_keys():
                self.table.setItem(r, cols.index(k),
                                   QtWidgets.QTableWidgetItem(e.meta.get(k, "")))
            self._number(r, cols.index("Events"), e.events)
        self._filling = False
        self._update_plan()
        self.table.resizeColumnsToContents()
        self.table.setColumnWidth(0, 30)

    def _number(self, row: int, col: int, value, muted: bool = False) -> None:
        item = QtWidgets.QTableWidgetItem("" if value is None else f"{value:,d}")
        item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight
                              | QtCore.Qt.AlignmentFlag.AlignVCenter)
        if muted:
            item.setForeground(QtGui.QBrush(QtGui.QColor(style.token("ink_faint"))))
        self.table.setItem(row, col, item)

    def _update_plan(self) -> None:
        """The In gate, Takes and Note columns, and the summary."""
        s = self.session
        cols = self._columns()
        c_gate, c_take, c_note = (cols.index(x) for x in ("In gate", "Takes", "Note"))
        plan = s.plan()
        self._filling = True
        for r, e in enumerate(s.entries):
            counted = s.available(e) and s.counted(e)
            self._number(r, c_gate, s.count(e) if counted else None)
            take = plan.take.get(e.name)
            self._number(r, c_take, take if take else None)
            if e.problem:
                note = e.problem
            elif not s.available(e):
                note = "no gates in the workspace"
            elif not counted:
                note = "counting…"
            elif e.name in plan.excluded:
                note = plan.excluded[e.name]
            elif plan.limited_by == e.name:
                note = "sets the limit"
            elif e.borrowed and s.population != ALL_EVENTS:
                note = "gates from the group"
            else:
                note = ""
            item = QtWidgets.QTableWidgetItem(note)
            quiet = note in ("", "counting…", "gates from the group", "sets the limit")
            item.setForeground(QtGui.QBrush(QtGui.QColor(style.token(
                "ink_muted" if quiet else "warning"))))
            self.table.setItem(r, c_note, item)
        self._filling = False

        if not s.entries:
            self.summary.setText("")
            self.notes.setText("")
            self.build_button.setEnabled(False)
            return
        included = [x for x in plan.take.values() if x > 0]
        text = (f"<b>{plan.total:,d} events</b> from {len(included)} samples"
                f" — {plan.used_fraction():.0%} of what they have.")
        totals = plan.group_totals()
        if len(totals) > 1:
            text += "<br>" + " · ".join(f"{g or '—'} {n:,d}" for g, n in sorted(totals.items()))
        if plan.limited_by:
            text += f"<br>Limited by <b>{plan.limited_by}</b>."
        self.summary.setText(text)
        notes = list(plan.notes)
        missing = [e.name for e in s.entries if e.include and not e.problem
                   and not s.available(e)]
        if missing:
            notes.append(f"{len(missing)} file(s) have no gates in the workspace "
                         "and are left out of this population.")
        for names in s.repeats():
            notes.append(" and ".join(names) + " look like the same animal acquired "
                         "twice: untick one, or it counts double.")
        self.notes.setText("\n".join(notes))
        self.build_button.setEnabled(plan.total > 0 and not self._busy)

    # --- from the person -------------------------------------------------------

    def _ticked(self, item) -> None:
        if self._filling or item.column() != 0:
            return
        self.session.entries[item.row()].include = (
            item.checkState() == QtCore.Qt.CheckState.Checked)
        self._update_plan()

    def _population_changed(self) -> None:
        if self._filling:
            return
        self.session.population = self.population.currentData() or ALL_EVENTS
        self._update_plan()
        self._count()

    def _extend_changed(self, on: bool) -> None:
        if self._filling:
            return
        self.session.set_extend_gates(on)
        self.refresh()

    def _pattern_changed(self) -> None:
        text = self.pattern.text()
        try:
            re.compile(text)
        except re.error as exc:
            self.sigStatus.emit(f"Name pattern: {exc}")
            return
        if text != self.session.pattern:
            self.session.annotate(text)
            self._fill_table()

    def _settings_changed(self) -> None:
        if self._filling:
            return
        s = self.session
        s.group_by = [k.strip() for k in self.group_by.text().split(",") if k.strip()]
        s.mode = self.mode.currentData()
        s.n = self.n.value() or None
        s.min_events = self.min_events.value()
        s.seed = self.seed.value()
        self._update_plan()

    # --- work ------------------------------------------------------------------

    def _count(self) -> None:
        s = self.session
        if s.population == ALL_EVENTS or not s.entries:
            return
        self.sigStatus.emit("Counting the gate in every sample…")

        def progress(done, total, name):
            self.sigStatus.emit(f"Counting {name} ({done}/{total})")
            self._update_plan()

        worker.run(s.count_all, on_done=lambda _: (self._update_plan(),
                                                   self.sigStatus.emit("")),
                   on_failed=self.sigStatus.emit, on_progress=progress)

    def concatenate(self) -> None:
        if self._busy:
            return
        self._busy = True
        self.build_button.setEnabled(False)
        self.sigStatus.emit("Concatenating…")

        def done(_):
            self._busy = False
            self._update_plan()
            ds = self.session.dataset
            self.sigStatus.emit(f"Concatenated {len(ds):,d} events from "
                                f"{len(ds.samples)} samples")
            self.sigConcatenated.emit()

        def failed(message):
            self._busy = False
            self._update_plan()
            self.sigStatus.emit(message)

        worker.run(self.session.build, on_done=done, on_failed=failed,
                   on_progress=lambda d, t, n: self.sigStatus.emit(
                       f"Reading {n} ({d}/{t})"))
