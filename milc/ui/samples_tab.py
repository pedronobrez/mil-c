"""
The first tab: which files, how they are grouped, how many events each
gives. The plan is on the table before anything is concatenated — what
each sample has, what it will give, and which one set the limit.
"""

from __future__ import annotations

import os
import re

from PyQt6 import QtCore, QtGui, QtWidgets

from .. import order, sampling
from ..session import ALL_EVENTS, Session, suggest_pattern
from . import style, worker
from .flow import ElidedLabel, flow, together

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
    #: a column was added, removed, renamed or chosen to group by
    sigColumnsChanged = QtCore.pyqtSignal()
    #: the order of a column's values was changed
    sigOrderChanged = QtCore.pyqtSignal()
    #: the workspace's gates were extended to the other files, or taken back
    sigGatesChanged = QtCore.pyqtSignal()
    #: the values the other tabs read may have changed under them
    sigNormalisationChanged = QtCore.pyqtSignal()

    def __init__(self, session: Session, parent=None):
        super().__init__(parent)
        self.session = session
        self._filling = False
        self._busy = False

        self.folder_label = ElidedLabel("No folder open")
        self.folder_label.setProperty("role", "strong")
        self.wsp_label = ElidedLabel("")
        self.wsp_label.setProperty("role", "caption")
        self.add_column_button = QtWidgets.QPushButton("Add column…")
        self.add_column_button.setToolTip(
            "A column of your own — day, tissue, donor, treatment — to balance, "
            "colour and split by. Type in its cells; with several rows selected, "
            "what you type goes into all of them.")
        self.order_button = QtWidgets.QPushButton("Order…")
        self.order_button.setToolTip("Which group comes first — Ctrl before Case, day 3 "
                                     "before day 14 — on every chart, legend and table")
        self.order_button.clicked.connect(self.ask_order)
        self.fill_button = QtWidgets.QPushButton("Fill from file names…")
        self.fill_button.setToolTip("Read columns out of the file names with a pattern")
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
        # in a row that wraps: on Windows its five buttons on one line were
        # wider than the window may be
        top.addWidget(flow(together(self.add_column_button, self.order_button),
                           together(self.fill_button, self.wsp_button, self.open_button),
                           spacing=6))

        self.table = QtWidgets.QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        # the columns of annotations are typed into; the rest is read
        self.table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.DoubleClicked
            | QtWidgets.QAbstractItemView.EditTrigger.EditKeyPressed
            | QtWidgets.QAbstractItemView.EditTrigger.AnyKeyPressed)
        header = self.table.horizontalHeader()
        header.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        header.customContextMenuRequested.connect(self._header_menu)
        header.sectionDoubleClicked.connect(self._rename_column)
        header.setToolTip("Double-click a column of yours to rename it; right-click to remove it")
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemChanged.connect(self._ticked)

        # --- the plan's settings
        self.population = QtWidgets.QComboBox()
        self.population.setSizeAdjustPolicy(
            QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.population.setMinimumContentsLength(24)
        self.population.setToolTip(
            "Which events of each sample are concatenated: every event, or those of one "
            "gate — the Tcrb population, say. It has to exist in every sample that takes part.")
        self.include = QtWidgets.QPushButton()
        self.include.setMenu(QtWidgets.QMenu(self.include))
        self.include.setToolTip(
            "Which samples take part: tick or untick every sample of a group, a day, an "
            "animal at once — to concatenate two groups only, untick the others. The same "
            "as the boxes in the first column of the table.")
        self.extend = QtWidgets.QCheckBox("Gate the files the workspace does not list")
        self.extend.setToolTip(
            "Every sample of the workspace has the same gates and the same "
            "compensation matrix. Ticked, the files it does not list are given "
            "those too — what FlowJo does when a file is dropped into the group. "
            "Nobody has looked at those gates on those files: check them.")
        self.group_by = QtWidgets.QPushButton()
        self.group_by.setMenu(QtWidgets.QMenu(self.group_by))
        self.group_by.setToolTip(
            "Which column — or columns together — make a group, for balancing "
            "\"the same total from every group\" and for comparing later. Tick "
            "two and a group is each combination: Case · 3, Case · 7…")
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
        self.norm_batch = QtWidgets.QComboBox()
        self.norm_batch.setToolTip(
            "The column that says which batch a sample was acquired in — a day, a run. "
            "Batches are brought onto one another before the samples are concatenated.")
        self.norm_column = QtWidgets.QComboBox()
        self.norm_value = QtWidgets.QComboBox()
        for w in (self.norm_column, self.norm_value):
            w.setToolTip(
                "The reference: samples that should look alike in every batch. Whatever "
                "separates them between batches is taken to be technical and is removed "
                "from every sample of the batch — including any real difference between "
                "them. Aliquots of one control are the safe choice; a control group is an "
                "assumption that it does not change between batches.")
        self.norm_clusters = QtWidgets.QSpinBox()
        self.norm_clusters.setRange(1, 30)
        self.norm_clusters.setSpecialValueText("each marker whole")
        self.norm_clusters.setToolTip(
            "Above 1, the correction is learnt within that many FlowSOM clusters of the "
            "reference, as CytoNorm does: for a shift that differs between cell types. "
            "It needs many reference events — at least 50 per batch in each cluster.")
        self.norm_goal = QtWidgets.QComboBox()
        self.norm_goal.setToolTip(
            "What the batches are brought onto. The mean of the batches moves every "
            "batch a little; one batch as the target leaves that batch exactly as it was "
            "acquired and brings the others to it.")
        self.norm_everywhere = QtWidgets.QCheckBox("Also in Gates and Statistics")
        self.norm_everywhere.setToolTip(
            "Normalise the samples as they are read, not only the concatenated events: "
            "plots, gates and the statistics table then show normalised values. Gates "
            "drawn before were drawn on the values as acquired and may need another look. "
            "The curves are learnt from the files as acquired, and are refreshed at each "
            "Concatenate.")
        self.norm_result = QtWidgets.QLabel("")
        self.norm_result.setWordWrap(True)
        self.norm_result.setProperty("role", "caption")
        self.norm_look = QtWidgets.QPushButton("Before and after…")
        self.norm_look.setToolTip("The reference samples of every batch drawn over one "
                                  "another, marker by marker, before and after")
        self.norm_look.setVisible(False)
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
        form.addRow("Include", self.include)
        form.addRow(self.extend)
        form.addRow("Group by", self.group_by)
        form.addRow(self._eyebrow("Downsampling"))
        form.addRow("Balance", self.mode)
        form.addRow("N", self.n)
        form.addRow("Minimum", self.min_events)
        form.addRow("Seed", self.seed)
        form.addRow(self._eyebrow("Batch normalisation"))
        form.addRow("Batches", self.norm_batch)
        ref = QtWidgets.QHBoxLayout()
        ref.addWidget(self.norm_column, 1)
        ref.addWidget(self.norm_value, 1)
        form.addRow("Reference", ref)
        form.addRow("Within", self.norm_clusters)
        form.addRow("Bring to", self.norm_goal)
        form.addRow("", self.norm_everywhere)
        # the settings scroll where the window is short; the button that
        # acts on them stays in sight below
        inner = QtWidgets.QVBoxLayout()
        inner.setContentsMargins(0, 0, 6, 0)
        inner.addLayout(form)
        inner.addSpacing(8)
        inner.addWidget(self.summary)
        inner.addWidget(self.notes)
        inner.addWidget(self.norm_result)
        inner.addWidget(self.norm_look)
        inner.addStretch(1)
        settings = QtWidgets.QWidget()
        settings.setLayout(inner)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(settings)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        side = QtWidgets.QVBoxLayout()
        side.setContentsMargins(0, 0, 0, 0)
        side.addWidget(scroll, 1)
        side.addSpacing(6)
        side.addWidget(self.build_button)
        panel = QtWidgets.QWidget()
        panel.setLayout(side)
        panel.setFixedWidth(396)

        body = QtWidgets.QHBoxLayout()
        body.addWidget(self.table, 1)
        body.addWidget(panel)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.addLayout(top)
        layout.addLayout(body, 1)

        self.population.currentIndexChanged.connect(self._population_changed)
        self.extend.toggled.connect(self._extend_changed)
        self.add_column_button.clicked.connect(self._ask_new_column)
        self.fill_button.clicked.connect(self._ask_pattern)
        self.mode.currentIndexChanged.connect(self._settings_changed)
        self.n.valueChanged.connect(self._settings_changed)
        self.min_events.valueChanged.connect(self._settings_changed)
        self.seed.valueChanged.connect(self._settings_changed)
        self.norm_batch.currentIndexChanged.connect(self._normalisation_changed)
        self.norm_column.currentIndexChanged.connect(self._reference_column_changed)
        self.norm_value.currentIndexChanged.connect(self._normalisation_changed)
        self.norm_clusters.valueChanged.connect(self._normalisation_changed)
        self.norm_goal.currentIndexChanged.connect(self._normalisation_changed)
        self.norm_everywhere.toggled.connect(self._normalisation_changed)
        self.norm_look.clicked.connect(self._look_at_normalisation)
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
        self._fill_group_menu()
        self._fill_normalisation()
        for b in (self.add_column_button, self.fill_button):
            b.setEnabled(bool(s.entries))
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
        self.table.setHorizontalHeaderLabels(cols)
        self.table.setRowCount(len(s.entries))
        for r, e in enumerate(s.entries):
            tick = QtWidgets.QTableWidgetItem()
            flags = QtCore.Qt.ItemFlag.ItemIsUserCheckable | QtCore.Qt.ItemFlag.ItemIsEnabled
            tick.setFlags(flags if not e.problem else QtCore.Qt.ItemFlag.NoItemFlags)
            tick.setCheckState(QtCore.Qt.CheckState.Checked if e.include and not e.problem
                               else QtCore.Qt.CheckState.Unchecked)
            self.table.setItem(r, 0, tick)
            self.table.setItem(r, 1, self._fixed(e.name))
            for k in s.columns:
                cell = QtWidgets.QTableWidgetItem(e.meta.get(k, ""))
                if e.problem:
                    cell.setFlags(QtCore.Qt.ItemFlag.NoItemFlags)
                self.table.setItem(r, cols.index(k), cell)
            self._number(r, cols.index("Events"), e.events)
        self._filling = False
        self._update_plan()
        self.table.resizeColumnsToContents()
        self.table.setColumnWidth(0, 30)

    @staticmethod
    def _fixed(text: str = "") -> QtWidgets.QTableWidgetItem:
        """A cell to read and to select, not to type into."""
        item = QtWidgets.QTableWidgetItem(text)
        item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
        return item

    def _number(self, row: int, col: int, value, muted: bool = False) -> None:
        item = self._fixed("" if value is None else f"{value:,d}")
        item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight
                              | QtCore.Qt.AlignmentFlag.AlignVCenter)
        if muted:
            item.setForeground(QtGui.QBrush(QtGui.QColor(style.token("ink_faint"))))
        self.table.setItem(row, col, item)

    def _update_plan(self) -> None:
        """The In gate, Takes and Note columns, and the summary."""
        s = self.session
        self._fill_include_menu()
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
            elif e.edited and s.population != ALL_EVENTS:
                note = "your gates"
            else:
                note = ""
            item = self._fixed(note)
            quiet = note in ("", "counting…", "gates from the group", "your gates",
                             "sets the limit")
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
        """A box was ticked, or a cell of a column was typed into."""
        if self._filling:
            return
        s = self.session
        col = self._columns()[item.column()]
        if item.column() == 0:
            s.entries[item.row()].include = (
                item.checkState() == QtCore.Qt.CheckState.Checked)
        elif col in s.columns:
            # with several rows selected, what is typed in one goes into
            # all of them: ten samples of one group are one entry, not ten
            rows = {i.row() for i in self.table.selectedIndexes()}
            rows = rows if item.row() in rows else {item.row()}
            targets = [s.entries[r] for r in rows if not s.entries[r].problem]
            s.set_value(targets, col, item.text())
            self._filling = True
            for r in rows:
                if not s.entries[r].problem:
                    self.table.item(r, item.column()).setText(item.text().strip())
            self._filling = False
        else:
            return
        s.save_project()
        self._update_plan()
        if item.column() != 0:
            self._fill_normalisation()   # the values a reference can take may have changed
        else:
            self._show_normalisation()

    # --- batch normalisation ----------------------------------------------------------

    def _fill_normalisation(self) -> None:
        s = self.session
        cfg = s.normalization or {}
        filling, self._filling = self._filling, True
        self.norm_batch.clear()
        self.norm_batch.addItem("None", "")
        self.norm_column.clear()
        for c in s.columns:
            self.norm_batch.addItem(c, c)
            self.norm_column.addItem(c, c)
        self.norm_batch.setCurrentIndex(max(self.norm_batch.findData(cfg.get("batch", "")), 0))
        column, value = cfg.get("reference", ["", ""])
        self.norm_column.setCurrentIndex(max(self.norm_column.findData(column), 0))
        self._fill_reference_values(value)
        self.norm_clusters.setValue(int(cfg.get("clusters", 1)))
        self._fill_goals(cfg.get("goal", "mean"))
        self.norm_everywhere.setChecked(bool(cfg.get("everywhere")))
        self._filling = filling
        self._show_normalisation()

    def _fill_reference_values(self, keep: str = "") -> None:
        column = self.norm_column.currentData()
        values = order.ordered({e.meta.get(column, "") for e in self.session.usable()} - {""})
        self.norm_value.clear()
        for v in values:
            self.norm_value.addItem(v, v)
        self.norm_value.setCurrentIndex(max(self.norm_value.findData(keep), 0))

    def _fill_goals(self, keep: str = "mean") -> None:
        batch = self.norm_batch.currentData()
        values = order.ordered({e.meta.get(batch, "") for e in self.session.usable()} - {""}) \
            if batch else []
        self.norm_goal.clear()
        self.norm_goal.addItem("The mean of the batches", "mean")
        for v in values:
            self.norm_goal.addItem(f"{batch} {v}, left as acquired", v)
        self.norm_goal.setCurrentIndex(max(self.norm_goal.findData(keep), 0))

    def _reference_column_changed(self) -> None:
        if self._filling:
            return
        self._filling = True
        self._fill_reference_values()
        self._filling = False
        self._normalisation_changed()

    def _normalisation_changed(self) -> None:
        if self._filling:
            return
        batch = self.norm_batch.currentData() or None
        if batch != (self.session.normalization or {}).get("batch"):
            # the targets on offer are the batches of this column
            self._filling = True
            self._fill_goals()
            self._filling = False
        self.session.set_normalization(
            batch,
            (self.norm_column.currentData() or "", self.norm_value.currentData() or ""),
            self.norm_clusters.value(), self.norm_goal.currentData() or "mean",
            self.norm_everywhere.isChecked())
        self._show_normalisation()
        self.sigNormalisationChanged.emit()

    def set_normalisation(self, batch: str | None, column: str = "", value: str = "",
                          clusters: int = 1, goal: str = "mean",
                          everywhere: bool = False) -> None:
        self.session.set_normalization(batch, (column, value), clusters, goal, everywhere)
        self._fill_normalisation()
        self.sigNormalisationChanged.emit()

    def _look_at_normalisation(self) -> None:
        if self.session.normalization_check:
            from .normalization_dialog import NormalizationDialog
            NormalizationDialog(self.session.normalization_check, self).exec()

    def _show_normalisation(self) -> None:
        """What will be done, before; what was done, after."""
        s = self.session
        on = bool(s.normalization)
        self.norm_look.setVisible(False)
        for w in (self.norm_column, self.norm_value, self.norm_clusters, self.norm_goal,
                  self.norm_everywhere):
            w.setEnabled(bool(self.norm_batch.currentData()))
        if not on:
            self.norm_result.setText("")
            return
        refs = s.reference_batches()
        column, value = s.normalization["reference"]
        # what stands in the way comes before what was done last time: the
        # table may have been edited since
        if len(refs) < 2:
            self.norm_result.setText(
                f"Needs samples with {column} = {value} in at least two batches; "
                f"there are such samples in {len(refs)}.")
            return
        goal = s.normalization.get("goal", "mean")
        if goal != "mean" and goal not in refs:
            self.norm_result.setText(
                f"Batch {goal}, which the others are to be brought onto, has no sample "
                f"with {column} = {value}.")
            return
        done = (s.dataset.provenance.get("normalization") if s.dataset is not None else None)
        if done and all(done.get(k) == s.normalization[k] for k in s.normalization) \
                and done.get("reference_samples") == {b: [e.name for e in es]
                                                      for b, es in refs.items()}:
            before, after = done["distance_before"], done["distance_after"]
            worst = max(range(len(before)), key=before.__getitem__)
            text = (f"Normalised. References apart by {sum(before) / len(before):.4f} before "
                    f"and {sum(after) / len(after):.4f} after, on average over the markers "
                    f"(share of the axis, on reference events held out). Most moved: "
                    f"{done['markers'][worst]}, {before[worst]:.4f} → {after[worst]:.4f}.")
            if done["fallbacks"]:
                text += (f" {len(done['fallbacks'])} batch-and-cluster pair(s) had too few "
                         "reference events and were left as they were.")
            if done["not_normalised"]:
                text += (" Not normalised, their batch having no reference: "
                         + ", ".join(done["not_normalised"]) + ".")
            self.norm_result.setText(text)
            self.norm_look.setVisible(s.normalization_check is not None)
            return
        self.norm_result.setText(
            "At the next Concatenate: " + "; ".join(
                f"{b or '—'}: {len(es)} reference sample(s)" for b, es in sorted(refs.items()))
            + f". Whatever separates the {value} samples between batches will be removed "
            "from every sample — assuming they should not differ."
            + (" Gates and Statistics already show normalised values."
               if s.normalised_everywhere() else ""))

    # --- columns --------------------------------------------------------------------

    def _fill_group_menu(self) -> None:
        s = self.session
        menu = self.group_by.menu()
        menu.clear()
        for c in s.columns:
            action = menu.addAction(c)
            action.setCheckable(True)
            action.setChecked(c in s.group_by)
            action.toggled.connect(lambda on, name=c: self.set_group_by(name, on))
        self.group_by.setText(" + ".join(s.group_by) if s.group_by else "Nothing: one group")

    def set_group_by(self, column: str, on: bool) -> None:
        s = self.session
        chosen = [c for c in s.group_by if c != column] + ([column] if on else [])
        s.group_by = [c for c in s.columns if c in chosen]   # in the table's order
        self.group_by.setText(" + ".join(s.group_by) if s.group_by else "Nothing: one group")
        s.save_project()
        self._update_plan()
        self.sigColumnsChanged.emit()

    def _fill_include_menu(self) -> None:
        """All, none, and each value of each column: ticked where every
        sample that has it takes part, unticked where none does."""
        s = self.session
        usable = s.usable()
        menu = self.include.menu()
        menu.clear()
        menu.addAction("Every sample", lambda: self.include_where(None, None, True))
        menu.addAction("None", lambda: self.include_where(None, None, False))
        for column in s.columns:
            values = s.values_of(column)
            if not 1 < len(values) <= 40:
                continue
            sub = menu.addMenu(column)
            for value in values:
                with_it = [e for e in usable if e.meta.get(column, "") == value]
                action = sub.addAction(f"{value}    ({len(with_it)})")
                action.setCheckable(True)
                action.setChecked(all(e.include for e in with_it))
                action.triggered.connect(
                    lambda on, c=column, v=value: self.include_where(c, v, on))
        n = sum(e.include for e in usable)
        self.include.setText(f"{n} of {len(usable)} samples" if usable else "No samples")

    def include_where(self, column: str | None, value: str | None, on: bool) -> None:
        """Tick or untick every sample whose `column` is `value` — or, with
        no column, every sample."""
        s = self.session
        for e in s.usable():
            if column is None or e.meta.get(column, "") == value:
                e.include = bool(on)
        s.save_project()
        self._fill_table()
        self._update_plan()

    def ask_order(self, column: str = "") -> None:
        from .order_dialog import ask
        if ask(self.session, self, column if isinstance(column, str) else ""):
            self.order_changed()

    def order_changed(self) -> None:
        self._fill_normalisation()
        self._update_plan()
        self.sigOrderChanged.emit()

    def _columns_changed(self) -> None:
        self.session.save_project()
        self._fill_group_menu()
        self._fill_normalisation()
        self._fill_table()
        self.sigColumnsChanged.emit()

    def add_column(self, name: str) -> None:
        try:
            self.session.add_column(name)
        except ValueError as exc:
            self.sigStatus.emit(str(exc).capitalize())
            return
        self._columns_changed()

    def _ask_new_column(self) -> None:
        name, ok = QtWidgets.QInputDialog.getText(
            self, "Add a column", "Name (Day, Tissue, Donor, Treatment…)")
        if ok and name.strip():
            self.add_column(name)

    def remove_column(self, name: str) -> None:
        self.session.remove_column(name)
        self._columns_changed()

    def rename_column(self, old: str, new: str) -> None:
        try:
            self.session.rename_column(old, new)
        except ValueError as exc:
            self.sigStatus.emit(str(exc).capitalize())
            return
        self._columns_changed()

    def _rename_column(self, index: int) -> None:
        name = self._columns()[index]
        if name not in self.session.columns:
            return
        new, ok = QtWidgets.QInputDialog.getText(self, "Rename the column", "Name", text=name)
        if ok and new.strip() and new.strip() != name:
            self.rename_column(name, new)

    def _header_menu(self, pos) -> None:
        index = self.table.horizontalHeader().logicalIndexAt(pos)
        name = self._columns()[index] if index >= 0 else ""
        if name not in self.session.columns:
            return
        menu = QtWidgets.QMenu(self)
        rename = menu.addAction(f"Rename \"{name}\"…")
        remove = menu.addAction(f"Remove \"{name}\"")
        chosen = menu.exec(self.table.horizontalHeader().mapToGlobal(pos))
        if chosen is rename:
            self._rename_column(index)
        elif chosen is remove:
            self.remove_column(name)

    def fill_from_names(self, pattern: str) -> list[str]:
        filled = self.session.annotate(pattern)
        self._columns_changed()
        return filled

    def _ask_pattern(self) -> None:
        dialog = PatternDialog(self.session, self)
        if dialog.exec():
            self.fill_from_names(dialog.pattern.text())

    def _population_changed(self) -> None:
        if self._filling:
            return
        self.session.population = self.population.currentData() or ALL_EVENTS
        self.session.save_project()
        self._update_plan()
        self._count()

    def _extend_changed(self, on: bool) -> None:
        if self._filling:
            return
        self.session.set_extend_gates(on)
        self.refresh()
        self.sigGatesChanged.emit()

    def _settings_changed(self) -> None:
        if self._filling:
            return
        s = self.session
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
            self._show_normalisation()
            self.sigStatus.emit(f"Concatenated {len(ds):,d} events from "
                                f"{len(ds.samples)} samples"
                                + (", normalised by batch" if "normalization" in ds.provenance
                                   else ""))
            self.sigConcatenated.emit()

        def failed(message):
            self._busy = False
            self._update_plan()
            self.sigStatus.emit(message)

        worker.run(self.session.build, on_done=done, on_failed=failed,
                   on_progress=lambda d, t, n: self.sigStatus.emit(
                       f"Reading {n} ({d}/{t})"))


class PatternDialog(QtWidgets.QDialog):
    """Columns read out of the file names: a pattern in which every named
    group is a column, and what it finds in the first few names, live."""

    def __init__(self, session: Session, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Fill columns from the file names")
        self.session = session
        names = [e.name for e in session.usable()]
        self.pattern = QtWidgets.QLineEdit(suggest_pattern(names[0]) if names else "")
        self.help = QtWidgets.QLabel(
            "Each <b>(?P&lt;Name&gt;…)</b> is a column called Name, filled with what the "
            "part in its brackets matches: <b>\\d+</b> a number, <b>[A-Za-z]+</b> letters. "
            "Rename the groups to what they are — Day, Group, Animal — and delete the "
            "ones you do not want as columns.")
        self.help.setWordWrap(True)
        self.help.setProperty("role", "caption")
        self.preview = QtWidgets.QTableWidget()
        self.preview.verticalHeader().setVisible(False)
        self.preview.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.state = QtWidgets.QLabel("")
        self.state.setProperty("role", "warning")
        self.box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                              | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.box.accepted.connect(self.accept)
        self.box.rejected.connect(self.reject)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.pattern)
        layout.addWidget(self.help)
        layout.addWidget(self.preview, 1)
        layout.addWidget(self.state)
        layout.addWidget(self.box)
        self.resize(620, 420)
        self._names = names
        self.pattern.textChanged.connect(self._show)
        self._show()

    def _show(self) -> None:
        ok = self.box.button(QtWidgets.QDialogButtonBox.StandardButton.Ok)
        try:
            rx = re.compile(self.pattern.text())
        except re.error as exc:
            self.state.setText(f"Not a pattern yet: {exc}")
            ok.setEnabled(False)
            return
        groups = list(rx.groupindex)
        shown = self._names[:8]
        self.preview.setColumnCount(1 + len(groups))
        self.preview.setHorizontalHeaderLabels(["File"] + groups)
        self.preview.setRowCount(len(shown))
        for r, name in enumerate(shown):
            m = rx.match(name)
            self.preview.setItem(r, 0, QtWidgets.QTableWidgetItem(name))
            for k, g in enumerate(groups):
                self.preview.setItem(r, 1 + k, QtWidgets.QTableWidgetItem(
                    (m.group(g) or "") if m else "—"))
        matched = sum(1 for n in self._names if rx.match(n))
        self.state.setText("" if matched == len(self._names) and groups else
                           "No named group: nothing would be filled" if not groups else
                           f"Matches {matched} of {len(self._names)} file names; "
                           "the others keep what they have")
        ok.setEnabled(bool(groups) and matched > 0)
