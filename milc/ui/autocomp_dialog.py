"""
A spillover matrix from single-stain controls, opened from the
compensation dialog: choose the control files, see and change which
detector each is for, compute it off the window's thread, and look at it
beside the matrix the samples carry. What is accepted goes back into the
compensation dialog's table, to be given to the samples selected there or
to all of them — the same road as a matrix typed by hand.
"""

from __future__ import annotations

import os

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets

from .help_window import describe
from .. import autocomp
from . import theme, worker

UNSTAINED = "Unstained"


def percent(value: float) -> str:
    """A spillover as the table shows it: in percent, in at most five
    characters — 100, 43.0, 8.70, 0.27, 0."""
    v = 100 * value
    if abs(v) < 0.005:
        return "0"
    if abs(v - 100) < 0.05:
        return "100"
    return f"{v:.1f}" if abs(v) >= 10 else f"{v:.2f}"


class ControlsDialog(QtWidgets.QDialog):
    def __init__(self, detectors: list[str], acquisition=None, acquisition_label: str = "",
                 start: str = "", parent=None):
        super().__init__(parent)
        describe(self, "autocompensation")
        self.setWindowTitle("Compensation from single-stain controls")
        self.detectors = list(detectors)
        self.acquisition = acquisition
        self.start = start
        self.paths: list[str] = []
        #: what the last run gave, and whether one is running
        self.result: autocomp.Result | None = None
        self.running = False
        #: what failed last, for tests and the status line
        self.error = ""

        head = QtWidgets.QLabel(
            "One FCS file per fluorochrome, each stained with that fluorochrome alone, and "
            "optionally an unstained one. Each control is given the detector it is brightest "
            "in; change it where that is wrong. The matrix found is shown beside the one the "
            "samples carry, and goes into the compensation table only when you say so.")
        head.setWordWrap(True)
        head.setProperty("role", "caption")

        add = QtWidgets.QPushButton("Add files…")
        add.setToolTip("Choose the control FCS files")
        add.clicked.connect(self._browse)
        self.remove = QtWidgets.QPushButton("Remove")
        self.remove.setToolTip("Take the selected controls out of the list")
        self.remove.clicked.connect(self._remove)
        self.method = QtWidgets.QComboBox()
        for key, text in ((autocomp.AUTOSPILL, "AutoSpill (robust, refined)"),
                          (autocomp.CLASSIC, "Classic (positive vs negative)")):
            self.method.addItem(text, key)
        self.method.setToolTip(
            "AutoSpill: a robust regression of every detector on the control's own, then "
            "corrected on the compensated controls until no spill is left (Roca et al. "
            "2021). Classic: the medians of the positive and the negative population.")
        self.gate = QtWidgets.QCheckBox("Gate on scatter")
        self.gate.setToolTip("Keep the main population on forward against side scatter "
                             "(a robust ellipse round it): debris and clumps out")
        self.universal = QtWidgets.QCheckBox("Unstained as negative")
        self.universal.setToolTip(
            "Classic only: the unstained control's events as every control's negative, in "
            "place of the control's own negative events")
        self.run_button = QtWidgets.QPushButton("Compute")
        self.run_button.setToolTip("Compute the matrix from the controls in the list")
        self.run_button.clicked.connect(self.run)
        for b in (add, self.remove, self.run_button):
            b.setAutoDefault(False)
        tools = QtWidgets.QHBoxLayout()
        tools.addWidget(add)
        tools.addWidget(self.remove)
        tools.addSpacing(12)
        tools.addWidget(self.method)
        tools.addWidget(self.gate)
        tools.addWidget(self.universal)
        tools.addStretch(1)
        tools.addWidget(self.run_button)

        self.files = QtWidgets.QTableWidget(0, 4)
        self.files.setHorizontalHeaderLabels(["Control", "Events", "Detector", "Brightest in"])
        self.files.verticalHeader().setVisible(False)
        self.files.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.files.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        h = self.files.horizontalHeader()
        h.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        for c in (1, 3):
            h.setSectionResizeMode(c, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        h.setSectionResizeMode(2, QtWidgets.QHeaderView.ResizeMode.Interactive)
        self.files.setMinimumHeight(140)

        self.status = QtWidgets.QLabel("Add the control files.")
        self.status.setWordWrap(True)
        self.status.setProperty("role", "caption")
        self.progress = QtWidgets.QProgressBar()
        self.progress.setVisible(False)
        self.progress.setMaximumHeight(14)
        self.progress.setTextVisible(False)

        self.found = self._matrix_table()
        self.carried = self._matrix_table()
        left = QtWidgets.QVBoxLayout()
        self.found_title = QtWidgets.QLabel("From the controls")
        self.found_title.setProperty("role", "strong")
        left.addWidget(self.found_title)
        left.addWidget(self.found, 1)
        right = QtWidgets.QVBoxLayout()
        title = QtWidgets.QLabel(acquisition_label or "The samples' matrix")
        title.setProperty("role", "strong")
        right.addWidget(title)
        right.addWidget(self.carried, 1)
        both = QtWidgets.QHBoxLayout()
        both.addLayout(left, 1)
        both.addLayout(right, 1)
        self.legend = QtWidgets.QLabel(
            "Percent of the row's fluorochrome seen in the column's detector. In colour: "
            "more than 1 point away from the samples' matrix.")
        self.legend.setWordWrap(True)
        self.legend.setProperty("role", "caption")
        self.notes = QtWidgets.QLabel("")
        self.notes.setWordWrap(True)
        self.notes.setProperty("role", "warning")

        self.use = QtWidgets.QPushButton("Put in the compensation table")
        self.use.setToolTip("Close this and put the matrix found in the compensation table, "
                            "to be given to the samples selected there, or to every sample")
        self.use.clicked.connect(self.accept)
        self.save = QtWidgets.QPushButton("Save as CSV…")
        self.save.setToolTip("The matrix found, as proportions, detectors in the first row "
                             "and column")
        self.save.clicked.connect(self._save)
        cancel = QtWidgets.QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        for b in (self.use, self.save, cancel):
            b.setAutoDefault(False)
        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(self.use)
        buttons.addWidget(self.save)
        buttons.addStretch(1)
        buttons.addWidget(cancel)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(head)
        layout.addLayout(tools)
        layout.addWidget(self.files, 1)
        layout.addWidget(self.status)
        layout.addWidget(self.progress)
        layout.addLayout(both, 2)
        layout.addWidget(self.legend)
        layout.addWidget(self.notes)
        layout.addLayout(buttons)

        self.method.currentIndexChanged.connect(self._enable)
        self.files.itemSelectionChanged.connect(self._enable)
        self._fill(self.carried, acquisition.matrix if acquisition is not None else None,
                   acquisition.channels if acquisition is not None else None)
        self._fill(self.found, None)
        self._enable()
        n = len(self.detectors)
        self.resize(min(1040, 420 + 2 * 62 * n), min(760, 460 + 24 * n))

    # --- the controls -----------------------------------------------------------------

    def _browse(self) -> None:
        chosen, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, "Single-stain controls", self.start, "FCS files (*.fcs *.FCS);;All files (*)")
        if chosen:
            self.add_files(chosen)

    def add_files(self, paths: list[str]) -> None:
        """Read the new files off the window's thread and give each its
        brightest detector; what was already set by hand stays."""
        new = [p for p in paths if p not in self.paths]
        if not new:
            return
        self.paths += new
        everything, detectors = list(self.paths), list(self.detectors)

        def work():
            controls, wrong = [], []
            for p in everything:
                try:
                    controls += autocomp.load_controls([p], detectors, max_events=20_000)
                except Exception as exc:    # a file that cannot be used is said, not fatal
                    wrong.append(str(exc))
            if controls:
                try:
                    autocomp.assign(controls, detectors)
                except ValueError as exc:
                    wrong.append(str(exc))
            b = autocomp.brightness(controls, detectors) if controls else []
            return [(c, detectors[int(np.argmax(row))]) for c, row in zip(controls, b)], wrong

        self._busy("Reading the controls…")
        worker.run(work, on_done=self._read, on_failed=self._read_failed)

    def _read(self, answer) -> None:
        found, wrong = answer
        self._idle()
        kept = {p: self.primary_of(p) for p in self._listed()}
        self.files.setRowCount(0)
        for c, brightest in found:
            r = self.files.rowCount()
            self.files.insertRow(r)
            item = QtWidgets.QTableWidgetItem(c.name)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, c.path)
            item.setToolTip(c.path)
            self.files.setItem(r, 0, item)
            n = QtWidgets.QTableWidgetItem(f"{c.total:,d}")
            n.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight
                               | QtCore.Qt.AlignmentFlag.AlignVCenter)
            self.files.setItem(r, 1, n)
            combo = QtWidgets.QComboBox()
            combo.addItem(UNSTAINED, None)
            for d in self.detectors:
                combo.addItem(d, d)
            primary = kept[c.path] if c.path in kept else c.primary
            combo.setCurrentIndex(max(combo.findData(primary), 0))
            combo.currentIndexChanged.connect(self._enable)
            combo.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToContents)
            self.files.setCellWidget(r, 2, combo)
            self.files.setColumnWidth(2, max(self.files.columnWidth(2),
                                             combo.sizeHint().width() + 24))
            self.files.setItem(r, 3, QtWidgets.QTableWidgetItem(
                "" if c.unstained else brightest))
        self.paths = self._listed()
        self.status.setText(f"{len(found)} control(s). Check the detector of each, then "
                            "compute.")
        self._enable()
        if wrong:
            self.error = "; ".join(wrong)
            self.status.setText("Left out: " + "; ".join(wrong) + ".")

    def _read_failed(self, message: str) -> None:
        self._idle()
        self.error = message
        # the files that could not be read are not kept in the list
        self.paths = self._listed()
        self.status.setText("Could not read the controls: " + message.split(": ", 1)[-1])
        self._enable()

    def _listed(self) -> list[str]:
        return [self.files.item(r, 0).data(QtCore.Qt.ItemDataRole.UserRole)
                for r in range(self.files.rowCount())]

    def primary_of(self, path: str):
        """The detector the list gives this control; None for unstained."""
        r = self._listed().index(path)
        return self.files.cellWidget(r, 2).currentData()

    def set_primary(self, name: str, detector: str | None) -> None:
        """For tests and scripts: the detector of the control called `name`."""
        for r in range(self.files.rowCount()):
            if self.files.item(r, 0).text() == name:
                combo = self.files.cellWidget(r, 2)
                combo.setCurrentIndex(combo.findData(detector))
                return
        raise KeyError(name)

    def _remove(self) -> None:
        rows = sorted({i.row() for i in self.files.selectedIndexes()}, reverse=True)
        for r in rows:
            path = self.files.item(r, 0).data(QtCore.Qt.ItemDataRole.UserRole)
            self.paths.remove(path)
            self.files.removeRow(r)
        self._enable()

    def problem(self) -> str:
        """Why the controls cannot be computed from as they are set; ""."""
        listed = self._listed()
        if not listed:
            return "Add the control files."
        given = [self.primary_of(p) for p in listed]
        stained = [d for d in given if d is not None]
        if not stained:
            return "No control has a detector."
        twice = sorted({d for d in stained if stained.count(d) > 1})
        if twice:
            return f"Two controls are given {', '.join(twice)}: each detector takes one."
        if given.count(None) > 1:
            return "Only one control can be the unstained one."
        return ""

    # --- running -----------------------------------------------------------------------

    def _enable(self, *_) -> None:
        why = self.problem()
        self.run_button.setEnabled(not why and not self.running)
        self.remove.setEnabled(bool(self.files.selectedIndexes()) and not self.running)
        self.universal.setEnabled(self.method.currentData() == autocomp.CLASSIC)
        has = self.result is not None
        self.use.setEnabled(has and not self.running)
        self.save.setEnabled(has and not self.running)
        if why and not self.running and self.files.rowCount():
            self.status.setText(why)

    def _busy(self, text: str) -> None:
        self.running = True
        self.status.setText(text)
        self.progress.setRange(0, 0)
        self.progress.setVisible(True)
        self._enable()

    def _idle(self) -> None:
        self.running = False
        self.progress.setVisible(False)
        self._enable()

    def run(self) -> None:
        if self.problem() or self.running:
            return
        paths = self._listed()
        primaries = {p: self.primary_of(p) for p in paths}
        method, gate = self.method.currentData(), self.gate.isChecked()
        universal = self.universal.isChecked()
        detectors, base = list(self.detectors), self.acquisition
        self.error = ""

        def work(progress):
            return autocomp.compute(paths, detectors, method, primaries, gate=gate, base=base,
                                    max_events=100_000, progress=progress,
                                    universal_negative=universal)

        self._busy("Computing…")
        worker.run(work, on_done=self._done, on_failed=self._failed,
                   on_progress=self._progress)

    def _progress(self, done: int, total: int, what: str) -> None:
        self.progress.setRange(0, max(total, 1))
        self.progress.setValue(done)
        self.status.setText(what[:1].upper() + what[1:] + "…")

    def _failed(self, message: str) -> None:
        self._idle()
        self.error = message
        self.status.setText("Could not compute the matrix: " + message.split(": ", 1)[-1])

    def _done(self, result) -> None:
        self.result = result
        self._idle()
        self._fill(self.found, result.matrix, result.channels)
        self.found_title.setText("From the controls — " + (
            "AutoSpill" if result.method == autocomp.AUTOSPILL else "classic"))
        parts = [f"{len(result.controls)} of {len(self.detectors)} detectors computed"]
        if result.history:
            parts.append(f"{len(result.history)} iterations, "
                         + ("converged" if result.converged else "did not converge"))
        if result.residual is not None:
            parts.append(f"largest spill left in the controls once compensated "
                         f"{100 * result.largest_residual():.2g}%")
        self.status.setText("; ".join(parts) + ".")
        self.notes.setText("\n".join(n[:1].upper() + n[1:] for n in result.notes))

    # --- the two matrices --------------------------------------------------------------

    def _matrix_table(self) -> QtWidgets.QTableWidget:
        t = QtWidgets.QTableWidget(len(self.detectors), len(self.detectors))
        short = [d[:-2] if d.endswith(("-A", "-H")) else d for d in self.detectors]
        # a long name on two lines, so that eight columns of each matrix fit side by side
        t.setHorizontalHeaderLabels([n.replace("-", "-\n", 1) if len(n) > 5 else n
                                     for n in short])
        t.setVerticalHeaderLabels(short)
        for j, d in enumerate(self.detectors):
            t.horizontalHeaderItem(j).setToolTip(d)
        t.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        t.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
        t.horizontalHeader().setMinimumSectionSize(36)
        t.verticalHeader().setDefaultSectionSize(22)
        font = t.font()
        font.setPointSizeF(max(font.pointSizeF() - 1, 8))
        t.setFont(font)
        return t

    def _fill(self, table, matrix, channels=None) -> None:
        n = len(self.detectors)
        other = self._carried_matrix()
        for i in range(n):
            for j in range(n):
                value = None
                if matrix is not None:
                    ch = channels or self.detectors
                    a, b = self.detectors[i], self.detectors[j]
                    if a in ch and b in ch:
                        value = matrix[ch.index(a), ch.index(b)]
                item = QtWidgets.QTableWidgetItem("" if value is None else percent(value))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight
                                      | QtCore.Qt.AlignmentFlag.AlignVCenter)
                if (table is self.found and value is not None and other is not None
                        and abs(value - other[i, j]) > 0.01):
                    item.setForeground(QtGui.QColor(theme.warning()))
                    item.setToolTip(f"{100 * other[i, j]:.2f} in the samples' matrix")
                table.setItem(i, j, item)

    def _carried_matrix(self) -> np.ndarray | None:
        a = self.acquisition
        if a is None:
            return None
        return np.array([[a.matrix[a.channels.index(x), a.channels.index(y)]
                          if x in a.channels and y in a.channels else float(x == y)
                          for y in self.detectors] for x in self.detectors])

    def matrix(self):
        """The matrix found, as a `Compensation`; None before a run."""
        return self.result.compensation() if self.result is not None else None

    def _save(self) -> None:
        if self.result is None:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save the matrix", os.path.join(self.start, "spillover.csv"), "CSV (*.csv)")
        if path:
            autocomp.write_csv(path, self.result)
            self.status.setText(f"Saved to {path}")
