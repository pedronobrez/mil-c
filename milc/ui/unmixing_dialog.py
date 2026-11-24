"""
Spectral unmixing, chosen and looked at: which files are the single-stain
reference controls and which fluorochrome each one is, the unstained
control, autofluorescence, OLS or WLS — then the spectra those give, how
alike they are (the similarity matrix and the complexity index), how well
the set explains each control, and Apply.

Applying unmixes every sample that has the spectra's detectors in place
of compensating it (`Session.set_unmixing`); Stop unmixing goes back to
compensation. The dialog is reached from the Gates tab (`button`).
"""

from __future__ import annotations

import dataclasses
import os
import re

import numpy as np
import pyqtgraph as pg
from PyQt6 import QtCore, QtGui, QtWidgets

from .. import unmixing as um
from . import palette, style, theme, worker

#: Cytek's warning: two spectra more alike than this are hard to unmix
SIMILAR = 0.98
CONTROL_COLUMNS = ["Use", "File", "Fluorochrome", "Marker", "Negative"]
CHECK_COLUMNS = ["Fluorochrome", "Peak", "Positives", "Unexplained", "Spill", "Into"]
NEGATIVES = (("Own dim", "internal"), ("Unstained", "unstained"))


def guess_fluor(path: str) -> str:
    """`A3 BUV395 (Beads)_Bead controls_1.fcs` → `BUV395`: the well, the
    kind of particle and SpectroFlo's group taken off the name."""
    name = os.path.splitext(os.path.basename(path))[0]
    name = re.sub(r"^[A-P]\d{1,2}\s+", "", name)
    name = re.sub(r"\s*\((Beads|Cells)\).*$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"_\d+$", "", name)
    return name.strip() or os.path.basename(path)


def looks_like_control(path: str) -> bool:
    name = os.path.basename(path).lower()
    return ("(beads)" in name or "(cells)" in name) and not looks_unstained(path)


def looks_unstained(path: str) -> bool:
    name = os.path.basename(path).lower()
    return "unstained" in name or "negative" in name


def button(tab) -> QtWidgets.QPushButton:
    """The Gates tab's way in: a button beside the others."""
    b = QtWidgets.QPushButton("Unmixing…")
    b.setToolTip("Spectral unmixing in place of compensation, for spectral cytometers: "
                 "reference spectra from single-stain controls, OLS or WLS")
    b.clicked.connect(lambda: UnmixingDialog(tab.session, lambda u: apply_to_tab(tab, u),
                                             tab).exec())
    return b


def apply_to_tab(tab, unmixing) -> None:
    """Unmix (or stop) and show the sample again with its new channels."""
    path = tab.node.path if getattr(tab, "node", None) is not None else ""
    tab.session.set_unmixing(unmixing)
    if tab.entry is not None:
        tab.sample = tab.session.load(tab.entry)
        tab._fill_axes()
        tab._fill_tree(select=path)
    tab.show_review()
    tab.sigGatesChanged.emit()


class UnmixingDialog(QtWidgets.QDialog):
    def __init__(self, session, on_apply, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Spectral unmixing")
        self.session, self.on_apply = session, on_apply
        #: what Compute or Read CSV gave, not yet applied
        self.unmixing: um.Unmixing | None = session.unmixing
        self.references: list[um.Reference] = []
        self.checks: list[dict] = []

        self.controls = QtWidgets.QTableWidget(0, len(CONTROL_COLUMNS))
        self.controls.setHorizontalHeaderLabels(CONTROL_COLUMNS)
        self.controls.verticalHeader().hide()
        head = self.controls.horizontalHeader()
        head.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Interactive)
        head.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        head.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.Stretch)
        head.setMinimumSectionSize(40)
        self.controls.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.controls.setColumnWidth(2, 110)
        self.controls.setColumnWidth(3, 60)
        self.controls.setColumnWidth(4, 118)
        self.controls.setToolTip("Tick each single-stain control and name its fluorochrome: "
                                 "the name becomes the unmixed channel's. The negative is "
                                 "the control's own dimmest events, or the unstained control")
        self.add_button = QtWidgets.QPushButton("Add files…")
        self.add_button.setToolTip("Reference controls kept outside the folder of samples")
        self.tick_all = QtWidgets.QPushButton("Tick all")
        self.untick_all = QtWidgets.QPushButton("Untick all")
        self.unstained = QtWidgets.QComboBox()
        self.unstained.setToolTip("Unstained cells or negative beads: the negative of the "
                                  "controls set to use it, the weights' floor for WLS, and "
                                  "the autofluorescence spectrum")
        self.autofluorescence = QtWidgets.QCheckBox("Autofluorescence as one more spectrum")
        self.autofluorescence.setToolTip("The unstained control's own spectrum, unmixed as "
                                         "a channel of its own, AF")
        self.method = QtWidgets.QComboBox()
        self.method.addItem("Ordinary least squares (OLS)", um.OLS)
        self.method.addItem("Weighted least squares (WLS)", um.WLS)
        self.method.setToolTip("WLS weighs each detector of each event by 1 / its value "
                               "(Poisson-like), as Cytek's WLS: slower, less spread")
        self.compute_button = QtWidgets.QPushButton("Compute spectra")
        self.compute_button.setDefault(True)
        self.read_button = QtWidgets.QPushButton("Read CSV…")
        self.read_button.setToolTip("Spectra from a table: one fluorochrome per row, one "
                                    "detector per column (or the other way round)")
        self.save_button = QtWidgets.QPushButton("Save CSV…")

        left = QtWidgets.QVBoxLayout()
        title = QtWidgets.QLabel("Reference controls")
        title.setFont(style.eyebrow_font())
        title.setProperty("role", "caption")
        left.addWidget(title)
        left.addWidget(self.controls, 1)
        row = QtWidgets.QHBoxLayout()
        for w in (self.add_button, self.tick_all, self.untick_all):
            row.addWidget(w)
        row.addStretch(1)
        left.addLayout(row)
        form = QtWidgets.QFormLayout()
        form.addRow("Unstained control", self.unstained)
        form.addRow("", self.autofluorescence)
        form.addRow("Method", self.method)
        left.addLayout(form)
        row = QtWidgets.QHBoxLayout()
        for w in (self.compute_button, self.read_button, self.save_button):
            row.addWidget(w)
        row.addStretch(1)
        left.addLayout(row)
        hint = QtWidgets.QLabel(
            "Controls in the folder of samples are also samples: untick them in Samples "
            "to keep them out of the concatenated events. While unmixing is on, the "
            "compensation matrix is not used.")
        hint.setWordWrap(True)
        hint.setProperty("role", "caption")
        left.addWidget(hint)
        left_box = QtWidgets.QWidget()
        left_box.setLayout(left)

        self.tabs = QtWidgets.QTabWidget()
        self.plot = pg.PlotWidget()
        self.plot.setBackground(theme.background())
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.plot.setMouseEnabled(False, False)
        self.highlight = QtWidgets.QComboBox()
        self.highlight.setToolTip("One spectrum drawn over the others, faint")
        spectra_box = QtWidgets.QWidget()
        spectra_layout = QtWidgets.QVBoxLayout(spectra_box)
        spectra_layout.setContentsMargins(0, 6, 0, 0)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Show"))
        row.addWidget(self.highlight, 1)
        spectra_layout.addLayout(row)
        spectra_layout.addWidget(self.plot, 1)
        self.tabs.addTab(spectra_box, "Spectra")
        sim_box = QtWidgets.QWidget()
        sim_layout = QtWidgets.QVBoxLayout(sim_box)
        sim_layout.setContentsMargins(0, 6, 0, 0)
        self.complexity = QtWidgets.QLabel("")
        self.complexity.setWordWrap(True)
        self.complexity.setProperty("role", "strong")
        self.similarity = QtWidgets.QTableWidget(0, 0)
        self.similarity.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        sim_layout.addWidget(self.complexity)
        sim_layout.addWidget(self.similarity, 1)
        self.tabs.addTab(sim_box, "Similarity")
        self.check = QtWidgets.QTableWidget(0, len(CHECK_COLUMNS))
        self.check.setHorizontalHeaderLabels(CHECK_COLUMNS)
        self.check.verticalHeader().hide()
        self.check.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.check.setToolTip(
            "Each control unmixed with the whole set, its brightest 5% of events: how much "
            "of them the spectra leave unexplained, and the largest share that lands in "
            "another channel. Both near 0 is a set that fits its controls.")
        self.tabs.addTab(self.check, "Check")

        split = QtWidgets.QSplitter()
        split.addWidget(left_box)
        split.addWidget(self.tabs)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 1)
        split.setSizes([540, 460])

        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        self.stop_button = QtWidgets.QPushButton("Stop unmixing")
        self.stop_button.setToolTip("Compensate the samples again, as before unmixing")
        self.apply_button = QtWidgets.QPushButton("Apply")
        self.apply_button.setToolTip("Unmix every sample that has these detectors")
        close = QtWidgets.QPushButton("Close")
        bottom = QtWidgets.QHBoxLayout()
        bottom.addWidget(self.status, 1)
        for w in (self.stop_button, self.apply_button, close):
            bottom.addWidget(w)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(split, 1)
        layout.addLayout(bottom)

        self.add_button.clicked.connect(self._add_files)
        self.tick_all.clicked.connect(lambda: self._tick(True))
        self.untick_all.clicked.connect(lambda: self._tick(False))
        self.compute_button.clicked.connect(self.compute)
        self.read_button.clicked.connect(self._read_csv)
        self.save_button.clicked.connect(self._save_csv)
        self.apply_button.clicked.connect(self.apply)
        self.stop_button.clicked.connect(self.stop)
        close.clicked.connect(self.reject)
        self.check.currentCellChanged.connect(self._check_row)
        self.highlight.currentIndexChanged.connect(lambda *_: self.draw())

        self._fill_files()
        if self.unmixing is not None:
            self.method.setCurrentIndex(self.method.findData(self.unmixing.method))
            self.autofluorescence.setChecked(um.AF in self.unmixing.fluors)
        self.show_result()
        self.resize(1020, 680)

    # --- the controls -----------------------------------------------------------

    def _fill_files(self) -> None:
        paths = [e.path for e in self.session.entries if not e.problem]
        self.unstained.addItem("None", "")
        for p in paths:
            self.unstained.addItem(os.path.basename(p), p)
        previous = {}
        if self.unmixing is not None:
            previous = dict(zip(self.unmixing.sources, zip(self.unmixing.fluors,
                                                           self.unmixing.markers)))
        for p in paths:
            self._add_row(p, previous.get(os.path.basename(p)))
        for i, p in enumerate(paths):
            if looks_unstained(p) and "cell" in os.path.basename(p).lower():
                self.unstained.setCurrentIndex(i + 1)
                break
        else:
            for i, p in enumerate(paths):
                if looks_unstained(p):
                    self.unstained.setCurrentIndex(i + 1)
                    break

    def _add_row(self, path: str, before: tuple | None = None) -> None:
        r = self.controls.rowCount()
        self.controls.insertRow(r)
        use = QtWidgets.QTableWidgetItem()
        use.setFlags(QtCore.Qt.ItemFlag.ItemIsUserCheckable | QtCore.Qt.ItemFlag.ItemIsEnabled)
        ticked = before is not None or (self.unmixing is None and looks_like_control(path))
        use.setCheckState(QtCore.Qt.CheckState.Checked if ticked
                          else QtCore.Qt.CheckState.Unchecked)
        self.controls.setItem(r, 0, use)
        name = QtWidgets.QTableWidgetItem(os.path.basename(path))
        name.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
        name.setData(QtCore.Qt.ItemDataRole.UserRole, path)
        name.setToolTip(path)
        self.controls.setItem(r, 1, name)
        self.controls.setItem(r, 2, QtWidgets.QTableWidgetItem(before[0] if before
                                                               else guess_fluor(path)))
        self.controls.setItem(r, 3, QtWidgets.QTableWidgetItem(before[1] if before else ""))
        negative = QtWidgets.QComboBox()
        for text, key in NEGATIVES:
            negative.addItem(text, key)
        self.controls.setCellWidget(r, 4, negative)

    def _add_files(self) -> None:
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, "Reference controls", self.session.folder, "FCS files (*.fcs)")
        for p in paths:
            self._add_row(p)
            self.unstained.addItem(os.path.basename(p), p)

    def _tick(self, on: bool) -> None:
        state = QtCore.Qt.CheckState.Checked if on else QtCore.Qt.CheckState.Unchecked
        for r in range(self.controls.rowCount()):
            self.controls.item(r, 0).setCheckState(state)

    def chosen(self) -> list[dict]:
        """The ticked controls, as `unmixing.read_references` takes them."""
        out = []
        for r in range(self.controls.rowCount()):
            if self.controls.item(r, 0).checkState() != QtCore.Qt.CheckState.Checked:
                continue
            out.append({"path": self.controls.item(r, 1).data(QtCore.Qt.ItemDataRole.UserRole),
                        "fluor": self.controls.item(r, 2).text().strip(),
                        "marker": self.controls.item(r, 3).text().strip(),
                        "negative": self.controls.cellWidget(r, 4).currentData()})
        return out

    # --- computing -----------------------------------------------------------------

    def compute(self) -> None:
        controls = self.chosen()
        unstained = self.unstained.currentData() or None
        af = self.autofluorescence.isChecked()
        method = self.method.currentData()
        if not controls:
            self.say("Tick the reference controls first.", warn=True)
            return
        if any(not c["fluor"] for c in controls):
            self.say("Every ticked control needs the name of its fluorochrome.", warn=True)
            return
        if (af or any(c["negative"] == "unstained" for c in controls)) and not unstained:
            self.say("Choose the unstained control.", warn=True)
            return
        controls = [c for c in controls if c["path"] != unstained]

        def job():
            dets, refs, noise = um.read_references(controls, unstained=unstained,
                                                   autofluorescence=af)
            u = um.Unmixing.from_references(dets, refs, method, noise)
            return u, refs, um.check_references(u, controls)

        self.compute_button.setEnabled(False)
        self.say(f"Reading {len(controls)} controls…")
        worker.run(job, self._computed, self._failed)

    def _computed(self, result) -> None:
        self.compute_button.setEnabled(True)
        self.unmixing, self.references, self.checks = result
        self.say(f"{len(self.unmixing.fluors)} spectra on {len(self.unmixing.detectors)} "
                 "detectors. Apply to unmix the samples with them.")
        self.show_result()

    def _failed(self, message: str) -> None:
        self.compute_button.setEnabled(True)
        self.say(message, warn=True)

    def _read_csv(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Spectra", self.session.folder, "CSV files (*.csv)")
        if path:
            self.read_csv(path)

    def read_csv(self, path: str) -> None:
        usable = self.session.usable()
        try:
            known = um.detector_channels(usable[0].channels) if usable else None
            try:
                dets, names, spectra = um.read_csv(path, known)
            except ValueError:
                dets, names, spectra = um.read_csv(path)
            self.unmixing = um.Unmixing(dets, names, spectra, method=self.method.currentData(),
                                        sources=[os.path.basename(path)] * len(names))
        except (OSError, ValueError) as exc:
            self.say(f"{os.path.basename(path)}: {exc}", warn=True)
            return
        self.references, self.checks = [], []
        self.say(f"{len(names)} spectra read from {os.path.basename(path)}.")
        self.show_result()

    def _save_csv(self) -> None:
        if self.unmixing is None:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Spectra", os.path.join(self.session.folder, "spectra.csv"),
            "CSV files (*.csv)")
        if path:
            um.write_csv(path, self.unmixing)
            self.say(f"Wrote {path}")

    # --- showing -----------------------------------------------------------------------

    def say(self, text: str, warn: bool = False) -> None:
        self.status.setText(text)
        self.status.setProperty("role", "warning" if warn else "")
        self.status.style().unpolish(self.status)
        self.status.style().polish(self.status)

    def show_result(self) -> None:
        u = self.unmixing
        has = u is not None
        self.apply_button.setEnabled(has)
        self.save_button.setEnabled(has)
        self.stop_button.setEnabled(self.session.unmixing is not None)
        if has and u is self.session.unmixing and not self.status.text():
            self.say(f"Unmixing is on: {len(u.fluors)} spectra, {u.method.upper()}.")
        self._fill_highlight()
        self._fill_similarity()
        self._fill_check()
        self.draw()

    def draw(self) -> None:
        self.plot.clear()
        theme.style_axes(self.plot.plotItem)
        u = self.unmixing
        if u is None:
            return
        x = np.arange(len(u.detectors))
        chosen = self.highlight.currentIndex() - 1     # 0 is every spectrum
        colours = palette.categorical(len(u.fluors))
        legend = self.plot.plotItem.legend or self.plot.plotItem.addLegend(offset=(-8, 8))
        legend.clear()
        legend.setBrush(theme.legend_brush())
        legend.setPen(theme.legend_pen())
        legend.setLabelTextColor(theme.foreground())
        for k, (name, row, colour) in enumerate(zip(u.fluors, u.spectra, colours)):
            bold = k == chosen
            if chosen >= 0 and not bold:
                faint = QtGui.QColor(theme.muted())
                faint.setAlpha(70)
                pen = pg.mkPen(faint, width=1.0)
            else:
                pen = pg.mkPen(colour, width=3.0 if bold else 1.4)
            curve = pg.PlotCurveItem(x, row, pen=pen)
            curve.setZValue(10 if bold else 0)
            self.plot.addItem(curve)
            if (chosen < 0 and len(u.fluors) <= 12) or bold:
                legend.addItem(curve, name)
        step = max(1, -(-len(u.detectors) // 8))
        self.plot.getAxis("bottom").setTicks(
            [[(i, u.detectors[i]) for i in range(0, len(u.detectors), step)]])
        self.plot.setLabel("bottom", "Detector")
        self.plot.setLabel("left", "Normalised emission")
        self.plot.setXRange(0, len(u.detectors) - 1, padding=0.02)
        self.plot.setYRange(-0.02, 1.05, padding=0)

    def _check_row(self, row: int, *_) -> None:
        item = self.check.item(row, 0) if row >= 0 else None
        if item is not None and self.unmixing is not None:
            self.highlight.setCurrentIndex(self.unmixing.fluors.index(item.text()) + 1)

    def _fill_highlight(self) -> None:
        keep = self.highlight.currentText()
        self.highlight.blockSignals(True)
        self.highlight.clear()
        self.highlight.addItem("Every spectrum")
        if self.unmixing is not None:
            self.highlight.addItems(self.unmixing.fluors)
        i = self.highlight.findText(keep)
        self.highlight.setCurrentIndex(max(i, 0))
        self.highlight.blockSignals(False)

    def _fill_similarity(self) -> None:
        u = self.unmixing
        if u is None:
            self.similarity.setRowCount(0)
            self.similarity.setColumnCount(0)
            self.complexity.setText("")
            return
        sim = u.similarity()
        k = len(u.fluors)
        self.similarity.setRowCount(k)
        self.similarity.setColumnCount(k)
        self.similarity.setHorizontalHeaderLabels(u.fluors)
        self.similarity.setVerticalHeaderLabels(u.fluors)
        warn = QtGui.QColor(theme.danger())
        for i in range(k):
            for j in range(k):
                item = QtWidgets.QTableWidgetItem(f"{sim[i, j] + 0.0:.2f}".replace("-0.00", "0.00"))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                if i != j and sim[i, j] > SIMILAR:
                    item.setForeground(warn)
                if i != j and sim[i, j] >= 0.9:
                    font = item.font()
                    font.setBold(True)
                    item.setFont(font)
                self.similarity.setItem(i, j, item)
        self.similarity.resizeColumnsToContents()
        iu = np.triu_indices(k, 1)
        if k > 1:
            t = int(np.argmax(sim[iu]))
            a, b = iu[0][t], iu[1][t]
            pair = f"; the most alike: {u.fluors[a]} and {u.fluors[b]}, {sim[a, b]:.2f}"
        else:
            pair = ""
        over = int((sim[iu] > SIMILAR).sum())
        self.complexity.setText(
            f"Complexity index {u.complexity():.2f}{pair}."
            + (f" {over} pair(s) above {SIMILAR}: hard to unmix." if over else ""))

    def _fill_check(self) -> None:
        refs = {r.fluor: r for r in self.references}
        rows = self.checks
        self.check.setRowCount(len(rows))
        for r, c in enumerate(rows):
            ref = refs.get(c["fluor"])
            values = [c["fluor"], ref.peak if ref else "",
                      f"{ref.positive_events:,d}" if ref else "",
                      f"{100 * c['residual']:.1f}%", f"{100 * c['leak']:.1f}%", c["leak_into"]]
            for col, v in enumerate(values):
                self.check.setItem(r, col, QtWidgets.QTableWidgetItem(v))
        self.check.resizeColumnsToContents()

    # --- applying ------------------------------------------------------------------------

    def apply(self) -> None:
        if self.unmixing is None:
            return
        if self.unmixing.method != self.method.currentData():
            self.unmixing = dataclasses.replace(self.unmixing, method=self.method.currentData())
        lost = um.lost_gate_channels(self.unmixing, self.session.usable())
        if lost:
            self.say("Gates are drawn on " + ", ".join(lost[:6])
                     + ", which unmixing replaces: remove those gates first.", warn=True)
            return
        try:
            self.on_apply(self.unmixing)
        except ValueError as exc:
            self.say(str(exc), warn=True)
            return
        self.say(f"Unmixing is on: {len(self.unmixing.fluors)} spectra, "
                 f"{self.unmixing.method.upper()}.")
        self.show_result()

    def stop(self) -> None:
        self.on_apply(None)
        self.say("Unmixing is off: the samples are compensated as before.")
        self.show_result()
