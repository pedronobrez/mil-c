"""
File ▸ Export gates…: the gates written for another program, and what
was left out of them.
"""

from __future__ import annotations

import os

from PyQt6 import QtWidgets

from .help_window import describe
from ..gate_export import FORMATS, GATINGML, WSP, Report

_EXTENSION = {WSP: ".wsp", GATINGML: ".xml"}
_ABOUT = {
    WSP: "One workspace with every sample: its file, its compensation, the scales and "
         "its gates. It has been read back here and by FlowKit, event for event; it has "
         "not been opened in FlowJo.",
    GATINGML: "The standard of ISAC, read by FlowKit, flowCore and others. One file where "
              "every sample is gated and compensated alike, else one for each sample, "
              "named after it.",
}


class GatesExportDialog(QtWidgets.QDialog):
    def __init__(self, folder: str, gated: int, ticked: int, parent=None):
        super().__init__(parent)
        describe(self, "workspaces-and-gating-ml")
        self.setWindowTitle("Export gates")
        self.format = QtWidgets.QComboBox()
        for key, name in FORMATS.items():
            self.format.addItem(f"{name} ({_EXTENSION[key]})", key)
        self.about = QtWidgets.QLabel()
        self.about.setWordWrap(True)
        self.about.setProperty("role", "hint")
        self.ticked = QtWidgets.QCheckBox("Only the samples ticked in Samples")
        self.ticked.setToolTip(f"Unticked: all {gated} with gates. Ticked: {ticked} of them")
        self.counts = QtWidgets.QCheckBox("Write each population's count")
        self.counts.setChecked(True)
        self.counts.setToolTip("Opens every file. The counts are this program's: another "
                               "has them to compare with its own")
        self.path = QtWidgets.QLineEdit(os.path.join(folder, "gates.wsp") if folder else "")
        browse = QtWidgets.QPushButton("Choose…")
        browse.clicked.connect(self._browse)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.path, 1)
        row.addWidget(browse)
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                         | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        box.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setText("Export")
        box.accepted.connect(self.accept)
        box.rejected.connect(self.reject)
        form = QtWidgets.QFormLayout(self)
        form.addRow("Format", self.format)
        form.addRow(self.about)
        form.addRow(self.ticked)
        form.addRow(self.counts)
        form.addRow("Into", row)
        form.addRow(box)
        self.format.currentIndexChanged.connect(self._format_changed)
        self._format_changed()
        self.resize(600, 260)

    def _format_changed(self) -> None:
        """The file's extension follows the format; the rest of its name stays."""
        key = self.format.currentData()
        self.about.setText(_ABOUT[key])
        stem, ext = os.path.splitext(self.path.text())
        if stem and ext.lower() in _EXTENSION.values():
            self.path.setText(stem + _EXTENSION[key])

    def _browse(self) -> None:
        key = self.format.currentData()
        chosen, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export the gates", self.path.text(),
            f"{FORMATS[key]} (*{_EXTENSION[key]})")
        if chosen:
            self.path.setText(chosen)


class GatesReportDialog(QtWidgets.QDialog):
    """What was written, what was not and why: the export's report."""

    def __init__(self, report: Report, parent=None):
        super().__init__(parent)
        describe(self, "workspaces-and-gating-ml")
        self.setWindowTitle("Gates exported")
        self.summary = QtWidgets.QLabel(report.summary())
        self.summary.setWordWrap(True)
        self.text = QtWidgets.QPlainTextEdit(report.text())
        self.text.setReadOnly(True)
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        box.accepted.connect(self.accept)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(self.summary)
        lay.addWidget(self.text, 1)
        if report.format == WSP:
            note = QtWidgets.QLabel("Not opened in FlowJo here: open it there and compare "
                                    "the counts before relying on it.")
            note.setWordWrap(True)
            note.setProperty("role", "warning")
            lay.addWidget(note)
        lay.addWidget(box)
        self.resize(720, 420)
