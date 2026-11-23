"""
A table of results — tests between groups — to read, sort, copy and save.
"""

from __future__ import annotations

import os

from PyQt6 import QtCore, QtWidgets

from ..stats import StatsTable
from .settings import settings
from .statistics_tab import _Cell

SETTING_FOLDER = "statistics/folder"


def _text(column: str, value) -> str:
    if value is None or (isinstance(value, float) and value != value):
        return ""
    if column.startswith("p") and isinstance(value, float):
        return "< 0.0001" if value < 1e-4 else f"{value:.4f}"
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


class ResultsDialog(QtWidgets.QDialog):
    def __init__(self, table: StatsTable, title: str, note: str = "", parent=None,
                 name: str = "comparison"):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.result = table
        self.name = name
        head = QtWidgets.QLabel(note)
        head.setWordWrap(True)
        head.setProperty("role", "caption")
        self.table = QtWidgets.QTableWidget(len(table.rows), len(table.columns))
        self.table.setHorizontalHeaderLabels(table.columns)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        for r, row in enumerate(table.rows):
            for k, value in enumerate(row):
                self.table.setItem(r, k, _Cell(value, _text(table.columns[k], value)))
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        export = QtWidgets.QPushButton("Export CSV…")
        copy = QtWidgets.QPushButton("Copy")
        close = QtWidgets.QPushButton("Close")
        export.clicked.connect(self.export)
        copy.clicked.connect(self.copy)
        close.clicked.connect(self.accept)
        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(export)
        buttons.addWidget(copy)
        buttons.addStretch(1)
        buttons.addWidget(close)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(head)
        layout.addWidget(self.table, 1)
        layout.addLayout(buttons)
        self.resize(1100, 560)

    def export(self) -> None:
        folder = settings().value(SETTING_FOLDER, "") or os.path.expanduser("~")
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export", os.path.join(folder, self.name + ".csv"),
            "Comma-separated values (*.csv)")
        if path:
            self.result.to_csv(path)
            settings().setValue(SETTING_FOLDER, os.path.dirname(path))

    def copy(self) -> None:
        QtWidgets.QApplication.clipboard().setText(self.result.to_text())


NOTE = ("Each sample is one observation: the n of a group is its number of samples, "
        "not of cells. {tests} test(s), p adjusted by {correction}. A blank p is a "
        "test that could not be done; the last column says why.")
