"""
The plate: the panels collected from the figure dialogs, put in order on
one page, and written as one figure — or as a report, behind a page that
says what the analysis was.

As in the figure dialog, the preview is drawn by the code that writes the
file (`plate.py`). Each panel remembers the plot it came from, so when
the analysis has changed the page can be drawn again from what those
plots show now. The plate lasts as long as the window: it is not written
to the project, because a panel is a plot's events and the plots are made
again from the analysis anyway.
"""

from __future__ import annotations

import os

from PyQt6 import QtCore, QtGui, QtWidgets

from .. import figures, plate
from . import style as ui_style
from .figure_dialog import SETTING_FOLDER
from .settings import settings

WIDTHS = [("Custom", None), ("Single column, 85 mm", 85), ("1.5 columns, 120 mm", 120),
          ("Double column, 180 mm", 180)]


class PlateWindow(QtWidgets.QDialog):
    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Plate")
        self.session = session
        self.plate = plate.Plate()
        self._filling = False
        self._timer = QtCore.QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(180)
        # a plate nobody is looking at is drawn when it is shown
        self._timer.timeout.connect(lambda: self.isVisible() and self.render())

        # --- preview
        self.preview = QtWidgets.QLabel()
        self.preview.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.preview.setMinimumSize(200, 200)
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidget(self.preview)
        self.scroll.setWidgetResizable(True)
        self.note = QtWidgets.QLabel("")
        self.note.setProperty("role", "caption")
        self.note.setWordWrap(True)

        # --- the panels
        self.list = QtWidgets.QListWidget()
        self.up_button = QtWidgets.QPushButton("Move up")
        self.down_button = QtWidgets.QPushButton("Move down")
        self.remove_button = QtWidgets.QPushButton("Remove")
        self.title = QtWidgets.QLineEdit()
        self.title.setPlaceholderText("none")
        self.colspan = self._count(1, 6, " column(s)")
        self.rowspan = self._count(1, 6, " row(s)")
        self.refresh_button = QtWidgets.QPushButton("Refresh from the current analysis")
        self.refresh_button.setToolTip(
            "Draw each panel again from the plot it came from, as that plot is now. "
            "A panel whose plot has nothing to show keeps what it had and is marked stale.")

        # --- the page
        self.columns = self._count(1, 6, "")
        self.width_preset = QtWidgets.QComboBox()
        for text, mm in WIDTHS:
            self.width_preset.addItem(text, mm)
        self.width = self._count(20, 600, " mm")
        self.height = self._count(0, 900, " mm")
        self.height.setSpecialValueText("Automatic")
        self.height.setToolTip("Automatic: three quarters of a cell's width for each row")
        self.letters = QtWidgets.QCheckBox("Letters: A, B, C…")
        self.font = QtWidgets.QComboBox()
        self.font.addItems(figures.fonts())
        self.font_size = QtWidgets.QDoubleSpinBox()
        self.font_size.setRange(4, 30)
        self.font_size.setSingleStep(0.5)
        self.font_size.setDecimals(1)
        self.font_size.setSuffix(" pt")
        self.dpi = self._count(72, 1200, " dpi")
        self.dpi.setSingleStep(50)

        moves = QtWidgets.QHBoxLayout()
        for b in (self.up_button, self.down_button, self.remove_button):
            moves.addWidget(b)
        spans = QtWidgets.QHBoxLayout()
        spans.addWidget(self.colspan)
        spans.addWidget(QtWidgets.QLabel("×"))
        spans.addWidget(self.rowspan)
        size = QtWidgets.QHBoxLayout()
        size.addWidget(self.width)
        size.addWidget(QtWidgets.QLabel("×"))
        size.addWidget(self.height)
        form = QtWidgets.QFormLayout()
        form.setFieldGrowthPolicy(
            QtWidgets.QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        form.addRow(self._eyebrow("The panel selected"))
        form.addRow("Title", self.title)
        form.addRow("Takes", spans)
        form.addRow(self._eyebrow("Page"))
        form.addRow("Panels per row", self.columns)
        form.addRow("Fit", self.width_preset)
        form.addRow("Width × height", size)
        form.addRow(self.letters)
        form.addRow("Font", self.font)
        form.addRow("Size", self.font_size)
        form.addRow("Resolution", self.dpi)

        self.export_button = QtWidgets.QPushButton("Export…")
        self.export_button.setProperty("primary", True)
        self.report_button = QtWidgets.QPushButton("Report…")
        self.report_button.setToolTip("A PDF: what the analysis was, then this page")
        self.close_button = QtWidgets.QPushButton("Close")
        buttons = QtWidgets.QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(self.close_button)
        buttons.addWidget(self.report_button)
        buttons.addWidget(self.export_button)

        side = QtWidgets.QVBoxLayout()
        side.setContentsMargins(0, 0, 0, 0)
        side.addWidget(self._eyebrow("Panels, in order"))
        side.addWidget(self.list, 1)
        side.addLayout(moves)
        side.addWidget(self.refresh_button)
        side.addLayout(form)
        side.addLayout(buttons)
        panel = QtWidgets.QWidget()
        panel.setLayout(side)
        panel.setFixedWidth(420)
        left = QtWidgets.QVBoxLayout()
        left.addWidget(self.scroll, 1)
        left.addWidget(self.note)
        layout = QtWidgets.QHBoxLayout(self)
        layout.addLayout(left, 1)
        layout.addWidget(panel)
        self.resize(1220, 800)

        self._fill()
        self.list.currentRowChanged.connect(self._selected)
        self.up_button.clicked.connect(lambda: self.move(-1))
        self.down_button.clicked.connect(lambda: self.move(1))
        self.remove_button.clicked.connect(self.remove)
        self.refresh_button.clicked.connect(self.refresh_from_analysis)
        self.title.textChanged.connect(self._panel_changed)
        for w in (self.colspan, self.rowspan):
            w.valueChanged.connect(self._panel_changed)
        for w in (self.columns, self.width, self.height, self.font_size, self.dpi):
            w.valueChanged.connect(self._page_changed)
        self.font.currentIndexChanged.connect(self._page_changed)
        self.letters.toggled.connect(self._page_changed)
        self.width_preset.currentIndexChanged.connect(self._width_chosen)
        self.export_button.clicked.connect(self.export)
        self.report_button.clicked.connect(self.export_report)
        self.close_button.clicked.connect(self.close)

    @staticmethod
    def _count(lo, hi, suffix):
        s = QtWidgets.QSpinBox()
        s.setRange(lo, hi)
        s.setSuffix(suffix)
        return s

    @staticmethod
    def _eyebrow(text):
        label = QtWidgets.QLabel(text)
        label.setFont(ui_style.eyebrow_font())
        label.setProperty("role", "caption")
        label.setContentsMargins(0, 8, 0, 0)
        return label

    # --- plate ↔ widgets ---------------------------------------------------------

    def _line(self, k: int) -> str:
        panel = self.plate.panels[k]
        text = f"{plate.letter(k)}   {panel.name or 'panel'}"
        if panel.title:
            text += f" — {panel.title}"
        if panel.stale:
            text += "   (stale: its plot has nothing to show)"
        return text

    def _fill(self, row: int | None = None) -> None:
        """Everything from the plate; `row` is the panel to leave selected."""
        p = self.plate
        row = self.list.currentRow() if row is None else row
        self._filling = True
        self.list.clear()
        for k in range(len(p.panels)):
            self.list.addItem(self._line(k))
        self.columns.setValue(p.columns)
        self.width.setValue(round(p.width_mm))
        self.width_preset.setCurrentIndex(max(self.width_preset.findData(round(p.width_mm)), 0))
        self.height.setValue(round(p.height_mm or 0))
        self.letters.setChecked(p.letters)
        i = self.font.findText(p.style.font)
        self.font.setCurrentIndex(i if i >= 0 else max(self.font.findText("DejaVu Sans"), 0))
        self.font_size.setValue(p.style.font_size)
        self.dpi.setValue(p.style.dpi)
        self._filling = False
        self.list.setCurrentRow(min(max(row, 0), len(p.panels) - 1))
        self._selected()

    def _selected(self, *_) -> None:
        if self._filling:
            return
        panel = self.selected()
        self._filling = True
        for w in (self.title, self.colspan, self.rowspan, self.up_button,
                  self.down_button, self.remove_button):
            w.setEnabled(panel is not None)
        self.title.setText(panel.title if panel else "")
        self.colspan.setValue(panel.colspan if panel else 1)
        self.rowspan.setValue(panel.rowspan if panel else 1)
        self._filling = False

    def selected(self) -> plate.Panel | None:
        row = self.list.currentRow()
        return self.plate.panels[row] if 0 <= row < len(self.plate.panels) else None

    def _panel_changed(self, *_) -> None:
        panel = self.selected()
        if self._filling or panel is None:
            return
        panel.title = self.title.text()
        panel.colspan, panel.rowspan = self.colspan.value(), self.rowspan.value()
        self.list.currentItem().setText(self._line(self.list.currentRow()))
        self._timer.start()

    def _page_changed(self, *_) -> None:
        if self._filling:
            return
        p = self.plate
        p.columns = self.columns.value()
        p.width_mm = float(self.width.value())
        p.height_mm = float(self.height.value()) or None
        p.letters = self.letters.isChecked()
        p.style.font, p.style.font_size = self.font.currentText(), self.font_size.value()
        p.style.dpi = self.dpi.value()
        self._filling = True
        self.width_preset.setCurrentIndex(max(self.width_preset.findData(self.width.value()), 0))
        self._filling = False
        self._timer.start()

    def _width_chosen(self) -> None:
        mm = self.width_preset.currentData()
        if mm and not self._filling:
            self.width.setValue(mm)

    # --- what is done to the plate -------------------------------------------------

    def add(self, content, style: figures.FigureStyle, source=None, name: str = "") -> str:
        """A panel as `style` draws it; the letter it was given. The first
        panel brings its type and resolution to the page: the ones after
        take the page's."""
        p = self.plate
        if not p.panels:
            p.style.font, p.style.font_size, p.style.dpi = style.font, style.font_size, style.dpi
        p.add(content, style, source=source, name=name)
        self._fill(len(p.panels) - 1)
        self._timer.start()
        return plate.letter(len(p.panels) - 1)

    def move(self, by: int) -> None:
        row = self.list.currentRow()
        if row >= 0:
            self._fill(self.plate.move(row, by))
            self._timer.start()

    def remove(self) -> None:
        row = self.list.currentRow()
        if row >= 0:
            del self.plate.panels[row]
            self._fill(row)
            self._timer.start()

    def refresh_from_analysis(self) -> None:
        stale = self.plate.refresh()
        self._fill()
        self.render()
        if stale:
            self.note.setText(f"{len(stale)} panel(s) kept as they were: their plots "
                              "have nothing to show now.")

    # --- drawing and writing ---------------------------------------------------------

    def render(self) -> None:
        if not self.plate.panels:
            self.preview.setPixmap(QtGui.QPixmap())
            self.note.setText("No panels yet: open a plot's Figure… and choose Add to plate.")
            return
        try:
            data = plate.preview_png(self.plate)
        except Exception as exc:
            self.note.setText(f"Cannot draw this: {exc}")
            return
        pix = QtGui.QPixmap()
        pix.loadFromData(data, "PNG")
        room = self.scroll.viewport().size() - QtCore.QSize(8, 8)
        # (a window not yet shown has no room to speak of: nothing to fit to)
        if room.isValid() and (pix.width() > room.width() or pix.height() > room.height()):
            pix = pix.scaled(room, QtCore.Qt.AspectRatioMode.KeepAspectRatio,
                             QtCore.Qt.TransformationMode.SmoothTransformation)
        self.preview.setPixmap(pix)
        w, h = self.plate.size_mm()
        dpi = self.plate.style.dpi
        self.note.setText(
            f"{w:g} × {h:g} mm — as PNG or TIFF at {dpi} dpi: "
            f"{round(w / 25.4 * dpi):,d} × {round(h / 25.4 * dpi):,d} px")

    def showEvent(self, ev) -> None:
        super().showEvent(ev)
        self._timer.start()

    def resizeEvent(self, ev) -> None:
        super().resizeEvent(ev)
        self._timer.start()

    def write(self, path: str) -> None:
        plate.save(self.plate, path)

    def write_report(self, path: str) -> int:
        """What the analysis was, where there is one, then the page."""
        return plate.report(path, [self.plate], self.session.dataset)

    def _ask(self, title: str, name: str, filters: str) -> tuple[str, str]:
        folder = settings().value(SETTING_FOLDER, "") or os.path.expanduser("~")
        return QtWidgets.QFileDialog.getSaveFileName(
            self, title, os.path.join(folder, name), filters)

    def _written(self, path: str, write) -> None:
        try:
            write(path)
        except Exception as exc:
            self.note.setText(f"Not written: {exc}")
            return
        settings().setValue(SETTING_FOLDER, os.path.dirname(path))
        self.note.setText(f"Wrote {path}")

    def export(self) -> None:
        path, chosen = self._ask("Export the plate", "plate.svg",
                                 ";;".join(figures.FORMATS.values()))
        if not path:
            return
        if "." not in os.path.basename(path):
            path += "." + next(k for k, v in figures.FORMATS.items() if v == chosen)
        self._written(path, self.write)

    def export_report(self) -> None:
        path, _ = self._ask("Write the report", "report.pdf", "PDF (*.pdf)")
        if path:
            self._written(path if path.lower().endswith(".pdf") else path + ".pdf",
                          self.write_report)
