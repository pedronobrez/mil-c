"""
The plate: the panels collected from the figure dialogs, put in order on
one page, and written as one figure — or as a report, behind a page that
says what the analysis was.

As in the figure dialog, the preview is drawn by the code that writes the
file (`plate.py`). Each panel remembers the plot it came from, so when
the analysis has changed the page can be drawn again from what those
plots show now.

A project has as many plates as the user makes, each with a name. They
are the session's (`session.plates`) and are written with the project
(`session.save_plates`) as they are: a panel read back from a project is
what it was when it was saved, and no longer follows a plot.
"""

from __future__ import annotations

import os
from dataclasses import replace

from PyQt6 import QtCore, QtGui, QtWidgets

from .. import figures, plate, plate_store
from . import style as ui_style
from .figure_dialog import SETTING_FOLDER
from .settings import settings

#: `add(to=NEW)`: on a plate made for it
NEW = "new"

WIDTHS = [("Custom", None), ("Single column, 85 mm", 85), ("1.5 columns, 120 mm", 120),
          ("Double column, 180 mm", 180)]


class PlateWindow(QtWidgets.QDialog):
    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Plate")
        self.session = session
        self._filling = False
        self._saver = QtCore.QTimer(self)
        self._saver.setSingleShot(True)
        self._saver.setInterval(1500)
        self._saver.timeout.connect(self.save)
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

        # --- the plates
        self.plates_box = QtWidgets.QComboBox()
        self.new_button = QtWidgets.QPushButton("New")
        self.rename_button = QtWidgets.QPushButton("Rename…")
        self.duplicate_button = QtWidgets.QPushButton("Duplicate")
        self.delete_button = QtWidgets.QPushButton("Delete")

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
        self.report_button.setToolTip(
            "A PDF: what the analysis was, then this page — or the plates chosen, "
            "a page each")
        self.close_button = QtWidgets.QPushButton("Close")
        buttons = QtWidgets.QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(self.close_button)
        buttons.addWidget(self.report_button)
        buttons.addWidget(self.export_button)

        side = QtWidgets.QVBoxLayout()
        side.setContentsMargins(0, 0, 0, 0)
        plates = QtWidgets.QHBoxLayout()
        for b in (self.new_button, self.rename_button, self.duplicate_button,
                  self.delete_button):
            plates.addWidget(b)
        side.addWidget(self._eyebrow("Plate"))
        side.addWidget(self.plates_box)
        side.addLayout(plates)
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
        self.plates_box.currentIndexChanged.connect(self._plate_chosen)
        self.new_button.clicked.connect(lambda: self.new_plate())
        self.rename_button.clicked.connect(self._ask_name)
        self.duplicate_button.clicked.connect(lambda: self.duplicate_plate())
        self.delete_button.clicked.connect(self._ask_delete)
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

    # --- the plates of the project ---------------------------------------------------

    @property
    def plates(self) -> list[plate.Plate]:
        return self.session.plates

    @property
    def plate(self) -> plate.Plate:
        """The plate being worked on. There is always one: the first is
        made when it is asked for."""
        s = self.session
        if not s.plates:
            s.plates.append(plate.Plate(name="Plate 1"))
        s.current_plate = min(max(s.current_plate, 0), len(s.plates) - 1)
        return s.plates[s.current_plate]

    def _free_name(self, wanted: str, but: plate.Plate | None = None) -> str:
        """`wanted`, or `wanted 2`, `wanted 3`…: two plates of one name
        cannot be told apart in a list."""
        taken = {p.name for p in self.plates if p is not but}
        name, k = wanted, 1
        while name in taken:
            k += 1
            name = f"{wanted} {k}"
        return name

    def _new_name(self) -> str:
        k = len(self.plates) + 1
        while f"Plate {k}" in {p.name for p in self.plates}:
            k += 1
        return f"Plate {k}"

    def switch_to(self, index: int) -> None:
        if 0 <= index < len(self.plates) and index != self.session.current_plate:
            self.session.current_plate = index
            self._changed()
        self._fill(0)
        self._timer.start()

    def _plate_chosen(self, index: int) -> None:
        if not self._filling:
            self.switch_to(index)

    def new_plate(self, name: str = "") -> plate.Plate:
        """An empty plate with the type and the size of the current one;
        it becomes the current one."""
        self.plate                                     # the first exists before a second
        name = self._free_name(name.strip()) if name.strip() else self._new_name()
        new = replace(self.plate.copy(name), panels=[])
        self.plates.append(new)
        self.switch_to(len(self.plates) - 1)
        return new

    def rename_plate(self, name: str) -> None:
        name = name.strip()
        if name:
            self.plate.name = self._free_name(name, self.plate)
            self._changed()
            self._fill()

    def duplicate_plate(self) -> plate.Plate:
        """A copy of the current plate, to change apart from it; it becomes
        the current one."""
        copy = self.plate.copy(self._free_name(f"{self.plate.name} copy"))
        self.plates.insert(self.session.current_plate + 1, copy)
        self.switch_to(self.session.current_plate + 1)
        return copy

    def delete_plate(self) -> None:
        """The current plate, with its panels. The one before it becomes
        the current one."""
        self.plate
        del self.plates[self.session.current_plate]
        self.session.current_plate = max(self.session.current_plate - 1, 0)
        self._changed()
        self._fill(0)
        self._timer.start()

    def _ask_name(self) -> None:
        name, ok = QtWidgets.QInputDialog.getText(
            self, "Rename the plate", "Name", text=self.plate.name)
        if ok:
            self.rename_plate(name)

    def _ask_delete(self) -> None:
        p = self.plate
        if p.panels:
            Box = QtWidgets.QMessageBox
            if Box.question(self, "Delete the plate",
                            f"Delete “{p.name}” and its {len(p.panels)} panel(s)? "
                            "This cannot be undone.") != Box.StandardButton.Yes:
                return
        self.delete_plate()

    def project_changed(self) -> None:
        """Another folder was opened: its plates are the ones to show."""
        self._saver.stop()
        self._fill(0)
        self._timer.start()

    # --- keeping them ------------------------------------------------------------------

    def _changed(self) -> None:
        """Something of a plate is no longer as the project has it: it is
        written once the changes pause, and at the latest when the window
        closes."""
        self._saver.start()

    def save(self) -> None:
        self._saver.stop()
        try:
            self.session.save_plates()
        except Exception as exc:
            self.note.setText(f"The plates were not saved with the project: {exc}")

    def flush(self) -> None:
        if self._saver.isActive():
            self.save()

    def hideEvent(self, ev) -> None:
        self.flush()
        super().hideEvent(ev)

    # --- plate ↔ widgets ---------------------------------------------------------

    def _line(self, k: int) -> str:
        panel = self.plate.panels[k]
        text = f"{plate.letter(k)}   {panel.name or 'panel'}"
        if panel.title:
            text += f" — {panel.title}"
        if isinstance(panel.content, plate_store.Missing):
            text += "   (not restored)"
        elif panel.stale:
            text += "   (stale: its plot has nothing to show)"
        return text

    def _fill(self, row: int | None = None) -> None:
        """Everything from the plate; `row` is the panel to leave selected."""
        p = self.plate
        row = self.list.currentRow() if row is None else row
        self._filling = True
        self.plates_box.clear()
        self.plates_box.addItems([q.name for q in self.plates])
        self.plates_box.setCurrentIndex(self.session.current_plate)
        self.setWindowTitle(f"Plate — {p.name}")
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
        self._changed()
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
        self._changed()
        self._timer.start()

    def _width_chosen(self) -> None:
        mm = self.width_preset.currentData()
        if mm and not self._filling:
            self.width.setValue(mm)

    # --- what is done to the plate -------------------------------------------------

    def add(self, content, style: figures.FigureStyle, source=None, name: str = "",
            to: int | str | None = None) -> str:
        """A panel as `style` draws it, on the current plate — or on the
        plate numbered `to`, or on a new one (`NEW`), which then becomes
        the current one; the letter it was given. The first panel brings
        its type and resolution to the page: the ones after take the
        page's. Written with the project at once."""
        if to == NEW:
            self.new_plate()
        elif to is not None and 0 <= to < len(self.plates):
            self.session.current_plate = to
        p = self.plate
        if not p.panels:
            p.style.font, p.style.font_size, p.style.dpi = style.font, style.font_size, style.dpi
        p.add(content, style, source=source, name=name)
        self._fill(len(p.panels) - 1)
        self._timer.start()
        self.save()
        return plate.letter(len(p.panels) - 1)

    def move(self, by: int) -> None:
        row = self.list.currentRow()
        if row >= 0:
            self._fill(self.plate.move(row, by))
            self._changed()
            self._timer.start()

    def remove(self) -> None:
        row = self.list.currentRow()
        if row >= 0:
            del self.plate.panels[row]
            self._fill(row)
            self._changed()
            self._timer.start()

    def refresh_from_analysis(self) -> None:
        stale = self.plate.refresh()
        fixed = sum(p.source is None for p in self.plate.panels)
        self._fill()
        self.render()
        self._changed()
        said = []
        if stale:
            said.append(f"{len(stale)} panel(s) kept as they were: their plots "
                        "have nothing to show now.")
        if fixed:
            said.append(f"{fixed} panel(s) came with the project and stay as they were "
                        "saved: add them again from their plots to follow the analysis.")
        if said:
            self.note.setText(" ".join(said))

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

    def write_report(self, path: str, plates: list[plate.Plate] | None = None) -> int:
        """What the analysis was, where there is one, then the page — or
        `plates`, a page each, in the order given."""
        pages = [self.plate] if plates is None else plates
        return plate.report(path, pages, self.session.dataset)

    def _ask_plates(self) -> list[plate.Plate] | None:
        """Which plates the report is to have; None when it was cancelled.
        Not asked when only one has anything on it."""
        full = [p for p in self.plates if p.panels]
        if len(full) < 2:
            return [self.plate]
        dialog = QtWidgets.QDialog(self)
        dialog.setWindowTitle("Report")
        boxes = []
        layout = QtWidgets.QVBoxLayout(dialog)
        layout.addWidget(QtWidgets.QLabel("Plates in the report, a page each:"))
        for p in full:
            box = QtWidgets.QCheckBox(f"{p.name} — {len(p.panels)} panel(s)")
            box.setChecked(True)
            boxes.append(box)
            layout.addWidget(box)
        ok = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                        | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        ok.accepted.connect(dialog.accept)
        ok.rejected.connect(dialog.reject)
        layout.addWidget(ok)
        if dialog.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return None
        return [p for p, box in zip(full, boxes) if box.isChecked()]

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
        plates = self._ask_plates()
        if not plates:
            return
        path, _ = self._ask("Write the report", "report.pdf", "PDF (*.pdf)")
        if path:
            self._written(path if path.lower().endswith(".pdf") else path + ".pdf",
                          lambda p: self.write_report(p, plates))
