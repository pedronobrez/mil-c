"""
A plot, made into a figure: its size in millimetres, its type, its words,
its colours and its legend — then SVG, PDF, PNG or TIFF.

The preview is drawn by the code that writes the file (`figures.py`), so
what is seen is what is written. The look is remembered, and can be kept
under a name: whoever exports one figure for a manuscript is about to
export six that must match.
"""

from __future__ import annotations

import json
import os

from PyQt6 import QtCore, QtGui, QtWidgets

from .. import figures
from ..figures import FigureStyle
from . import palette
from . import style as ui_style
from .ramp_picker import RampPicker
from .settings import settings

SETTING_LAST = "figure/last"
SETTING_PRESETS = "figure/presets"
SETTING_FOLDER = "figure/folder"

FRAMES = [("open", "Left and bottom"), ("box", "Box"), ("none", "None")]
LEGENDS = [("right", "Right"), ("bottom", "Below"), ("inside", "Inside"), ("none", "Hidden")]
SIZES = [("Custom", None), ("Single column, 85 × 65 mm", (85, 65)),
         ("1.5 columns, 120 × 90 mm", (120, 90)),
         ("Double column, 180 × 110 mm", (180, 110)),
         ("Square, 85 × 85 mm", (85, 85))]


def load_presets() -> dict[str, dict]:
    try:
        return json.loads(settings().value(SETTING_PRESETS, "{}") or "{}")
    except (TypeError, ValueError):
        return {}


def last_style() -> FigureStyle:
    try:
        return FigureStyle.from_dict(json.loads(settings().value(SETTING_LAST, "{}") or "{}"))
    except (TypeError, ValueError):
        return FigureStyle()


class FigureDialog(QtWidgets.QDialog):
    def __init__(self, content, parent=None, name: str = "figure"):
        super().__init__(parent)
        self.setWindowTitle("Figure")
        self.content = content
        self.name = name
        self.style = last_style()
        # the look is remembered; this figure's own words are not
        self.style.title, self.style.xlabel, self.style.ylabel = "", None, None
        self.style.legend_title = None
        # a figure starts from the ramp on the screen: what was being looked
        # at is what is expected on paper
        self.style.colormap = palette.ramp_name()
        self.style.colormap_colors = palette.custom_colors()
        self._filling = False
        self._timer = QtCore.QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(180)
        self._timer.timeout.connect(self.render)

        # --- preview
        self.preview = QtWidgets.QLabel()
        self.preview.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.preview.setMinimumSize(200, 200)
        scroll = self.scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.preview)
        scroll.setWidgetResizable(True)
        self.note = QtWidgets.QLabel("")
        self.note.setProperty("role", "caption")

        # --- controls
        self.preset = QtWidgets.QComboBox()
        self.save_preset = QtWidgets.QPushButton("Save as…")
        self.delete_preset = QtWidgets.QPushButton("Delete")
        self.size_preset = QtWidgets.QComboBox()
        for text, wh in SIZES:
            self.size_preset.addItem(text, wh)
        self.width = self._spin(20, 600, " mm")
        self.height = self._spin(20, 600, " mm")
        self.dpi = QtWidgets.QSpinBox()
        self.dpi.setRange(72, 1200)
        self.dpi.setSingleStep(50)
        self.dpi.setSuffix(" dpi")
        self.font = QtWidgets.QComboBox()
        self.font.addItems(figures.fonts())
        self.font_size = self._spin(4, 30, " pt", 0.5)
        self.title = QtWidgets.QLineEdit()
        self.xlabel = QtWidgets.QLineEdit()
        self.ylabel = QtWidgets.QLineEdit()
        for edit in (self.xlabel, self.ylabel):
            edit.setPlaceholderText("as on the plot")
        self.frame = self._combo(FRAMES)
        self.ticks = QtWidgets.QCheckBox("Ticks")
        self.grid = QtWidgets.QCheckBox("Grid")
        self.point_size = self._spin(0.2, 20, " pt", 0.25)
        self.opacity = self._spin(0.05, 1.0, "", 0.05)
        self.rasterize = QtWidgets.QCheckBox("Points as a bitmap inside vector files")
        self.rasterize.setToolTip("For maps of very many events: text and axes stay "
                                  "vector, the cloud of points becomes an image.")
        self.legend = self._combo(LEGENDS)
        self.legend_title = QtWidgets.QLineEdit()
        self.legend_title.setPlaceholderText("as on the plot")
        self.legend_counts = QtWidgets.QCheckBox("Counts")
        self.colormap = RampPicker(follow=False)
        self.annotate = QtWidgets.QCheckBox("Numbers in the cells")
        self.transparent = QtWidgets.QCheckBox("Transparent background")
        self.categories = QtWidgets.QTableWidget(0, 2)
        self.categories.setHorizontalHeaderLabels(["Colour", "Label"])
        self.categories.verticalHeader().setVisible(False)
        self.categories.horizontalHeader().setStretchLastSection(True)
        self.categories.setColumnWidth(0, 70)
        self.categories.setIconSize(QtCore.QSize(44, 18))
        self.categories.setToolTip("Double-click a colour to change it; type over a label")
        self.categories_title = self._eyebrow("Categories — double-click a colour, type over a label")

        form = QtWidgets.QFormLayout()
        form.setFieldGrowthPolicy(
            QtWidgets.QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        form.setContentsMargins(0, 0, 14, 0)   # clear of the scroll bar
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.preset, 1)
        row.addWidget(self.save_preset)
        row.addWidget(self.delete_preset)
        form.addRow("Preset", row)
        form.addRow(self._eyebrow("Size"))
        form.addRow("Fit", self.size_preset)
        wh = QtWidgets.QHBoxLayout()
        wh.addWidget(self.width)
        wh.addWidget(QtWidgets.QLabel("×"))
        wh.addWidget(self.height)
        form.addRow("Width × height", wh)
        form.addRow("Resolution", self.dpi)
        form.addRow(self._eyebrow("Text"))
        form.addRow("Font", self.font)
        form.addRow("Size", self.font_size)
        form.addRow("Title", self.title)
        self._row_x = self._add(form, "X label", self.xlabel)
        self._row_y = self._add(form, "Y label", self.ylabel)
        form.addRow(self._eyebrow("Axes and points"))
        fr = QtWidgets.QHBoxLayout()
        fr.addWidget(self.frame, 1)
        fr.addWidget(self.ticks)
        fr.addWidget(self.grid)
        form.addRow("Frame", fr)
        self._row_point = self._add(form, "Point size", self.point_size)
        self._row_opacity = self._add(form, "Opacity", self.opacity)
        form.addRow(self.rasterize)
        form.addRow(self._eyebrow("Legend and colours"))
        lg = QtWidgets.QHBoxLayout()
        lg.addWidget(self.legend, 1)
        lg.addWidget(self.legend_counts)
        form.addRow("Legend", lg)
        form.addRow("Legend title", self.legend_title)
        self._row_cmap = self._add(form, "Colour ramp", self.colormap)
        form.addRow(self.annotate)
        form.addRow(self.transparent)

        self.export_button = QtWidgets.QPushButton("Export…")
        self.export_button.setProperty("primary", True)
        self.close_button = QtWidgets.QPushButton("Close")
        buttons = QtWidgets.QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(self.close_button)
        buttons.addWidget(self.export_button)

        # everything but the buttons scrolls: on a laptop screen the list of
        # categories was squeezed to one row under the form and could not be
        # reached at all
        inner = QtWidgets.QVBoxLayout()
        inner.addLayout(form)
        inner.addWidget(self.categories_title)
        inner.addWidget(self.categories)
        inner.addStretch(1)
        body = QtWidgets.QWidget()
        body.setLayout(inner)
        side_scroll = QtWidgets.QScrollArea()
        side_scroll.setWidget(body)
        side_scroll.setWidgetResizable(True)
        side_scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        side_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        side = QtWidgets.QVBoxLayout()
        side.setContentsMargins(0, 0, 0, 0)
        side.addWidget(side_scroll, 1)
        side.addLayout(buttons)
        panel = QtWidgets.QWidget()
        panel.setLayout(side)
        panel.setFixedWidth(480)
        left = QtWidgets.QVBoxLayout()
        left.addWidget(scroll, 1)
        left.addWidget(self.note)
        layout = QtWidgets.QHBoxLayout(self)
        layout.addLayout(left, 1)
        layout.addWidget(panel)
        self.resize(1220, 800)

        # what this content has no use for is not offered
        is_map = isinstance(content, figures.MapContent)
        is_tree = isinstance(content, figures.TreeContent)
        is_heat = isinstance(content, figures.HeatmapContent)
        self._visible(self._row_point, not is_heat)
        self._visible(self._row_opacity, is_map)
        self.rasterize.setVisible(is_map)
        self._visible(self._row_cmap, figures.uses_colormap(content))
        self.annotate.setVisible(is_heat)
        self._visible(self._row_x, not is_tree)
        self._visible(self._row_y, not is_tree)
        self.grid.setVisible(not (is_tree or is_heat))
        self.ticks.setVisible(not (is_tree or is_heat))
        has_categories = bool(figures.category_names(content))
        self.categories.setVisible(has_categories)
        self.categories_title.setVisible(has_categories)

        self._fill_presets()
        self._to_widgets()
        for w in (self.width, self.height, self.font_size, self.point_size, self.opacity):
            w.valueChanged.connect(self._changed)
        self.dpi.valueChanged.connect(self._changed)
        for w in (self.font, self.frame, self.legend):
            w.currentIndexChanged.connect(self._changed)
        self.colormap.sigChosen.connect(self._changed)
        for w in (self.title, self.xlabel, self.ylabel, self.legend_title):
            w.textEdited.connect(self._changed)
        for w in (self.ticks, self.grid, self.rasterize, self.legend_counts,
                  self.annotate, self.transparent):
            w.toggled.connect(self._changed)
        self.size_preset.currentIndexChanged.connect(self._size_chosen)
        self.preset.activated.connect(self._preset_chosen)
        self.save_preset.clicked.connect(self._save_preset)
        self.delete_preset.clicked.connect(self._delete_preset)
        self.categories.itemChanged.connect(self._category_edited)
        self.categories.cellDoubleClicked.connect(self._pick_colour)
        self.export_button.clicked.connect(self.export)
        self.close_button.clicked.connect(self.accept)
        self.render()

    # --- small builders ----------------------------------------------------------

    @staticmethod
    def _spin(lo, hi, suffix, step=1.0):
        s = QtWidgets.QDoubleSpinBox()
        s.setRange(lo, hi)
        s.setSingleStep(step)
        s.setDecimals(2 if step < 1 else 0)
        s.setSuffix(suffix)
        return s

    @staticmethod
    def _combo(items):
        c = QtWidgets.QComboBox()
        for key, text in items:
            c.addItem(text, key)
        return c

    @staticmethod
    def _eyebrow(text):
        label = QtWidgets.QLabel(text)
        label.setFont(ui_style.eyebrow_font())
        label.setProperty("role", "caption")
        label.setContentsMargins(0, 8, 0, 0)
        return label

    @staticmethod
    def _add(form, text, widget):
        label = QtWidgets.QLabel(text)
        form.addRow(label, widget)
        return label, widget

    @staticmethod
    def _visible(row, on: bool):
        for w in row:
            w.setVisible(on)

    # --- style ↔ widgets ---------------------------------------------------------

    def _to_widgets(self) -> None:
        s = self.style
        self._filling = True
        self.width.setValue(s.width_mm)
        self.height.setValue(s.height_mm)
        self.dpi.setValue(s.dpi)
        i = self.font.findText(s.font)
        self.font.setCurrentIndex(i if i >= 0 else max(self.font.findText("DejaVu Sans"), 0))
        self.font_size.setValue(s.font_size)
        self.title.setText(s.title)
        self.xlabel.setText(s.xlabel or "")
        self.ylabel.setText(s.ylabel or "")
        self.frame.setCurrentIndex(max(self.frame.findData(s.frame), 0))
        self.ticks.setChecked(s.ticks)
        self.grid.setChecked(s.grid)
        self.point_size.setValue(s.point_size)
        self.opacity.setValue(s.opacity)
        self.rasterize.setChecked(s.rasterize)
        self.legend.setCurrentIndex(max(self.legend.findData(s.legend), 0))
        self.legend_title.setText(s.legend_title or "")
        self.legend_counts.setChecked(s.legend_counts)
        self.colormap.set_choice(s.colormap, s.colormap_colors)
        self.annotate.setChecked(s.annotate)
        self.transparent.setChecked(s.transparent)
        names = figures.category_names(self.content)
        self.categories.setRowCount(len(names))
        for r, (name, default) in enumerate(names):
            swatch = QtWidgets.QTableWidgetItem()
            swatch.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
            # as a picture, not as the cell's background: the stylesheet
            # paints every cell itself and the colours never showed
            colour = s.color(name, default)
            pix = QtGui.QPixmap(44, 18)
            pix.fill(QtGui.QColor(colour))
            swatch.setIcon(QtGui.QIcon(pix))
            swatch.setData(QtCore.Qt.ItemDataRole.UserRole, name)
            swatch.setData(QtCore.Qt.ItemDataRole.UserRole + 1, colour)
            self.categories.setItem(r, 0, swatch)
            self.categories.setItem(r, 1, QtWidgets.QTableWidgetItem(s.label(name)))
        # as tall as its rows, up to a point: every category in sight
        rows = self.categories.rowCount()
        self.categories.setFixedHeight(min(34 + 30 * max(rows, 1), 34 + 30 * 12))
        self._filling = False

    def _from_widgets(self) -> None:
        s = self.style
        s.width_mm, s.height_mm = self.width.value(), self.height.value()
        s.dpi = self.dpi.value()
        s.font, s.font_size = self.font.currentText(), self.font_size.value()
        s.title = self.title.text()
        s.xlabel = self.xlabel.text() or None
        s.ylabel = self.ylabel.text() or None
        s.frame = self.frame.currentData()
        s.ticks, s.grid = self.ticks.isChecked(), self.grid.isChecked()
        s.point_size, s.opacity = self.point_size.value(), self.opacity.value()
        s.rasterize = self.rasterize.isChecked()
        s.legend = self.legend.currentData()
        s.legend_title = self.legend_title.text() or None
        s.legend_counts = self.legend_counts.isChecked()
        s.colormap = self.colormap.name()
        s.colormap_colors = list(self.colormap.colors)
        s.annotate = self.annotate.isChecked()
        s.transparent = self.transparent.isChecked()

    def _changed(self, *_):
        if self._filling:
            return
        self._from_widgets()
        self._timer.start()

    def _size_chosen(self) -> None:
        wh = self.size_preset.currentData()
        if wh:
            self._filling = True
            self.width.setValue(wh[0])
            self.height.setValue(wh[1])
            self._filling = False
            self._changed()

    # --- categories ----------------------------------------------------------------

    def _category_edited(self, item) -> None:
        if self._filling or item.column() != 1:
            return
        name = self.categories.item(item.row(), 0).data(QtCore.Qt.ItemDataRole.UserRole)
        text = item.text().strip()
        if text and text != name:
            self.style.labels[name] = text
        else:
            self.style.labels.pop(name, None)
        self._timer.start()

    def _pick_colour(self, row: int, col: int) -> None:
        if col != 0:
            return
        item = self.categories.item(row, 0)
        current = QtGui.QColor(item.data(QtCore.Qt.ItemDataRole.UserRole + 1))
        chosen = QtWidgets.QColorDialog.getColor(current, self, "Colour")
        if chosen.isValid():
            self.set_colour(item.data(QtCore.Qt.ItemDataRole.UserRole), chosen.name())

    def set_colour(self, name: str, colour: str) -> None:
        self.style.colors[name] = colour
        self._to_widgets()
        self._timer.start()

    # --- presets ---------------------------------------------------------------------

    def _fill_presets(self) -> None:
        self.preset.clear()
        self.preset.addItem("Last used", "")
        self.preset.addItem("Defaults", "*")
        for name in sorted(load_presets()):
            self.preset.addItem(name, name)

    def _preset_chosen(self) -> None:
        key = self.preset.currentData()
        if key == "":
            return
        look = {} if key == "*" else load_presets().get(key, {})
        self.apply_look(look)

    def apply_look(self, look: dict) -> None:
        """A preset changes the look and leaves this figure's words alone."""
        fresh = FigureStyle.from_dict(look)
        for name in figures.PRESET_FIELDS:
            setattr(self.style, name, getattr(fresh, name))
        self._to_widgets()
        self.render()

    def store_preset(self, name: str) -> None:
        presets = load_presets()
        full = self.style.to_dict()
        presets[name] = {k: full[k] for k in figures.PRESET_FIELDS}
        settings().setValue(SETTING_PRESETS, json.dumps(presets))
        self._fill_presets()
        self.preset.setCurrentIndex(self.preset.findData(name))

    def _save_preset(self) -> None:
        name, ok = QtWidgets.QInputDialog.getText(self, "Save the look", "Name")
        if ok and name.strip():
            self.store_preset(name.strip())

    def _delete_preset(self) -> None:
        key = self.preset.currentData()
        presets = load_presets()
        if key in presets:
            del presets[key]
            settings().setValue(SETTING_PRESETS, json.dumps(presets))
            self._fill_presets()

    # --- drawing and writing ---------------------------------------------------------

    def render(self) -> None:
        try:
            data = figures.preview_png(self.content, self.style)
        except Exception as exc:
            self.note.setText(f"Cannot draw this: {exc}")
            return
        pix = QtGui.QPixmap()
        pix.loadFromData(data, "PNG")
        # the whole figure, always: a preview cut at the edge hides the
        # legend, which is half of what is being adjusted
        room = self.scroll.viewport().size() - QtCore.QSize(8, 8)
        if pix.width() > room.width() or pix.height() > room.height():
            pix = pix.scaled(room, QtCore.Qt.AspectRatioMode.KeepAspectRatio,
                             QtCore.Qt.TransformationMode.SmoothTransformation)
        self.preview.setPixmap(pix)
        s = self.style
        self.note.setText(
            f"{s.width_mm:g} × {s.height_mm:g} mm — as PNG or TIFF at {s.dpi} dpi: "
            f"{round(s.width_mm / 25.4 * s.dpi):,d} × {round(s.height_mm / 25.4 * s.dpi):,d} px")

    def resizeEvent(self, ev) -> None:
        super().resizeEvent(ev)
        self._timer.start()

    def remember(self) -> None:
        settings().setValue(SETTING_LAST, json.dumps(self.style.to_dict()))

    def write(self, path: str) -> None:
        self._from_widgets()
        figures.save(self.content, self.style, path)
        self.remember()

    def export(self) -> None:
        folder = settings().value(SETTING_FOLDER, "") or os.path.expanduser("~")
        path, chosen = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export the figure", os.path.join(folder, self.name + ".svg"),
            ";;".join(figures.FORMATS.values()))
        if not path:
            return
        if "." not in os.path.basename(path):
            ext = next(k for k, v in figures.FORMATS.items() if v == chosen)
            path += "." + ext
        try:
            self.write(path)
        except Exception as exc:
            self.note.setText(f"Not written: {exc}")
            return
        settings().setValue(SETTING_FOLDER, os.path.dirname(path))
        self.note.setText(f"Wrote {path}")

    def accept(self) -> None:
        self._from_widgets()
        self.remember()
        super().accept()


def open_for(view, parent, name: str = "figure") -> FigureDialog | None:
    """The dialog for whatever `view` (anything with `content()`) shows
    now; None when it shows nothing."""
    content = view.content()
    if content is None:
        return None
    dialog = FigureDialog(content, parent, name)
    dialog.exec()
    return dialog
