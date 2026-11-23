"""
Choosing the colour ramp: the named ones, each shown as a strip of
itself, and "Custom…" — two or three colours of the user's own.
"""

from __future__ import annotations

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets

from .. import figures
from . import palette


def _strip(name: str, colors=None, width: int = 72, height: int = 12) -> QtGui.QIcon:
    lut = figures.lookup_table(name, colors, width)
    row = np.ascontiguousarray(np.repeat(lut[None, :, :], height, axis=0))
    image = QtGui.QImage(row.data, width, height, width * 4,
                         QtGui.QImage.Format.Format_RGBA8888).copy()
    return QtGui.QIcon(QtGui.QPixmap.fromImage(image))


class CustomRampDialog(QtWidgets.QDialog):
    """Low, middle (optional) and high."""

    def __init__(self, colors: list[str], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Your own colour ramp")
        colors = list(colors) or list(figures.DEFAULT_CUSTOM)
        self.low = colors[0]
        self.high = colors[-1]
        self.middle = colors[1] if len(colors) > 2 else ""
        self.buttons = {}
        self.use_middle = QtWidgets.QCheckBox("Through a middle colour")
        self.use_middle.setChecked(bool(self.middle))
        self.preview = QtWidgets.QLabel()
        form = QtWidgets.QFormLayout(self)
        for key, text in (("low", "Low"), ("middle", "Middle"), ("high", "High")):
            b = QtWidgets.QPushButton()
            b.setFixedWidth(120)
            b.clicked.connect(lambda _=False, k=key: self._pick(k))
            self.buttons[key] = b
            form.addRow(text, b)
        form.addRow(self.use_middle)
        form.addRow(self.preview)
        box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                         | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        box.accepted.connect(self.accept)
        box.rejected.connect(self.reject)
        form.addRow(box)
        self.use_middle.toggled.connect(self._show)
        self._show()

    def colors(self) -> list[str]:
        middle = [self.middle or "#888888"] if self.use_middle.isChecked() else []
        return [self.low, *middle, self.high]

    def _pick(self, key: str) -> None:
        current = getattr(self, key) or "#888888"
        chosen = QtWidgets.QColorDialog.getColor(QtGui.QColor(current), self, "Colour")
        if chosen.isValid():
            setattr(self, key, chosen.name())
            if key == "middle":
                self.use_middle.setChecked(True)
            self._show()

    def _show(self) -> None:
        for key, b in self.buttons.items():
            value = getattr(self, key) or "#888888"
            b.setStyleSheet(f"background: {value}; border: 1px solid #888; min-height: 18px;")
        self.buttons["middle"].setEnabled(self.use_middle.isChecked())
        self.preview.setPixmap(_strip(figures.CUSTOM, self.colors(), 260, 16).pixmap(260, 16))


class RampPicker(QtWidgets.QComboBox):
    """Follows and sets the ramp every plot on the screen uses. With
    `follow=False` it only holds a choice (`name`, `colors`) of its own:
    a figure may differ from the screen."""

    sigChosen = QtCore.pyqtSignal()

    def __init__(self, parent=None, follow: bool = True):
        super().__init__(parent)
        self.follow = follow
        self.colors = palette.custom_colors()
        self.setIconSize(QtCore.QSize(72, 12))
        for name in figures.COLORMAPS:
            self.addItem(_strip(name), name, name)
        self.addItem(_strip(figures.CUSTOM, self.colors), "Custom…", figures.CUSTOM)
        self.set_choice(palette.ramp_name(), self.colors)
        self.activated.connect(self._activated)
        if follow:
            palette.notifier.changed.connect(self._sync)

    def name(self) -> str:
        return self.currentData()

    def set_choice(self, name: str, colors: list[str] | None = None) -> None:
        if colors:
            self.colors = list(colors)
            self.setItemIcon(self.count() - 1, _strip(figures.CUSTOM, self.colors))
        block = self.blockSignals(True)
        self.setCurrentIndex(max(self.findData(name), 0))
        self.blockSignals(block)

    def _sync(self) -> None:
        self.set_choice(palette.ramp_name(), palette.custom_colors())

    def choose_custom(self, colors: list[str]) -> None:
        """What the dialog does when it is accepted."""
        self.set_choice(figures.CUSTOM, colors)
        self._publish()

    def _activated(self, _index: int) -> None:
        if self.name() == figures.CUSTOM:
            dialog = CustomRampDialog(self.colors, self)
            if dialog.exec():
                self.choose_custom(dialog.colors())
                return
            if self.follow:
                self._sync()
                return
        self._publish()

    def _publish(self) -> None:
        if self.follow:
            palette.set_ramp(self.name(), self.colors)
        self.sigChosen.emit()
