"""A row of controls that goes on to the next line when the window is
narrower than the row, and a label that shortens itself: what lets the
window fit a laptop's screen."""

from __future__ import annotations

from PyQt6 import QtCore, QtGui, QtWidgets


class FlowLayout(QtWidgets.QLayout):
    """Widgets left to right, wrapping; each line as tall as its tallest,
    and everything on it centred on the line."""

    #: the row does not ask to be narrower than this, if it is wider on one
    #: line: a layout sizes a window for the height a row has at its
    #: narrowest, and a row a button wide is as tall as the screen
    FLOOR = 640

    def __init__(self, parent=None, spacing: int = 8, line_spacing: int = 6):
        super().__init__(parent)
        self._items: list[QtWidgets.QLayoutItem] = []
        self._gap, self._line_gap = spacing, line_spacing
        self.setContentsMargins(0, 0, 0, 0)

    def addItem(self, item) -> None:
        self._items.append(item)

    def count(self) -> int:
        return len(self._items)

    def itemAt(self, index: int):
        return self._items[index] if 0 <= index < len(self._items) else None

    def takeAt(self, index: int):
        return self._items.pop(index) if 0 <= index < len(self._items) else None

    def expandingDirections(self):
        return QtCore.Qt.Orientation(0)

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, width: int) -> int:
        return self._place(QtCore.QRect(0, 0, width, 0), move=False)

    def setGeometry(self, rect) -> None:
        super().setGeometry(rect)
        self._place(rect, move=True)

    def sizeHint(self):
        return self.minimumSize()

    def minimumSize(self):
        """As wide as the widest single item, or as `FLOOR` if all of it
        on one line is wider than that: the rest can wrap."""
        size = QtCore.QSize()
        shown = self._shown()
        for item in shown:
            size = size.expandedTo(item.minimumSize())
        line = sum(i.sizeHint().width() for i in shown) + self._gap * max(len(shown) - 1, 0)
        size.setWidth(max(size.width(), min(line, self.FLOOR)))
        m = self.contentsMargins()
        return size + QtCore.QSize(m.left() + m.right(), m.top() + m.bottom())

    def _shown(self):
        return [i for i in self._items if not i.isEmpty()]

    def _place(self, rect, move: bool) -> int:
        m = self.contentsMargins()
        area = rect.adjusted(m.left(), m.top(), -m.right(), -m.bottom())
        lines: list[list] = [[]]
        x = 0
        for item in self._shown():
            w = item.sizeHint().width()
            if lines[-1] and x + w > area.width():
                lines.append([])
                x = 0
            lines[-1].append(item)
            x += w + self._gap
        y = area.y()
        for line in lines:
            if not line:
                continue
            height = max(i.sizeHint().height() for i in line)
            x = area.x()
            for item in line:
                hint = item.sizeHint()
                if move:
                    item.setGeometry(QtCore.QRect(
                        x, y + (height - hint.height()) // 2,
                        min(hint.width(), area.width()), hint.height()))
                x += hint.width() + self._gap
            y += height + self._line_gap
        return y - self._line_gap - rect.y() + m.bottom()


def together(*widgets, spacing: int = 6) -> QtWidgets.QWidget:
    """Widgets that stay on one line of a flow: a caption and its box."""
    box = QtWidgets.QWidget()
    lay = QtWidgets.QHBoxLayout(box)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(spacing)
    for w in widgets:
        lay.addWidget(w)
    return box


def captioned(text: str, *widgets) -> QtWidgets.QWidget:
    cap = QtWidgets.QLabel(text)
    cap.setProperty("role", "caption")
    return together(cap, *widgets)


def flow(*widgets, spacing: int = 14) -> QtWidgets.QWidget:
    """A widget holding `widgets` in a flow: its height follows its width."""
    box = QtWidgets.QWidget()
    lay = FlowLayout(box, spacing=spacing)
    for w in widgets:
        lay.addWidget(w)
    policy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred,
                                   QtWidgets.QSizePolicy.Policy.Minimum)
    policy.setHeightForWidth(True)
    box.setSizePolicy(policy)
    return box


class ElidedLabel(QtWidgets.QLabel):
    """One line of text that gives way: where it does not fit, its middle
    becomes an ellipsis and the whole of it is in the tooltip."""

    def minimumSizeHint(self):
        return QtCore.QSize(60, super().minimumSizeHint().height())

    def paintEvent(self, ev) -> None:
        painter = QtGui.QPainter(self)
        text = self.fontMetrics().elidedText(
            self.text(), QtCore.Qt.TextElideMode.ElideMiddle, self.width())
        self.setToolTip(self.text() if text != self.text() else "")
        painter.drawText(self.rect(), int(self.alignment()), text)
