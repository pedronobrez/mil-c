"""Which group comes first: the order of a column's values on every chart,
legend and table."""

from __future__ import annotations

from PyQt6 import QtCore, QtWidgets
from .help_window import describe


class OrderDialog(QtWidgets.QDialog):
    """A column, and its values to drag — or move with the buttons — into
    the order they should be shown in. Each column keeps its own; what is
    arranged for one stays while another is looked at, and all of it is
    applied on OK."""

    def __init__(self, session, parent=None, column: str = ""):
        super().__init__(parent)
        describe(self, "sample-order")
        self.setWindowTitle("Order of the groups")
        self.session = session
        #: column → its values as arranged here, or None for counting order
        self.arranged: dict[str, list[str] | None] = {}
        self._column = ""

        head = QtWidgets.QLabel(
            "The order of a column's values wherever they are shown: bars, legends, panels, "
            "tables, the Prism file. Drag a value, or select it and use the buttons.")
        head.setWordWrap(True)
        head.setProperty("role", "caption")
        self.column = QtWidgets.QComboBox()
        for c in session.columns:
            self.column.addItem(c, c)
        self.list = QtWidgets.QListWidget()
        self.list.setDragDropMode(QtWidgets.QAbstractItemView.DragDropMode.InternalMove)
        self.list.setDefaultDropAction(QtCore.Qt.DropAction.MoveAction)
        self.up = QtWidgets.QPushButton("Move up")
        self.down = QtWidgets.QPushButton("Move down")
        self.reverse = QtWidgets.QPushButton("Reverse")
        self.counting = QtWidgets.QPushButton("Counting order")
        self.counting.setToolTip("Back to the order nobody chose: 3, 7, 14; A, B, C")
        self.state = QtWidgets.QLabel("")
        self.state.setProperty("role", "caption")
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        side = QtWidgets.QVBoxLayout()
        for b in (self.up, self.down, self.reverse, self.counting):
            side.addWidget(b)
        side.addStretch(1)
        middle = QtWidgets.QHBoxLayout()
        middle.addWidget(self.list, 1)
        middle.addLayout(side)
        top = QtWidgets.QHBoxLayout()
        cap = QtWidgets.QLabel("Column")
        cap.setProperty("role", "caption")
        top.addWidget(cap)
        top.addWidget(self.column, 1)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(head)
        layout.addLayout(top)
        layout.addLayout(middle, 1)
        layout.addWidget(self.state)
        layout.addWidget(buttons)
        self.resize(420, 400)

        self.up.clicked.connect(lambda: self.move(-1))
        self.down.clicked.connect(lambda: self.move(1))
        self.reverse.clicked.connect(self.reverse_all)
        self.counting.clicked.connect(self.counting_order)
        self.column.currentIndexChanged.connect(self._column_chosen)
        self.list.model().rowsMoved.connect(self._moved)
        self.list.currentRowChanged.connect(self._sync)
        self.set_column(column if column in session.columns else self.column.currentData() or "")

    # --- what is shown ---------------------------------------------------------

    def values(self) -> list[str]:
        return [self.list.item(i).text() for i in range(self.list.count())]

    def set_column(self, column: str) -> None:
        self._keep()
        self._column = column
        self.column.blockSignals(True)
        self.column.setCurrentIndex(max(self.column.findData(column), 0))
        self.column.blockSignals(False)
        self._fill(self.arranged[column] if self.arranged.get(column)
                   else self.session.values_of(column))

    def _column_chosen(self, *_) -> None:
        self.set_column(self.column.currentData() or "")

    def _fill(self, values: list[str]) -> None:
        self.list.clear()
        self.list.addItems(values)
        if values:
            self.list.setCurrentRow(0)
        self._sync()

    def _keep(self) -> None:
        """What is arranged for the column shown, before another is."""
        if self._column and self.list.count():
            self.arranged[self._column] = self.values()

    def _given(self, column: str) -> bool:
        """Whether this column will have an order of its own on OK."""
        if column in self.arranged:
            return self.arranged[column] is not None
        return column in self.session.orders

    def _sync(self, *_) -> None:
        row, n = self.list.currentRow(), self.list.count()
        self.up.setEnabled(row > 0)
        self.down.setEnabled(0 <= row < n - 1)
        self.reverse.setEnabled(n > 1)
        self.state.setText("No values in this column yet." if n == 0 else
                           "In the order given here." if self._given(self._column) else
                           "In counting order: nobody chose another.")

    # --- changing it -----------------------------------------------------------

    def _moved(self, *_) -> None:
        self.arranged[self._column] = self.values()
        self._sync()

    def move(self, by: int) -> None:
        row = self.list.currentRow()
        to = row + by
        if row < 0 or not 0 <= to < self.list.count():
            return
        self.list.insertItem(to, self.list.takeItem(row))
        self.list.setCurrentRow(to)
        self._moved()

    def set_values(self, values: list[str]) -> None:
        """The column shown, in this order; values not named keep their
        place after those that are."""
        rest = [v for v in self.values() if v not in values]
        self._fill([v for v in values if v in self.values()] + rest)
        self._moved()

    def reverse_all(self) -> None:
        self.set_values(self.values()[::-1])

    def counting_order(self) -> None:
        from .. import order
        self.arranged[self._column] = None
        # as it would be with nothing given for this column
        kept = dict(self.session.orders)
        try:
            order.set_order({c: v for c, v in kept.items() if c != self._column})
            values = order.ordered(self.values())
        finally:
            order.set_order(kept)
        self._fill(values)

    def apply(self) -> bool:
        """Give the session what was arranged; whether anything changed."""
        if self.arranged.get(self._column) is not None:
            self._keep()
        changed = False
        for column, values in self.arranged.items():
            if (values or None) != self.session.orders.get(column):
                self.session.set_order(column, values)
                changed = True
        return changed


def ask(session, parent=None, column: str = "") -> bool:
    """Show the dialog; whether the order changed."""
    dialog = OrderDialog(session, parent, column)
    return bool(dialog.exec()) and dialog.apply()
