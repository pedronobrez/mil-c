"""
What can be done to a cluster, as a row of buttons under the Clusters tab:
name it, annotate it, merge it with others, have names proposed from what
the clusters express, turn it into a population of the files, keep a
lasso as a cluster of its own.

Every button asks a question and then calls a plain method that does the
work — `rename`, `merge`, `make_population`… — so the work can be done
without the question. The edits themselves are `milc.clusters`'; nothing
here knows how they are kept.
"""

from __future__ import annotations

from PyQt6 import QtCore, QtWidgets

from .. import clusters
from . import palette
from .scatter import Categories

ONLY_HERE = ("It exists in this dataset only. It cannot become a population of the full "
             "samples: the map it was outlined on has coordinates for the events drawn "
             "here and for no others.")


def cluster_categories(now: clusters.Clustering) -> Categories:
    """The clusters as the user has them, for a plot to colour by."""
    return Categories(now.codes(), list(now.names), palette.categorical(len(now.ids)),
                      f"{now.key} cluster")


class ClusterTools(QtWidgets.QWidget):
    def __init__(self, tab):
        super().__init__(tab)
        self.tab = tab
        row = QtWidgets.QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        self.buttons: dict[str, QtWidgets.QPushButton] = {}
        for text, slot, tip in (
                ("Rename…", self._ask_rename,
                 "Call the selected cluster something: the name is shown wherever "
                 "its number was"),
                ("Note…", self._ask_note, "A note kept with the selected cluster"),
                ("Merge with…", self._ask_merge,
                 "Make one cluster of the selected one and others. The algorithm's "
                 "clusters are kept underneath: Undo merges brings them back"),
                ("Undo merges", self.undo_merges,
                 "Back to the clusters the algorithm gave"),
                ("Suggest names…", self._ask_names,
                 "Propose a name for each cluster from the markers it is positive "
                 "for — CD4+ Foxp3+ — to accept, edit or refuse"),
                ("Make population", self._make_population,
                 "The selected cluster as a population in Gates and Statistics, on "
                 "every event of every file: each event goes to its nearest node of "
                 "the FlowSOM map"),
                ("Save selection…", self._ask_selection,
                 "Keep the events outlined on the map (Shift+drag) as a cluster of "
                 "their own; the other events stay in their clusters. " + ONLY_HERE),
                ("Remove", self.remove,
                 "Remove the selected hand-drawn cluster: its events go back to "
                 "the clusters they were in")):
            button = QtWidgets.QPushButton(text)
            button.setToolTip(tip)
            button.clicked.connect(lambda _=False, s=slot: self._try(s))
            row.addWidget(button)
            self.buttons[text.rstrip("…")] = button
        row.addStretch(1)

    @property
    def ds(self):
        return self.tab.session.dataset

    def now(self) -> clusters.Clustering | None:
        return self.tab.session.clusters(self.tab.key)

    def sync(self) -> None:
        """Each button is live when there is something for it to act on."""
        now, chosen = self.now(), self.tab.selected
        manual = chosen is not None and now is not None and now.is_manual(chosen)
        algorithm = chosen is not None and not manual
        live = {"Rename": chosen is not None, "Note": chosen is not None,
                "Merge with": algorithm and sum(not now.is_manual(c) for c in now.ids) > 1,
                "Undo merges": now is not None and now.merged(),
                "Suggest names": now is not None,
                "Make population": algorithm,
                "Save selection": self.tab._lasso is not None,
                "Remove": manual}
        for name, on in live.items():
            self.buttons[name].setEnabled(bool(on))

    def _try(self, action) -> None:
        try:
            action()
        except ValueError as exc:
            self.tab.sigStatus.emit(str(exc))

    # --- what the buttons do ---------------------------------------------------

    def rename(self, name: str, cluster: int | None = None) -> None:
        cluster = self.tab.selected if cluster is None else cluster
        clusters.rename(self.ds, {cluster: name}, self.tab.key)
        self.tab.clusters_changed(cluster)

    def set_note(self, note: str, cluster: int | None = None) -> None:
        cluster = self.tab.selected if cluster is None else cluster
        clusters.set_note(self.ds, cluster, note, self.tab.key)
        self.tab.clusters_changed(cluster)

    def merge(self, others: list[int]) -> int:
        """The selected cluster and `others` become one, which stays selected."""
        merged = clusters.merge(self.ds, [self.tab.selected, *others], self.tab.key)
        self.tab.clusters_changed(merged)
        return merged

    def undo_merges(self) -> None:
        clusters.undo_merges(self.ds, self.tab.key)
        self.tab.clusters_changed()

    def suggestions(self) -> dict[int, str]:
        return self.tab.session.suggest_cluster_names(self.tab.key)

    def apply_names(self, names: dict[int, str]) -> None:
        clusters.rename(self.ds, names, self.tab.key)
        self.tab.clusters_changed(self.tab.selected)

    def make_population(self, name: str | None = None) -> str:
        path = self.tab.session.make_population(self.tab.selected, self.tab.key, name)
        self.tab.sigGatesChanged.emit()
        return path

    def save_selection(self, name: str) -> int:
        """The lasso as a cluster, which becomes the selected one."""
        if self.tab._lasso is None:
            raise ValueError("outline some events on the map first: Shift+drag")
        cluster = clusters.add_manual(self.ds, self.tab._lasso, name, self.tab.key)
        self.tab.clusters_changed(cluster)
        return cluster

    def remove(self) -> None:
        clusters.remove_manual(self.ds, self.tab.selected, self.tab.key)
        self.tab.clusters_changed()

    # --- the questions -----------------------------------------------------------

    def _ask_rename(self) -> None:
        now = self.now()
        name, ok = QtWidgets.QInputDialog.getText(
            self, "Rename cluster", "Name (empty: back to its number)",
            text=now.name(self.tab.selected))
        if ok:
            self.rename(name)

    def _ask_note(self) -> None:
        now = self.now()
        note, ok = QtWidgets.QInputDialog.getMultiLineText(
            self, "Note", f"A note on cluster {now.name(self.tab.selected)}",
            now.note(self.tab.selected))
        if ok:
            self.set_note(note)

    def _ask_merge(self) -> None:
        now = self.now()
        dialog = MergeDialog(now, self.tab.selected, self)
        if dialog.exec() and dialog.chosen():
            self.merge(dialog.chosen())

    def _ask_names(self) -> None:
        found = self.suggestions()
        if not found:
            self.tab.sigStatus.emit("No marker separates the clusters: there is no name "
                                    "to suggest")
            return
        dialog = NamesDialog(self.now(), found, self)
        if dialog.exec():
            self.apply_names(dialog.names())

    def _make_population(self) -> None:
        path = self.make_population()
        self.tab.sigStatus.emit(
            f"{path} is now a population of every file: see Gates and Statistics")

    def _ask_selection(self) -> None:
        name, ok = QtWidgets.QInputDialog.getText(
            self, "Save selection as a cluster",
            "A name for the outlined events.\n" + ONLY_HERE.replace(". ", ".\n", 1))
        if ok:
            self.save_selection(name)


class MergeDialog(QtWidgets.QDialog):
    """Which clusters to merge the selected one with."""

    def __init__(self, now: clusters.Clustering, selected: int, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Merge clusters")
        head = QtWidgets.QLabel(f"Merge cluster {now.name(selected)} with:")
        self.list = QtWidgets.QListWidget()
        for cluster, name in zip(now.ids, now.names):
            if cluster == selected or now.is_manual(cluster):
                continue
            item = QtWidgets.QListWidgetItem(name)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, cluster)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.list.addItem(item)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(head)
        layout.addWidget(self.list, 1)
        layout.addWidget(_ok_cancel(self, "Merge"))

    def chosen(self) -> list[int]:
        items = [self.list.item(i) for i in range(self.list.count())]
        return [it.data(QtCore.Qt.ItemDataRole.UserRole) for it in items
                if it.checkState() == QtCore.Qt.CheckState.Checked]


class NamesDialog(QtWidgets.QDialog):
    """The names proposed, beside the ones in use, to edit before any is
    applied. A row left empty is a cluster left as it is."""

    def __init__(self, now: clusters.Clustering, found: dict[int, str], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Suggested names")
        head = QtWidgets.QLabel(
            "For each marker the clusters are cut in two at the widest gap between "
            "their values: those above are positive, those below negative. A marker "
            "whose widest gap is narrow separates nothing and is left out. Type over "
            "any name; clear one to keep the cluster as it is.")
        head.setWordWrap(True)
        head.setProperty("role", "caption")
        self._ids = list(found)
        self.table = QtWidgets.QTableWidget(len(found), 2)
        self.table.setHorizontalHeaderLabels(["Now", "Suggested"])
        self.table.verticalHeader().setVisible(False)
        for r, (cluster, name) in enumerate(found.items()):
            current = QtWidgets.QTableWidgetItem(now.name(cluster))
            current.setFlags(current.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(r, 0, current)
            self.table.setItem(r, 1, QtWidgets.QTableWidgetItem(name))
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.resizeColumnToContents(0)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(head)
        layout.addWidget(self.table, 1)
        layout.addWidget(_ok_cancel(self, "Apply"))
        self.resize(460, 420)

    def names(self) -> dict[int, str]:
        typed = {c: self.table.item(r, 1).text().strip() for r, c in enumerate(self._ids)}
        return {c: n for c, n in typed.items() if n}


def _ok_cancel(dialog: QtWidgets.QDialog, accept: str) -> QtWidgets.QDialogButtonBox:
    box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                     | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
    box.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setText(accept)
    box.accepted.connect(dialog.accept)
    box.rejected.connect(dialog.reject)
    return box
