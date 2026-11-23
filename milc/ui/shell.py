"""
The application shell: one window holding the tabs, which all read the
same session. The window owns the files, the menus and the status bar.
"""

from __future__ import annotations

import os

from PyQt6 import QtCore, QtGui, QtWidgets

from ..dataset import Dataset
from ..session import Session
from . import style, theme
from .clusters_tab import ClustersTab
from .explore_tab import ExploreTab
from .gates_tab import GatesTab
from .samples_tab import SamplesTab
from .settings import settings
from .statistics_tab import StatisticsTab

SETTING_FOLDER = "files/last_folder"
SETTING_APPEARANCE = "view/appearance"


class MainShell(QtWidgets.QMainWindow):
    WANTED_SIZE = (1500, 920)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("MIL-C")
        self.session = Session()
        self.samples = SamplesTab(self.session)
        self.explore = ExploreTab(self.session)
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setDocumentMode(True)
        self.gates = GatesTab(self.session)
        self.tabs.addTab(self.samples, "Samples")
        self.tabs.addTab(self.gates, "Gates")
        self.statistics = StatisticsTab(self.session)
        self.tabs.addTab(self.statistics, "Statistics")
        self.clusters = ClustersTab(self.session)
        self.tabs.addTab(self.explore, "Explore")
        self.tabs.addTab(self.clusters, "Clusters")
        self.setCentralWidget(self.tabs)
        self.statusBar().setSizeGripEnabled(False)

        self.samples.open_button.clicked.connect(self.choose_folder)
        self.samples.wsp_button.clicked.connect(self.choose_workspace)
        self.samples.sigStatus.connect(self.status)
        self.gates.sigStatus.connect(self.status)
        self.gates.sigGatesChanged.connect(self.samples.refresh)
        self.gates.sigGatesChanged.connect(self.statistics.refresh)
        self.samples.sigGatesChanged.connect(self.statistics.refresh)
        self.samples.sigColumnsChanged.connect(self.statistics.mark_stale)
        self.statistics.sigStatus.connect(self.status)
        self.samples.sigGatesChanged.connect(self.gates.refresh)
        self.explore.sigStatus.connect(self.status)
        self.explore.sigResults.connect(self.clusters.refresh)
        self.clusters.sigStatus.connect(self.status)
        self.samples.sigConcatenated.connect(self._concatenated)

        menu = self.menuBar().addMenu("&File")
        self._action(menu, "Open folder…", self.choose_folder, QtGui.QKeySequence.StandardKey.Open)
        self._action(menu, "Open workspace…", self.choose_workspace)
        menu.addSeparator()
        self._action(menu, "Open concatenated file…", self.choose_dataset)
        self._action(menu, "Save concatenated file…", self.save_dataset,
                     QtGui.QKeySequence.StandardKey.Save)
        self._action(menu, "Export tables…", self.export_tables)

        view = self.menuBar().addMenu("&View")
        looks = view.addMenu("Appearance")
        group = QtGui.QActionGroup(self)
        chosen = settings().value(SETTING_APPEARANCE, "") or ""
        self._appearance = {}
        for key, text in (("", "As the system"), ("light", "Light"), ("dark", "Dark")):
            action = looks.addAction(text)
            action.setCheckable(True)
            action.setChecked(key == chosen)
            action.triggered.connect(lambda _=False, k=key: self.set_appearance(k))
            group.addAction(action)
            self._appearance[key] = action

        screen = QtWidgets.QApplication.primaryScreen()
        room = screen.availableGeometry() if screen else QtCore.QRect(0, 0, 1280, 800)
        self.resize(min(self.WANTED_SIZE[0], room.width() - 40),
                    min(self.WANTED_SIZE[1], room.height() - 60))

        app = QtWidgets.QApplication.instance()
        self._watcher = style.ThemeWatcher(app, self)
        self._watcher.sigThemeChanged.connect(self._retheme)
        if chosen:
            self.set_appearance(chosen)

    def _action(self, menu, text, slot, shortcut=None):
        action = menu.addAction(text)
        action.triggered.connect(lambda _=False: slot())
        if shortcut is not None:
            action.setShortcut(shortcut)
        return action

    def status(self, text: str) -> None:
        self.statusBar().showMessage(text)

    def set_appearance(self, scheme: str) -> None:
        """Light, dark, or "" for whatever the system is in."""
        style.set_appearance(scheme or None)
        settings().setValue(SETTING_APPEARANCE, scheme)
        self._appearance[scheme if scheme in self._appearance else ""].setChecked(True)
        style.apply(QtWidgets.QApplication.instance())
        if hasattr(self, "_watcher"):
            self._retheme()

    def _retheme(self) -> None:
        theme.apply_defaults()
        self.explore.view.retheme()
        self.clusters.map.setBackground(theme.background())
        self.gates._draw()
        self.clusters.refresh()
        self.samples.refresh()

    # --- files -----------------------------------------------------------------

    def open_folder(self, folder: str) -> None:
        self.session.open_folder(folder)
        settings().setValue(SETTING_FOLDER, folder)
        self.samples.refresh()
        self.gates.refresh()
        self.statistics.result = None
        self.statistics.refresh()
        self.explore.dataset_changed()
        self.clusters.refresh()
        self.tabs.setCurrentWidget(self.samples)
        n = len(self.session.entries)
        self.status(f"{n} file(s) in {folder}" if n else f"No .fcs files in {folder}")

    def choose_folder(self) -> None:
        start = settings().value(SETTING_FOLDER, "") or os.path.expanduser("~")
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self, "A folder of .fcs files", start)
        if folder:
            self.open_folder(folder)

    def choose_workspace(self) -> None:
        if not self.session.folder:
            self.status("Open a folder first")
            return
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "FlowJo workspace", self.session.folder, "FlowJo workspace (*.wsp)")
        if path:
            try:
                self.session.set_workspace(path)
            except Exception as exc:
                self.status(f"{os.path.basename(path)}: {exc}")
                return
            self.samples.refresh()
            self.gates.refresh()
            self.statistics.refresh()

    def open_dataset(self, path: str) -> None:
        self.session.dataset = (Dataset.load(path) if path.endswith(".npz")
                                else Dataset.from_fcs(path))
        self.session.flowsom = None
        self._concatenated()
        self.status(f"{len(self.session.dataset):,d} events from {os.path.basename(path)}")

    def choose_dataset(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Concatenated file", self.session.folder,
            "Concatenated (*.npz *.fcs)")
        if path:
            try:
                self.open_dataset(path)
            except Exception as exc:
                self.status(f"{os.path.basename(path)}: {exc}")

    def save_dataset(self) -> None:
        ds = self.session.dataset
        if ds is None:
            self.status("Nothing concatenated yet")
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save the concatenated events",
            os.path.join(self.session.folder, "concat.npz"),
            "MIL-C dataset, with results (*.npz);;FCS, for other programs (*.fcs)")
        if not path:
            return
        ds.to_fcs(path) if path.endswith(".fcs") else ds.save(path)
        self.status(f"Wrote {path}")

    def export_tables(self) -> None:
        tables = self.session.cluster_tables()
        if tables is None:
            self.status("Run FlowSOM first: the tables are per cluster")
            return
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Where to write the tables", self.session.folder)
        if not folder:
            return
        expression, freq = tables
        expression.to_csv(os.path.join(folder, "cluster_expression.csv"), "cluster")
        freq.to_csv(os.path.join(folder, "cluster_frequency.csv"), "cluster")
        self.status(f"Wrote cluster_expression.csv and cluster_frequency.csv to {folder}")

    def _concatenated(self) -> None:
        self.explore.dataset_changed()
        self.clusters.refresh()
        self.tabs.setCurrentWidget(self.explore)
