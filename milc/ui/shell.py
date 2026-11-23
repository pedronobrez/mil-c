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
        self.stale = QtWidgets.QLabel("")
        self.stale.setProperty("role", "warning")
        self.stale.setWordWrap(True)
        self.stale.setContentsMargins(16, 6, 16, 0)
        self.stale.hide()
        centre = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(centre)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self.stale)
        lay.addWidget(self.tabs, 1)
        self.setCentralWidget(centre)
        self.statusBar().setSizeGripEnabled(False)

        self.samples.open_button.clicked.connect(self.choose_folder)
        self.samples.wsp_button.clicked.connect(self.choose_workspace)
        self.samples.sigStatus.connect(self.status)
        self.gates.sigStatus.connect(self.status)
        self.gates.sigGatesChanged.connect(self.samples.refresh)
        self.gates.sigScalesChanged.connect(self.explore._draw)
        self.gates.sigGatesChanged.connect(self.statistics.refresh)
        self.samples.sigGatesChanged.connect(self.statistics.refresh)
        self.samples.sigColumnsChanged.connect(self.statistics.mark_stale)
        self.statistics.sigStatus.connect(self.status)
        self.samples.sigGatesChanged.connect(self.gates.refresh)
        self.samples.sigNormalisationChanged.connect(self.gates.refresh)
        self.samples.sigNormalisationChanged.connect(self.statistics.mark_stale)
        self.explore.sigStatus.connect(self.status)
        self.explore.sigResults.connect(self.clusters.refresh)
        self.clusters.sigStatus.connect(self.status)
        self.clusters.sigClustersChanged.connect(self.explore._draw)
        for tab in (self.samples, self.gates, self.statistics):
            self.clusters.sigGatesChanged.connect(tab.refresh)
        self.samples.sigConcatenated.connect(self._concatenated)
        self.samples.sigOrderChanged.connect(self.order_changed)
        self.clusters.sigOrderAsked.connect(self.samples.ask_order)
        for changed in (self.gates.sigGatesChanged, self.samples.sigNormalisationChanged,
                        self.samples.sigConcatenated):
            changed.connect(self.show_stale)

        menu = self.menuBar().addMenu("&File")
        self._action(menu, "Open folder…", self.choose_folder, QtGui.QKeySequence.StandardKey.Open)
        self._action(menu, "Open workspace…", self.choose_workspace)
        menu.addSeparator()
        self._action(menu, "Open concatenated file…", self.choose_dataset)
        self._action(menu, "Save concatenated file…", self.save_dataset,
                     QtGui.QKeySequence.StandardKey.Save)
        self._action(menu, "Export tables…", self.export_tables)
        self._action(menu, "Export gates…", self.choose_gates_export)

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

        view.addSeparator()
        self._action(view, "Plate…", self.show_plate)
        self._plate = None

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

    # --- the plate ---------------------------------------------------------------

    def plate_window(self):
        """The page of several panels: one for the session, made when first
        asked for."""
        if self._plate is None:
            from .plate_window import PlateWindow
            self._plate = PlateWindow(self.session, self)
        return self._plate

    def add_to_plate(self, content, style, source=None, name: str = "", to=None) -> str:
        """`to`: the number of a plate, or "new"; None is the current one."""
        window = self.plate_window()
        letter = window.add(content, style, source, name, to)
        self.status(f"Panel {letter} added to “{window.plate.name}” (View ▸ Plate…)")
        return letter

    def show_plate(self) -> None:
        window = self.plate_window()
        window.show()
        window.raise_()

    # --- files -----------------------------------------------------------------

    def closeEvent(self, ev) -> None:
        if self._plate is not None:
            self._plate.flush()
        super().closeEvent(ev)

    def open_folder(self, folder: str) -> None:
        if self._plate is not None:
            self._plate.flush()                        # into the folder they belong to
        self.session.open_folder(folder)
        if self._plate is not None:
            self._plate.project_changed()
        settings().setValue(SETTING_FOLDER, folder)
        self.samples.refresh()
        self.gates.refresh()
        self.statistics.result = None
        self.statistics.refresh()
        self.explore.dataset_changed()
        self.clusters.refresh()
        self.show_stale()
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
            report = self.session.workspace_report()
            if report:
                self.status(f"{os.path.basename(path)}: {report}")

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

    # --- the gates, for another program ------------------------------------------------

    def export_gates(self, path: str, fmt: str | None = None, ticked_only: bool = False,
                     counts: bool = True):
        """Write the gates of every sample that has any — or of those
        ticked in Samples — as a FlowJo workspace or Gating-ML, by the
        file's name where `fmt` is not given. Returns the report."""
        from ..gate_export import export_gates
        entries = [e for e in self.session.usable()
                   if e.gated and (e.include or not ticked_only)]
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            report = export_gates(self.session, path, fmt, entries, counts)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.status(report.summary() if report.files else "No gates were written")
        return report

    def choose_gates_export(self) -> None:
        from .gates_export_dialog import GatesExportDialog, GatesReportDialog
        gated = [e for e in self.session.usable() if e.gated]
        if not gated:
            self.status("No sample has gates to export")
            return
        dialog = GatesExportDialog(self.session.folder, len(gated),
                                   sum(e.include for e in gated), self)
        if not dialog.exec() or not dialog.path.text():
            return
        try:
            report = self.export_gates(dialog.path.text(), dialog.format.currentData(),
                                       dialog.ticked.isChecked(), dialog.counts.isChecked())
        except OSError as exc:
            self.status(f"{os.path.basename(dialog.path.text())}: {exc}")
            return
        GatesReportDialog(report, self).exec()

    def order_changed(self) -> None:
        """The groups come in another order: everything that shows them is
        drawn again. Nothing is computed again — the order is of what is
        shown, not of what was measured."""
        self.gates._draw()
        self.statistics.refresh()
        if self.session.dataset is not None:
            self.explore._draw()
            self.clusters.refresh()
        self.status("The order of the groups changed")

    def show_stale(self) -> None:
        """Above the tabs, for as long as the concatenated events are older
        than the values they were taken from."""
        why = self.session.dataset_stale
        self.stale.setText(
            f"The concatenated events were taken before {why}: Explore and Clusters still "
            "show the old values. Concatenate again in Samples to bring them up to date."
            if why else "")
        self.stale.setVisible(bool(why))

    def _concatenated(self) -> None:
        self.explore.dataset_changed()
        self.clusters.refresh()
        self.tabs.setCurrentWidget(self.explore)
