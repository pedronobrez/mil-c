"""
The second tab: run the algorithms on the concatenated events and look at
the result — any two axes, coloured by group, sample, cluster or the
expression of a marker, in one panel or one per group.
"""

from __future__ import annotations

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets

from .. import figures
from ..session import Session
from ..transforms import Biex, decade_ticks
from . import palette, style, worker
from .ramp_picker import RampPicker
from .scatter import Categories, MapView, Values, categories_of
from .cluster_tools import cluster_categories



def _spec(kind: str, key: str = "") -> str:
    """What a choice in a combo stands for, as the string Qt can keep and
    find again: a tuple comes back as a list and `findData` misses it."""
    return f"{kind}|{key}"


def _parse(data) -> tuple[str, str] | None:
    if not data:
        return None
    kind, _, key = str(data).partition("|")
    return kind, key


class ExploreTab(QtWidgets.QWidget):
    sigStatus = QtCore.pyqtSignal(str)
    #: an algorithm finished: there is something new to show elsewhere
    sigResults = QtCore.pyqtSignal()

    def __init__(self, session: Session, parent=None):
        super().__init__(parent)
        self.session = session
        self._filling = False
        self._running: set[str] = set()

        # --- left: what to compute -------------------------------------------
        self.markers = QtWidgets.QListWidget()
        self.markers.setToolTip("The markers the algorithms see. Leave out "
                                "what was used to gate: viability, lineage.")
        self.transform = QtWidgets.QComboBox()
        self.transform.addItem("Biexponential", "biex")
        self.transform.setToolTip("Biexponential is what FlowJo feeds its plugins")
        self.transform.addItem("Arcsinh, cofactor 150", "arcsinh")
        self.seed = QtWidgets.QSpinBox()
        self.seed.setRange(0, 999_999)

        self.neighbors = QtWidgets.QSpinBox()
        self.neighbors.setRange(2, 500)
        self.neighbors.setValue(15)
        self.min_dist = QtWidgets.QDoubleSpinBox()
        self.min_dist.setRange(0.0, 0.99)
        self.min_dist.setSingleStep(0.05)
        self.min_dist.setValue(0.5)
        self.perplexity = QtWidgets.QDoubleSpinBox()
        self.perplexity.setRange(2, 500)
        self.perplexity.setValue(30)
        self.clusters = QtWidgets.QSpinBox()
        self.clusters.setRange(2, 60)
        self.clusters.setValue(5)
        self.grid = QtWidgets.QSpinBox()
        self.grid.setRange(2, 30)
        self.grid.setValue(10)
        self.grid.setToolTip("The map is this many nodes on each side")
        self.k = QtWidgets.QSpinBox()
        self.k.setRange(2, 500)
        self.k.setValue(30)
        self.k.setToolTip("Each event is joined to this many nearest "
                          "neighbours: fewer finds smaller clusters")
        self.resolution = QtWidgets.QDoubleSpinBox()
        self.resolution.setRange(0.05, 10.0)
        self.resolution.setSingleStep(0.1)
        self.resolution.setValue(1.0)
        self.resolution.setToolTip("Above 1, more and smaller clusters; below, "
                                   "fewer and larger — lower it when one "
                                   "population comes out in pieces")

        self.run_buttons = {
            "UMAP": QtWidgets.QPushButton("Run UMAP"),
            "tSNE": QtWidgets.QPushButton("Run tSNE"),
            "FlowSOM": QtWidgets.QPushButton("Run FlowSOM"),
            "PhenoGraph": QtWidgets.QPushButton("Run PhenoGraph"),
            "PCA": QtWidgets.QPushButton("Run PCA"),
        }
        self.run_all = QtWidgets.QPushButton("Run all three")
        self.run_all.setProperty("primary", True)
        self.result = QtWidgets.QLabel("")
        self.result.setWordWrap(True)
        self.result.setProperty("role", "caption")

        left = QtWidgets.QVBoxLayout()
        left.addWidget(self._eyebrow("Markers"))
        left.addWidget(self.markers, 1)
        form = QtWidgets.QFormLayout()
        form.addRow("Scale", self.transform)
        form.addRow("Seed", self.seed)
        form.addRow(self._eyebrow("UMAP"))
        form.addRow("Neighbours", self.neighbors)
        form.addRow("Min. distance", self.min_dist)
        form.addRow(self.run_buttons["UMAP"])
        form.addRow(self._eyebrow("tSNE"))
        form.addRow("Perplexity", self.perplexity)
        form.addRow(self.run_buttons["tSNE"])
        form.addRow(self._eyebrow("FlowSOM"))
        form.addRow("Clusters", self.clusters)
        form.addRow("Grid", self.grid)
        form.addRow(self.run_buttons["FlowSOM"])
        form.addRow(self._eyebrow("PhenoGraph"))
        form.addRow("Neighbours", self.k)
        form.addRow("Resolution", self.resolution)
        form.addRow(self.run_buttons["PhenoGraph"])
        form.addRow(self._eyebrow("PCA"))
        form.addRow(self.run_buttons["PCA"])
        left.addLayout(form)
        # the settings are taller than a laptop's screen now: they scroll,
        # sooner than squeeze the list of markers to nothing
        self.markers.setMinimumHeight(120)
        left.setContentsMargins(0, 0, 0, 0)
        settings = QtWidgets.QWidget()
        settings.setLayout(left)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(settings)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        # room for the bar beside the fields, not over them
        scroll.setMinimumWidth(settings.minimumSizeHint().width()
                               + scroll.verticalScrollBar().sizeHint().width())
        left = QtWidgets.QVBoxLayout()
        left.addWidget(scroll, 1)
        left.addSpacing(6)
        left.addWidget(self.run_all)
        left.addWidget(self.result)
        left_panel = QtWidgets.QWidget()
        left_panel.setLayout(left)
        left_panel.setFixedWidth(max(280, left_panel.minimumSizeHint().width()))

        # --- middle: the plot --------------------------------------------------
        self.view = MapView()
        self.empty = QtWidgets.QLabel(
            "Concatenate the samples in the first tab, then run an algorithm here.")
        self.empty.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.empty.setProperty("role", "hint")
        self.stack = QtWidgets.QStackedWidget()
        self.stack.addWidget(self.empty)
        self.stack.addWidget(self.view)

        # --- right: how to look at it --------------------------------------------
        self.x_axis = QtWidgets.QComboBox()
        self.y_axis = QtWidgets.QComboBox()
        self.color_by = QtWidgets.QComboBox()
        self.split_by = QtWidgets.QComboBox()
        self.split_by2 = QtWidgets.QComboBox()
        self.split_by2.setToolTip("A second subdivision: one panel for every "
                                  "combination, Case · 3, Case · 7…")
        self.size = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.size.setRange(1, 14)
        self.size.setValue(4)
        self.alpha = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.alpha.setRange(5, 100)
        self.alpha.setValue(80)
        self.legend = QtWidgets.QListWidget()
        self.legend.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        self.legend.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        self.legend.setWordWrap(True)
        self.legend.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        right_form = QtWidgets.QFormLayout()
        right_form.addRow(self._eyebrow("Axes"))
        right_form.addRow("X", self.x_axis)
        right_form.addRow("Y", self.y_axis)
        right_form.addRow(self._eyebrow("Colour"))
        right_form.addRow("Colour by", self.color_by)
        right_form.addRow("Split by", self.split_by)
        right_form.addRow("and by", self.split_by2)
        self.ramp = RampPicker()
        right_form.addRow("Ramp", self.ramp)
        right_form.addRow("Point size", self.size)
        right_form.addRow("Opacity", self.alpha)
        right = QtWidgets.QVBoxLayout()
        right.addLayout(right_form)
        right.addWidget(self._eyebrow("Legend"))
        right.addWidget(self.legend, 1)
        self.figure_button = QtWidgets.QPushButton("Figure…")
        self.figure_button.setToolTip("This plot as a figure: size, type, colours, "
                                      "legend — SVG, PDF, PNG or TIFF")
        right.addWidget(self.figure_button)
        right_panel = QtWidgets.QWidget()
        right_panel.setLayout(right)
        right_panel.setFixedWidth(250)

        layout = QtWidgets.QHBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.addWidget(left_panel)
        layout.addWidget(self.stack, 1)
        layout.addWidget(right_panel)

        for name, button in self.run_buttons.items():
            button.clicked.connect(lambda _=False, n=name: self.run([n]))
        self.run_all.clicked.connect(lambda: self.run(["UMAP", "tSNE", "FlowSOM"]))
        for combo in (self.x_axis, self.y_axis, self.color_by, self.split_by, self.split_by2):
            combo.currentIndexChanged.connect(self._draw)
        self.size.valueChanged.connect(lambda v: self.view.set_point_size(v))
        self.figure_button.clicked.connect(self.make_figure)
        palette.notifier.changed.connect(self._draw)
        self.alpha.valueChanged.connect(lambda v: self.view.set_opacity(v / 100))
        self.dataset_changed()

    @staticmethod
    def _eyebrow(text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setFont(style.eyebrow_font())
        label.setProperty("role", "caption")
        label.setContentsMargins(0, 10, 0, 2)
        return label

    # --- the dataset ---------------------------------------------------------

    def dataset_changed(self) -> None:
        """A new concatenation: the markers, and nothing computed yet."""
        ds = self.session.dataset
        self._filling = True
        self.markers.clear()
        if ds is not None:
            chosen = set(self.session.default_channels())
            for c, m, comp in zip(ds.channels, ds.markers, ds.compensated):
                if not comp:
                    continue  # scatter and time are not markers
                item = QtWidgets.QListWidgetItem(m or c)
                item.setData(QtCore.Qt.ItemDataRole.UserRole, c)
                item.setToolTip(c)
                item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(QtCore.Qt.CheckState.Checked if c in chosen
                                   else QtCore.Qt.CheckState.Unchecked)
                self.markers.addItem(item)
        self._filling = False
        for w in (*self.run_buttons.values(), self.run_all):
            w.setEnabled(ds is not None)
        self.results_changed(prefer=None)

    def chosen_markers(self) -> list[str]:
        return [self.markers.item(i).data(QtCore.Qt.ItemDataRole.UserRole)
                for i in range(self.markers.count())
                if self.markers.item(i).checkState() == QtCore.Qt.CheckState.Checked]

    def set_markers(self, channels: list[str]) -> None:
        for i in range(self.markers.count()):
            item = self.markers.item(i)
            item.setCheckState(
                QtCore.Qt.CheckState.Checked
                if item.data(QtCore.Qt.ItemDataRole.UserRole) in channels
                else QtCore.Qt.CheckState.Unchecked)

    def results_changed(self, prefer: str | None) -> None:
        """Refill the choices, keeping what was chosen where it still
        exists; `prefer` is the map that has just been computed."""
        ds = self.session.dataset
        keep = {c: c.currentData() for c in
                (self.x_axis, self.y_axis, self.color_by, self.split_by, self.split_by2)}
        self._filling = True
        for c in keep:
            c.clear()
        if ds is not None:
            for combo in (self.x_axis, self.y_axis):
                for key in ds.derived:
                    if key[-1] in "12":
                        combo.addItem(key, _spec("derived", key))
                for ch in ds.channels:
                    combo.addItem(ds.label(ch), _spec("channel", ch))
            self.color_by.addItem("One colour", "")
            self.split_by.addItem("One panel", "")
            self.split_by2.addItem("—", "")
            keys = self._annotation_keys()
            for combo in (self.color_by, self.split_by, self.split_by2):
                if len(self.session.group_by) > 1 and combo is not self.split_by2:
                    combo.addItem(" + ".join(self.session.group_by), _spec("group", ""))
                for k in keys:
                    combo.addItem(k, _spec("meta", k))
                combo.addItem("Sample", _spec("sample", ""))
                if "FlowSOM" in ds.derived:
                    combo.addItem("FlowSOM cluster", _spec("cluster", "FlowSOM"))
                if "PhenoGraph" in ds.derived:
                    combo.addItem("PhenoGraph cluster", _spec("cluster", "PhenoGraph"))
            self.color_by.insertSeparator(self.color_by.count())
            for ch, m, comp in zip(ds.channels, ds.markers, ds.compensated):
                if comp:
                    self.color_by.addItem(m or ch, _spec("expression", ch))
            for combo, value in keep.items():
                i = combo.findData(value) if value else -1
                if i >= 0:
                    combo.setCurrentIndex(i)
            if self.y_axis.currentIndex() == self.x_axis.currentIndex() and self.y_axis.count() > 1:
                self.y_axis.setCurrentIndex(self.x_axis.currentIndex() + 1)
            if prefer in ("UMAP", "tSNE", "PCA"):
                self.x_axis.setCurrentIndex(self.x_axis.findData(_spec("derived", prefer + "1")))
                self.y_axis.setCurrentIndex(self.y_axis.findData(_spec("derived", prefer + "2")))
            if prefer in ("FlowSOM", "PhenoGraph"):
                self.color_by.setCurrentIndex(self.color_by.findData(_spec("cluster", prefer)))
            elif not self.color_by.currentData() and keys:
                i = self.color_by.findData(_spec("meta", self.session.group_by[0]
                                            if self.session.group_by else keys[0]))
                self.color_by.setCurrentIndex(max(i, 0))
        self._filling = False
        self._describe()
        self._draw()

    def _annotation_keys(self) -> list[str]:
        ds = self.session.dataset
        keys: list[str] = []
        for m in ds.meta:
            for k, v in m.items():
                if k not in keys:
                    keys.append(k)
        # a key every sample shares, or none fills in, says nothing
        return [k for k in keys if len({m.get(k, "") for m in ds.meta}) > 1]

    def _describe(self) -> None:
        ds = self.session.dataset
        if ds is None:
            self.result.setText("")
            return
        lines = [f"{len(ds):,d} events, {len(ds.samples)} samples"]
        for name in ("UMAP", "tSNE"):
            run = ds.runs.get(name)
            if run and "neighbours_kept" in run:
                lines.append(f"{name} keeps {run['neighbours_kept']:.0%} of each "
                             "event's 15 nearest neighbours")
        run = ds.runs.get("PhenoGraph")
        if run:
            lines.append(f"PhenoGraph: {run['n_clusters']} clusters, "
                         f"modularity {run['modularity']:.2f}")
        lines += self._components()
        self.result.setText("\n".join(lines))

    def _components(self) -> list[str]:
        """What each principal component is made of: its share of the
        variance and the three markers that weigh most in it."""
        ds = self.session.dataset
        run = ds.runs.get("PCA")
        if not run:
            return []
        names = [ds.label(c) for c in run["channels"]]
        lines = []
        for i, (share, row) in enumerate(zip(run["explained"], run["loadings"]), 1):
            heaviest = sorted(zip(row, names), key=lambda t: -abs(t[0]))[:3]
            lines.append(f"PCA{i}, {share:.0%}: " + ", ".join(
                f"{'+' if w >= 0 else '−'}{abs(w):.2f} {name}" for w, name in heaviest))
        return lines

    def _derived_label(self, key: str) -> str:
        """A principal component says how much of the variance it holds."""
        run = self.session.dataset.runs.get("PCA")
        if run and key.startswith("PCA"):
            return f"{key} ({run['explained'][int(key[3:]) - 1]:.0%})"
        return key

    # --- drawing -------------------------------------------------------------

    def _axis(self, spec):
        """(values, label, ticks). A marker is drawn on FlowJo's
        biexponential scale and ticked in decades; anything else is a
        plain number."""
        ds = self.session.dataset
        kind, key = _parse(spec)
        if kind == "derived":
            return ds.derived[key], self._derived_label(key), None
        x = ds.column(key).astype(np.float64)
        t = self._scale(key)
        if t is not None:
            return t(x), ds.label(key), decade_ticks(t)
        return x, ds.label(key), None

    def _scale(self, channel: str):
        """The scale chosen for the parameter in Gates; without a folder
        open — a dataset read from a file — the usual one."""
        ds = self.session.dataset
        if channel in self.session.scales or self.session.usable():
            return self.session.scale_of(channel)
        return Biex() if ds.compensated[ds.channels.index(channel)] else None

    def _labels(self, spec) -> Categories | None:
        ds = self.session.dataset
        kind, key = _parse(spec)
        if kind == "sample":
            return categories_of(ds.obs("sample"), "Sample", order=ds.samples)
        if kind == "meta":
            return categories_of(ds.obs(key), key)
        if kind == "group":
            groups = np.array([self.session.group_of(m) for m in ds.meta],
                              dtype=object)[ds.sample_index]
            return categories_of(groups, "Group")
        if kind == "cluster":
            return cluster_categories(self.session.clusters(key))
        return None

    def _draw(self) -> None:
        if self._filling:
            return
        ds = self.session.dataset
        xs, ys = self.x_axis.currentData(), self.y_axis.currentData()
        self.legend.clear()
        if ds is None or not xs or not ys:
            self.stack.setCurrentWidget(self.empty)
            self.view.clear_points()
            return
        self.stack.setCurrentWidget(self.view)
        x, xlabel, xticks = self._axis(xs)
        y, ylabel, yticks = self._axis(ys)
        spec = self.color_by.currentData()
        if not spec:
            color = None
        elif _parse(spec)[0] == "expression":
            ch = _parse(spec)[1]
            t = self._scale(ch) or Biex()
            color = Values(t(ds.column(ch).astype(np.float64)),
                           title=ds.label(ch), ticks=decade_ticks(t))
        else:
            color = self._labels(spec)
        split_spec, second = self.split_by.currentData(), self.split_by2.currentData()
        first = self._labels(split_spec) if split_spec else None
        then = self._labels(second) if second and second != split_spec else None
        split = figures.combine(first, then)
        # split twice, a row is one value of the first column and the
        # columns line up the values of the second
        columns = len(then.names) if first is not None and then is not None else None
        self.view.show_points(x, y, xlabel, ylabel, color, split, xticks, yticks, columns)
        if isinstance(color, Categories):
            for name, hexa, n in zip(color.names, color.colors, color.counts()):
                pix = QtGui.QPixmap(12, 12)
                pix.fill(QtGui.QColor(hexa))
                item = QtWidgets.QListWidgetItem(
                    QtGui.QIcon(pix), f"{name or '—'}    {int(n):,d}")
                self.legend.addItem(item)
        elif isinstance(color, Values):
            lo, hi = color.limits()
            self.legend.addItem(f"{color.title}: {lo:.2f} to {hi:.2f}")
            self.legend.addItem("Biexponential scale; the colours span "
                                "the 1st to the 99th percentile.")

    def make_figure(self) -> None:
        from .figure_dialog import open_for
        if open_for(self.view, self, "map") is None:
            self.sigStatus.emit("Nothing on the plot yet")

    # --- running -------------------------------------------------------------

    def run(self, names: list[str]) -> None:
        s = self.session
        channels = self.chosen_markers()
        if s.dataset is None:
            return
        if len(channels) < 2:
            self.sigStatus.emit("Tick at least two markers")
            return
        common = dict(channels=channels, transform=self.transform.currentData(),
                      seed=self.seed.value())
        jobs = {
            "UMAP": lambda: s.run_umap(n_neighbors=self.neighbors.value(),
                                       min_dist=self.min_dist.value(), **common),
            "tSNE": lambda: s.run_tsne(perplexity=self.perplexity.value(), **common),
            "FlowSOM": lambda: s.run_flowsom(n_clusters=self.clusters.value(),
                                             grid=self.grid.value(), **common),
            "PhenoGraph": lambda: s.run_phenograph(k=self.k.value(),
                                                   resolution=self.resolution.value(),
                                                   **common),
            "PCA": lambda: s.run_pca(channels, common["transform"]),
        }
        # the settings are read now, on this thread, not when the work starts
        for name in names:
            if name in self._running:
                continue
            self._running.add(name)
            self.run_buttons[name].setEnabled(False)
            self.run_buttons[name].setText(f"Running {name}…")
            self.sigStatus.emit(f"Running {', '.join(sorted(self._running))}…")
            worker.run(jobs[name], compute=True,
                       on_done=lambda _, n=name: self._finished(n, None),
                       on_failed=lambda msg, n=name: self._finished(n, msg))

    def _finished(self, name: str, error: str | None) -> None:
        self._running.discard(name)
        self.run_buttons[name].setEnabled(True)
        self.run_buttons[name].setText(f"Run {name}")
        if error:
            self.sigStatus.emit(f"{name}: {error}")
            return
        self.sigStatus.emit(f"Running {', '.join(sorted(self._running))}…"
                            if self._running else f"{name} done")
        self.results_changed(prefer=name)
        self.sigResults.emit()
