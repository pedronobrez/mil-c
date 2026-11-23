"""
The third tab: the clusters, from four sides at once. The map says where
a cluster is, the heatmap what it expresses, the bars in which samples
and groups it is found, the tree how the map's nodes relate. A click on a
cluster in any of them marks it in all four; a lasso on the map
(Shift+drag) asks the same questions of any set of events.
"""

from __future__ import annotations

import html

import numpy as np
from PyQt6 import QtCore, QtWidgets

from .. import cluster_stats, clusters, compare
from ..session import Session
from . import palette, style, theme
from .ramp_picker import RampPicker
from .cluster_tools import ClusterTools, cluster_categories
from .cluster_views import FrequencyView, HeatmapView, TreeView
from .scatter import Categories, MapView

METACLUSTER = "Metacluster"


class ClustersTab(QtWidgets.QWidget):
    sigStatus = QtCore.pyqtSignal(str)
    #: clusters were renamed, merged or drawn: whatever shows them is out of date
    sigClustersChanged = QtCore.pyqtSignal()
    #: a cluster became a population of the files: the gating trees changed
    sigGatesChanged = QtCore.pyqtSignal()

    def __init__(self, session: Session, parent=None):
        super().__init__(parent)
        self.session = session
        #: the clustering shown, a key of the dataset's results
        self.key = clusters.FLOWSOM
        self.selected: int | None = None
        self._lasso: np.ndarray | None = None
        self._filling = False
        self._tree = None

        self.embedding = QtWidgets.QComboBox()
        self.all_markers = QtWidgets.QCheckBox("Other markers too")
        self.ramp = RampPicker()
        self.ramp.setToolTip("The colours of the heatmap, of a tree coloured by a marker "
                             "and of expression on the maps")
        self.stat = QtWidgets.QComboBox()
        self.stat.addItem("Median", "median")
        self.stat.addItem("Mean", "mean")
        self.per_marker = QtWidgets.QCheckBox("Colour per marker")
        self.per_marker.setChecked(True)
        self.per_marker.setToolTip(
            "Each column runs over the whole ramp by itself, so a dim marker "
            "shows its own differences. Unticked, one ramp for the table.")
        self.group_by = QtWidgets.QComboBox()
        self.group_by.setToolTip("What the bars compare, by colour")
        self.facet_by = QtWidgets.QComboBox()
        self.facet_by.setToolTip(
            "A second subdivision: one panel for each of its values, so the groups "
            "are compared within each — the groups within each day, say")
        self.tree_color = QtWidgets.QComboBox()
        self.test = QtWidgets.QComboBox()
        self.test.addItem("No test", "")
        for key, text in compare.KINDS.items():
            self.test.addItem(text, key)
        self.test.setToolTip(
            "Compare the groups' share of each cluster, one value per sample. Two groups: "
            "Welch's t or Mann-Whitney. Three or more: ANOVA or Kruskal-Wallis. The p "
            "written above each cluster is adjusted for the number of clusters tested.")
        self.correction = QtWidgets.QComboBox()
        for name in compare.CORRECTIONS:
            self.correction.addItem(name, name)
        self.paired_by = QtWidgets.QComboBox()
        self.paired_by.setToolTip(
            "Where the same animals or donors were measured in every group: the column "
            "that says which one a sample is of. The test is then of each one's change — "
            "paired t or Wilcoxon for two groups, repeated-measures ANOVA or Friedman for "
            "more — and those not found in every group are left out.")
        self.pairing = QtWidgets.QLabel("")
        self.pairing.setProperty("role", "warning")
        # a row of its own: the bar above is already wider than a laptop's screen
        cap = QtWidgets.QLabel("Paired by")
        cap.setProperty("role", "caption")
        pairing_row = QtWidgets.QHBoxLayout()
        pairing_row.addWidget(cap)
        pairing_row.addWidget(self.paired_by)
        pairing_row.addSpacing(10)
        pairing_row.addWidget(self.pairing, 1)
        self.results_button = QtWidgets.QPushButton("Results…")
        self.results_button.setToolTip("Every comparison as a table, pair by pair where "
                                       "there are more than two groups")
        self.prism_button = QtWidgets.QPushButton("Prism…")
        self.prism_button.setToolTip(
            "The share of each cluster in each sample as a GraphPad Prism file: a grouped "
            "table, clusters down, groups across, a subcolumn for each sample — one table "
            "for each panel of the chart")
        self.abundance_button = QtWidgets.QPushButton("Differential abundance…")
        self.abundance_button.setToolTip(
            "Two of the groups compared on the clusters' cell counts, by a negative-binomial "
            "model in the manner of diffcyt and edgeR, with a volcano")
        self.clear_button = QtWidgets.QPushButton("Clear selection")

        bar = QtWidgets.QHBoxLayout()
        for label, widget in (("Map", self.embedding), ("Heatmap", self.stat),
                              ("Colours", self.ramp),
                              (None, self.per_marker), (None, self.all_markers),
                              ("Compare", self.group_by), ("within each", self.facet_by),
                              ("Tree", self.tree_color), ("Test", self.test),
                              (None, self.correction)):
            if label:
                cap = QtWidgets.QLabel(label)
                cap.setProperty("role", "caption")
                bar.addWidget(cap)
            bar.addWidget(widget)
            bar.addSpacing(10)
        bar.addWidget(self.results_button)
        bar.addWidget(self.prism_button)
        bar.addWidget(self.abundance_button)
        bar.addStretch(1)
        bar.addWidget(self.clear_button)

        self.map = MapView()
        self.heatmap = HeatmapView()
        self.freq = FrequencyView()
        self.tree = TreeView()
        top = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        top.addWidget(self._titled("Map", self.map, "cluster_map"))
        top.addWidget(self._titled("Expression", self.heatmap, "cluster_expression"))
        bottom = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        bottom.addWidget(self._titled("FlowSOM tree", self.tree, "flowsom_tree"))
        bottom.addWidget(self._titled("Frequency", self.freq, "cluster_frequency"))
        grid = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        grid.addWidget(top)
        grid.addWidget(bottom)
        for sp in (top, bottom, grid):
            sp.setChildrenCollapsible(False)
            sp.setSizes([1000, 1000])

        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)
        self.summary.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.empty = QtWidgets.QLabel("Run FlowSOM in the Explore tab: the clusters appear here.")
        self.empty.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.empty.setProperty("role", "hint")
        body = QtWidgets.QWidget()
        inner = QtWidgets.QVBoxLayout(body)
        inner.setContentsMargins(0, 0, 0, 0)
        inner.addLayout(bar)
        inner.addLayout(pairing_row)
        inner.addWidget(grid, 1)
        inner.addWidget(self.summary)
        self.tools = ClusterTools(self)
        inner.addWidget(self.tools)
        self.stack = QtWidgets.QStackedWidget()
        self.stack.addWidget(self.empty)
        self.stack.addWidget(body)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(16, 14, 16, 14)
        layout.addWidget(self.stack)

        for view in (self.heatmap, self.freq, self.tree):
            view.sigCluster.connect(self.select)
        self.map.sigSelected.connect(self._lassoed)
        self.clear_button.clicked.connect(self.clear_selection)
        self.embedding.currentIndexChanged.connect(self._draw_map)
        self.stat.currentIndexChanged.connect(self._draw_tables)
        self.per_marker.toggled.connect(self._draw_tables)
        self.all_markers.toggled.connect(self._draw_tables)
        self.group_by.currentIndexChanged.connect(self._draw_tables)
        self.group_by.currentIndexChanged.connect(self._summarise)
        self.facet_by.currentIndexChanged.connect(self._draw_tables)
        self.facet_by.currentIndexChanged.connect(self._summarise)
        self.tree_color.currentIndexChanged.connect(self._draw_tree)
        self.test.currentIndexChanged.connect(self._draw_tables)
        self.correction.currentIndexChanged.connect(self._draw_tables)
        self.paired_by.currentIndexChanged.connect(self._draw_tables)
        self.results_button.clicked.connect(self.show_results)
        self.prism_button.clicked.connect(self._ask_prism)
        self.abundance_button.clicked.connect(lambda: self.abundance_dialog().exec())
        palette.notifier.changed.connect(self._recolour)
        self.refresh()

    def _titled(self, title: str, widget: QtWidgets.QWidget, name: str) -> QtWidgets.QWidget:
        """A panel under its title, with the button that makes it a figure."""
        box = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(box)
        lay.setContentsMargins(0, 6, 0, 0)
        lay.setSpacing(4)
        label = QtWidgets.QLabel(title)
        label.setFont(style.eyebrow_font())
        label.setProperty("role", "caption")
        button = QtWidgets.QToolButton()
        button.setText("Figure…")
        button.setAutoRaise(True)
        button.clicked.connect(lambda _=False, v=widget, n=name: self.make_figure(v, n))
        head = QtWidgets.QHBoxLayout()
        head.addWidget(label)
        head.addStretch(1)
        head.addWidget(button)
        lay.addLayout(head)
        lay.addWidget(widget, 1)
        return box

    def make_figure(self, view, name: str) -> None:
        from .figure_dialog import open_for
        open_for(view, self, name)

    # --- from the session ------------------------------------------------------

    def refresh(self) -> None:
        """After a run, or a new dataset."""
        ds = self.session.dataset
        now = self.session.clusters(self.key)
        if now is None:
            self.stack.setCurrentIndex(0)
            self.selected, self._lasso, self._tree = None, None, None
            return
        self.stack.setCurrentIndex(1)
        if self.selected not in now.ids:
            self.selected = None
        if self._lasso is not None and self._lasso.size != len(ds):
            self._lasso = None
        self._tree = self.session.flowsom_tree()

        # the box used to do nothing whenever FlowSOM had been run on every
        # marker, which is what the Explore tab ticks to begin with; now it
        # says how many there are to add, or that there are none
        others = self._other_markers()
        self.all_markers.setEnabled(bool(others))
        self.all_markers.setText(f"Other markers too ({len(others)})" if others
                                 else "Other markers too (none)")
        self.all_markers.setToolTip(
            "Also show " + ", ".join(ds.label(c).split()[0] for c in others)
            + ": stained, but not used to build the clusters" if others else
            "FlowSOM was run on every stained marker, so there is nothing to add. "
            "Untick the lineage and viability markers in Explore and run it again "
            "to cluster on fewer.")
        self._filling = True
        keep = {c: c.currentData() for c in (self.embedding, self.group_by, self.facet_by,
                                             self.tree_color)}
        self.embedding.clear()
        for name in ("UMAP", "tSNE"):
            if name + "1" in ds.derived:
                self.embedding.addItem(name, name)
        run = ds.runs.get(self.key, {})
        pair = (run.get("channels") or self.session.default_channels())[:2]
        self.embedding.addItem(" × ".join(ds.label(c).split()[0] for c in pair),
                               "channels")
        self.group_by.clear()
        self.facet_by.clear()
        keys = [k for k in dict.fromkeys(k for m in ds.meta for k in m)
                if len({m.get(k, "") for m in ds.meta}) > 1]
        if len(self.session.group_by) > 1:
            self.group_by.addItem(" + ".join(self.session.group_by), "*")
        for k in keys:
            self.group_by.addItem(k, k)
        if self.group_by.count() == 0:
            self.group_by.addItem("All samples", "")
        self.facet_by.addItem("—", "")
        for k in keys:
            self.facet_by.addItem(k, k)
        self._offer_pairing(keys)
        self.tree_color.clear()
        self.tree_color.addItem(METACLUSTER, "")
        for m in (self._tree or {}).get("markers", []):
            self.tree_color.addItem(m.split()[0], m)
        # the grouping chosen for the sampling is the one asked about first
        first = (self.session.group_by or [""])[0]
        self.group_by.setCurrentIndex(max(self.group_by.findData(first), 0))
        # two columns chosen to group by: compare the first within each of the second
        second = self.session.group_by[1] if len(self.session.group_by) > 1 else ""
        self.facet_by.setCurrentIndex(max(self.facet_by.findData(second), 0))
        for combo, value in keep.items():
            i = combo.findData(value) if value else -1
            if i >= 0:
                combo.setCurrentIndex(i)
        self._filling = False
        self._draw_all()

    def _recolour(self) -> None:
        if self.stack.currentIndex() == 1:
            self._draw_all()

    def _other_markers(self) -> list[str]:
        """Stained markers FlowSOM was not run on."""
        ds = self.session.dataset
        used = set(ds.runs.get(self.key, {}).get("channels") or [])
        return [c for c, m, comp in zip(ds.channels, ds.markers, ds.compensated)
                if comp and m and c not in used]

    def _offer_pairing(self, keys: list[str]) -> None:
        """The columns the samples can be paired by, keeping the one chosen."""
        chosen = self.paired_by.currentData()
        self.paired_by.clear()
        self.paired_by.addItem("—", "")
        for k in keys:
            self.paired_by.addItem(k, k)
        self.paired_by.setCurrentIndex(max(self.paired_by.findData(chosen), 0))

    def _facets(self) -> list[str] | None:
        key = self.facet_by.currentData()
        if not key or key == self.group_by.currentData():
            return None
        return [m.get(key, "") for m in self.session.dataset.meta]

    def _summary_groups(self) -> list[str]:
        """The groups a selection is described by: every combination of
        what is compared and what it is compared within."""
        groups, facets = self._groups(), self._facets()
        if not facets:
            return groups
        return [f"{g or '—'} · {f or '—'}" for g, f in zip(groups, facets)]

    def _groups(self) -> list[str]:
        ds = self.session.dataset
        key = self.group_by.currentData()
        if key == "*":
            return self.session.sample_groups()
        if not key:
            return ["All samples"] * len(ds.samples)
        return [m.get(key, "") for m in ds.meta]

    # --- selection ---------------------------------------------------------------

    def select(self, cluster: int) -> None:
        """A second click on the selected cluster lets it go."""
        self.selected = None if cluster == self.selected else int(cluster)
        self._lasso = None
        self._draw_all()

    def clear_selection(self) -> None:
        self.selected, self._lasso = None, None
        self._draw_all()

    def _lassoed(self, mask: np.ndarray) -> None:
        self._lasso = mask if mask.any() else None
        self.selected = None
        self._draw_all()

    def selection_mask(self) -> np.ndarray | None:
        if self._lasso is not None:
            return self._lasso
        if self.selected is not None:
            return self.session.clusters(self.key).labels == self.selected
        return None

    def clusters_changed(self, select: int | None = None) -> None:
        """The clusters were edited: show them as they now are, here and
        wherever else they are shown."""
        self.selected, self._lasso = select, None
        self.refresh()
        self.sigClustersChanged.emit()

    # --- drawing -------------------------------------------------------------------

    def _draw_all(self) -> None:
        if self._filling:
            return
        self._draw_map()
        self._draw_tables()
        self._draw_tree()
        self._summarise()
        self.tools.sync()

    def _draw_map(self) -> None:
        ds = self.session.dataset
        now = self.session.clusters(self.key)
        if self._filling or now is None:
            return
        which = self.embedding.currentData()
        if which in ("UMAP", "tSNE"):
            x, y = ds.derived[which + "1"], ds.derived[which + "2"]
            xl, yl = which + "1", which + "2"
        else:
            run = ds.runs.get(self.key, {})
            a, b = (run.get("channels") or self.session.default_channels())[:2]
            X = self.session.features([a, b], run.get("transform", "biex"))
            x, y, xl, yl = X[:, 0], X[:, 1], ds.label(a), ds.label(b)
        cats = cluster_categories(now)
        if self.selected is not None:
            # the others stay where they are, greyed: the cluster is read
            # against the whole map
            cats = Categories(cats.codes, cats.names,
                              [c if i == self.selected else theme.faint_axis()
                               for c, i in zip(cats.colors, now.ids)], cats.title)
        self.map.selected = self._lasso
        self.map.show_points(x, y, xl, yl, cats)

    def _draw_tables(self) -> None:
        if self._filling or self.session.dataset is None:
            return
        ds = self.session.dataset
        used = ds.runs.get(self.key, {}).get("channels") or self.session.default_channels()
        show_others = self.all_markers.isChecked() and self.all_markers.isEnabled()
        channels = list(used) + self._other_markers() if show_others else None
        tables = self.session.cluster_tables(channels, self.stat.currentData(), self.key)
        if tables is None:
            return
        expression, freq = tables
        self.heatmap.show_table(expression, self.per_marker.isChecked(), self.selected)
        self._frequency = freq
        self.correction.setEnabled(bool(self.test.currentData()))
        self.paired_by.setEnabled(bool(self.test.currentData()))
        self.pairing.setText(self.pairing_problem())
        self.results_button.setEnabled(bool(self.test.currentData()))
        # on the plot, what a reader looks for: the p where there is a
        # difference, "ns" where there is not. Every number is in Results.
        notes = {c.key: compare.label(c.p_adjusted) if c.p_adjusted < 0.05 else "ns"
                 for c in self.comparisons() if c.p_adjusted == c.p_adjusted}
        self.freq.show_frequencies(freq, self._groups(), self.selected, self._facets(), notes)

    def _values(self) -> dict:
        """(panel, cluster) → group → one share per sample."""
        freq = self._frequency
        groups, facets = self._groups(), self._facets()
        out: dict = {}
        for i, cluster in enumerate(freq.rows):
            for j, group in enumerate(groups):
                value = freq.values[i, j]
                if value == value:
                    panel = facets[j] if facets else None
                    out.setdefault((panel, cluster), {}).setdefault(group, []).append(float(value))
        return out

    def _subjects(self) -> dict | None:
        """As `_values`, whose each share is: the animal or donor of its
        sample, in the column chosen to pair by. None where none is."""
        column = self.paired_by.currentData()
        if not column:
            return None
        freq, facets = self._frequency, self._facets()
        whose = [m.get(column, "") for m in self.session.dataset.meta]
        out: dict = {}
        for i, cluster in enumerate(freq.rows):
            for j, group in enumerate(self._groups()):
                if freq.values[i, j] == freq.values[i, j]:
                    panel = facets[j] if facets else None
                    out.setdefault((panel, cluster), {}).setdefault(group, []).append(whose[j])
        return out

    def pairing_problem(self) -> str:
        """Why the samples cannot be paired by the column chosen; "" where
        they can, or where no test or no pairing is asked for."""
        if not self.test.currentData() or getattr(self, "_frequency", None) is None:
            return ""
        subjects = self._subjects()
        problem = compare.pairing_problem(self._values(), subjects) if subjects else ""
        return (f"No test: the samples cannot be paired by {self.paired_by.currentData()} — "
                f"{problem}.") if problem else ""

    def comparisons(self, pairwise: bool = False) -> list:
        """The groups compared in each cluster — in each panel, where the
        chart has panels — with the test and the correction chosen."""
        kind = self.test.currentData()
        if not kind or getattr(self, "_frequency", None) is None:
            return []
        if self.pairing_problem():   # nothing is computed from a pairing that is not one
            return []
        return compare.compare(self._values(), kind, self.correction.currentData(), pairwise,
                               self._subjects())

    def prism_tables(self) -> list:
        from .. import prism
        return prism.from_frequencies(self._frequency, self._groups(), self._facets(),
                                      "Cluster frequency (% of sample)")

    def save_prism(self, path: str) -> int:
        from .. import prism
        tables = self.prism_tables()
        prism.write_pzfx(path, tables)
        self.sigStatus.emit(f"Wrote {len(tables)} Prism table(s) to {path}")
        return len(tables)

    def _ask_prism(self) -> None:
        if getattr(self, "_frequency", None) is None:
            return
        import os
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save for GraphPad Prism",
            os.path.join(self.session.folder, "cluster_frequency.pzfx"), "GraphPad Prism (*.pzfx)")
        if path:
            self.save_prism(path if path.endswith(".pzfx") else path + ".pzfx")

    def results_table(self):
        found = self.comparisons(pairwise=True)
        panels = self._facets() is not None
        table = compare.to_table(found, [self.facet_by.currentText(), "Cluster"] if panels
                                 else ["Cluster"], self.correction.currentData(),
                                 self.paired_by.currentData() or None)
        if not panels:
            for row in table.rows:   # the empty panel is not a column
                del row[0]
        return table, found

    def show_results(self) -> None:
        from .results_dialog import NOTE, ResultsDialog
        table, found = self.results_table()
        if not found:
            return
        note = NOTE.format(tests=sum(1 for c in found if c.p == c.p),
                           correction=self.correction.currentData())
        ResultsDialog(table, "Clusters between groups",
                      "The share of each sample's events in each cluster, compared between "
                      f"{self.group_by.currentText().lower()}s. " + note, self,
                      "cluster_comparison").exec()

    def abundance_dialog(self):
        """Differential abundance of what the tab compares now: the
        dialog, not yet shown."""
        from .abundance_dialog import AbundanceDialog
        return AbundanceDialog.for_tab(self)

    def _draw_tree(self) -> None:
        if self._filling:
            return
        self.tree.show_tree(self._tree, self.tree_color.currentData() or None,
                            self.selected)

    def _summarise(self, *_) -> None:
        if self._filling:
            return
        ds = self.session.dataset
        mask = self.selection_mask()
        if ds is None or mask is None:
            self.summary.setText(
                "<span style='color:%s'>Click a cluster in any panel to follow it "
                "through the others. Shift+drag on the map to outline any set "
                "of events.</span>" % theme.muted())
            return
        info = cluster_stats.selection(mask, ds.sample_index, ds.samples,
                                       self._summary_groups())
        now = self.session.clusters(self.key)
        what = (f"Cluster {html.escape(now.name(self.selected))}" if self._lasso is None
                else "Selection")
        parts = [f"<b>{what}</b>: {info['n']:,d} events, {info['share']:.1%} of all"]
        if self._lasso is None and now.note(self.selected):
            parts.append(f"<i>{html.escape(now.note(self.selected))}</i>")
        for g, (mean, sd, n) in info["by_group"].items():
            spread = f" ± {sd:.1f}" if n > 1 else ""
            parts.append(f"{g or '—'} {mean:.1f}{spread}% (n={n})")
        if self._lasso is not None:
            labels = now.labels[mask]
            ids, counts = np.unique(labels, return_counts=True)
            top = np.argsort(-counts)[:3]
            parts.append("clusters " + ", ".join(
                f"{html.escape(now.name(int(ids[i])))} ({counts[i] / labels.size:.0%})"
                for i in top))
        self.summary.setText(" &nbsp;·&nbsp; ".join(parts))
