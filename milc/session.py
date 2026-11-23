"""
The open analysis: a folder of files, the workspace that gates them, how
they are grouped and sampled, and the concatenated dataset with what has
been computed on it. The window reads and writes this and nothing else;
it never imports Qt.
"""

from __future__ import annotations

import glob
import os
from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from . import cluster_stats, sampling
from .dataset import Dataset, concatenate
from .fcs import read_header
from .sample import Sample, annotate
from .transforms import Arcsinh, Biex
from .wsp import Workspace, WspSample, read_wsp

#: `14Case4` → day 14, group Case, animal 4; a trailing ` marks a re-run
DEFAULT_PATTERN = r"(?P<day>\d+)(?P<group>[A-Za-z]+)(?P<animal>\d+)(?P<rerun>`)?$"

#: every event of the file, no gate
ALL_EVENTS = ""


@dataclass
class Entry:
    """One file of the folder, known from its header until it is needed."""
    path: str
    events: int
    channels: list[str]
    markers: list[str]
    meta: dict[str, str] = field(default_factory=dict)
    include: bool = True
    wsp: WspSample | None = None
    #: the workspace does not list this file; it is given the gates and
    #: the matrix all the listed ones share
    borrowed: bool = False
    #: set when the file cannot go in with the others, and why
    problem: str = ""
    sample: Sample | None = field(default=None, repr=False)

    @property
    def name(self) -> str:
        return os.path.splitext(os.path.basename(self.path))[0]

    @property
    def gated(self) -> bool:
        return self.borrowed or (self.wsp is not None and bool(self.wsp.gates.children))


class Session:
    def __init__(self):
        self.folder = ""
        self.workspace: Workspace | None = None
        self.entries: list[Entry] = []
        self.pattern = DEFAULT_PATTERN
        self.group_by: list[str] = ["group"]
        self.population = ALL_EVENTS
        self.mode = sampling.EQUAL_PER_SAMPLE
        self.n: int | None = None
        self.min_events = 1
        self.seed = 0
        self.dataset: Dataset | None = None
        #: give the files the workspace does not list the gates and the
        #: matrix the listed ones share, where they do share them
        self.extend_gates = True
        self._shared = None
        self._counts: dict[tuple[str, str], int] = {}
        #: the last FlowSOM, kept for its map and tree
        self.flowsom = None

    # --- files -----------------------------------------------------------

    def open_folder(self, folder: str) -> None:
        self.folder = folder
        self.entries, self._counts, self.dataset = [], {}, None
        for path in sorted(glob.glob(os.path.join(folder, "*.fcs"))):
            try:
                kv = read_header(path)
            except Exception as exc:  # a broken file is listed, not fatal
                self.entries.append(Entry(path, 0, [], [], include=False,
                                          problem=f"cannot be read: {exc}"))
                continue
            n = int(kv["$PAR"])
            self.entries.append(Entry(
                path, int(kv.get("$TOT", 0) or 0),
                [kv[f"$P{i}N"].strip() for i in range(1, n + 1)],
                [kv.get(f"$P{i}S", "").strip() for i in range(1, n + 1)]))
        # the files of an experiment share their channels; one that does
        # not (another program's concatenated export, say) stays out
        usual = Counter(tuple(e.channels) for e in self.entries if e.channels)
        if usual:
            common = usual.most_common(1)[0][0]
            for e in self.entries:
                if e.channels and tuple(e.channels) != common:
                    e.problem, e.include = "different channels from the others", False
        self.annotate()
        self.workspace = None
        found = sorted(glob.glob(os.path.join(folder, "*.wsp")))
        if found:
            self.set_workspace(found[0])

    def set_workspace(self, path: str | None) -> None:
        self.workspace = read_wsp(path) if path else None
        by_file = ({w.filename: w for w in self.workspace.samples}
                   if self.workspace else {})
        self._shared = self.workspace.shared() if self.workspace else None
        for e in self.entries:
            e.wsp = by_file.get(os.path.basename(e.path))
            e.sample = None
        self._borrow()
        paths = self.populations()
        # the deepest gate every gated sample has is what was meant
        self.population = paths[-1] if paths else ALL_EVENTS

    def _borrow(self) -> None:
        for e in self.entries:
            e.borrowed = bool(self.extend_gates and self._shared is not None
                              and e.wsp is None and not e.problem)
            if e.wsp is None:
                e.sample = None
        self._counts.clear()

    def can_extend(self) -> bool:
        """Is there one set of gates and one matrix to give the others?"""
        return self._shared is not None

    def set_extend_gates(self, on: bool) -> None:
        self.extend_gates = bool(on)
        self._borrow()
        paths = self.populations()
        if self.population not in paths:
            self.population = paths[-1] if paths else ALL_EVENTS

    def gates_of(self, e: Entry):
        if e.borrowed:
            return self._shared[0]
        return e.wsp.gates if e.wsp is not None else None

    def annotate(self, pattern: str | None = None) -> None:
        if pattern is not None:
            self.pattern = pattern
        for e in self.entries:
            try:
                e.meta = annotate(e.name, self.pattern)
            except Exception:
                e.meta = {}

    def repeats(self) -> list[list[str]]:
        """Included files that carry the same annotations apart from the
        re-run mark: one animal acquired twice, which would be counted as
        two."""
        seen: dict[tuple, list[str]] = {}
        for e in self.entries:
            if e.include and not e.problem and e.meta:
                key = tuple(sorted((k, v) for k, v in e.meta.items() if k != "rerun"))
                seen.setdefault(key, []).append(e.name)
        return [names for names in seen.values() if len(names) > 1]

    def meta_keys(self) -> list[str]:
        keys: list[str] = []
        for e in self.entries:
            for k in e.meta:
                if k not in keys:
                    keys.append(k)
        return keys

    def usable(self) -> list[Entry]:
        return [e for e in self.entries if not e.problem]

    def populations(self) -> list[str]:
        """Gate paths every gated file has, shallowest first; without the
        ungated root."""
        gated = [e for e in self.usable() if e.gated]
        if not gated:
            return []
        sets = [[n.path for n in self.gates_of(e).walk() if n.gate is not None]
                for e in gated]
        return [p for p in sets[0] if all(p in s for s in sets[1:])]

    def available(self, e: Entry) -> bool:
        """Can this file give events of the chosen population?"""
        return not e.problem and (self.population == ALL_EVENTS or e.gated)

    def group_of(self, meta: dict[str, str]) -> str:
        return " ".join(meta.get(k, "") for k in self.group_by).strip()

    # --- events ----------------------------------------------------------

    def load(self, e: Entry) -> Sample:
        if e.sample is None:
            if e.borrowed:
                gates, matrix = self._shared
                e.sample = Sample.open(e.path, compensation=matrix, gates=gates)
            elif e.wsp is not None and e.wsp.compensation is not None:
                e.sample = Sample.open(e.path, compensation=e.wsp.compensation,
                                       gates=e.wsp.gates, transforms=e.wsp.transforms)
            else:
                e.sample = Sample.open(e.path)
            e.sample.meta = dict(e.meta)
        return e.sample

    def count(self, e: Entry) -> int | None:
        """Events of the chosen population; None where it does not exist
        for this file. Reads the file the first time."""
        if not self.available(e):
            return None
        if self.population == ALL_EVENTS:
            return e.events
        key = (e.path, self.population)
        if key not in self._counts:
            self._counts[key] = self.load(e).count(self.population)
        return self._counts[key]

    def counted(self, e: Entry) -> bool:
        return (self.population == ALL_EVENTS
                or (e.path, self.population) in self._counts)

    def count_all(self, progress=None) -> None:
        todo = [e for e in self.entries if self.available(e)]
        for i, e in enumerate(todo):
            self.count(e)
            if progress:
                progress(i + 1, len(todo), e.name)

    def plan(self) -> sampling.Plan:
        """From what has been counted so far."""
        chosen = [e for e in self.entries
                  if e.include and self.available(e) and self.counted(e)]
        return sampling.plan(
            {e.name: self.count(e) for e in chosen},
            {e.name: self.group_of(e.meta) for e in chosen},
            self.mode, self.n, self.min_events)

    def build(self, progress=None) -> Dataset:
        self.count_all(progress)
        p = self.plan()
        samples = []
        for e in self.entries:
            if p.take.get(e.name, 0) > 0:
                s = self.load(e)
                s.meta = dict(e.meta)
                samples.append(s)
        self.dataset = concatenate(samples, p, self.population, self.seed)
        self.flowsom = None
        # the events are in the dataset now; the files need not stay in memory
        for e in self.entries:
            e.sample = None
        return self.dataset

    # --- analysis --------------------------------------------------------

    def default_channels(self) -> list[str]:
        """The stained fluorescence channels: what a first analysis is
        usually run on. The user unticks lineage and viability dyes."""
        ds = self.dataset
        return [c for c, m, comp in zip(ds.channels, ds.markers, ds.compensated)
                if comp and m]

    def features(self, channels: list[str], transform: str = "biex",
                 cofactor: float = 150.0) -> np.ndarray:
        if transform == "biex":
            return self.dataset.matrix(channels, Biex(), 1023.0)
        return self.dataset.matrix(channels, Arcsinh(cofactor))

    def run_umap(self, channels, transform="biex", cofactor=150.0,
                 n_neighbors=15, min_dist=0.5, seed=0) -> None:
        from . import embed
        X = self.features(channels, transform, cofactor)
        Y = embed.umap(X, n_neighbors, min_dist, seed=seed)
        ds = self.dataset
        ds.derived["UMAP1"], ds.derived["UMAP2"] = Y[:, 0], Y[:, 1]
        ds.runs["UMAP"] = {
            "channels": list(channels), "transform": transform, "seed": seed,
            "n_neighbors": n_neighbors, "min_dist": min_dist,
            "neighbours_kept": embed.neighbour_preservation(X, Y)}

    def run_tsne(self, channels, transform="biex", cofactor=150.0,
                 perplexity=30.0, seed=0) -> None:
        from . import embed
        X = self.features(channels, transform, cofactor)
        Y = embed.tsne(X, perplexity, seed=seed)
        ds = self.dataset
        ds.derived["tSNE1"], ds.derived["tSNE2"] = Y[:, 0], Y[:, 1]
        ds.runs["tSNE"] = {
            "channels": list(channels), "transform": transform, "seed": seed,
            "perplexity": perplexity,
            "neighbours_kept": embed.neighbour_preservation(X, Y)}

    def run_flowsom(self, channels, transform="biex", cofactor=150.0,
                    n_clusters=5, grid=10, seed=0) -> None:
        from . import flowsom
        X = self.features(channels, transform, cofactor)
        fs = flowsom.flowsom(X, n_clusters, grid, grid, seed=seed)
        ds = self.dataset
        ds.derived["FlowSOM"] = fs.labels.astype(np.int32)
        ds.derived["FlowSOM_node"] = fs.node.astype(np.int32)
        ds.runs["FlowSOM"] = {
            "channels": list(channels), "transform": transform, "seed": seed,
            "n_clusters": n_clusters, "grid": grid}
        self.flowsom = fs

    def cluster_tables(self, channels: list[str] | None = None,
                       stat: str = "median"):
        """(expression, frequency per sample) of the FlowSOM clusters, or
        None before there are any. The markers are the ones FlowSOM ran
        on unless others are asked for."""
        ds = self.dataset
        if ds is None or "FlowSOM" not in ds.derived:
            return None
        run = ds.runs.get("FlowSOM", {})
        channels = channels or run.get("channels") or self.default_channels()
        X = self.features(channels, run.get("transform", "biex"))
        labels = ds.derived["FlowSOM"]
        return (cluster_stats.expression(X, labels, [ds.label(c) for c in channels], stat),
                cluster_stats.frequencies(labels, ds.sample_index, ds.samples))

    def sample_groups(self) -> list[str]:
        """The group of each sample of the dataset, in its order."""
        return [self.group_of(m) for m in self.dataset.meta]

    def flowsom_tree(self):
        """What the tree is drawn from, or None: positions and edges, each
        node's events and metacluster, and its median of every marker the
        map was built on (nodes x markers, with their names)."""
        from . import flowsom
        ds = self.dataset
        if ds is None or "FlowSOM_node" not in ds.derived:
            return None
        run = ds.runs.get("FlowSOM", {})
        channels = run.get("channels") or self.default_channels()
        X = self.features(channels, run.get("transform", "biex"))
        node, labels = ds.derived["FlowSOM_node"], ds.derived["FlowSOM"]
        n_nodes = (self.flowsom.codes.shape[0] if self.flowsom is not None
                   else int(node.max()) + 1)
        sizes = np.bincount(node, minlength=n_nodes)
        medians = np.full((n_nodes, len(channels)), np.nan)
        cluster = np.zeros(n_nodes, dtype=np.int64)
        for i in np.flatnonzero(sizes):
            rows = node == i
            medians[i] = np.median(X[rows], axis=0)
            cluster[i] = labels[rows][0]
        # a dataset opened from a file has no trained map: the nodes'
        # own medians stand in for its codes
        codes = (self.flowsom.codes if self.flowsom is not None
                 else np.nan_to_num(medians))
        pos, edges = flowsom.tree(codes, sizes)
        return {"pos": pos, "edges": edges, "sizes": sizes, "cluster": cluster,
                "medians": medians, "markers": [ds.label(c) for c in channels]}
