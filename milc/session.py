"""
The open analysis: a folder of files, the workspace that gates them, how
they are grouped and sampled, and the concatenated dataset with what has
been computed on it. The window reads and writes this and nothing else;
it never imports Qt.
"""

from __future__ import annotations

import glob
import json
import os
import re
from collections import Counter
from dataclasses import dataclass, field

import numpy as np

from . import cluster_stats, sampling, stats
from .dataset import Dataset, concatenate
from .compensation import Compensation, strip_prefix
from .fcs import read_header, write_fcs
from .gating import GateNode, copy_tree, follow_scale, tree_from_dict, tree_to_dict
from .sample import Sample, annotate
from .transforms import Arcsinh, Biex, from_spec, to_spec
from .wsp import Workspace, WspSample, read_wsp

#: the one column a new folder starts with; every other is the user's
GROUP = "Group"
#: the table's own columns, which a column of the user's cannot be called
RESERVED = ("Sample", "Events", "In gate", "Takes", "Note")

#: every event of the file, no gate
ALL_EVENTS = ""

#: what the user did to a folder — columns and their values, which files
#: are left out, gates drawn — kept beside the files it is about
PROJECT_FILE = "MIL-C.json"
#: where gates alone were kept, for a day, before there was more to keep
OLD_GATES_FILE = "MIL-C gates.json"


def guess_groups(names: list[str]) -> list[str]:
    """A first guess at the group of each file from its name: the run of
    letters that takes few values across the files (`14Case4`, `3Ctrl1` →
    Case, Ctrl). Nothing where no run of letters does; it is a starting
    point for a column the user can type over, not a rule."""
    tokens = [re.findall(r"[A-Za-z]+", n) for n in names]
    if not tokens or not all(tokens):
        return [""] * len(names)
    best = None
    for i in range(min(len(t) for t in tokens)):
        values = [t[i] for t in tokens]
        k = len(set(values))
        if 1 < k <= max(2, len(names) // 2) and (best is None or k < best[0]):
            best = (k, values)
    return best[1] if best else [""] * len(names)


def suggest_pattern(name: str) -> str:
    """A pattern with a named group for every run of digits or letters in
    a file name, to be edited: `14Case4` → number, Group, number."""
    parts, n, letters = [], 0, 0
    for token in re.findall(r"\d+|[A-Za-z]+|[^A-Za-z\d]+", name):
        if token.isdigit():
            n += 1
            parts.append(f"(?P<Number{n}>\\d+)")
        elif token.isalpha():
            letters += 1
            parts.append("(?P<%s>[A-Za-z]+)" % (GROUP if letters == 1 else f"Text{letters}"))
        else:
            parts.append(re.escape(token))
    # a mark at the end of one name (a re-run's backtick) is not in the others
    while parts and not parts[-1].startswith("(?P<"):
        parts.pop()
    return "".join(parts)


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
    #: this file's own gating tree, to look at and to edit: a copy of the
    #: workspace's, of the group's, of what was saved — or empty
    gates: GateNode = field(default_factory=lambda: GateNode(""))
    #: the workspace does not list this file; it was given the gates and
    #: the matrix all the listed ones share
    borrowed: bool = False
    #: the gates are the user's own, saved beside the files
    edited: bool = False
    #: set when the file cannot go in with the others, and why
    problem: str = ""
    sample: Sample | None = field(default=None, repr=False)

    @property
    def name(self) -> str:
        return os.path.splitext(os.path.basename(self.path))[0]

    @property
    def gated(self) -> bool:
        return bool(self.gates.children)


class Session:
    def __init__(self):
        self.folder = ""
        self.workspace: Workspace | None = None
        self.entries: list[Entry] = []
        #: the annotation columns, in order: free text per sample, to group,
        #: balance, colour and split by
        self.columns: list[str] = [GROUP]
        self.group_by: list[str] = [GROUP]
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
        #: channel → the scale its axis is drawn on, where the user chose
        #: one; a channel not here has the usual one (see `scale_of`)
        self.scales: dict[str, object] = {}
        #: a compensation matrix edited here, in place of every other
        self.matrix: object = None
        #: batch normalisation, or None: {"batch": column, "reference":
        #: [column, value], "clusters": n}
        self.normalization: dict | None = None
        #: the last normalisation, to look at: the held-out reference events
        #: of each batch before and after, on each marker's scale
        self.normalization_check: dict | None = None
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
        self.columns, self.group_by = [GROUP], [GROUP]
        self.scales, self.matrix, self.normalization = {}, None, None
        usable = [e for e in self.entries if not e.problem]
        for e, group in zip(usable, guess_groups([e.name for e in usable])):
            e.meta = {GROUP: group}
        self._project = self._read_project()
        self._restore_annotations()
        found = sorted(glob.glob(os.path.join(folder, "*.wsp")))
        self.set_workspace(found[0] if found else None)

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
        """Every file its own gating tree: what was saved beside the files
        if anything was, else the workspace's, else the group's, else none."""
        saved = self._saved_gates()
        for e in self.entries:
            e.sample = None
            e.borrowed = e.edited = False
            if e.problem:
                e.gates = GateNode(e.name)
            elif os.path.basename(e.path) in saved:
                e.gates, e.edited = saved[os.path.basename(e.path)], True
            elif e.wsp is not None and e.wsp.gates.children:
                e.gates = copy_tree(e.wsp.gates)
            elif self.extend_gates and self._shared is not None and e.wsp is None:
                e.gates, e.borrowed = copy_tree(self._shared[0], keep_counts=False), True
            else:
                e.gates = GateNode(e.name)
        self._counts.clear()

    def can_extend(self) -> bool:
        """Is there one set of gates and one matrix to give the others?"""
        return self._shared is not None

    def set_extend_gates(self, on: bool) -> None:
        self.extend_gates = bool(on)
        self._borrow()
        self._keep_population()

    def _keep_population(self) -> None:
        paths = self.populations()
        if self.population != ALL_EVENTS and self.population not in paths:
            self.population = paths[-1] if paths else ALL_EVENTS

    def gates_of(self, e: Entry) -> GateNode:
        return e.gates

    def matrix_of(self, e: Entry):
        """The compensation this file is given: the workspace's for it,
        the one the workspace's samples share where it lists none, and
        otherwise None — the file's own."""
        if self.matrix is not None:
            return self.matrix
        if e.wsp is not None and e.wsp.compensation is not None:
            return e.wsp.compensation
        if self._shared is not None and self.extend_gates:
            return self._shared[1]
        return None

    # --- compensation ------------------------------------------------------------

    def current_matrix(self):
        """(matrix, where it comes from) as the first usable file has it:
        `edited`, `workspace` or `file`. None where there is none at all."""
        entries = self.usable()
        if not entries:
            return None, ""
        first = entries[0]
        given = self.matrix_of(first)
        if given is not None:
            return given, "edited" if self.matrix is not None else "workspace"
        from . import compensation
        own = compensation.from_keywords(read_header(first.path))
        return own, "file" if own is not None else ""

    def set_matrix(self, matrix) -> None:
        """Compensate every file with this matrix — or, with None, go
        back to the workspace's or the files' own. Everything computed
        from compensated values is forgotten."""
        self.matrix = matrix
        for e in self.entries:
            e.sample = None
        self._counts.clear()
        self.save_project()

    # --- scales ----------------------------------------------------------------------

    def fluorescence(self) -> set[str]:
        """The channels a compensation matrix covers."""
        matrix, _ = self.current_matrix()
        return set(matrix.channels) if matrix is not None else set()

    def scale_of(self, name: str):
        """The scale a parameter's axis is drawn on: the one chosen for
        it, else FlowJo's biexponential for fluorescence and none — the
        values as they are — for scatter and time. `name` may carry the
        `Comp-` of a gate's dimension."""
        channel = strip_prefix(name)[0]
        if channel in self.scales:
            return self.scales[channel]
        return Biex() if channel in self.fluorescence() else None

    def set_scale(self, channels: list[str], scale) -> None:
        """Give these parameters this scale (None: linear, as stored).
        Polygons drawn on them follow; counts are forgotten."""
        for channel in channels:
            channel = strip_prefix(channel)[0]
            self.scales[channel] = scale
            for e in self.usable():
                for dim in (channel, "Comp-" + channel):
                    if follow_scale(e.gates, dim, scale) and e.sample is not None:
                        e.sample._masks = None
        self._counts.clear()
        self.save_project()

    def reset_scale(self, channels: list[str]) -> None:
        """Back to the usual scale for these parameters."""
        for channel in channels:
            channel = strip_prefix(channel)[0]
            self.scales.pop(channel, None)
            self.set_scale([channel], self.scale_of(channel))
            self.scales.pop(channel, None)   # the usual one is not a choice to keep
        self.save_project()

    def project_file(self) -> str:
        return os.path.join(self.folder, PROJECT_FILE)

    def _read_project(self) -> dict:
        try:
            with open(self.project_file(), encoding="utf-8") as f:
                data = json.load(f)
            return data if isinstance(data, dict) else {}
        except (OSError, ValueError):
            pass
        try:  # gates drawn before the file held anything else
            with open(os.path.join(self.folder, OLD_GATES_FILE), encoding="utf-8") as f:
                return {"gates": json.load(f).get("samples", {})}
        except (OSError, ValueError, AttributeError):
            return {}

    def _saved_gates(self) -> dict[str, GateNode]:
        try:
            return {name: tree_from_dict(tree)
                    for name, tree in getattr(self, "_project", {}).get("gates", {}).items()}
        except (KeyError, TypeError, ValueError):
            return {}

    def _restore_annotations(self) -> None:
        data = self._project
        self.scales = {c: from_spec(spec) for c, spec in data.get("scales", {}).items()}
        self.normalization = data.get("normalization") or None
        saved = data.get("compensation")
        if saved:
            try:
                self.matrix = Compensation(list(saved["channels"]),
                                           np.array(saved["matrix"], dtype=float), "edited")
            except (KeyError, TypeError, ValueError):
                self.matrix = None
        if not data.get("columns"):
            return
        self.columns = [str(c) for c in data["columns"]]
        self.group_by = [c for c in data.get("group_by", []) if c in self.columns] \
            or self.columns[:1]
        notes, left_out = data.get("annotations", {}), set(data.get("excluded", []))
        for e in self.entries:
            file = os.path.basename(e.path)
            if file in notes:
                e.meta = {c: str(notes[file].get(c, "")) for c in self.columns}
            if file in left_out:
                e.include = False

    def save_project(self) -> None:
        """Written when the user changes something, never just for opening
        a folder: a folder that was only looked at is left as it was."""
        usable = self.usable()
        data = {
            "version": 2, "columns": self.columns, "group_by": self.group_by,
            "annotations": {os.path.basename(e.path): {c: e.meta.get(c, "") for c in self.columns}
                            for e in usable},
            "excluded": [os.path.basename(e.path) for e in usable if not e.include],
            "gates": {os.path.basename(e.path): tree_to_dict(e.gates)
                      for e in usable if e.edited},
            "scales": {c: to_spec(t) for c, t in self.scales.items()},
            "normalization": self.normalization,
            "compensation": ({"channels": self.matrix.channels,
                              "matrix": np.asarray(self.matrix.matrix).tolist()}
                             if self.matrix is not None else None),
        }
        self._project = data
        with open(self.project_file(), "w", encoding="utf-8") as f:
            json.dump(data, f, indent=1)

    def gates_changed(self, e: Entry, everyone: bool = True) -> None:
        """`e`'s tree was edited. With `everyone`, the other files get a
        copy — FlowJo's group, where one gate is every sample's gate.
        Counts are forgotten and the trees saved beside the files."""
        changed = [e]
        if everyone:
            for other in self.usable():
                if other is not e:
                    other.gates = copy_tree(e.gates, keep_counts=False)
                    changed.append(other)
        for node in e.gates.walk():
            node.reference_count = None   # no longer what another program counted
        for c in changed:
            c.edited, c.borrowed = True, False
            if c.sample is not None:
                c.sample.gates, c.sample._masks = c.gates, None
        self._counts.clear()
        self._keep_population()
        self.save_project()

    def discard_gates(self) -> None:
        """Back to the workspace's gates: the saved ones are deleted. The
        columns and what was typed in them stay."""
        try:
            os.remove(os.path.join(self.folder, OLD_GATES_FILE))
        except OSError:
            pass
        self._project = {**getattr(self, "_project", {}), "gates": {}}
        self._borrow()
        self._keep_population()
        if os.path.exists(self.project_file()):
            self.save_project()

    # --- annotations -------------------------------------------------------------

    def annotate(self, pattern: str) -> list[str]:
        """Fill columns from the file names: every named group of the
        pattern, `(?P<Day>\\d+)`, is a column, added if it is new. A file
        the pattern does not match keeps what it had. Returns the columns
        it filled."""
        names = list(re.compile(pattern).groupindex)
        for name in names:
            if name not in self.columns:
                self.columns.append(name)
        for e in self.usable():
            found = annotate(e.name, pattern)
            for name in names:
                if name in found:
                    e.meta[name] = found[name]
        return names

    def add_column(self, name: str) -> str:
        name = name.strip()
        if not name or name in self.columns or name in RESERVED:
            raise ValueError(f"there is already a column {name!r}" if name
                             else "a column needs a name")
        self.columns.append(name)
        for e in self.entries:
            e.meta.setdefault(name, "")
        return name

    def remove_column(self, name: str) -> None:
        self.columns.remove(name)
        self.group_by = [c for c in self.group_by if c != name]
        for e in self.entries:
            e.meta.pop(name, None)

    def rename_column(self, old: str, new: str) -> None:
        new = new.strip()
        if not new or new in RESERVED or (new in self.columns and new != old):
            raise ValueError(f"there is already a column {new!r}" if new
                             else "a column needs a name")
        self.columns[self.columns.index(old)] = new
        self.group_by = [new if c == old else c for c in self.group_by]
        for e in self.entries:
            e.meta[new] = e.meta.pop(old, "")

    def set_value(self, entries: list[Entry], column: str, value: str) -> None:
        for e in entries:
            e.meta[column] = value.strip()

    def repeats(self) -> list[list[str]]:
        """Included files whose names differ only by a mark at the end —
        `14Case4` and ``14Case4` `` — which is how a re-acquisition of one
        tube is usually named: one animal that would be counted as two."""
        seen: dict[str, list[str]] = {}
        for e in self.entries:
            if e.include and not e.problem:
                seen.setdefault(re.sub(r"[^A-Za-z0-9]+$", "", e.name), []).append(e.name)
        return [names for names in seen.values() if len(names) > 1]

    def meta_keys(self) -> list[str]:
        return list(self.columns)

    def usable(self) -> list[Entry]:
        return [e for e in self.entries if not e.problem]

    def populations(self) -> list[str]:
        """Gate paths every gated file has, shallowest first; without the
        ungated root."""
        gated = [e for e in self.usable() if e.gated]
        if not gated:
            return []
        sets = [[n.path for n in e.gates.walk() if n.gate is not None]
                for e in gated]
        return [p for p in sets[0] if all(p in s for s in sets[1:])]

    def available(self, e: Entry) -> bool:
        """Can this file give events of the chosen population?"""
        return not e.problem and (self.population == ALL_EVENTS or e.gated)

    def group_of(self, meta: dict[str, str]) -> str:
        return " · ".join(meta.get(k, "") or "—" for k in self.group_by) if self.group_by else ""

    # --- events ----------------------------------------------------------

    def load(self, e: Entry) -> Sample:
        if e.sample is None:
            matrix = self.matrix_of(e)
            e.sample = (Sample.open(e.path, compensation=matrix, gates=e.gates)
                        if matrix is not None else Sample.open(e.path, gates=e.gates))
            e.sample.meta = {c: e.meta.get(c, "") for c in self.columns}
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
                s.meta = {c: e.meta.get(c, "") for c in self.columns}
                samples.append(s)
        self.dataset = concatenate(samples, p, self.population, self.seed)
        self.normalization_check = None
        if self.normalization:
            self._normalise(self.dataset, {s.name: s for s in samples})
        self.flowsom = None
        # the events are in the dataset now; the files need not stay in memory
        for e in self.entries:
            e.sample = None
        return self.dataset

    # --- taking a population elsewhere ---------------------------------------

    def export_population(self, entries: list[Entry], population: str, folder: str,
                          fmt: str = "fcs", compensated: bool = True,
                          progress=None) -> list[str]:
        """One file per sample holding only the events of `population`.

        Compensated, the fluorescence columns carry the values this
        program works with and are named `Comp-…`, so nothing downstream
        compensates them again. Not compensated, the events are as
        acquired and the FCS file carries the matrix in use as its
        `$SPILLOVER`, for whatever reads it next to apply."""
        os.makedirs(folder, exist_ok=True)
        tail = re.sub(r"[^A-Za-z0-9+._ -]+", "_", population.split("/")[-1]).strip() or "all"
        written = []
        for i, e in enumerate(entries):
            sample = self.load(e)
            masks = sample.masks()
            if population not in masks:
                continue
            rows = masks[population]
            f = sample.fcs
            touched = set(sample.compensation.channels) if sample.compensation else set()
            data = (sample.compensated if compensated else f.events)[rows]
            names = [("Comp-" + c) if compensated and c in touched else c for c in f.channels]
            path = os.path.join(folder, f"{e.name}_{tail}.{fmt}")
            if fmt == "csv":
                heads = [m or c for c, m in zip(names, f.markers)]
                with open(path, "w", encoding=cluster_stats.CSV_ENCODING) as out:
                    out.write(",".join('"' + h.replace('"', '""') + '"' if "," in h else h
                                       for h in heads) + "\n")
                    np.savetxt(out, data, delimiter=",", fmt="%.6g")
            else:
                keep = {k: f.keywords[k] for k in ("$CYT", "$DATE", "$BTIM", "$ETIM", "TUBE NAME")
                        if f.keywords.get(k)}
                keep["MILC POPULATION"] = population or "All events"
                keep["MILC SOURCE"] = f.name
                if not compensated and sample.compensation is not None:
                    m = sample.compensation
                    keep["$SPILLOVER"] = ",".join(
                        [str(len(m.channels)), *m.channels,
                         *(f"{v:.10g}" for v in np.asarray(m.matrix).ravel())])
                write_fcs(path, data, names, list(f.markers), list(f.ranges), keep)
            written.append(path)
            if progress:
                progress(i + 1, len(entries), e.name)
        return written

    # --- statistics ------------------------------------------------------

    def stat_channels(self) -> list[tuple[str, str]]:
        """(label, column) for every parameter a statistic can be asked
        of: the stain's name where there is one, and `Comp-…` as the column
        where the parameter is compensated — the value FlowJo reports."""
        entries = self.usable()
        if not entries:
            return []
        first = entries[0]
        matrix = self.matrix_of(first)
        spill = None
        if matrix is None:
            from . import compensation
            spill = compensation.from_keywords(read_header(first.path))
        compensated = set((matrix or spill).channels) if (matrix or spill) else set()
        return [(m or c, ("Comp-" + c) if c in compensated else c)
                for c, m in zip(first.channels, first.markers)]

    def statistics(self, populations: list[str] | None = None,
                   channels: list[tuple[str, str]] | None = None,
                   statistics: list[str] | None = None, progress=None,
                   only_included: bool = True) -> stats.StatsTable:
        """The table: a row for each population of each sample, with the
        annotation columns beside the name so it can be sorted, filtered
        and compared by group wherever it is taken."""
        chosen = [e for e in self.usable() if e.include or not only_included]
        table = stats.StatsTable(stats.headers(self.columns, channels, statistics))
        for i, e in enumerate(chosen):
            sample = self.load(e)
            notes = [e.meta.get(c, "") for c in self.columns]
            for row in stats.population_rows(sample, populations, channels, statistics):
                table.rows.append([e.name, *notes, row[1], row[0], *row[2:]])
            if progress:
                progress(i + 1, len(chosen), e.name)
        return table

    # --- batch normalisation -----------------------------------------------------

    def set_normalization(self, batch: str | None, reference: tuple[str, str] | None = None,
                          clusters: int = 1) -> None:
        """Normalise at the next concatenation: `batch` is the column that
        says which batch a sample was acquired in, `reference` the
        (column, value) of the samples that should look alike in every
        batch. None for `batch` turns it off."""
        self.normalization = ({"batch": batch, "reference": list(reference),
                               "clusters": int(clusters)}
                              if batch and reference and reference[0] else None)
        self.save_project()

    def reference_batches(self) -> dict[str, list[Entry]]:
        """Batch → its reference samples, as the settings stand: included,
        with the population, and carrying the reference value."""
        cfg = self.normalization
        if not cfg:
            return {}
        column, value = cfg["reference"]
        out: dict[str, list[Entry]] = {}
        for e in self.entries:
            if e.include and self.available(e) and e.meta.get(column, "") == value:
                out.setdefault(e.meta.get(cfg["batch"], ""), []).append(e)
        return out

    def _normalise(self, ds: Dataset, samples: dict) -> None:
        """Bring every batch of `ds` onto the mean of the batches, as
        learnt from the reference samples — from all the events of their
        population, not only the ones drawn — on each marker's own scale.
        What was done, and how far apart the batches' references were
        before and after, is written into the dataset."""
        from . import normalize
        cfg = self.normalization
        refs = self.reference_batches()
        if len(refs) < 2:
            raise ValueError(
                "batch normalisation needs reference samples in at least two batches; "
                f"found them in {len(refs)}")
        markers = [c for c, comp in zip(ds.channels, ds.compensated) if comp]
        cols = [ds.channels.index(c) for c in markers]
        scales = [self.scale_of(c) if c in self.scales else Biex() for c in markers]

        def forward(values):
            return np.column_stack([s(values[:, k]) if s is not None else values[:, k]
                                    for k, s in enumerate(scales)])

        def back(values):
            return np.column_stack([s.inverse(values[:, k]) if s is not None else values[:, k]
                                    for k, s in enumerate(scales)])

        reference = {}
        for batch, entries in refs.items():
            parts = []
            for e in entries:
                sample = samples.get(e.name) or self.load(e)
                parts.append(sample.compensated[sample.mask(self.population)][:, cols])
            reference[batch] = forward(np.vstack(parts).astype(np.float64))
        model = normalize.fit(reference, n_clusters=cfg.get("clusters", 1), seed=self.seed)
        batch_of = [m.get(cfg["batch"], "") for m in ds.meta]
        skipped = []
        for i, batch in enumerate(batch_of):
            rows = ds.sample_index == i
            if batch not in model.batches:
                skipped.append(ds.samples[i])   # a batch with no reference stays as it is
                continue
            fixed = normalize.apply(model, forward(ds.X[rows][:, cols].astype(np.float64)), batch)
            ds.X[np.ix_(np.flatnonzero(rows), cols)] = back(fixed).astype(ds.X.dtype)
        # How far apart the references are afterwards is measured on events
        # the curves were not fitted on: half of each batch's is set aside,
        # a second model is fitted on the other half and judged on these.
        # Measured on the events it was fitted on, the answer is zero by
        # construction and says nothing.
        rng = np.random.default_rng(self.seed)
        held, rest = {}, {}
        for b, x in reference.items():
            order = rng.permutation(len(x))
            held[b], rest[b] = x[order[: len(x) // 2]], x[order[len(x) // 2:]]
        judge = normalize.fit(rest, n_clusters=cfg.get("clusters", 1), seed=self.seed)
        after = {b: normalize.apply(judge, x, b) for b, x in held.items()}
        self.normalization_check = {
            "markers": [ds.label(c) for c in markers], "scales": scales,
            "before": held, "after": after}
        ds.provenance["normalization"] = {
            **cfg, "markers": [ds.label(c) for c in markers], "batches": model.batches,
            "reference_samples": {b: [e.name for e in es] for b, es in refs.items()},
            "reference_events": {b: int(len(x)) for b, x in reference.items()},
            "clusters_used": model.n_clusters,
            "distance_before": normalize.distance(held).tolist(),
            "distance_after": normalize.distance(after).tolist(),
            "fallbacks": [list(f) for f in model.fallbacks], "not_normalised": skipped,
        }

    # --- analysis --------------------------------------------------------

    def default_channels(self) -> list[str]:
        """The stained fluorescence channels: what a first analysis is
        usually run on. The user unticks lineage and viability dyes."""
        ds = self.dataset
        return [c for c, m, comp in zip(ds.channels, ds.markers, ds.compensated)
                if comp and m]

    def features(self, channels: list[str], transform: str = "biex",
                 cofactor: float = 150.0) -> np.ndarray:
        """What the algorithms are fed. `biex` is each marker on the scale
        its axis has — FlowJo's biexponential unless another was chosen —
        on 0 to 1023, which is what FlowJo hands its plugins."""
        if transform != "biex":
            return self.dataset.matrix(channels, Arcsinh(cofactor))
        columns = []
        for c in channels:
            x = self.dataset.column(c).astype(np.float64)
            scale = self.scale_of(c) if c in self.scales else Biex()
            if scale is None:
                columns.append(x / 262144.0 * 1023.0)
            elif isinstance(scale, Arcsinh):
                columns.append(scale(x))
            else:
                columns.append(scale(x) * 1023.0)
        return np.column_stack(columns)

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
