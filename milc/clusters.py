"""
Clusters as a person works with them: named, annotated, merged, with a
hand-drawn one beside them.

An algorithm leaves a number per event in `ds.derived[key]` and that is
never written over. What the user makes of it is kept apart, as edits on
top — which clusters were merged into which, what each is called, a note —
so a merge can always be undone and nothing the algorithm said is lost.
`view` puts the two together and is the one place the rest of the program
asks what the clusters now are.

The edits live in `ds.runs[key]["edits"]`: they are saved in the dataset's
header with the run they belong to, and a new run of the algorithm, which
replaces that entry, starts clean — names given to the old clusters would
be wrong on the new ones. Their keys are text because JSON's are.
A cluster drawn by hand is per event and so sits in `ds.derived`, under
`key + "_manual"`: 0, or the cluster's number.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cluster_stats import Table

#: the clustering asked about where none is named
FLOWSOM = "FlowSOM"
MANUAL = "_manual"

#: how far apart the two groups of clusters must be for a marker to be
#: called positive in one and negative in the other, by the scale the
#: table is in: a tenth of the 0-1023 the biexponential scale is fed as
#: (about half a decade where that scale is logarithmic), and 1 on an
#: arcsinh scale (a factor of e, 0.43 of a decade, above the cofactor)
MIN_GAP = {"biex": 100.0, "logicle": 100.0, "arcsinh": 1.0}


@dataclass
class Clustering:
    """The clusters as the user now has them."""
    key: str
    #: the clusters that hold events, in the order they are shown: the
    #: algorithm's (a merged one under its lowest number), then the manual
    ids: list[int]
    #: what each is called, in the same order; no two alike
    names: list[str]
    notes: list[str]
    #: per event: its cluster, one of `ids`
    labels: np.ndarray
    #: cluster → the algorithm's clusters it is made of; none for a manual one
    members: dict[int, list[int]]
    #: the colour the user gave each, in the same order; "" where none was
    #: given and the palette's is used
    colors: list[str] = field(default_factory=list)

    def palette(self, defaults: list[str]) -> list[str]:
        """A colour for each cluster: the user's where one was given, else
        the one of `defaults` in its place."""
        given = self.colors or [""] * len(self.ids)
        return [c or d for c, d in zip(given, defaults)]

    def color_of(self, defaults: list[str]) -> dict[int, str]:
        return dict(zip(self.ids, self.palette(defaults)))

    def name(self, cluster: int) -> str:
        return self.names[self.ids.index(cluster)]

    def note(self, cluster: int) -> str:
        return self.notes[self.ids.index(cluster)]

    def codes(self) -> np.ndarray:
        """Per event: the position of its cluster in `ids`."""
        return np.searchsorted(np.asarray(self.ids), self.labels)

    def is_manual(self, cluster: int) -> bool:
        return not self.members[cluster]

    def merged(self) -> bool:
        return any(len(m) > 1 for m in self.members.values())


def _edits(ds, key: str, write: bool = False) -> dict:
    if not write:
        found = ds.runs.get(key, {}).get("edits", {})
        return {k: found.get(k, {} if k != "manual" else [])
                for k in ("names", "notes", "merged", "manual", "colors")}
    edits = ds.runs.setdefault(key, {}).setdefault(
        "edits", {"names": {}, "notes": {}, "merged": {}, "manual": []})
    edits.setdefault("colors", {})        # runs kept before clusters had colours
    return edits


def merged_labels(ds, key: str = FLOWSOM) -> np.ndarray:
    """Per event: the algorithm's cluster after the merges, without the
    manual clusters — what a node of the map belongs to."""
    raw = ds.derived[key]
    lookup = np.arange(int(raw.max()) + 1)
    for old, new in _edits(ds, key)["merged"].items():
        if int(old) < lookup.size:
            lookup[int(old)] = int(new)
    return lookup[raw]


def view(ds, key: str = FLOWSOM) -> Clustering | None:
    """The clustering `key` of the dataset with the user's edits applied,
    or None where there is no such clustering."""
    if ds is None or key not in ds.derived:
        return None
    e = _edits(ds, key)
    raw = ds.derived[key]
    labels = merged_labels(ds, key)
    members: dict[int, list[int]] = {}
    for c in np.unique(raw):
        target = int(e["merged"].get(str(int(c)), c))
        members.setdefault(target, []).append(int(c))
    drawn = ds.derived.get(key + MANUAL)
    if drawn is not None and e["manual"]:
        inside = np.isin(drawn, e["manual"])
        labels = np.where(inside, drawn, labels)
    ids = [int(c) for c in np.unique(labels)]
    members = {c: members.get(c, []) if c not in e["manual"] else [] for c in ids}
    names = [e["names"].get(str(c)) or "+".join(str(m) for m in members[c]) or str(c)
             for c in ids]
    return Clustering(key, ids, names, [e["notes"].get(str(c), "") for c in ids],
                      labels, members, [e["colors"].get(str(c), "") for c in ids])


def named(table: Table, clustering: Clustering) -> Table:
    """A table whose rows are cluster numbers, with their names instead
    and the numbers beside them."""
    ids = [int(r) for r in table.rows]
    return Table([clustering.name(c) for c in ids], table.columns, table.values, ids)


# --- edits -------------------------------------------------------------------

def rename(ds, names: dict[int, str], key: str = FLOWSOM) -> None:
    """Call these clusters this. An empty name gives a cluster back the
    one it had; two clusters cannot share a name, because a name is what
    a row of every table is known by."""
    now = view(ds, key)
    after = dict(zip(now.ids, now.names))
    for cluster, name in names.items():
        default = "+".join(str(m) for m in now.members[cluster]) or str(cluster)
        after[cluster] = name.strip() or default
    taken = [n for n in set(after.values()) if list(after.values()).count(n) > 1]
    if taken:
        raise ValueError(f"there is already a cluster called {taken[0]!r}")
    kept = _edits(ds, key, write=True)["names"]
    for cluster, name in names.items():
        if name.strip():
            kept[str(cluster)] = name.strip()
        else:
            kept.pop(str(cluster), None)


def set_color(ds, cluster: int, color: str, key: str = FLOWSOM) -> None:
    """Draw this cluster in this colour (`#rrggbb`) wherever clusters are
    coloured. An empty colour gives it back the palette's."""
    colors = _edits(ds, key, write=True)["colors"]
    if color:
        colors[str(int(cluster))] = str(color)
    else:
        colors.pop(str(int(cluster)), None)


def set_note(ds, cluster: int, note: str, key: str = FLOWSOM) -> None:
    notes = _edits(ds, key, write=True)["notes"]
    if note.strip():
        notes[str(cluster)] = note.strip()
    else:
        notes.pop(str(cluster), None)


def merge(ds, clusters: list[int], key: str = FLOWSOM) -> int:
    """Make one cluster of these — its events are theirs together — and
    return its number, the lowest of them. It keeps the name and the note
    of that one; the others' are dropped with them."""
    now = view(ds, key)
    chosen = sorted(set(int(c) for c in clusters))
    if len(chosen) < 2:
        raise ValueError("choose two or more clusters to merge")
    if any(c not in now.ids or now.is_manual(c) for c in chosen):
        raise ValueError("only the algorithm's clusters can be merged")
    e = _edits(ds, key, write=True)
    target = chosen[0]
    for c in chosen[1:]:
        for m in now.members[c]:
            e["merged"][str(m)] = target
        e["names"].pop(str(c), None)
        e["notes"].pop(str(c), None)
    return target


def undo_merges(ds, key: str = FLOWSOM) -> None:
    """Back to the clusters the algorithm gave. What a merged cluster was
    called described the whole of it, so that name and note go too."""
    now = view(ds, key)
    e = _edits(ds, key, write=True)
    for cluster, members in now.members.items():
        if len(members) > 1:
            e["names"].pop(str(cluster), None)
            e["notes"].pop(str(cluster), None)
    e["merged"].clear()


def add_manual(ds, mask: np.ndarray, name: str, key: str = FLOWSOM) -> int:
    """These events as a cluster of their own, called `name`; every other
    event stays in its cluster. Returns the new cluster's number."""
    now = view(ds, key)
    name = name.strip()
    if not name or name in now.names:
        raise ValueError(f"there is already a cluster called {name!r}" if name
                         else "a cluster drawn by hand needs a name")
    if not mask.any():
        raise ValueError("no events are selected")
    e = _edits(ds, key, write=True)
    drawn = ds.derived.get(key + MANUAL)
    if drawn is None or not e["manual"]:   # what an earlier run left is not this run's
        drawn = np.zeros(len(ds), dtype=np.int32)
    cluster = max(int(ds.derived[key].max()), int(drawn.max())) + 1
    drawn = drawn.copy()
    drawn[mask] = cluster
    ds.derived[key + MANUAL] = drawn
    # one drawn over entirely by this one is gone
    e["manual"][:] = [c for c in e["manual"] if (drawn == c).any()] + [cluster]
    e["names"][str(cluster)] = name
    return cluster


def remove_manual(ds, cluster: int, key: str = FLOWSOM) -> None:
    """Its events go back to the clusters they were in."""
    e = _edits(ds, key, write=True)
    if cluster not in e["manual"]:
        raise ValueError("only a cluster drawn by hand can be removed")
    e["manual"].remove(cluster)
    e["names"].pop(str(cluster), None)
    e["notes"].pop(str(cluster), None)
    drawn = ds.derived[key + MANUAL].copy()
    drawn[drawn == cluster] = 0
    ds.derived[key + MANUAL] = drawn


# --- names from what the clusters express --------------------------------------

def calls(expression: Table, min_gap: float) -> np.ndarray:
    """Cluster x marker: +1 where the cluster is positive for the marker,
    -1 where it is negative, 0 where the marker says nothing.

    The rule, per marker: sort the clusters' values and cut them in two at
    the widest gap between neighbours. The clusters above the cut are
    positive, those below negative — unless that gap is narrower than
    `min_gap` (in the table's own units), in which case the marker does
    not tell the clusters apart and is left out for all of them."""
    v = np.asarray(expression.values, dtype=np.float64)
    out = np.zeros(v.shape, dtype=np.int8)
    if v.shape[0] < 2:
        return out
    for j in range(v.shape[1]):
        ordered = np.sort(v[:, j])
        gaps = np.diff(ordered)
        k = int(np.argmax(gaps))
        if gaps[k] >= min_gap:
            out[:, j] = np.where(v[:, j] > (ordered[k] + ordered[k + 1]) / 2, 1, -1)
    return out


def suggest_names(expression: Table, min_gap: float) -> list[str]:
    """A name for each row of the table in the form `CD4+ Foxp3+`: the
    markers the cluster is positive for by `calls`, in the table's order;
    for a cluster positive for none, the ones it is negative for
    (`CD4- CD8-`). Empty where no marker separates the clusters. Two
    clusters the markers do not tell apart get the same name with
    ` (2)`, ` (3)` after it."""
    sign = calls(expression, min_gap)
    short = [c.split()[0] if c else c for c in expression.columns]
    names, seen = [], {}
    for row in sign:
        marks = ([m + "+" for m, s in zip(short, row) if s > 0]
                 or [m + "-" for m, s in zip(short, row) if s < 0])
        name = " ".join(marks)
        seen[name] = seen.get(name, 0) + 1
        names.append(f"{name} ({seen[name]})" if name and seen[name] > 1 else name)
    return names
