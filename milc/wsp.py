"""
A FlowJo 10 workspace (.wsp), read for what the analysis needs from it:
which files, the compensation matrix each was given, the axis transforms,
and the gating tree with the counts FlowJo reported.

Nothing is written back. Plugin nodes (UMAP, FlowSOM, ClusterExplorer) are
skipped: their results are recomputed here, not imported.
"""

from __future__ import annotations

import os
import re
import urllib.parse
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

import numpy as np

from .compensation import Compensation
from .gating import GateNode, PolygonGate, RectangleGate
from .transforms import Biex, Linear


def _tag(e) -> str:
    return e.tag.rsplit("}", 1)[-1]


def _attrs(e) -> dict[str, str]:
    return {k.rsplit("}", 1)[-1]: v for k, v in e.attrib.items()}


def _children(e, name: str):
    return [c for c in e if _tag(c) == name]


def _child(e, name: str):
    found = _children(e, name)
    return found[0] if found else None


@dataclass
class WspSample:
    name: str
    #: where FlowJo found the file — usually a path on another machine
    uri: str
    sample_id: str
    compensation: Compensation | None
    transforms: dict[str, object]
    gates: GateNode
    keywords: dict[str, str] = field(default_factory=dict)

    @property
    def filename(self) -> str:
        path = urllib.parse.unquote(urllib.parse.urlparse(self.uri).path)
        return re.split(r"[\\/]", path)[-1] or self.name

    def locate(self, folder: str) -> str | None:
        path = os.path.join(folder, self.filename)
        return path if os.path.exists(path) else None


@dataclass
class Workspace:
    path: str
    flowjo_version: str
    samples: list[WspSample]
    #: the gates FlowJo keeps on a group, by the group's name: what it
    #: applies to every sample put in that group
    groups: dict[str, GateNode] = field(default_factory=dict)

    def gated(self) -> list[WspSample]:
        return [s for s in self.samples
                if s.compensation is not None and s.gates.children]

    def shared(self) -> tuple[GateNode, Compensation] | None:
        """The gates and the matrix every gated sample of the workspace
        has in common — what a file the workspace does not list would get
        if it were dropped into the group in FlowJo. None unless they
        really are in common: a sample gated or compensated its own way
        means there is no one answer to give the others."""
        gated = self.gated()
        if not gated:
            return None
        first = gated[0]
        for s in gated[1:]:
            if (s.compensation.channels != first.compensation.channels
                    or not np.array_equal(s.compensation.matrix, first.compensation.matrix)
                    or signature(s.gates) != signature(first.gates)):
                return None
        for tree in self.groups.values():
            if signature(tree) == signature(first.gates):
                return tree, first.compensation
        return None

    def sample(self, name: str) -> WspSample:
        for s in self.samples:
            if name in (s.name, s.filename):
                return s
        raise KeyError(name)


def read_wsp(path: str) -> Workspace:
    root = ET.parse(path).getroot()
    samples = []
    sample_list = _child(root, "SampleList")
    for s in (_children(sample_list, "Sample") if sample_list is not None else []):
        dataset = _child(s, "DataSet")
        node = _child(s, "SampleNode")
        if dataset is None or node is None:
            continue
        kw_el = _child(s, "Keywords")
        keywords = {k.get("name"): k.get("value", "") for k in kw_el} if kw_el is not None else {}
        tree = GateNode(name=node.get("name", ""),
                        reference_count=_int(node.get("count")))
        _read_populations(node, tree)
        samples.append(WspSample(
            name=node.get("name", ""),
            uri=dataset.get("uri", ""),
            sample_id=dataset.get("sampleID", ""),
            compensation=_read_matrix(_child(s, "spilloverMatrix")),
            transforms=_read_transforms(_child(s, "Transformations")),
            gates=tree,
            keywords=keywords,
        ))
    groups = {}
    groups_el = _child(root, "Groups")
    for g in (_children(groups_el, "GroupNode") if groups_el is not None else []):
        tree = GateNode(name=g.get("name", ""))
        _read_populations(g, tree)
        if tree.children:
            # the counts stored on a group's gates are one sample's, and
            # would be read as FlowJo's count for whichever file gets them
            for node in tree.walk():
                node.reference_count = None
            groups[tree.name] = tree
    return Workspace(path=path, flowjo_version=root.get("flowJoVersion", ""),
                     samples=samples, groups=groups)


def signature(tree: GateNode) -> list:
    """A gating tree as plain values, to tell whether two are the same."""
    out = []
    for node in tree.walk():
        g = node.gate
        if g is None:
            continue
        shape = (np.asarray(g.vertices).round(6).tolist() if isinstance(g, PolygonGate)
                 else [g.mins, g.maxs])
        out.append((node.path, type(g).__name__, list(g.dims), g.inside, shape))
    return out


def _int(text) -> int | None:
    try:
        v = int(text)
        return v if v >= 0 else None
    except (TypeError, ValueError):
        return None


def _float(text) -> float | None:
    return None if text in (None, "") else float(text)


def _read_matrix(el) -> Compensation | None:
    if el is None:
        return None
    params = _child(el, "parameters")
    names = [_attrs(p)["name"] for p in params]
    matrix = np.zeros((len(names), len(names)))
    for row in _children(el, "spillover"):
        i = names.index(_attrs(row)["parameter"])
        for c in _children(row, "coefficient"):
            a = _attrs(c)
            matrix[i, names.index(a["parameter"])] = float(a["value"])
    return Compensation(names, matrix, name=el.get("name", ""))


def _read_transforms(el) -> dict[str, object]:
    out: dict[str, object] = {}
    for t in (el if el is not None else []):
        a = _attrs(t)
        param = _child(t, "parameter")
        if param is None:
            continue
        name = _attrs(param)["name"]
        kind = _tag(t)
        if kind == "linear":
            out[name] = Linear(float(a.get("minRange", 0)), float(a.get("maxRange", 262144)))
        elif kind == "biex":
            out[name] = Biex(int(float(a.get("length", 256))),
                             float(a.get("maxRange", 262144)),
                             float(a.get("neg", 0)), float(a.get("width", -100)),
                             float(a.get("pos", 4.5)))
        # log, logicle, hyperlog, fasinh: not met yet, left out rather
        # than guessed
    return out


def _read_gate(el):
    shape = next(iter(el), None)
    if shape is None:
        return None
    a = _attrs(shape)
    inside = a.get("eventsInside", "1") != "0"
    dims, mins, maxs = [], [], []
    for d in _children(shape, "dimension"):
        da = _attrs(d)
        dims.append(_attrs(next(iter(d)))["name"])
        mins.append(_float(da.get("min")))
        maxs.append(_float(da.get("max")))
    kind = _tag(shape)
    if kind == "RectangleGate":
        return RectangleGate(dims, mins, maxs, inside)
    if kind == "PolygonGate":
        verts = [[float(_attrs(c)["value"]) for c in _children(v, "coordinate")]
                 for v in _children(shape, "vertex")]
        return PolygonGate(dims, np.array(verts), inside)
    return None  # ellipses, curly quadrants, booleans: not met yet


def _read_populations(el, parent: GateNode) -> None:
    for sub in _children(el, "Subpopulations"):
        for p in _children(sub, "Population"):
            gate_el = _child(p, "Gate")
            gate = _read_gate(gate_el) if gate_el is not None else None
            if gate is None:
                continue
            node = parent.add(GateNode(name=p.get("name", "").strip(), gate=gate,
                                       reference_count=_int(p.get("count"))))
            _read_populations(p, node)
