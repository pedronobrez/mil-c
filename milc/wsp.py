"""
A FlowJo 10 workspace (.wsp), read for what the analysis needs from it:
which files, the compensation matrix each was given, the axis transforms,
and the gating tree with the counts FlowJo reported.

Writing one is `wsp_write.py`. Plugin nodes (UMAP, FlowSOM, ClusterExplorer) are
skipped: their results are recomputed here, not imported.

Rectangles and polygons are read by rules checked against FlowJo's own
counts on a real experiment (tests/test_real.py). Quadrants, combinations
(and, or, not), ellipses and curly quadrants are not: no workspace of ours
holds one. They
are read the way two other readers of the format do — FlowKit (Python)
and CytoML with cytolib (Bioconductor) — where the two agree, and the
workspaces FlowKit publishes with its tests were looked at to see what
FlowJo 10.6 to 10.8 really writes. So three things hold:

  - a population read that way says so in its `note`;
  - what could not be read with confidence is not read, and is listed in
    `Workspace.skipped` with the reason;
  - FlowJo's count is beside every population, and `count_mismatches`
    lists those whose count here is not FlowJo's. That is the check.
"""

from __future__ import annotations

import os
import re
import urllib.parse
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

import numpy as np

from .compensation import Compensation
from .gating import (AND, CURL_TOP, NOT, OR, QUADRANTS, BooleanGate, CurlyQuadrantGate,
                     EllipseGate, GateNode, PolygonGate, QuadrantGate, RectangleGate,
                     tree_to_dict)
from .transforms import Biex, Linear, Logicle


def _tag(e) -> str:
    return e.tag.rsplit("}", 1)[-1]


def _attrs(e) -> dict[str, str]:
    return {k.rsplit("}", 1)[-1]: v for k, v in e.attrib.items()}


def _children(e, name: str):
    return [c for c in e if _tag(c) == name]


def _child(e, name: str):
    found = _children(e, name)
    return found[0] if found else None


#: begins the `note` of every population read by a rule taken from other
#: readers of the format, not checked here against FlowJo's counts
UNVERIFIED = "read from the workspace by a rule not yet checked against FlowJo"

#: how far a count may be from FlowJo's before it is worth a look: what was
#: measured on 72 rectangles and polygons, with its margin. FlowJo counts on
#: a 256-bin grid by a rule not reproduced, so equality is not expected.
COUNT_TOLERANCE = 0.012
COUNT_SLACK = 3


@dataclass
class Skipped:
    """A population of the workspace that was not read."""
    #: the sample or the group whose tree it is in
    owner: str
    path: str
    #: the workspace's own word for it: `EllipsoidGate`, `OrNode`…
    kind: str
    reason: str
    #: how many populations under it went with it
    below: int = 0


@dataclass
class Mismatch:
    """A population counted differently here and in FlowJo."""
    path: str
    count: int
    flowjo: int
    #: the population's `note`: whether its gate was read by an unchecked rule
    note: str = ""


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
    #: every population that was left unread, and why
    skipped: list[Skipped] = field(default_factory=list)

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
    skipped: list[Skipped] = []
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
        transforms = _read_transforms(_child(s, "Transformations"))
        seconds = _time_unit(keywords)
        if seconds and isinstance(transforms.get(seconds[0]), Linear):
            t = transforms[seconds[0]]
            transforms[seconds[0]] = Linear(t.min / seconds[1], t.max / seconds[1])
        _read_populations(node, tree, transforms, skipped, tree.name)
        if seconds:
            _time_as_stored(tree, *seconds)
        samples.append(WspSample(
            name=node.get("name", ""),
            uri=dataset.get("uri", ""),
            sample_id=dataset.get("sampleID", ""),
            compensation=_read_matrix(_child(s, "spilloverMatrix")),
            transforms=transforms,
            gates=tree,
            keywords=keywords,
        ))
    groups = {}
    groups_el = _child(root, "Groups")
    for g in (_children(groups_el, "GroupNode") if groups_el is not None else []):
        tree = GateNode(name=g.get("name", ""))
        _read_populations(g, tree, None, skipped, tree.name)
        if tree.children:
            # the counts stored on a group's gates are one sample's, and
            # would be read as FlowJo's count for whichever file gets them
            for node in tree.walk():
                node.reference_count = None
            groups[tree.name] = tree
    return Workspace(path=path, flowjo_version=root.get("flowJoVersion", ""),
                     samples=samples, groups=groups, skipped=skipped)


def _time_unit(keywords: dict[str, str]) -> tuple[str, float] | None:
    """(the Time parameter, seconds in one of its units), where the file
    says. FlowJo shows Time in seconds — the stored values times
    `$TIMESTEP`: the range it writes for Time in a workspace is that of
    the file times the step — so that is what a gate on Time is in."""
    try:
        step = float(keywords.get("$TIMESTEP") or 0)
    except ValueError:
        return None
    name = next((v.strip() for k, v in keywords.items()
                 if k and re.fullmatch(r"\$P\d+N", k) and v.strip().lower() == "time"), None)
    return (name, step) if name and step > 0 and step != 1 else None


def _time_as_stored(tree: GateNode, name: str, step: float) -> None:
    """Gates on Time, from FlowJo's seconds to the units of the file."""
    for node in tree.walk():
        g = node.gate
        if g is None or name not in g.dims:
            continue
        k = g.dims.index(name)
        if isinstance(g, RectangleGate):
            g.mins[k] = None if g.mins[k] is None else g.mins[k] / step
            g.maxs[k] = None if g.maxs[k] is None else g.maxs[k] / step
        elif isinstance(g, PolygonGate):
            g.vertices[:, k] /= step
        elif isinstance(g, QuadrantGate):
            g.x, g.y = (g.x / step, g.y) if k == 0 else (g.x, g.y / step)
        elif isinstance(g, CurlyQuadrantGate):
            g.center[k] /= step


def signature(tree: GateNode) -> list:
    """A gating tree as plain values, to tell whether two are the same."""
    def plain(d: dict) -> dict:
        gate = dict(d.get("gate") or {})
        if "vertices" in gate:
            gate["vertices"] = np.round(gate["vertices"], 6).tolist()
        return {"name": d["name"], "gate": gate,
                "children": [plain(c) for c in d["children"]]}
    return [plain(tree_to_dict(c)) for c in tree.children]


def unverified(tree: GateNode) -> list[GateNode]:
    """The populations of a tree whose gates were read by a rule not
    checked against FlowJo: the ones to warn about."""
    return [node for node in tree.walk() if node.note.startswith(UNVERIFIED)]


def count_mismatches(tree: GateNode, count, tolerance: float = COUNT_TOLERANCE,
                     slack: int = COUNT_SLACK) -> list[Mismatch]:
    """The populations whose count here is not the one FlowJo wrote in the
    workspace, by more than the tolerance. `count` gives ours for a path
    (`Sample.count`). For the gates no FlowJo file was at hand to test,
    this is the only check there is."""
    out = []
    for node in tree.walk():
        ref = node.reference_count
        if node.gate is None or ref is None:
            continue
        ours = int(count(node.path))
        if abs(ours - ref) > max(tolerance * ref, slack):
            out.append(Mismatch(node.path, ours, ref, node.note))
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
        elif kind == "logicle":
            # T, W, M and A under their own names, as FlowKit reads them
            # (wsp_utils.py, _parse_wsp_transforms); a scale that cannot
            # be is left out like one that is not known
            try:
                out[name] = Logicle(float(a.get("T", 262144)), float(a.get("W", 0.5)),
                                    float(a.get("M", 4.5)), float(a.get("A", 0)))
            except ValueError:
                pass
        # log, hyperlog, fasinh: not met yet, left out rather than guessed
    return out


class _Unread(Exception):
    """A gate left alone; the message is the reason."""


def _vertices(el) -> np.ndarray:
    return np.array([[float(_attrs(c)["value"]) for c in _children(v, "coordinate")]
                     for v in (_children(el, "vertex") if el is not None else [])])


def _read_gate(el, transforms=None):
    shape = next(iter(el), None) if el is not None else None
    if shape is None:
        raise _Unread("the workspace holds no gate for it")
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
        # FlowJo stores the vertices in the file's units and tests the
        # events where the polygon was drawn: on the axis' own scale
        scales = [t if isinstance(t, (Biex, Logicle)) else None
                  for t in ((transforms or {}).get(d) for d in dims)]
        return PolygonGate(dims, _vertices(shape), inside,
                           scales if any(s is not None for s in scales) else None)
    if kind == "EllipsoidGate":
        return _read_ellipse(shape, dims, inside, transforms)
    if kind == "CurlyQuad":
        return _read_curly(dims, mins, maxs, inside, transforms)
    raise _Unread(f"a kind of gate ({kind}) that is not read")


def curl_scale(transform):
    """The scale a curly quadrant's curve is measured on, from the axis'
    transform in the workspace: a biexponential or a Logicle as it is; a
    linear axis from 0 — cytolib scales it by its top alone, whatever its
    bottom — to its top; an axis with none, linear to 262,144."""
    if isinstance(transform, (Biex, Logicle)):
        return transform
    return Linear(0.0, transform.max if isinstance(transform, Linear) else CURL_TOP)


def _read_curly(dims, mins, maxs, inside, transforms) -> CurlyQuadrantGate:
    """FlowJo's curly quadrant: a corner of a crossing whose arms above
    and to the right of the centre bend.

    What the workspace holds is a quadrant's, and nothing else: `gating:
    CurlyQuad` with one `gating:dimension` for each axis, bounded on one
    side only, `min` where the corner is the high side of the axis and
    `max` where it is the low one — so the crossing point and which corner
    (CytoML, inst/include/CytoML/winFlowJoWorkspace.hpp,
    getGate(wsCurlyQuadGateNode&); and so in the FlowJo 10.10 workspace
    allotropy publishes, tests/parsers/flowjo/testdata/
    flowjo_example_2.wsp). No curve, no constant: FlowJo draws the curve
    itself, from a model of photon-counting spread (docs.flowjo.com,
    Curly Quadrant: p2 = a·√p1 + b, `a` its "curliness") whose `a` is not
    written. So the curve FlowJo draws is not known, and the one given is
    CytoML's (cytolib src/gate.cpp, CurlyQuadGate::interpolate; see
    `gating.CURL_BEND`), checked against CytoML itself (tests/
    test_curly_quadrant.py) and not against FlowJo. FlowJo's own count
    is beside each corner: `count_mismatches` says where it is not met."""
    if transforms is None:
        raise _Unread("a curly quadrant bends on the scale of each axis, and a group's "
                      "gates come without the samples' scales")
    if len(dims) != 2 or any((lo is None) == (hi is None) for lo, hi in zip(mins, maxs)):
        raise _Unread("a curly quadrant without one bound on each of two axes")
    high = tuple(lo is not None for lo in mins)
    at = [lo if lo is not None else hi for lo, hi in zip(mins, maxs)]
    corner = {v: k for k, v in QUADRANTS.items()}[high]
    return CurlyQuadrantGate(dims, at, corner, "", inside,
                             [curl_scale(transforms.get(d)) for d in dims])


def _join_curly(parent: GateNode) -> None:
    """The corners of one curly crossing under a parent share a `group`,
    as the four of a plain quadrant do, so that moving one moves them
    all: the same two axes, the same centre."""
    for node in parent.children:
        g = node.gate
        if isinstance(g, CurlyQuadrantGate):
            g.group = "|".join([parent.path, "curly", *g.dims, *map(repr, g.center)])


def _read_ellipse(shape, dims, inside, transforms) -> EllipseGate:
    """FlowJo's ellipse is not Gating-ML's (a mean and a covariance): it is
    two foci and four points of the edge, and unlike every other gate they
    are not in the file's units but on the 256 channels of the plot — each
    axis' own scale, 0 at the bottom and 256 at the top.

      - FlowKit, src/flowkit/_models/gates/_wsp_gates.py,
        WSPEllipsoidGate.convert_to_polygon_gate: divides by 256 into the
        transform's own range; centre and angle from the foci, the long
        half-axis from the first pair of edge points, the short one from
        that and the foci.
      - cytolib, src/gate.cpp, ellipsoidGate::transforming, with CytoML,
        inst/include/CytoML/winFlowJoWorkspace.hpp, getGate(wsEllipseGateNode&)
        and cytolib src/ellipse2points.cpp, parseEllipse: the same 256
        channels; centre, both half-axes and angle from the four edge
        points alone, the longer pair being the long axis.

    The ellipse is read when both readings give the same one to within a
    channel, and left alone otherwise. In the two workspaces of FlowJo's
    looked at they give the same one exactly: the edge points are in whole
    channels and the foci are those of the two half-lengths — though the
    second pair of points is not at right angles to the first. Measured once, 2026-11-23, on the workspace FlowKit
    publishes with its file (data/simple_line_example/
    single_ellipse_51_events.wsp, FlowJo 10.6.1, linear axes, 100 events):
    this gives FlowJo's 51; half a channel more or less on the long axis
    gives 52 or 49; taking the `distance` attribute for the sum of the
    focal distances, which neither reader does, gives 66. Nothing was
    measured on a biexponential axis, where the two readers also build
    FlowJo's table at different resolutions (4096 and 256 steps) and can
    differ by most of a channel near zero."""
    if transforms is None:
        raise _Unread("an ellipse is stored on the scale of each axis, and a group's "
                      "gates come without the samples' scales")
    scales = []
    for d in dims:
        t = transforms.get(d)
        # FlowKit takes only 256-channel biexponentials (wsp_utils.py,
        # _parse_wsp_transforms); on a linear axis it adds minRange where
        # CytoML ignores it, so the two agree only when it is zero
        if not (isinstance(t, Biex) and t.length == 256
                or isinstance(t, Linear) and t.min == 0):
            raise _Unread(f"an ellipse is stored on the scale of each axis, and that of "
                          f"{d} is not one its conversion is known for")
        scales.append(t)
    foci, edge = _vertices(_child(shape, "foci")), _vertices(_child(shape, "edge"))
    if len(dims) != 2 or foci.shape != (2, 2) or edge.shape != (4, 2):
        raise _Unread("an ellipse without two foci and four edge points on two axes")
    first, second = edge[:2], edge[2:]
    if np.hypot(*(first[0] - first[1])) < np.hypot(*(second[0] - second[1])):
        first, second = second, first
    centre = foci.mean(axis=0)
    along = foci[1] - foci[0]
    c = np.hypot(*along) / 2
    a = np.hypot(*(first[0] - first[1])) / 2
    if c == 0:
        along = first[1] - first[0]
    if a <= c or a == 0:
        raise _Unread("an ellipse whose foci do not lie inside its edge points")
    b = np.sqrt(a * a - c * c)
    along = along / np.hypot(*along)
    arm = first[0] - centre
    if (np.abs(edge.mean(axis=0) - centre).max() > 1
            or abs(arm[0] * along[1] - arm[1] * along[0]) > 1
            or abs(np.hypot(*(second[0] - second[1])) / 2 - b) > 1):
        raise _Unread("its foci and its edge points do not describe the same ellipse "
                      "to within a channel")
    return EllipseGate(dims, [float(v) / 256 for v in centre], [float(a) / 256, float(b) / 256],
                       float(np.arctan2(along[1], along[0])), inside, scales)


def _join_quadrants(parent: GateNode) -> None:
    """FlowJo writes a quadrant gate as four populations, each a rectangle
    with one bound on each axis and the other side open; nothing else ties
    them together (no `quadId` in the workspaces of FlowJo 10.6 and 10.8
    that FlowKit publishes: data/simple_diamond_example/
    simple_diamond_example_quad_gate.wsp, data/8_color_data_set/
    reused_quad_gate_with_child.wsp).

      - CytoML, inst/include/CytoML/winFlowJoWorkspace.hpp,
        getGate(wsRectGateNode&): a rectangle on two axes with one bound
        missing on each is a quadGate, the missing sides giving the
        corner — no min on X and no max on Y is Q1, and so on round.
      - FlowKit keeps them as the rectangles they are (tests/
        workspace_tests.py, test_load_wsp_single_quad) and gets FlowJo's
        four counts.

    So the events of each corner are the rectangle's, as before. What is
    added is knowing the four belong together, and only when they do:
    four sisters on the same two axes, bounded at the same two values,
    one in each corner."""
    corners = {v: k for k, v in QUADRANTS.items()}
    sets: dict[tuple, dict[tuple, GateNode | None]] = {}
    for node in parent.children:
        g = node.gate
        if not (isinstance(g, RectangleGate) and len(g.dims) == 2 and g.inside
                and all((lo is None) != (hi is None) for lo, hi in zip(g.mins, g.maxs))):
            continue
        high = tuple(lo is not None for lo in g.mins)
        at = tuple(lo if lo is not None else hi for lo, hi in zip(g.mins, g.maxs))
        members = sets.setdefault((tuple(g.dims), at), {})
        # a corner taken twice is not a quadrant gate: None spoils the set
        members[high] = None if high in members else node
    for (dims, at), members in sets.items():
        if len(members) != 4 or None in members.values():
            continue
        group = "|".join([parent.path, *dims, repr(at[0]), repr(at[1])])
        for high, node in members.items():
            node.gate = QuadrantGate(list(dims), at[0], at[1], corners[high], group)
            node.note = f"{UNVERIFIED}: four open rectangles taken as one quadrant gate"


_COMBINATIONS = {"AndNode": AND, "OrNode": OR, "NotNode": NOT}


def _read_combination(p, parent: GateNode, sisters: list[str], known: list[str]) -> BooleanGate:
    """A combination is a node of its own beside the populations, naming
    what it is made of by path from the sample.

      - FlowKit, src/flowkit/_utils/wsp_utils.py, _parse_boolean_node and
        _recurse_wsp_sub_populations: AndNode, OrNode, NotNode; each
        `Dependents/Dependent name=` a path split on `/`.
      - CytoML, inst/include/CytoML/winFlowJoWorkspace.hpp,
        getGate(wsPopNode&): the same; a NotNode has exactly one
        Dependent, the others at least two.

    The counts in FlowKit's data/8_color_data_set/
    8_color_ICS_boolean_gate_testing.wsp (FlowJo 10.8.1) fit: of 81,992
    CD4+ with 74 IL2+, the NotNode of IL2+ has 81,918; 62 and 607 with 6 in
    both have 663 in either. A `BooleanGate` combines sisters, so that is
    what is read; anything made of populations elsewhere is left."""
    deps_el = _child(p, "Dependents")
    deps = [d.get("name", "") for d in (_children(deps_el, "Dependent")
                                        if deps_el is not None else [])]
    op = _COMBINATIONS[_tag(p)]
    if (len(deps) != 1) if op == NOT else (len(deps) < 2):
        raise _Unread(f"it names {len(deps)} populations, which is not what "
                      f"{'a negation' if op == NOT else 'a combination'} is made of")
    prefix = parent.path + "/" if parent.path else ""
    operands = []
    for dep in deps:
        name = dep[len(prefix):] if dep.startswith(prefix) else None
        if name not in sisters:
            raise _Unread(f"it is made of {dep}, which is not one of its sisters")
        if sisters.count(name) > 1:
            raise _Unread(f"it is made of {dep}, and two of its sisters have that name")
        if name not in known:
            raise _Unread(f"it is made of {dep}, which was not read, or is a "
                          f"combination listed after it")
        operands.append(name)
    return BooleanGate(op, operands)


def _read_populations(el, parent: GateNode, transforms=None,
                      skipped: list[Skipped] | None = None, owner: str = "") -> None:
    def is_node(e) -> bool:
        return _tag(e) == "Population" or _tag(e) in _COMBINATIONS

    nodes = [p for sub in _children(el, "Subpopulations") for p in sub if is_node(p)]
    sisters = [p.get("name", "").strip() for p in nodes]
    # the drawn gates first, as `GateNode.apply` takes them, so that a
    # combination can be told whether what it names was read
    gates: list = []
    for p in nodes:
        try:
            gates.append(_read_gate(_child(p, "Gate"), transforms)
                         if _tag(p) == "Population" else None)
        except _Unread as why:
            gates.append(why)
    known = [name for name, g in zip(sisters, gates)
             if g is not None and not isinstance(g, _Unread)]
    read = []
    for p, name, gate in zip(nodes, sisters, gates):
        note = ""
        if gate is None:
            try:
                gate = _read_combination(p, parent, sisters, known)
                known.append(name)
                note = f"{UNVERIFIED}: a combination of its sister populations"
            except _Unread as why:
                gate = why
        if isinstance(gate, _Unread):
            if skipped is not None:
                gate_el = _child(p, "Gate")
                shape = next(iter(gate_el), p) if gate_el is not None else p
                skipped.append(Skipped(owner, "/".join(filter(None, [parent.path, name])),
                                       _tag(shape), str(gate),
                                       sum(is_node(d) for d in p.iter()) - 1))
            continue
        if isinstance(gate, EllipseGate):
            note = f"{UNVERIFIED}: an ellipse converted from FlowJo's 256-channel grid"
        elif isinstance(gate, CurlyQuadrantGate):
            note = (f"{UNVERIFIED}: a curly quadrant, whose curve the workspace does not hold; "
                    "it is drawn as CytoML draws it")
        read.append((p, parent.add(GateNode(name=name, gate=gate, note=note,
                                            reference_count=_int(p.get("count"))))))
    _join_quadrants(parent)
    _join_curly(parent)
    for p, node in read:
        _read_populations(p, node, transforms, skipped, owner)
