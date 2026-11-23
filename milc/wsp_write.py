"""
The gates written as a FlowJo 10 workspace (.wsp): the files, the matrix
each is compensated with, the scale of each parameter and the gating tree
of each sample, in the form `wsp.py` reads — which is the form FlowJo
10.10 wrote in the one workspace of its own this was modelled on.

It has NOT been opened in FlowJo. What is checked is that it comes back:
read by `wsp.read_wsp`, every population of every sample holds the events
it holds here; and FlowKit (another reader of the format) agrees.

How each thing is written:

  - a rectangle, its bounds as they are; a quadrant, as FlowJo does: a
    rectangle open on two sides, one for each corner.
  - a polygon, its vertices in the file's units. FlowJo tests it on the
    scale of each axis, so the scale it was drawn on is the one given to
    the parameter; drawn on another (or on one FlowJo is not known to
    test on: arcsinh, log), it gets vertices added along its sides until
    they are the same line on any scale.
  - an ellipse, on FlowJo's 256 channels: two foci and four points of
    the edge. That is defined only for linear and biexponential axes; on
    others it becomes a polygon.
  - a combination, as an AndNode, OrNode or NotNode naming its sisters.
  - Time, in seconds: FlowJo multiplies the file's Time by `$TIMESTEP`
    (its own range for Time in a workspace says so), and so does this.
"""

from __future__ import annotations

import re
import time
import uuid
import xml.etree.ElementTree as ET

import numpy as np

from .gate_export import (WSP, ExportSample, Prepared, Report, Unwritable, affine, axis_scales,
                          bend, counted, ellipse_axes, ellipse_covariance, ellipse_points,
                          is_helper, off, on, prepare, scales_of, walk)
from .gating import (AND, NOT, OR, QUADRANTS, BooleanGate, EllipseGate, GateNode, PolygonGate,
                     QuadrantGate, RectangleGate)
from .transforms import Arcsinh, Biex, Linear, Log, Logicle, same

#: the version this claims to be written by: old enough for any FlowJo 10
#: in use to open, and what other readers look for
FLOWJO_VERSION = "10.8.1"
GROUP = "All Samples"

NS = {"gating": "http://www.isac-net.org/std/Gating-ML/v2.0/gating",
      "transforms": "http://www.isac-net.org/std/Gating-ML/v2.0/transformations",
      "data-type": "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes",
      "xsi": "http://www.w3.org/2001/XMLSchema-instance"}
SCHEMAS = " ".join(f"{NS[k]} http://www.isac-net.org/std/Gating-ML/v2.0/gating/{f}.v2.0.xsd"
                   for k, f in (("gating", "Gating-ML"), ("transforms", "Transformations"),
                                ("data-type", "DataTypes"))) + " "
_NODES = {AND: "AndNode", OR: "OrNode", NOT: "NotNode"}
#: what FlowJo writes on every gate it draws
_DRAWN = {"annoOffsetX": "0", "annoOffsetY": "0", "tint": "#000000", "isTinted": "0",
          "lineWeight": "Normal", "userDefined": "1"}
_UNSAFE = re.compile("[^\t\n\r\x20-\ud7ff\ue000-\ufffd\U00010000-\U0010ffff]")


def q(prefix: str, name: str) -> str:
    return f"{{{NS[prefix]}}}{name}"


def num(value) -> str:
    return repr(float(value))


def text(value) -> str:
    """Without the characters XML cannot hold (a keyword may have any)."""
    return _UNSAFE.sub("", str(value))


def file_uri(path: str) -> str:
    """`file:/C:/folder/a.fcs`, as Java — FlowJo — writes a file's URI."""
    import urllib.parse
    path = path.replace("\\", "/")
    return "file:" + urllib.parse.quote(path if path.startswith("/") else "/" + path, safe="/:")


class _Sample:
    """One sample being written: its names, scales and identifiers."""

    def __init__(self, prepared: Prepared, sample_id: int, report: Report, ids):
        self.sample, self.tree, self.counts = prepared.sample, prepared.tree, prepared.counts
        self.id, self.report, self.ids = str(sample_id), report, ids
        self.scales = axis_scales(self.sample, self.tree)
        self.written = 0

    # --- parameters ----------------------------------------------------------------

    def name(self, dim: str) -> str:
        channel, compensated = self.sample.dimension(dim)
        return "Comp-" + channel if compensated else channel

    def scale(self, dim: str):
        """The scale the parameter is given in the file."""
        if dim in self.scales:
            return self.scales[dim]
        if "Comp-" + dim in self.scales:     # as FlowJo: one scale, compensated or not
            return self.scales["Comp-" + dim]
        return self.sample.scales.get(self.sample.dimension(dim)[0])

    def unit(self, dim: str) -> float:
        """What FlowJo multiplies the file's values by: seconds for Time."""
        return self.sample.timestep if self.sample.dimension(dim)[0] == self.sample.time else 1.0

    def top(self, dim: str) -> float:
        return self.sample.range_of(self.sample.dimension(dim)[0])

    def transform(self, dim: str, name: str) -> ET.Element:
        scale, unit = self.scale(dim), self.unit(dim)
        if isinstance(scale, Biex):
            el = ET.Element(q("transforms", "biex"), {
                q("transforms", "length"): str(scale.length),
                q("transforms", "maxRange"): num(scale.max_range),
                q("transforms", "neg"): num(scale.neg), q("transforms", "width"): num(scale.width),
                q("transforms", "pos"): num(scale.pos)})
        elif isinstance(scale, Logicle):
            el = ET.Element(q("transforms", "logicle"), {
                q("transforms", "length"): "256", q("transforms", "T"): num(scale.T),
                q("transforms", "A"): num(scale.A), q("transforms", "W"): num(scale.W),
                q("transforms", "M"): num(scale.M)})
        elif isinstance(scale, Arcsinh):
            # T / sinh(M ln 10) is the cofactor
            top = self.top(dim)
            el = ET.Element(q("transforms", "fasinh"), {
                q("transforms", "length"): "256", q("transforms", "maxRange"): num(top),
                q("transforms", "T"): num(top), q("transforms", "A"): "0.0",
                q("transforms", "M"): num(np.arcsinh(top / scale.cofactor) / np.log(10.0)),
                q("transforms", "W"): "0.0"})
        elif isinstance(scale, Log):
            el = ET.Element(q("transforms", "log"), {
                q("transforms", "offset"): num(scale.min),
                q("transforms", "decades"): num(np.log10(scale.max / scale.min))})
        else:
            lo, hi = ((scale.min, scale.max) if isinstance(scale, Linear)
                      else (0.0, self.top(dim)))
            el = ET.Element(q("transforms", "linear"), {
                q("transforms", "minRange"): num(lo * unit),
                q("transforms", "maxRange"): num(hi * unit), "gain": "1"})
        ET.SubElement(el, q("data-type", "parameter"), {q("data-type", "name"): name})
        return el

    def tested_on(self, dim: str):
        """The scale a polygon on this parameter is tested on by whoever
        reads the workspace: the parameter's, where it is a biexponential
        or a Logicle; the file's units otherwise."""
        scale = self.scale(dim)
        return scale if isinstance(scale, (Biex, Logicle)) else None

    # --- gates -----------------------------------------------------------------------

    def dimension(self, parent: ET.Element, dim: str, **bounds) -> None:
        el = ET.SubElement(parent, q("gating", "dimension"),
                           {q("gating", k): num(v * self.unit(dim))
                            for k, v in bounds.items() if v is not None})
        ET.SubElement(el, q("data-type", "fcs-dimension"), {q("data-type", "name"): self.name(dim)})

    def rectangle(self, dims, mins, maxs, inside: bool) -> ET.Element:
        attrs = {"eventsInside": "1" if inside else "0", **_DRAWN}
        if len(dims) == 2:
            attrs.update(percentX="0", percentY="0")
        el = ET.Element(q("gating", "RectangleGate"), attrs)
        for dim, lo, hi in zip(dims, mins, maxs):
            self.dimension(el, dim, min=lo, max=hi)
        return el

    @staticmethod
    def vertex(parent: ET.Element, point) -> None:
        v = ET.SubElement(parent, q("gating", "vertex"))
        for value in point:
            ET.SubElement(v, q("gating", "coordinate"), {q("data-type", "value"): num(value)})

    def polygon(self, node: GateNode, dims, vertices, scales, inside: bool) -> ET.Element:
        """`vertices` in the file's units, straight on `scales`."""
        vertices = np.asarray(vertices, dtype=np.float64)
        tested = [self.tested_on(d) for d in dims]
        forward, back = [], []
        for drawn, target in zip(scales, tested):
            alike = (affine(drawn) is not None if target is None else same(drawn, target))
            forward.append(None if alike else (lambda s, a=drawn, b=target: on(b, off(a, s))))
            back.append(None if alike else (lambda s, a=drawn, b=target: on(a, off(b, s))))
        if any(f is not None for f in forward):
            where = np.column_stack([on(s, vertices[:, k]) for k, s in enumerate(scales)])
            moved = bend(where, forward, back)
            out = np.column_stack([off(drawn if f is None else target, moved[:, k])
                                   for k, (drawn, target, f)
                                   in enumerate(zip(scales, tested, forward))])
            if len(out) > len(vertices):
                self.report.note(self.sample.name, getattr(node, "source", None) or node.path,
                                 "its sides are straight on a scale the workspace does not "
                                 f"test it on: vertices were added along them "
                                 f"({len(vertices)} became {len(out)})")
            vertices = out
        el = ET.Element(q("gating", "PolygonGate"),
                        {"eventsInside": "1" if inside else "0", **_DRAWN,
                         "quadId": "-1", "gateResolution": "256"})
        for dim in dims:
            self.dimension(el, dim)
        units = np.array([self.unit(d) for d in dims])
        for point in vertices * units:
            self.vertex(el, point)
        return el

    def ellipse(self, node: GateNode) -> ET.Element:
        g = node.gate
        scales = scales_of(g)
        stretch, shift = [], []
        for dim, drawn in zip(g.dims, scales):
            given = self.scale(dim)
            line = affine(drawn)
            if isinstance(given, Biex) and given.length == 256 and same(drawn, given):
                stretch.append(256.0)
                shift.append(0.0)
            elif line is not None and (given is None or isinstance(given, Linear)
                                       and given.min == 0):
                top = given.max if isinstance(given, Linear) else self.top(dim)
                stretch.append(256.0 * line[0] / top)
                shift.append(256.0 * line[1] / top)
            else:
                break
        else:
            return self._ellipse_256(g, np.array(stretch), np.array(shift))
        # FlowJo's ellipse lives on 256 channels of a linear or a
        # biexponential axis; on any other it is not known what it means
        self.report.note(self.sample.name, getattr(node, "source", None) or node.path,
                         "an ellipse on an axis that is neither linear nor FlowJo's "
                         "biexponential is written as a polygon of "
                         f"{len(ellipse_points(g.center, g.radii, g.angle))} sides")
        where = ellipse_points(g.center, g.radii, g.angle)
        vertices = np.column_stack([off(s, where[:, k]) for k, s in enumerate(scales)])
        return self.polygon(node, g.dims, vertices, scales, g.inside)

    def _ellipse_256(self, g: EllipseGate, stretch, shift) -> ET.Element:
        centre = np.asarray(g.center, dtype=np.float64) * stretch + shift
        a, b, angle = ellipse_axes(ellipse_covariance(g.radii, g.angle, stretch))
        along = np.array([np.cos(angle), np.sin(angle)])
        across = np.array([-along[1], along[0]])
        c = np.sqrt(max(a * a - b * b, 0.0))
        # as FlowJo orders them: the ends of the long axis, that by the
        # second focus first, then the ends of the short one
        el = ET.Element(q("gating", "EllipsoidGate"),
                        {"eventsInside": "1" if g.inside else "0", **_DRAWN,
                         q("gating", "distance"): num(2 * a)})
        for dim in g.dims:
            self.dimension(el, dim)
        foci = ET.SubElement(el, q("gating", "foci"))
        for point in (centre - c * along, centre + c * along):
            self.vertex(foci, point)
        edge = ET.SubElement(el, q("gating", "edge"))
        for point in (centre + a * along, centre - a * along,
                      centre - b * across, centre + b * across):
            self.vertex(edge, point)
        return el

    def shape(self, node: GateNode) -> ET.Element:
        g = node.gate
        if isinstance(g, RectangleGate):
            return self.rectangle(g.dims, g.mins, g.maxs, g.inside)
        if isinstance(g, QuadrantGate):
            high = QUADRANTS[g.quadrant]
            at = (g.x, g.y)
            return self.rectangle(g.dims, [v if h else None for v, h in zip(at, high)],
                                  [None if h else v for v, h in zip(at, high)], g.inside)
        if isinstance(g, PolygonGate):
            return self.polygon(node, g.dims, g.vertices, scales_of(g), g.inside)
        if isinstance(g, EllipseGate):
            return self.ellipse(node)
        raise Unwritable(f"a kind of gate ({type(g).__name__}) a workspace does not hold")

    # --- populations -------------------------------------------------------------------

    def graph(self, node: GateNode) -> ET.Element:
        """The plot FlowJo opens on the population: the axes of the first
        gate drawn in it, or of its own."""
        inner = next((c.gate for c in node.children
                      if c.gate is not None and not isinstance(c.gate, BooleanGate)), None)
        gate = inner or (None if isinstance(node.gate, BooleanGate) else node.gate)
        try:
            names = [self.name(d) for d in (gate.dims if gate is not None else [])]
        except Unwritable:
            names = []
        names = (names + ["", ""])[:2]
        el = ET.Element("Graph", {"smoothing": "0", "backColor": "#ffffff", "foreColor": "#000000",
                                  "type": "Pseudocolor" if names[1] else "Histogram", "fast": "1"})
        for axis, name in zip("xy", names):
            ET.SubElement(el, "Axis", {"dimension": axis, "name": name, "label": "",
                                       "auto": "auto"})
        return el

    def count(self, node: GateNode) -> str:
        return str(self.counts[node.path]) if self.counts and node.path in self.counts else "-1"

    def population(self, node: GateNode, parent_gate: str | None, owner: str):
        """(the population's element, the identifier of its gate)."""
        g = node.gate
        attrs = {"name": text(node.name), "annotation": "", "owningGroup": owner,
                 "expanded": "1", "sortPriority": "10", "count": self.count(node)}
        if isinstance(g, BooleanGate):
            el = ET.Element(_NODES[g.op], attrs)
            el.append(self.graph(node))
            deps = ET.SubElement(el, "Dependents")
            before = node.parent.path + "/" if node.parent is not None and node.parent.path else ""
            for name in g.operands:
                ET.SubElement(deps, "Dependent", {"name": text(before + name)})
            gate_id = None
        else:
            shape = self.shape(node)     # first: a gate that fails leaves nothing behind
            el = ET.Element("Population", attrs)
            el.append(self.graph(node))
            gate_id = next(self.ids)
            gate = ET.SubElement(el, "Gate", {q("gating", "id"): gate_id})
            if parent_gate:
                gate.set(q("gating", "parent_id"), parent_gate)
            gate.append(shape)
        return el, gate_id

    def populations(self, parent_el: ET.Element, parent: GateNode, parent_gate: str | None,
                    owner: str = "") -> None:
        made = walk(parent, lambda n: self.population(n, parent_gate, owner),
                    self.sample.name, self.report)
        if not made:
            return
        sub = ET.SubElement(parent_el, "Subpopulations")
        for node, (el, gate_id) in made:
            sub.append(el)
            if not is_helper(node):
                self.written += 1
            self.populations(el, node, gate_id, owner)

    # --- the sample ----------------------------------------------------------------------

    def matrix(self, inside_sample: bool = True) -> ET.Element | None:
        comp = self.sample.compensation
        if comp is None:
            return None
        own = self.sample.matrix_source == "file"
        attrs = {"spectral": "0", "gainNormalized": "0", "weightOptAlgorithmType": "OLS",
                 "prefix": "Comp-", "name": text(matrix_name(self.sample)),
                 "editable": "0" if own else "1",
                 "matrixType": "acquisitionDefined" if own else "userDefined",
                 "color": "#c0c0c0" if own else "#00ccff",
                 "version": "FlowJo-" + FLOWJO_VERSION, "status": "FINALIZED",
                 q("transforms", "id"): matrix_id(comp), "suffix": ""}
        if inside_sample:
            attrs = {"sampleID": file_uri(self.sample.path), **attrs}
        el = ET.Element(q("transforms", "spilloverMatrix"), attrs)
        params = ET.SubElement(el, q("data-type", "parameters"))
        for c in comp.channels:
            ET.SubElement(params, q("data-type", "parameter"),
                          {q("data-type", "name"): c, "userProvidedCompInfix": "Comp-" + c})
        for c, row in zip(comp.channels, np.asarray(comp.matrix, dtype=np.float64)):
            spill = ET.SubElement(el, q("transforms", "spillover"),
                                  {q("data-type", "parameter"): c,
                                   "userProvidedCompInfix": "Comp-" + c})
            for other, value in zip(comp.channels, row):
                ET.SubElement(spill, q("transforms", "coefficient"),
                              {q("data-type", "parameter"): other,
                               q("transforms", "value"): num(value)})
        return el

    def element(self, owner: str = "") -> ET.Element:
        s = self.sample
        el = ET.Element("Sample")
        ET.SubElement(el, "DataSet", {"uri": file_uri(s.path), "sampleID": self.id})
        matrix = self.matrix()
        if matrix is not None:
            el.append(matrix)
        scales = ET.SubElement(el, "Transformations")
        for c in s.channels:
            scales.append(self.transform(c, c))
        for c in (s.compensation.channels if s.compensation is not None else []):
            if c in s.channels and "Comp-" + c not in s.channels:
                scales.append(self.transform("Comp-" + c, "Comp-" + c))
        keywords = ET.SubElement(el, "Keywords")
        for key, value in s.keywords.items():
            if text(key):
                ET.SubElement(keywords, "Keyword", {"name": text(key), "value": text(value)})
        node = ET.SubElement(el, "SampleNode", {
            "name": text(s.name), "annotation": "", "owningGroup": "", "expanded": "1",
            "sortPriority": "10", "count": str(s.events), "sampleID": self.id})
        node.append(self.graph(self.tree))
        self.populations(node, self.tree, None, owner)
        return el


def matrix_name(sample: ExportSample) -> str:
    if sample.matrix_source == "file":
        return "Acquisition-defined"
    return sample.compensation.name or "MIL-C"


def matrix_id(comp) -> str:
    """One identifier for one matrix, whichever sample it is written for."""
    key = ",".join(comp.channels) + "|" + ",".join(
        num(v) for v in np.asarray(comp.matrix, dtype=np.float64).ravel())
    return str(uuid.uuid5(uuid.NAMESPACE_OID, key))


def _gate_ids():
    n = 100000000
    while True:
        n += 1
        yield f"ID{n}"


def _plain(el: ET.Element | None) -> str:
    """A sample's populations as text, without what differs between two
    samples gated alike: the identifiers and the counts."""
    if el is None:
        return ""
    return re.sub(r' (?:count|\w+:id|\w+:parent_id)="[^"]*"', "", ET.tostring(el, encoding="unicode"))


def write_wsp(samples: list[ExportSample], path: str) -> Report:
    """Write a workspace holding these samples. Never raises for a gate
    or a sample that cannot be written: they are in the report."""
    for prefix, uri in NS.items():
        ET.register_namespace(prefix, uri)
    report = Report(format=WSP)
    ids = _gate_ids()
    writers = []
    for sample in samples:
        try:
            tree = prepare(sample, report)
            counts = counted(sample, tree, report)
        except Exception as exc:
            report.note(sample.name, "", f"its gates could not be written: {exc}")
            tree, counts = GateNode(sample.gates.name), None
        writers.append(_Sample(Prepared(sample, tree, counts), len(writers) + 1, report, ids))

    root = ET.Element("Workspace", {
        "version": "20.0", "modDate": time.strftime("%a %b %d %H:%M:%S %Z %Y"),
        "flowJoVersion": FLOWJO_VERSION, "curGroup": GROUP,
        q("xsi", "schemaLocation"): SCHEMAS})
    matrices = ET.SubElement(root, "Matrices")
    groups = ET.SubElement(root, "Groups")
    sample_list = ET.SubElement(root, "SampleList")

    elements = []
    for w in writers:
        try:
            elements.append(w.element())
        except Exception as exc:
            report.note(w.sample.name, "", f"the sample could not be written: {exc}")
            elements.append(None)
    kept = [(w, el) for w, el in zip(writers, elements) if el is not None]

    # the gates every sample has in common are the group's too, as FlowJo
    # keeps them — except with an ellipse, which a group cannot be read to
    # hold (it is stored on the scales of a sample)
    trees = [_plain(el.find("SampleNode/Subpopulations")) for _, el in kept]
    shared = bool(trees) and bool(trees[0]) and all(t == trees[0] for t in trees) \
        and "EllipsoidGate" not in trees[0]
    seen = set()
    for w, el in kept:
        if shared:
            for p in el.iter():
                if "owningGroup" in p.attrib and p.tag != "SampleNode":
                    p.set("owningGroup", GROUP)
        sample_list.append(el)
        report.written += w.written
        comp = w.sample.compensation
        if comp is not None and matrix_id(comp) not in seen:
            seen.add(matrix_id(comp))
            matrices.append(w.matrix(inside_sample=False))

    group = ET.SubElement(groups, "GroupNode", {
        "name": GROUP, "annotation": "", "owningGroup": GROUP, "expanded": "1",
        "sortPriority": "10", "count": "-1"})
    if shared:
        first = kept[0][0]
        quiet = Report()
        _Sample(Prepared(first.sample, first.tree, first.counts), 0, quiet, ids).populations(
            group, first.tree, None, GROUP)
    members = ET.SubElement(group, "Group", {
        "name": GROUP, "live": "1", "role": "ws.group.dlog.test", "key": "", "synchronized": "0",
        "foreground": "#000000", "fontStyle": "bold"})
    ET.SubElement(members, "Criteria")
    refs = ET.SubElement(members, "SampleRefs")
    for w, _ in kept:
        ET.SubElement(refs, "SampleRef", {"sampleID": w.id})

    report.samples = len(kept)
    ET.indent(root, space="  ")
    ET.ElementTree(root).write(path, encoding="UTF-8", xml_declaration=True)
    report.files.append(path)
    return report
