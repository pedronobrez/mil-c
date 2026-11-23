"""
The gates written as Gating-ML 2.0, the standard of ISAC for describing
gates (Spidlen et al., Cytometry A 2015): the compensation matrix, the
scales and every population, each naming its parent.

Gating-ML has no samples in it. So there is one file where every sample
has the same gates and the same matrix, and one file for each sample
where they differ.

How each thing is written:

  - compensation: a `spectrumMatrix`, its fluorochromes named as its
    detectors; a dimension compensated here refers to it, any other is
    `uncompensated`.
  - a rectangle: `RectangleGate`, its bounds in the file's units, no
    transformation. A range of time likewise.
  - a polygon: `PolygonGate`, its vertices on the scale it was drawn on,
    which is written as a transformation: `logicle` for a Logicle,
    `fasinh` for an arcsinh (T / sinh(M ln 10) is the cofactor), `flog`
    for a logarithm; none for a linear axis, where the vertices are in the
    file's units.
  - FlowJo's biexponential has no form in Gating-ML. It is written as the
    Logicle nearest to it (`logicle_like`: a few 1e-6 of the scale away),
    and the vertices are placed — more of them where a side would bend —
    so that the polygon holds the events it holds here, not nearly those.
  - an ellipse: `EllipsoidGate`, a centre and a covariance matrix with a
    distance of 1 — except on a biexponential axis, where it is not an
    ellipse on any scale Gating-ML has, and becomes a polygon.
  - quadrants: one `QuadrantGate` for the corners made together, each a
    `Quadrant` that the populations under it name as their parent.
  - a combination: `BooleanGate`. A gate that keeps what is outside it is
    the gate and the `not` of it.

Each gate carries its name, its path and its count here in a
`custom_info`, since a Gating-ML identifier cannot hold a name as it is.
Checked against FlowKit, which reads the file and gates the same events;
valid against the schema of the standard.
"""

from __future__ import annotations

import os
import re
import time
import xml.etree.ElementTree as ET

import numpy as np

from . import __version__
from .gate_export import (GATINGML, ExportSample, Prepared, Report, Unwritable, affine, bend,
                          counted, ellipse_covariance, ellipse_points, is_helper, logicle_like,
                          off, on, prepare, scales_of, source_of, walk)
from .gating import (AND, NOT, OR, QUADRANTS, BooleanGate, EllipseGate, GateNode, PolygonGate,
                     QuadrantGate, RectangleGate)
from .transforms import Arcsinh, Biex, Log, Logicle, to_spec

NS = {"gating": "http://www.isac-net.org/std/Gating-ML/v2.0/gating",
      "transforms": "http://www.isac-net.org/std/Gating-ML/v2.0/transformations",
      "data-type": "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes",
      "xsi": "http://www.w3.org/2001/XMLSchema-instance"}
SCHEMAS = " ".join(f"{NS[k]} http://www.isac-net.org/std/Gating-ML/v2.0/gating/{f}.v2.0.xsd"
                   for k, f in (("gating", "Gating-ML"), ("transforms", "Transformations"),
                                ("data-type", "DataTypes")))
UNCOMPENSATED = "uncompensated"
MATRIX = "Compensation"
_UNSAFE = re.compile("[^\t\n\r\x20-퟿-�\U00010000-\U0010ffff]")
#: the top given to an arcsinh: any will do, the cofactor is what counts
_ASINH_TOP = 262144.0


def q(prefix: str, name: str) -> str:
    return f"{{{NS[prefix]}}}{name}"


def num(value) -> str:
    return repr(float(value))


def text(value) -> str:
    return _UNSAFE.sub("", str(value))


class _Axis:
    """How one scale of MIL-C is written: the transformation a dimension
    refers to (None: the file's units) and where a position on the scale
    falls there."""

    def __init__(self, ref: str | None, line=None, forward=None, back=None):
        self.ref = ref
        #: (k, c): the coordinate in the file is k·position + c — a
        #: straight side stays straight, an ellipse an ellipse
        self.line = line
        self.forward, self.back = forward, back

    def place(self, positions):
        positions = np.asarray(positions, dtype=np.float64)
        if self.line is not None:
            return positions * self.line[0] + self.line[1]
        return self.forward(positions)


class _Document:
    """The Gating-ML of one sample's tree."""

    def __init__(self, prepared: Prepared, report: Report):
        self.sample, self.tree, self.counts = prepared.sample, prepared.tree, prepared.counts
        self.report = report
        self.taken: set[str] = {MATRIX}
        self.transformations: dict[str, tuple[_Axis, ET.Element | None]] = {}
        self.gates: list[ET.Element] = []
        #: gate identifier → its `custom_info`, to put the counts in
        self.info: dict[str, ET.Element] = {}
        self.count_of: dict[str, int] = {}
        self.quadrants: dict[tuple, tuple[ET.Element, list[str]]] = {}
        self.fits: list[str] = []
        #: identifier of a Quadrant → the population it is
        self.named: dict[str, GateNode] = {}
        self.written = 0
        comp = self.sample.compensation
        self.matrix = comp if comp is not None and len(comp.channels) >= 2 else None
        if comp is not None and self.matrix is None:
            report.note(self.sample.name, "", "a matrix of one channel cannot be written in "
                                              "Gating-ML: its gates are left uncompensated")

    # --- names ---------------------------------------------------------------------

    def ident(self, name: str) -> str:
        """An XML identifier nobody has yet, as near the name as one can be."""
        base = re.sub(r"[^A-Za-z0-9_.-]", "_", name.strip()) or "gate"
        if not re.match(r"[A-Za-z_]", base):
            base = "_" + base
        out, k = base, 2
        while out in self.taken:
            out, k = f"{base}.{k}", k + 1
        self.taken.add(out)
        return out

    def describe(self, gate: ET.Element, gate_id: str, node: GateNode | None,
                 helper: str = "") -> None:
        info = ET.SubElement(ET.SubElement(gate, q("data-type", "custom_info")), "milc")
        if helper or node is None or is_helper(node):
            ET.SubElement(info, "helper").text = text(
                helper or "written to say the population that names it; not one in MIL-C")
        if node is not None and not helper:
            ET.SubElement(info, "name").text = text(node.name)
            ET.SubElement(info, "path").text = text(source_of(node))
            self.info[gate_id] = info
            if self.counts and node.path in self.counts:
                self.count_of[gate_id] = self.counts[node.path]

    # --- parameters and scales -------------------------------------------------------

    def dimension(self, parent: ET.Element, dim: str, axis: _Axis | None = None,
                  tag: str = "dimension", **attrs) -> ET.Element:
        channel, compensated = self.sample.dimension(dim)
        compensated = compensated and self.matrix is not None
        el = ET.SubElement(parent, q("gating", tag), {
            **{q("gating", k): v for k, v in attrs.items() if v is not None},
            q("gating", "compensation-ref"): MATRIX if compensated else UNCOMPENSATED})
        if axis is not None and axis.ref is not None:
            el.set(q("gating", "transformation-ref"), axis.ref)
        ET.SubElement(el, q("data-type", "fcs-dimension"), {q("data-type", "name"): channel})
        return el

    def axis(self, scale) -> _Axis:
        key = repr(to_spec(scale)) if scale is None or to_spec(scale) is not None \
            else f"{type(scale).__name__}@{id(scale)}"
        if key not in self.transformations:
            self.transformations[key] = self._axis(scale)
        return self.transformations[key][0]

    def _transformation(self, name: str, kind: str, **values) -> tuple[str, ET.Element]:
        ref = self.ident(name)
        el = ET.Element(q("transforms", "transformation"), {q("transforms", "id"): ref})
        ET.SubElement(el, q("transforms", kind),
                      {q("transforms", k): num(v) for k, v in values.items()})
        return ref, el

    def _axis(self, scale) -> tuple[_Axis, ET.Element | None]:
        line = affine(scale)
        if line is not None:
            return _Axis(None, line), None
        if isinstance(scale, Logicle):
            ref, el = self._transformation("logicle", "logicle", T=scale.T, W=scale.W,
                                           M=scale.M, A=scale.A)
            return _Axis(ref, (1.0, 0.0)), el
        if isinstance(scale, Arcsinh):
            m = float(np.arcsinh(_ASINH_TOP / scale.cofactor) / np.log(10.0))
            ref, el = self._transformation("arcsinh", "fasinh", T=_ASINH_TOP, M=m, A=0.0)
            return _Axis(ref, (1.0 / (m * np.log(10.0)), 0.0)), el
        if isinstance(scale, Log):
            ref, el = self._transformation("log", "flog", T=scale.max,
                                           M=np.log10(scale.max / scale.min))
            self.report.note(self.sample.name, "", f"values at or below {scale.min:g} sit on the "
                             "edge of a logarithmic axis here; in Gating-ML they are off it, "
                             "outside every gate drawn on that axis")
            return _Axis(ref, (1.0, 0.0)), el
        if isinstance(scale, Biex):
            fit = logicle_like(scale)
            if fit is not None:
                near, apart = fit
                ref, el = self._transformation("biex_as_logicle", "logicle", T=near.T, W=near.W,
                                               M=near.M, A=near.A)
                said = (f"FlowJo's biexponential (width {scale.width:g}, positive {scale.pos:g}, "
                        f"negative {scale.neg:g}) is written as the Logicle nearest to it "
                        f"(T={near.T:.6g}, W={near.W:.6g}, M={near.M:.6g}, A={near.A:.6g}: "
                        f"{apart:.1e} of the scale away at most, {apart * 256:.4f} of FlowJo's "
                        "256 channels), with the vertices placed to hold the same events")
                self.fits.append(said)
                self.report.note(self.sample.name, "", said)
                return _Axis(ref, None, lambda s: near(scale.inverse(s)),
                             lambda c: scale(near.inverse(c))), el
        if not (callable(scale) and hasattr(scale, "inverse")):
            raise Unwritable(f"it is drawn on a scale ({type(scale).__name__}) that cannot be "
                             "written")
        # no transformation to name: in the file's units, with vertices added
        return _Axis(None, None, scale.inverse, scale), None

    # --- gates -----------------------------------------------------------------------

    def gate(self, kind: str, gate_id: str, parent: str | None) -> ET.Element:
        el = ET.Element(q("gating", kind), {q("gating", "id"): gate_id})
        if parent:
            el.set(q("gating", "parent_id"), parent)
        return el

    def polygon(self, el: ET.Element, node: GateNode, dims, positions, scales) -> None:
        """`positions`: the vertices on `scales`, where the sides are straight."""
        axes = [self.axis(s) for s in scales]
        bent = [a.line is None for a in axes]
        if any(bent):
            moved = bend(positions, [a.forward if b else None for a, b in zip(axes, bent)],
                         [a.back if b else None for a, b in zip(axes, bent)])
            points = np.column_stack([moved[:, k] if b else a.place(moved[:, k])
                                      for k, (a, b) in enumerate(zip(axes, bent))])
        else:
            points = np.column_stack([a.place(positions[:, k]) for k, a in enumerate(axes)])
        for dim, axis in zip(dims, axes):
            self.dimension(el, dim, axis)
        for point in points:
            v = ET.SubElement(el, q("gating", "vertex"))
            for value in point:
                ET.SubElement(v, q("gating", "coordinate"), {q("data-type", "value"): num(value)})

    def drawn(self, node: GateNode, gate_id: str, parent: str | None) -> list[ET.Element]:
        """A rectangle, polygon or ellipse — and the `not` of it where the
        population is what lies outside."""
        g = node.gate
        own = gate_id if g.inside else self.ident(gate_id + ".inside")
        if isinstance(g, RectangleGate):
            el = self.gate("RectangleGate", own, parent)
            fill = lambda: [self.dimension(el, d, min=None if lo is None else num(lo),   # noqa: E731
                                           max=None if hi is None else num(hi))
                            for d, lo, hi in zip(g.dims, g.mins, g.maxs)]
        elif isinstance(g, PolygonGate):
            el = self.gate("PolygonGate", own, parent)
            scales = scales_of(g)
            vertices = np.asarray(g.vertices, dtype=np.float64)
            where = np.column_stack([on(s, vertices[:, k]) for k, s in enumerate(scales)])
            fill = lambda: self.polygon(el, node, g.dims, where, scales)   # noqa: E731
        elif isinstance(g, EllipseGate):
            scales = scales_of(g)
            axes = [self.axis(s) for s in scales]
            if all(a.line is not None for a in axes):
                el = self.gate("EllipsoidGate", own, parent)
                fill = lambda: self.ellipse(el, g, axes)   # noqa: E731
            else:
                el = self.gate("PolygonGate", own, parent)
                where = ellipse_points(g.center, g.radii, g.angle)
                self.report.note(self.sample.name, source_of(node),
                                 "an ellipse drawn on a biexponential axis is not an ellipse on "
                                 f"any scale Gating-ML has: written as a polygon of {len(where)} "
                                 "sides")
                fill = lambda: self.polygon(el, node, g.dims, where, scales)   # noqa: E731
        else:
            raise Unwritable(f"a kind of gate ({type(g).__name__}) Gating-ML does not hold")
        if g.inside:
            self.describe(el, own, node)
            fill()
            return [el]
        self.describe(el, own, None, f"the gate of {source_of(node)}, which is what lies "
                                     "outside it")
        fill()
        return [el, self.negation(node, gate_id, parent, own)]

    def ellipse(self, el: ET.Element, g: EllipseGate, axes: list[_Axis]) -> None:
        stretch = [a.line[0] for a in axes]
        mean = [a.place(c) for a, c in zip(axes, g.center)]
        cov = ellipse_covariance(g.radii, g.angle, stretch)
        for dim, axis in zip(g.dims, axes):
            self.dimension(el, dim, axis)
        point = ET.SubElement(el, q("gating", "mean"))
        for value in mean:
            ET.SubElement(point, q("gating", "coordinate"), {q("data-type", "value"): num(value)})
        matrix = ET.SubElement(el, q("gating", "covarianceMatrix"))
        for row in cov:
            r = ET.SubElement(matrix, q("gating", "row"))
            for value in row:
                ET.SubElement(r, q("gating", "entry"), {q("data-type", "value"): num(value)})
        ET.SubElement(el, q("gating", "distanceSquare"), {q("data-type", "value"): "1.0"})

    def negation(self, node: GateNode, gate_id: str, parent: str | None, of: str) -> ET.Element:
        el = self.gate("BooleanGate", gate_id, parent)
        self.describe(el, gate_id, node)
        ET.SubElement(ET.SubElement(el, q("gating", "not")), q("gating", "gateReference"),
                      {q("gating", "ref"): of})
        return el

    def quadrant(self, node: GateNode, gate_id: str, parent: str | None) -> list[ET.Element]:
        g = node.gate
        own = gate_id if g.inside else self.ident(gate_id + ".inside")
        key = (parent, g.group, tuple(g.dims), float(g.x), float(g.y)) if g.group and g.inside \
            else (own,)
        new = key not in self.quadrants
        if new:
            el = self.gate("QuadrantGate", self.ident("Quadrants"), parent)
            self.describe(el, "", None, "the corners of one quadrant gate; each is a population")
            dividers = []
            for dim, at in zip(g.dims, (g.x, g.y)):
                dividers.append(self.ident("divider"))
                d = self.dimension(el, dim, tag="divider", id=dividers[-1])
                ET.SubElement(d, q("gating", "value")).text = num(at)
            self.quadrants[key] = (el, dividers)
        el, dividers = self.quadrants[key]
        corner = ET.SubElement(el, q("gating", "Quadrant"), {q("gating", "id"): own})
        for divider, at, high in zip(dividers, (g.x, g.y), QUADRANTS[g.quadrant]):
            aside = max(1.0, abs(float(at)))
            ET.SubElement(corner, q("gating", "position"), {
                q("gating", "divider_ref"): divider,
                q("gating", "location"): num(at + aside if high else at - aside)})
        # a Quadrant has no room for a custom_info: what it is called is said here
        made = [el] if new else []
        if g.inside:
            self.named[own] = node
            return made
        return made + [self.negation(node, gate_id, parent, own)]

    def combination(self, node: GateNode, gate_id: str, parent: str | None,
                    sisters: dict[str, str]) -> list[ET.Element]:
        g = node.gate
        el = self.gate("BooleanGate", gate_id, parent)
        self.describe(el, gate_id, node)
        inner = ET.SubElement(el, q("gating", {AND: "and", OR: "or", NOT: "not"}[g.op]))
        for name in g.operands:
            ET.SubElement(inner, q("gating", "gateReference"), {q("gating", "ref"): sisters[name]})
        return [el]

    def populations(self, parent: GateNode, parent_id: str | None) -> None:
        sisters: dict[str, str] = {}

        def build(node: GateNode):
            taken = set(self.taken)
            try:
                gate_id = self.ident(node.name)
                if isinstance(node.gate, BooleanGate):
                    made = self.combination(node, gate_id, parent_id, sisters)
                elif isinstance(node.gate, QuadrantGate):
                    made = self.quadrant(node, gate_id, parent_id)
                else:
                    made = self.drawn(node, gate_id, parent_id)
            except Exception:
                self.taken = taken       # a gate that fails leaves nothing behind
                raise
            sisters[node.name] = gate_id
            return gate_id, made

        for node, (gate_id, made) in walk(parent, build, self.sample.name, self.report):
            self.gates += made
            if not is_helper(node):
                self.written += 1
            self.populations(node, gate_id)

    # --- the file ------------------------------------------------------------------------

    def build(self) -> None:
        self.populations(self.tree, None)
        # the quadrants, which cannot describe themselves, in their gate's custom_info
        for el, _ in self.quadrants.values():
            info = el.find(f"{q('data-type', 'custom_info')}/milc")
            for corner in el.findall(q("gating", "Quadrant")):
                node = self.named.get(corner.get(q("gating", "id")))
                if node is None:
                    continue
                entry = ET.SubElement(info, "quadrant", {"id": corner.get(q("gating", "id"))})
                ET.SubElement(entry, "name").text = text(node.name)
                ET.SubElement(entry, "path").text = text(source_of(node))
                self.info[corner.get(q("gating", "id"))] = entry
                if self.counts and node.path in self.counts:
                    self.count_of[corner.get(q("gating", "id"))] = self.counts[node.path]

    def spectrum(self) -> ET.Element | None:
        if self.matrix is None:
            return None
        el = ET.Element(q("transforms", "spectrumMatrix"), {
            q("transforms", "id"): MATRIX, q("transforms", "matrix-inverted-already"): "false"})
        for tag in ("fluorochromes", "detectors"):
            names = ET.SubElement(el, q("transforms", tag))
            for c in self.matrix.channels:
                ET.SubElement(names, q("data-type", "fcs-dimension"), {q("data-type", "name"): c})
        for row in np.asarray(self.matrix.matrix, dtype=np.float64):
            line = ET.SubElement(el, q("transforms", "spectrum"))
            for value in row:
                ET.SubElement(line, q("transforms", "coefficient"),
                              {q("transforms", "value"): num(value)})
        return el

    def body(self) -> list[ET.Element]:
        """Everything but the heading, in the order of the file."""
        uses_matrix = any(e.get(q("gating", "compensation-ref")) == MATRIX
                          for g in self.gates for e in g.iter())
        out = [self.spectrum()] if uses_matrix else []
        out += [el for _, el in self.transformations.values() if el is not None]
        return out + self.gates

    def plain(self) -> str:
        return "".join(ET.tostring(el, encoding="unicode") for el in self.body())


def _write(path: str, documents: list[_Document]) -> None:
    """One file for these samples, which are gated alike: the gates of the
    first, the counts of each."""
    first = documents[0]
    root = ET.Element(q("gating", "Gating-ML"), {q("xsi", "schemaLocation"): SCHEMAS})
    head = ET.SubElement(ET.SubElement(root, q("data-type", "custom_info")), "milc", {
        "version": __version__, "written": time.strftime("%Y-%m-%dT%H:%M:%S")})
    for d in documents:
        ET.SubElement(head, "sample", {"name": text(d.sample.name), "events": str(d.sample.events)})
    for said in dict.fromkeys(first.fits):
        ET.SubElement(head, "note").text = said
    for gate_id, info in first.info.items():
        for d in documents:
            if gate_id in d.count_of:
                ET.SubElement(info, "count", {"sample": text(d.sample.name)}).text = str(
                    d.count_of[gate_id])
    root.extend(first.body())
    ET.indent(root, space="  ")
    ET.ElementTree(root).write(path, encoding="UTF-8", xml_declaration=True)


def write_gatingml(samples: list[ExportSample], path: str) -> Report:
    """Write the gates of these samples: to `path` if they are all gated
    and compensated alike, else to a file for each beside it
    (`gates.xml` → `gates_<sample>.xml`). Never raises for a gate or a
    sample that cannot be written: they are in the report."""
    for prefix, uri in NS.items():
        ET.register_namespace(prefix, uri)
    report = Report(format=GATINGML)
    documents = []
    for sample in samples:
        try:
            tree = prepare(sample, report)
            document = _Document(Prepared(sample, tree, counted(sample, tree, report)), report)
            document.build()
        except Exception as exc:
            report.note(sample.name, "", f"its gates could not be written: {exc}")
            continue
        if not document.gates:
            report.note(sample.name, "", "it has no gate to write: no file for it")
            continue
        documents.append(document)
    if not documents:
        return report
    report.samples = len(documents)
    report.written = sum(d.written for d in documents)
    if all(d.plain() == documents[0].plain() for d in documents):
        _write(path, documents)
        report.files.append(path)
        return report
    stem, ext = os.path.splitext(path)
    for d in documents:
        name = re.sub(r"[^\w.-]+", "_", os.path.splitext(d.sample.name)[0])
        report.files.append(f"{stem}_{name}{ext or '.xml'}")
        _write(report.files[-1], [d])
    return report
