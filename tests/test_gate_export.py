"""The gates written for another program: a FlowJo workspace and Gating-ML.

What is asserted is membership: a population read back from what was
written holds the events it holds here. The workspace is read by MIL-C's
own reader; Gating-ML by a small reader written here from the standard —
and by FlowKit where it is installed (it is not a dependency).
"""

import os
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from milc import gate_export, gatingml, wsp, wsp_write
from milc.compensation import Compensation
from milc.fcs import write_fcs
from milc.gating import (AND, NOT, OR, QUADRANTS, BooleanGate, ClusterGate, EllipseGate, GateNode,
                         PolygonGate, QuadrantGate, RectangleGate, TimeGate, points_in_polygon,
                         wrap)
from milc.sample import Sample
from milc.session import Session
from milc.transforms import Arcsinh, Biex, Linear, Log, Logicle

CHANNELS = ["FSC-A", "SSC-A", "FITC-A", "PE-A", "APC-A", "BV421-A", "Time"]
SPILL = "4,FITC-A,PE-A,APC-A,BV421-A,1,0.08,0,0.01,0.03,1,0.02,0,0,0.05,1,0,0.02,0,0,1"
L, A, G, B = Logicle(262144, 0.8, 4.5, 0), Arcsinh(150), Log(10, 262144), Biex()
#: the populations whose shape is redrawn with more vertices: an event on
#: the very edge may fall on the other side
REDRAWN = 2


@pytest.fixture
def folder(tmp_path):
    rng = np.random.default_rng(7)
    for name, n in (("a", 3000), ("b", 2500), ("c", 4000)):
        scatter = rng.normal([90000, 50000], [25000, 15000], (n, 2))
        fluo = rng.lognormal(6.5, 1.3, (n, 4)) - 150
        time = np.sort(rng.uniform(0, 60000, n)).round()
        write_fcs(str(tmp_path / f"{name}.fcs"), np.column_stack([scatter, fluo, time]), CHANNELS,
                  markers=["", "", "CD8", "CD4", "Foxp3", "TIGIT", ""],
                  ranges=[262144.0] * 7, keywords={"SPILL": SPILL, "$TIMESTEP": "0.01"})
    return str(tmp_path)


class Curly:
    """A kind of gate the export has never heard of."""
    dims = ["FSC-A", "SSC-A"]
    inside = True

    def contains(self, column):
        return column("FSC-A") > 60000


def polygon(dims, scales, points):
    v = np.column_stack([np.asarray(p, float) if s is None else s.inverse(np.asarray(p, float))
                         for s, p in zip(scales, zip(*points))])
    return PolygonGate(dims, v, True, scales if any(s is not None for s in scales) else None)


def every_gate(session: Session) -> dict[str, int]:
    """Every kind of gate on every scale, in every sample's tree; each
    sample cleaned of its own stretches of time, the first compensated
    its own way. Returns how many events a population may differ by."""
    session.set_scale(["PE-A"], L)
    session.set_scale(["APC-A"], A)
    session.set_scale(["BV421-A"], G)
    entries = session.usable()
    m, _ = session.matrix_in_use(entries[0])
    edited = np.array(m.matrix, dtype=float)
    edited[0, 1] += 0.05
    session.set_matrix(Compensation(list(m.channels), edited, "edited"), [entries[0]])
    slack = {}
    for k, e in enumerate(entries):
        t = GateNode(e.name)
        cells = t.add(GateNode("Cells", polygon(["FSC-A", "SSC-A"], [None, None], [
            (30000, 15000), (150000, 20000), (160000, 90000), (40000, 85000)])))
        cells.add(GateNode("big", RectangleGate(["FSC-A"], [100000.0], [None])))
        cells.add(GateNode("not small", RectangleGate(["FSC-A", "SSC-A"], [0.0, 0.0],
                                                      [70000.0, 45000.0], inside=False)))
        cells.add(GateNode("linear ellipse", EllipseGate(
            ["FSC-A", "SSC-A"], [90000.0, 50000.0], [40000.0, 15000.0], 0.4)))
        cells.add(GateNode("round", EllipseGate(["FSC-A", "SSC-A"], [90000.0, 50000.0],
                                                [20000.0, 20000.0], 0.0)))
        cells.add(GateNode("on a range", EllipseGate(
            ["FSC-A", "SSC-A"], [0.35, 0.2], [0.1, 0.05], 0.3, True,
            [Linear(0, 262144), Linear(0, 262144)])))
        for q in QUADRANTS:
            cells.add(GateNode(q, QuadrantGate(["Comp-FITC-A", "Comp-PE-A"], 500.0, 700.0, q, "one")))
        cells.find("Cells/Q2").add(GateNode("biex in Q2", polygon(
            ["Comp-FITC-A", "Comp-BV421-A"], [B, G], [(0.3, 0.2), (0.9, 0.25), (0.8, 0.8), (0.35, 0.7)])))
        cells.add(GateNode("outside Q3", QuadrantGate(["Comp-FITC-A", "Comp-PE-A"], 500.0, 700.0,
                                                      "Q3", "", inside=False)))
        cells.add(GateNode("biex", polygon(["Comp-FITC-A", "FSC-A"], [B, None], [
            (0.2, 40000), (0.7, 50000), (0.75, 130000), (0.3, 120000)])))
        cells.add(GateNode("logicle", polygon(["Comp-PE-A", "Comp-FITC-A"], [L, B], [
            (0.25, 0.2), (0.7, 0.25), (0.8, 0.7), (0.3, 0.75)])))
        cells.add(GateNode("arcsinh", polygon(["Comp-APC-A", "Comp-PE-A"], [A, L], [
            (0.3, 0.2), (4.0, 0.3), (5.0, 0.7), (1.0, 0.8)])))
        cells.add(GateNode("log", polygon(["Comp-BV421-A", "FSC-A"], [G, None], [
            (0.1, 40000), (0.6, 50000), (0.7, 130000), (0.2, 120000)])))
        cells.add(GateNode("biex ellipse", EllipseGate(
            ["Comp-FITC-A", "Comp-BV421-A"], [0.5, 0.45], [0.2, 0.12], 0.5, True, [B, B])))
        cells.add(GateNode("mixed ellipse", EllipseGate(
            ["FSC-A", "Comp-FITC-A"], [90000.0, 0.5], [35000.0, 0.15], 0.0, True, [None, B])))
        cells.add(GateNode("logicle ellipse", EllipseGate(
            ["Comp-PE-A", "Comp-APC-A"], [0.5, 2.5], [0.2, 1.5], 0.2, True, [L, A])))
        cells.add(GateNode("outside an ellipse", EllipseGate(
            ["Comp-FITC-A", "Comp-BV421-A"], [0.5, 0.45], [0.2, 0.12], 0.5, False, [B, B])))
        cells.add(GateNode("Q1 or Q3", BooleanGate(OR, ["Q1", "Q3"]))).add(
            GateNode("under a combination", RectangleGate(["SSC-A"], [45000.0], [None])))
        cells.add(GateNode("both", BooleanGate(AND, ["biex", "linear ellipse"])))
        cells.add(GateNode("not Q2", BooleanGate(NOT, ["Q2"])))
        cells.add(GateNode("neither", BooleanGate(NOT, ["biex", "logicle"])))
        cells.add(GateNode("just one", BooleanGate(AND, ["big"])))
        cells.add(GateNode("of a combination", BooleanGate(AND, ["Q1 or Q3", "not small"])))
        time = Sample.open(e.path).column("Time")
        at = lambda a, b: [float(np.quantile(time, a)), float(np.quantile(time, b))]  # noqa: E731
        wrap(t, "Clean", (TimeGate(["Time"], [at(0.1, 0.2), at(0.5, 0.55), at(0.9, 0.95)], False)
                          if k == 0 else TimeGate(["Time"], [at(0.05, 0.6)], True) if k == 1
                          else TimeGate(["Time"], [], False)))
        e.gates, e.sample = t, None
    session._counts.clear()
    for name in ("arcsinh", "log", "logicle ellipse", "biex in Q2", "biex", "logicle",
                 "biex ellipse", "mixed ellipse", "outside an ellipse", "both", "neither"):
        slack["Clean/Cells/" + name] = REDRAWN
    return slack


def read_back(session, path):
    """(sample name → masks as gated here, masks as gated from the
    workspace read back), and the workspace."""
    w = wsp.read_wsp(path)
    out = {}
    for e in session.usable():
        back = w.sample(os.path.basename(e.path))
        matrix = session.matrix_in_use(e)[0]
        assert back.compensation.channels == matrix.channels
        assert np.array_equal(back.compensation.matrix, matrix.matrix)
        again = Sample.open(e.path, compensation=back.compensation, gates=back.gates)
        assert wsp.count_mismatches(back.gates, again.count, tolerance=0, slack=0) == []
        out[e.name] = (session.load(e).masks(), again.masks())
    return out, w


def test_every_gate_comes_back_from_a_workspace(folder, tmp_path):
    s = Session()
    s.open_folder(folder)
    slack = every_gate(s)
    out = str(tmp_path / "out.wsp")
    report = gate_export.export_gates(s, out)
    assert report.files == [out] and report.samples == 3 and report.skipped == []
    assert report.written == sum(1 for e in s.usable() for n in e.gates.walk() if n.gate is not None)
    masks, w = read_back(s, out)
    assert w.skipped == [] and [x.filename for x in w.samples] == ["a.fcs", "b.fcs", "c.fcs"]
    for name, (own, again) in masks.items():
        assert len(own) > 25
        for path, mask in own.items():
            apart = int((mask != again[path]).sum())
            assert apart <= slack.get(path, 0), (name, path, apart)
        # every drawn gate holds some events and leaves some out: the test means something
        inner = [p for p in own if p.startswith("Clean/Cells/")]
        assert all(0 < own[p].sum() < own["Clean/Cells"].sum() for p in inner
                   if p.split("/")[-1] not in ("Q1", "just one")), name
    # what was added to say a population is beside it, and nothing else is new
    extra = {p for p in masks["a"][1] if p not in masks["a"][0]}
    assert extra == {"Clean [1]", "Clean [2]", "Clean [3]", "Clean [4]", "Clean/Cells/neither [any]"}
    assert {p for p in masks["b"][1] if p not in masks["b"][0]} == {"Clean/Cells/neither [any]"}
    said = report.text()
    assert "4 ranges beside it" in said and "negation of their union" in said
    assert "polygon of 512 sides" in said and "vertices were added" in said


def test_a_workspace_is_written_as_flowjo_writes_it(folder, tmp_path):
    s = Session()
    s.open_folder(folder)
    every_gate(s)
    out = str(tmp_path / "out.wsp")
    gate_export.export_gates(s, out)
    root = ET.parse(out).getroot()
    g, t, d = (f"{{{wsp_write.NS[k]}}}" for k in ("gating", "transforms", "data-type"))
    assert root.tag == "Workspace" and root.get("flowJoVersion").startswith("10.")
    assert [c.tag for c in root] == ["Matrices", "Groups", "SampleList"]
    assert len(root.find("Matrices")) == 2          # the edited one and the files' own
    refs = root.findall("Groups/GroupNode/Group/SampleRefs/SampleRef")
    assert [r.get("sampleID") for r in refs] == ["1", "2", "3"]
    first = root.find("SampleList/Sample")
    assert [c.tag.split("}")[-1] for c in first] == ["DataSet", "spilloverMatrix", "Transformations",
                                                    "Keywords", "SampleNode"]
    assert first.find("DataSet").get("uri").startswith("file:/") and \
        first.find("DataSet").get("uri").endswith("/a.fcs")
    kinds = {x.find(d + "parameter").get(d + "name"): x.tag.split("}")[-1]
             for x in first.find("Transformations")}
    assert kinds == {"FSC-A": "linear", "SSC-A": "linear", "FITC-A": "biex", "PE-A": "logicle",
                     "APC-A": "fasinh", "BV421-A": "log", "Time": "linear",
                     "Comp-FITC-A": "biex", "Comp-PE-A": "logicle", "Comp-APC-A": "fasinh",
                     "Comp-BV421-A": "log"}
    keys = {k.get("name"): k.get("value") for k in first.find("Keywords")}
    assert keys["$TOT"] == "3000" and keys["$P3N"] == "FITC-A" and keys["SPILL"] == SPILL
    node = first.find("SampleNode")
    assert node.get("name") == "a.fcs" and node.get("count") == "3000" and node.get("sampleID") == "1"

    def pop(path, start=node):
        for name in path.split("/"):
            start = next(p for p in start.find("Subpopulations") if p.get("name") == name)
        return start

    # Time in seconds: the file's values times $TIMESTEP, as FlowJo keeps it
    time = Sample.open(os.path.join(folder, "a.fcs")).column("Time")
    lo, hi = (float(v) for v in next(iter(pop("Clean [2]").find("Gate")))[0].attrib.values())
    assert np.quantile(time, 0.2) * 0.01 < lo < hi <= np.quantile(time, 0.5) * 0.01
    assert pop("Clean").tag == "OrNode" and [x.get("name") for x in pop("Clean").find("Dependents")] == [
        "Clean [1]", "Clean [2]", "Clean [3]", "Clean [4]"]
    span = next(x for x in first.find("Transformations") if x.find(d + "parameter").get(d + "name") == "Time")
    assert float(span.get(t + "maxRange")) == pytest.approx(2621.44)

    # a quadrant is a rectangle open on two sides; the gates say who their parent is
    q1 = pop("Clean/Cells/Q1").find("Gate")
    assert [dict(x.attrib) for x in q1[0]] == [{g + "max": "500.0"}, {g + "min": "700.0"}]
    assert q1.get(g + "parent_id") == pop("Clean/Cells").find("Gate").get(g + "id")
    assert pop("Clean/Cells/not small").find("Gate")[0].get("eventsInside") == "0"
    assert pop("Clean/Cells/not Q2").tag == "NotNode"
    assert pop("Clean/Cells/not Q2").find("Dependents")[0].get("name") == "Clean/Cells/Q2"
    # an ellipse: two foci, four points of the edge, on 256 channels
    ell = pop("Clean/Cells/linear ellipse").find("Gate")[0]
    assert ell.tag == g + "EllipsoidGate"
    assert len(ell.find(g + "foci")) == 2 and len(ell.find(g + "edge")) == 4
    points = np.array([[float(c.get(d + "value")) for c in v] for v in ell.find(g + "edge")])
    assert points.mean(axis=0) == pytest.approx([90000 / 1024, 50000 / 1024])
    assert float(ell.get(g + "distance")) == pytest.approx(np.hypot(*(points[0] - points[1])))
    # a polygon on its own scales keeps its vertices; on one FlowJo does not test on, it gets more
    assert len(pop("Clean/Cells/logicle").find("Gate")[0].findall(g + "vertex")) == 4
    assert len(pop("Clean/Cells/arcsinh").find("Gate")[0].findall(g + "vertex")) > 50
    assert pop("Clean/Cells/logicle ellipse").find("Gate")[0].tag == g + "PolygonGate"
    # a plot to open on each population: the axes of what is gated in it
    assert [x.get("name") for x in pop("Clean/Cells").find("Graph")] == ["FSC-A", ""]
    assert [x.get("name") for x in pop("Clean/Cells/Q2").find("Graph")] == ["Comp-FITC-A", "Comp-BV421-A"]


def test_gates_every_sample_shares_are_the_groups_too(folder, tmp_path):
    s = Session()
    s.open_folder(folder)
    for e in s.usable():
        cells = e.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [60000.0], [None])))
        cells.add(GateNode("CD4", RectangleGate(["Comp-PE-A", "Comp-FITC-A"], [300.0, None],
                                                [20000.0, 900.0])))
        e.sample = None
    out = str(tmp_path / "shared.wsp")
    report = gate_export.export_gates(s, out)
    assert report.written == 6 and not report.skipped and not report.notes
    masks, w = read_back(s, out)
    assert all(np.array_equal(own[p], again[p]) for own, again in masks.values() for p in own)
    # as a workspace of FlowJo's: one tree and one matrix, for a file it does not list too
    assert list(w.groups) == ["All Samples"] and w.shared() is not None
    assert [n.path for n in w.groups["All Samples"].walk()] == ["", "Cells", "Cells/CD4"]
    root = ET.parse(out).getroot()
    assert {p.get("owningGroup") for p in root.iter("Population")} == {"All Samples"}
    # one sample gated another way: nothing is the group's any more
    s.usable()[1].gates.children[0].gate.mins[0] = 65000.0
    gate_export.export_gates(s, out)
    w = wsp.read_wsp(out)
    assert w.groups == {} and w.shared() is None
    assert {p.get("owningGroup") for p in ET.parse(out).getroot().iter("Population")} == {""}
    # opened by the session as any workspace is
    again = Session()
    again.open_folder(folder)
    again.set_workspace(out)
    assert [e.gates.find("Cells").gate.mins[0] for e in again.usable()] == [60000.0, 65000.0, 60000.0]
    assert again.workspace_report() == ""


def test_what_cannot_be_written_is_left_out_and_said(folder, tmp_path):
    s = Session()
    s.open_folder(folder)
    for e in s.usable():
        cells = e.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [60000.0], [None])))
        cells.add(GateNode("cluster", ClusterGate(["Comp-FITC-A"], [B], [1.0], np.zeros((2, 1)), [0])))
        curly = cells.add(GateNode("curly", Curly()))
        curly.add(GateNode("under", RectangleGate(["FSC-A"], [0.0], [1e5]))).add(
            GateNode("deeper", RectangleGate(["FSC-A"], [0.0], [9e4])))
        cells.add(GateNode("kept", RectangleGate(["SSC-A"], [None], [50000.0])))
        cells.add(GateNode("of the curly", BooleanGate(OR, ["curly", "kept"])))
        cells.add(GateNode("of a later one", BooleanGate(AND, ["kept", "last"])))
        cells.add(GateNode("last", BooleanGate(NOT, ["kept"])))
        cells.add(GateNode("no such", RectangleGate(["Comp-Nope-A"], [0.0], [1.0])))
        cells.add(GateNode("broken", PolygonGate(["FSC-A", "SSC-A"], np.array(
            [[0.0, 0.0], [1.0, np.nan], [2.0, 2.0]]))))
        cells.add(GateNode("kept", RectangleGate(["SSC-A"], [None], [1.0])))
        e.sample = None
    s.usable()[2].gates = GateNode("c")            # one sample with no gates at all
    for fmt, name in (("wsp", "out.wsp"), ("gatingml", "out.xml")):
        report = gate_export.export_gates(s, str(tmp_path / name), fmt, counts=False)
        lost = {(k.path, k.kind): k for k in report.skipped if k.sample == "a.fcs"}
        assert set(lost) == {("Cells/cluster", "ClusterGate"), ("Cells/curly", "Curly"),
                             ("Cells/of the curly", "BooleanGate"), ("Cells/of a later one", "BooleanGate"),
                             ("Cells/no such", "RectangleGate"), ("Cells/broken", "PolygonGate"),
                             ("Cells/kept", "RectangleGate")}
        assert lost[("Cells/curly", "Curly")].below == 2
        assert "does not know" in lost[("Cells/curly", "Curly")].reason
        assert "no other program has the map" in lost[("Cells/cluster", "ClusterGate")].reason
        assert "made of curly" in lost[("Cells/of the curly", "BooleanGate")].reason
        assert "no channel Comp-Nope-A" in lost[("Cells/no such", "RectangleGate")].reason
        assert "not a number" in lost[("Cells/broken", "PolygonGate")].reason
        assert "has its name" in lost[("Cells/kept", "RectangleGate")].reason
        assert report.written == 2 * 3 and "(and 2 under it)" in report.text()
        assert "18 not written" in report.summary()
    # the rest is there, and reads
    w = wsp.read_wsp(str(tmp_path / "out.wsp"))
    assert [n.path for n in w.sample("a.fcs").gates.walk()] == ["", "Cells", "Cells/kept", "Cells/last"]
    assert w.sample("c.fcs").gates.children == [] and len(w.samples) == 3
    # Gating-ML has nothing to say of a sample without gates
    assert sorted(os.listdir(tmp_path))[-1] == "out.xml" and "no gate to write" in report.text()
    assert report.samples == 2

    # a file that cannot be opened for its counts costs the counts, not the export
    def cannot():
        raise OSError("gone")
    samples = gate_export.samples_of(s)
    samples[0].load = cannot
    report = wsp_write.write_wsp(samples, str(tmp_path / "again.wsp"))
    assert report.samples == 3 and "the counts are not written: gone" in report.text()
    again = wsp.read_wsp(str(tmp_path / "again.wsp"))
    assert again.sample("a.fcs").gates.find("Cells").reference_count is None


# --- Gating-ML, read here from the standard ------------------------------------------

def read_gatingml(path: str, sample: Sample) -> dict[str, np.ndarray]:
    """MIL-C path → mask, by the rules of Gating-ML 2.0 and nothing of
    the writer's: the transformations by their formulas, a gate under its
    parent, a quadrant by its dividers."""
    ns = gatingml.NS
    g, t, d = (f"{{{ns[k]}}}" for k in ("gating", "transforms", "data-type"))
    root = ET.parse(path).getroot()
    scales = {}
    for el in root.findall(t + "transformation"):
        kind = el[0]
        p = {k.split("}")[-1]: float(v) for k, v in kind.attrib.items()}
        ln10 = np.log(10.0)
        if kind.tag == t + "logicle":
            scales[el.get(t + "id")] = Logicle(p["T"], p["W"], p["M"], p["A"])
        elif kind.tag == t + "fasinh":
            scales[el.get(t + "id")] = lambda x, p=p: (
                np.arcsinh(x * np.sinh(p["M"] * ln10) / p["T"]) + p["A"] * ln10) / ((p["M"] + p["A"]) * ln10)
        elif kind.tag == t + "flog":
            def flog(x, p=p):
                with np.errstate(invalid="ignore", divide="ignore"):
                    return np.log10(x / p["T"]) / p["M"] + 1
            scales[el.get(t + "id")] = flog
    matrix = root.find(t + "spectrumMatrix")
    if matrix is not None:
        names = [x.get(d + "name") for x in matrix.find(t + "detectors")]
        assert names == [x.get(d + "name") for x in matrix.find(t + "fluorochromes")]
        values = np.array([[float(c.get(t + "value")) for c in row]
                           for row in matrix.findall(t + "spectrum")])
        assert names == sample.compensation.channels and np.array_equal(values, sample.compensation.matrix)

    def column(dim):
        name = dim.find(d + "fcs-dimension").get(d + "name")
        ref = dim.get(g + "compensation-ref")
        assert ref in ("uncompensated", matrix.get(t + "id") if matrix is not None else None)
        x = sample.column(name if ref == "uncompensated" else "Comp-" + name).astype(np.float64)
        scale = dim.get(g + "transformation-ref")
        return scales[scale](x) if scale else x

    gates = {el.get(g + "id"): el for el in root if el.tag.startswith(g)}
    corners = {q.get(g + "id"): (el, q) for el in gates.values() for q in el.findall(g + "Quadrant")}
    done: dict[str, np.ndarray] = {}

    def events(gate_id):
        if gate_id in done:
            return done[gate_id]
        if gate_id in corners:
            el, corner = corners[gate_id]
            mask = np.ones(len(sample), dtype=bool)
            dividers = {x.get(g + "id"): x for x in el.findall(g + "divider")}
            for position in corner:
                divider = dividers[position.get(g + "divider_ref")]
                at, x = float(divider.find(g + "value").text), column(divider)
                mask &= (x >= at) if float(position.get(g + "location")) >= at else (x < at)
        else:
            el = gates[gate_id]
            kind = el.tag.split("}")[-1]
            dims = el.findall(g + "dimension")
            if kind == "RectangleGate":
                mask = np.ones(len(sample), dtype=bool)
                for dim in dims:
                    x = column(dim)
                    if dim.get(g + "min") is not None:
                        mask &= x >= np.float32(dim.get(g + "min"))
                    if dim.get(g + "max") is not None:
                        mask &= x < np.float32(dim.get(g + "max"))
            elif kind == "PolygonGate":
                poly = np.array([[float(c.get(d + "value")) for c in v] for v in el.findall(g + "vertex")])
                mask = points_in_polygon(column(dims[0]), column(dims[1]), poly)
            elif kind == "EllipsoidGate":
                mean = np.array([float(c.get(d + "value")) for c in el.find(g + "mean")])
                cov = np.array([[float(c.get(d + "value")) for c in row]
                                for row in el.find(g + "covarianceMatrix")])
                off = np.column_stack([column(dims[0]), column(dims[1])]) - mean
                mask = np.einsum("ij,jk,ik->i", off, np.linalg.inv(cov), off) <= float(
                    el.find(g + "distanceSquare").get(d + "value"))
            else:
                op = next(x for x in el if x.tag in (g + "and", g + "or", g + "not"))
                parts = [events(r.get(g + "ref")) for r in op]
                mask = (np.logical_and.reduce(parts) if op.tag == g + "and"
                        else np.logical_or.reduce(parts) if op.tag == g + "or" else ~parts[0])
        parent = el.get(g + "parent_id")
        done[gate_id] = mask & events(parent) if parent else mask
        return done[gate_id]

    out = {}
    for gate_id, el in gates.items():
        info = el.find(d + "custom_info/milc")
        if info.find("helper") is None:
            out[info.findtext("path")] = events(gate_id)
        for corner in info.findall("quadrant"):
            out[corner.findtext("path")] = events(corner.get("id"))
    return out


def test_every_gate_comes_back_from_gatingml(folder, tmp_path):
    s = Session()
    s.open_folder(folder)
    slack = every_gate(s)
    out = str(tmp_path / "gates.xml")
    report = gate_export.export_gates(s, out)
    # each sample has its own stretches of time (and one its own matrix): a file for each
    assert [os.path.basename(f) for f in report.files] == ["gates_a.xml", "gates_b.xml", "gates_c.xml"]
    assert report.skipped == [] and report.samples == 3
    for e, path in zip(s.usable(), report.files):
        own = s.load(e).masks()
        again = read_gatingml(path, s.load(e))
        assert set(again) == set(own) - {""}
        for p, mask in again.items():
            assert int((mask != own[p]).sum()) <= slack.get(p, 0), (e.name, p)
    said = report.text()
    assert "written as the Logicle nearest to it" in said and "of the scale away" in said
    assert "edge of a logarithmic axis" in said and "polygon of 512 sides" in said
    root = ET.parse(report.files[0]).getroot()
    g, t, d = (f"{{{gatingml.NS[k]}}}" for k in ("gating", "transforms", "data-type"))
    kinds = sorted({el[0].tag.split("}")[-1] for el in root.findall(t + "transformation")})
    assert kinds == ["fasinh", "flog", "logicle"]
    quads = root.findall(g + "QuadrantGate")
    assert len(quads) == 2 and len(quads[0].findall(g + "Quadrant")) == 4
    # an ellipse stays one where its axes allow; the counts are beside the names
    ellipses = {el.find(d + "custom_info/milc").findtext("name") for el in root.findall(g + "EllipsoidGate")}
    assert ellipses == {"linear ellipse", "round", "on a range", "logicle ellipse"}
    cells = next(el for el in root.findall(g + "PolygonGate") if el.get(g + "id") == "Cells")
    assert cells.get(g + "parent_id") == "Clean"
    assert cells.find(d + "custom_info/milc/count").text == str(s.load(s.usable()[0]).count("Clean/Cells"))
    # identifiers are XML names, unique in the file
    ids = [el.get(g + "id") for el in root.iter() if el.get(g + "id")] + [
        el.get(t + "id") for el in root.iter() if el.get(t + "id")]
    assert len(ids) == len(set(ids)) and all(i[0].isalpha() or i[0] == "_" for i in ids)
    assert not any(" " in i for i in ids)


def test_one_gatingml_file_where_the_samples_are_gated_alike(folder, tmp_path):
    s = Session()
    s.open_folder(folder)
    for e in s.usable():
        cells = e.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [60000.0], [None])))
        cells.add(GateNode("Single Cells", polygon(["Comp-PE-A", "Comp-FITC-A"], [B, B], [
            (0.3, 0.1), (0.9, 0.1), (0.9, 0.6), (0.3, 0.6)]))).add(
            GateNode("Single Cells", RectangleGate(["Comp-APC-A"], [100.0], [5000.0])))
        e.sample = None
    out = str(tmp_path / "gates.xml")
    report = gate_export.export_gates(s, out)
    assert report.files == [out] and report.written == 9
    root = ET.parse(out).getroot()
    g, d = (f"{{{gatingml.NS[k]}}}" for k in ("gating", "data-type"))
    assert [x.get("name") for x in root.iter("sample")] == ["a.fcs", "b.fcs", "c.fcs"]
    # a name twice in the tree is two identifiers; each sample's count is with it
    assert [el.get(g + "id") for el in root if el.get(g + "id")] == [
        "Cells", "Single_Cells", "Single_Cells.2"]
    counts = root.find(g + "RectangleGate").findall(d + "custom_info/milc/count")
    assert [(c.get("sample"), int(c.text)) for c in counts] == [
        (os.path.basename(e.path), s.load(e).count("Cells")) for e in s.usable()]
    for e in s.usable():
        own, again = s.load(e).masks(), read_gatingml(out, s.load(e))
        assert all(np.array_equal(again[p], own[p]) for p in again) and len(again) == 3
    # compensated another way, a sample needs a file of its own
    m, _ = s.matrix_in_use(s.usable()[0])
    s.set_matrix(Compensation(list(m.channels), m.matrix * 1.01, "other"), [s.usable()[0]])
    assert len(gate_export.export_gates(s, out).files) == 3


def test_gatingml_is_read_by_flowkit(folder, tmp_path):
    fk = pytest.importorskip("flowkit")
    s = Session()
    s.open_folder(folder)
    every_gate(s)
    report = gate_export.export_gates(s, str(tmp_path / "gates.xml"))
    for e, path in zip(s.usable(), report.files):
        own = s.load(e).masks()
        strategy = fk.parse_gating_xml(path)            # validates against the schema
        sample = fk.Sample(e.path, sample_id=e.name, preprocess=False)
        session = fk.Session(gating_strategy=strategy, fcs_samples=[sample])
        session.analyze_samples(sample_id=e.name, use_mp=False)
        again = read_gatingml(path, s.load(e))
        root = ET.parse(path).getroot()
        d = f"{{{gatingml.NS['data-type']}}}"
        named = {}
        for el in root:
            info = el.find(d + "custom_info/milc")
            if info is not None and info.find("helper") is None and info.find("path") is not None:
                named[el.get(f"{{{gatingml.NS['gating']}}}id")] = info.findtext("path")
            for corner in (info.findall("quadrant") if info is not None else []):
                named[corner.get("id")] = corner.findtext("path")
        assert set(named.values()) == set(again)
        for gate_id, gate_path in strategy.get_gate_ids():
            if gate_id in named:
                theirs = session.get_gate_membership(e.name, gate_id, gate_path)
                assert int((theirs != own[named[gate_id]]).sum()) <= REDRAWN, named[gate_id]


def test_the_workspace_is_read_by_flowkit(folder, tmp_path):
    fk = pytest.importorskip("flowkit")
    s = Session()
    s.open_folder(folder)
    every_gate(s)
    out = str(tmp_path / "gates.wsp")
    gate_export.export_gates(s, out)
    w = fk.Workspace(out, find_fcs_files_from_wsp=True)      # by the URIs written
    assert sorted(w.get_sample_ids()) == ["a.fcs", "b.fcs", "c.fcs"]
    # FlowKit makes a polygon of 128 sides of FlowJo's ellipse, inside it: it
    # loses the events of the sliver between
    ellipses = ("ellipse", "on a range", "round", "both", "of a combination")
    for e in s.usable():
        own = s.load(e).masks()
        name = os.path.basename(e.path)
        w.analyze_samples(sample_id=name, use_mp=False)
        for gate_id, gate_path in w.get_gate_ids(name):
            path = "/".join(gate_path[1:] + (gate_id,))
            if path not in own:
                continue
            theirs = w.get_gate_membership(name, gate_id, gate_path)
            apart = int((theirs != own[path]).sum())
            allowed = max(REDRAWN, 0.01 * own[path].sum()) if any(x in gate_id for x in ellipses) \
                else REDRAWN
            assert apart <= allowed, (name, path, apart)


# --- the pieces ------------------------------------------------------------------------

@pytest.mark.parametrize("biex", [Biex(), Biex(width=-10), Biex(width=-800), Biex(neg=1.0),
                                  Biex(256, 10000, 0.5, -30, 4.0)])
def test_the_logicle_nearest_a_biexponential(biex):
    near, apart = gate_export.logicle_like(biex)
    assert isinstance(near, Logicle) and apart < 1e-5
    x = biex.inverse(np.linspace(0, 1, 4001))
    assert np.abs(near(x) - biex(x)).max() <= apart * 1.5


def test_a_side_redrawn_on_another_scale_stays_where_it_was():
    # straight on an arcsinh axis, written in the file's units: the pieces are within the tolerance
    drawn = np.array([[0.3, 100.0], [5.0, 200.0], [6.0, 900.0], [1.0, 800.0]])
    out = gate_export.bend(drawn, [A.inverse, None], [A, None])
    assert len(out) > 40 and np.allclose(out[0], [A.inverse(0.3), 100.0])
    back = np.column_stack([A(out[:, 0]), out[:, 1]])
    x = A(np.random.default_rng(0).uniform(A.inverse(0.0), A.inverse(6.5), 200000))
    y = np.random.default_rng(1).uniform(0, 1000, 200000)
    apart = points_in_polygon(x, y, drawn) != points_in_polygon(A.inverse(x), y, out)
    assert apart.sum() <= 2 and len(back) == len(out)
    # nothing to do where the axes only stretch
    assert np.array_equal(gate_export.bend(drawn, [None, None], [None, None]), drawn)


def test_time_ranges_hold_the_events_of_the_intervals():
    t = np.sort(np.random.default_rng(3).uniform(0, 50000, 5000)).astype(np.float32)
    column = lambda _: t   # noqa: E731
    for gate in (TimeGate(["Time"], [[float(t[100]), float(t[200])], [float(t[900]), 30000.4]], True),
                 TimeGate(["Time"], [[float(t[100]), float(t[200])], [float(t[150]), float(t[400])],
                                     [40000.7, 60000.0]], False),
                 TimeGate(["Time"], [], False), TimeGate(["Time"], [[5.0, 1.0]], True)):
        ranges = gate_export.time_ranges(gate)
        mask = np.zeros(len(t), dtype=bool)
        for lo, hi in ranges:
            mask |= RectangleGate(["Time"], [lo], [hi]).contains(column)
            # in seconds and back, as a workspace keeps them
            assert RectangleGate(["Time"], [None if lo is None else lo * 0.001 / 0.001],
                                 [None if hi is None else hi * 0.001 / 0.001]).contains(column).sum() == \
                RectangleGate(["Time"], [lo], [hi]).contains(column).sum()
        assert np.array_equal(mask, gate.contains(column))


def test_a_workspace_of_flowjo_with_a_gate_on_time_is_read_in_the_files_units(tmp_path):
    # FlowJo keeps Time in seconds: a range of 2 to 5 s is 200 to 500 of a file with a step of 0.01
    text = """<?xml version="1.0"?>
<Workspace flowJoVersion="10.8.1" xmlns:g="http://g" xmlns:t="http://t" xmlns:d="http://d">
 <Groups/><SampleList><Sample><DataSet uri="file:/C:/x/a.fcs" sampleID="1"/>
 <Transformations><t:linear t:minRange="0.5" t:maxRange="600" gain="1"><d:parameter d:name="Time"/></t:linear>
 </Transformations>
 <Keywords><Keyword name="$P1N" value="FSC-A"/><Keyword name="$P2N" value="Time"/>
  <Keyword name="$TIMESTEP" value="0.01"/></Keywords>
 <SampleNode name="a.fcs" count="10"><Subpopulations>
  <Population name="steady" count="3"><Gate><g:RectangleGate eventsInside="1">
   <g:dimension g:min="2" g:max="5"><d:fcs-dimension d:name="Time"/></g:dimension>
   <g:dimension g:min="10"><d:fcs-dimension d:name="FSC-A"/></g:dimension>
  </g:RectangleGate></Gate></Population></Subpopulations></SampleNode>
 </Sample></SampleList></Workspace>"""
    path = tmp_path / "t.wsp"
    path.write_text(text)
    s = wsp.read_wsp(str(path)).samples[0]
    gate = s.gates.find("steady").gate
    assert gate.mins == [pytest.approx(200.0), 10.0] and gate.maxs == [pytest.approx(500.0), None]
    span = s.transforms["Time"]
    assert span.min == pytest.approx(50.0) and span.max == pytest.approx(60000.0)


# --- the command line ---------------------------------------------------------------------

def test_cli_writes_the_gates(folder, tmp_path, capsys):
    from milc.cli import main
    assert main(["gates", folder, "--out", str(tmp_path / "none.wsp")]) == 2
    assert "no sample here has gates" in capsys.readouterr().err
    s = Session()
    s.open_folder(folder)
    e = s.usable()[0]
    e.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [60000.0], [None]))).add(
        GateNode("cluster", ClusterGate(["Comp-FITC-A"], [B], [1.0], np.zeros((2, 1)), [0])))
    s.gates_changed(e)                              # saved beside the files, for every sample
    out = str(tmp_path / "cli.wsp")
    assert main(["gates", folder, "--out", out]) == 0
    said = capsys.readouterr().out
    assert "Wrote cli.wsp: 3 population(s) of 3 sample(s); 3 not written" in said
    assert "Cells/cluster (3 samples): a cluster is" in said and "not been opened in FlowJo" in said
    w = wsp.read_wsp(out)
    counts = [x.gates.find("Cells").reference_count for x in w.samples]
    assert counts == [s.load(x).count("Cells") for x in s.usable()] and min(counts) > 100
    # Gating-ML by the name of the file, or by --format; --qc puts Clean above the gates
    xml = str(tmp_path / "cli.xml")
    assert main(["gates", folder, "--out", xml, "--no-counts"]) == 0
    assert "Wrote cli.xml" in capsys.readouterr().out and os.path.exists(xml)
    assert ET.parse(xml).getroot().find(".//count") is None
    odd = str(tmp_path / "odd.txt")
    assert main(["gates", folder, "--out", odd, "--format", "wsp", "--qc"]) == 0
    assert wsp.read_wsp(odd).samples[0].gates.children[0].name == "Clean"


def test_a_quadrant_with_offset_arms_goes_out_as_four_polygons(folder, tmp_path):
    """Neither format has such a quadrant: each corner is written as the
    polygon its two arms cut out, and holds the same events."""
    from milc.gating import ARMS, OffsetQuadrantGate
    s = Session()
    s.open_folder(folder)
    for e in s.usable():
        # on scatter, as stored; and on two fluorescence axes, biexponential
        for dims, centre, scales, tilt in (
                (["FSC-A", "SSC-A"], [90000.0, 50000.0], None, 30000.0),
                (["Comp-FITC-A", "Comp-PE-A"], [600.0, 500.0], [B, B], 0.15)):
            sc = scales or [None, None]
            at = [float(c) if t is None else float(t(np.float64(c))) for c, t in zip(centre, sc)]
            step = {"right": (1, 0.3), "up": (-0.2, 1), "left": (-1, 0.1), "down": (0.3, -1)}
            arms = {}
            for arm in ARMS:
                p = [at[0] + step[arm][0] * tilt, at[1] + step[arm][1] * tilt]
                arms[arm] = [float(v) if t is None else float(t.inverse(np.float64(v)))
                             for v, t in zip(p, sc)]
            for q in QUADRANTS:
                e.gates.add(GateNode(f"{dims[0][:4]} {q}",
                                     OffsetQuadrantGate(dims, centre, arms, q, "g", True, scales)))
        e.edited = True
    for path, write in ((str(tmp_path / "out.wsp"), None), (str(tmp_path / "out.xml"), None)):
        report = gate_export.export_gates(s, path)
        assert report.skipped == [] and "written as a polygon" in report.text()
    masks, w = read_back(s, str(tmp_path / "out.wsp"))
    assert w.skipped == []
    for name, (own, again) in masks.items():
        for prefix in ("FSC-", "Comp"):
            corners = [p for p in own if p.startswith(prefix)]
            assert len(corners) == 4 and all(own[p].sum() > 20 for p in corners), (name, prefix)
            for p in corners:
                assert int((own[p] != again[p]).sum()) <= REDRAWN, (name, p)
