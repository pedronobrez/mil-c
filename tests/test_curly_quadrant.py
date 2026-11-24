"""
FlowJo's curly quadrant (`gating.CurlyQuadrantGate`): read from a
workspace, drawn, saved, written out.

The workspace holds a curly corner as it holds a quadrant's — the crossing
point and which corner — and no curve: whoever reads it draws the curve
(milc/wsp.py, `_read_curly`). The curve here is CytoML's, and what is
checked is that it is CytoML's: the events of every corner against a
port of cytolib's own construction (CurlyQuadGate::interpolate), and the
counts CytoML gave on the workspace `curly_files` writes, measured once
with `tests/reference/curly_compare.py`. Against FlowJo nothing has been
measured: no workspace with curly quadrants and its files was at hand.
"""

import os
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from milc import gate_export, gatingml, wsp, wsp_write
from milc.compensation import Compensation
from milc.fcs import read_header, write_fcs
from milc.gating import (CURL_BEND, CURL_CHANNELS, CURL_POWER, QUADRANTS, CurlyQuadrantGate,
                         GateNode, PolygonGate, QuadrantGate, copy_tree, follow_scale,
                         points_in_polygon, tree_from_dict, tree_to_dict)
from milc.sample import Sample
from milc.transforms import Arcsinh, Biex, Linear

# --- a file with spread, and curly quadrants on it -------------------------------

CHANNELS = ["FSC-A", "SSC-A", "FITC-A", "PE-A", "BV786-A"]
FLUO = ["FITC-A", "PE-A", "BV786-A"]
#: row = fluorochrome, column = detector, in the order of FLUO
SPILL = np.array([[1.0, 0.22, 0.004], [0.012, 1.0, 0.06], [0.0, 0.03, 1.0]])
#: the scales FlowJo gave PE and BV786 in the FlowJo 10.10 workspace with
#: curly quadrants that allotropy publishes (see `test_flowjo_writes_...`)
PE_SCALE = Biex(256, 262144, 0, -15.8489319246, 4.418539922)
BV_SCALE = Biex(256, 262144, 0, -63.0957344480193, 4.418539922)
FITC_SCALE = Linear(0.0, 262144.0)
#: the crossings: (name, dims, centre, scales, parent)
CROSSINGS = [
    ("CD8 CD4", ["Comp-PE-A", "Comp-BV786-A"], [292.6491374776836, 434.6014486060488],
     [PE_SCALE, BV_SCALE], ""),
    ("FITC PE", ["Comp-FITC-A", "Comp-PE-A"], [3000.0, 700.0], [FITC_SCALE, PE_SCALE], ""),
    ("high", ["Comp-PE-A", "Comp-BV786-A"], [2500.0, 1800.0], [PE_SCALE, BV_SCALE], "Cells"),
]


def curly_events(n: int = 40000, seed: int = 11) -> np.ndarray:
    """Events as acquired: a negative population, PE+, BV786+ and both,
    each detector's signal spread as photons are counted (its variance
    grows with the light) — so that, compensated, a negative population
    spreads wider the brighter it is on the other axis, which is what a
    curly quadrant is for."""
    rng = np.random.default_rng(seed)
    kind = rng.choice(4, n, p=[0.35, 0.25, 0.25, 0.15])
    pe = np.where((kind == 1) | (kind == 3), rng.lognormal(np.log(20000), 0.6, n), 0.0)
    bv = np.where((kind == 2) | (kind == 3), rng.lognormal(np.log(8000), 0.6, n), 0.0)
    fitc = rng.lognormal(np.log(1500), 1.2, n)
    light = np.column_stack([fitc, pe, bv]) @ SPILL
    detected = light + rng.normal(0.0, 1.0, light.shape) * np.sqrt(40.0 * light + 60.0 ** 2)
    scatter = rng.normal([90000, 50000], [25000, 15000], (n, 2))
    return np.column_stack([scatter, detected])


def curly_tree(name: str = "curly.fcs") -> GateNode:
    tree = GateNode(name)
    cells = tree.add(GateNode("Cells", PolygonGate(
        ["FSC-A", "SSC-A"], np.array([[30000.0, 15000], [150000, 20000], [160000, 90000],
                                      [40000, 85000]]))))
    for label, dims, centre, scales, parent in CROSSINGS:
        under = cells if parent else tree
        for q in QUADRANTS:
            under.add(GateNode(f"{label} {q}", CurlyQuadrantGate(
                list(dims), list(centre), q, label, True, list(scales))))
    return tree


def curly_files(folder: str) -> tuple[str, str]:
    """(the .fcs, the workspace MIL-C writes for it with `curly_tree`):
    the workspace `curly_compare.py` gives CytoML. NOT written by FlowJo."""
    fcs = os.path.join(folder, "curly.fcs")
    spill = ",".join([str(len(FLUO)), *FLUO, *(repr(float(v)) for v in SPILL.ravel())])
    write_fcs(fcs, curly_events(), CHANNELS, ranges=[262144.0] * len(CHANNELS),
              keywords={"SPILL": spill, "$FIL": "curly.fcs"})
    comp = Compensation(list(FLUO), SPILL.copy(), "Acquisition-defined")
    tree = curly_tree()
    sample = gate_export.ExportSample(
        name="curly.fcs", path=fcs, events=40000, channels=list(CHANNELS),
        keywords=read_header(fcs), compensation=comp, matrix_source="file",
        scales={"PE-A": PE_SCALE, "BV786-A": BV_SCALE, "FITC-A": FITC_SCALE}, gates=tree,
        load=lambda: Sample.open(fcs, compensation=comp, gates=copy_tree(tree)))
    path = os.path.join(folder, "curly.wsp")
    wsp_write.write_wsp([sample], path)
    return fcs, path


def opened(fcs: str, tree: GateNode) -> Sample:
    comp = Compensation(list(FLUO), SPILL.copy())
    return Sample.open(fcs, compensation=comp, gates=tree)


# --- what FlowJo writes ----------------------------------------------------------

#: a curly corner as FlowJo 10.10.0 wrote it, in the workspace allotropy
#: (Benchling's instrument-data converters) publishes with its tests:
#: tests/parsers/flowjo/testdata/flowjo_example_2.wsp, the population
#: "Q1: CD8- , CD4+" of "CD8, CD4 subset (Curly Gate)". The attributes,
#: the bounds and the scales are FlowJo's; the rest of the file is cut.
CURLY = """<Population name="{q}: CD8{sx} , CD4{sy}" count="{n}"><Gate gating:id="ID{k}">
 <gating:CurlyQuad eventsInside="1"  annoOffsetX="0"  annoOffsetY="0"  tint="#000000"  isTinted="0"  lineWeight="Normal"  userDefined="1"  percentX="0"  percentY="0" >
  <gating:dimension gating:{bx}="292.6491374776836" >
   <data-type:fcs-dimension data-type:name="Comp-PE-A" />
  </gating:dimension>
  <gating:dimension gating:{by}="434.6014486060488" >
   <data-type:fcs-dimension data-type:name="Comp-BV786-A" />
  </gating:dimension>
 </gating:CurlyQuad></Gate></Population>"""
FLOWJO = """<?xml version="1.0" encoding="UTF-8"?>
<Workspace version="20.0" flowJoVersion="10.10.0"
 xmlns:gating="http://www.isac-net.org/std/Gating-ML/v2.0/gating"
 xmlns:transforms="http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
 xmlns:data-type="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes">
 <Groups><GroupNode name="All Samples"><Subpopulations>{corners}</Subpopulations></GroupNode></Groups>
 <SampleList><Sample><DataSet uri="file:/Y:/x/Saline_004.fcs" sampleID="1"/>
  <Transformations>
   <transforms:biex transforms:length="256"  transforms:maxRange="262144"  transforms:neg="0"  transforms:width="-15.8489319246"  transforms:pos="4.418539922" >
    <data-type:parameter data-type:name="Comp-PE-A" /></transforms:biex>
   <transforms:biex transforms:length="256"  transforms:maxRange="262144"  transforms:neg="0"  transforms:width="-63.0957344480193"  transforms:pos="4.418539922" >
    <data-type:parameter data-type:name="Comp-BV786-A" /></transforms:biex>
  </Transformations>
  <SampleNode name="Saline_004.fcs" count="40"><Subpopulations>{corners}</Subpopulations></SampleNode>
 </Sample></SampleList>
</Workspace>"""
#: FlowJo's corners: which side each axis is bounded on, and the count it
#: wrote in that workspace
FLOWJO_CORNERS = {"Q1": ("max", "min", 5), "Q2": ("min", "min", 2), "Q3": ("min", "max", 6),
                  "Q4": ("max", "max", 0)}


def flowjo_workspace(tmp_path) -> str:
    corners = "".join(
        CURLY.format(q=q, sx="+" if bx == "min" else "-", sy="+" if by == "min" else "-",
                     n=n, k=k, bx=bx, by=by)
        for k, (q, (bx, by, n)) in enumerate(FLOWJO_CORNERS.items()))
    path = tmp_path / "flowjo.wsp"
    path.write_text(FLOWJO.format(corners=corners), encoding="utf-8")
    return str(path)


def test_flowjo_writes_a_curly_quadrant_as_a_quadrant(tmp_path):
    """A corner and a crossing point, nothing about the curve: read as a
    curly corner on the scales FlowJo gave the two axes, said to be read
    by a rule not checked against FlowJo."""
    w = wsp.read_wsp(flowjo_workspace(tmp_path))
    tree = w.samples[0].gates
    assert [c.name[:2] for c in tree.children] == list(FLOWJO_CORNERS)
    for node in tree.children:
        g = node.gate
        assert isinstance(g, CurlyQuadrantGate) and g.quadrant == node.name[:2]
        assert g.dims == ["Comp-PE-A", "Comp-BV786-A"]
        assert g.center == [292.6491374776836, 434.6014486060488]
        assert g.scales[0].width == pytest.approx(-15.8489319246)
        assert g.scales[1].width == pytest.approx(-63.0957344480193)
        assert (g.bend, g.power, g.channels) == (CURL_BEND, CURL_POWER, CURL_CHANNELS)
        assert node.reference_count == FLOWJO_CORNERS[g.quadrant][2]
        assert node.note.startswith(wsp.UNVERIFIED) and "curly" in node.note
    assert len({c.gate.group for c in tree.children}) == 1
    # a group's gates come without scales: left, and said why
    assert "All Samples" not in w.groups
    assert [(s.owner, s.kind) for s in w.skipped] == [("All Samples", "CurlyQuad")] * 4
    assert all("scale" in s.reason for s in w.skipped)


def test_a_linear_axis_and_one_with_no_scale():
    """cytolib measures a linear axis by its top alone, and one with no
    scale as linear to 262,144; so does this."""
    assert wsp.curl_scale(Linear(-200.0, 100000.0)) == Linear(0.0, 100000.0)
    assert wsp.curl_scale(None) == Linear(0.0, 262144.0)
    b = Biex()
    assert wsp.curl_scale(b) is b


# --- the geometry ------------------------------------------------------------------

def plane(n=30000, seed=3):
    rng = np.random.default_rng(seed)
    pe = np.concatenate([rng.normal(0, 300, n // 2), rng.lognormal(8, 1.5, n - n // 2)])
    bv = np.concatenate([rng.lognormal(7, 1.7, n // 2), rng.normal(100, 400, n - n // 2)])
    return {"Comp-PE-A": pe, "Comp-BV786-A": bv}


def four(centre=(292.6491374776836, 434.6014486060488), **curve):
    return {q: CurlyQuadrantGate(["Comp-PE-A", "Comp-BV786-A"], list(centre), q, "g", True,
                                 [PE_SCALE, BV_SCALE], **curve) for q in QUADRANTS}


def test_every_event_in_one_corner_and_straight_without_a_bend():
    data = plane()
    masks = {q: g.contains(data.__getitem__) for q, g in four().items()}
    assert np.array_equal(sum(m.astype(int) for m in masks.values()), np.ones(30000, int))
    # the arms bend away from the centre: Q2 loses to Q1 and Q3, Q4 keeps its own
    straight = {q: QuadrantGate(["Comp-PE-A", "Comp-BV786-A"], 292.6491374776836,
                                434.6014486060488, q).contains(data.__getitem__)
                for q in QUADRANTS}
    assert np.array_equal(masks["Q4"], straight["Q4"])
    assert masks["Q2"].sum() < straight["Q2"].sum()
    assert not (masks["Q2"] & ~straight["Q2"]).any()
    assert masks["Q1"].sum() > straight["Q1"].sum() and masks["Q3"].sum() > straight["Q3"].sum()
    for q, g in four(bend=0.0).items():
        assert np.array_equal(g.contains(data.__getitem__), straight[q]), q
    outside = four()["Q2"]
    outside.inside = False
    assert np.array_equal(outside.contains(data.__getitem__), ~masks["Q2"])


def cytolib_polygon(gate: CurlyQuadrantGate) -> np.ndarray:
    """cytolib's polygon for the corner, in channels — its construction,
    ported as it is (src/gate.cpp, CurlyQuadGate::interpolate): 40 points
    on each curve, the edges of the plot at 255 and -4000."""
    cx, cy = gate.centre_channels()
    top, low, n = 255.0, -4e3, 40
    d1 = (top - cx) / n * np.arange(n)
    curve1 = np.column_stack([cx + d1, 0.001 * d1 ** 2 + cy])
    d2 = (top - cy) / n * np.arange(n)
    curve2 = np.column_stack([0.001 * d2 ** 2 + cx, cy + d2])
    q = gate.quadrant
    if q == "Q1":
        return np.vstack([curve2, [[low, top], [low, cy], curve2[0]]])
    if q == "Q2":
        return np.vstack([curve1, [[top, top], curve2[-1]], curve2[::-1]])
    if q == "Q3":
        return np.vstack([curve1, [[top, low], [cx, low], [cx, cy]]])
    return np.array([[cx, cy], [cx, low], [low, low], [low, cy], [cx, cy]])


@pytest.mark.parametrize("centre", [(292.6491374776836, 434.6014486060488), (2500.0, 1800.0),
                                    (-150.0, 30.0)])
def test_the_curve_is_cytolibs(centre):
    """On the events inside cytolib's plot (below 255 channels on both
    axes, short of the strip its polygons leave between their last points
    and the corners of the plot), each corner holds what cytolib's
    polygon holds."""
    data = plane()
    for q, g in four(centre).items():
        x = g.to_channels(data["Comp-PE-A"], 0)
        y = g.to_channels(data["Comp-BV786-A"], 1)
        cx, cy = g.centre_channels()
        within = (x < cx + 0.97 * (255 - cx)) & (y < cy + 0.97 * (255 - cy))
        theirs = points_in_polygon(x, y, cytolib_polygon(g))
        ours = g.contains(data.__getitem__)
        assert within.sum() > 25000
        # 40 straight pieces against the curve: an event or two in between
        assert int((ours != theirs)[within].sum()) <= 2, q


def test_a_corner_as_a_polygon():
    """For a format that has no curly quadrant: the polygon of each corner
    holds the corner's events, inside the box it is cut from."""
    data = plane(seed=8)
    box = ((-0.1, 1.1), (-0.1, 1.1))
    total = np.zeros(30000, int)
    for q, g in four().items():
        poly = g.corner(*box)
        assert poly[0] == pytest.approx([PE_SCALE(292.6491374776836), BV_SCALE(434.6014486060488)])
        inside = points_in_polygon(PE_SCALE(data["Comp-PE-A"]), BV_SCALE(data["Comp-BV786-A"]),
                                   poly)
        assert int((inside != g.contains(data.__getitem__)).sum()) <= 2, q
        total += inside
    assert int((total != 1).sum()) <= 4
    assert four()["Q2"].corner((0.5, 1.0), (0.0, 1.0)) is None        # centre outside


def test_saved_copied_and_scaled():
    tree = curly_tree()
    data = {**plane(40000), "FSC-A": np.full(40000, 90000.0), "SSC-A": np.full(40000, 50000.0),
            "Comp-FITC-A": np.random.default_rng(1).lognormal(7, 1.5, 40000)}
    again = tree_from_dict(tree_to_dict(tree))
    assert tree_to_dict(again) == tree_to_dict(tree)
    copied = copy_tree(tree)
    for t in (again, copied):
        a, b = tree.apply(data.__getitem__, 40000), t.apply(data.__getitem__, 40000)
        assert all(np.array_equal(a[k], b[k]) for k in a)
    copied.children[1].gate.center[0] = 0.0
    assert tree.children[1].gate.center[0] != 0.0
    other = Arcsinh(150)
    assert follow_scale(copied, "Comp-PE-A", other)
    assert all(n.gate.scales[n.gate.dims.index("Comp-PE-A")] is other
               for n in copied.walk() if n.gate is not None and "Comp-PE-A" in n.gate.dims)
    # a curve of its own is kept
    odd = GateNode("", CurlyQuadrantGate(["a", "b"], [1.0, 2.0], "Q3", bend=0.012, power=1.5,
                                         channels=256.0))
    back = tree_from_dict(tree_to_dict(odd)).gate
    assert (back.bend, back.power, back.channels, back.scales) == (0.012, 1.5, 256.0, None)


def test_outline_is_four_arms_two_of_them_curved():
    from milc import views
    g = four()["Q1"]
    lines = views.outline(g, lambda dim: {"Comp-PE-A": PE_SCALE, "Comp-BV786-A": BV_SCALE}[dim],
                          [(0.0, 1.0), (0.0, 1.0)])
    assert len(lines) == 4
    centre = [PE_SCALE(292.6491374776836), BV_SCALE(434.6014486060488)]
    for line in lines:
        assert line[0] == pytest.approx(centre)
    right, up, left, down = lines
    assert np.all(np.diff(right[:, 1]) >= 0) and right[-1, 1] > centre[1] + 0.05
    assert np.all(np.diff(up[:, 0]) >= 0) and up[-1, 0] > centre[0] + 0.05
    assert np.allclose(left[:, 1], centre[1]) and np.allclose(down[:, 0], centre[0])


# --- against CytoML ------------------------------------------------------------------

#: CytoML 2.24.0 (flowWorkspace 4.24.0, cytolib 2.24.0) on the workspace
#: `curly_files` writes, read with flowjo_to_gatingset — measured on
#: 2026-11-23 with tests/reference/curly_compare.py
CYTOML = {"Cells": 38389, "Cells/high Q1": 9517, "Cells/high Q2": 5628, "Cells/high Q3": 9599,
          "Cells/high Q4": 13644, "CD8 CD4 Q1": 9547, "CD8 CD4 Q2": 6402, "CD8 CD4 Q3": 10607,
          "CD8 CD4 Q4": 13444, "FITC PE Q1": 15033, "FITC PE Q2": 956, "FITC PE Q3": 6670,
          "FITC PE Q4": 17339}


def test_the_workspace_and_cytoml(tmp_path):
    """Written as FlowJo writes curly corners, read back the same; counted
    as CytoML counts them, to the events its biexponential puts elsewhere.

    CytoML draws the curve on cytolib's own table of the biexponential,
    built on 255 channels and with a root found to the nearest whole step
    (transformation.cpp, biexpTrans::logRoot), where MIL-C uses FlowJo's
    scale as it draws every other gate: up to 65 events of 40,000 change
    corner (1% of a corner), all of them on the axes with that scale. With
    cytolib's table in its place the two put at most 1 event differently
    (curly_compare.py)."""
    fcs, path = curly_files(str(tmp_path))
    text = open(path, encoding="utf-8").read()
    assert text.count("<gating:CurlyQuad ") == 12 and "PolygonGate" in text
    w = wsp.read_wsp(path)
    assert w.skipped == []
    own = opened(fcs, curly_tree()).masks()
    back = opened(fcs, w.samples[0].gates).masks()
    for p, n in CYTOML.items():
        assert np.array_equal(own[p], back[p]), p
        assert w.samples[0].gates.find(p).reference_count == int(own[p].sum())
        assert abs(int(own[p].sum()) - n) <= max(70, 0.01 * n), p
    for p in ("Cells/high Q1", "CD8 CD4 Q4", "FITC PE Q3", "FITC PE Q4"):
        assert abs(int(own[p].sum()) - CYTOML[p]) <= 1, p


def test_written_as_polygons_where_a_workspace_would_draw_another_curve(tmp_path):
    """A curly corner whose curve a reader would not give back from the
    workspace — another bend, here — goes out as the polygon of its corner,
    and so does every curly corner in Gating-ML, which has none."""
    fcs, _ = curly_files(str(tmp_path))
    tree = curly_tree()
    for node in tree.walk():
        if isinstance(node.gate, CurlyQuadrantGate) and node.name.startswith("CD8"):
            node.gate.bend, node.gate.power, node.gate.channels = 0.012, 1.5, 256.0
    comp = Compensation(list(FLUO), SPILL.copy(), "Acquisition-defined")
    sample = gate_export.ExportSample(
        name="curly.fcs", path=fcs, events=40000, channels=list(CHANNELS),
        keywords=read_header(fcs), compensation=comp, matrix_source="file",
        scales={"PE-A": PE_SCALE, "BV786-A": BV_SCALE, "FITC-A": FITC_SCALE}, gates=tree,
        load=lambda: Sample.open(fcs, compensation=comp, gates=copy_tree(tree)))
    own = opened(fcs, copy_tree(tree)).masks()
    out = str(tmp_path / "bent.wsp")
    report = wsp_write.write_wsp([sample], out)
    assert report.skipped == [] and "written as a polygon" in report.text()
    text = open(out, encoding="utf-8").read()
    assert text.count("<gating:CurlyQuad ") == 8
    w = wsp.read_wsp(out)
    back = opened(fcs, w.samples[0].gates).masks()
    for p in CYTOML:
        kind = type(w.samples[0].gates.find(p).gate)
        assert kind is (CurlyQuadrantGate if " Q" in p and not p.startswith("CD8")
                        else PolygonGate), p
        assert int((own[p] != back[p]).sum()) <= 2, p

    xml = str(tmp_path / "curly.xml")
    report = gatingml.write_gatingml([sample], xml)
    assert report.skipped == [] and "curly quadrant is written as a polygon" in report.text()
    prepared = gate_export.prepare(sample, gate_export.Report())
    assert not any(isinstance(n.gate, CurlyQuadrantGate) for n in prepared.walk())
    again = opened(fcs, prepared).masks()
    for p in CYTOML:
        assert int((own[p] != again[p]).sum()) <= 2, p
    root = ET.parse(xml).getroot()
    assert len([e for e in root.iter() if e.tag.endswith("PolygonGate")]) == 13


# --- in the window -------------------------------------------------------------------

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    QtWidgets = pytest.importorskip("PyQt6.QtWidgets")
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


def test_a_curly_quadrant_in_the_window(app, tmp_path):
    """Read with the folder's workspace, drawn as its arms — two of them
    curves — with a handle at the centre that moves the four corners."""
    from milc.ui import worker
    from milc.ui.shell import MainShell

    curly_files(str(tmp_path))
    w = MainShell()
    w.open_folder(str(tmp_path))
    w.show()
    g = w.gates
    w.tabs.setCurrentWidget(g)
    for _ in range(3):
        worker.wait_all()
        app.processEvents()
    corners = [c for c in g.sample.gates.children if c.name.startswith("CD8 CD4")]
    assert all(isinstance(c.gate, CurlyQuadrantGate) for c in corners) and len(corners) == 4
    g.set_axes("Comp-PE-A", "Comp-BV786-A")
    [curl] = g._curls
    _, _, _, labels = g._crossings[-1]
    assert len(labels) == 4
    shown = curl.shape(curl.centre)
    right, up = shown["right"], shown["up"]
    assert right[0] == pytest.approx(curl.centre) and right[-1, 1] > curl.centre[1]
    assert up[-1, 0] > curl.centre[0]
    assert np.allclose(shown["left"][:, 1], curl.centre[1])

    # the centre moves all four, and the counts follow
    scale = g._scale("Comp-PE-A"), g._scale("Comp-BV786-A")
    target = (float(scale[0](np.float64(800.0))), float(scale[1](np.float64(600.0))))
    curl.move_centre(*target)
    for c in corners:
        assert c.gate.center == pytest.approx([800.0, 600.0], rel=1e-3)
    data = g.sample.column
    for c in corners:
        assert g.sample.count(c.path) == int((c.gate.contains(data)
                                               & g.sample.mask(c.parent.path)).sum())
    assert sum(g.sample.count(c.path) for c in corners) == len(g.sample)

    # seen from the other side: the same arms, mirrored
    g.set_axes("Comp-BV786-A", "Comp-PE-A")
    [mirrored] = g._curls
    assert mirrored.centre == pytest.approx(target[::-1])
    assert mirrored.shape(mirrored.centre)["right"][-1, 1] > mirrored.centre[1]
    w.close()
