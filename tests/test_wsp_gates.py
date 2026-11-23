"""
Quadrants, combinations and ellipses out of a FlowJo workspace.

None of this is FlowJo's word: the workspaces below are written by hand
the way FlowKit and CytoML read them and the way the workspaces FlowKit
publishes look (see milc/wsp.py for which file says what). What is tested
is that each construct becomes the gate those sources describe — every
gate against a calculation of its own in numpy — that what is read this
way says so, and that what is not read is listed with the reason.
"""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

from milc import wsp
from milc.gating import (BooleanGate, EllipseGate, QuadrantGate, RectangleGate, copy_tree,
                         tree_from_dict, tree_to_dict)
from milc.transforms import Biex, Linear

WSP = """<?xml version="1.0"?>
<Workspace flowJoVersion="10.8.1" xmlns:g="http://g" xmlns:t="http://t" xmlns:d="http://d">
 <Groups><GroupNode name="All Samples"><Subpopulations>{group}</Subpopulations></GroupNode></Groups>
 <SampleList><Sample><DataSet uri="file:/C:/x/a.fcs" sampleID="1"/>
 <Transformations>{scales}</Transformations>
 <SampleNode name="a.fcs" count="2000"><Subpopulations>{gates}</Subpopulations></SampleNode>
 </Sample></SampleList>
</Workspace>"""
LINEAR = """<t:linear t:minRange="{lo}" t:maxRange="262144" gain="1"><d:parameter d:name="{name}"/></t:linear>"""
BIEX = """<t:biex t:length="{length}" t:maxRange="262144" t:neg="0" t:width="-100" t:pos="4.42"><d:parameter d:name="{name}"/></t:biex>"""
POP = """<Population name="{name}" count="{n}"><Gate>{gate}</Gate><Subpopulations>{inner}</Subpopulations></Population>"""
RECT = """<g:RectangleGate eventsInside="1" percentX="0" percentY="0">
 <g:dimension {x}><d:fcs-dimension d:name="A"/></g:dimension>
 <g:dimension {y}><d:fcs-dimension d:name="B"/></g:dimension></g:RectangleGate>"""
NODE = """<{kind} name="{name}" count="{n}"><Subpopulations>{inner}</Subpopulations>
 <Dependents>{deps}</Dependents></{kind}>"""
ELLIPSE = """<g:EllipsoidGate eventsInside="{inside}" g:distance="89.2">
 <g:dimension><d:fcs-dimension d:name="A"/></g:dimension>
 <g:dimension><d:fcs-dimension d:name="B"/></g:dimension>
 <g:foci>{foci}</g:foci><g:edge>{edge}</g:edge></g:EllipsoidGate>"""

SCALES = LINEAR.format(lo=0, name="A") + BIEX.format(length=256, name="B")
#: the corners as FlowJo names and bounds them: Q1 is low on X, high on Y
CORNERS = {"Q1": ("max", "min"), "Q2": ("min", "min"), "Q3": ("min", "max"), "Q4": ("max", "max")}


def pop(name, gate, n=0, inner=""):
    return POP.format(name=name, gate=gate, n=n, inner=inner)


def rect(x="", y=""):
    return RECT.format(x=x, y=y)


def quadrants(x=30000, y=500, corners=CORNERS, counts=(0, 0, 0, 0)):
    return "".join(pop(f"{q}: A{'+' if sx == 'min' else '-'} , B{'+' if sy == 'min' else '-'}",
                       rect(f'g:{sx}="{x}"', f'g:{sy}="{y}"'), n)
                   for (q, (sx, sy)), n in zip(corners.items(), counts))


def node(kind, name, deps, n=0, inner=""):
    return NODE.format(kind=kind, name=name, n=n, inner=inner,
                       deps="".join(f'<Dependent name="{d}"/>' for d in deps))


def vertices(points):
    return "".join("<g:vertex>" + "".join(f'<g:coordinate d:value="{v!r}"/>' for v in p)
                   + "</g:vertex>" for p in np.asarray(points, dtype=float).tolist())


def ellipse(centre, a, b, angle, inside=1, whole=False, short_first=False):
    """As FlowJo writes one, on the plot's 256 channels: the two foci, the
    ends of the long axis, the ends of the short one. In its own files
    the ends are in whole channels and the foci are worked out from them
    (`whole`)."""
    centre = np.asarray(centre, dtype=float)
    u, v = np.array([np.cos(angle), np.sin(angle)]), np.array([-np.sin(angle), np.cos(angle)])
    edge = np.array([centre + a * u, centre - a * u, centre + b * v, centre - b * v])
    if whole:
        edge = np.round(edge)
        centre, u = edge[:2].mean(axis=0), (edge[0] - edge[1]) / np.hypot(*(edge[0] - edge[1]))
        a, b = np.hypot(*(edge[0] - edge[1])) / 2, np.hypot(*(edge[2] - edge[3])) / 2
    c = np.sqrt(a * a - b * b)
    return ELLIPSE.format(inside=inside, foci=vertices([centre - c * u, centre + c * u]),
                          edge=vertices(edge[[2, 3, 0, 1]] if short_first else edge))


def mismatched(**other):
    """The foci of one ellipse with the edge points of another."""
    one = dict(centre=(120, 110), a=70, b=30, angle=0.6)
    return (ellipse(**one).split("<g:edge>")[0] + "<g:edge>"
            + ellipse(**{**one, **other}).split("<g:edge>")[1])


def read(tmp_path, gates, scales=SCALES, group=""):
    path = tmp_path / "a.wsp"
    path.write_text(WSP.format(gates=gates, scales=scales, group=group))
    return wsp.read_wsp(str(path))


@pytest.fixture
def events():
    rng = np.random.default_rng(3)
    data = {"A": rng.uniform(0, 262144, 2000),
            "B": np.concatenate([rng.normal(0, 150, 1000), rng.lognormal(7, 1.5, 1000)])}
    return data.__getitem__


def masks(tree, events):
    return tree.apply(events, 2000)


# --- quadrants ---------------------------------------------------------------

def test_four_open_rectangles_are_one_quadrant_gate(tmp_path, events):
    w = read(tmp_path, pop("Cells", rect('g:min="1000"'), inner=quadrants()))
    cells = w.samples[0].gates.find("Cells")
    assert isinstance(cells.gate, RectangleGate) and cells.note == ""
    assert [type(c.gate) for c in cells.children] == [QuadrantGate] * 4
    assert [c.gate.quadrant for c in cells.children] == [c.name[:2] for c in cells.children]
    assert {(c.gate.x, c.gate.y, tuple(c.gate.dims)) for c in cells.children} == {(30000, 500, ("A", "B"))}
    assert len({c.gate.group for c in cells.children}) == 1 and cells.children[0].gate.group
    assert all(c.note.startswith(wsp.UNVERIFIED) for c in cells.children)
    assert w.skipped == []

    got = masks(w.samples[0].gates, events)
    a, b, parent = events("A"), events("B"), events("A") >= 1000
    assert np.array_equal(got["Cells/Q1: A- , B+"], parent & (a < 30000) & (b >= 500))
    assert np.array_equal(got["Cells/Q2: A+ , B+"], parent & (a >= 30000) & (b >= 500))
    assert np.array_equal(got["Cells/Q3: A+ , B-"], parent & (a >= 30000) & (b < 500))
    assert np.array_equal(got["Cells/Q4: A- , B-"], parent & (a < 30000) & (b < 500))
    # every event of the parent is in one corner and only one
    assert np.array_equal(sum(got[c.path].astype(int) for c in cells.children), parent.astype(int))


def test_two_crossings_under_one_parent_are_two_groups(tmp_path):
    w = read(tmp_path, quadrants(30000, 500) + quadrants(60000, 500))
    groups = [c.gate.group for c in w.samples[0].gates.children]
    assert len(set(groups)) == 2 and groups[:4] == [groups[0]] * 4


@pytest.mark.parametrize("gates", [
    # one corner missing
    quadrants(corners={q: CORNERS[q] for q in ("Q1", "Q2", "Q3")}),
    # one corner cut somewhere else
    quadrants(corners={q: CORNERS[q] for q in ("Q1", "Q2", "Q3")})
    + pop("Q4", rect('g:max="30000"', 'g:max="501"')),
    # a corner taken twice
    quadrants() + pop("again", rect('g:max="30000"', 'g:min="500"')),
    # a corner that keeps what is outside it
    quadrants(corners={q: CORNERS[q] for q in ("Q1", "Q2", "Q3")})
    + pop("Q4", rect('g:max="30000"', 'g:max="500"').replace('eventsInside="1"', 'eventsInside="0"')),
])
def test_rectangles_that_do_not_make_a_quadrant_stay_rectangles(tmp_path, events, gates):
    w = read(tmp_path, gates)
    tree = w.samples[0].gates
    assert all(isinstance(c.gate, RectangleGate) and c.note == "" for c in tree.children)
    assert wsp.unverified(tree) == [] and w.skipped == []
    a, b = events("A"), events("B")
    assert np.array_equal(masks(tree, events)["Q1: A- , B+"], (a < 30000) & (b >= 500))


def test_quadrants_of_group_and_samples_are_the_same_tree(tmp_path):
    """`Workspace.shared` compares the group's gates with the samples':
    the four corners have to be joined under the same name in both."""
    gates = pop("Cells", rect('g:min="1000"'), inner=quadrants())
    w = read(tmp_path, gates, group=gates)
    assert wsp.signature(w.groups["All Samples"]) == wsp.signature(w.samples[0].gates)


# --- combinations ------------------------------------------------------------

def test_combinations_of_sisters(tmp_path, events):
    inner = (node("OrNode", "either", ["Cells/high", "Cells/bright"])  # before what it names
             + pop("high", rect('g:min="100000"')) + pop("bright", rect(y='g:min="1000"'))
             + node("AndNode", "both", ["Cells/high", "Cells/bright"])
             + node("NotNode", "not high", ["Cells/high"],
                    inner=pop("dim", rect(y='g:max="1000"')))
             + node("NotNode", "neither", ["Cells/either"]))
    w = read(tmp_path, pop("Cells", rect('g:min="1000"'), inner=inner))
    cells = w.samples[0].gates.find("Cells")
    assert [c.name for c in cells.children] == ["either", "high", "bright", "both", "not high", "neither"]
    found = {c.name: (c.gate.op, c.gate.operands) for c in cells.children
             if isinstance(c.gate, BooleanGate)}
    assert found == {"either": ("or", ["high", "bright"]), "both": ("and", ["high", "bright"]),
                     "not high": ("not", ["high"]), "neither": ("not", ["either"])}
    assert [n.name for n in wsp.unverified(w.samples[0].gates)] == list(found)
    assert w.skipped == []

    got = masks(w.samples[0].gates, events)
    parent, high, bright = events("A") >= 1000, events("A") >= 100000, events("B") >= 1000
    assert np.array_equal(got["Cells/either"], parent & (high | bright))
    assert np.array_equal(got["Cells/both"], parent & high & bright)
    assert np.array_equal(got["Cells/not high"], parent & ~high)
    assert np.array_equal(got["Cells/not high/dim"], parent & ~high & ~bright)
    assert np.array_equal(got["Cells/neither"], parent & ~(high | bright))


def test_combinations_at_the_top_of_the_tree(tmp_path, events):
    w = read(tmp_path, pop("high", rect('g:min="100000"')) + node("NotNode", "low", ["high"]))
    assert np.array_equal(masks(w.samples[0].gates, events)["low"], events("A") < 100000)


@pytest.mark.parametrize("extra, reason", [
    # made of a population that is not beside it
    (node("AndNode", "x", ["Cells/high", "Elsewhere/bright"]), "not one of its sisters"),
    (node("OrNode", "x", ["high", "Cells/bright"]), "not one of its sisters"),
    (node("NotNode", "x", ["Cells/gone"]), "not one of its sisters"),
    # made of one that could not be read
    (node("OrNode", "x", ["Cells/high", "Cells/curly"]), "was not read"),
    # made of a combination that comes after it
    (node("NotNode", "x", ["Cells/later"]) + node("NotNode", "later", ["Cells/high"]),
     "listed after it"),
    # two sisters by the name it gives
    (node("NotNode", "x", ["Cells/high"]) + pop("high", rect('g:min="5"')), "two of its sisters"),
    # not the number of populations its kind is made of
    (node("NotNode", "x", ["Cells/high", "Cells/bright"]), "names 2 populations"),
    (node("AndNode", "x", ["Cells/high"]), "names 1 populations"),
    (node("OrNode", "x", []), "names 0 populations"),
])
def test_combinations_that_are_left_alone(tmp_path, extra, reason):
    inner = (pop("high", rect('g:min="100000"')) + pop("bright", rect(y='g:min="1000"'))
             + pop("curly", "<g:CurlyQuad/>") + extra)
    w = read(tmp_path, pop("Cells", rect('g:min="1000"'), inner=inner))
    cells = w.samples[0].gates.find("Cells")
    assert "x" not in [c.name for c in cells.children]
    [skipped] = [s for s in w.skipped if s.path == "Cells/x"]
    assert reason in skipped.reason and skipped.owner == "a.fcs"
    assert skipped.kind in ("AndNode", "OrNode", "NotNode")


# --- ellipses ----------------------------------------------------------------

def on_the_plot(events):
    """Where each event is drawn on FlowJo's 256 channels: A linear to
    262144, B biexponential."""
    return events("A") / 262144 * 256, Biex(256, 262144, 0, -100, 4.42)(events("B")) * 256


def by_the_foci(x, y, centre, a, b, angle):
    """An ellipse as the points whose distances to the two foci add up to
    no more than its long axis — not the way `EllipseGate` works it out."""
    c = np.sqrt(a * a - b * b) * np.array([np.cos(angle), np.sin(angle)])
    (x1, y1), (x2, y2) = np.asarray(centre) - c, np.asarray(centre) + c
    return np.hypot(x - x1, y - y1) + np.hypot(x - x2, y - y2) <= 2 * a


@pytest.mark.parametrize("centre, a, b, angle, short_first", [
    ((120.0, 110.0), 70.0, 30.0, 0.6, False),
    ((90.0, 130.0), 60.0, 45.0, -1.1, False),
    ((100.0, 90.0), 80.0, 25.0, np.pi / 2, False),
    ((128.0, 100.0), 40.0, 40.0, 0.0, False),  # a circle: the foci are one point
    ((120.0, 110.0), 70.0, 30.0, 0.6, True),
])
def test_ellipse_is_the_one_drawn_on_the_256_channels(tmp_path, events, centre, a, b, angle,
                                                      short_first):
    drawn = ellipse(centre, a, b, angle, short_first=short_first)
    w = read(tmp_path, pop("Cells", rect('g:min="1000"'),
                           inner=pop("round", drawn, inner=pop("in", rect()))))
    tree = w.samples[0].gates
    found = tree.find("Cells/round")
    assert isinstance(found.gate, EllipseGate) and found.note.startswith(wsp.UNVERIFIED)
    assert wsp.unverified(tree) == [found] and w.skipped == []
    assert found.gate.center == pytest.approx([centre[0] / 256, centre[1] / 256])
    assert found.gate.radii == pytest.approx([a / 256, b / 256])

    x, y = on_the_plot(events)
    expected = (events("A") >= 1000) & by_the_foci(x, y, centre, a, b, angle)
    got = masks(tree, events)["Cells/round"]
    assert 100 < expected.sum() < 1900
    # the two sums differ in the last place for a point on the edge itself
    assert (got != expected).sum() <= 1


def test_ellipse_that_keeps_the_outside(tmp_path, events):
    w = read(tmp_path, pop("out", ellipse((120, 110), 70, 30, 0.6, inside=0)))
    x, y = on_the_plot(events)
    got = masks(w.samples[0].gates, events)["out"]
    assert (got != ~by_the_foci(x, y, (120, 110), 70, 30, 0.6)).sum() <= 1


def test_ellipse_with_its_edge_in_whole_channels(tmp_path):
    # the ends fall on (178, 150), (62, 71) and (103, 135), (138, 85)
    w = read(tmp_path, pop("round", ellipse((120.4, 110.3), 70.2, 30.4, 0.6, whole=True)))
    gate = w.samples[0].gates.find("round").gate
    assert np.array(gate.center) * 256 == pytest.approx([120, 110.5])
    assert np.array(gate.radii) * 256 == pytest.approx([np.hypot(116, 79) / 2, np.hypot(35, 50) / 2])
    assert gate.angle == pytest.approx(np.arctan2(79, 116))


def test_ellipse_of_flowkits_published_workspace():
    """The numbers of `single_ellipse_51_events.wsp` (FlowJo 10.6.1, in
    FlowKit's data): the third and fourth edge points are not at right
    angles to the first two, and the ellipse is still the one both
    readers make of it. FlowJo counted 51 of that file's 100 events in
    it, and so did this, measured by hand on 2026-11-23 — the events are
    not here to count again."""
    xml = ELLIPSE.format(inside=1,
                         foci=vertices([[62.7724519002, 157.4044547167], [94.2275480998, 93.5955452833]]),
                         edge=vertices([[96, 90], [61, 161], [67, 113], [90, 139]]))
    shape = ET.fromstring(f'<Gate xmlns:g="http://g" xmlns:d="http://d">{xml}</Gate>')[0]
    gate = wsp._read_ellipse(shape, ["A", "B"], True, {"A": Linear(), "B": Linear()})
    assert np.array(gate.center) * 256 == pytest.approx([78.5, 125.5])
    assert np.array(gate.radii) * 256 == pytest.approx([39.579, 17.357], abs=1e-3)
    assert gate.angle == pytest.approx(np.arctan2(-63.8089094334, 31.4550961996))


@pytest.mark.parametrize("gate, scales, kind, reason", [
    # an axis whose scale is not read (FlowJo's logicle), or is not given
    (ellipse((120, 110), 70, 30, 0.6), LINEAR.format(lo=0, name="A"), "EllipsoidGate", "that of B"),
    # a linear axis that does not start at zero: the two readers differ
    (ellipse((120, 110), 70, 30, 0.6),
     LINEAR.format(lo=-100, name="A") + BIEX.format(length=256, name="B"),
     "EllipsoidGate", "that of A"),
    # a biexponential that is not on 256 channels
    (ellipse((120, 110), 70, 30, 0.6),
     LINEAR.format(lo=0, name="A") + BIEX.format(length=4096, name="B"),
     "EllipsoidGate", "that of B"),
    # edge points that are not on the ellipse the foci give: wider, moved, turned
    (mismatched(b=32), SCALES, "EllipsoidGate", "same ellipse"),
    (mismatched(centre=(122, 110)), SCALES, "EllipsoidGate", "same ellipse"),
    (mismatched(angle=0.62), SCALES, "EllipsoidGate", "same ellipse"),
    (mismatched(a=20, b=10), SCALES, "EllipsoidGate", "foci do not lie inside"),
    # not two foci and four points
    (ELLIPSE.format(inside=1, foci=vertices([[100, 100]]), edge=vertices([[1, 2]] * 4)),
     SCALES, "EllipsoidGate", "two foci and four edge points"),
    # a curly quadrant, a gate of a kind never met, a population without a gate
    ("<g:CurlyQuad/>", SCALES, "CurlyQuad", "curly quadrant"),
    ("<g:SplineGate/>", SCALES, "SplineGate", "SplineGate"),
    ("", SCALES, "Population", "holds no gate"),
])
def test_what_is_not_read_is_listed_with_the_reason(tmp_path, gate, scales, kind, reason):
    inner = pop("odd", gate, inner=pop("under", rect(), inner=node("NotNode", "n", ["x"])))
    w = read(tmp_path, pop("Cells", rect('g:min="1000"'), inner=inner + pop("kept", rect())), scales)
    assert [n.path for n in w.samples[0].gates.walk()][1:] == ["Cells", "Cells/kept"]
    [skipped] = w.skipped
    assert (skipped.owner, skipped.path, skipped.below) == ("a.fcs", "Cells/odd", 2)
    assert skipped.kind == kind and reason in skipped.reason


def test_ellipse_on_a_group_is_left_alone(tmp_path):
    """A group's gates come without scales, and an ellipse cannot be
    placed without them; the sample's own copy is read."""
    gates = pop("round", ellipse((120, 110), 70, 30, 0.6))
    w = read(tmp_path, gates, group=pop("Cells", rect('g:min="1000"')) + gates)
    assert isinstance(w.samples[0].gates.find("round").gate, EllipseGate)
    assert [n.name for n in w.groups["All Samples"].children] == ["Cells"]
    [skipped] = w.skipped
    assert skipped.owner == "All Samples" and "group" in skipped.reason


# --- what the window is given ------------------------------------------------

def test_the_note_stays_with_the_gate(tmp_path):
    w = read(tmp_path, quadrants() + pop("plain", rect('g:min="5"'))
             + pop("round", ellipse((120, 110), 70, 30, 0.6))
             + node("NotNode", "not plain", ["plain"]))
    tree = w.samples[0].gates
    marked = [n.name for n in wsp.unverified(tree)]
    assert len(marked) == 6 and "plain" not in marked and tree.find("plain").note == ""
    for other in (copy_tree(tree), copy_tree(tree, keep_counts=False),
                  tree_from_dict(tree_to_dict(tree))):
        assert [n.name for n in wsp.unverified(other)] == marked
        assert other.find("round").note == tree.find("round").note
    assert "note" not in tree_to_dict(tree.find("plain"))


def test_counts_that_are_not_flowjos_are_listed(tmp_path, events):
    a, b = events("A"), events("B")
    high = int((a >= 100000).sum())
    corners = [int((x & y).sum()) for x, y in ((a < 30000, b >= 500), (a >= 30000, b >= 500),
                                               (a >= 30000, b < 500), (a < 30000, b < 500))]
    gates = (pop("high", rect('g:min="100000"'), high)
             # within what was measured on polygons: 1.2%, or three events
             + pop("near", rect('g:min="100000"'), high + int(0.012 * high))
             + pop("few", rect('g:min="262140"'), int((a >= 262140).sum()) + 3)
             + pop("far", rect('g:min="100000"'), high + int(0.013 * high) + 2)
             + quadrants(counts=[corners[0] // 2] + corners[1:]))
    w = read(tmp_path, gates)
    tree = w.samples[0].gates
    got = masks(tree, events)
    found = wsp.count_mismatches(tree, lambda path: got[path].sum())
    by_path = {m.path: m for m in found}
    assert set(by_path) == {"far", "Q1: A- , B+"}
    assert (by_path["far"].count, by_path["far"].flowjo) == (high, high + int(0.013 * high) + 2)
    assert by_path["far"].note == "" and by_path["Q1: A- , B+"].note.startswith(wsp.UNVERIFIED)
    assert (by_path["Q1: A- , B+"].count, by_path["Q1: A- , B+"].flowjo) == (corners[0], corners[0] // 2)
    assert wsp.count_mismatches(tree, lambda path: got[path].sum(), tolerance=2.0) == []
    # a tree without FlowJo's counts has nothing to be compared with
    assert wsp.count_mismatches(copy_tree(tree, keep_counts=False), lambda path: 0) == []
