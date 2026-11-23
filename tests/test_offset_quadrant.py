"""A quadrant whose arms are moved one by one: the geometry, what is kept
of it, and its handles in the window."""

import json
import os

import numpy as np
import pytest

from milc.gating import (ARMS, QUADRANTS, GateNode, OffsetQuadrantGate, QuadrantGate, copy_tree,
                         follow_scale, offset_quadrant, tree_from_dict, tree_to_dict)
from milc.transforms import Biex


def events(n=20000, seed=0):
    rng = np.random.default_rng(seed)
    x, y = rng.normal(1.0, 3.0, n), rng.normal(2.0, 3.0, n)
    # on the centre, on each arm, and on the lines past the centre
    x[:9] = [1, 4, 1, -3, 1, 2, 0, 5, -5]
    y[:9] = [2, 2, 6, 2, -4, 3, 1, 2, 2]
    return {"a": x, "b": y}


def corners(center, arms, scales=None):
    return {q: OffsetQuadrantGate(["a", "b"], list(center), arms, q, "g", scales=scales)
            for q in QUADRANTS}


def by_angle(data, center, arms):
    """Each event's corner from the angle it is seen at from the centre:
    another way to the same answer, away from the arms themselves."""
    theta = np.arctan2(data["b"] - center[1], data["a"] - center[0])
    at = {a: np.arctan2(arms[a][1] - center[1], arms[a][0] - center[0]) for a in ARMS}

    def turn(angle):
        return (angle - at["right"]) % (2 * np.pi)
    t, up, left, down = turn(theta), turn(at["up"]), turn(at["left"]), turn(at["down"])
    return {"Q2": t < up, "Q1": (t > up) & (t < left), "Q4": (t > left) & (t < down),
            "Q3": t > down}


def test_arms_along_the_axes_are_the_plain_quadrant():
    data = events()
    for scales in (None, [Biex(), None]):
        for q in QUADRANTS:
            plain = QuadrantGate(["a", "b"], 1.0, 2.0, q, "g")
            moved = offset_quadrant(plain, scales, reach=(0.2, 5.0))
            assert moved.ordered() and moved.center == [1.0, 2.0] and moved.group == "g"
            assert moved.arms["up"][0] == 1.0 and moved.arms["left"][1] == 2.0
            assert np.array_equal(plain.contains(data.__getitem__),
                                  moved.contains(data.__getitem__)), (q, scales)
            assert np.array_equal(~plain.contains(data.__getitem__),
                                  OffsetQuadrantGate(moved.dims, moved.center, moved.arms, q,
                                                     inside=False).contains(data.__getitem__))


def test_tilted_arms_cut_where_the_arms_are():
    data = events()
    center = (1.0, 2.0)
    arms = {"right": [6.0, 1.0], "up": [2.5, 7.0], "left": [-4.0, 4.0], "down": [0.0, -5.0]}
    gates = corners(center, arms)
    assert gates["Q1"].ordered()
    masks = {q: g.contains(data.__getitem__) for q, g in gates.items()}
    # every event in exactly one corner
    assert np.array_equal(sum(m.astype(int) for m in masks.values()), np.ones(20000, int))
    expected = by_angle(data, center, arms)
    off_the_arms = np.ones(20000, bool)
    off_the_arms[0] = False       # the centre itself has no angle
    for q in QUADRANTS:
        assert np.array_equal(masks[q][off_the_arms], expected[q][off_the_arms]), q
    # known points: just either side of the `up` arm, which leans right
    probe = {"a": np.array([2.4, 2.6, 1.0, 1.0]), "b": np.array([7.0, 7.0, 2.0, 9.0])}
    assert gates["Q1"].contains(probe.__getitem__).tolist() == [True, False, False, True]
    assert gates["Q2"].contains(probe.__getitem__).tolist() == [False, True, True, False]


def test_an_event_on_an_arm_is_in_one_corner_and_the_plain_quadrant_s():
    center = (0.0, 0.0)
    arms = {"right": [2.0, 1.0], "up": [-1.0, 2.0], "left": [-2.0, -1.0], "down": [1.0, -2.0]}
    on = {"a": np.array([0.0, 4.0, -2.0, -4.0, 2.0]), "b": np.array([0.0, 2.0, 4.0, -2.0, -4.0])}
    found = {q: g.contains(on.__getitem__).tolist() for q, g in corners(center, arms).items()}
    # the centre, then a point on right, up, left, down
    assert found == {"Q2": [True, True, True, False, False],
                     "Q1": [False, False, False, True, False],
                     "Q3": [False, False, False, False, True],
                     "Q4": [False, False, False, False, False]}


def test_a_corner_wider_than_half_a_turn():
    data = events(seed=3)
    center = (1.0, 2.0)
    # right and up both lean down-right and up-left: Q2 is narrow, Q4 holds 200 degrees
    arms = {"right": [5.0, 6.0], "up": [0.0, 7.0], "left": [-1.0, 5.0], "down": [6.0, 4.0]}
    gates = corners(center, arms)
    assert gates["Q4"].ordered()
    masks = {q: g.contains(data.__getitem__) for q, g in gates.items()}
    assert np.array_equal(sum(m.astype(int) for m in masks.values()), np.ones(20000, int))
    expected = by_angle(data, center, arms)
    for q in QUADRANTS:
        assert np.array_equal(masks[q][9:], expected[q][9:]), q
    assert masks["Q4"].mean() > 0.5


def test_degenerate_arms():
    data = events()
    center = (1.0, 2.0)
    straight = {"right": [5.0, 2.0], "up": [1.0, 6.0], "left": [-3.0, 2.0], "down": [1.0, -2.0]}
    plain = {q: QuadrantGate(["a", "b"], 1.0, 2.0, q).contains(data.__getitem__)
             for q in QUADRANTS}
    # an arm whose point is the centre has no direction: it is taken along its axis
    for arm in ARMS:
        arms = {**straight, arm: list(center)}
        gates = corners(center, arms)
        assert gates["Q1"].ordered()
        for q in QUADRANTS:
            assert np.array_equal(gates[q].contains(data.__getitem__), plain[q]), (arm, q)
    # one arm on top of its neighbour: the corner between them is empty but for the
    # arm, the arms are no longer in order, and the events are still shared out
    arms = {**straight, "up": [9.0, 2.0]}
    gates = corners(center, arms)
    assert not gates["Q2"].ordered()
    masks = {q: g.contains(data.__getitem__) for q, g in gates.items()}
    assert np.array_equal(sum(m.astype(int) for m in masks.values()), np.ones(20000, int))
    assert np.array_equal(masks["Q2"], (data["b"] == 2.0) & (data["a"] >= 1.0))
    # an arm past its neighbour, arms in line with each other, nothing finite
    for arms in ({**straight, "up": [5.0, 1.0]}, {**straight, "left": [9.0, 2.0]},
                 {**straight, "down": [float("nan"), 1.0]}):
        masks = [g.contains(data.__getitem__) for g in corners(center, arms).values()]
        assert np.array_equal(sum(m.astype(int) for m in masks), np.ones(20000, int))
    assert not corners(center, {**straight, "up": [5.0, 1.0]})["Q1"].ordered()
    # opposite arms in one line are an ordinary crossing turned
    turned = {"right": [2.0, 3.0], "up": [0.0, 3.0], "left": [0.0, 1.0], "down": [2.0, 1.0]}
    assert corners(center, turned)["Q1"].ordered()


def test_tested_on_the_scale_it_was_drawn_on():
    biex = Biex()
    rng = np.random.default_rng(1)
    data = {"a": rng.lognormal(6, 2, 5000) - 200, "b": rng.lognormal(6, 2, 5000) - 200}
    center = [500.0, 800.0]
    arms = {"right": [50000.0, 2000.0], "up": [900.0, 60000.0], "left": [-150.0, 800.0],
            "down": [500.0, -150.0]}
    on_scale = corners(center, arms, [biex, biex])["Q2"].contains(data.__getitem__)
    shown = {k: biex(v) for k, v in data.items()}
    same = OffsetQuadrantGate(["a", "b"], biex(np.array(center)).tolist(),
                              {a: biex(np.array(p)).tolist() for a, p in arms.items()},
                              "Q2").contains(shown.__getitem__)
    assert np.array_equal(on_scale, same)
    # a straight arm on the biexponential scale is not straight in the file's units
    assert not np.array_equal(on_scale, corners(center, arms)["Q2"].contains(data.__getitem__))
    # given another scale, the arms are straight there: as a polygon follows its axis
    tree = GateNode("")
    node = tree.add(GateNode("Q2", corners(center, arms, [biex, biex])["Q2"]))
    assert follow_scale(tree, "a", None) and node.gate.scales == [None, biex]
    assert follow_scale(tree, "b", None) and node.gate.scales is None
    assert not follow_scale(tree, "c", None)


def test_kept_in_the_project_and_copied():
    center = [1.0, 2.0]
    arms = {"right": [6.0, 1.0], "up": [2.5, 7.0], "left": [-4.0, 4.0], "down": [0.0, -5.0]}
    tree = GateNode("")
    for q, gate in corners(center, arms, [Biex(), None]).items():
        tree.add(GateNode(q, gate))
    plain = tree.add(GateNode("plain", QuadrantGate(["a", "b"], 3.0, 4.0, "Q3", "old")))
    again = tree_from_dict(json.loads(json.dumps(tree_to_dict(tree))))
    data = events()
    for before, after in zip(tree.children[:4], again.children):
        g = after.gate
        assert isinstance(g, OffsetQuadrantGate) and g.quadrant == before.name
        assert (g.center, g.arms, g.group) == (center, arms, "g") and g.scales[1] is None
        assert np.array_equal(g.contains(data.__getitem__),
                              before.gate.contains(data.__getitem__))
    assert again.children[4].gate == plain.gate      # and the plain one stays plain

    copy = copy_tree(tree)
    copy.children[0].gate.arms["up"][0] = 99.0
    copy.children[0].gate.center[0] = 5.0
    assert tree.children[0].gate.arms["up"] == [2.5, 7.0]
    assert tree.children[0].gate.center == [1.0, 2.0]


def test_a_project_saved_before_there_were_arms():
    old = {"name": "", "children": [
        {"name": "CD4+ CD8-", "children": [],
         "gate": {"type": "quadrant", "dims": ["a", "b"], "inside": True, "x": 1.0, "y": 2.0,
                  "quadrant": "Q3", "group": "a|b|0"}}]}
    tree = tree_from_dict(old)
    assert tree.children[0].gate == QuadrantGate(["a", "b"], 1.0, 2.0, "Q3", "a|b|0")
    assert tree.children[0].fit is None and tree_to_dict(tree) == old


def test_outline_is_the_four_arms():
    from milc import views
    gate = OffsetQuadrantGate(["a", "b"], [1.0, 2.0],
                              {"right": [6.0, 1.0], "up": [2.5, 7.0], "left": [-4.0, 4.0],
                               "down": [0.0, -5.0]}, "Q1")
    lines = views.outline(gate, lambda dim: None, [(-10.0, 10.0), (-10.0, 10.0)])
    assert len(lines) == 4
    for line, arm in zip(lines, ARMS):
        assert line[0] == pytest.approx([1.0, 2.0])
        d, along = line[-1] - line[0], np.array(gate.arms[arm]) - [1.0, 2.0]
        assert d[0] * along[1] - d[1] * along[0] == pytest.approx(0, abs=1e-9)
        assert d @ along > 0 and max(abs(line[-1])) > 10      # out of the plot


def test_a_corner_as_a_polygon():
    """For a program that has polygons and no such quadrant: the polygon
    of a corner holds the events the corner does, inside the box."""
    from milc.gating import points_in_polygon
    data = events(seed=5)
    box = ((-12.0, 14.0), (-11.0, 15.0))
    within = ((data["a"] > box[0][0]) & (data["a"] < box[0][1])
              & (data["b"] > box[1][0]) & (data["b"] < box[1][1]))
    within[:9] = False                  # not the ones put on the arms
    assert within.sum() > 19000
    for arms in ({"right": [6.0, 1.0], "up": [2.5, 7.0], "left": [-4.0, 4.0],
                  "down": [0.0, -5.0]},
                 {"right": [5.0, 6.0], "up": [0.0, 7.0], "left": [-1.0, 5.0],
                  "down": [6.0, 4.0]}):          # one corner wider than half a turn
        total = np.zeros(20000, int)
        for q, gate in corners((1.0, 2.0), arms).items():
            polygon = gate.corner(*box)
            assert polygon[0] == pytest.approx([1.0, 2.0])
            inside = points_in_polygon(data["a"], data["b"], polygon)
            assert np.array_equal(inside[within], gate.contains(data.__getitem__)[within]), q
            total += inside
        assert np.array_equal(total[within], np.ones(int(within.sum()), int))
    # along the axes, Q2 is the rectangle from the centre to the box's corner
    plain = offset_quadrant(QuadrantGate(["a", "b"], 1.0, 2.0, "Q2"))
    assert plain.corner(*box) == pytest.approx(
        np.array([[1.0, 2.0], [14.0, 2.0], [14.0, 15.0], [1.0, 15.0]]))
    assert plain.corner((2.0, 5.0), (0.0, 5.0)) is None      # the centre is outside


# --- in the window ----------------------------------------------------------------

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    QtWidgets = pytest.importorskip("PyQt6.QtWidgets")
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


@pytest.fixture
def folder(tmp_path):
    from milc.fcs import write_fcs
    rng = np.random.default_rng(0)
    centres = np.array([[300, 3000, 100], [3000, 200, 100], [200, 300, 4000.0]])
    for name, n in (("3Ctrl1", 600), ("3Case1", 900)):
        pop = rng.integers(0, 3, n)
        fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (n, 3)))
        ev = np.column_stack([rng.normal(90000, 9000, n), fluo]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), ev, ["FSC-A", "FITC-A", "PE-A", "APC-A"],
                  markers=["", "CD8", "CD4", "Foxp3"],
                  keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
    return str(tmp_path)


def test_the_arms_of_a_quadrant_in_the_window(app, folder):
    from milc.session import Session
    from milc.ui import worker
    from milc.ui.shell import MainShell

    w = MainShell()
    w.open_folder(folder)
    w.show()
    g = w.gates
    w.tabs.setCurrentWidget(g)
    for _ in range(3):
        worker.wait_all()
        app.processEvents()
    said = []
    g.sigStatus.connect(said.append)
    s = g.sample
    t = Biex()
    pe, fitc = t(s.column("Comp-PE-A")), t(s.column("Comp-FITC-A"))
    shown = {"a": pe, "b": fitc}      # what `corners` calls the two axes
    g.set_axes("Comp-PE-A", "Comp-FITC-A")
    made = g.add_quadrants(0.4, 0.4)

    # as drawn: two lines, a handle on each arm inside the plot, none at the centre
    [arms] = g._arms
    assert not arms.offset and arms.vertical.isVisible() and not arms.middle.isVisible()
    assert all(isinstance(n.gate, QuadrantGate) for n in made)
    assert arms.ends["right"] == pytest.approx([0.4 + 0.3 * (arms.xr[1] - arms.xr[0]), 0.4])
    assert arms.ends["up"][0] == pytest.approx(0.4) and arms.ends["up"][1] > 0.4
    # dragging a line takes the handles with it, and the gates stay plain quadrants
    arms.vertical.setValue(0.42)
    assert arms.ends["down"][0] == pytest.approx(0.42)
    assert arms.handles["down"].pos().x() == pytest.approx(0.42)
    g._crossing_moved(g._crossings[0][0], 0.42, 0.4)
    assert all(isinstance(n.gate, QuadrantGate) for n in made)
    assert g.sample.count(made[1].path) == ((pe >= 0.42) & (fitc >= 0.4)).sum()

    # dragging one arm's handle tilts that arm alone: through the CD8+ cells
    before = {n.name: g.sample.count(n.path) for n in made}
    arms.move_arm("up", 0.24, 0.53)
    assert arms.offset and not arms.vertical.isVisible() and arms.middle.isVisible()
    assert all(isinstance(n.gate, OffsetQuadrantGate) for n in made)
    gate = made[0].gate
    assert gate.scales is not None and gate.dims == ["Comp-PE-A", "Comp-FITC-A"]
    assert t(np.array(gate.center)) == pytest.approx([0.42, 0.4])
    assert t(np.array(gate.arms["up"])) == pytest.approx([0.24, 0.53])
    expected = corners((0.42, 0.4), {a: arms.ends[a].tolist() for a in ARMS})
    after = {n.name: g.sample.count(n.path) for n in made}
    for n in made:
        assert after[n.name] == expected[n.gate.quadrant].contains(shown.__getitem__).sum()
    assert sum(after.values()) == len(s) and after != before
    assert 0 < after["CD4+ CD8+"] < before["CD4- CD8+"]
    assert after["CD4- CD8-"] == before["CD4- CD8-"]     # that corner has neither arm
    # the label on the plot and the tree follow
    label = next(lab for child, lab in g._crossings[0][3] if child is made[0])
    assert label.toPlainText() == f"{made[0].name}  {100 * after[made[0].name] / len(s):.1f}%"

    # an arm cannot go past its neighbour: it is put back, and nothing changes
    arms.move_arm("up", 0.9, 0.3)
    assert "past" in said[-1] and arms.ends["up"] == pytest.approx([0.24, 0.53])
    assert arms.handles["up"].pos().y() == pytest.approx(0.53)
    assert {n.name: g.sample.count(n.path) for n in made} == after

    # the centre moves them all
    arms.move_centre(0.45, 0.38)
    assert arms.ends["up"] == pytest.approx([0.27, 0.51])
    moved = corners((0.45, 0.38), {a: arms.ends[a].tolist() for a in ARMS})
    assert g.sample.count(made[1].path) == moved["Q2"].contains(shown.__getitem__).sum()

    # seen from the other side, the same arms mirrored — and still to be moved
    g.set_axes("Comp-FITC-A", "Comp-PE-A")
    [mirrored] = g._arms
    assert mirrored.offset and mirrored.centre == pytest.approx([0.38, 0.45])
    assert mirrored.ends["right"] == pytest.approx([0.51, 0.27])
    mirrored.move_arm("right", 0.55, 0.30)
    g.set_axes("Comp-PE-A", "Comp-FITC-A")
    assert g._arms[0].ends["up"] == pytest.approx([0.30, 0.55])
    final = corners((0.45, 0.38), {a: g._arms[0].ends[a].tolist() for a in ARMS})
    counts = {n.gate.quadrant: g.sample.count(n.path) for n in made}
    assert counts == {q: int(gate.contains(shown.__getitem__).sum())
                      for q, gate in final.items()}

    # every sample has it, and it is there when the folder is opened again
    assert all(isinstance(e.gates.children[0].gate, OffsetQuadrantGate)
               for e in w.session.usable())
    w.close()
    again = Session()
    again.open_folder(folder)
    first = again.usable()[0]
    assert isinstance(first.gates.children[0].gate, OffsetQuadrantGate)
    assert {n.gate.quadrant: again.load(first).count(n.path)
            for n in first.gates.children} == counts
