"""Gates the program proposes, on events made up so that the right answer
is known: which events are the population, where the valley is."""

import os

import numpy as np
import pytest

from milc import autogate as ag
from milc.gating import (EllipseGate, GateNode, PolygonGate, RectangleGate, copy_tree,
                         keep_fits, tree_from_dict, tree_to_dict)
from milc.transforms import Biex

UNIT = (0.0, 1.0)


# --- the valley of a histogram ----------------------------------------------------

def two_peaks(n=20000, share=0.3, seed=0, low=(0.25, 0.03), high=(0.62, 0.06)):
    rng = np.random.default_rng(seed)
    positive = rng.random(n) < share
    x = np.where(positive, rng.normal(*high, n), rng.normal(*low, n))
    return x, positive


def test_valley_between_two_peaks():
    x, positive = two_peaks()
    found = ag.valley(x, *UNIT)
    assert 0.36 < found.threshold < 0.46            # where the two densities cross: 0.40
    assert found.modes == pytest.approx((0.25, 0.62), abs=0.02)
    assert found.depth < 0.05 and found.above == pytest.approx(0.3, abs=0.01)
    assert ((x >= found.threshold) == positive).mean() > 0.995
    assert found.why.startswith("valley at 0% of the lower peak") and "above it" in found.why
    # the same events, the same answer; and on an axis ten times as long
    assert ag.valley(x, *UNIT).threshold == found.threshold
    assert ag.valley(x * 10, 0.0, 10.0).threshold == pytest.approx(10 * found.threshold)


def test_valley_finds_a_small_population_beside_a_tall_peak():
    # 3% positive and wide: its peak is a hundredth of the negative one's height
    x, positive = two_peaks(n=30000, share=0.03, low=(0.25, 0.01), high=(0.65, 0.08))
    found = ag.valley(x, *UNIT)
    assert found.threshold is not None and 0.3 < found.threshold < 0.5
    assert ((x >= found.threshold) == positive).mean() > 0.995


@pytest.mark.parametrize("x, reason", [
    # one population
    (np.random.default_rng(0).normal(0.4, 0.05, 20000), "one peak only"),
    # one population with a shoulder that never dips: what a dim marker looks like
    (np.concatenate([np.random.default_rng(1).normal(0.30, 0.02, 16000),
                     np.random.default_rng(2).normal(0.38, 0.06, 4000)]), "one peak only"),
    # a long tail of noise is not a second population
    (np.concatenate([np.random.default_rng(3).normal(0.3, 0.03, 20000),
                     np.random.default_rng(4).uniform(0.3, 0.9, 400)]), "one peak only"),
    # nothing to count
    (np.random.default_rng(5).normal(0.4, 0.05, 60), "too few events"),
    (np.full(500, 0.4), "same value"),
])
def test_no_valley_is_said_not_guessed(x, reason):
    found = ag.valley(x, *UNIT)
    assert found.threshold is None and reason in found.why


def test_valley_that_is_shallow_or_one_sided_is_not_proposed():
    # two peaks that overlap: a dip, but most of the way up
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(0.40, 0.05, 10000), rng.normal(0.53, 0.05, 10000)])
    found = ag.valley(x, *UNIT)
    assert found.threshold is None and "shallow" in found.why and found.depth > 0.7
    assert found.modes == pytest.approx((0.40, 0.53), abs=0.03)
    # and with the bar lowered it is where it should be
    assert ag.valley(x, *UNIT, max_depth=0.99).threshold == pytest.approx(0.465, abs=0.02)
    # a real second peak holding less than is asked for
    x, _ = two_peaks(n=40000, share=0.02, low=(0.25, 0.01), high=(0.65, 0.03))
    found = ag.valley(x, *UNIT, min_share=0.05)
    assert found.threshold is None and "would hold only 2.0%" in found.why


def test_valley_ignores_the_piles_at_the_ends_of_the_axis():
    x, positive = two_peaks()
    piled = np.concatenate([x, np.zeros(5000), np.ones(5000)])    # off scale, both ends
    assert ag.valley(piled, *UNIT).threshold == pytest.approx(ag.valley(x, *UNIT).threshold)


# --- a population on two parameters ----------------------------------------------------

def scatter(seed=0, n=(6000, 12000, 1500)):
    """Debris by the origin, cells, a smaller cluster — and where each is."""
    rng = np.random.default_rng(seed)
    debris = rng.multivariate_normal([0.10, 0.08], [[0.0010, 0], [0, 0.0006]], n[0])
    cells = rng.multivariate_normal([0.45, 0.30], [[0.0100, 0.0030], [0.0030, 0.0040]], n[1])
    other = rng.multivariate_normal([0.25, 0.80], [[0.0008, 0], [0, 0.0040]], n[2])
    noise = rng.uniform(0, 1, (800, 2))
    P = np.vstack([debris, cells, other, noise])
    label = np.repeat([0, 1, 2, 3], [n[0], n[1], n[2], 800])
    return P[:, 0], P[:, 1], label


def test_landscape_finds_the_populations():
    x, y, label = scatter()
    land = ag.Landscape(x, y, UNIT, UNIT)
    assert len(land.modes) == 3
    cells, debris, other = land.modes             # by how many events each holds
    assert (cells.x, cells.y) == pytest.approx((0.45, 0.30), abs=0.04)
    assert (debris.x, debris.y) == pytest.approx((0.10, 0.08), abs=0.02)
    assert (other.x, other.y) == pytest.approx((0.25, 0.80), abs=0.04)
    assert cells.share == pytest.approx(12000 / 20300, abs=0.05)
    assert sum(m.share for m in land.modes) == pytest.approx(1.0)
    # a point on a population chooses it; one on nothing, the nearest peak
    assert land.choose(None) == 0 and land.choose((0.5, 0.33)) == 0
    assert land.choose((0.11, 0.09)) == 1 and land.choose((0.26, 0.85)) == 2
    assert land.choose((0.02, 0.99)) == 2


def test_contour_of_the_largest_and_of_the_one_clicked():
    from milc.gating import points_in_polygon
    x, y, label = scatter()
    largest = ag.population(x, y, UNIT, UNIT)
    assert 6 <= len(largest.vertices) <= 32 and largest.populations == 3
    inside = points_in_polygon(x, y, largest.vertices)
    cells = label == 1
    # nine cells in ten, nearly none of the rest
    assert inside[cells].mean() == pytest.approx(0.90, abs=0.04)
    assert inside[label == 0].mean() < 0.01 and inside[label == 2].mean() < 0.01
    assert ag.jaccard(inside, cells) > 0.85
    assert "90% contour" in largest.why and "2 other populations" in largest.why
    assert ag.population(x, y, UNIT, UNIT, near=(0.5, 0.35)).vertices == pytest.approx(
        largest.vertices)

    debris = ag.population(x, y, UNIT, UNIT, near=(0.1, 0.1))
    inside = points_in_polygon(x, y, debris.vertices)
    assert ag.jaccard(inside, label == 0) > 0.85 and inside[cells].mean() < 0.01
    assert (debris.mode.x, debris.mode.y) == pytest.approx((0.10, 0.08), abs=0.02)

    # a smaller share is a tighter contour inside the first
    tight = ag.population(x, y, UNIT, UNIT, share=0.5)
    held = points_in_polygon(x, y, tight.vertices)
    assert held[cells].mean() == pytest.approx(0.5, abs=0.05)
    assert (held & ~points_in_polygon(x, y, largest.vertices)).sum() < 0.01 * held.sum()

    # the corner count is a limit, not a number to reach
    assert len(ag.population(x, y, UNIT, UNIT, max_vertices=8).vertices) <= 8


def test_contour_stops_at_the_valley_to_a_neighbour():
    """Two populations that touch: the contour of one does not run round both."""
    from milc.gating import points_in_polygon
    rng = np.random.default_rng(0)
    a = rng.multivariate_normal([0.35, 0.5], [[0.004, 0], [0, 0.004]], 8000)
    b = rng.multivariate_normal([0.60, 0.5], [[0.004, 0], [0, 0.004]], 8000)
    x, y = np.concatenate([a[:, 0], b[:, 0]]), np.concatenate([a[:, 1], b[:, 1]])
    found = ag.population(x, y, UNIT, UNIT, near=(0.33, 0.5))
    assert found.populations == 2
    assert found.vertices[:, 0].max() < 0.50          # the valley is at 0.475
    inside = points_in_polygon(x, y, found.vertices)
    assert inside[:8000].mean() > 0.8 and inside[8000:].mean() < 0.05


def test_nothing_to_outline():
    assert ag.population(np.array([]), np.array([]), UNIT, UNIT).vertices is None
    few = np.random.default_rng(0).uniform(0.2, 0.8, (50, 2))
    for found in (ag.population(few[:, 0], few[:, 1], UNIT, UNIT),
                  ag.gaussian(few[:, 0], few[:, 1], UNIT, UNIT)):
        assert "too few events" in found.why
    assert ag.gaussian(few[:, 0], few[:, 1], UNIT, UNIT).center is None
    # events on one line have no width for an ellipse
    line = np.linspace(0.2, 0.8, 5000)
    assert ag.gaussian(line, np.full(5000, 0.5), UNIT, UNIT).center is None


def test_ellipse_is_the_population_s_gaussian():
    x, y, label = scatter()
    found = ag.gaussian(x, y, UNIT, UNIT)
    assert found.center == pytest.approx((0.45, 0.30), abs=0.01)
    # the half-axes of the 95% ellipse of the covariance the cells were drawn with
    values, vectors = np.linalg.eigh(np.array([[0.0100, 0.0030], [0.0030, 0.0040]]))
    reach = np.sqrt(-2 * np.log(0.05))
    assert found.radii == pytest.approx((reach * values[1] ** 0.5, reach * values[0] ** 0.5),
                                        rel=0.08)
    expected = np.arctan2(vectors[1, 1], vectors[0, 1]) % np.pi
    assert found.angle % np.pi == pytest.approx(expected, abs=0.05)
    gate = EllipseGate(["x", "y"], list(found.center), list(found.radii), found.angle)
    inside = gate.contains({"x": x, "y": y}.__getitem__)
    assert inside[label == 1].mean() == pytest.approx(0.95, abs=0.02)
    assert inside[label == 0].mean() < 0.01 and inside[label == 2].mean() < 0.01
    # the debris, when it is the one clicked
    debris = ag.gaussian(x, y, UNIT, UNIT, near=(0.1, 0.1))
    assert debris.center == pytest.approx((0.10, 0.08), abs=0.01)
    assert debris.radii[0] == pytest.approx(reach * 0.0010 ** 0.5, rel=0.12)


# --- singlets -----------------------------------------------------------------------------

def pulses(seed=0, n=20000, doublets=0.12):
    """Height and area of a pulse: area 1.4 times height for a singlet, half
    as much again or more for a doublet."""
    rng = np.random.default_rng(seed)
    height = rng.uniform(30000, 120000, n)
    double = rng.random(n) < doublets
    area = 1.4 * height + rng.normal(0, 2500, n)
    area[double] = (1.4 * height * rng.uniform(1.25, 1.9, n))[double]
    return height, area, double


def test_band_holds_the_singlets():
    from milc.gating import points_in_polygon
    height, area, double = pulses()
    r = (0.0, 262144.0)
    found = ag.band(height, area, r, r)
    assert found.slope == pytest.approx(1.4, rel=0.01)
    assert found.half_width == pytest.approx(3 * 2500, rel=0.1)
    assert abs(found.intercept) < 1500 and found.vertices.shape == (4, 2)
    inside = points_in_polygon(height, area, found.vertices)
    assert inside[~double].mean() > 0.99 and inside[double].mean() < 0.03
    assert ag.jaccard(inside, ~double) > 0.98
    assert "3 standard deviations" in found.why
    # the other way round — height against area — the doublets fall below the line
    turned = ag.band(area, height, r, r)
    assert turned.slope == pytest.approx(1 / 1.4, rel=0.02)
    assert ag.jaccard(points_in_polygon(area, height, turned.vertices), ~double) > 0.97
    # narrower when asked
    narrow = ag.band(height, area, r, r, width=1.0)
    assert narrow.half_width == pytest.approx(found.half_width / 3, rel=1e-6)
    assert points_in_polygon(height, area, narrow.vertices)[~double].mean() == pytest.approx(
        0.68, abs=0.03)


def test_tight_band_for_a_line_with_a_fan_on_one_side():
    """Half the events spread above the line: the spread of both sides is
    then the fan's, and the line's own is on the side below."""
    from milc.gating import points_in_polygon
    rng = np.random.default_rng(1)
    n = 20000
    height = rng.uniform(30000, 120000, n)
    area = 1.4 * height + rng.normal(0, 2000, n)
    fanned = rng.random(n) < 0.45
    area[fanned] += np.abs(rng.normal(0, 15000, n))[fanned]
    r = (0.0, 262144.0)
    loose, tight = ag.band(height, area, r, r), ag.band(height, area, r, r, tight=True)
    assert tight.half_width == pytest.approx(3 * 2000, rel=0.25)
    assert tight.half_width < 0.6 * loose.half_width
    assert tight.slope == pytest.approx(1.4, rel=0.03)
    inside = points_in_polygon(height, area, tight.vertices)
    assert inside[~fanned].mean() > 0.97 and "cleaner side" in tight.why


def test_no_band_where_there_is_no_line():
    rng = np.random.default_rng(0)
    cloud = rng.normal(0.5, 0.1, (5000, 2))
    found = ag.band(cloud[:, 0], cloud[:, 1], UNIT, UNIT)
    assert found.vertices is None and "do not follow a line" in found.why
    assert "too few" in ag.band(cloud[:50, 0], cloud[:50, 1], UNIT, UNIT).why
    assert "same X" in ag.band(np.full(500, 0.5), cloud[:500, 1], UNIT, UNIT).why


# --- from a description to a gate ------------------------------------------------------------

def test_a_proposal_is_an_ordinary_gate_where_a_hand_drawn_one_would_be():
    biex = Biex()
    rng = np.random.default_rng(0)
    stored = {"M": np.concatenate([rng.normal(50, 60, 14000), rng.lognormal(8.5, 0.4, 6000)]),
              "F": rng.normal(90000, 15000, 20000)}
    stored["N"] = np.concatenate([rng.lognormal(8.0, 0.3, 14000), rng.normal(0, 80, 6000)])

    def shown(dim):
        return biex(stored[dim]) if dim != "F" else stored[dim]

    # a valley: a range open on the side kept, in the file's units
    spec = ag.make_spec(ag.VALLEY, ["M"], side=ag.ABOVE)
    assert spec == {"method": "valley", "dims": ["M"], "side": "above"}
    p = ag.propose(spec, shown, [(biex(-200.0), 1.0)], [biex])
    assert p.found and isinstance(p.gate, RectangleGate) and p.gate.dims == ["M"]
    assert p.gate.maxs == [None] and 400 < p.gate.mins[0] < 2500
    assert p.gate.contains(stored.__getitem__).sum() == pytest.approx(6000, abs=40)
    below = ag.propose(ag.make_spec(ag.VALLEY, ["M"], side=ag.BELOW), shown,
                       [(biex(-200.0), 1.0)], [biex]).gate
    assert below.mins == [None] and below.maxs == p.gate.mins
    assert p.shown[0] == pytest.approx(float(biex(p.gate.mins[0])))

    # a contour: a polygon in the file's units, tested on the scales of the plot
    ranges = [(float(biex(-200.0)), 1.0), (float(biex(-200.0)), 1.0)]
    p = ag.propose(ag.make_spec(ag.POPULATION, ["M", "N"]), shown, ranges, [biex, biex])
    assert isinstance(p.gate, PolygonGate) and p.gate.scales == [biex, biex]
    assert p.gate.vertices[:, 0].max() < 1000 and p.gate.vertices[:, 1].min() > 500
    assert p.gate.contains(stored.__getitem__)[:14000].mean() == pytest.approx(0.9, abs=0.04)
    assert p.gate.contains(stored.__getitem__)[14000:].mean() < 0.01
    assert p.near == pytest.approx((float(biex(50.0)), float(biex(np.exp(8.0)))), abs=0.03)
    assert np.allclose(p.shown[0], p.shown[-1])              # closed, to draw

    # an ellipse stays where it was drawn; a band is a polygon of four corners
    p = ag.propose(ag.make_spec(ag.ELLIPSE, ["M", "N"], near=(0.9, 0.2)), shown, ranges,
                   [biex, biex])
    assert isinstance(p.gate, EllipseGate) and p.gate.scales == [biex, biex]
    assert p.gate.contains(stored.__getitem__)[14000:].mean() == pytest.approx(0.95, abs=0.03)
    height, area, double = pulses()
    p = ag.propose(ag.make_spec(ag.BAND, ["H", "A"], tight=True),
                   {"H": height, "A": area}.__getitem__, [(0, 262144.0)] * 2, [None, None])
    assert isinstance(p.gate, PolygonGate) and p.gate.scales is None and p.near is None
    assert len(p.gate.vertices) == 4

    # nothing found: no gate, and the reason
    p = ag.propose(ag.make_spec(ag.VALLEY, ["F"]), shown, [(0.0, 262144.0)], [None])
    assert not p.found and p.gate is None and "one peak only" in p.why
    with pytest.raises(ValueError):
        ag.propose({"method": "magic", "dims": ["M"]}, shown, [UNIT], [None])


def test_jaccard():
    a = np.array([1, 1, 0, 0], bool)
    assert ag.jaccard(a, a) == 1 and ag.jaccard(a, ~a) == 0
    assert ag.jaccard(a, np.array([1, 0, 1, 0], bool)) == pytest.approx(1 / 3)
    assert ag.jaccard(np.zeros(3, bool), np.zeros(3, bool)) == 1


def test_how_a_gate_was_fitted_is_kept_with_it():
    spec = ag.make_spec(ag.POPULATION, ["A", "B"], near=np.array([0.5, 0.25]), share=0.9)
    tree = GateNode("")
    node = tree.add(GateNode("Cells", PolygonGate(["A", "B"], np.eye(2).repeat(2, 0)), fit=spec,
                             note=ag.note(spec)))
    node.add(GateNode("In", RectangleGate(["A"], [0.0], [1.0])))
    again = tree_from_dict(tree_to_dict(tree))
    assert again.children[0].fit == spec == {"method": "population", "dims": ["A", "B"],
                                             "near": [0.5, 0.25], "share": 0.9}
    assert again.children[0].children[0].fit is None
    assert "fit" not in tree_to_dict(tree)["children"][0]["children"][0]
    copy = copy_tree(tree)
    copy.children[0].fit["near"][0] = 9.0
    assert node.fit["near"] == [0.5, 0.25]

    # another sample's tree given to this one: its own fitted gate stays its own
    own = copy_tree(tree)
    own.children[0].gate.vertices += 5.0
    own.children[0].note = "its own note"
    given = copy_tree(tree)
    given.add(GateNode("New", RectangleGate(["A"], [0.0], [1.0]),
                       fit=ag.make_spec(ag.VALLEY, ["A"])))
    missing = keep_fits(given, own)
    assert [n.name for n in missing] == ["New"]
    assert np.array_equal(given.children[0].gate.vertices, own.children[0].gate.vertices)
    assert given.children[0].gate is not own.children[0].gate
    assert given.children[0].note == "its own note"
    # fitted another way since: that is not its own any more
    given = copy_tree(tree)
    given.children[0].fit["share"] = 0.5
    assert [n.name for n in keep_fits(given, own)] == ["Cells"]


# --- every sample its own ---------------------------------------------------------------------

def shifting(tmp_path, seed=0):
    """Three samples whose positive cells sit at 3000, 6000 and 12000, and a
    fourth with none; cells and debris on scatter, further out on the third."""
    from milc.fcs import write_fcs
    rng = np.random.default_rng(seed)
    truth = {}
    for name, bright, size in (("a1", 3000, 1.0), ("a2", 6000, 1.0), ("a3", 12000, 1.3),
                               ("a4", None, 1.0)):
        n = 6000
        cells = rng.random(n) < 0.7
        fsc = np.where(cells, rng.normal(90000 * size, 9000, n), rng.normal(20000, 4000, n))
        ssc = np.where(cells, rng.normal(40000 * size, 6000, n), rng.normal(9000, 3000, n))
        positive = cells & (rng.random(n) < 0.4) & (bright is not None)
        pe = np.where(positive, (bright or 1) * rng.lognormal(0, 0.25, n), rng.normal(60, 70, n))
        ev = np.column_stack([fsc, ssc, pe]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), ev, ["FSC-A", "SSC-A", "PE-A"],
                  markers=["", "", "CD4"], keywords={"SPILL": "1,PE-A,1"})
        truth[name] = (cells, positive)
    return str(tmp_path), truth


def test_a_gate_fitted_on_each_sample_separately(tmp_path):
    from milc.session import Session
    folder, truth = shifting(tmp_path)
    s = Session()
    s.open_folder(folder)
    first, second, third, fourth = s.usable()
    assert s.axis_range(first, "FSC-A") == (0.0, s.load(first).fcs.ranges[0])
    assert s.axis_range(first, "Comp-PE-A")[1] == pytest.approx(1.0, abs=0.01)

    # cells: the population at a click, on each sample its own
    spec = ag.make_spec(ag.POPULATION, ["FSC-A", "SSC-A"], near=(90000, 40000))
    proposal = s.fit_gate(first, "", spec)
    assert proposal.found and proposal.near == pytest.approx((90000, 40000), rel=0.1)
    cells = first.gates.add(GateNode("Cells", proposal.gate, fit=spec, note=ag.note(spec)))
    s.gates_changed(first)
    assert s.unfitted == []
    for e in s.usable():
        sample, (is_cell, _) = s.load(e), truth[e.name]
        node = e.gates.find("Cells")
        assert node.fit == spec and node.note.startswith("fitted to this sample's own events")
        assert ag.jaccard(sample.mask("Cells"), is_cell) > 0.85, e.name
    # the third sample's cells are further out, and so is its gate; copied
    # coordinates would have missed a third of them
    own, copied = third.gates.find("Cells").gate, first.gates.find("Cells").gate
    assert own.vertices[:, 0].mean() > 1.2 * copied.vertices[:, 0].mean()
    assert ag.jaccard(copied.contains(s.load(third).column), truth["a3"][0]) < 0.75

    # inside it, the valley of the marker: where each sample's own is
    spec = ag.make_spec(ag.VALLEY, ["Comp-PE-A"], side=ag.ABOVE)
    proposal = s.fit_gate(first, "Cells", spec)
    cells.add(GateNode("CD4+", proposal.gate, fit=spec, note=ag.note(spec)))
    s.gates_changed(first)
    cuts = [e.gates.find("Cells/CD4+").gate.mins[0] for e in (first, second, third)]
    assert cuts[0] < cuts[1] < cuts[2] and 300 < cuts[0] < 2000
    for e in (first, second, third):
        sample, (is_cell, positive) = s.load(e), truth[e.name]
        inside = sample.mask("Cells")
        assert ag.jaccard(sample.mask("Cells/CD4+"), positive & inside) > 0.98, e.name
    # the fourth has no positive cells: nothing is fitted there, it is said, and
    # the gate it is left with is the first sample's
    assert s.unfitted == [("a4", "Cells/CD4+", "one peak only: there is no valley to split at")]
    left = fourth.gates.find("Cells/CD4+")
    assert left.gate.mins == [cuts[0]] and left.note.startswith(ag.NOT_FITTED)
    assert "one peak only" in left.note and "a1" in left.note

    # a gate drawn by hand afterwards goes to every sample; the fitted ones stay
    cells.add(GateNode("Low", RectangleGate(["Comp-PE-A"], [None], [100.0])))
    s.gates_changed(first)
    assert s.unfitted == []
    assert [e.gates.find("Cells/CD4+").gate.mins[0] for e in (first, second, third)] == cuts
    assert all(e.gates.find("Cells/Low").gate.maxs == [100.0] for e in s.usable())
    assert fourth.gates.find("Cells/CD4+").note.startswith(ag.NOT_FITTED)
    # moving a fitted gate on one sample moves it there only
    second.gates.find("Cells/CD4+").gate.mins[0] = 5000.0
    s.gates_changed(second)
    assert first.gates.find("Cells/CD4+").gate.mins[0] == cuts[0]
    assert third.gates.find("Cells/CD4+").gate.mins[0] == cuts[2]
    # only this sample: nobody else is touched, nothing is fitted
    first.gates.find("Cells").add(GateNode("Here", RectangleGate(["FSC-A"], [0.0], [1.0])))
    s.gates_changed(first, everyone=False)
    assert s.unfitted == [] and "Cells/Here" not in [n.path for n in second.gates.walk()]

    # and all of it is there when the folder is opened again
    again = Session()
    again.open_folder(folder)
    a1, a2, a3, a4 = again.usable()
    assert a2.gates.find("Cells/CD4+").gate.mins == [5000.0]
    assert a3.gates.find("Cells/CD4+").gate.mins == [cuts[2]]
    assert a3.gates.find("Cells/CD4+").fit == spec
    assert a4.gates.find("Cells/CD4+").note.startswith(ag.NOT_FITTED)
    assert again.load(a3).count("Cells") == s.load(third).count("Cells")


def fitted_under_a_hand_drawn_gate(tmp_path):
    """Every sample: `Big`, drawn by hand around every event on FSC, and in
    it `Cells`, the largest population fitted on each sample's own events,
    and in that the valley of the marker, fitted on each sample too."""
    from milc.session import Session
    folder, truth = shifting(tmp_path)
    s = Session()
    s.open_folder(folder)
    first = s.usable()[0]
    first.gates.add(GateNode("Big", RectangleGate(["FSC-A"], [0.0], [1e6])))
    s.gates_changed(first)
    spec = ag.make_spec(ag.POPULATION, ["FSC-A", "SSC-A"])
    found = s.fit_gate(first, "Big", spec)
    spec["near"] = list(found.near)
    cells = first.gates.find("Big").add(GateNode("Cells", found.gate, fit=spec,
                                                 note=ag.note(spec)))
    marker = ag.make_spec(ag.VALLEY, ["Comp-PE-A"], side=ag.ABOVE)
    cells.add(GateNode("CD4+", s.fit_gate(first, "Big", marker).gate, fit=marker))
    s.gates_changed(first)
    return s, folder, truth


def centre(e, path="Big/Cells"):
    return float(e.gates.find(path).gate.vertices[:, 0].mean())


def test_fitted_gates_follow_their_parent_when_it_is_moved(tmp_path):
    from milc.session import Session
    s, folder, truth = fitted_under_a_hand_drawn_gate(tmp_path)
    first, second, third, fourth = s.usable()
    # on each sample the cells: where they are on that sample
    assert [round(centre(e), -4) for e in s.usable()] == [90000, 90000, 120000, 90000]
    before = {e.name: s.load(e).count("Big/Cells") for e in s.usable()}
    assert all(3500 < n < 4000 for n in before.values())
    # moved by hand on the second sample: that one is its own now
    second.gates.find("Big/Cells").gate.vertices[:, 0] += 5000.0
    s.gates_changed(second)
    moved = centre(second)
    # what is fitted inside it follows it, there only
    assert s.refitted == [("a2", "Big/Cells/CD4+")] and s.outdated() == []

    # the hand-drawn parent is narrowed to the debris: the cells fitted on each
    # sample are fitted again there — the largest population is now the debris
    first.gates.find("Big").gate.maxs = [45000.0]
    s.gates_changed(first)
    assert all(e.gates.find("Big").gate.maxs == [45000.0] for e in s.usable())
    assert sorted(s.refitted) == sorted(
        [(e.name, "Big/Cells") for e in (first, third, fourth)]
        + [(e.name, "Big/Cells/CD4+") for e in s.usable()])
    for e in (first, third, fourth):
        assert 15000 < centre(e) < 25000, e.name
        assert s.load(e).count("Big/Cells") < 0.5 * before[e.name]
        assert s.load(e).count("Big/Cells") > 0.8 * (s.load(e).count("Big"))
    # but not the one moved by hand: it is where it was, and said to be out of date
    assert centre(second) == moved
    late = s.outdated()
    assert [(e.name, n.path) for e, n in late] == [("a2", "Big/Cells")]

    # it stays out of date in the project
    again = Session()
    again.open_folder(folder)
    assert [(e.name, n.path) for e, n in again.outdated()] == [("a2", "Big/Cells")]
    assert centre(again.usable()[1]) == moved
    assert centre(again.usable()[0]) == centre(first)
    # and nothing else is fitted again for opening the folder and saving
    again.gates_changed(again.usable()[0])
    assert again.refitted == []

    # one click fits it again, giving up the hand's work, and its child with it
    again.refit_outdated()
    assert sorted(again.refitted) == [("a2", "Big/Cells"), ("a2", "Big/Cells/CD4+")]
    assert 15000 < centre(again.usable()[1]) < 25000 and again.outdated() == []
    # or it is kept as it is: no longer out of date, but still its own
    s.keep_outdated()
    assert s.outdated() == [] and centre(second) == moved
    third.gates.find("Big").gate.maxs = [1e6]
    s.gates_changed(third)
    assert centre(second) == moved and [n.path for _, n in s.outdated()] == ["Big/Cells"]
    assert 80000 < centre(first) < 100000


def test_fitted_gates_follow_the_values_under_them(tmp_path):
    from milc.compensation import Compensation
    from milc.transforms import Arcsinh
    s, folder, truth = fitted_under_a_hand_drawn_gate(tmp_path)
    first = s.usable()[0]
    cut = first.gates.find("Big/Cells/CD4+").gate.mins[0]
    s.load(first)
    # values halved by another compensation: the valley is fitted again on them
    s.set_matrix(Compensation(["PE-A"], np.array([[2.0]]), "edited"))
    assert ("a1", "Big/Cells/CD4+") in s.refitted and ("a1", "Big/Cells") not in s.refitted
    assert first.gates.find("Big/Cells/CD4+").gate.mins[0] == pytest.approx(cut / 2, rel=0.15)
    assert s.review == "the compensation changed"   # hand-drawn gates are still to be looked at
    # another scale on the axis: the same
    s.set_scale(["PE-A"], Arcsinh(150.0))
    assert ("a1", "Big/Cells/CD4+") in s.refitted
    s.set_scale(["PE-A"], Arcsinh(150.0))
    assert s.refitted == []


def test_a_gate_fitted_on_one_sample_and_given_to_all(tmp_path):
    from milc.session import Session
    folder, truth = shifting(tmp_path)
    s = Session()
    s.open_folder(folder)
    first, second, third, fourth = s.usable()
    first.gates.add(GateNode("Big", RectangleGate(["FSC-A"], [0.0], [1e6])))
    s.gates_changed(first)
    spec = {**ag.make_spec(ag.POPULATION, ["FSC-A", "SSC-A"]), "sample": "a3.fcs"}
    first.gates.find("Big").add(GateNode("Cells", s.fit_gate(third, "Big", spec).gate,
                                         fit=spec, note=ag.note(spec)))
    s.gates_changed(first)
    assert s.unfitted == [] and "fitted to the events of a3.fcs" in ag.note(spec)
    gates = [e.gates.find("Big/Cells").gate.vertices for e in s.usable()]
    assert all(np.array_equal(v, gates[0]) for v in gates)    # the same gate everywhere
    # the parent moved on the first sample: fitted again on the third's events, given to all
    first.gates.find("Big").gate.maxs = [45000.0]
    s.gates_changed(first)
    assert {name for name, _ in s.refitted} == {"a1", "a2", "a3", "a4"}
    gates = [e.gates.find("Big/Cells").gate.vertices for e in s.usable()]
    assert all(np.array_equal(v, gates[0]) for v in gates)
    assert np.array_equal(gates[0], s.fit_gate(third, "Big", spec).gate.vertices)
    assert 15000 < centre(first) < 25000
    # moved on one sample only: there it is its own, fitted again on that one's
    fourth.gates.find("Big").gate.maxs = [1e6]
    s.gates_changed(fourth, everyone=False)
    assert s.refitted == [("a4", "Big/Cells")] and 80000 < centre(fourth) < 100000
    assert 15000 < centre(third) < 25000


def test_renaming_keeps_what_each_sample_has_of_its_own(tmp_path):
    from milc.gating import BooleanGate
    from milc.session import Session
    s, folder, truth = fitted_under_a_hand_drawn_gate(tmp_path)
    first, second, third, fourth = s.usable()
    cells = first.gates.find("Big/Cells")
    cells.add(GateNode("Low", RectangleGate(["Comp-PE-A"], [None], [100.0])))
    cells.add(GateNode("Not CD4", BooleanGate("NOT", ["CD4+"])))
    s.gates_changed(first)
    own = {e.name: e.gates.find("Big/Cells/CD4+").gate.mins[0] for e in s.usable()}
    assert own["a1"] < own["a2"] < own["a3"]
    second.gates.find("Big/Cells/CD4+").gate.mins = [5000.0]      # by hand, there only
    s.gates_changed(second)
    shapes = {e.name: centre(e) for e in s.usable()}

    s.rename_gate(first, "Big/Cells/CD4+", "CD4 pos")
    assert s.refitted == [] and s.outdated() == []
    for e in s.usable():
        node = e.gates.find("Big/Cells/CD4 pos")
        assert node.fit["method"] == "valley" and node.fitted is not None
        assert node.gate.mins[0] == (5000.0 if e.name == "a2" else own[e.name]), e.name
        assert e.gates.find("Big/Cells/Not CD4").gate.operands == ["CD4 pos"]
        assert "Big/Cells/CD4+" not in [n.path for n in e.gates.walk()]
    # a fitted parent renamed: each sample's own outline stays, and so do its children
    s.rename_gate(second, "Big/Cells", "Lymphocytes")
    assert s.refitted == []
    assert {e.name: centre(e, "Big/Lymphocytes") for e in s.usable()} == shapes
    assert second.gates.find("Big/Lymphocytes/CD4 pos").gate.mins == [5000.0]

    # and after the folder is opened again
    again = Session()
    again.open_folder(folder)
    for e in again.usable():
        assert e.gates.find("Big/Lymphocytes/CD4 pos").gate.mins[0] == \
            (5000.0 if e.name == "a2" else own[e.name])
        assert centre(e, "Big/Lymphocytes") == shapes[e.name]
    assert again.outdated() == []
    # the moved one still counts as moved: a parent change leaves it alone
    a1 = again.usable()[0]
    a1.gates.find("Big").gate.maxs = [2e5]
    again.gates_changed(a1)
    assert again.usable()[1].gates.find("Big/Lymphocytes/CD4 pos").gate.mins == [5000.0]


# --- in the window ----------------------------------------------------------------------------

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    QtWidgets = pytest.importorskip("PyQt6.QtWidgets")
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


def test_automatic_gates_in_the_window(app, tmp_path, monkeypatch):
    from PyQt6 import QtWidgets
    from milc.ui import worker
    from milc.ui.gates_tab import HISTOGRAM
    from milc.ui.shell import MainShell

    folder, truth = shifting(tmp_path)
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
    assert g.entry.name == "a1"
    # one button, a menu of what can be proposed, and how it goes to the others
    texts = [a.text() for a in g.auto_button.menu().actions() if a.text()]
    assert len(texts) == 9 and texts[0].startswith("Valley") and texts[-1] == \
        "Fit on each sample separately"
    assert g.per_sample.isCheckable() and not g.per_sample.isChecked()

    # a population: proposed first — an outline, not yet a gate — then accepted
    g.set_axes("FSC-A", "SSC-A")
    found = g.propose_gate(ag.POPULATION, near=(88000, 42000))
    assert found.found and g.proposal is not None and g.box._shape is not None
    assert g.node.children == [] and w.session.populations() == []
    cells = g.accept_proposal("Cells")
    assert g.proposal is None and g.box._shape is None
    # the same gate for every sample, but how it was fitted, and on which, is kept
    assert isinstance(cells.gate, PolygonGate) and cells.fit["sample"] == "a1.fcs"
    assert cells.note.startswith("fitted to the events of a1.fcs")
    assert ag.jaccard(s.mask("Cells"), truth["a1"][0]) > 0.85
    # like any polygon: on the plot with its share, on every sample the same, to be dragged
    child, roi, label = next(o for o in g._overlays if o[0] is cells)
    assert label.toPlainText() == f"Cells  {100 * s.count('Cells') / len(s):.1f}%"
    others = [e.gates.find("Cells").gate.vertices for e in w.session.usable()]
    assert all(np.array_equal(v, cells.gate.vertices) for v in others)
    before = g.sample.count("Cells")
    roi.setPos(roi.pos().x() + 30000, roi.pos().y())
    g._gate_moved(child, roi)
    assert g.sample.count("Cells") < before
    assert "90% contour" in " ".join(said)

    # the largest population without a click, as an ellipse; and the band of a diagonal
    ellipse = g.auto_gate(ag.ELLIPSE)
    assert ellipse.name == "Population" and isinstance(ellipse.gate, EllipseGate)
    assert ag.jaccard(g.sample.mask("Population"), truth["a1"][0]) > 0.85
    assert g.auto_gate(ag.BAND) is None and "do not follow a line" in g.auto_said

    # a valley, on a histogram: the positive side, fitted on each sample separately
    g.select("Population")
    g.set_axes("Comp-PE-A", HISTOGRAM)
    g.per_sample.setChecked(True)
    cd4 = g.auto_gate(ag.VALLEY, side=ag.ABOVE)
    assert cd4.name == "CD4+" and cd4.gate.dims == ["Comp-PE-A"] and cd4.gate.maxs == [None]
    assert cd4.fit["method"] == "valley" and cd4.note.startswith("fitted to this sample")
    assert "except a4" in said[-1] and "one peak only" in said[-1]
    cuts = [e.gates.find("Population/CD4+").gate.mins[0] for e in w.session.usable()]
    assert cuts[0] < cuts[1] < cuts[2] and cuts[3] == cuts[0]
    assert ag.jaccard(g.sample.mask("Population/CD4+"),
                      truth["a1"][1] & g.sample.mask("Population")) > 0.98
    # in the tree it is set apart, and says why
    item = g.tree.currentItem().child(0)
    assert item.text(0) == "CD4+" and item.font(0).italic()
    assert "each sample has its own" in item.toolTip(0)
    # the negative side takes the other name; the same coordinates everywhere when not asked
    g.per_sample.setChecked(False)
    negative = g.auto_gate(ag.VALLEY, side=ag.BELOW)
    assert negative.name == "CD4-" and negative.fit["sample"] == "a1.fcs"
    assert negative.gate.maxs == cd4.gate.mins
    assert len({e.gates.find("Population/CD4-").gate.maxs[0] for e in w.session.usable()}) == 1
    assert [e.gates.find("Population/CD4+").gate.mins[0] for e in w.session.usable()] == cuts
    # only on this sample: nothing to fit elsewhere, whatever the menu says
    g.per_sample.setChecked(True)
    g.everyone.setChecked(False)
    here = g.auto_gate(ag.VALLEY, "Here", side=ag.ABOVE)
    assert here.fit["sample"] == "a1.fcs" and "Population/Here" not in w.session.populations()
    g.everyone.setChecked(True)
    g.per_sample.setChecked(False)

    # where nothing can be proposed it is said, and no gate is made
    g.select("Population/CD4+")
    count = len(g.node.children)
    assert g.auto_gate(ag.VALLEY, side=ag.ABOVE) is None
    assert g.auto_said == "one peak only: there is no valley to split at"
    assert len(g.node.children) == count and g.proposal is None
    assert g.auto_gate(ag.POPULATION) is None and "two parameters" in g.auto_said

    # from the menu: the reason in a box; a proposal with a name to accept or not
    shown = []
    monkeypatch.setattr(QtWidgets.QMessageBox, "information",
                        lambda *a, **k: shown.append(a[2]))
    g._ask_auto(ag.VALLEY, side=ag.ABOVE)
    assert shown == ["One peak only: there is no valley to split at."]
    g.select("")
    g.set_axes("FSC-A", "SSC-A")
    asked = []
    monkeypatch.setattr(QtWidgets.QInputDialog, "getText",
                        lambda *a, **k: (asked.append((a[2], k.get("text"))), ("", False))[1])
    count = len(g.node.children)
    g._ask_auto(ag.POPULATION)
    assert asked[0][0].startswith("Proposed: the 90% contour") and asked[0][1] == "Population 2"
    assert len(g.node.children) == count and g.box._shape is None     # turned down
    # "at a click": the next click on the plot says which population
    g._seed(ag.POPULATION)
    assert g.box.mode == "seed" and "Click on the population" in said[-1]
    monkeypatch.setattr(QtWidgets.QInputDialog, "getText", lambda *a, **k: ("Debris", True))
    g.box.sigPoint.emit(20000.0, 9000.0)
    assert g.box.mode is None
    debris = g.entry.gates.find("Debris")
    assert ag.jaccard(g.sample.mask("Debris"), ~truth["a1"][0]) > 0.85
    assert debris.gate.vertices[:, 0].max() < 60000
    # a click of the mouse itself, on the events: their dots take the clicks on
    # them, and the tool is given it all the same
    from PyQt6 import QtCore
    from PyQt6.QtTest import QTest
    g._seed(ag.ELLIPSE)
    spot = g.plot.mapFromScene(g.box.mapViewToScene(QtCore.QPointF(90000.0, 40000.0)))
    QTest.mouseClick(g.plot.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=spot)
    app.processEvents()
    clicked = g.entry.gates.find("Debris 2")
    assert g.box.mode is None and isinstance(clicked.gate, EllipseGate)
    assert clicked.gate.center == pytest.approx([90000, 40000], rel=0.05)
    # a click is not asked for on a histogram
    g.set_axes("Comp-PE-A", HISTOGRAM)
    g._seed(ag.ELLIPSE)
    assert g.box.mode is None and "two parameters" in said[-1]
    # and the quadrant tool still takes its click
    g.set_axes("FSC-A", "SSC-A")
    g.quad_button.setChecked(True)
    g.box.sigPoint.emit(50000.0, 20000.0)
    assert len(g._crossings) == 1 and not g.quad_button.isChecked()
    assert w.minimumSizeHint().width() <= 1280
    w.close()


def test_fitted_gates_refit_and_rename_in_the_window(app, tmp_path):
    from milc.ui import worker
    from milc.ui.shell import MainShell

    folder, truth = shifting(tmp_path)
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
    s = w.session
    a1, a2, a3, a4 = s.usable()
    g.set_axes("FSC-A", "SSC-A")
    g.add_rectangle(0.0, 0.0, 250000.0, 250000.0, "Big")
    g.select("Big")
    g.per_sample.setChecked(True)
    cells = g.auto_gate(ag.POPULATION, "Cells")
    assert cells.fit and "sample" not in cells.fit
    assert not g.outdated_bar.isVisibleTo(g)

    # moved by hand on the second sample
    g.samples.setCurrentRow(1)
    app.processEvents()
    assert g.entry is a2
    g.select("Big")
    child, roi, _ = next(o for o in g._overlays if o[0].name == "Cells")
    roi.setPos(roi.pos().x() + 5000, roi.pos().y())
    g._gate_moved(child, roi)
    moved = centre(a2)
    # renamed on the first: the second keeps its own
    g.samples.setCurrentRow(0)
    app.processEvents()
    g.select("Big/Cells")
    g.rename("Lymphocytes")
    assert centre(a2, "Big/Lymphocytes") == moved and s.refitted == []

    # the parent narrowed: fitted again where not moved by hand, and said
    g.select("Big")
    big = a1.gates.find("Big")
    big.gate.mins, big.gate.maxs = [0.0, 0.0], [45000.0, 250000.0]
    g._commit()
    assert "Fitted again on the population as it now is: Big/Lymphocytes, on 3 samples" \
        in said[-1]
    assert 15000 < centre(a1, "Big/Lymphocytes") < 25000
    assert centre(a2, "Big/Lymphocytes") == moved
    assert g.outdated_bar.isVisibleTo(g)
    assert "Big/Lymphocytes was fitted automatically and then moved by hand on a2" \
        in g.outdated.text()
    w.resize(1052, 760)
    app.processEvents()
    g.grab().save(str(tmp_path / "outdated.png"))
    if os.environ.get("MILC_SHOTS"):
        g.grab().save(os.path.join(os.environ["MILC_SHOTS"], "gates_outdated.png"))
    assert w.minimumSizeHint().width() <= 1052, _asking(w)
    # one click: fitted again there too
    g.refit_button.click()
    assert 15000 < centre(a2, "Big/Lymphocytes") < 25000
    assert not g.outdated_bar.isVisibleTo(g) and s.outdated() == []
    w.close()


def _asking(widget) -> list:
    """For the message of a failure: the visible widgets that ask for the
    most width, widest first, with their text."""
    from PyQt6 import QtWidgets
    found = []
    for child in widget.findChildren(QtWidgets.QWidget):
        if child.isVisible() and child.minimumSizeHint().width() > 200:
            text = getattr(child, "text", lambda: "")()
            found.append((child.minimumSizeHint().width(), type(child).__name__,
                          str(text)[:30]))
    return sorted(found, reverse=True)[:15]
