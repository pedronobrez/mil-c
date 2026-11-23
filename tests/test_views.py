"""Violins and boxes, the matrix of every pair, backgating: the numbers
against an independent calculation, the figures as files, and the Gates
tab driven without a screen."""

import os

import numpy as np
import pytest

from milc import density, figures, views
from milc.gating import (NOT, BooleanGate, EllipseGate, GateNode, PolygonGate, QuadrantGate,
                         RectangleGate)
from milc.transforms import Arcsinh, Biex, decade_ticks

STYLE = figures.FigureStyle(font="DejaVu Sans")


def _svg(content, tmp_path, style=STYLE) -> str:
    path = str(tmp_path / "figure.svg")
    figures.save(content, style, path)
    with open(path, encoding="utf-8") as f:
        return f.read()


# --- violins and boxes -----------------------------------------------------------

def test_summary_is_numpy_percentiles_and_tukeys_whiskers():
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(0.4, 0.05, 5000), rng.uniform(0, 1, 300)])
    s = views.summary(x)
    q1, median, q3 = np.percentile(x, 25), np.median(x), np.percentile(x, 75)
    assert (s.n, s.median, s.q1, s.q3) == (5300, pytest.approx(median), pytest.approx(q1),
                                           pytest.approx(q3))
    # the whiskers, the slow way: the furthest events within 1.5 IQR of the box
    ordered = sorted(x.tolist())
    low = next(v for v in ordered if v >= q1 - 1.5 * (q3 - q1))
    high = next(v for v in reversed(ordered) if v <= q3 + 1.5 * (q3 - q1))
    assert (s.low, s.high) == (low, high)
    assert s.outliers == sum(v < low or v > high for v in ordered) > 200

    # an event exactly at the limit of the rule is inside it
    edge = views.summary(np.array([0.0, 1.0, 2.0, 3.0, 4.0, 7.0]))
    assert (edge.q1, edge.q3, edge.high, edge.outliers) == (1.25, 3.75, 7.0, 0)
    assert views.summary(np.array([0.0, 1.0, 2.0, 3.0, 4.0, 7.6])).outliers == 1

    one = views.summary(np.array([3.0]))
    assert (one.median, one.low, one.high, one.outliers) == (3.0, 3.0, 3.0, 0)
    empty = views.summary(np.array([]))
    assert empty.n == 0 and np.isnan(empty.median) and empty.outliers == 0


def test_violin_is_the_overlay_histogram_and_outliers_are_capped():
    rng = np.random.default_rng(1)
    x = np.concatenate([rng.normal(0.3, 0.03, 20000), rng.uniform(0, 1, 2000)])
    d = views.distribution(x, 0.0, 1.0, "s1", "Case")
    centres, heights = density.histogram(x, 0.0, 1.0)
    assert np.array_equal(d.x, centres) and np.array_equal(d.density, heights)
    at, half = views.violin(d)
    keep = np.isin(centres, at)
    assert np.allclose(half, 0.42 * heights[keep] / 100) and half.max() == pytest.approx(0.42)
    assert np.all(heights[~keep] <= 0.5)       # only the hair-thin tails are left out

    s = d.summary
    far = x[(x < s.low) | (x > s.high)]
    assert s.outliers == far.size > views.MAX_OUTLIERS == d.outliers.size
    assert d.outliers.min() == far.min() and d.outliers.max() == far.max()
    assert np.isin(d.outliers, far).all()
    few = views.distribution(x[:300], 0.0, 1.0, "s2")
    assert few.group == "s2" and few.outliers.size == few.summary.outliers

    nothing = views.distribution(np.array([]), 0.0, 1.0, "none")
    assert nothing.summary.n == 0 and views.violin(nothing)[0].size == 0

    content = views.DistributionContent([d, few], views.BOX, "CD4")
    drawn = d.outliers.size + few.outliers.size
    assert f"{drawn:,d} of {s.outliers + few.summary.outliers:,d} outliers drawn" \
        in content.note()
    assert "1.5 IQR" in content.note()
    assert "outliers" not in views.DistributionContent([few], views.BOX).note()
    assert "median" in views.DistributionContent([d], views.VIOLIN).note()


def test_distributions_as_figures(tmp_path):
    rng = np.random.default_rng(2)
    b = Biex()
    series = [views.distribution(b(rng.lognormal(6 + k, 1.0, 4000)), 0.0, 1.0, name, group)
              for k, (name, group) in enumerate((("3Ctrl1", "Ctrl"), ("3Ctrl2", "Ctrl"),
                                                 ("3Case1", "Case"), ("nobody", "Case")))]
    series[-1] = views.distribution(np.array([]), 0.0, 1.0, "nobody", "Case")
    for kind in (views.VIOLIN, views.BOX):
        c = views.DistributionContent(series, kind, "CD4 PE", decade_ticks(b), (0.0, 1.0))
        assert figures.category_names(c) == [("Ctrl", figures.CATEGORICAL[0]),
                                             ("Case", figures.CATEGORICAL[1])]
        assert not figures.uses_colormap(c)
        assert figures.preview_png(c, STYLE, dpi=40)[:4] == b"\x89PNG"
        look = figures.FigureStyle(font="DejaVu Sans", labels={"Case": "Treated"},
                                   colors={"Ctrl": "#123456"})
        svg = _svg(c, tmp_path, look)
        for text in ("3Ctrl1", "3Case1", "nobody", "CD4 PE", "Treated", "Ctrl"):
            assert f">{text}<" in svg
        assert "#123456" in svg and "<image" not in svg
    # a colour for each row: the axis names them, and there is no legend
    alone = views.DistributionContent(series[:1], views.VIOLIN, "CD4 PE")
    assert _svg(alone, tmp_path).count(">3Ctrl1<") == 1


# --- every pair -------------------------------------------------------------------

def test_pairs_are_symmetric_and_count_every_event():
    rng = np.random.default_rng(3)
    n, k, bins = 5000, 4, 16
    columns = [rng.normal(0.5, 0.25, n) for _ in range(k)]   # some beyond both ends
    columns[1] = 0.5 * columns[0] + 0.5 * columns[1]
    ranges = [(0.0, 1.0), (0.0, 1.0), (-0.5, 2.0), (0.2, 0.9)]
    counts = views.pairs(columns, ranges, bins=bins)
    assert counts.shape == (k, k, bins, bins)
    assert (counts.sum(axis=(2, 3)) == n).all()
    clipped = [np.clip(x, lo, hi) for x, (lo, hi) in zip(columns, ranges)]
    for i in range(k):
        h, _ = np.histogram(clipped[i], bins=bins, range=ranges[i])
        assert np.array_equal(np.diagonal(counts[i, i]), h)
        assert counts[i, i].sum() == h.sum()               # and nothing off its diagonal
        for j in range(k):
            assert np.array_equal(counts[i, j], counts[j, i].T)
            H, _, _ = np.histogram2d(clipped[i], clipped[j], bins=bins,
                                     range=[ranges[i], ranges[j]])
            assert np.array_equal(counts[i, j], H)
    assert views.pairs(columns, ranges).shape[-1] == views.MATRIX_BINS == 64


def test_matrix_as_a_figure(tmp_path):
    rng = np.random.default_rng(4)
    b = Biex()
    raw = [rng.lognormal(6, 1.2, 3000) for _ in range(3)]
    scatter = rng.normal(90000, 9000, 3000)
    ranges = [(0.0, 1.0)] * 3 + [(0.0, 262144.0)]
    ticks = [decade_ticks(b)] * 3 + [views.plain_ticks(0.0, 262144.0)]
    assert ticks[-1] == [(0.0, "0"), (100000.0, "100K"), (200000.0, "200K")]
    c = views.MatrixContent(["CD8 FITC", "CD4 PE", "Foxp3 APC", "FSC-A"],
                            views.pairs([b(x) for x in raw] + [scatter], ranges),
                            ranges, ticks)
    assert figures.uses_colormap(c) and figures.category_names(c) == []
    assert figures.preview_png(c, STYLE, dpi=40)[:4] == b"\x89PNG"
    svg = _svg(c, tmp_path)
    for text in ("CD8", "CD4", "Foxp3", "FSC-A", "100K"):
        assert f">{text}<" in svg
    # every density cell is in one bitmap, whatever the number of events
    assert svg.count("<image") == 1
    more = views.MatrixContent(c.labels, c.counts * 1000, ranges, ticks)
    assert abs(len(_svg(more, tmp_path)) - len(svg)) < 0.05 * len(svg)
    # cells too small for the decades to be written keep their marks only
    small = figures.FigureStyle(font="DejaVu Sans", width_mm=60, height_mm=60)
    assert ">100K<" not in _svg(c, tmp_path, small)


# --- backgating -------------------------------------------------------------------

def _tree():
    """Every kind of gate, one inside the other, on events laid out so
    that each keeps some and loses some."""
    rng = np.random.default_rng(5)
    n = 6000
    events = {"FSC": rng.uniform(0, 1000, n), "SSC": rng.uniform(0, 1000, n),
              "A": rng.lognormal(5, 1.5, n), "B": rng.lognormal(5, 1.5, n),
              "C": rng.uniform(0, 100, n), "D": rng.uniform(0, 100, n)}
    scale = Arcsinh(150.0)
    tree = GateNode("")
    cells = tree.add(GateNode("Cells", RectangleGate(["FSC", "SSC"], [200.0, None],
                                                     [None, 800.0])))
    square = np.array([[20.0, 20.0], [2000.0, 20.0], [2000.0, 2000.0], [20.0, 2000.0]])
    poly = cells.add(GateNode("Poly", PolygonGate(["A", "B"], square, scales=[scale, scale])))
    cells.add(GateNode("Other", RectangleGate(["A", "B"], [0.0, 0.0], [10.0, 10.0])))
    round_ = poly.add(GateNode("Round", EllipseGate(["C", "D"], [50.0, 50.0], [40.0, 25.0], 0.5)))
    corner = round_.add(GateNode("Q2", QuadrantGate(["C", "D"], 45.0, 50.0, "Q2", "g")))
    corner.add(GateNode("Low", RectangleGate(["C"], [None], [70.0])))
    high = corner.add(GateNode("High", RectangleGate(["C"], [70.0], [None])))
    corner.add(GateNode("Not high", BooleanGate(NOT, ["High"])))
    high.add(GateNode("Nothing", BooleanGate(NOT, ["Missing"])))
    return tree, events, scale


def test_backgate_steps_on_a_tree_of_every_gate():
    tree, events, _ = _tree()
    masks = tree.apply(lambda d: events[d], 6000)
    path = "Cells/Poly/Round/Q2/Not high"
    found = views.steps(tree, path, masks)
    assert [s.name for s in found] == ["Cells", "Poly", "Round", "Q2", "Not high"]
    assert [s.path for s in found] == ["Cells", "Cells/Poly", "Cells/Poly/Round",
                                       "Cells/Poly/Round/Q2", path]
    # a combination is looked at on what it is made of
    assert [s.dims for s in found] == [["FSC", "SSC"], ["A", "B"], ["C", "D"], ["C", "D"], ["C"]]
    parents = ["", "Cells", "Cells/Poly", "Cells/Poly/Round", "Cells/Poly/Round/Q2"]
    for step, parent in zip(found, parents):
        assert np.array_equal(step.parent, masks[parent])
        assert 0 < masks[path].sum() < step.parent.sum()
        assert not (masks[path] & ~step.parent).any()   # it is inside every one of them
    assert views.steps(tree, "", masks) == []
    assert [s.name for s in views.steps(tree, "Cells/Other", masks)] == ["Cells", "Other"]
    # made of nothing that is there: a step all the same, with nothing to draw it on
    assert views.steps(tree, "Cells/Poly/Round/Q2/High/Nothing", masks)[-1].dims == []
    with pytest.raises(KeyError):
        views.steps(tree, "Cells/Nobody", masks)


def test_gate_outlines():
    tree, _, drawn = _tree()
    shown = Biex()
    plain = lambda dim: None                                        # noqa: E731
    cells = tree.find("Cells").gate
    (box,) = views.outline(cells, plain, [(0.0, 1000.0), (-50.0, 1200.0)])
    # the open sides are at the ends of the axes
    assert box.tolist() == [[200, -50], [1000, -50], [1000, 800], [200, 800], [200, -50]]

    poly = tree.find("Cells/Poly").gate
    (same,) = views.outline(poly, lambda dim: drawn, [(0, 1), (0, 1)])
    assert np.allclose(same, drawn(np.vstack([poly.vertices, poly.vertices[:1]])))
    # on another scale its sides are curves: every point of the outline is
    # still on a side of the polygon where it is tested
    (curve,) = views.outline(poly, lambda dim: shown, [(0, 1), (0, 1)])
    back = drawn(shown.inverse(curve))
    lo, hi = drawn(20.0), drawn(2000.0)
    on_side = (np.isclose(back, lo, atol=1e-6) | np.isclose(back, hi, atol=1e-6)).any(axis=1)
    assert len(curve) > 50 and on_side.all()
    assert (back >= lo - 1e-6).all() and (back <= hi + 1e-6).all()
    # a slanted side, drawn on plain axes and shown with one of them
    # arcsinh, is the curve y = sinh(x): not the chord between its ends
    slanted = PolygonGate(["A", "B"], np.array([[0.0, 0.0], [900.0, 900.0], [900.0, 0.0]]))
    bent = Arcsinh(150.0)
    (curve,) = views.outline(slanted, {"A": bent, "B": None}.get, [(0, 3), (0, 1000)])
    x, y = bent.inverse(curve[:, 0]), curve[:, 1]
    assert (np.isclose(x, y) | np.isclose(x, 900.0) | np.isclose(y, 0.0)).all()
    assert np.isclose(x, 450.0).any() and np.isclose(x, y)[:24].all()

    round_ = tree.find("Cells/Poly/Round").gate
    for scale in (None, Arcsinh(5.0)):
        (line,) = views.outline(round_, lambda dim: scale, [(0, 1), (0, 1)])
        stored = line if scale is None else scale.inverse(line)
        dx, dy = stored[:, 0] - 50.0, stored[:, 1] - 50.0
        u = dx * np.cos(0.5) + dy * np.sin(0.5)
        v = -dx * np.sin(0.5) + dy * np.cos(0.5)
        assert np.allclose((u / 40.0) ** 2 + (v / 25.0) ** 2, 1.0)

    quadrant = tree.find("Cells/Poly/Round/Q2").gate
    vertical, horizontal = views.outline(quadrant, plain, [(0.0, 100.0), (-10.0, 110.0)])
    assert vertical.tolist() == [[45, -10], [45, 110]]
    assert horizontal.tolist() == [[0, 50], [100, 50]]

    low = tree.find("Cells/Poly/Round/Q2/Low").gate
    (bound,) = views.outline(low, plain, [(0.0, 100.0)])
    assert bound.tolist() == [[70, 0], [70, 1]]
    both = views.outline(RectangleGate(["C"], [10.0], [70.0]), lambda dim: Arcsinh(5.0),
                         [(0.0, 5.0)])
    assert [line[0, 0] for line in both] == pytest.approx([np.arcsinh(2), np.arcsinh(14)])
    assert views.outline(tree.find("Cells/Poly/Round/Q2/Not high").gate, plain, []) == []


def test_backgate_content_and_figure(tmp_path):
    tree, events, drawn = _tree()
    masks = tree.apply(lambda d: events[d], 6000)
    scales = {"A": drawn, "B": drawn}
    ranges = {"FSC": (0.0, 1000.0), "SSC": (0.0, 1000.0), "C": (0.0, 100.0), "D": (0.0, 100.0)}

    def display(dim):
        return scales[dim](events[dim]) if dim in scales else events[dim]

    path = "Cells/Poly/Round/Q2/Not high"
    c = views.backgate(tree, path, masks, display, scales.get,
                       lambda d: ranges.get(d, (0.0, 8.0)), lambda d: f"{d} marker")
    chosen = int(masks[path].sum())
    assert (c.population, c.events, c.drawn) == ("Not high", chosen, chosen)
    assert [p.title for p in c.panels] == ["Cells", "Poly", "Round", "Q2", "Not high"]
    for panel, parent in zip(c.panels[:4], ["", "Cells", "Cells/Poly", "Cells/Poly/Round"]):
        assert panel.counts.sum() == masks[parent].sum()       # every event of the parent
        assert panel.points.shape == (chosen, 2) and len(panel.outline) >= 1
    first = c.panels[0]
    H, _, _ = np.histogram2d(events["SSC"][masks[""]], events["FSC"][masks[""]],
                             bins=views.BACKGATE_BINS, range=[(0, 1000), (0, 1000)])
    assert np.array_equal(first.counts, H)                     # [y, x]
    assert (first.xlabel, first.ylabel, first.xticks) == ("FSC marker", "SSC marker", None)
    assert np.array_equal(first.points[:, 0], events["FSC"][masks[path]])
    assert c.panels[1].xticks == decade_ticks(drawn) and c.panels[1].yrange == (0.0, 8.0)
    last = c.panels[-1]
    assert last.counts is None and last.ylabel == "Events"
    # one parameter: the population's histogram sits inside its parent's
    assert last.parent.sum() == pytest.approx(masks["Cells/Poly/Round/Q2"].sum())
    assert last.chosen.sum() == pytest.approx(chosen) and (last.chosen <= last.parent + 1e-9).all()
    assert last.outline == []                                   # a combination has none

    assert figures.category_names(c) == [(views.EVERYONE, "#6f7682"),
                                         ("Not high", figures.CATEGORICAL[1])]
    assert not figures.uses_colormap(c)
    assert figures.preview_png(c, STYLE, dpi=40)[:4] == b"\x89PNG"
    svg = _svg(c, tmp_path, figures.FigureStyle(font="DejaVu Sans", legend_counts=True,
                                                labels={"Not high": "Mine"}))
    for text in ("Cells", "Poly", "Round", "Q2", "FSC marker", views.EVERYONE,
                 f"Mine ({chosen:,d})"):
        assert f">{text}<" in svg
    # four panels of density and their dots are bitmaps; the last is a histogram
    assert svg.count("<image") == 8

    # nothing to draw a step on: its panel says so and the figure is still made
    lost = views.backgate(tree, "Cells/Poly/Round/Q2/High/Nothing", masks, display,
                          scales.get, lambda d: ranges.get(d, (0.0, 8.0)), str)
    assert lost.panels[-1].counts is None and lost.panels[-1].x is None
    assert "a combination" in _svg(lost, tmp_path)


def test_backgate_draws_a_capped_number_of_events():
    rng = np.random.default_rng(6)
    events = {"X": rng.uniform(0, 1, 20000), "Y": rng.uniform(0, 1, 20000)}
    tree = GateNode("")
    tree.add(GateNode("Most", RectangleGate(["X", "Y"], [0.1, 0.1], [None, None])))
    masks = tree.apply(lambda d: events[d], 20000)
    c = views.backgate(tree, "Most", masks, lambda d: events[d], lambda d: None,
                       lambda d: (0.0, 1.0), str)
    assert c.events == masks["Most"].sum() > views.MAX_HIGHLIGHTED == c.drawn
    assert c.panels[0].points.shape == (views.MAX_HIGHLIGHTED, 2)
    assert (c.panels[0].points >= 0.1).all()


# --- the window -------------------------------------------------------------------

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
    for name, n in (("3Ctrl1", 240), ("3Ctrl2", 400), ("3Case1", 900), ("3Case2", 1500)):
        pop = rng.integers(0, 3, n)
        fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (n, 3)))
        ev = np.column_stack([rng.normal(90000, 9000, n), fluo]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), ev, ["FSC-A", "FITC-A", "PE-A", "APC-A"],
                  markers=["", "CD8", "CD4", ""],
                  keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
    return str(tmp_path)


@pytest.fixture
def gates(app, folder):
    from milc.ui import worker
    from milc.ui.shell import MainShell
    w = MainShell()
    w.open_folder(folder)
    w.show()
    w.tabs.setCurrentWidget(w.gates)
    for _ in range(3):
        worker.wait_all()
        app.processEvents()
    yield w.gates
    w.close()


def test_violins_and_boxes_in_gates(gates):
    from milc.ui.gates_tab import HISTOGRAM
    g = gates
    g.set_axes("Comp-PE-A", HISTOGRAM)
    assert not g.overlay.isEnabled() and not g.offset.isEnabled()   # one sample: no choice
    names = ["3Case1", "3Case2", "3Ctrl1", "3Ctrl2"]
    g.select_samples(names)
    assert g.overlay.isEnabled() and g.overlay.text() == "Overlaid"
    overlay = g.content()
    assert isinstance(overlay, figures.OverlayContent)

    g.set_overlay("violin")
    c = g.content()
    assert isinstance(c, views.DistributionContent) and c.kind == views.VIOLIN
    assert g.overlay.text() == "Violins" and not g.offset.isChecked()
    assert [s.label for s in c.series] == names == c.groups()
    assert c.xlabel == "CD4" and c.xrange == (0.0, 1.0) and c.xticks[0][1] == "0"
    # a violin is the overlaid histogram of the same sample, and its numbers
    # are those of the events on the axis
    for violin, curve in zip(c.series, overlay.series):
        assert np.array_equal(violin.density, curve.y) and np.array_equal(violin.x, curve.x)
    sample = g.session.load(g.entries[[e.name for e in g.entries].index("3Case2")])
    shown = g._display("Comp-PE-A", sample)
    trt2 = c.series[1].summary
    assert trt2.n == 1500 and trt2.median == pytest.approx(np.median(shown))
    assert (trt2.q1, trt2.q3) == pytest.approx(tuple(np.percentile(shown, [25, 75])))
    assert "median" in g.caption.text() and "4 samples" in g.caption.text()

    g.colour_by.setCurrentIndex(g.colour_by.findData("Group"))
    g.set_overlay("box")
    c = g.content()
    assert c.kind == views.BOX and [s.label for s in c.series] == names
    assert c.groups() == ["Case", "Ctrl"] and "1.5 IQR" in g.caption.text()
    assert figures.preview_png(c, STYLE, dpi=40)[:4] == b"\x89PNG"

    # one for each group: its samples' events together
    g.set_overlay("box", pooled=True)
    c = g.content()
    assert [(s.label, s.group, s.summary.n) for s in c.series] == [("Case", "Case", 2400),
                                                                   ("Ctrl", "Ctrl", 640)]
    both = np.concatenate([g._display("Comp-PE-A", g.session.load(e))
                           for e in g.entries if e.name.startswith("3Case")])
    assert c.series[0].summary.median == pytest.approx(np.median(both))

    # within a population, and with a range gate still on the plot
    g.add_rectangle(0.5, 0, 0.9, 0, "CD4+")
    assert isinstance(g.content(), views.DistributionContent) and len(g._overlays) == 1
    g.select("CD4+")
    inside = g.content().series[0].summary
    assert 0 < inside.n < 2400 and inside.low >= 0.5 - 1e-9
    g.select("")

    # the choice the Staggered box used to be is still one, and excludes the others
    g.offset.setChecked(True)
    assert g.overlay.mode() == "staggered" and g.content().offset
    g.set_overlay("violin", pooled=True)
    assert not g.offset.isChecked() and len(g.content().series) == 2
    # on two parameters, or back to one sample, it is not asked
    g.set_axes("Comp-PE-A", "Comp-FITC-A")
    assert isinstance(g.content(), figures.OverlayContent) and not g.overlay.isEnabled()
    g.set_axes("Comp-PE-A", HISTOGRAM)
    g.select_samples([g.entry.name])
    assert isinstance(g.content(), figures.OverlayContent) and not g.overlay.isEnabled()

    from milc.ui.figure_dialog import FigureDialog
    g.select_samples(names)
    d = FigureDialog(g.content(), g, "violins")
    assert d.categories.rowCount() == 2 and not d.colormap.isVisibleTo(d)
    d.close()


def test_matrix_from_gates(gates, tmp_path):
    g = gates
    assert g.matrix_button.text() == "Matrix…"
    # the stained ones: APC has no marker in these files
    assert g.matrix_default() == ["Comp-FITC-A", "Comp-PE-A"]
    g.set_axes("FSC-A", "Comp-PE-A")
    g.add_rectangle(80000, 0.0, 100000, 1.0, "Cells")
    g.select("Cells")
    n = g.sample.count("Cells")
    c = g.matrix_content()
    assert c.labels == ["CD8", "CD4"] and c.counts.shape == (2, 2, 64, 64)
    assert (c.counts.sum(axis=(2, 3)) == n).all() and 0 < n < len(g.sample)
    assert np.array_equal(c.counts[0, 1], c.counts[1, 0].T)

    d = g.matrix_dialog()
    assert "Cells" in d.windowTitle() and g.entry.name in d.windowTitle()
    assert d.ticked() == ["Comp-FITC-A", "Comp-PE-A"] and d.list.count() == 4
    assert d.content.labels == ["CD8", "CD4"] and not d.preview.pixmap().isNull()
    d.set_parameters(["FSC-A", "Comp-PE-A", "Comp-APC-A"])
    assert d.content.labels == ["FSC-A", "CD4", "APC-A"]
    assert d.content.ticks[0][1][1].endswith("K") and d.content.ticks[1][0] == (pytest.approx(
        float(g.session.scale_of("PE-A")(0.0))), "0")
    assert (d.content.counts.sum(axis=(2, 3)) == n).all()
    figures.save(d.content, STYLE, str(tmp_path / "matrix.svg"))
    d.set_parameters(["Comp-PE-A"])
    assert d.content is None and not d.figure_button.isEnabled()
    assert "two parameters" in d.note.text() and d.make_figure() is None
    d.close()

    from milc.ui.view_dialogs import MatrixDialog
    seen = []

    def build(dims):
        seen.append(dims)
        return g.matrix_content()

    many = [(f"P{i}", f"P{i}") for i in range(15)]
    d = MatrixDialog("Every pair", many, [p for _, p in many], build, g)
    assert len(seen[-1]) == views.MAX_PARAMETERS == 12
    assert "first 12 of the 15" in d.note.text()
    d.close()


def test_backgate_from_gates(gates, tmp_path):
    from milc.ui.gates_tab import HISTOGRAM
    g = gates
    assert g.backgate_button.text() == "Backgate…"
    said = []
    g.sigStatus.connect(said.append)
    assert g.backgate_dialog() is None and "gated population" in said[-1]

    g.set_axes("FSC-A", "Comp-PE-A")
    g.add_rectangle(80000, 0.0, 100000, 1.0, "Cells")
    g.select("Cells")
    g.set_axes("Comp-PE-A", "Comp-FITC-A")
    g.add_polygon(np.array([[0.45, 0.05], [0.95, 0.05], [0.95, 0.5], [0.45, 0.5]]), "CD4")
    g.add_ellipse(0.35, 0.6, 0.15, 0.15, "CD8")
    g.add_quadrants(0.5, 0.5)
    g.select("Cells/CD4")
    g.set_axes("Comp-APC-A", HISTOGRAM)
    g.add_rectangle(0.0, 0, 0.45, 0, "Dim")
    g.select("Cells/CD4/Dim")
    c = g.backgate_content()
    chosen = g.sample.count("Cells/CD4/Dim")
    assert (c.population, c.events, c.drawn) == ("Dim", chosen, chosen) and chosen > 0
    assert [(p.title, p.xlabel, p.ylabel) for p in c.panels] == [
        ("Cells", "FSC-A", "CD4"), ("CD4", "CD4", "CD8"), ("Dim", "APC-A", "Events")]
    assert c.panels[0].counts.sum() == len(g.sample)
    assert c.panels[1].counts.sum() == g.sample.count("Cells")
    assert c.panels[0].xticks is None and c.panels[1].xticks[0][1] == "0"
    # the polygon where it was drawn, the range as its two bounds
    assert np.allclose(c.panels[1].outline[0][:4],
                       [[0.45, 0.05], [0.95, 0.05], [0.95, 0.5], [0.45, 0.5]])
    assert [line[0, 0] for line in c.panels[2].outline] == pytest.approx([0.0, 0.45])

    d = g.backgate_dialog()
    assert "Dim" in d.windowTitle() and not d.preview.pixmap().isNull()
    assert "events of Dim" in d.note.text() and d.content.population == "Dim"
    path = str(tmp_path / "backgate.svg")
    figures.save(d.content, STYLE, path)
    with open(path, encoding="utf-8") as f:
        svg = f.read()
    assert ">Dim<" in svg and ">CD4<" in svg and ">APC-A<" in svg
    d.close()

    # an ellipse and a corner of a quadrant are steps like any other
    for path, lines in (("Cells/CD8", 1), ("Cells/CD4+ CD8-", 2)):
        g.select(path)
        last = g.backgate_content().panels[-1]
        assert len(last.outline) == lines and last.counts.sum() == g.sample.count("Cells")
    from milc.ui.figure_dialog import FigureDialog
    f = FigureDialog(g.backgate_content(), g, "backgate")
    assert f.categories.rowCount() == 2 and not f.colormap.isVisibleTo(f)
    f.close()
