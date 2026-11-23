import os

import numpy as np
import pytest

from milc import sampling
from milc.compensation import Compensation, from_keywords, strip_prefix
from milc.dataset import Dataset, concatenate
from milc.fcs import read_fcs, read_header, write_fcs
from milc.gating import GateNode, PolygonGate, RectangleGate, points_in_polygon
from milc.sample import Sample, annotate
from milc.transforms import Arcsinh, Biex, Linear


def _write(tmp_path, name, events, channels=("FSC-A", "FITC-A", "PE-A"), **kw):
    path = str(tmp_path / name)
    write_fcs(path, events, list(channels), **kw)
    return path


# --- fcs -----------------------------------------------------------------

def test_fcs_round_trip(tmp_path):
    rng = np.random.default_rng(1)
    ev = rng.normal(1000, 300, (500, 3)).astype(np.float32)
    path = _write(tmp_path, "a.fcs", ev, markers=["", "CD8", "CD4"],
                  keywords={"SPILL": "2,FITC-A,PE-A,1,0.1,0.02,1", "NOTE": "a|b"})
    f = read_fcs(path)
    assert np.array_equal(f.events, ev)
    assert f.channels == ["FSC-A", "FITC-A", "PE-A"]
    assert f.markers == ["", "CD8", "CD4"]
    assert f.keywords["NOTE"] == "a|b"  # the delimiter survives in a value
    assert f.label("FITC-A") == "CD8 (FITC-A)" and f.label("FSC-A") == "FSC-A"
    assert read_header(path)["$TOT"] == "500"


def test_fcs_integer_data(tmp_path):
    vals = np.array([[0, 1023], [512, 5]], dtype="<u2")
    text = ("|$BEGINANALYSIS|0|$ENDANALYSIS|0|$BEGINSTEXT|0|$ENDSTEXT|0"
            "|$BYTEORD|1,2|$DATATYPE|I|$MODE|L|$NEXTDATA|0|$PAR|2|$TOT|2"
            "|$P1B|16|$P1E|0,0|$P1N|A|$P1R|1024"
            "|$P2B|16|$P2E|4,1|$P2N|B|$P2R|1024|").encode()
    t0, d0 = 58, 58 + len(text)
    head = b"FCS3.0    " + b"%8d%8d%8d%8d%8d%8d" % (
        t0, d0 - 1, d0, d0 + vals.nbytes - 1, 0, 0)
    p = tmp_path / "i.fcs"
    p.write_bytes(head + text + vals.tobytes())
    f = read_fcs(str(p))
    assert f.events[:, 0].tolist() == [0, 512]
    # four decades over 1024 channels
    assert f.events[0, 1] == pytest.approx(10 ** (4 * 1023 / 1024), rel=1e-5)


# --- compensation --------------------------------------------------------

def test_compensation_recovers_the_true_signal():
    S = np.array([[1, 0.2], [0.05, 1]])
    true = np.array([[1000.0, 0], [0, 500], [300, 300]])
    observed = np.column_stack([np.full(3, 7.0), true @ S])
    c = Compensation(["FITC-A", "PE-A"], S)
    out = c.apply(observed.astype(np.float32), ["FSC-A", "FITC-A", "PE-A"])
    assert np.allclose(out[:, 1:], true, atol=1e-3)
    assert np.all(out[:, 0] == 7)  # scatter untouched


def test_spill_keyword_and_prefixes():
    c = from_keywords({"SPILL": "2,FITC-A,PE-A,1,0.1,0.02,1"})
    assert c.channels == ["FITC-A", "PE-A"] and c.matrix[0, 1] == 0.1
    assert from_keywords({}) is None
    assert strip_prefix("Comp-PE-A") == ("PE-A", True)
    assert strip_prefix("FJComp-PE-A") == ("PE-A", True)
    assert strip_prefix("FSC-A") == ("FSC-A", False)


# --- transforms ----------------------------------------------------------

@pytest.mark.parametrize("t", [Linear(0, 262144), Arcsinh(150), Biex()])
def test_transforms_invert(t):
    x = np.array([-500.0, -10, 0, 10, 1000, 50000, 262144])
    assert np.allclose(t.inverse(t(x)), x, rtol=1e-6, atol=1e-3)
    assert np.all(np.diff(t(x)) > 0)


def test_biex_shape():
    b = Biex()
    # the routine's table stops a hair short of the top channel
    assert b(262144) == pytest.approx(1.0, abs=1e-3)
    assert 0 < b(0) < 0.5
    # symmetric about zero, as a biexponential is
    assert b(0) - b(-200) == pytest.approx(b(200) - b(0), rel=1e-6)


# --- gates ---------------------------------------------------------------

def test_polygon():
    square = np.array([[0, 0], [10, 0], [10, 10], [0, 10.0]])
    x = np.array([5, 15, -1, 9.99, 5])
    y = np.array([5, 5, 5, 9.99, 11])
    assert points_in_polygon(x, y, square).tolist() == [True, False, False, True, False]
    # concave: the notch is outside
    c = np.array([[0, 0], [10, 0], [10, 10], [5, 2], [0, 10.0]])
    assert points_in_polygon(np.array([5.0, 1]), np.array([5.0, 1]), c).tolist() == [False, True]


def test_tree_paths_and_masks():
    cols = {"a": np.array([1.0, 5, 9, 5]), "b": np.array([1.0, 5, 5, 50])}
    root = GateNode("s")
    cells = root.add(GateNode("Cells", RectangleGate(["a"], [2], [None])))
    one = cells.add(GateNode("Single", RectangleGate(["b"], [None], [10])))
    two = one.add(GateNode("Single", RectangleGate(["a", "b"], [0, 0], [6, 6])))
    cells.add(GateNode("Not", RectangleGate(["b"], [None], [10], inside=False)))
    assert two.path == "Cells/Single/Single"
    m = root.apply(cols.__getitem__, 4)
    assert m[""].sum() == 4
    assert m["Cells"].tolist() == [False, True, True, True]
    assert m["Cells/Single"].tolist() == [False, True, True, False]
    assert m["Cells/Single/Single"].tolist() == [False, True, False, False]
    assert m["Cells/Not"].tolist() == [False, False, False, True]
    assert [n.path for n in root.leaves()] == ["Cells/Single/Single", "Cells/Not"]
    assert root.find("Cells/Single") is one


# --- sampling ------------------------------------------------------------

AVAIL = {"m1": 200, "m2": 300, "h1": 1000, "h2": 5000}
GROUPS = {"m1": "Ctrl", "m2": "Ctrl", "h1": "Case", "h2": "Case"}


def test_equal_per_sample_is_set_by_the_smallest():
    p = sampling.plan(AVAIL, GROUPS)
    assert set(p.take.values()) == {200} and p.limited_by == "m1"
    assert p.group_totals() == {"Ctrl": 400, "Case": 400}
    p = sampling.plan(AVAIL, GROUPS, n=50)
    assert set(p.take.values()) == {50} and p.limited_by is None
    p = sampling.plan(AVAIL, GROUPS, n=900)
    assert set(p.take.values()) == {200} and p.notes


def test_minimum_excludes_instead_of_dragging_down():
    p = sampling.plan(AVAIL, GROUPS, min_events=250)
    assert p.take == {"m1": 0, "m2": 300, "h1": 300, "h2": 300}
    assert "m1" in p.excluded and p.limited_by == "m2"
    # a measured zero is a sample with nothing, not a missing setting
    p = sampling.plan({"a": 0, "b": 10})
    assert p.take == {"a": 0, "b": 10} and "a" in p.excluded


def test_equal_per_group():
    p = sampling.plan(AVAIL, GROUPS, sampling.EQUAL_PER_GROUP)
    assert p.group_totals() == {"Ctrl": 500, "Case": 500}
    assert p.take == {"m1": 200, "m2": 300, "h1": 250, "h2": 250}
    assert p.limited_by == "Ctrl"
    p = sampling.plan(AVAIL, GROUPS, sampling.EQUAL_PER_GROUP, n=450)
    assert p.take["m1"] == 200 and p.take["m2"] == 250
    assert p.group_totals() == {"Ctrl": 450, "Case": 450} and p.notes
    p = sampling.plan({"a": 3, "b": 100, "c": 100}, {}, sampling.EQUAL_PER_GROUP, n=10)
    assert sum(p.take.values()) == 10 and p.take["a"] == 3


def test_ceiling():
    p = sampling.plan(AVAIL, GROUPS, sampling.CEILING, n=400)
    assert p.take == {"m1": 200, "m2": 300, "h1": 400, "h2": 400} and p.notes
    with pytest.raises(ValueError):
        sampling.plan(AVAIL, GROUPS, sampling.CEILING)


def test_draw_is_reproducible_and_independent_of_the_batch():
    mask = np.zeros(1000, dtype=bool)
    mask[::2] = True
    a = sampling.draw(mask, 50, seed=3, key="s1")
    assert np.array_equal(a, sampling.draw(mask, 50, seed=3, key="s1"))
    assert not np.array_equal(a, sampling.draw(mask, 50, seed=4, key="s1"))
    assert not np.array_equal(a, sampling.draw(mask, 50, seed=3, key="s2"))
    assert len(set(a)) == 50 and mask[a].all() and np.all(np.diff(a) > 0)
    assert sampling.draw(mask, 9999, 0, "x").size == 500


# --- end to end ----------------------------------------------------------

def test_concatenate_and_files(tmp_path):
    rng = np.random.default_rng(0)
    samples = []
    for name, n in (("3Ctrl1", 300), ("3Case1", 2000)):
        ev = np.abs(rng.normal(1000, 400, (n, 3))).astype(np.float32)
        path = _write(tmp_path, name + ".fcs", ev,
                      keywords={"SPILL": "2,FITC-A,PE-A,1,0.1,0.02,1"})
        s = Sample.open(path)
        s.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [800], [None])))
        s.meta = annotate(s.name, r"(?P<day>\d+)(?P<group>[A-Za-z]+)(?P<animal>\d+)")
        samples.append(s)
    assert samples[1].meta == {"day": "3", "group": "Case", "animal": "1"}
    avail = {s.name: s.count("Cells") for s in samples}
    p = sampling.plan(avail, {s.name: s.meta["group"] for s in samples})
    ds = concatenate(samples, p, "Cells", seed=7)
    k = avail["3Ctrl1"]
    assert ds.counts() == {"3Ctrl1": k, "3Case1": k}
    assert ds.compensated == [False, True, True]
    assert (ds.column("FSC-A") >= 800).all()
    # every row is the compensated event it says it is
    trt = ds.sample_index == 1
    assert np.array_equal(ds.X[trt], samples[1].compensated[ds.event_index[trt]])
    assert set(ds.obs("group")[trt]) == {"Case"}
    assert np.array_equal(ds.X, concatenate(samples, p, "Cells", seed=7).X)

    ds.save(str(tmp_path / "c.npz"))
    back = Dataset.load(str(tmp_path / "c.npz"))
    assert np.array_equal(back.X, ds.X) and back.meta == ds.meta
    assert back.provenance["seed"] == 7 and back.provenance["limited_by"] == "3Ctrl1"

    ds.to_fcs(str(tmp_path / "c.fcs"))
    f = read_fcs(str(tmp_path / "c.fcs"))
    assert f.channels == ["FSC-A", "Comp-FITC-A", "Comp-PE-A", "SampleID"]
    assert np.array_equal(f.events[:, :3], ds.X)
    assert f.events[:, 3].tolist() == (ds.sample_index + 1).tolist()


# --- analysis ------------------------------------------------------------

def _blobs(n=300, seed=0):
    """Three populations in four markers, well apart."""
    rng = np.random.default_rng(seed)
    centres = np.array([[200, 200, 600, 200], [600, 200, 200, 200], [200, 600, 200, 600.0]])
    truth = np.repeat(np.arange(3), n)
    return centres[truth] + rng.normal(0, 25, (3 * n, 4)), truth


def test_flowsom_finds_the_populations():
    from milc import flowsom as fs
    X, truth = _blobs()
    r = fs.flowsom(X, n_clusters=3, xdim=6, ydim=6, seed=1)
    assert fs.adjusted_rand(r.labels, truth) > 0.98
    assert r.n_clusters == 3 and r.codes.shape == (36, 4)
    sizes = np.bincount(r.labels)[1:]
    assert np.all(np.diff(sizes) <= 0)  # 1 is the largest
    assert r.mst.shape == (35, 2)  # a tree over the nodes
    assert np.array_equal(r.labels, fs.flowsom(X, 3, 6, 6, seed=1).labels)
    assert np.array_equal(r.map(X), r.labels)
    assert r.node_sizes().sum() == len(X)


def test_embeddings_keep_neighbours():
    from milc import embed
    X, truth = _blobs(150)
    for Y in (embed.umap(X, seed=2), embed.tsne(X, seed=2)):
        assert Y.shape == (450, 2) and np.isfinite(Y).all()
        assert embed.neighbour_preservation(X, Y) > 0.25  # chance is 0.03
        # the three populations are apart on the map too
        c = np.array([Y[truth == i].mean(axis=0) for i in range(3)])
        nearest = np.argmin(((Y[:, None, :] - c[None]) ** 2).sum(axis=2), axis=1)
        assert (nearest == truth).mean() > 0.97
    assert embed.neighbour_preservation(X, X) == 1.0
    assert np.array_equal(embed.umap(X, seed=2), embed.umap(X, seed=2))


def test_cluster_tables(tmp_path):
    from milc import cluster_stats as cs
    labels = np.array([1, 1, 2, 2, 2, 1])
    sidx = np.array([0, 0, 0, 1, 1, 1])
    X = np.array([[1.0, 10], [3, 10], [5, 20], [5, 40], [7, 30], [100, 10]])
    e = cs.expression(X, labels, ["a", "b"])
    assert e.rows == ["1", "2"] and e.values.tolist() == [[3, 10], [5, 30]]
    assert cs.expression(X, labels, ["a", "b"], "mean").values[0, 0] == pytest.approx(104 / 3)
    f = cs.frequencies(labels, sidx, ["s1", "s2", "empty"])
    assert np.allclose(f.values[:, :2], [[200 / 3, 100 / 3], [100 / 3, 200 / 3]])
    assert np.isnan(f.values[:, 2]).all()  # no events is not 0%
    mean, sd, n = cs.by_group(f, ["A", "A", "B"])
    assert mean.columns == ["A", "B"] and mean.values[0, 0] == pytest.approx(50)
    assert n.values[0].tolist() == [2, 0] and np.isnan(mean.values[0, 1])
    assert sd.values[0, 0] == pytest.approx(np.std([200 / 3, 100 / 3], ddof=1))
    f.to_csv(str(tmp_path / "f.csv"), "cluster")
    assert (tmp_path / "f.csv").read_text(encoding="utf-8-sig").splitlines()[0] == "cluster,s1,s2,empty"


def test_analyze_command(tmp_path):
    from milc.cli import main
    X, _ = _blobs(120)
    ev = np.column_stack([X, np.repeat([1, 2, 3, 4], 90)]).astype(np.float32)
    path = _write(tmp_path, "c.fcs", ev,
                  channels=("FJComp-FITC-A", "FJComp-PE-A", "APC-A", "BV421-A", "SampleID"),
                  markers=["CD8 FITC-A", "CD4 PE-A", "Foxp3", "", ""])
    out = str(tmp_path / "run")
    assert main(["analyze", path, "--channels", "cd8,CD4,Foxp3,BV421-A", "--umap",
                 "--flowsom", "--clusters", "3", "--grid", "5", "--out", out]) == 0
    ds = Dataset.load(out + ".npz")
    assert ds.channels[:2] == ["FITC-A", "PE-A"] and ds.compensated[:2] == [True, True]
    assert set(ds.derived) == {"UMAP1", "UMAP2", "FlowSOM", "FlowSOM_node"}
    assert ds.runs["FlowSOM"]["channels"] == ["FITC-A", "PE-A", "APC-A", "BV421-A"]
    assert len(ds.samples) == 4 and set(np.unique(ds.derived["FlowSOM"])) == {1, 2, 3}
    lines = (tmp_path / "run_events.csv").read_text(encoding="utf-8-sig").splitlines()
    assert lines[0] == "sample,event,UMAP1,UMAP2,FlowSOM,FlowSOM_node" and len(lines) == 361
    assert (tmp_path / "run_mfi.csv").read_text(encoding="utf-8-sig").startswith("cluster,CD8 FITC-A,CD4 PE-A,Foxp3,BV421-A")
    with pytest.raises(KeyError):
        ds.resolve("CD")


def test_tree_and_selection():
    from milc import cluster_stats as cs
    from milc import flowsom as fs
    X, _ = _blobs()
    r = fs.flowsom(X, 3, 5, 5, seed=0)
    sizes = r.node_sizes()
    pos, edges = fs.tree(r.codes, sizes)
    used = sizes > 0
    assert np.isfinite(pos[used]).all() and np.isnan(pos[~used]).all()
    assert len(edges) == used.sum() - 1 and used[edges].all()  # a tree over what is populated
    assert np.array_equal(pos[used], fs.tree(r.codes, sizes)[0][used])
    # nodes joined by an edge sit closer than the picture is wide
    span = np.ptp(pos[used], axis=0).max()
    assert np.hypot(*(pos[edges[:, 0]] - pos[edges[:, 1]]).T).max() < span / 2

    mask = np.array([True, True, False, False, True, False])
    info = cs.selection(mask, np.array([0, 0, 0, 1, 1, 2]), ["a", "b", "c", "d"],
                        ["7", "7", "14", "14"])
    assert info["n"] == 3 and info["share"] == 0.5
    assert list(info["by_group"]) == ["7", "14"]  # days in order, not as text
    assert info["by_group"]["7"][0] == pytest.approx((200 / 3 + 50) / 2)
    assert info["by_group"]["14"] == (0.0, pytest.approx(float("nan"), nan_ok=True), 1)
    assert cs.ordered(["Ctrl", "Case", "Ctrl"]) == ["Case", "Ctrl"]
    t = cs.Table(["1", "2"], ["a", "b"], np.array([[0.0, 5], [10, 5]]))
    assert cs.scaled(t).tolist() == [[0, 0], [1, 0]]
    assert cs.scaled(t, per_column=False).tolist() == [[0, 0.5], [1, 0.5]]


# --- figures ---------------------------------------------------------------

def _map_content(n=400):
    from milc import figures
    rng = np.random.default_rng(0)
    group = np.where(np.arange(n) % 2 == 0, "Case", "Ctrl")
    return figures.MapContent(rng.normal(size=n), rng.normal(size=n), "UMAP1", "UMAP2",
                              color=figures.categories_of(group, "Group"))


def test_figure_files(tmp_path):
    from PIL import Image

    from milc import figures
    c = _map_content()
    style = figures.FigureStyle(width_mm=85, height_mm=65, dpi=300, font="DejaVu Sans",
                                colors={"Case": "#c0392b"}, labels={"Case": "Case-1"},
                                title="T cells")
    for ext in ("svg", "pdf", "png", "tiff"):
        figures.save(c, style, str(tmp_path / f"f.{ext}"))
    # the size asked for, at the resolution asked for
    for ext in ("png", "tiff"):
        im = Image.open(tmp_path / f"f.{ext}")
        # to the pixel, give or take the one matplotlib rounds away
        assert abs(im.size[0] - 85 / 25.4 * 300) <= 1 and abs(im.size[1] - 65 / 25.4 * 300) <= 1
        assert round(im.info["dpi"][0]) == 300
    svg = (tmp_path / "f.svg").read_text(encoding="utf-8")
    # text is text — editable — and says what was asked
    assert "<text" in svg and "Case-1" in svg and "T cells" in svg and "UMAP1" in svg
    assert "#c0392b" in svg.lower() and "Ctrl" in svg
    assert 'width="240.9' in svg  # 85 mm in points
    assert (tmp_path / "f.pdf").read_bytes().startswith(b"%PDF")
    with pytest.raises(ValueError):
        figures.save(c, style, str(tmp_path / "f.bmp"))

    # as a bitmap inside the vector file, the points are no longer paths
    style.rasterize = True
    figures.save(c, style, str(tmp_path / "r.svg"))
    assert "<image" in (tmp_path / "r.svg").read_text(encoding="utf-8") and "<image" not in svg


def test_every_kind_of_figure_draws(tmp_path):
    from milc import cluster_stats as cs
    from milc import figures
    from milc import flowsom as fs
    from milc.transforms import Biex, decade_ticks
    X, _ = _blobs(100)
    r = fs.flowsom(X, 3, 4, 4, seed=0)
    sidx = np.arange(len(X)) % 4
    names = ["CD8 FITC", "CD4 PE", "Foxp3", "IFNg"]
    pos, edges = fs.tree(r.codes, r.node_sizes())
    tree = {"pos": pos, "edges": edges, "sizes": r.node_sizes(),
            "cluster": r.node_cluster, "medians": r.codes, "markers": names}
    ticks = decade_ticks(Biex())
    assert ticks[0][1] == "0" and [t for _, t in ticks][-2:] == ["10⁴", "10⁵"]
    assert all(b[0] - a[0] > 0.04 for a, b in zip(ticks, ticks[1:]))  # none on top of another
    split = figures.categories_of(sidx.astype(str), "Sample")
    contents = {
        "values": figures.MapContent(X[:, 0], X[:, 1], "a", "b", xticks=ticks,
                                     color=figures.Values(X[:, 2], "Foxp3", ticks=ticks)),
        "split": figures.MapContent(X[:, 0], X[:, 1], "a", "b", split=split,
                                    color=figures.categories_of(r.labels, "Cluster")),
        "heat": figures.HeatmapContent(cs.expression(X, r.labels, names)),
        "freq": figures.FrequencyContent(
            cs.frequencies(r.labels, sidx, ["a", "b", "c", "d"]), ["7", "7", "14", "14"]),
        "tree": figures.TreeContent(tree),
        "tree_marker": figures.TreeContent(tree, "Foxp3"),
    }
    for legend in ("right", "bottom", "inside", "none"):
        for frame in ("open", "box", "none"):
            style = figures.FigureStyle(font="DejaVu Sans", legend=legend, frame=frame,
                                        legend_counts=True)
            for name, c in contents.items():
                assert figures.preview_png(c, style, dpi=40)[:4] == b"\x89PNG", name
    assert [n for n, _ in figures.category_names(contents["freq"])] == ["7", "14"]
    assert figures.uses_colormap(contents["heat"]) and not figures.uses_colormap(contents["tree"])
    figures.save(contents["split"], figures.FigureStyle(font="DejaVu Sans"),
                 str(tmp_path / "s.svg"))
    assert (tmp_path / "s.svg").read_text(encoding="utf-8").count("Cluster") >= 1

    # a preset carries the look and none of a figure's own words
    s = figures.FigureStyle(title="x", colormap="magma")
    assert "title" not in figures.PRESET_FIELDS and "colormap" in figures.PRESET_FIELDS
    assert figures.FigureStyle.from_dict({**s.to_dict(), "unknown": 1}).colormap == "magma"


WSP = """<?xml version="1.0"?>
<Workspace flowJoVersion="10" xmlns:g="http://g" xmlns:t="http://t" xmlns:d="http://d">
 <Groups><GroupNode name="All Samples"><Subpopulations>{gate}</Subpopulations></GroupNode></Groups>
 <SampleList>{samples}</SampleList>
</Workspace>"""
GATE = """<Population name="Cells" count="{n}"><Gate><g:RectangleGate eventsInside="1">
 <g:dimension g:min="{lo}"><d:fcs-dimension d:name="FSC-A"/></g:dimension>
 </g:RectangleGate></Gate></Population>"""
SAMPLE = """<Sample><DataSet uri="file:/C:/x/{name}.fcs" sampleID="1"/>
 <t:spilloverMatrix name="edited"><d:parameters><d:parameter d:name="FITC-A"/><d:parameter d:name="PE-A"/></d:parameters>
 <t:spillover d:parameter="FITC-A"><t:coefficient d:parameter="FITC-A" t:value="1"/><t:coefficient d:parameter="PE-A" t:value="{spill}"/></t:spillover>
 <t:spillover d:parameter="PE-A"><t:coefficient d:parameter="FITC-A" t:value="0"/><t:coefficient d:parameter="PE-A" t:value="1"/></t:spillover>
 </t:spilloverMatrix>
 <SampleNode name="{name}.fcs" count="100"><Subpopulations>{gate}</Subpopulations></SampleNode></Sample>"""


def test_shared_gates_reach_unlisted_files(tmp_path):
    from milc.session import Session
    rng = np.random.default_rng(0)
    for name in ("3Ctrl1", "3Case1", "3Case2"):
        ev = np.column_stack([np.linspace(0, 999, 100), rng.random((100, 2)) * 1000])
        _write(tmp_path, name + ".fcs", ev.astype(np.float32))

    def workspace(spills, lows):
        gate = GATE.format(n=50, lo=lows[0])
        (tmp_path / "a.wsp").write_text(WSP.format(gate=gate, samples="".join(
            SAMPLE.format(name=n, spill=sp, gate=GATE.format(n=50, lo=lo))
            for n, sp, lo in zip(("3Ctrl1", "3Case1"), spills, lows))))
        s = Session()
        s.open_folder(str(tmp_path))
        return s

    s = workspace((0.2, 0.2), (500, 500))
    extra = s.entries[1]
    assert extra.name == "3Case2" and extra.borrowed and s.can_extend()
    assert s.population == "Cells" and s.count(extra) == 50
    assert s.load(extra).compensation.matrix[0, 1] == 0.2  # the edited matrix, not its own
    assert s.load(extra).gates.find("Cells").reference_count is None
    s.set_extend_gates(False)
    assert not extra.borrowed and not s.available(extra)

    # one sample compensated or gated its own way: nothing to hand on
    assert not workspace((0.2, 0.3), (500, 500)).can_extend()
    s = workspace((0.2, 0.2), (500, 600))
    assert not s.can_extend() and not s.entries[1].borrowed


def test_selftest_never_raises_when_logging(tmp_path, monkeypatch):
    """A windowed build on Windows cannot show an error: it has to end up
    in the file and in the exit code."""
    from milc import cli
    log = tmp_path / "s.log"
    assert cli.main(["selftest", "--no-window", "--log", str(log)]) == 0
    assert "self-test passed" in log.read_text(encoding="utf-8")
    monkeypatch.setattr("milc.session.Session.run_umap",
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("no UMAP today")))
    assert cli.main(["selftest", "--no-window", "--log", str(log)]) == 1
    assert "no UMAP today" in log.read_text(encoding="utf-8") and "Traceback" in log.read_text(encoding="utf-8")


# --- statistics --------------------------------------------------------------

def test_statistic_definitions():
    from milc import stats
    x = np.array([1.0, 2, 3, 4, 100])
    assert stats.compute("Median", x) == 3 and stats.compute("Mean", x) == 22
    assert stats.compute("SD", x) == pytest.approx(np.std(x, ddof=1))
    assert stats.compute("CV", x) == pytest.approx(100 * np.std(x, ddof=1) / 22)
    assert stats.compute("MAD", x) == 1
    lo, hi = np.percentile(x, [15.87, 84.13])
    assert stats.compute("Robust SD", x) == pytest.approx((hi - lo) / 2)
    assert stats.compute("Robust CV", x) == pytest.approx(100 * (hi - lo) / 2 / 3)
    # no logarithm at or below zero: those events are left out, not an error
    assert stats.compute("Geometric mean", np.array([-5.0, 0, 10, 1000])) == pytest.approx(100)
    assert np.isnan(stats.compute("Geometric mean", np.array([-1.0, 0])))
    # nothing to compute on is no value, for every statistic
    assert all(np.isnan(stats.compute(name, np.array([]))) for name in stats.STATISTICS)
    assert np.isnan(stats.compute("SD", np.array([5.0]))) and np.isnan(stats.compute("CV", np.zeros(3)))
    # a robust SD of normal data is its SD
    normal = np.random.default_rng(0).normal(50, 7, 200_000)
    assert stats.compute("Robust SD", normal) == pytest.approx(7, rel=0.01)


def test_statistics_table(tmp_path):
    from milc import stats
    from milc.cli import main
    from milc.session import Session
    rng = np.random.default_rng(1)
    for name, n in (("3Ctrl1", 400), ("3Case1", 900)):
        ev = np.column_stack([np.linspace(0, 999, n), np.abs(rng.normal(500, 100, (n, 2)))])
        _write(tmp_path, name + ".fcs", ev.astype(np.float32),
               markers=["", "CD8", "CD4"], keywords={"SPILL": "2,FITC-A,PE-A,1,0.1,0.02,1"})
    s = Session()
    s.open_folder(str(tmp_path))
    assert s.stat_channels() == [("FSC-A", "FSC-A"), ("CD8", "Comp-FITC-A"), ("CD4", "Comp-PE-A")]
    for e in s.usable():
        cells = e.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [500], [None])))
        cells.add(GateNode("Top", RectangleGate(["FSC-A"], [750], [None])))
    t = s.statistics(channels=[("CD8", "Comp-FITC-A")], statistics=["Median", "CV"])
    assert t.columns == ["Sample", "Group", "Population", "Path", "Count", "% of parent",
                         "% of grandparent", "% of total", "CD8 · Median", "CD8 · CV"]
    assert len(t.rows) == 6 and t.column("Population")[:3] == ["All events", "Cells", "Top"]
    trt = {r[3]: r for r in t.rows if r[0] == "3Case1"}
    sample = s.load(s.entries[0])
    fsc, cd8 = sample.column("FSC-A"), sample.column("Comp-FITC-A")
    assert trt["Cells"][4] == (fsc >= 500).sum() and trt["Cells/Top"][4] == (fsc >= 750).sum()
    assert trt["Cells"][5] == pytest.approx(100 * (fsc >= 500).mean())
    assert np.isnan(trt[""][5]) and np.isnan(trt["Cells"][6])       # no parent, no grandparent
    assert trt["Cells/Top"][5] == pytest.approx(100 * (fsc >= 750).sum() / (fsc >= 500).sum())
    assert trt["Cells/Top"][6] == pytest.approx(trt["Cells/Top"][7])  # grandparent is everything
    assert trt["Cells/Top"][8] == pytest.approx(np.median(cd8[fsc >= 750]))   # compensated values
    assert trt["Cells"][1] == "Case"                                  # the sample's own column
    # only what was asked for
    one = s.statistics(populations=["Cells"])
    assert len(one.rows) == 2 and one.columns[-1] == "% of total"

    path = tmp_path / "out.csv"
    t.to_csv(str(path))
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    assert lines[0].startswith("Sample,Group,Population,Path,Count") and len(lines) == 7
    assert lines[1].split(",")[5] == ""                              # no value is an empty cell
    assert path.read_bytes()[:3] == b"\xef\xbb\xbf"                   # the mark Excel needs
    assert t.to_text().splitlines()[0].split("\t")[-1] == "CD8 · CV"

    s.save_project()                                                 # the gates, for the command line
    out = tmp_path / "cli.csv"
    assert main(["stats", str(tmp_path), "--channels", "cd8,FSC-A", "--stats", "median,robust cv",
                 "--out", str(out)]) == 0
    head = out.read_text(encoding="utf-8-sig").splitlines()[0]
    assert head.endswith("CD8 · Median,CD8 · Robust CV,FSC-A · Median,FSC-A · Robust CV")
    assert main(["stats", str(tmp_path), "--channels", "nope"]) == 2
    assert main(["stats", str(tmp_path), "--stats", "mode"]) == 2


# --- density: histograms and contours ---------------------------------------------

def test_histogram_and_contours():
    from milc import density
    from milc.gating import points_in_polygon
    rng = np.random.default_rng(0)
    x, y = rng.normal(0.5, 0.08, 20000), rng.normal(0.4, 0.05, 20000)
    centres, heights = density.histogram(x, 0, 1)
    assert heights.max() == pytest.approx(100) and abs(centres[np.argmax(heights)] - 0.5) < 0.03
    _, counts = density.histogram(x, 0, 1, smooth=0, normalise=density.COUNT)
    assert counts.sum() == 20000
    # events beyond the axis pile up at its end; they are not lost
    _, edge = density.histogram(np.array([5.0, 6.0, 0.5]), 0, 1, bins=10, smooth=0,
                                normalise=density.COUNT)
    assert edge[-1] == 2 and edge.sum() == 3
    # a contour line encloses the share of the events it is named for
    levels = density.contours(x, y, (0, 1), (0, 1))
    assert [share for share, _ in levels] == list(density.LEVELS)
    for share, lines in levels:
        inside = sum(points_in_polygon(x, y, line).sum() for line in lines) / x.size
        assert inside == pytest.approx(share, abs=0.04)
    far = density.outside(x, y, (0, 1), (0, 1))
    assert far.mean() == pytest.approx(0.10, abs=0.03)
    assert density.contours(x[:5], y[:5], (0, 1), (0, 1)) == []   # too few to have a shape


def test_overlay_figures(tmp_path):
    from milc import density, figures
    rng = np.random.default_rng(1)
    hist, lines = [], []
    for k, (name, group) in enumerate((("a", "Ctrl"), ("b", "Case"), ("c", "Case"))):
        x, y = rng.normal(0.4 + 0.1 * k, 0.05, 3000), rng.normal(0.5, 0.06, 3000)
        hist.append(figures.Series(name, group, *density.histogram(x, 0, 1)))
        lines.append(figures.Series(name, group, levels=density.contours(x, y, (0, 1), (0, 1))))
    overlay = figures.OverlayContent(hist, "CD4", "% of the tallest bin", xrange=(0, 1))
    assert overlay.is_histogram and overlay.groups() == ["Ctrl", "Case"]
    assert [n for n, _ in figures.category_names(overlay)] == ["Ctrl", "Case"]
    style = figures.FigureStyle(font="DejaVu Sans", colors={"Case": "#c0392b"})
    figures.save(overlay, style, str(tmp_path / "h.svg"))
    svg = (tmp_path / "h.svg").read_text(encoding="utf-8")
    assert "#c0392b" in svg.lower() and "CD4" in svg and "Ctrl" in svg
    overlay.offset = True
    figures.save(overlay, style, str(tmp_path / "o.svg"))
    assert ">a<" in (tmp_path / "o.svg").read_text(encoding="utf-8").replace("\n", "")
    contour = figures.OverlayContent(lines, "CD4", "CD8", xrange=(0, 1), yrange=(0, 1))
    assert not contour.is_histogram
    assert figures.preview_png(contour, style, dpi=40)[:4] == b"\x89PNG"


# --- ellipse, quadrant and combined gates -----------------------------------------

def test_ellipse_quadrant_and_boolean_gates():
    from milc.gating import (AND, NOT, OR, BooleanGate, EllipseGate, QuadrantGate,
                             copy_tree, tree_from_dict, tree_to_dict)
    from milc.transforms import Biex
    rng = np.random.default_rng(0)
    cols = {"a": rng.uniform(-10, 10, 4000), "b": rng.uniform(-10, 10, 4000)}
    col = cols.__getitem__
    a, b = cols["a"], cols["b"]

    # a circle, an ellipse along an axis, and the same ellipse turned
    assert np.array_equal(EllipseGate(["a", "b"], [0, 0], [5, 5]).contains(col),
                          a ** 2 + b ** 2 <= 25)
    flat = EllipseGate(["a", "b"], [1, -2], [6, 2]).contains(col)
    assert np.array_equal(flat, ((a - 1) / 6) ** 2 + ((b + 2) / 2) ** 2 <= 1)
    upright = EllipseGate(["a", "b"], [1, -2], [6, 2], angle=np.pi / 2).contains(col)
    assert np.array_equal(upright, ((a - 1) / 2) ** 2 + ((b + 2) / 6) ** 2 <= 1)
    assert flat.sum() == pytest.approx(upright.sum(), rel=0.1)   # same area, turned
    assert np.array_equal(EllipseGate(["a", "b"], [0, 0], [5, 5], inside=False).contains(col),
                          a ** 2 + b ** 2 > 25)
    # on a scale, it is the drawn ellipse that counts
    t = Biex()
    raw = {"a": rng.uniform(0, 20000, 2000), "b": rng.uniform(0, 20000, 2000)}
    drawn = EllipseGate(["a", "b"], [0.6, 0.6], [0.1, 0.05], scales=[t, t]).contains(raw.__getitem__)
    assert np.array_equal(drawn, ((t(raw["a"]) - 0.6) / 0.1) ** 2 + ((t(raw["b"]) - 0.6) / 0.05) ** 2 <= 1)

    # the four corners of a crossing: every event in exactly one of them
    corners = {q: QuadrantGate(["a", "b"], 2.0, -1.0, q, "g").contains(col)
               for q in ("Q1", "Q2", "Q3", "Q4")}
    assert np.array_equal(sum(m.astype(int) for m in corners.values()), np.ones(4000))
    assert np.array_equal(corners["Q2"], (a >= 2) & (b >= -1))    # upper right: both high
    assert np.array_equal(corners["Q1"], (a < 2) & (b >= -1))     # upper left
    assert np.array_equal(corners["Q4"], (a < 2) & (b < -1))      # lower left

    # combinations of sisters, inside their parent
    root = GateNode("s")
    parent = root.add(GateNode("P", RectangleGate(["a"], [-8], [8])))
    parent.add(GateNode("both", BooleanGate(AND, ["X", "Y"])))    # listed before what it uses
    parent.add(GateNode("X", RectangleGate(["a"], [0], [None])))
    parent.add(GateNode("Y", RectangleGate(["b"], [0], [None])))
    parent.add(GateNode("either", BooleanGate(OR, ["X", "Y"])))
    parent.add(GateNode("neither", BooleanGate(NOT, ["X", "Y"])))
    parent.add(GateNode("lost", BooleanGate(AND, ["X", "gone"])))
    m = root.apply(col, 4000)
    p = (a >= -8) & (a < 8)
    assert np.array_equal(m["P/both"], p & (a >= 0) & (b >= 0))
    assert np.array_equal(m["P/either"], p & ((a >= 0) | (b >= 0)))
    assert np.array_equal(m["P/neither"], p & ~((a >= 0) | (b >= 0)))   # never outside the parent
    assert m["P/either"].sum() + m["P/neither"].sum() == p.sum()
    assert not m["P/lost"].any()                                 # made of something that is gone

    # all of it survives a copy and a file
    parent.add(GateNode("E", EllipseGate(["a", "b"], [0.5, 0.5], [0.2, 0.1], 0.3, scales=[t, None])))
    parent.add(GateNode("Q", QuadrantGate(["a", "b"], 1.0, 2.0, "Q3", "g1")))
    for other in (copy_tree(root), tree_from_dict(tree_to_dict(root))):
        assert tree_to_dict(other) == tree_to_dict(root)
        again = other.apply(col, 4000)
        assert all(np.array_equal(again[k], v) for k, v in root.apply(col, 4000).items())
    e = tree_from_dict(tree_to_dict(root)).find("P/E").gate
    assert isinstance(e.scales[0], Biex) and e.scales[1] is None and e.angle == 0.3


# --- scales, compensation, export ------------------------------------------------

def test_scales_as_values():
    from milc.transforms import (Arcsinh, Biex, Linear, Log, decade_ticks, display_range,
                                 from_spec, same, to_spec)
    log = Log(10)
    x = np.array([10.0, 100, 1000, 262144])
    assert log(x) == pytest.approx(np.log10(x / 10) / np.log10(26214.4))
    assert np.allclose(log.inverse(log(x)), x) and log(-50) == 0 and log(0) == 0
    for t in (Biex(), Biex(256, 262144, 1.0, -400, 5.0), Arcsinh(25), Log(3), Linear(0, 1000)):
        assert same(from_spec(to_spec(t)), t)
    assert to_spec(None) is None and from_spec(None) is None and from_spec({"what": 1}) is None
    assert not same(Biex(), Biex(width=-400)) and not same(Biex(), Arcsinh(150))
    assert display_range(Biex()) == (0.0, 1.0) and display_range(None, 1000) == (0.0, 1000.0)
    lo, hi = display_range(Arcsinh(150))
    assert lo < 0 < hi == pytest.approx(np.arcsinh(262144 / 150))
    # a logarithmic axis has decades and no zero; the others begin at zero
    assert [t for _, t in decade_ticks(Log(10))] == ["10²", "10³", "10⁴", "10⁵"]
    assert decade_ticks(Arcsinh(150))[0] == (0.0, "0")
    # a wider linear part pushes the decades up the axis
    assert Biex(width=-1000)(0.0) > Biex(width=-100)(0.0)


def _two_samples(tmp_path):
    from milc.session import Session
    rng = np.random.default_rng(3)
    for name, n in (("3Ctrl1", 500), ("3Case1", 800)):
        ev = np.column_stack([rng.uniform(0, 1000, n), np.abs(rng.lognormal(6, 1.2, (n, 2)))])
        _write(tmp_path, name + ".fcs", ev.astype(np.float32), markers=["", "CD8", "CD4"],
               keywords={"SPILL": "2,FITC-A,PE-A,1,0.1,0.02,1"})
    s = Session()
    s.open_folder(str(tmp_path))
    return s


def test_scale_of_a_parameter(tmp_path):
    from milc.session import Session
    from milc.transforms import Arcsinh, Biex, same
    s = _two_samples(tmp_path)
    assert s.fluorescence() == {"FITC-A", "PE-A"}
    assert same(s.scale_of("Comp-PE-A"), Biex()) and s.scale_of("FSC-A") is None
    # a polygon drawn on the usual scale…
    t = Biex()
    box = np.array([[0.4, 0.4], [0.9, 0.4], [0.9, 0.9], [0.4, 0.9]])
    raw = np.column_stack([t.inverse(box[:, 0]), t.inverse(box[:, 1])])
    for e in s.usable():
        e.gates.add(GateNode("P", PolygonGate(["Comp-PE-A", "Comp-FITC-A"], raw, scales=[t, t])))
    first = s.usable()[0]
    s.gates_changed(first)                      # drawn, so kept beside the files
    before = s.load(first).count("P")
    # …follows its axis onto another: same corners, sides straight somewhere else
    wide = Biex(width=-2000)
    s.set_scale(["Comp-PE-A"], wide)
    gate = first.gates.find("P").gate
    assert same(gate.scales[0], wide) and same(gate.scales[1], t) and same(s.scale_of("PE-A"), wide)
    sample = s.load(first)
    pe, fitc = wide(sample.column("Comp-PE-A")), t(sample.column("Comp-FITC-A"))
    lo, hi = wide(raw[0, 0]), wide(raw[1, 0])
    assert sample.count("P") == ((pe > lo) & (pe < hi) & (fitc > 0.4) & (fitc < 0.9)).sum()
    assert sample.count("P") == before          # a box with sides along the axes is the same box
    s.set_scale(["PE-A"], None)                 # linear: no scale at all on that side
    assert first.gates.find("P").gate.scales[0] is None and s.scale_of("PE-A") is None
    s.set_scale(["FITC-A"], Arcsinh(150))
    # it is all there next time, and the usual one comes back when asked
    again = Session()
    again.open_folder(str(tmp_path))
    assert again.scale_of("PE-A") is None and same(again.scale_of("FITC-A"), Arcsinh(150))
    assert "PE-A" in again.scales
    again.reset_scale(["PE-A", "FITC-A"])
    assert again.scales == {} and same(again.scale_of("Comp-PE-A"), Biex())
    assert same(again.usable()[0].gates.find("P").gate.scales[0], Biex())


def test_features_follow_the_scale(tmp_path):
    from milc import sampling
    from milc.transforms import Arcsinh, Biex
    s = _two_samples(tmp_path)
    s.build()
    ds = s.dataset
    x = ds.column("PE-A").astype(np.float64)
    assert np.allclose(s.features(["PE-A"])[:, 0], Biex()(x) * 1023)      # as FlowJo feeds it
    s.set_scale(["PE-A"], Biex(width=-1000))
    assert np.allclose(s.features(["PE-A"])[:, 0], Biex(width=-1000)(x) * 1023)
    s.set_scale(["PE-A"], Arcsinh(150))
    assert np.allclose(s.features(["PE-A"])[:, 0], np.arcsinh(x / 150))
    assert np.allclose(s.features(["PE-A"], "arcsinh", 5)[:, 0], np.arcsinh(x / 5))


def test_compensation_can_be_edited(tmp_path):
    from milc.compensation import Compensation
    from milc.session import Session
    s = _two_samples(tmp_path)
    matrix, source = s.current_matrix()
    assert source == "file" and matrix.matrix[0, 1] == 0.1
    first = s.usable()[0]
    raw = s.load(first).fcs.events.astype(np.float64)
    usual = s.load(first).compensated.copy()
    edited = Compensation(["FITC-A", "PE-A"], np.array([[1, 0.3], [0.02, 1]]), "edited")
    s.set_matrix(edited)
    assert s.current_matrix()[1] == "edited"
    now = s.load(first).compensated
    assert np.allclose(now[:, 1:], raw[:, 1:] @ np.linalg.inv(edited.matrix), rtol=1e-4)
    assert not np.allclose(now, usual) and np.array_equal(now[:, 0], usual[:, 0])   # scatter untouched
    again = Session()
    again.open_folder(str(tmp_path))
    assert again.current_matrix()[1] == "edited" and again.matrix.matrix[0, 1] == 0.3
    again.set_matrix(None)
    assert again.current_matrix()[1] == "file"
    assert np.allclose(again.load(again.usable()[0]).compensated, usual)


def test_export_a_population(tmp_path):
    from milc.cli import main
    s = _two_samples(tmp_path)
    for e in s.usable():
        e.gates.add(GateNode("Big cells", RectangleGate(["FSC-A"], [600], [None])))
    out = tmp_path / "out"
    written = s.export_population(s.usable(), "Big cells", str(out))
    assert [os.path.basename(p) for p in written] == ["3Case1_Big cells.fcs", "3Ctrl1_Big cells.fcs"]
    sample = s.load(s.usable()[0])
    mask = sample.mask("Big cells")
    f = read_fcs(written[0])
    assert f.channels == ["FSC-A", "Comp-FITC-A", "Comp-PE-A"] and f.markers == ["", "CD8", "CD4"]
    assert np.allclose(f.events, sample.compensated[mask]) and len(f) == mask.sum()
    assert f.keywords["MILC POPULATION"] == "Big cells" and "$SPILLOVER" not in f.keywords
    # as acquired: the events untouched and the matrix beside them
    raw = s.export_population(s.usable()[:1], "Big cells", str(tmp_path / "raw"), compensated=False)
    g = read_fcs(raw[0])
    assert g.channels == ["FSC-A", "FITC-A", "PE-A"] and np.allclose(g.events, sample.fcs.events[mask])
    from milc.compensation import from_keywords
    assert np.allclose(from_keywords(g.keywords).matrix, [[1, 0.1], [0.02, 1]])
    # as a table
    table = s.export_population(s.usable()[:1], "Big cells", str(tmp_path / "csv"), fmt="csv")
    lines = open(table[0], encoding="utf-8-sig").read().splitlines()
    assert lines[0] == "FSC-A,CD8,CD4" and len(lines) == mask.sum() + 1
    assert float(lines[1].split(",")[1]) == pytest.approx(sample.compensated[mask][0, 1], rel=1e-5)
    # every event when no population is named, and from the command line
    assert len(read_fcs(s.export_population(s.usable()[:1], "", str(tmp_path / "all"))[0])) == 800
    s.gates_changed(s.usable()[0])
    assert main(["export", str(tmp_path), "--population", "Big cells", "--format", "csv",
                 "--out", str(tmp_path / "cli")]) == 0
    assert sorted(os.listdir(tmp_path / "cli")) == ["3Case1_Big cells.csv", "3Ctrl1_Big cells.csv"]
    assert main(["export", str(tmp_path), "--population", "Nothing"]) == 2


# --- comparing groups -------------------------------------------------------------

def test_tests_between_groups():
    from scipy import stats as st

    from milc import compare as c
    a, b, d = [1.0, 2, 3, 4, 5], [4.0, 5, 6, 7, 9], [9.0, 10, 11]
    name, stat, p, note = c.test_groups([a, b])
    assert name == "Welch t" and p == pytest.approx(st.ttest_ind(a, b, equal_var=False).pvalue)
    assert c.test_groups([a, b], c.RANK)[2] == pytest.approx(
        st.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    assert c.test_groups([a, b, d])[0] == "ANOVA"
    assert c.test_groups([a, b, d])[2] == pytest.approx(st.f_oneway(a, b, d).pvalue)
    assert c.test_groups([a, b, d], c.RANK)[2] == pytest.approx(st.kruskal(a, b, d).pvalue)
    # three against three cannot go below 0.1 by ranks, and it says so
    name, _, p, note = c.test_groups([[1.0, 2, 3], [7.0, 8, 9]], c.RANK)
    assert p == pytest.approx(0.1) and "0.1" in note
    # what cannot be tested is no p value, with the reason
    for samples, why in (([a], "one group"), ([a, [3.0]], "fewer than 2"),
                         ([[2.0, 2], [2.0, 2]], "same value"),
                         ([[1.0, 1], [2.0, 2]], "no variation")):
        _, _, p, note = c.test_groups(samples)
        assert np.isnan(p) and why in note
    assert c.test_groups([[1.0, np.nan, 3, 5], [7.0, 8, 9]])[2] == pytest.approx(
        st.ttest_ind([1.0, 3, 5], [7.0, 8, 9], equal_var=False).pvalue)   # a missing value is skipped

    # corrections: against scipy where it has one, by hand where not
    p = [0.01, 0.04, 0.03, 0.005]
    assert np.allclose(c.adjust(p, c.BH), st.false_discovery_control(p))
    assert np.allclose(c.adjust(p, c.BONFERRONI), [0.04, 0.16, 0.12, 0.02])
    assert np.allclose(c.adjust(p, c.HOLM), [0.03, 0.06, 0.06, 0.02])
    assert np.allclose(c.adjust(p, c.NONE), p)
    assert np.allclose(c.adjust([0.5, 0.9], c.BONFERRONI), [1.0, 1.0])     # never above 1
    mixed = c.adjust([0.01, np.nan, 0.04], c.BONFERRONI)                   # a test not done
    assert np.isnan(mixed[1]) and np.allclose(mixed[[0, 2]], [0.02, 0.08])  # is not counted
    assert all(np.all(np.diff(c.adjust(sorted(p), m)) >= 0) for m in c.CORRECTIONS)

    # many things at once, and every pair
    found = c.compare({("1",): {"7": a, "14": b, "3": d}, ("2",): {"7": a, "14": a, "3": a}},
                      pairwise=True)
    assert len(found) == 8 and found[0].groups == ["3", "7", "14"]        # days in order
    assert [f.test for f in found[:4]] == ["ANOVA", "Welch t", "Welch t", "Welch t"]
    assert found[0].n == [3, 5, 5] and found[0].mean[1] == 3
    assert np.allclose([f.p_adjusted for f in found], c.adjust([f.p for f in found], c.BH))
    assert all(f.p_adjusted >= f.p for f in found)
    table = c.to_table(found, ["Cluster"])
    assert table.columns[:3] == ["Cluster", "Groups", "Mean ± SD (n)"]
    assert table.rows[0][1] == "3 vs 7 vs 14" and "(3)" in table.rows[0][2]
    assert c.stars(0.03) == "*" and c.stars(0.2) == "ns" and c.stars(float("nan")) == ""
    assert c.label(0.034) == "p = 0.03" and c.label(0.0041) == "p = 0.0041"
    assert c.label(1e-6) == "p < 0.0001"


def test_comparing_a_statistics_table():
    from milc import compare as c
    from milc.stats import StatsTable
    t = StatsTable(["Sample", "Group", "Day", "Population", "Path", "Count", "% of parent"])
    rng = np.random.default_rng(0)
    for day in ("3", "7"):
        for group, centre in (("Ctrl", 10), ("Case", 30 if day == "7" else 10)):
            for animal in range(5):
                value = float(rng.normal(centre, 2))
                t.rows.append([f"{day}{group}{animal}", group, day, "T", "Cells/T", 100, value])
                t.rows.append([f"{day}{group}{animal}", group, day, "B", "Cells/B", 100, 50.0])
    values, names = c.from_table(t, "% of parent", ["Group"])
    assert names == ["Population"] and set(values) == {("T",), ("B",)}
    assert len(values[("T",)]["Case"]) == 10
    values, names = c.from_table(t, "% of parent", ["Group"], within="Day")
    assert names == ["Day", "Population"] and len(values[("7", "T")]["Ctrl"]) == 5
    found = {f.key: f for f in c.compare(values)}
    assert found[("7", "T")].p_adjusted < 0.001 and found[("3", "T")].p_adjusted > 0.05
    assert np.isnan(found[("3", "B")].p) and "same value" in found[("3", "B")].note
    # two columns together are the groups
    values, _ = c.from_table(t, "% of parent", ["Group", "Day"])
    assert set(values[("T",)]) == {"Ctrl · 3", "Ctrl · 7", "Case · 3", "Case · 7"}


# --- Logicle ---------------------------------------------------------------------

#: FlowUtils 1.2.2, `flowutils.transforms.logicle` — the C code of Moore &
#: Parks (Cytometry A 2012) — at these values, for (T, W, M, A)
LOGICLE_AT = [-1000, -100, -10, -1, 0, 1, 10, 100, 1000, 10000, 100000, 262144, 1e6]
LOGICLE_REFERENCE = {
    (262144, 0.5, 4.5, 0): [
        -0.23211535395017646, 0.009041134692025388, 0.099917946544774, 0.1099906900730379,
        0.1111111111111111, 0.11223153214918431, 0.12230427567744821, 0.21318108753019682,
        0.45433757617239867, 0.6838326572265573, 0.9069275914814589, 1.0, 1.1292427085676322],
    (262144, 1.0, 4.5, 0): [
        0.03813323848718597, 0.1996489705135676, 0.2199597782467301, 0.22199597288068995,
        0.2222222222222222, 0.22244847156375447, 0.22448466619771432, 0.24479547393087683,
        0.40631120595725845, 0.6762208182193373, 0.9063558264967823, 1.0, 1.129519464719908],
    (262144, 0.5, 4.5, 0.5): [
        -0.10890381855515885, 0.10813702122282293, 0.18992615189029663, 0.19899162106573412,
        0.2, 0.2010083789342659, 0.2100738481097034, 0.2918629787771771, 0.5089038185551589,
        0.7154493915039016, 0.9162348323333129, 1.0, 1.116318437710869],
    (10000, 0, 4, 0): [
        -0.7500001074877236, -0.500010854648286, -0.25106980177441446, -0.052246909576938755,
        0.0, 0.052246909576938755, 0.25106980177441446, 0.500010854648286, 0.7500001074877236,
        1.0, 1.2499999989251211, 1.354634979403759, 1.4999999989143724],
    (4194304, 0.8, 5.5, 1.0): [
        0.1391642638587352, 0.2584957251415965, 0.27507300715869754, 0.2767380628539924,
        0.27692307692307694, 0.2771080909921615, 0.27877314668745634, 0.29535042870455736,
        0.4146818899874187, 0.5922874485753206, 0.7498848736020712, 0.814570792022158,
        0.9041666523419553],
    (1000, 1, 4, 1): [
        -0.19999999999999996, 0.008361879241263415, 0.25405935894582443, 0.3830012926556203,
        0.4, 0.41699870734437977, 0.5459406410541756, 0.7916381207587366, 1.0,
        1.2010029960319952, 1.4011140722863589, 1.4848300871320008, 1.601125824523326],
}


@pytest.mark.parametrize("parameters", list(LOGICLE_REFERENCE))
def test_logicle_reference(parameters):
    from milc.transforms import Logicle
    t = Logicle(*parameters)
    x = np.array(LOGICLE_AT, dtype=np.float64)
    assert np.abs(t(x) - LOGICLE_REFERENCE[parameters]).max() < 1e-12
    assert np.allclose(t.inverse(LOGICLE_REFERENCE[parameters]), x, rtol=1e-9, atol=1e-9)


def test_logicle_shape_and_round_trip():
    from milc.transforms import Logicle
    t = Logicle()
    T, W, M = 262144.0, 0.5, 4.5
    # the published properties: T at the top, zero at W / M, mirrored about zero
    assert float(t(T)) == pytest.approx(1.0, abs=1e-14) and float(t(0.0)) == pytest.approx(W / M)
    x = np.concatenate([-np.logspace(-3, 6, 400), [0.0], np.logspace(-3, 8, 600)])
    y = t(x)
    assert np.all(np.diff(y[np.argsort(x)]) > 0)
    assert np.allclose(t(-x) + y, 2 * W / M, atol=1e-13)
    assert np.allclose(t.inverse(y), x, rtol=1e-11, atol=1e-9)
    ys = np.linspace(-0.3, 1.4, 5001)
    assert np.abs(t(t.inverse(ys)) - ys).max() < 1e-12
    # far above zero it is the logarithm it replaces: a decade is 1 / M
    assert float(t(1e5) - t(1e4)) == pytest.approx(1 / M, abs=2e-3)
    # the shape of what goes in comes out, single precision and scalars too
    assert t(np.zeros((3, 2), np.float32)).shape == (3, 2) and float(t(5.0)) > W / M
    assert np.isnan(t(np.array([np.nan, 1.0]))[0])
    for bad in ((0, 0.5, 4.5, 0), (262144, 3, 4.5, 0), (262144, 0.5, 4.5, -1), (262144, 1, 4, 3)):
        with pytest.raises(ValueError):
            Logicle(*bad)


def test_logicle_is_a_scale_like_the_others(tmp_path):
    from milc.session import Session
    from milc.transforms import (Biex, Logicle, decade_ticks, display_range, from_spec, same,
                                 to_spec)
    t = Logicle(262144, 1.0, 4.5, 0.5)
    assert to_spec(t) == {"logicle": [262144.0, 1.0, 4.5, 0.5]} and same(from_spec(to_spec(t)), t)
    assert not same(t, Logicle()) and not same(Logicle(), Biex())
    assert display_range(t) == (0.0, 1.0)
    ticks = decade_ticks(Logicle())
    assert ticks[0] == (pytest.approx(0.5 / 4.5), "0") and ticks[-1][1] == "10⁵"
    assert [label for _, label in ticks[1:]] == ["10²", "10³", "10⁴", "10⁵"]

    # a polygon drawn on it is tested on it, and the scale is kept with the project
    s = _two_samples(tmp_path)
    s.set_scale(["PE-A", "FITC-A"], t)
    box = np.array([[0.4, 0.4], [0.9, 0.4], [0.9, 0.9], [0.4, 0.9]])
    raw = np.column_stack([t.inverse(box[:, 0]), t.inverse(box[:, 1])])
    first = s.usable()[0]
    for e in s.usable():
        e.gates.add(GateNode("P", PolygonGate(["Comp-PE-A", "Comp-FITC-A"], raw, scales=[t, t])))
    s.gates_changed(first)
    sample = s.load(first)
    pe, fitc = t(sample.column("Comp-PE-A")), t(sample.column("Comp-FITC-A"))
    assert sample.count("P") == ((pe > 0.4) & (pe < 0.9) & (fitc > 0.4) & (fitc < 0.9)).sum() > 0
    again = Session()
    again.open_folder(str(tmp_path))
    assert same(again.scale_of("Comp-PE-A"), t)
    assert same(again.usable()[0].gates.find("P").gate.scales[0], t)
    assert again.load(again.usable()[0]).count("P") == sample.count("P")

    # the algorithms are fed it: the marker's own, or the standard one for every marker
    again.build()
    x = again.dataset.column("PE-A").astype(np.float64)
    assert np.allclose(again.features(["PE-A"])[:, 0], t(x) * 1023)
    assert np.allclose(again.features(["PE-A"], "logicle")[:, 0], Logicle()(x) * 1023)


def test_analyze_takes_logicle(tmp_path):
    from milc.cli import main
    from milc.dataset import Dataset
    from milc.transforms import Logicle
    rng = np.random.default_rng(2)
    ev = np.column_stack([np.abs(rng.lognormal(6, 1.5, (600, 2))),
                          np.tile([1, 2, 3], 200)]).astype(np.float32)
    path = str(tmp_path / "c.fcs")
    write_fcs(path, ev, ["FITC-A", "PE-A", "SampleID"], markers=["CD8", "CD4", ""])
    out = str(tmp_path / "run")
    assert main(["analyze", path, "--channels", "FITC-A,PE-A", "--transform", "logicle",
                 "--pca", "--out", out]) == 0
    ds = Dataset.load(out + ".npz")
    assert ds.runs["PCA"]["transform"] == "logicle" and "cofactor" not in ds.runs["PCA"]
    X = ds.matrix(["FITC-A", "PE-A"], Logicle(), 1023)
    centred = X - X.mean(axis=0)
    first = np.linalg.svd(centred, full_matrices=False)[2][0]
    assert np.allclose(np.abs(ds.derived["PCA1"]), np.abs(centred @ first), atol=1e-3 * 1023)


def test_a_workspace_logicle_is_read(tmp_path):
    from milc.transforms import Logicle, same
    from milc.wsp import read_wsp
    text = """<Workspace version="20.0" flowJoVersion="10.8.1"
 xmlns:transforms="http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
 xmlns:data-type="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes"
 xmlns:gating="http://www.isac-net.org/std/Gating-ML/v2.0/gating">
 <SampleList><Sample><DataSet uri="file:/x/a.fcs" sampleID="1"/>
 <Transformations>
  <transforms:logicle transforms:length="256" transforms:T="262144" transforms:A="0"
     transforms:W="1" transforms:M="4.5"><data-type:parameter data-type:name="Comp-PE-A"/>
  </transforms:logicle>
  <transforms:logicle transforms:T="262144" transforms:A="0" transforms:W="9"
     transforms:M="4.5"><data-type:parameter data-type:name="Comp-APC-A"/></transforms:logicle>
  <transforms:hyperlog transforms:T="262144"><data-type:parameter data-type:name="Comp-X-A"/>
  </transforms:hyperlog>
 </Transformations>
 <SampleNode name="a.fcs" sampleID="1" count="10"><Subpopulations>
  <Population name="P" count="4"><Gate><gating:PolygonGate>
   <gating:dimension><data-type:fcs-dimension data-type:name="Comp-PE-A"/></gating:dimension>
   <gating:dimension><data-type:fcs-dimension data-type:name="FSC-A"/></gating:dimension>
   <gating:vertex><gating:coordinate data-type:value="0"/><gating:coordinate data-type:value="0"/></gating:vertex>
   <gating:vertex><gating:coordinate data-type:value="1000"/><gating:coordinate data-type:value="0"/></gating:vertex>
   <gating:vertex><gating:coordinate data-type:value="1000"/><gating:coordinate data-type:value="9"/></gating:vertex>
  </gating:PolygonGate></Gate></Population>
 </Subpopulations></SampleNode></Sample></SampleList></Workspace>"""
    path = tmp_path / "w.wsp"
    path.write_text(text, encoding="utf-8")
    sample = read_wsp(str(path)).samples[0]
    assert same(sample.transforms["Comp-PE-A"], Logicle(262144, 1, 4.5, 0))
    assert set(sample.transforms) == {"Comp-PE-A"}        # an impossible W, an unknown scale: left out
    gate = sample.gates.find("P").gate
    assert same(gate.scales[0], Logicle(262144, 1, 4.5, 0)) and gate.scales[1] is None


# --- a compensation matrix per sample ----------------------------------------------

def _own_spill(tmp_path):
    """Two files that carry different matrices, and a third like the first."""
    from milc.session import Session
    rng = np.random.default_rng(5)
    for name, spill in (("3Ctrl1", 0.1), ("3Case1", 0.3), ("3Ctrl2", 0.1)):
        ev = np.column_stack([rng.uniform(0, 1000, 300), np.abs(rng.lognormal(6, 1.2, (300, 2)))])
        _write(tmp_path, name + ".fcs", ev.astype(np.float32), markers=["", "CD8", "CD4"],
               keywords={"SPILL": f"2,FITC-A,PE-A,1,{spill},0.02,1"})
    s = Session()
    s.open_folder(str(tmp_path))
    return s


def _by_name(s):
    return {e.name: e for e in s.usable()}


def test_each_sample_is_compensated_with_its_own_matrix(tmp_path):
    from milc.compensation import Compensation
    from milc.session import Session
    s = _own_spill(tmp_path)
    e = _by_name(s)
    # as acquired: each file's own, and they differ
    assert {n: s.matrix_in_use(x)[1] for n, x in e.items()} == dict.fromkeys(e, "file")
    assert s.matrix_in_use(e["3Case1"])[0].matrix[0, 1] == 0.3 and s.distinct_matrices() == 2
    assert s.load(e["3Case1"]).compensation.matrix[0, 1] == 0.3
    assert s.load(e["3Ctrl1"]).compensation.matrix[0, 1] == 0.1
    usual = {n: s.load(x).compensated.copy() for n, x in e.items()}
    assert not os.path.exists(s.project_file())            # looking writes nothing

    # one sample edited: only that one moves
    edited = Compensation(["FITC-A", "PE-A"], np.array([[1, 0.2], [0.02, 1]]), "edited")
    s.set_matrix(edited, [e["3Ctrl1"]])
    assert s.matrix_in_use(e["3Ctrl1"]) == (edited, "edited")
    assert s.matrix_in_use(e["3Case1"])[1] == "file" and s.distinct_matrices() == 3
    raw = s.load(e["3Ctrl1"]).fcs.events.astype(np.float64)
    assert np.allclose(s.load(e["3Ctrl1"]).compensated[:, 1:],
                       raw[:, 1:] @ np.linalg.inv(edited.matrix), rtol=1e-4)
    assert np.array_equal(s.load(e["3Case1"]).compensated, usual["3Case1"])
    assert s.matrix is None and list(s.own_matrices) == ["3Ctrl1.fcs"]

    # kept with the project
    again = Session()
    again.open_folder(str(tmp_path))
    a = _by_name(again)
    assert again.matrix_in_use(a["3Ctrl1"])[1] == "edited"
    assert again.matrix_of(a["3Ctrl1"]).matrix[0, 1] == 0.2
    assert again.matrix_in_use(a["3Ctrl2"])[1] == "file"
    assert np.allclose(again.load(a["3Ctrl1"]).compensated, s.load(e["3Ctrl1"]).compensated)

    # to all: one shared matrix in place of every other, the edited one included
    shared = Compensation(["FITC-A", "PE-A"], np.array([[1, 0.5], [0, 1]]), "edited")
    again.set_matrix(shared)
    assert {again.matrix_in_use(x)[1] for x in a.values()} == {"shared"}
    assert again.distinct_matrices() == 1 and again.own_matrices == {}
    # one of them back on its file's, another edited on its own: the rest still share
    again.set_matrix(None, [a["3Case1"]])
    again.set_matrix(edited, [a["3Ctrl2"]])
    assert [again.matrix_in_use(a[n])[1] for n in ("3Ctrl1", "3Case1", "3Ctrl2")] == [
        "shared", "file", "edited"]
    assert np.array_equal(again.load(a["3Case1"]).compensated, usual["3Case1"])
    third = Session()
    third.open_folder(str(tmp_path))
    t = _by_name(third)
    assert [third.matrix_in_use(t[n])[1] for n in ("3Ctrl1", "3Case1", "3Ctrl2")] == [
        "shared", "file", "edited"]
    # back to the original is the file's, not the shared one; and once the
    # last one sharing has gone back too, nothing is shared any more
    third.set_matrix(None, [t["3Ctrl2"]])
    assert third.matrix_in_use(t["3Ctrl2"])[1] == "file" and third.matrix is not None
    third.set_matrix(None, [t["3Ctrl1"]])
    assert third.matrix is None and third.own_matrices == {}
    assert {third.matrix_in_use(x)[1] for x in t.values()} == {"file"}
    for n, x in t.items():
        assert np.array_equal(third.load(x).compensated, usual[n])

    # a matrix for channels a file does not have is not given to it
    with pytest.raises(ValueError, match="BV421-A"):
        third.set_matrix(Compensation(["FITC-A", "BV421-A"], np.eye(2)), [t["3Ctrl1"]])
    assert third.own_matrices == {}


def test_a_project_saved_before_matrices_were_per_sample_still_loads(tmp_path):
    """The project file as it was written before: one `compensation`
    for every sample, and no `compensation_samples`."""
    import json
    from milc.session import PROJECT_FILE, Session
    s = _own_spill(tmp_path)
    usual = s.load(s.usable()[0]).compensated.copy()
    old = {"version": 2, "columns": ["Group"], "group_by": ["Group"],
           "annotations": {}, "excluded": [], "gates": {}, "scales": {"PE-A": {"arcsinh": 150.0}},
           "normalization": None,
           "compensation": {"channels": ["FITC-A", "PE-A"], "matrix": [[1, 0.25], [0.02, 1]]}}
    with open(tmp_path / PROJECT_FILE, "w", encoding="utf-8") as f:
        json.dump(old, f)
    s = Session()
    s.open_folder(str(tmp_path))
    assert s.matrix.matrix[0, 1] == 0.25 and s.own_matrices == {}
    assert s.current_matrix()[1] == "edited"                     # what it said before
    for e in s.usable():                                        # every sample, as before
        assert s.matrix_in_use(e)[1] == "shared"
        sample = s.load(e)
        assert np.allclose(sample.compensated[:, 1:], sample.fcs.events[:, 1:].astype(np.float64)
                           @ np.linalg.inv(s.matrix.matrix), rtol=1e-4)
    # and one without any edited matrix: each file's own
    old["compensation"] = None
    with open(tmp_path / PROJECT_FILE, "w", encoding="utf-8") as f:
        json.dump(old, f)
    s = Session()
    s.open_folder(str(tmp_path))
    assert s.matrix is None and {s.matrix_in_use(e)[1] for e in s.usable()} == {"file"}
    assert np.array_equal(s.load(s.usable()[0]).compensated, usual)
    # a mangled entry for one sample does not stop the folder opening
    old["compensation_samples"] = {"3Ctrl1.fcs": {"channels": ["FITC-A"]}, "3Case1.fcs": 7}
    with open(tmp_path / PROJECT_FILE, "w", encoding="utf-8") as f:
        json.dump(old, f)
    s = Session()
    s.open_folder(str(tmp_path))
    assert s.own_matrices == {} and len(s.usable()) == 3
