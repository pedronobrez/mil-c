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
