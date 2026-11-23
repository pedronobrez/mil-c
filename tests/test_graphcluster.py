"""PhenoGraph-style clustering and PCA, against what is known: mixtures
whose populations were drawn, a graph whose best partition is in the
literature, modularity written out from its definition, scikit-learn's
PCA. Then the same through the session, the command line and the window."""

import os

import numpy as np
import pytest
from sklearn.metrics import adjusted_rand_score

from milc import embed
from milc import graphcluster as gc
from milc.fcs import write_fcs

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def _mixture(shares, apart, n=6000, p=8, seed=0):
    """Gaussian populations of unit spread, of these shares, their
    centres `apart` from one another along different markers' axes."""
    rng = np.random.default_rng(seed)
    truth = rng.choice(len(shares), n, p=shares)
    centres = np.zeros((len(shares), p))
    for i in range(len(shares)):
        centres[i, i % p] = apart / np.sqrt(2) * (1 + i // p)
    return centres[truth] + rng.normal(0, 1, (n, p)), truth


SHARES = (0.40, 0.25, 0.18, 0.10, 0.06, 0.01)


def _dense_modularity(A, labels, resolution=1.0):
    """The definition, pair by pair: (1/2m) sum (A_ij - g k_i k_j / 2m)
    over the pairs that share a community."""
    k = A.sum(axis=1)
    two_m = k.sum()
    same = labels[:, None] == labels[None, :]
    return float(((A - resolution * np.outer(k, k) / two_m) * same).sum() / two_m)


def _ring_of_cliques(cliques, size):
    from scipy import sparse
    n = cliques * size
    A = np.zeros((n, n))
    for c in range(cliques):
        A[c * size:(c + 1) * size, c * size:(c + 1) * size] = 1.0
        a, b = c * size, ((c + 1) % cliques) * size + 1
        A[a, b] = A[b, a] = 1.0
    np.fill_diagonal(A, 0.0)
    return sparse.csr_matrix(A), np.repeat(np.arange(cliques), size)


# --- the graph ---------------------------------------------------------------

def test_neighbours_are_the_nearest():
    rng = np.random.default_rng(0)
    X = rng.normal(0, 1, (400, 5))
    X[10] = X[11] = X[12]  # identical events
    near = gc.neighbours(X, 7)
    d = ((X[:, None, :] - X[None, :, :]) ** 2).sum(axis=2)
    np.fill_diagonal(d, np.inf)
    assert near.shape == (400, 7)
    assert not (near == np.arange(400)[:, None]).any()
    kth = np.sort(d, axis=1)[:, 6]
    assert np.all(d[np.arange(400)[:, None], near] <= kth[:, None] + 1e-12)
    assert gc.neighbours(X[:5], 30).shape == (5, 4)  # never more than there are


def test_approximate_neighbours(monkeypatch):
    """Beyond 50,000 events on more than ten markers the search is
    approximate: nearly all of the true neighbours, and the same ones
    for the same seed."""
    rng = np.random.default_rng(6)
    X = rng.normal(0, 1, (3000, 12)) + 6 * rng.integers(0, 2, (3000, 12))
    exact = gc.neighbours(X, 15)
    monkeypatch.setattr(gc, "EXACT_EVENTS", 0)
    approximate = gc.neighbours(X, 15, seed=1)
    assert approximate.shape == (3000, 15)
    assert not (approximate == np.arange(3000)[:, None]).any()
    found = np.mean([len(np.intersect1d(a, b)) for a, b in zip(exact, approximate)]) / 15
    assert found > 0.97
    assert np.array_equal(approximate, gc.neighbours(X, 15, seed=1))


def test_jaccard_weights():
    rng = np.random.default_rng(1)
    near = gc.neighbours(rng.normal(0, 1, (300, 4)), 10)
    G = gc.jaccard_graph(near)
    assert (G != G.T).nnz == 0 and G.diagonal().sum() == 0
    hoods = [set(row) | {i} for i, row in enumerate(near)]
    dense = G.toarray()
    for i in range(300):
        for j in range(300):
            joined = i != j and (j in hoods[i] or i in hoods[j])
            expected = len(hoods[i] & hoods[j]) / len(hoods[i] | hoods[j]) if joined else 0.0
            assert dense[i, j] == pytest.approx(expected)
    # an event is in its own neighbourhood: no join weighs nothing
    assert dense[np.arange(300)[:, None], near].min() > 0


# --- the communities ---------------------------------------------------------

def test_modularity_is_the_definition():
    from scipy import sparse
    rng = np.random.default_rng(2)
    A = np.triu(rng.random((60, 60)) * (rng.random((60, 60)) < 0.2), 1)
    A = A + A.T
    labels = rng.integers(0, 5, 60)
    for resolution in (1.0, 0.5, 2.0):
        assert gc.modularity(sparse.csr_matrix(A), labels, resolution) == \
            pytest.approx(_dense_modularity(A, labels, resolution), abs=1e-12)
    assert gc.modularity(sparse.csr_matrix(A), np.zeros(60, int)) == pytest.approx(0.0)


def test_ring_of_cliques():
    """Ten cliques of five in a ring: the cliques are the optimum, at
    Q = 1 - 2/(m(m-1)+2) - 1/n (Fortunato & Barthelemy, PNAS 2007)."""
    G, cliques = _ring_of_cliques(10, 5)
    for seed in range(5):
        member, q = gc.louvain(G, seed=seed)
        assert adjusted_rand_score(cliques, member) == 1.0
        assert q == pytest.approx(1 - 2 / 22 - 1 / 10)
        assert q == pytest.approx(_dense_modularity(G.toarray(), member))


def test_resolution_limit_and_the_way_past_it():
    """Thirty such cliques: modularity is higher with the cliques paired
    than apart — the resolution limit of the same paper — and Louvain
    should find that; a higher resolution brings the cliques back."""
    G, cliques = _ring_of_cliques(30, 5)
    apart, paired = 1 - 2 / 22 - 1 / 30, 1 - 1 / 22 - 2 / 30
    assert gc.modularity(G, cliques) == pytest.approx(apart)
    for seed in range(5):
        member, q = gc.louvain(G, seed=seed)
        assert apart < q <= paired + 1e-12
        # no clique is split
        assert all(len(set(member[cliques == c])) == 1 for c in range(30))
        assert len(set(member)) < 30
        sharp, _ = gc.louvain(G, resolution=2.0, seed=seed)
        assert adjusted_rand_score(cliques, sharp) == 1.0


def test_graph_without_joins():
    from scipy import sparse
    member, q = gc.louvain(sparse.csr_matrix((4, 4)))
    assert member.tolist() == [0, 1, 2, 3] and q == 0.0


# --- mixtures ----------------------------------------------------------------

def test_well_separated_populations():
    X, truth = _mixture(SHARES, apart=10)
    r = gc.phenograph(X, k=30, seed=0)
    assert adjusted_rand_score(truth, r.labels) > 0.99
    assert r.n_clusters == 6 and r.k == 30
    sizes = np.bincount(r.labels)[1:]
    assert r.labels.min() == 1 and np.all(np.diff(sizes) <= 0)  # 1 is the largest
    # the rare one, 1% of the events, is a cluster of its own
    assert set(r.labels[truth == 5]) == {6} and sizes[5] == (truth == 5).sum()
    graph = gc.jaccard_graph(gc.neighbours(X, 30))
    assert r.modularity == pytest.approx(gc.modularity(graph, truth), abs=1e-3)
    assert r.modularity == pytest.approx(
        _dense_modularity(graph[:1500, :1500].toarray(), r.labels[:1500]), abs=0.05)


def test_same_seed_same_clusters_and_other_seeds_agree():
    X, truth = _mixture(SHARES, apart=10, n=3000)
    first = gc.phenograph(X, seed=3)
    assert np.array_equal(first.labels, gc.phenograph(X, seed=3).labels)
    Xo, _ = _mixture(SHARES, apart=4.5, n=3000)
    runs = [gc.phenograph(Xo, seed=s).labels for s in range(4)]
    for other in runs[1:]:
        assert adjusted_rand_score(runs[0], other) > 0.85


def _purity(truth, labels):
    """The share of events whose cluster is mostly their population."""
    table = np.zeros((labels.max() + 1, truth.max() + 1))
    np.add.at(table, (labels, truth), 1)
    return table.max(axis=1).sum() / len(truth)


def test_overlapping_populations():
    """Centres 6 SD apart, where neighbouring populations begin to share
    events: the truth is still found, the rare population (1%, 60
    events) with it. At 4.5 SD the truth is no longer the best partition
    of the graph: the rare population is taken into a neighbour and the
    largest is cut. What can be asked there is that the partition is at
    least as good as the truth by the measure being optimised, and that
    it follows the truth as far as it goes."""
    X, truth = _mixture(SHARES, apart=6.0)
    r = gc.phenograph(X, k=30, seed=0)
    assert r.n_clusters == 6 and adjusted_rand_score(truth, r.labels) > 0.97
    rare = truth == 5
    found = r.labels == np.bincount(r.labels[rare]).argmax()
    assert (found & rare).sum() / rare.sum() > 0.9       # recall
    assert (found & rare).sum() / found.sum() > 0.9      # precision

    X, truth = _mixture(SHARES, apart=4.5)
    r = gc.phenograph(X, k=30, seed=0)
    graph = gc.jaccard_graph(gc.neighbours(X, 30))
    assert r.modularity >= gc.modularity(graph, truth)
    assert adjusted_rand_score(truth, r.labels) > 0.6
    # fewer, larger clusters at a lower resolution; more at a higher
    low = gc.phenograph(X, k=30, resolution=0.5, seed=0)
    high = gc.phenograph(X, k=30, resolution=3.0, seed=0)
    assert low.n_clusters < r.n_clusters < high.n_clusters
    assert adjusted_rand_score(truth, low.labels) > 0.85


def test_a_large_population_is_split_at_resolution_one():
    """What modularity does and PhenoGraph is known for: with enough
    events, cutting a large uniform population in two scores higher
    than leaving it whole. The pieces are pure, and a lower resolution
    puts them together."""
    X, truth = _mixture(SHARES, apart=10, n=20000)
    r = gc.phenograph(X, seed=0)
    graph = gc.jaccard_graph(gc.neighbours(X, 30))
    assert r.n_clusters > 6 and r.modularity > gc.modularity(graph, truth)
    assert _purity(truth, r.labels) > 0.999
    whole = gc.phenograph(X, resolution=0.5, seed=0)
    assert whole.n_clusters == 6 and adjusted_rand_score(truth, whole.labels) > 0.999


# --- PCA ---------------------------------------------------------------------

def test_pca_against_scikit_learn():
    from sklearn.decomposition import PCA
    rng = np.random.default_rng(4)
    X = rng.normal(0, 1, (500, 6)) @ rng.normal(0, 1, (6, 6)) + rng.normal(5, 1, 6)
    ours, theirs = embed.pca(X, 3), PCA(3, svd_solver="full").fit(X)
    assert ours.scores.shape == (500, 3) and ours.loadings.shape == (3, 6)
    assert np.allclose(ours.explained, theirs.explained_variance_ratio_)
    scores = theirs.transform(X)
    for i in range(3):
        sign = np.sign(ours.loadings[i] @ theirs.components_[i])
        assert np.allclose(ours.loadings[i], sign * theirs.components_[i], atol=1e-8)
        assert np.allclose(ours.scores[:, i], sign * scores[:, i], atol=1e-3)
        assert ours.loadings[i, np.abs(ours.loadings[i]).argmax()] > 0
    assert np.allclose(ours.loadings @ ours.loadings.T, np.eye(3), atol=1e-10)
    assert embed.pca(X, 99).scores.shape == (500, 6)


def test_pca_finds_the_direction_that_was_put_there():
    rng = np.random.default_rng(5)
    t = rng.normal(0, 10, 2000)
    X = np.outer(t, [0.6, 0.8, 0.0]) + rng.normal(0, 1, (2000, 3))
    r = embed.pca(X)
    assert np.allclose(r.loadings[0], [0.6, 0.8, 0.0], atol=0.02)
    assert r.explained[0] > 0.95


# --- session, command line, window ---------------------------------------------

@pytest.fixture
def folder(tmp_path):
    """Four samples of three populations; the last column says which."""
    rng = np.random.default_rng(0)
    centres = np.array([[300, 3000, 100], [3000, 200, 100], [200, 300, 4000.0]])
    for name, n in (("3Ctrl1", 600), ("3Ctrl2", 700), ("3Case1", 900), ("3Case2", 1200)):
        pop = rng.choice(3, n, p=(0.6, 0.3, 0.1))
        fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (n, 3)))
        ev = np.column_stack([rng.normal(90000, 9000, n), fluo, pop]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), ev,
                  ["FSC-A", "FITC-A", "PE-A", "APC-A", "Truth"],
                  markers=["", "CD8", "CD4", "Foxp3", ""],
                  keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
    return str(tmp_path)


def test_session_runs(folder, tmp_path):
    from milc.dataset import Dataset
    from milc.session import Session
    s = Session()
    s.open_folder(folder)
    ds = s.build()
    channels = s.default_channels()
    assert len(ds) == 2400 and channels == ["FITC-A", "PE-A", "APC-A"]
    truth = ds.column("Truth").astype(int)

    # three markers and one population of 60%: at resolution 1 they come
    # out in pieces, each of them of one population; far below 1, whole
    s.run_phenograph(channels)
    labels, run = ds.derived["PhenoGraph"], ds.runs["PhenoGraph"]
    assert (run["k"], run["resolution"], run["seed"]) == (30, 1.0, 0)
    assert run["n_clusters"] == labels.max() > 3 and _purity(truth, labels) > 0.99

    s.run_phenograph(channels, k=40, resolution=0.05, seed=1)
    labels, run = ds.derived["PhenoGraph"], ds.runs["PhenoGraph"]
    assert labels.dtype == np.int32 and labels.min() == 1
    assert run["channels"] == channels and run["transform"] == "biex"
    assert (run["k"], run["resolution"], run["seed"]) == (40, 0.05, 1)
    assert run["n_clusters"] == 3 and 0.3 < run["modularity"] < 1
    assert adjusted_rand_score(truth, labels) > 0.99

    s.run_pca(channels)
    run = ds.runs["PCA"]
    X = s.features(channels)
    assert np.allclose(ds.derived["PCA1"], embed.pca(X).scores[:, 0])
    assert len(run["explained"]) == 2 and np.array(run["loadings"]).shape == (2, 3)
    assert run["explained"][0] > run["explained"][1] > 0
    # the populations are apart along the components too
    means = [ds.derived["PCA1"][truth == i].mean() for i in range(3)]
    assert np.ptp(means) > 3 * ds.derived["PCA1"][truth == 0].std()

    path = str(tmp_path / "kept.npz")
    ds.save(path)
    again = Dataset.load(path)
    assert again.runs["PhenoGraph"] == ds.runs["PhenoGraph"]
    assert again.runs["PCA"] == ds.runs["PCA"]
    assert np.array_equal(again.derived["PhenoGraph"], labels)


def test_analyze_command(tmp_path, capsys):
    from milc.cli import main
    from milc.dataset import Dataset
    rng = np.random.default_rng(0)
    centres = np.array([[200, 200, 600, 200], [600, 200, 200, 200], [200, 600, 200, 600.0]])
    truth = np.repeat(np.arange(3), 200)
    X = centres[truth] + rng.normal(0, 25, (600, 4))
    ev = np.column_stack([X, np.tile([1, 2, 3, 4], 150)]).astype(np.float32)
    path = str(tmp_path / "c.fcs")
    write_fcs(path, ev, ["FITC-A", "PE-A", "APC-A", "BV421-A", "SampleID"],
              markers=["CD8", "CD4", "Foxp3", "", ""])
    out = str(tmp_path / "run")
    assert main(["analyze", path, "--channels", "cd8,CD4,Foxp3,BV421-A", "--transform",
                 "arcsinh", "--phenograph", "--k", "20", "--resolution", "0.5", "--pca",
                 "--out", out]) == 0
    said = capsys.readouterr().out
    assert "PhenoGraph: 3 clusters, modularity" in said and "PCA1 holds" in said
    ds = Dataset.load(out + ".npz")
    assert set(ds.derived) == {"PhenoGraph", "PCA1", "PCA2"}
    assert adjusted_rand_score(truth, ds.derived["PhenoGraph"]) == 1.0
    run = ds.runs["PhenoGraph"]
    assert (run["k"], run["resolution"], run["n_clusters"]) == (20, 0.5, 3)
    assert run["cofactor"] == 150.0 and len(ds.runs["PCA"]["explained"]) == 2

    def lines(name):
        return (tmp_path / name).read_text(encoding="utf-8-sig").splitlines()
    assert lines("run_events.csv")[0] == "sample,event,PhenoGraph,PCA1,PCA2"
    assert lines("run_phenograph_mfi.csv")[0] == "cluster,CD8,CD4,Foxp3,BV421-A"
    assert len(lines("run_phenograph_mfi.csv")) == 4
    assert len(lines("run_phenograph_freq.csv")) == 4
    loadings = lines("run_pca_loadings.csv")
    assert loadings[0] == "component,CD8,CD4,Foxp3,BV421-A" and loadings[1].startswith("PCA1,")
    assert not (tmp_path / "run_mfi.csv").exists()  # that one is FlowSOM's


def test_explore_controls(folder, tmp_path):
    QtWidgets = pytest.importorskip("PyQt6.QtWidgets")
    os.environ["MILC_SETTINGS"] = str(tmp_path / "s.ini")
    from milc.ui import prepare, worker
    from milc.ui.shell import MainShell
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(folder)
    w.samples.concatenate()
    settle()
    ds, e = w.session.dataset, w.explore
    assert e.k.value() == 30 and e.resolution.value() == 1.0
    assert e.color_by.findData("cluster|PhenoGraph") == -1
    assert e.x_axis.findData("derived|PCA1") == -1

    e.k.setValue(40)
    e.resolution.setValue(0.05)
    e.seed.setValue(2)
    e.run(["PhenoGraph"])
    assert not e.run_buttons["PhenoGraph"].isEnabled()
    settle()
    assert e.run_buttons["PhenoGraph"].isEnabled()
    run = ds.runs["PhenoGraph"]
    assert (run["k"], run["resolution"], run["seed"]) == (40, 0.05, 2)
    assert adjusted_rand_score(ds.column("Truth").astype(int), ds.derived["PhenoGraph"]) > 0.99
    assert e.color_by.currentData() == "cluster|PhenoGraph"
    assert e.color_by.currentText() == "PhenoGraph cluster"
    assert e.legend.count() == run["n_clusters"]
    assert f"PhenoGraph: {run['n_clusters']} clusters, modularity" in e.result.text()
    e.split_by.setCurrentIndex(e.split_by.findData("cluster|PhenoGraph"))
    assert len(e.view.panels) == run["n_clusters"]
    assert e.view.content().split.title == "PhenoGraph cluster"
    e.split_by.setCurrentIndex(0)

    e.run(["PCA"])
    settle()
    assert e.x_axis.currentData() == "derived|PCA1" and e.y_axis.currentData() == "derived|PCA2"
    shares = ds.runs["PCA"]["explained"]
    content = e.view.content()
    assert content.xlabel == f"PCA1 ({shares[0]:.0%})"
    assert content.ylabel == f"PCA2 ({shares[1]:.0%})"
    # what the components are made of, by the markers' names
    text = e.result.text()
    assert f"PCA1, {shares[0]:.0%}: " in text and "PCA2, " in text
    assert all(name in text for name in ("CD8", "CD4", "Foxp3"))
    assert e.color_by.currentData() == "cluster|PhenoGraph"  # a new map keeps the colours

    # FlowSOM beside it: both clusterings are offered, each under its name
    e.clusters.setValue(3)
    e.grid.setValue(5)
    e.run(["FlowSOM"])
    settle()
    assert e.color_by.currentData() == "cluster|FlowSOM"
    assert e.split_by.findData("cluster|PhenoGraph") >= 0
    assert e.view.content().color.title == "FlowSOM cluster"

    # a new concatenation takes both away again
    w.samples.concatenate()
    settle()
    assert e.color_by.findData("cluster|PhenoGraph") == -1
    assert e.x_axis.findData("derived|PCA1") == -1 and "PCA" not in e.result.text()
    w.close()
