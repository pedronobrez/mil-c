"""Clusters a person works with: names, merges and hand-drawn clusters on
top of the algorithm's labels, names proposed from expression, and a
cluster as a gate on the full files. Numbers are checked against a
calculation that does not go through the code under test."""

import json
import os

import numpy as np
import pytest

from milc import clusters
from milc.cluster_stats import Table
from milc.dataset import Dataset
from milc.fcs import write_fcs
from milc.gating import (ClusterGate, GateNode, RectangleGate, copy_tree, tree_from_dict,
                         tree_to_dict)
from milc.session import Session
from milc.transforms import Arcsinh, Biex, Linear, Logicle


def dataset(labels, samples=2) -> Dataset:
    labels = np.asarray(labels, dtype=np.int32)
    n = labels.size
    rng = np.random.default_rng(1)
    return Dataset(X=rng.normal(size=(n, 2)).astype(np.float32), channels=["A", "B"],
                   markers=["CD4", "CD8"], compensated=[True, True],
                   sample_index=(np.arange(n) % samples).astype(np.int32),
                   event_index=np.arange(n, dtype=np.int64),
                   samples=[f"s{i}" for i in range(samples)],
                   meta=[{} for _ in range(samples)],
                   derived={"FlowSOM": labels}, runs={"FlowSOM": {"channels": ["A", "B"]}})


# --- names, notes, merges ------------------------------------------------------

def test_untouched_clusters_are_the_algorithms():
    ds = dataset([1, 1, 2, 3, 3, 3])
    now = clusters.view(ds)
    assert now.ids == [1, 2, 3] and now.names == ["1", "2", "3"]
    assert np.array_equal(now.labels, ds.derived["FlowSOM"])
    assert clusters.view(ds, "Other") is None and clusters.view(None) is None


def test_rename_and_note():
    ds = dataset([1, 1, 2, 3, 3, 3])
    clusters.rename(ds, {2: " CD8+ "})
    clusters.set_note(ds, 2, "cytotoxic")
    now = clusters.view(ds)
    assert now.names == ["1", "CD8+", "3"] and now.note(2) == "cytotoxic"
    with pytest.raises(ValueError):
        clusters.rename(ds, {1: "CD8+"})       # a name is one cluster's
    with pytest.raises(ValueError):
        clusters.rename(ds, {1: "3"})          # and so is another's number
    clusters.rename(ds, {1: "CD8+", 2: "CD4+"})  # together, a swap is no clash
    assert clusters.view(ds).names == ["CD8+", "CD4+", "3"]
    clusters.rename(ds, {1: ""})
    assert clusters.view(ds).names == ["1", "CD4+", "3"]


def test_merge_is_a_mapping_and_can_be_undone():
    raw = [1, 1, 2, 3, 3, 3, 4, 2]
    ds = dataset(raw)
    clusters.rename(ds, {3: "gone"})
    assert clusters.merge(ds, [3, 2]) == 2
    now = clusters.view(ds)
    assert now.ids == [1, 2, 4] and now.names == ["1", "2+3", "4"]
    assert now.members == {1: [1], 2: [2, 3], 4: [4]}
    # the union, counted by hand
    assert np.array_equal(now.labels, [1, 1, 2, 2, 2, 2, 4, 2])
    assert np.array_equal(ds.derived["FlowSOM"], raw)     # the algorithm's are untouched
    # a merged cluster merges on: everything it was made of follows
    clusters.rename(ds, {2: "T cells"})
    assert clusters.merge(ds, [4, 2]) == 2
    now = clusters.view(ds)
    assert now.members == {1: [1], 2: [2, 3, 4]} and now.names == ["1", "T cells"]
    clusters.undo_merges(ds)
    now = clusters.view(ds)
    assert now.ids == [1, 2, 3, 4] and now.names == ["1", "2", "3", "4"]
    assert np.array_equal(now.labels, raw)
    with pytest.raises(ValueError):
        clusters.merge(ds, [1])
    with pytest.raises(ValueError):
        clusters.merge(ds, [1, 9])


def test_manual_cluster_takes_its_events_and_gives_them_back():
    raw = np.array([1, 1, 2, 3, 3, 3])
    ds = dataset(raw)
    mask = np.array([False, True, True, False, False, False])
    c = clusters.add_manual(ds, mask, "odd ones")
    now = clusters.view(ds)
    assert c == 4 and now.ids == [1, 3, 4]   # 2 lost its only event
    assert now.names == ["1", "3", "odd ones"] and now.is_manual(4) and not now.is_manual(1)
    assert np.array_equal(now.labels, [1, 4, 4, 3, 3, 3])
    assert np.array_equal(clusters.merged_labels(ds), raw)   # the map's nodes do not move
    with pytest.raises(ValueError):
        clusters.merge(ds, [1, 4])
    with pytest.raises(ValueError):
        clusters.add_manual(ds, mask, "odd ones")
    with pytest.raises(ValueError):
        clusters.add_manual(ds, np.zeros(6, bool), "nothing")
    clusters.remove_manual(ds, 4)
    assert np.array_equal(clusters.view(ds).labels, raw)
    with pytest.raises(ValueError):
        clusters.remove_manual(ds, 1)


def test_edits_are_saved_with_the_dataset_and_a_new_run_starts_clean(tmp_path):
    ds = dataset([1, 1, 2, 3, 3, 3, 4, 4])
    clusters.merge(ds, [1, 2])
    clusters.rename(ds, {1: "CD4+", 3: "CD8+"})
    clusters.set_note(ds, 3, "after the merge")
    clusters.add_manual(ds, np.arange(8) == 7, "one event")
    before = clusters.view(ds)
    path = str(tmp_path / "d.npz")
    ds.save(path)
    again = clusters.view(Dataset.load(path))
    assert (again.ids, again.names, again.notes, again.members) == \
        (before.ids, before.names, before.notes, before.members)
    assert np.array_equal(again.labels, before.labels)
    # what the algorithm writes when it runs again: labels and its own entry
    back = Dataset.load(path)
    back.derived["FlowSOM"] = np.array([2, 2, 1, 1, 1, 1, 1, 1], dtype=np.int32)
    back.runs["FlowSOM"] = {"channels": ["A", "B"]}
    fresh = clusters.view(back)
    assert fresh.ids == [1, 2] and fresh.names == ["1", "2"]
    assert np.array_equal(fresh.labels, back.derived["FlowSOM"])
    # and a cluster drawn after it is not numbered into the old one's place
    c = clusters.add_manual(back, np.arange(8) == 0, "new")
    assert (clusters.view(back).labels == c).sum() == 1


def test_tables_carry_names_and_merged_frequencies_add_up():
    labels = np.array([1, 2, 2, 3, 1, 3, 3, 2, 1, 1])
    ds = dataset(labels, samples=2)
    s = Session()
    s.dataset = ds
    plain = s.cluster_tables()[1]
    assert plain.rows == ["1", "2", "3"] and plain.ids == [1, 2, 3]
    clusters.merge(ds, [2, 3])
    clusters.rename(ds, {1: "CD4+"})
    expression, freq = s.cluster_tables()
    assert freq.rows == expression.rows == ["CD4+", "2+3"] and freq.ids == [1, 2]
    assert freq.row_of(2) == 1 and freq.cluster_at(0) == 1 and freq.row_of(None) is None
    assert np.allclose(freq.values[1], plain.values[1] + plain.values[2])
    X = s.features(["A", "B"])
    assert np.allclose(expression.values[1], np.median(X[labels != 1], axis=0))
    assert s.clusters().names == ["CD4+", "2+3"]


# --- names from expression -------------------------------------------------------

def test_suggested_names_on_known_profiles():
    #                 CD4  CD8  Foxp3  CD45 (on everything: separates nothing)
    table = Table(["1", "2", "3", "4"], ["CD4 PE-A", "CD8 FITC-A", "Foxp3 APC-A", "CD45"],
                  np.array([[700, 150, 160, 800],
                            [690, 140, 600, 790],
                            [120, 720, 150, 810],
                            [130, 130, 170, 805.0]]))
    assert np.array_equal(clusters.calls(table, 100.0),
                          [[1, -1, -1, 0], [1, -1, 1, 0], [-1, 1, -1, 0], [-1, -1, -1, 0]])
    assert clusters.suggest_names(table, 100.0) == \
        ["CD4+", "CD4+ Foxp3+", "CD8+", "CD4- CD8- Foxp3-"]


def test_suggested_names_cut_at_the_widest_gap_and_number_twins():
    # three levels: the cut falls in the wider of the two gaps
    table = Table(["1", "2", "3"], ["CD4"], np.array([[100.0], [250.0], [800.0]]))
    assert clusters.calls(table, 100.0).ravel().tolist() == [-1, -1, 1]
    assert clusters.suggest_names(table, 100.0) == ["CD4-", "CD4- (2)", "CD4+"]
    # the threshold is in the table's units
    assert clusters.suggest_names(table, 600.0) == ["", "", ""]


def test_no_marker_separates_no_name_is_proposed():
    flat = Table(["1", "2", "3"], ["CD4", "CD8"],
                 np.array([[400, 300], [430, 310], [460, 290.0]]))
    assert not clusters.calls(flat, 100.0).any()
    assert clusters.suggest_names(flat, 100.0) == ["", "", ""]
    one = Table(["1"], ["CD4"], np.array([[900.0]]))
    assert clusters.suggest_names(one, 100.0) == [""]
    s = Session()
    assert s.suggest_cluster_names() == {}


# --- a cluster as a gate ---------------------------------------------------------

def nearest(X, codes):
    """Each row's nearest code, by the whole distance matrix."""
    d = ((X[:, None, :] - codes[None, :, :]) ** 2).sum(axis=2)
    return d.argmin(axis=1)


def test_cluster_gate_keeps_the_events_of_its_nodes():
    rng = np.random.default_rng(3)
    raw = np.abs(rng.lognormal(6, 1.5, (500, 2))).astype(np.float32)
    codes = rng.uniform(0, 1, (9, 2)) * [1023.0, 6.0]   # in the map's units: channels, asinh
    scales = [Biex(), Arcsinh(150.0)]
    gate = ClusterGate(["Comp-A", "B"], scales, [1023.0, 1.0], codes, [0, 4, 7])
    column = {"Comp-A": raw[:, 0], "B": raw[:, 1]}.__getitem__
    X = np.column_stack([Biex()(raw[:, 0].astype(np.float64)) * 1023,
                         np.arcsinh(raw[:, 1].astype(np.float64) / 150)])
    expected = np.isin(nearest(X, codes), [0, 4, 7])
    assert 0 < expected.sum() < 500
    assert np.array_equal(gate.contains(column), expected)
    gate.inside = False
    assert np.array_equal(gate.contains(column), ~expected)


def test_cluster_gate_is_copied_and_saved_whole():
    rng = np.random.default_rng(4)
    codes = rng.normal(size=(6, 3)) * 1e3 / 7        # nothing a short decimal holds
    gate = ClusterGate(["Comp-A", "Comp-B", "C"], [Biex(), Arcsinh(5.0), Linear()],
                       [1023.0, 1.0, 1023.0], codes, [1, 5])
    tree = GateNode("")
    tree.add(GateNode("Cells", RectangleGate(["C"], [0.0], [None]))).add(GateNode("CD8+", gate))
    copy = copy_tree(tree)
    saved = tree_from_dict(json.loads(json.dumps(tree_to_dict(tree))))
    for other in (copy, saved):
        g = other.find("Cells/CD8+").gate
        assert isinstance(g, ClusterGate) and g is not gate and g.codes is not gate.codes
        assert np.array_equal(g.codes, codes)        # to the last bit
        assert g.channels == gate.channels and g.nodes == [1, 5] and g.factors == gate.factors
        assert [type(s) for s in g.scales] == [Biex, Arcsinh, Linear]
        assert g.scales[1].cofactor == 5.0
    copy.find("Cells/CD8+").gate.nodes.append(2)
    assert gate.nodes == [1, 5]


@pytest.fixture
def folder(tmp_path):
    """Four files of three populations and, below 70,000 on FSC, debris
    stained like nothing in particular."""
    rng = np.random.default_rng(0)
    centres = np.array([[300, 3000, 100], [3000, 200, 100], [200, 300, 4000.0]])
    for name, n in (("3Ctrl1", 700), ("3Ctrl2", 900), ("3Case1", 1200), ("3Case2", 1500)):
        pop = rng.integers(0, 3, n)
        fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (n, 3)))
        fsc = rng.normal(90000, 9000, n)
        ev = np.column_stack([fsc, fluo]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), ev, ["FSC-A", "FITC-A", "PE-A", "APC-A"],
                  markers=["", "CD8", "CD4", "Foxp3"],
                  keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
    return str(tmp_path)


def gated_session(folder, transform="biex") -> Session:
    s = Session()
    s.open_folder(folder)
    first = s.usable()[0]
    first.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [85000.0], [None])))
    s.gates_changed(first)
    s.population, s.n = "Cells", 300
    s.build()
    s.run_flowsom(s.default_channels(), transform=transform, n_clusters=3, grid=5)
    return s


def agreement(s: Session, path: str, cluster: int) -> int:
    """The population's events in every file, having checked that each
    event of the dataset is in it exactly when the dataset says so."""
    ds, now = s.dataset, s.clusters()
    total = 0
    for i, name in enumerate(ds.samples):
        entry = next(e for e in s.usable() if e.name == name)
        mask = s.load(entry).mask(path)
        rows = ds.sample_index == i
        assert np.array_equal(mask[ds.event_index[rows]], now.labels[rows] == cluster), name
        total += int(mask.sum())
    return total


@pytest.mark.parametrize("transform", ["biex", "arcsinh", "logicle"])
def test_population_agrees_with_the_dataset_on_its_events(folder, transform):
    s = gated_session(folder, transform)
    ds = s.dataset
    assert len(ds) == 4 * 300 and s.clusters().ids == [1, 2, 3]
    paths = {c: s.make_population(c) for c in (1, 2, 3)}
    assert paths[1] == "Cells/FlowSOM 1"       # under what the map was trained on
    assert all(p in s.populations() for p in paths.values())
    totals = {c: agreement(s, p, c) for c, p in paths.items()}
    # every cell goes to some node: the clusters share out the parent, in
    # all the events of the files and not only the 300 drawn from each
    cells = sum(s.load(e).count("Cells") for e in s.usable())
    assert cells > 4 * 300 and sum(totals.values()) == cells
    assert sum(s.load(e).count("") for e in s.usable()) > cells
    # the same events by a calculation of our own, on one whole file
    entry = s.usable()[0]
    sample = s.load(entry)
    cols = [sample.fcs.index(c) for c in ds.runs["FlowSOM"]["channels"]]
    values = sample.compensated[:, cols].astype(np.float32).astype(np.float64)
    X = (Biex()(values) * 1023.0 if transform == "biex"
         else Logicle()(values) * 1023.0 if transform == "logicle"
         else np.arcsinh(values / 150.0))
    mine = s.flowsom.node_cluster[nearest(X, s.flowsom.codes)] == 2
    assert np.array_equal(sample.mask(paths[2]), mine & sample.mask("Cells"))


def test_population_of_a_merged_and_renamed_cluster_and_its_project(folder):
    s = gated_session(folder)
    s.set_scale(["FITC-A"], Arcsinh(150.0))     # the user's own scale on one marker
    s.run_flowsom(s.default_channels(), n_clusters=3, grid=5)
    ds = s.dataset
    merged = clusters.merge(ds, [2, 3])
    clusters.rename(ds, {merged: "CD4+ / Treg"})
    path = s.make_population(merged)
    assert path == "Cells/CD4+ - Treg"          # a slash would be a level of the tree
    total = agreement(s, path, merged)
    # made again after a change, it is the same population with the new gate
    clusters.undo_merges(ds)
    assert s.make_population(2, name="CD4+ - Treg") == path
    assert [c.name for c in s.usable()[0].gates.find("Cells").children] == ["CD4+ - Treg"]
    assert agreement(s, path, 2) < total
    # another session, opening the folder: the map came with the project
    with open(os.path.join(folder, "MIL-C.json"), encoding="utf-8") as f:
        assert '"type": "cluster"' in f.read()
    again = Session()
    again.open_folder(folder)
    for e, other in zip(s.usable(), again.usable()):
        assert np.array_equal(again.load(other).mask(path), s.load(e).mask(path))
    table = again.statistics(populations=[path])
    assert [row[table.columns.index("Count")] for row in table.rows] == \
        [int(s.load(e).mask(path).sum()) for e in s.usable()]


def test_what_cannot_become_a_population_says_why(folder, tmp_path):
    s = gated_session(folder)
    ds = s.dataset
    drawn = clusters.add_manual(ds, np.arange(len(ds)) < 50, "by hand")
    with pytest.raises(ValueError, match="only in this dataset"):
        s.make_population(drawn)
    with pytest.raises(ValueError, match="no such cluster"):
        s.make_population(99)
    path = str(tmp_path / "d.npz")
    ds.save(path)
    s.dataset, s.flowsom = Dataset.load(path), None
    with pytest.raises(ValueError, match="run FlowSOM again"):
        s.make_population(1)
    assert "Cells/FlowSOM 1" not in s.populations()


# --- the window --------------------------------------------------------------------

@pytest.fixture(scope="module")
def app(tmp_path_factory):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    QtWidgets = pytest.importorskip("PyQt6.QtWidgets")
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


@pytest.fixture
def window(app, folder):
    from milc.ui import worker
    from milc.ui.shell import MainShell

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(folder)
    w.samples.concatenate()
    settle()
    w.explore.clusters.setValue(3)
    w.explore.grid.setValue(5)
    w.explore.run(["UMAP", "FlowSOM"])
    settle()
    return w


def legend(w) -> list[str]:
    return [w.explore.legend.item(i).text().split("    ")[0]
            for i in range(w.explore.legend.count())]


def test_names_reach_every_view(window):
    w, c, s = window, window.clusters, window.session
    tools = c.tools
    assert not tools.buttons["Rename"].isEnabled()      # nothing selected yet
    assert tools.buttons["Suggest names"].isEnabled()
    c.select(2)
    assert tools.buttons["Rename"].isEnabled() and tools.buttons["Make population"].isEnabled()
    assert not tools.buttons["Undo merges"].isEnabled()
    assert not tools.buttons["Remove"].isEnabled()
    tools.rename("Tregs")
    tools.set_note("Foxp3 high")
    assert c.selected == 2
    assert c.heatmap.table.rows == c.freq.freq.rows == ["1", "Tregs", "3"]
    assert c.heatmap.cluster_at(0.5, 1.5) == 2 and c.freq.cluster_at(1.0) == 2
    assert c.map.content().color.names == ["1", "Tregs", "3"]
    assert c.tree.tree["names"] == {1: "1", 2: "Tregs", 3: "3"}
    assert "Cluster Tregs" in c.summary.text() and "Foxp3 high" in c.summary.text()
    assert legend(w) == ["1", "Tregs", "3"]
    assert [t.rows for t in c.prism_tables()] == [["1", "Tregs", "3"]]
    w.explore.split_by.setCurrentIndex(w.explore.split_by.findData("cluster|FlowSOM"))
    assert w.explore.view.content().split.names == ["1", "Tregs", "3"]
    c.test.setCurrentIndex(c.test.findData("rank"))
    assert {row[0] for row in c.results_table()[0].rows} == {"1", "Tregs", "3"}
    from milc import figures
    assert [n for n, _ in figures.category_names(c.tree.content())] == ["1", "Tregs", "3"]
    assert figures.preview_png(c.tree.content(), figures.FigureStyle())[:4] == b"\x89PNG"
    # a name in use is refused, and said, not raised at the user
    said = []
    c.sigStatus.connect(said.append)
    c.select(1)
    tools._try(lambda: tools.rename("Tregs"))
    assert "already" in said[-1] and s.clusters().names == ["1", "Tregs", "3"]


def test_merge_undo_and_suggested_names_in_the_window(window):
    w, c, s = window, window.clusters, window.session
    ds = s.dataset
    found = c.tools.suggestions()
    assert sorted(found.values()) == ["CD4+", "CD8+", "Foxp3+"]
    assert s.clusters().names == ["1", "2", "3"]         # proposed, not applied
    c.tools.apply_names(found)
    assert set(c.heatmap.table.rows) == {"CD4+", "CD8+", "Foxp3+"}
    raw = ds.derived["FlowSOM"].copy()
    c.select(3)
    assert c.tools.merge([2]) == 2 and c.selected == 2
    now = s.clusters()
    assert now.ids == [1, 2] and len(c.heatmap.table.rows) == 2 and len(legend(w)) == 2
    assert int(c.selection_mask().sum()) == int(np.isin(raw, [2, 3]).sum())
    assert set(int(k) for k in c.tree.tree["cluster"][c.tree.tree["sizes"] > 0]) == {1, 2}
    assert c.tools.buttons["Undo merges"].isEnabled()
    assert c.tools.buttons["Merge with"].isEnabled()     # there is still cluster 1
    c.tools.undo_merges()
    assert s.clusters().ids == [1, 2, 3] and np.array_equal(ds.derived["FlowSOM"], raw)
    assert c.selected is None and len(legend(w)) == 3


def test_lasso_becomes_a_cluster_but_not_a_population(window):
    w, c, s = window, window.clusters, window.session
    ds = s.dataset
    assert not c.tools.buttons["Save selection"].isEnabled()
    with pytest.raises(ValueError):
        c.tools.save_selection("nothing outlined")
    mask = np.zeros(len(ds), dtype=bool)
    mask[:40] = True
    c._lassoed(mask)
    assert c.tools.buttons["Save selection"].isEnabled()
    assert "cannot become a population" in c.tools.buttons["Save selection"].toolTip()
    drawn = c.tools.save_selection("outlined")
    assert c.selected == drawn and int(c.selection_mask().sum()) == 40
    assert c.heatmap.table.rows[-1] == "outlined" and legend(w)[-1] == "outlined"
    assert c.tools.buttons["Remove"].isEnabled()
    assert not c.tools.buttons["Make population"].isEnabled()
    with pytest.raises(ValueError):
        c.tools.make_population()
    c.tools.remove()
    assert s.clusters().ids == [1, 2, 3] and c.selected is None


def test_cluster_becomes_a_population_the_gates_tab_lists(window):
    w, c, s = window, window.clusters, window.session
    c.select(1)
    path = c.tools.make_population()
    assert path == "FlowSOM 1"                 # no gate was the dataset's: at the root
    assert agreement(s, path, 1) > 0
    w.gates.samples.setCurrentRow(0)           # the tab was told: it lists what is there now
    names = []
    it = w.gates.tree.invisibleRootItem()
    stack = [it.child(i) for i in range(it.childCount())]
    while stack:
        item = stack.pop()
        names.append(item.text(0))
        stack.extend(item.child(i) for i in range(item.childCount()))
    assert "FlowSOM 1" in names
    assert path in s.populations()
    w.gates.select(path)                       # a gate with no outline: nothing to draw, no error
    assert w.gates.node.name == "FlowSOM 1"

