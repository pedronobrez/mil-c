"""The window, driven without a screen: open a folder, concatenate, run,
look. What is asserted is what the session and the widgets hold, never
pixels."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

QtWidgets = pytest.importorskip("PyQt6.QtWidgets")

from milc.fcs import write_fcs  # noqa: E402


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import style
    from milc.ui.theme import apply_defaults
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    style.apply(app)
    apply_defaults()
    return app


@pytest.fixture
def folder(tmp_path):
    rng = np.random.default_rng(0)
    centres = np.array([[300, 3000, 100], [3000, 200, 100], [200, 300, 4000.0]])
    for name, n in (("3Ctrl1", 240), ("3Ctrl2", 400), ("3Case1", 900), ("3Case2", 1500)):
        pop = rng.integers(0, 3, n)
        fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (n, 3)))
        ev = np.column_stack([rng.normal(90000, 9000, n), fluo]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), ev, ["FSC-A", "FITC-A", "PE-A", "APC-A"],
                  markers=["", "CD8", "CD4", "Foxp3"],
                  keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
    write_fcs(str(tmp_path / "concat_other.fcs"), np.zeros((5, 2), np.float32), ["A", "B"])
    return str(tmp_path)


def test_folder_to_map(app, folder):
    from milc.ui import worker
    from milc.ui.shell import MainShell

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(folder)
    s = w.session
    assert [e.name for e in s.entries if e.problem] == ["concat_other"]
    assert s.meta_keys() == ["day", "group", "animal"]
    plan = s.plan()
    assert plan.limited_by == "3Ctrl1" and plan.total == 4 * 240
    assert "Limited by" in w.samples.summary.text()

    # unticking the smallest sample lifts the limit to the next one
    w.samples.table.item(2, 0).setCheckState(w.samples.table.item(2, 0).checkState().Unchecked)
    assert s.entries[2].name == "3Ctrl1" and s.plan().limited_by == "3Ctrl2"
    w.samples.table.item(2, 0).setCheckState(w.samples.table.item(2, 0).checkState().Checked)

    w.samples.concatenate()
    settle()
    ds = s.dataset
    assert len(ds) == 960 and w.tabs.currentWidget() is w.explore
    e = w.explore
    assert e.chosen_markers() == ["FITC-A", "PE-A", "APC-A"]
    assert e.x_axis.currentData() == "channel|FSC-A"  # something to look at already

    e.clusters.setValue(3)
    e.grid.setValue(5)
    e.run(["UMAP", "FlowSOM"])
    settle()
    assert {"UMAP1", "UMAP2", "FlowSOM"} <= set(ds.derived)
    assert e.x_axis.currentData() == "derived|UMAP1" and e.y_axis.currentData() == "derived|UMAP2"
    assert e.color_by.currentData() == "cluster|FlowSOM"
    assert e.legend.count() == 3 and len(e.view.panels) == 1

    e.split_by.setCurrentIndex(e.split_by.findData("meta|group"))
    assert len(e.view.panels) == 2
    e.color_by.setCurrentIndex(e.color_by.findData("expression|PE-A"))
    assert "CD4" in e.legend.item(0).text()
    expression, freq = s.cluster_tables()
    assert expression.values.shape == (3, 3) and freq.columns == ds.samples

    # the clusters tab: one click marks a cluster everywhere
    c = w.clusters
    assert c.stack.currentIndex() == 1 and c.heatmap.table.values.shape == (3, 3)
    assert c.group_by.currentData() == "group" and c.embedding.currentData() == "UMAP"
    c.heatmap.sigCluster.emit(2)
    assert c.selected == 2 and "Cluster 2" in c.summary.text()
    assert int(c.selection_mask().sum()) == int((ds.derived["FlowSOM"] == 2).sum())
    assert "Case" in c.summary.text() and "Ctrl" in c.summary.text()
    tree = c.tree.tree
    node = int(np.flatnonzero(tree["sizes"] > 0)[0])
    assert c.tree.cluster_at(*tree["pos"][node]) == tree["cluster"][node]
    assert c.tree.cluster_at(1e6, 1e6) is None
    c.heatmap.sigCluster.emit(2)  # again: let go
    assert c.selected is None

    # a lasso asks the same of any events
    x, y = ds.derived["UMAP1"], ds.derived["UMAP2"]
    box = np.array([[x.min() - 1, y.min() - 1], [np.median(x), y.min() - 1],
                    [np.median(x), y.max() + 1], [x.min() - 1, y.max() + 1]])
    c.map.lasso(box)
    assert c.selection_mask().sum() == (x < np.median(x)).sum()
    assert "Selection" in c.summary.text() and "clusters" in c.summary.text()
    c.clear_selection()
    assert c.selection_mask() is None

    # any plot as a figure; the look is remembered, the words are not
    from milc.ui.figure_dialog import FigureDialog, last_style, load_presets
    d = FigureDialog(e.view.content(), w, "map")
    assert not d.categories.isVisibleTo(d)  # coloured by expression: a ramp, no categories
    assert d.colormap.isVisibleTo(d)
    d.close()
    e.color_by.setCurrentIndex(e.color_by.findData("meta|group"))
    d = FigureDialog(e.view.content(), w, "map")
    assert d.categories.rowCount() == 2
    d.width.setValue(85)
    d.title.setText("Mine")
    d.title.textEdited.emit("Mine")
    d.set_colour("Case", "#112233")
    d.categories.item(1, 1).setText("Control")
    assert d.style.labels == {"Ctrl": "Control"}
    path = os.path.join(folder, "out.svg")
    d.write(path)
    svg = open(path).read()
    assert "Mine" in svg and "Control" in svg and "#112233" in svg
    d.store_preset("Paper")
    d.accept()
    assert load_presets()["Paper"]["width_mm"] == 85 and "title" not in load_presets()["Paper"]
    assert last_style().colors == {"Case": "#112233"}
    d = FigureDialog(c.heatmap.content(), w, "heat")
    assert d.style.width_mm == 85 and d.style.title == ""  # look kept, title not
    assert d.annotate.isVisibleTo(d) and not d.opacity.isVisibleTo(d)
    d.apply_look({})
    assert d.style.width_mm == 120
    d.close()
    for view in (c.freq, c.tree, c.map):
        assert view.content() is not None

    # a new concatenation takes the clusters away again
    w.samples.concatenate()
    settle()
    assert c.stack.currentIndex() == 0
    w.close()
