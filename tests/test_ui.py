"""The window, driven without a screen: open a folder, concatenate, run,
look. What is asserted is what the session and the widgets hold, never
pixels."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

QtWidgets = pytest.importorskip("PyQt6.QtWidgets")
from PyQt6 import QtCore  # noqa: E402

from milc.fcs import write_fcs  # noqa: E402


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
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
    assert s.meta_keys() == ["Group"] and s.entries[0].meta == {"Group": "Case"}
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

    e.split_by.setCurrentIndex(e.split_by.findData("meta|Group"))
    assert len(e.view.panels) == 2
    e.color_by.setCurrentIndex(e.color_by.findData("expression|PE-A"))
    assert "CD4" in e.legend.item(0).text()
    expression, freq = s.cluster_tables()
    assert expression.values.shape == (3, 3) and freq.columns == ds.samples

    # the clusters tab: one click marks a cluster everywhere
    c = w.clusters
    assert c.stack.currentIndex() == 1 and c.heatmap.table.values.shape == (3, 3)
    assert c.group_by.currentData() == "Group" and c.embedding.currentData() == "UMAP"
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
    e.color_by.setCurrentIndex(e.color_by.findData("meta|Group"))
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


def test_gates_by_hand(app, folder):
    """Open a sample, look at it, gate it: a rectangle, a polygon on a
    biexponential axis, a range on a histogram — and every sample gets
    them, and they are still there when the folder is opened again."""
    from milc.session import PROJECT_FILE as GATES_FILE
    from milc.transforms import Biex
    from milc.ui import worker
    from milc.ui.gates_tab import HISTOGRAM
    from milc.ui.shell import MainShell

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(folder)
    w.show()
    g = w.gates
    w.tabs.setCurrentWidget(g)
    settle()
    assert g.entry.name == "3Case1" and g.node.path == "" and g.tree.topLevelItemCount() == 1
    s = g.sample
    assert g.x_axis.findData("Comp-FITC-A") >= 0 and g.x_axis.findData("FSC-A") >= 0

    # a rectangle on scatter, in the file's own units
    g.set_axes("FSC-A", "Comp-PE-A")
    cells = g.add_rectangle(80000, 0.0, 100000, 1.0, "Cells")
    fsc = s.column("FSC-A")
    assert s.count("Cells") == int(((fsc >= 80000) & (fsc < 100000)).sum())
    assert cells.gate.mins[1] < 0 and cells.gate.maxs[1] == pytest.approx(262144, rel=1e-2)

    # inside it, a polygon drawn on two biexponential axes: tested where drawn
    g.select("Cells")
    g.set_axes("Comp-PE-A", "Comp-FITC-A")
    assert g.node.path == "Cells"
    box = np.array([[0.5, 0.0], [1.0, 0.0], [1.0, 0.45], [0.5, 0.45]])
    cd4 = g.add_polygon(box, "CD4")
    t = Biex()
    pe, fitc = t(s.column("Comp-PE-A")), t(s.column("Comp-FITC-A"))
    expected = (s.mask("Cells") & (pe > 0.5) & (pe < 1.0) & (fitc > 0.0) & (fitc < 0.45)).sum()
    assert s.count("Cells/CD4") == expected > 20
    assert cd4.gate.scales is not None and len(g._overlays) == 1
    assert g.x_axis.currentData() == "Comp-PE-A"   # the view did not jump away

    # a range on a histogram is a gate of one parameter
    g.set_axes("Comp-APC-A", HISTOGRAM)
    foxp3 = g.add_rectangle(0.6, 0, 1.0, 0, "Foxp3+")
    assert foxp3.gate.dims == ["Comp-APC-A"] and 0 < s.count("Cells/Foxp3+") < s.count("Cells")
    assert g.add_rectangle(0.6, 0, 1.0, 0, "Foxp3+").name == "Foxp3+ 2"   # names stay apart

    # every sample has them, the Samples tab can now ask for them, and they are on disk
    assert all(e.edited and e.gates.find("Cells/CD4") for e in w.session.usable())
    assert "Cells/CD4" in w.session.populations()
    assert w.samples.population.findData("Cells/CD4") >= 0
    assert os.path.exists(os.path.join(folder, GATES_FILE))

    # dragging an outline moves the gate: the count follows
    g.select("Cells")
    g.set_axes("Comp-PE-A", "Comp-FITC-A")
    before = s.count("Cells/CD4")
    child, roi, label = next(o for o in g._overlays if o[0].name == "CD4")
    roi.setPos(roi.pos().x() + 0.3, roi.pos().y())
    g._gate_moved(child, roi)
    assert g.sample.count("Cells/CD4") < before and label.toPlainText().startswith("CD4")

    # renaming and deleting, with what is gated inside
    g.select("Cells/Foxp3+ 2")
    g.delete_selected()
    assert g.node.path == "Cells" and "Cells/Foxp3+ 2" not in w.session.populations()
    g.select("Cells/CD4")
    g.rename("T helper")
    assert "Cells/T helper" in w.session.populations()

    # only this sample, when asked
    g.everyone.setChecked(False)
    g.select("Cells")
    g.add_rectangle(0.0, 0.0, 0.2, 0.2, "Just here")
    assert "Cells/Just here" not in w.session.populations()   # not common to all
    g.everyone.setChecked(True)
    w.close()

    # opened again, the gates are the ones that were drawn
    from milc.session import Session
    again = Session()
    again.open_folder(folder)
    assert again.populations() == ["Cells", "Cells/T helper", "Cells/Foxp3+"]
    again.population = "Cells/T helper"
    first = again.usable()[0]
    assert again.count(first) == again.load(first).count("Cells/T helper")
    assert again.load(first).count("Cells") > 0 and first.edited
    again.discard_gates()
    assert again.populations() == [] and not any(e.edited for e in again.usable())
    third = Session()
    third.open_folder(folder)
    assert third.populations() == []   # and they stay gone


def test_ramp_and_appearance(app, folder):
    from milc import figures
    from milc.ui import palette, style
    from milc.ui.shell import MainShell
    w = MainShell()
    w.clusters.ramp.set_choice("magma")
    w.clusters.ramp._publish()
    assert palette.ramp_name() == "magma" and w.explore.ramp.name() == "magma"
    assert palette.ramp(2)[0].tolist() == figures.lookup_table("magma", steps=2)[0].tolist()
    w.explore.ramp.choose_custom(["#000000", "#ff0000"])
    assert palette.ramp_name() == "custom" and palette.ramp(2)[-1].tolist() == [255, 0, 0, 255]
    assert w.clusters.ramp.name() == "custom"
    # a figure starts from what is on the screen
    from milc.ui.figure_dialog import FigureDialog
    d = FigureDialog(figures.HeatmapContent(
        __import__("milc.cluster_stats", fromlist=["Table"]).Table(
            ["1", "2"], ["a", "b"], np.array([[0.0, 1], [2, 3]]))), w)
    assert d.style.colormap == "custom" and d.style.colormap_colors == ["#000000", "#ff0000"]
    assert figures.preview_png(d.content, d.style, dpi=30)[:4] == b"\x89PNG"
    d.close()
    w.set_appearance("dark")
    assert style.is_dark() and w._appearance["dark"].isChecked()
    w.set_appearance("light")
    assert not style.is_dark()
    w.set_appearance("")
    palette.set_ramp("viridis")
    w.close()


def test_columns_are_the_users(app, folder):
    """The table is typed into: a column of one's own, a value for several
    rows at once, columns out of the file names — and two of them together
    make the groups, in the plan, on the map and in the frequency chart."""
    from milc.session import PROJECT_FILE, Session
    from milc.ui import worker
    from milc.ui.samples_tab import PatternDialog
    from milc.ui.shell import MainShell

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(folder)
    t, s = w.samples, w.session
    cols = t._columns()
    assert cols[:3] == ["", "Sample", "Group"] and not os.path.exists(os.path.join(folder, PROJECT_FILE))
    group = cols.index("Group")
    flags = t.table.item(0, group).flags()
    assert flags & flags.ItemIsEditable and not (t.table.item(0, 1).flags() & flags.ItemIsEditable)

    # typing in one cell
    t.table.item(0, group).setText(" Treated ")
    assert s.entries[0].meta["Group"] == "Treated" and t.table.item(0, group).text() == "Treated"
    # typing with several rows selected fills them all
    t.add_column("Day")
    day = t._columns().index("Day")
    t.table.clearSelection()
    for r in (0, 1, 2):
        t.table.selectRow(r) if r == 0 else t.table.selectionModel().select(
            t.table.model().index(r, day),
            t.table.selectionModel().SelectionFlag.Select | t.table.selectionModel().SelectionFlag.Rows)
    t.table.item(1, day).setText("3")
    assert [e.meta["Day"] for e in s.entries[:4]] == ["3", "3", "3", ""]
    t.add_column("Day")                      # a second one of the same name is refused
    assert s.columns == ["Group", "Day"]
    t.rename_column("Day", "Time point")
    assert s.columns == ["Group", "Time point"] and s.entries[0].meta["Time point"] == "3"
    t.remove_column("Time point")
    assert s.columns == ["Group"] and "Time point" not in s.entries[0].meta

    # columns out of the names, with a pattern the dialog proposes
    d = PatternDialog(s, t)
    assert d.pattern.text() == r"(?P<Number1>\d+)(?P<Group>[A-Za-z]+)(?P<Number2>\d+)"
    assert d.preview.item(0, 0).text() == "3Case1" and d.preview.item(0, 2).text() == "Case"
    d.pattern.setText(r"(?P<Day>\d+)(?P<Group>")
    assert not d.box.button(d.box.StandardButton.Ok).isEnabled()   # not a pattern yet
    d.close()
    assert t.fill_from_names(r"(?P<Day>\d+)(?P<Group>[A-Za-z]+)(?P<Animal>\d+)") == ["Day", "Group", "Animal"]
    assert s.columns == ["Group", "Day", "Animal"] and s.entries[0].meta["Group"] == "Case"
    s.set_value(s.entries[1:2], "Day", "7")
    s.set_value(s.entries[3:4], "Day", "7")

    # two columns together are the groups
    t.set_group_by("Day", True)
    assert s.group_by == ["Group", "Day"] and t.group_by.text() == "Group + Day"
    assert s.group_of(s.entries[1].meta) == "Case · 7"
    assert set(s.plan().group_totals()) == {"Case · 3", "Case · 7", "Ctrl · 3", "Ctrl · 7"}

    # it is all still there when the folder is opened again
    again = Session()
    again.open_folder(folder)
    assert again.columns == ["Group", "Day", "Animal"] and again.group_by == ["Group", "Day"]
    assert again.entries[1].meta == {"Group": "Case", "Day": "7", "Animal": "2"}

    # on the map: split by one column and then by the other
    t.concatenate()
    settle()
    e = w.explore
    e.grid.setValue(4)
    e.clusters.setValue(3)
    e.run(["FlowSOM"])
    settle()
    e.split_by.setCurrentIndex(e.split_by.findData("meta|Group"))
    assert len(e.view.panels) == 2
    e.split_by2.setCurrentIndex(e.split_by2.findData("meta|Day"))
    assert len(e.view.panels) == 4
    assert [n for n in e.view.content().split.names] == ["Case · 3", "Case · 7", "Ctrl · 3", "Ctrl · 7"]

    # in the clusters: the groups compared within each day
    c = w.clusters
    assert c.group_by.currentData() == "Group" and c.facet_by.currentData() == "Day"
    assert len(c.freq.panels) == 2 and c.freq.content().facets == ["3", "7", "3", "7"]
    c.freq.sigCluster.emit(1)
    assert "Case · 3" in c.summary.text() and "Ctrl · 7" in c.summary.text()
    from milc import figures
    assert figures.preview_png(c.freq.content(), figures.FigureStyle(font="DejaVu Sans"),
                               dpi=40)[:4] == b"\x89PNG"
    c.facet_by.setCurrentIndex(0)
    assert len(c.freq.panels) == 1 and c.freq.content().facets is None
    w.close()


def test_statistics_tab(app, folder, tmp_path):
    from milc.ui import worker
    from milc.ui.gates_tab import HISTOGRAM
    from milc.ui.shell import MainShell

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(folder)
    w.show()
    st = w.statistics
    assert st.populations.topLevelItemCount() == 1 and st.parameters.count() == 4
    assert st.chosen_statistics() == ["Median"] and not st.export_button.isEnabled()

    # a gate drawn in Gates is offered here at once
    w.tabs.setCurrentWidget(w.gates)
    settle()
    w.gates.set_axes("FSC-A", HISTOGRAM)
    w.gates.add_rectangle(85000, 0, 95000, 0, "Cells")
    assert st.chosen_populations() == ["", "Cells"]

    st.tick(parameters=["CD4", "Comp-FITC-A"], statistics=["Median", "Robust CV"])
    st.compute()
    settle()
    t = st.result
    assert len(t.rows) == 8 and t.columns[-4:] == ["CD8 · Median", "CD8 · Robust CV",
                                                    "CD4 · Median", "CD4 · Robust CV"]
    assert st.table.rowCount() == 8 and st.table.columnCount() == len(t.columns) - 1  # no Path
    assert "4 samples × 2 populations" in st.caption.text() and st.export_button.isEnabled()

    # sorting is by the number, not by its text
    count = [st.table.horizontalHeaderItem(i).text() for i in range(st.table.columnCount())].index("Count")
    st.table.sortItems(count, QtCore.Qt.SortOrder.DescendingOrder)
    shown = [st.table.item(r, count).value for r in range(8)]
    assert shown == sorted(shown, reverse=True) and st.table.item(0, count).text() == "1,500"

    # leaving a sample out, and knowing when the table is behind
    w.session.entries[0].include = False
    st.compute()
    settle()
    assert len(st.result.rows) == 6 and st.stale.text() == ""
    w.samples.add_column("Day")
    assert "again" in st.stale.text()
    st.tick(populations=["Cells"])
    st.compute()
    settle()
    assert len(st.result.rows) == 3 and st.result.columns[:3] == ["Sample", "Group", "Day"]
    st.copy()
    assert QtWidgets.QApplication.clipboard().text().splitlines()[0].startswith("Sample\tGroup\tDay")
    w.close()
