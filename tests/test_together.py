"""Where the pieces meet: a second clustering in the Clusters tab, named
clusters in differential abundance, a workspace's unverified gates."""

import numpy as np

from milc import clusters
from test_ui import app, folder  # noqa: F401  (fixtures)


def _analysed(app, folder):
    from milc.ui import worker
    from milc.ui.shell import MainShell
    w = MainShell()
    w.open_folder(folder)
    def settle():
        for _ in range(4):
            worker.wait_all()
            app.processEvents()

    w.samples.concatenate()
    settle()
    for what in (["FlowSOM"], ["PhenoGraph"]):
        w.explore.run(what)
        settle()
    return w


def test_the_clusters_tab_shows_either_clustering(app, folder):
    w = _analysed(app, folder)
    c, ds = w.clusters, w.session.dataset
    assert [c.clustering.itemData(i) for i in range(c.clustering.count())] == [
        "FlowSOM", "PhenoGraph"]
    assert c.key == "FlowSOM" and c._tree is not None
    flowsom_rows = list(c._frequency.rows)
    c.select(w.session.clusters("FlowSOM").ids[0])

    c.clustering.setCurrentIndex(1)
    assert c.key == "PhenoGraph" and c.selected is None
    found = w.session.clusters("PhenoGraph")
    assert list(c._frequency.rows) == found.names
    assert c._frequency.values.shape[0] == len(np.unique(ds.derived["PhenoGraph"]))
    # the tree is FlowSOM's, and is not drawn for another clustering
    assert c._tree is None and not c.tree_color.isEnabled()
    # the same tools work on it; FlowSOM's own clusters are untouched
    c.select(found.ids[0])
    c.tools.rename("Big")
    assert w.session.clusters("PhenoGraph").names[0] == "Big"
    assert "Big" in c._frequency.rows
    c.set_clustering("FlowSOM")
    assert list(c._frequency.rows) == flowsom_rows and c._tree is not None
    w.close()


def test_differential_abundance_counts_the_clusters_as_named_and_merged(app, folder):
    w = _analysed(app, folder)
    c, ds = w.clusters, w.session.dataset
    ids = w.session.clusters("FlowSOM").ids
    clusters.rename(ds, {ids[0]: "Tregs"})
    if len(ids) > 2:
        clusters.merge(ds, [ids[1], ids[2]])
    c.clusters_changed()
    now = w.session.clusters("FlowSOM")
    d = c.abundance_dialog()
    assert list(d.counts.rows) == now.names and "Tregs" in d.counts.rows
    assert d.counts.values.sum() == len(ds)
    for name, row in zip(now.names, d.counts.values):
        assert row.sum() == (now.labels == now.ids[now.names.index(name)]).sum()
    w.close()


def test_a_workspace_says_what_it_left_unread(tmp_path):
    from milc.session import Session
    from milc.wsp import Skipped, UNVERIFIED
    s = Session()
    assert s.workspace_report() == ""

    class Fake:
        skipped = [Skipped("a.fcs", "Cells/Odd", "CurlyQuad", "not read"),
                   Skipped("b.fcs", "Cells/Odd", "CurlyQuad", "not read")]
        samples = []
    s.workspace = Fake()
    assert s.workspace_report() == "2 population(s) not read (2 CurlyQuad)"
    assert UNVERIFIED
