"""Several panels on a page, and the report: every kind of content as a
panel, the page's size, what the files say, and the window that collects
the panels. What is asserted is text, artists and sizes, never pixels."""

import os
import re

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from milc import cluster_stats as cs  # noqa: E402
from milc import density, figures, plate  # noqa: E402
from milc.dataset import Dataset  # noqa: E402
from milc.fcs import write_fcs  # noqa: E402

FONT = "DejaVu Sans"


def _contents() -> dict:
    from milc import flowsom as fs
    rng = np.random.default_rng(0)
    centres = np.array([[0.2, 0.8, 0.2, 0.5], [0.8, 0.2, 0.2, 0.5], [0.2, 0.2, 0.9, 0.5]])
    X = np.vstack([c + rng.normal(0, 0.04, (100, 4)) for c in centres])
    r = fs.flowsom(X, 3, 4, 4, seed=0)
    sidx = np.arange(len(X)) % 4
    names = ["CD8 FITC", "CD4 PE", "Foxp3", "IFNg"]
    pos, edges = fs.tree(r.codes, r.node_sizes())
    tree = {"pos": pos, "edges": edges, "sizes": r.node_sizes(),
            "cluster": r.node_cluster, "medians": r.codes, "markers": names}
    group = np.where(sidx % 2 == 0, "Case", "Ctrl")
    series = [figures.Series(name, g, *density.histogram(X[sidx == k, 0], 0, 1))
              for k, (name, g) in enumerate((("s1", "Case"), ("s2", "Ctrl")))]
    return {
        "map": figures.MapContent(X[:, 0], X[:, 1], "UMAP1", "UMAP2",
                                  color=figures.categories_of(group, "Group")),
        "values": figures.MapContent(X[:, 0], X[:, 1], "tSNE1", "tSNE2",
                                     color=figures.Values(X[:, 2], "Foxp3"),
                                     split=figures.categories_of(group, "Group")),
        "heat": figures.HeatmapContent(cs.expression(X, r.labels, names)),
        "freq": figures.FrequencyContent(
            cs.frequencies(r.labels, sidx, ["a", "b", "c", "d"]), ["7", "7", "14", "14"]),
        "tree": figures.TreeContent(tree),
        "tree_marker": figures.TreeContent(tree, "Foxp3"),
        "overlay": figures.OverlayContent(series, "CD8", "% of the tallest bin",
                                          xrange=(0, 1)),
    }


def _drawn_texts(fig) -> list[list[str]]:
    """Per panel, every string matplotlib will write in it."""
    from matplotlib.text import Text
    return [[t.get_text() for t in sub.findobj(Text) if t.get_text()]
            for sub in fig.subfigs]


def _pdf_pages(path: str) -> int:
    try:
        from PyQt6 import QtPdf
    except ImportError:
        with open(path, "rb") as f:
            return len(re.findall(rb"/Type\s*/Page\b", f.read()))
    doc = QtPdf.QPdfDocument(None)
    doc.load(path)
    return doc.pageCount()


def test_each_content_alone_and_mixed(tmp_path):
    contents = _contents()
    style = figures.FigureStyle(font=FONT)
    for name, content in contents.items():
        for legend in ("right", "bottom", "inside", "none"):
            p = plate.Plate(columns=1, width_mm=85, style=style)
            p.add(content, figures.FigureStyle(font=FONT, legend=legend), title=name)
            assert plate.preview_png(p, dpi=40)[:4] == b"\x89PNG", (name, legend)

    p = plate.Plate(columns=3, width_mm=180, style=style)
    for name, content in contents.items():
        p.add(content, title=f"The {name}", name=name)
    assert p.rows() == 3 and p.size_mm() == (180, 135)
    with plate._drawn(p) as fig:
        subs = dict(zip(contents, fig.subfigs))
        texts = dict(zip(contents, _drawn_texts(fig)))
        # each panel has its own legend or colour bar, inside its own cell
        assert [t.get_text() for t in subs["map"].legends[0].get_texts()] == ["Case", "Ctrl"]
        assert subs["map"].legends[0].get_title().get_text() == "Group"
        assert [t.get_text() for t in subs["freq"].legends[0].get_texts()] == ["7", "14"]
        assert [t.get_text() for t in subs["overlay"].legends[0].get_texts()] == ["Case", "Ctrl"]
        assert len(subs["tree"].legends[0].get_texts()) == 3
        bars = {n: [a for a in s.axes if a.get_label() == "<colorbar>"]
                for n, s in subs.items()}
        assert {n for n, b in bars.items() if b} == {"values", "heat", "tree_marker"}
        assert bars["values"][0].get_ylabel() == "Foxp3"
        assert bars["heat"][0].get_ylabel() == "Per marker"
        fig.canvas.draw()
        for name, sub in subs.items():
            cell = sub.bbox
            for artist in sub.legends + sub.axes:
                box = artist.get_window_extent()
                assert (box.x0 >= cell.x0 - 1 and box.x1 <= cell.x1 + 1
                        and box.y0 >= cell.y0 - 1 and box.y1 <= cell.y1 + 1), name
            # the letter's band is above the plot, not on it
            top = max(a.get_position().y1 for a in sub.axes)
            assert sub._suptitle.get_window_extent().y0 >= cell.y0 + top * cell.height - 1
    for k, (name, said) in enumerate(texts.items()):
        assert plate.letter(k) in said and f"The {name}" in said
    assert "UMAP1" in texts["map"] and "tSNE2" in texts["values"]
    assert "CD4" in texts["heat"] and "CD8" in texts["overlay"]


def test_files_size_and_text(tmp_path):
    from PIL import Image
    contents = _contents()
    p = plate.Plate(columns=2, width_mm=120, height_mm=95,
                    style=figures.FigureStyle(font=FONT, dpi=200, colors={"Case": "#c0392b"}))
    p.add(contents["map"], title="Map of it")
    p.add(contents["heat"], figures.FigureStyle(font=FONT, legend_title="Scaled"))
    p.add(contents["freq"], figures.FigureStyle(font=FONT, ylabel="Share"), title="Clusters")
    p.panels[2].colspan = 2
    assert p.cells() == [(0, 0, 1, 1), (0, 1, 1, 1), (1, 0, 1, 2)]
    for ext in ("svg", "pdf", "png", "tiff"):
        plate.save(p, str(tmp_path / f"p.{ext}"))
    for ext in ("png", "tiff"):
        im = Image.open(tmp_path / f"p.{ext}")
        assert abs(im.size[0] - 120 / 25.4 * 200) <= 1 and abs(im.size[1] - 95 / 25.4 * 200) <= 1
        assert round(im.info["dpi"][0]) == 200
    svg = (tmp_path / "p.svg").read_text(encoding="utf-8")
    assert 'width="340.1' in svg                       # 120 mm in points
    written = re.findall(r"<text[^>]*>([^<]*)</text>", svg)
    for word in ("A", "B", "C", "Map of it", "Clusters", "UMAP1", "Ctrl", "Group",
                 "Scaled", "min", "max", "Share", "14", "Foxp3"):
        assert word in written, word
    # one colour for a name across the page, from the page's style
    assert "#c0392b" in svg.lower()
    assert (tmp_path / "p.pdf").read_bytes().startswith(b"%PDF")
    assert _pdf_pages(str(tmp_path / "p.pdf")) == 1

    p.letters = False
    plate.save(p, str(tmp_path / "bare.svg"))
    bare = re.findall(r"<text[^>]*>([^<]*)</text>",
                      (tmp_path / "bare.svg").read_text(encoding="utf-8"))
    assert "Map of it" in bare and not {"A", "B", "C"} & set(bare)

    with pytest.raises(ValueError):
        plate.save(p, str(tmp_path / "p.bmp"))
    with pytest.raises(ValueError):
        plate.preview_png(plate.Plate())


def test_grid_order_and_styles():
    c = _contents()["map"]
    p = plate.Plate(columns=3, style=figures.FigureStyle(font=FONT, font_size=7))
    wide = p.add(c, name="wide")
    tall = p.add(c, name="tall")
    for name in ("a", "b", "c"):
        p.add(c, name=name)
    wide.colspan, tall.rowspan = 2, 2
    # in order, each in the first free place after the one before
    assert p.cells() == [(0, 0, 1, 2), (0, 2, 2, 1), (1, 0, 1, 1), (1, 1, 1, 1), (2, 0, 1, 1)]
    wide.colspan = 9                                    # no wider than the page
    assert p.cells()[0] == (0, 0, 1, 3) and p.rows() == 3
    assert p.move(4, -9) == 0 and [x.name for x in p.panels] == ["c", "wide", "tall", "a", "b"]
    assert p.move(0, 1) == 1 and p.move(4, 3) == 4
    assert [plate.letter(k) for k in (0, 25, 26, 27, 52)] == ["A", "Z", "AA", "AB", "BA"]

    # a panel keeps its own look and takes the page's type
    theirs = figures.FigureStyle(font="Courier", font_size=12, dpi=72, legend="bottom",
                                 title="T", colors={"Case": "#111111"}, width_mm=40)
    p = plate.Plate(style=figures.FigureStyle(font=FONT, font_size=7,
                                              colors={"Case": "#222222", "Ctrl": "#333333"}))
    panel = p.add(c, theirs)
    assert panel.title == "T" and panel.style == {"legend": "bottom",
                                                  "colors": {"Case": "#111111"}}
    s = p.style_of(panel)
    assert (s.font, s.font_size, s.dpi, s.legend, s.title) == (FONT, 7, 300, "bottom", "")
    assert s.colors == {"Case": "#111111", "Ctrl": "#333333"}
    assert p.style.colors["Case"] == "#222222"           # the page's is not written over


def test_refresh_from_the_sources():
    c = _contents()

    class View:
        def __init__(self, content):
            self.shown = content

        def content(self):
            return self.shown

    a, b = View(c["map"]), View(c["heat"])
    p = plate.Plate(style=figures.FigureStyle(font=FONT))
    p.add(a.content(), source=a)
    p.add(b.content(), source=b)
    p.add(c["tree"])                                    # from nowhere: left alone
    a.shown, b.shown = c["values"], None
    assert p.refresh() == [p.panels[1]]
    assert p.panels[0].content is c["values"] and not p.panels[0].stale
    assert p.panels[1].content is c["heat"] and p.panels[1].stale
    assert p.panels[2].content is c["tree"] and not p.panels[2].stale
    b.shown = c["freq"]
    assert p.refresh() == [] and p.panels[1].content is c["freq"]


def _dataset(samples: int = 3) -> Dataset:
    names = [f"s{k}" for k in range(samples)]
    return Dataset(
        X=np.zeros((10 * samples, 2), np.float32), channels=["FITC-A", "PE-A"],
        markers=["CD8", "CD4"], compensated=[True, True],
        sample_index=np.repeat(np.arange(samples), 10).astype(np.int32),
        event_index=np.arange(10 * samples), samples=names, meta=[{} for _ in names],
        provenance={
            "population": "Lymphocytes/CD4", "seed": 7, "mode": "equal_per_sample",
            "take": {n: 10 for n in names}, "available": {n: 10 + k for k, n in enumerate(names)},
            "excluded": {"tiny": "below the minimum of 5"}, "limited_by": "s0",
            "normalization": {
                "batch": "Batch", "reference": ["Type", "ref"], "clusters": 1,
                "markers": ["CD8", "CD4"], "batches": ["1", "2"],
                "reference_samples": {"1": ["r1"], "2": ["r2"]},
                "reference_events": {"1": 1000, "2": 1200}, "clusters_used": 1,
                "distance_before": [0.05, 0.03], "distance_after": [0.01, 0.01],
                "fallbacks": [], "not_normalised": ["s2"]}},
        runs={"FlowSOM": {"channels": ["FITC-A", "PE-A"], "transform": "biex", "seed": 3,
                          "n_clusters": 5, "grid": 10},
              "UMAP": {"channels": ["FITC-A", "PE-A"], "transform": "arcsinh", "cofactor": 150.0,
                       "seed": 3, "n_neighbors": 15, "min_dist": 0.5,
                       "neighbours_kept": 0.4321}})


def test_summary_says_what_is_recorded():
    lines = plate.summary(_dataset())
    text = "\n".join(lines)
    assert lines[0] == "30 events from 3 sample(s)."
    assert "Population: Lymphocytes/CD4." in lines
    assert "The same number of events from every sample; seed 7." in lines
    assert "The number was set by s0." in lines and "  s1: 10 of 11" in lines
    assert "  tiny: left out (below the minimum of 5)" in lines
    assert "Batches by the column 'Batch'; reference samples: Type = ref." in lines
    assert "0.0400 before, 0.0100 after" in text and "s2" in lines[-5]
    assert lines[-2:] == [
        "UMAP: on CD8, CD4; scale: arcsinh (cofactor 150); n_neighbors 15; min_dist 0.5; "
        "seed 3. Keeps 43.2% of each event's 15 nearest neighbours.",
        "FlowSOM: on CD8, CD4; scale: biex; 5 clusters on a 10 × 10 grid; seed 3."]
    # a file somebody else concatenated records only where it came from
    bare = Dataset(np.zeros((4, 1), np.float32), ["A"], [""], [False], np.zeros(4, np.int32),
                   np.arange(4), ["x"], [{}], provenance={"source": "x"})
    assert plate.summary(bare) == ["4 events from 1 sample(s).",
                                   "Read from an already concatenated file: x."]


def test_report_pages(tmp_path):
    contents = _contents()
    style = figures.FigureStyle(font=FONT, width_mm=85, height_mm=65)
    p = plate.Plate(style=style)
    p.add(contents["map"])
    p.add(contents["heat"])
    pages = [p, (contents["freq"], style), (contents["tree"], style)]
    path = str(tmp_path / "r.pdf")
    assert plate.report(path, pages, _dataset()) == 4 and _pdf_pages(path) == 4
    assert plate.report(path, pages) == 3 and _pdf_pages(path) == 3   # nothing to say first
    # more samples than a page of text holds: the text goes on
    assert plate.report(path, [p], _dataset(samples=120)) == 4 and _pdf_pages(path) == 4
    with open(path, "rb") as f:
        assert b"/Type /Font" in f.read()                # text, not outlines


# --- the window ----------------------------------------------------------------------

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
    rng = np.random.default_rng(0)
    centres = np.array([[300, 3000, 100], [3000, 200, 100], [200, 300, 4000.0]])
    for name, n in (("3Ctrl1", 240), ("3Ctrl2", 400), ("3Case1", 900), ("3Case2", 1500)):
        pop = rng.integers(0, 3, n)
        fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (n, 3)))
        ev = np.column_stack([rng.normal(90000, 9000, n), fluo]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), ev, ["FSC-A", "FITC-A", "PE-A", "APC-A"],
                  markers=["", "CD8", "CD4", "Foxp3"],
                  keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
    return str(tmp_path)


def test_the_window_collects_orders_refreshes_and_exports(app, folder, tmp_path):
    from milc.ui import worker
    from milc.ui.figure_dialog import FigureDialog
    from milc.ui.shell import MainShell

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(folder)
    w.samples.concatenate()
    settle()
    e, c, g = w.explore, w.clusters, w.gates
    e.clusters.setValue(3)
    e.grid.setValue(5)
    e.run(["UMAP", "FlowSOM"])
    settle()
    e.color_by.setCurrentIndex(e.color_by.findData("meta|Group"))
    w.show()
    w.tabs.setCurrentWidget(g)                          # the gates draw a sample when looked at
    settle()
    g.set_axes("Comp-PE-A", "Comp-FITC-A")

    # a figure dialog adds what it shows, as it shows it
    d = FigureDialog(e.view.content(), w, "map")
    d.source = e.view                                   # as `open_for` does
    assert d.plate_button.isVisibleTo(d)
    d.font.setCurrentIndex(d.font.findText(FONT))
    d.legend.setCurrentIndex(d.legend.findData("bottom"))
    d.title.setText("The map")
    d.add_to_plate()
    assert "panel A" in d.note.text()
    d.close()
    pw = w.plate_window()
    assert pw is w.plate_window() and len(pw.plate.panels) == 1
    first = pw.plate.panels[0]
    assert first.title == "The map" and first.style["legend"] == "bottom"
    assert first.source is e.view and pw.plate.style.font == FONT   # the first sets the type
    for view, name in ((c.heatmap, "heatmap"), (c.freq, "frequencies"), (c.tree, "tree"),
                       (g, "gates")):
        d = FigureDialog(view.content(), w, name)
        d.source = view
        d.add_to_plate()
        d.close()
    assert not FigureDialog(e.view.content(), None).plate_button.isVisibleTo(d)
    assert pw.list.count() == 5 and pw.list.item(0).text().startswith("A   map — The map")
    assert [p.name for p in pw.plate.panels] == ["map", "heatmap", "frequencies", "tree", "gates"]

    # order, titles, spans, the page
    pw.list.setCurrentRow(4)
    pw.move(-1)
    pw.move(-1)
    assert [p.name for p in pw.plate.panels][2:] == ["gates", "frequencies", "tree"]
    assert pw.list.currentRow() == 2 and pw.list.item(2).text().startswith("C   gates")
    pw.title.setText("One sample")
    pw.colspan.setValue(2)
    assert pw.plate.panels[2].title == "One sample" and pw.plate.panels[2].colspan == 2
    assert "One sample" in pw.list.item(2).text()
    pw.list.setCurrentRow(4)
    pw.remove()
    assert [p.name for p in pw.plate.panels] == ["map", "heatmap", "gates", "frequencies"]
    assert pw.list.currentRow() == 3
    pw.columns.setValue(3)
    pw.width_preset.setCurrentIndex(pw.width_preset.findData(120))
    assert pw.plate.columns == 3 and pw.plate.width_mm == 120 and pw.plate.height_mm is None
    assert pw.plate.cells() == [(0, 0, 1, 1), (0, 1, 1, 1), (1, 0, 1, 2), (1, 2, 1, 1)]
    pw.width.setValue(130)
    assert pw.width_preset.currentData() is None        # no longer one of the presets
    pw.width_preset.setCurrentIndex(pw.width_preset.findData(180))
    pw.height.setValue(140)
    pw.render()
    assert "180 × 140 mm" in pw.note.text() and not pw.preview.pixmap().isNull()

    svg_path = str(tmp_path / "plate.svg")
    pw.write(svg_path)
    written = re.findall(r"<text[^>]*>([^<]*)</text>",
                         open(svg_path, encoding="utf-8").read())
    for word in ("A", "B", "C", "D", "The map", "One sample", "Case", "Ctrl", "UMAP1", "CD4"):
        assert word in written, word
    pw.letters.setChecked(False)
    pw.write(svg_path)
    assert "D" not in re.findall(r"<text[^>]*>([^<]*)</text>",
                                 open(svg_path, encoding="utf-8").read())
    pw.letters.setChecked(True)

    # the report: what the dataset records, then the page
    report = str(tmp_path / "report.pdf")
    assert pw.write_report(report) == 2 and _pdf_pages(report) == 2
    said = "\n".join(plate.summary(w.session.dataset))
    assert "960 events from 4 sample(s)" in said and "set by 3Ctrl1" in said
    assert "UMAP: on CD8, CD4, Foxp3" in said and "3 clusters on a 5 × 5 grid" in said

    # the analysis moves on: the panels are asked again
    e.color_by.setCurrentIndex(e.color_by.findData("cluster|FlowSOM"))
    assert pw.plate.panels[0].content.color.title == "Group"
    pw.refresh_from_analysis()
    assert pw.plate.panels[0].content.color.names == ["1", "2", "3"]
    assert pw.plate.panels[0].style["legend"] == "bottom"   # its look is kept
    assert not any(p.stale for p in pw.plate.panels)
    kept = pw.plate.panels[1].content
    c.heatmap.table = None                              # a plot with nothing to show
    pw.refresh_from_analysis()
    assert pw.plate.panels[1].stale and pw.plate.panels[1].content is kept
    assert "stale" in pw.list.item(1).text() and "1 panel(s) kept" in pw.note.text()
    pw.write(svg_path)                                  # and it is still drawn

    w.show_plate()
    assert pw.isVisible()
    pw.close()
