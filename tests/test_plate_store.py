"""Plates kept with the project: every kind of content comes back as it
was, a kind nobody told the store about too; what cannot come back is a
panel that says so; a project from before the plates opens as it did; and
the window makes, names, copies and deletes plates that are there again
when the folder is opened."""

import dataclasses
import enum
import json
import os
from dataclasses import dataclass, field

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from milc import abundance, figures, plate, plate_store, views  # noqa: E402
from milc.gating import GateNode, PolygonGate, RectangleGate  # noqa: E402
from milc.session import PROJECT_FILE, Session  # noqa: E402
from milc.transforms import Biex, decade_ticks  # noqa: E402
from milc.volcano import VolcanoContent  # noqa: E402
from test_plate import FONT, _contents, _pdf_pages, app, folder  # noqa: E402,F401


def _all_contents() -> dict:
    """Every kind of content that can be put on a plate today."""
    rng = np.random.default_rng(1)
    c = _contents()
    b = Biex()
    c["freq"].facets = ["x", "x", "y", "y"]
    c["freq"].notes = {(None, "1"): "*", ("x", "2"): "ns"}     # keys that are tuples
    c["tree"].tree["names"] = {1: "Tregs"}                      # and numbers
    lines = [np.column_stack([np.cos(t), np.sin(t)]) * r
             for r in (1.0, 0.5) for t in [np.linspace(0, 6.28, 90)]]
    c["contours"] = figures.OverlayContent(
        [figures.Series("s1", "Case", levels=[(0.9, lines[:1]), (0.5, lines[1:])],
                        outliers=rng.normal(0, 1, (30, 2)))],
        "CD8", "CD4", xticks=decade_ticks(b), xrange=(-1.5, 1.5), yrange=(-1.5, 1.5))
    series = [views.distribution(b(rng.lognormal(6 + k, 1.0, 3000)), 0.0, 1.0, name, group)
              for k, (name, group) in enumerate((("3Ctrl1", "Ctrl"), ("3Case1", "Case")))]
    series.append(views.distribution(np.array([]), 0.0, 1.0, "nobody", "Case"))   # NaN inside
    c["violin"] = views.DistributionContent(series, views.VIOLIN, "CD4", decade_ticks(b),
                                            (0.0, 1.0))
    c["box"] = views.DistributionContent(series, views.BOX, "CD4", decade_ticks(b), (0.0, 1.0))
    ranges = [(0.0, 1.0)] * 3
    c["matrix"] = views.MatrixContent(
        ["CD8", "CD4", "Foxp3"],
        views.pairs([b(rng.lognormal(6, 1.2, 2000)) for _ in range(3)], ranges),
        ranges, [decade_ticks(b)] * 3)
    events = {"FSC": rng.uniform(0, 1000, 4000), "SSC": rng.uniform(0, 1000, 4000),
              "A": rng.uniform(0, 100, 4000), "B": rng.uniform(0, 100, 4000)}
    tree = GateNode("")
    cells = tree.add(GateNode("Cells", RectangleGate(["FSC", "SSC"], [200.0, None],
                                                     [None, 800.0])))
    square = np.array([[20.0, 20.0], [80.0, 20.0], [80.0, 80.0], [20.0, 80.0]])
    cells.add(GateNode("Poly", PolygonGate(["A", "B"], square)))
    c["backgate"] = views.backgate(
        tree, "Cells/Poly", tree.apply(lambda d: events[d], 4000), lambda d: events[d],
        lambda d: None, lambda d: (0.0, 1000.0 if d in ("FSC", "SSC") else 100.0),
        lambda d: f"{d} marker")
    counts = rng.poisson(200, (8, 10)).astype(float)
    counts[0, 5:] *= 4
    c["volcano"] = VolcanoContent(abundance.differential(
        counts, [str(i + 1) for i in range(8)], ["A"] * 5 + ["B"] * 5, ("A", "B")))
    return c


def same(a, b, where="") -> None:
    """`a` and `b` alike all the way down: the kinds, the numbers, the
    types of the arrays, a NaN where there was one."""
    assert type(a) is type(b), (where, type(a), type(b))
    if isinstance(a, np.ndarray):
        assert a.dtype == b.dtype and a.shape == b.shape, (where, a.dtype, b.dtype)
        if a.dtype.kind == "O":
            for k, (x, y) in enumerate(zip(a.ravel(), b.ravel())):
                same(x, y, f"{where}[{k}]")
        else:
            assert np.array_equal(a, b, equal_nan=a.dtype.kind == "f"), where
    elif dataclasses.is_dataclass(a):
        for f in dataclasses.fields(a):
            same(getattr(a, f.name), getattr(b, f.name), f"{where}.{f.name}")
    elif isinstance(a, dict):
        assert list(a) == list(b), where
        for k in a:
            same(a[k], b[k], f"{where}[{k!r}]")
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b), where
        for k, (x, y) in enumerate(zip(a, b)):
            same(x, y, f"{where}[{k}]")
    elif isinstance(a, float) and a != a:
        assert b != b, where
    else:
        assert a == b, (where, a, b)


def _through_json(plates, current=0):
    record, arrays = plate_store.dump(plates, current)
    return plate_store.restore(json.loads(json.dumps(record, allow_nan=False)), arrays)


def test_every_content_round_trips_and_draws_the_same():
    contents = _all_contents()
    assert {type(c) for c in contents.values()} == set(figures._DRAW) - {
        plate_store.Missing} - {t for t in figures._DRAW if t.__module__ == __name__}
    page = figures.FigureStyle(font=FONT, font_size=7.0, dpi=200)
    first = plate.Plate(columns=3, width_mm=180, style=page, name="Figure 1")
    second = plate.Plate(columns=4, width_mm=240, height_mm=200.0, letters=False,
                         style=figures.FigureStyle(font=FONT, colors={"Case": "#000000"}),
                         name="Supplement")
    for k, (name, content) in enumerate(contents.items()):
        p = (first, second)[k % 2]
        look = figures.FigureStyle(font=FONT, legend=("bottom", "right", "none")[k % 3],
                                   colormap="magma", point_size=1.0 + k,
                                   labels={"Case": "Treated"}, title=f"The {name}")
        panel = p.add(content, look, name=name)
        panel.colspan, panel.rowspan = 1 + k % 2, 1 + (k == 3)
    first.panels[1].stale = True

    plates, current = _through_json([first, second], current=1)
    assert current == 1 and [p.name for p in plates] == ["Figure 1", "Supplement"]
    for before, after in zip((first, second), plates):
        for f in ("columns", "width_mm", "height_mm", "letters", "name", "style"):
            assert getattr(before, f) == getattr(after, f), f
        assert after.cells() == before.cells() and after.size_mm() == before.size_mm()
        for p, q in zip(before.panels, after.panels, strict=True):
            assert (p.title, p.name, p.style, p.colspan, p.rowspan) == (
                q.title, q.name, q.style, q.colspan, q.rowspan)
            # what it is, not where it came from: no plot behind it to be stale against
            assert q.source is None and not q.stale
            same(p.content, q.content, p.name)
        assert plate.preview_png(before, dpi=50) == plate.preview_png(after, dpi=50)
    assert plates[1].style.colors == {"Case": "#000000"}


# --- a kind of content the store was never told about ---------------------------------

class Shade(enum.Enum):
    LIGHT = 1
    DARK = 2


@dataclass
class Part:
    weights: np.ndarray
    shade: Shade = Shade.LIGHT


@dataclass(frozen=True)
class Later:
    """As another module may add one day: registered in `figures._DRAW`
    and nothing more."""
    parts: list
    big: np.ndarray
    flags: np.ndarray
    words: np.ndarray
    mixed: np.ndarray
    one: np.float32
    yes: np.bool_
    limits: tuple = (0.0, float("inf"))
    missing: float = float("nan")
    chosen: set = field(default_factory=set)
    by_number: dict = field(default_factory=dict)
    shade: Shade | None = None
    nothing: None = None
    made: int = field(init=False, default=0)


def _draw_later(fig, c, style) -> None:
    fig.subplots().plot(c.big[:, 0], c.big[:, 1])


def _later() -> Later:
    mixed = np.empty(3, dtype=object)
    mixed[:] = ["a", 2, None]
    c = Later([Part(np.arange(3.0)), Part(np.full(200, np.nan), Shade.DARK)],
              np.random.default_rng(0).normal(0, 1, (500, 2)).astype(np.float32),
              np.array([True, False]), np.array(["CD4", "Foxp3"]), mixed,
              np.float32(0.1), np.bool_(True), chosen={"a", "b"},
              by_number={3: "three", (1, 2): [1, 2]}, shade=Shade.DARK)
    object.__setattr__(c, "made", 7)
    return c


@pytest.fixture
def later():
    figures._DRAW[Later] = _draw_later
    yield _later()
    figures._DRAW.pop(Later, None)


def test_a_content_added_later_round_trips_unannounced(later):
    p = plate.Plate(style=figures.FigureStyle(font=FONT), name="x")
    p.add(later, title="New")
    record, arrays = plate_store.dump([p])
    # the big arrays are beside the JSON, the small ones in it
    assert len(arrays) == 2 and sorted(a.size for a in arrays.values()) == [200, 1000]
    assert len(json.dumps(record)) < 3000
    (back,), _ = plate_store.restore(json.loads(json.dumps(record, allow_nan=False)), arrays)
    got = back.panels[0].content
    same(later, got)
    assert got.made == 7 and got.by_number[(1, 2)] == [1, 2] and got.shade is Shade.DARK
    assert plate.preview_png(back, dpi=40) == plate.preview_png(p, dpi=40)


# --- what cannot come back -------------------------------------------------------------

def _texts(p: plate.Plate) -> list[str]:
    from matplotlib.text import Text
    with plate._drawn(p) as fig:
        return [t.get_text() for t in fig.findobj(Text) if t.get_text()]


def test_a_panel_that_cannot_be_restored_is_a_placeholder(later, tmp_path):
    folder = str(tmp_path)
    contents = _contents()
    p = plate.Plate(style=figures.FigureStyle(font=FONT), name="Figure")
    p.add(contents["map"], title="The map", name="map")
    p.add(later, title="The new one", name="later")
    p.add(contents["heat"], title="The heatmap", name="heatmap")
    p.panels[1].colspan = 2
    record = json.loads(json.dumps(plate_store.write(folder, [p])))

    # its kind is not known to this program (an older one, or a module not there)
    del figures._DRAW[Later]
    (back,), _ = plate_store.read(folder, record)
    gone = back.panels[1]
    assert isinstance(gone.content, plate_store.Missing) and "Later" in gone.content.note
    assert (gone.title, gone.name, gone.colspan) == ("The new one", "later", 2)
    assert not isinstance(back.panels[0].content, plate_store.Missing)
    texts = _texts(back)
    assert "The new one" in texts and "B" in texts
    assert any("could not be restored" in t for t in texts)
    assert figures.category_names(gone.content) == [] and not figures.uses_colormap(gone.content)
    # saved again as it is, nothing of it is lost: it is there when its kind is
    again = json.loads(json.dumps(plate_store.write(folder, [back])))
    assert again["plates"][0]["panels"][1] == record["plates"][0]["panels"][1]
    figures._DRAW[Later] = _draw_later
    (whole,), _ = plate_store.read(folder, again)
    same(later, whole.panels[1].content)

    # the arrays are gone: every panel that had some says so, by its title
    os.remove(os.path.join(folder, plate_store.FILE))
    (bare,), _ = plate_store.read(folder, record)
    assert [isinstance(q.content, plate_store.Missing) for q in bare.panels] == [
        True, True, False]                              # a small table is in the JSON itself
    assert [q.title for q in bare.panels] == ["The map", "The new one", "The heatmap"]
    assert plate_store.FILE in bare.panels[0].content.note
    assert plate.preview_png(bare, dpi=40)[:4] == b"\x89PNG"

    # a content that cannot be written costs its own panel, not the save
    @dataclass
    class Odd:
        thing: object

    figures._DRAW[Odd] = _draw_later
    try:
        q = plate.Plate(name="odd")
        q.add(Odd(lambda: 0), title="Odd one")
        q.add(contents["heat"], title="Fine")
        (back,), _ = plate_store.read(folder, json.loads(json.dumps(
            plate_store.write(folder, [q]))))
    finally:
        del figures._DRAW[Odd]
    assert isinstance(back.panels[0].content, plate_store.Missing)
    assert "Odd could not be stored" in back.panels[0].content.note
    assert back.panels[0].title == "Odd one"
    same(contents["heat"], back.panels[1].content)

    # and a project somebody's hand has been in is not an error
    for record in (None, "x", 3, {"plates": 4}, {"plates": [3, {"panels": [5, {}]}]},
                   {"plates": [{"columns": "many"}, {"name": "ok", "panels": [
                       {"title": "T", "kind": "MapContent",
                        "content": {"$": "data", "type": "os:system", "v": {}}},
                       {"content": {"$": "data", "type": "milc.figures:Nope", "v": {}}},
                       {"content": {"$": "what"}, "colspan": "x"}]}], "current": "z"}):
        plates, current = plate_store.read(folder, record)
        assert current == 0
        for q in plates:
            assert all(isinstance(x.content, plate_store.Missing) for x in q.panels)
            plate.preview_png(q, dpi=30) if q.panels else None
    assert [len(q.panels) for q in plates] == [3] and plates[0].panels[0].title == "T"


def test_the_arrays_are_written_once_and_only_when_they_change(tmp_path, monkeypatch):
    folder, path = str(tmp_path), str(tmp_path / plate_store.FILE)
    contents = _contents()
    writes = []
    real = plate_store._write_arrays
    monkeypatch.setattr(plate_store, "_write_arrays",
                        lambda f, arrays: (writes.append(len(arrays)), real(f, arrays))[1])
    assert plate_store.write(folder, []) is None and not os.path.exists(path)
    p = plate.Plate(name="a")
    p.add(contents["map"])
    p.add(contents["values"])                           # the same x and y, coloured otherwise
    plate_store.write(folder, [p])
    with np.load(path) as held:
        assert len(held.files) == 4                     # x, y, the group once, the marker
    assert writes == [4]
    # the page changes, a copy is made: the arrays are the ones already there
    p.columns, p.panels[0].title = 3, "Renamed"
    record = plate_store.write(folder, [p, p.copy("b")], current=1)
    assert writes == [4] and record["plates"][1]["name"] == "b" and record["current"] == 1
    p.add(figures.MapContent(np.arange(300.0), np.arange(300.0)[::-1]))
    plate_store.write(folder, [p])
    assert writes == [4, 6] and not os.path.exists(path + ".part")
    del p.panels[:]
    plate_store.write(folder, [p])                      # an empty plate has no arrays
    assert not os.path.exists(path)
    (back,), _ = plate_store.read(folder, plate_store.write(folder, [p]))
    assert back.name == "a" and back.panels == [] and back.columns == 3


def test_a_copy_is_changed_apart():
    p = plate.Plate(name="a", style=figures.FigureStyle(font=FONT))
    first = p.add(_contents()["map"], figures.FigureStyle(font=FONT, legend="bottom"))
    q = p.copy("b")
    q.panels[0].title, q.panels[0].style["legend"], q.style.font_size = "Other", "none", 12.0
    q.style.colors["Case"] = "#000000"
    q.add(first.content)
    assert (first.title, first.style["legend"], p.style.font_size) == ("", "bottom", 8.0)
    assert p.style.colors == {}
    assert len(p.panels) == 1 and q.panels[0].content is first.content and q.name == "b"


# --- with the project --------------------------------------------------------------------

def test_projects_without_plates_and_what_else_a_project_holds(folder):
    old = {"version": 2, "columns": ["Group", "Day"], "group_by": ["Group"],
           "annotations": {"3Ctrl1.fcs": {"Group": "Ctrl", "Day": "3"}},
           "excluded": ["3Case2.fcs"], "gates": {}, "scales": {}, "normalization": None,
           "compensation": None}
    with open(os.path.join(folder, PROJECT_FILE), "w", encoding="utf-8") as f:
        json.dump(old, f)
    s = Session()
    assert s.plates == [] and not s.save_plates()       # no folder: nowhere to write
    s.open_folder(folder)
    assert s.plates == [] and s.current_plate == 0
    assert s.columns == ["Group", "Day"] and not s.entries[1].include
    s.save_project()
    saved = json.load(open(s.project_file(), encoding="utf-8"))
    assert "plates" not in saved and not os.path.exists(os.path.join(folder, plate_store.FILE))
    assert set(saved) == set(old) and saved["columns"] == old["columns"]

    p = plate.Plate(name="Figure 1", style=figures.FigureStyle(font=FONT))
    p.add(_contents()["map"], title="Map")
    s.plates.append(p)
    assert s.save_plates()
    saved = json.load(open(s.project_file(), encoding="utf-8"))
    # what was there is as it was, with the plates beside it
    assert set(saved) == set(old) | {"plates"} and saved["version"] == 2
    assert {k: saved[k] for k in ("columns", "group_by", "excluded")} == {
        k: old[k] for k in ("columns", "group_by", "excluded")}
    # whatever else writes the project keeps them
    s.set_scale(["FITC-A"], None)
    assert json.load(open(s.project_file(), encoding="utf-8"))["plates"] == saved["plates"]

    again = Session()
    again.open_folder(folder)
    assert [q.name for q in again.plates] == ["Figure 1"] and again.columns == ["Group", "Day"]
    same(p.panels[0].content, again.plates[0].panels[0].content)
    del again.plates[:]
    again.save_plates()
    assert "plates" not in json.load(open(s.project_file(), encoding="utf-8"))
    assert not os.path.exists(os.path.join(folder, plate_store.FILE))


def test_plates_in_the_window_and_back_from_the_project(app, folder, tmp_path):  # noqa: F811
    from milc.ui import worker
    from milc.ui.figure_dialog import FigureDialog
    from milc.ui.shell import MainShell

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    def add(w, view, name, to=None, choose=None):
        d = FigureDialog(view.content(), w, name)
        d.source = view
        d.font.setCurrentIndex(d.font.findText(FONT))
        if choose is not None:
            d.plate_choice.setCurrentIndex(d.plate_choice.findText(choose))
        d.add_to_plate(to)
        return d

    w = MainShell()
    w.open_folder(folder)
    w.samples.concatenate()
    settle()
    e, c = w.explore, w.clusters
    e.clusters.setValue(3)
    e.grid.setValue(5)
    e.run(["UMAP", "FlowSOM"])
    settle()
    assert w.session.plates == []

    # the first panel makes the first plate; the dialog offers it and a new one
    d = FigureDialog(e.view.content(), w, "map")
    assert [d.plate_choice.itemText(i) for i in range(d.plate_choice.count())] == [
        "Plate 1", "New plate"]
    assert d.plate_choice.isVisibleTo(d) and d.plate_choice.currentIndex() == 0
    d = add(w, e.view, "map")
    assert "“Plate 1” as panel A" in d.note.text()
    pw = w.plate_window()
    assert [p.name for p in pw.plates] == ["Plate 1"] and pw.plates is w.session.plates
    add(w, c.heatmap, "heatmap")
    pw.rename_plate("Clusters")
    assert pw.plates_box.currentText() == "Clusters" and "Clusters" in pw.windowTitle()

    # a new plate from the dialog, which becomes the one shown
    d = add(w, c.freq, "frequencies", choose="New plate")
    assert "“Plate 2” as panel A" in d.note.text() and "Plate 2" in w.statusBar().currentMessage()
    assert [d.plate_choice.itemText(i) for i in range(d.plate_choice.count())] == [
        "Clusters", "Plate 2", "New plate"] and d.plate_choice.currentText() == "Plate 2"
    assert w.session.current_plate == 1 and pw.list.count() == 1
    assert pw.plate.style.font == FONT                  # the type of the plate it was made from
    # and a plate that is not the current one, by its place or by its name
    add(w, c.tree, "tree", to=0)
    add(w, c.freq, "frequencies", choose="Plate 2")
    assert [[p.name for p in q.panels] for q in pw.plates] == [
        ["map", "heatmap", "tree"], ["frequencies", "frequencies"]]

    # switch, duplicate, rename, delete
    pw.plates_box.setCurrentIndex(0)
    assert w.session.current_plate == 0 and pw.list.count() == 3
    pw.columns.setValue(3)
    pw.list.setCurrentRow(1)
    pw.title.setText("Markers")
    pw.colspan.setValue(2)
    copy = pw.duplicate_plate()
    assert [p.name for p in pw.plates] == ["Clusters", "Clusters copy", "Plate 2"]
    assert pw.plate is copy and pw.plates_box.currentIndex() == 1 and pw.list.count() == 3
    pw.list.setCurrentRow(0)
    pw.remove()
    pw.letters.setChecked(False)
    assert len(pw.plates[0].panels) == 3 and pw.plates[0].letters   # the original is as it was
    pw.rename_plate("Plate 2")                          # taken: two of a name cannot be told apart
    assert pw.plate.name == "Plate 2 2"
    pw.rename_plate("  ")
    assert pw.plate.name == "Plate 2 2"
    pw.rename_plate("Short")
    pw.new_plate("Empty")
    assert pw.plate.panels == [] and pw.plate.columns == 3 and not pw.plate.letters
    pw.render()
    assert "No panels yet" in pw.note.text()
    pw.new_plate()
    assert [p.name for p in pw.plates] == ["Clusters", "Short", "Plate 2", "Empty", "Plate 5"]
    pw.delete_plate()
    assert [p.name for p in pw.plates] == ["Clusters", "Short", "Plate 2", "Empty"]
    assert pw.plate.name == "Empty" and pw.plates_box.count() == 4

    # the report: this plate, or the ones chosen, a page each behind the text
    pw.switch_to(0)
    report = str(tmp_path / "report.pdf")
    assert pw.write_report(report) == 2
    assert pw.write_report(report, [p for p in pw.plates if p.panels]) == 4
    assert _pdf_pages(report) == 4
    assert pw.write_report(report, pw.plates[1:2]) == 2
    before = [plate.preview_png(p, dpi=50) for p in pw.plates if p.panels]

    # what is on the page is written with the project: at once when a panel
    # is added, and when the changes pause or the window goes
    assert pw._saver.isActive()
    w.close()
    assert not pw._saver.isActive()
    saved = json.load(open(w.session.project_file(), encoding="utf-8"))["plates"]
    assert [p["name"] for p in saved["plates"]] == ["Clusters", "Short", "Plate 2", "Empty"]
    assert saved["current"] == 0 and os.path.exists(os.path.join(folder, plate_store.FILE))

    # another session opens the folder: no analysis, and the plates as they were
    w2 = MainShell()
    w2.open_folder(folder)
    assert w2.session.dataset is None
    pw2 = w2.plate_window()
    assert [p.name for p in pw2.plates] == ["Clusters", "Short", "Plate 2", "Empty"]
    assert pw2.plates_box.currentText() == "Clusters" and pw2.list.count() == 3
    assert pw2.list.item(1).text().startswith("B   heatmap — Markers")
    assert pw2.plates[0].panels[1].colspan == 2 and not pw2.plates[1].letters
    assert [plate.preview_png(p, dpi=50) for p in pw2.plates if p.panels] == before
    w2.show_plate()
    pw2.render()
    assert not pw2.preview.pixmap().isNull()
    # they are what they were, and no longer follow a plot
    pw2.refresh_from_analysis()
    assert "3 panel(s) came with the project" in pw2.note.text()
    d = FigureDialog(pw2.plates[0].panels[0].content, w2, "map")
    assert d.plate_choice.currentText() == "Clusters" and d.plate_choice.count() == 5
    d.add_to_plate()
    assert len(pw2.plates[0].panels) == 4 and pw2.list.count() == 4

    # a folder without plates has none, and this one's are not written into it
    other = tmp_path / "other"
    other.mkdir()
    w2.open_folder(str(other))
    assert [p.panels for p in pw2.plates] == [[]] and pw2.list.count() == 0
    assert not os.path.exists(other / PROJECT_FILE)
    w2.open_folder(folder)
    assert len(pw2.plates) == 4 and len(pw2.plate.panels) == 4
    pw2.close()
    w2.close()


def test_a_long_axis_name_in_a_short_cell_keeps_clear_of_the_letter():
    contents = _contents()
    style = figures.FigureStyle(font=FONT)
    for letters, title in ((True, ""), (True, "Frequencies"), (False, "Frequencies")):
        p = plate.Plate(columns=3, width_mm=180, style=style, letters=letters)
        for name in ("map", "heat", "freq"):
            p.add(contents[name], title=title)
        with plate._drawn(p) as fig:
            fig.draw_without_rendering()
            sub = fig.subfigs[2]
            label = sub.get_axes()[0].yaxis.label
            assert label.get_text().replace("\n", " ") == "% of the sample's events"
            # too long for its axis, and in the way of the letter: on several lines
            assert ("\n" in label.get_text()) == letters
            box = label.get_window_extent()
            heads = [sub._suptitle] + list(sub.texts)
            assert heads and not any(box.overlaps(h.get_window_extent()) for h in heads)
            assert sub.bbox.y0 <= box.y0 and box.y1 <= sub.bbox.y1
            # a name that fits is left alone
            assert fig.subfigs[0].get_axes()[0].yaxis.label.get_text() == "UMAP2"
    # and where the cell is tall enough nothing is broken
    p = plate.Plate(columns=1, width_mm=120, style=style)
    p.add(contents["freq"])
    with plate._drawn(p) as fig:
        assert fig.subfigs[0].get_axes()[0].yaxis.label.get_text() == "% of the sample's events"
