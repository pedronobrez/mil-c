"""Quality control of the acquisition: synthetic acquisitions, each with one
defect whose interval is known, and the gate that removes what is found.

The acquisitions run for 60 s at 1,000 events/s — bins of 150 events,
0.15 s each — with two scatter parameters and two fluorescence ones, the
latter with two populations each (40% of the cells in one). What was
measured when this was written is said beside each assertion; the
assertions leave a margin."""

import json
import os

import numpy as np
import pytest

from milc import qc
from milc.fcs import write_fcs
from milc.gating import (GateNode, RectangleGate, TimeGate, copy_tree, keep_times, time_node,
                         tree_from_dict, tree_to_dict, unwrap, wrap)
from milc.session import Session
from milc.transforms import Biex

SECONDS, RATE = 60.0, 1000.0
CHANNELS = ["FSC-A", "SSC-A", "FITC-A", "PE-A"]
SCALES = [None, None, Biex(), Biex()]
BIN = 0.15   # seconds a bin of 150 events lasts


def arrivals(rng, seconds=SECONDS, rate=RATE) -> np.ndarray:
    return np.sort(rng.uniform(0, seconds, rng.poisson(seconds * rate)))


def cells(rng, n: int, bright: float = 0.4) -> np.ndarray:
    bright = rng.random(n) < bright
    return np.column_stack([
        rng.normal(90000, 9000, n), rng.normal(40000, 6000, n),
        np.where(bright, 3000, 150) * rng.lognormal(0, 0.3, n) + rng.normal(0, 40, n),
        np.where(bright, 200, 5000) * rng.lognormal(0, 0.3, n) + rng.normal(0, 40, n)])


def acquisition(seed: int):
    rng = np.random.default_rng(seed)
    t = arrivals(rng)
    return rng, t, cells(rng, len(t))


def check(x, t, **kw) -> qc.Result:
    return qc.check(x, t, CHANNELS, SCALES, **kw)


def removed_within(r: qc.Result, t, lo, hi) -> float:
    """Share of the events acquired in [lo, hi] that the check removes."""
    there = (t >= lo) & (t <= hi)
    return float((~r.keep[there]).mean())


def removed_outside(r: qc.Result, t, lo, hi) -> float:
    there = (t < lo) | (t > hi)
    return float((~r.keep[there]).mean())


# --- clean data ------------------------------------------------------------------

def test_clean_acquisitions_are_left_almost_untouched():
    """Whatever share of the cells is in each population of a marker —
    half and half is where the median of all of them would jump about."""
    for bright in (0.4, 0.48, 0.25, 0.1):
        lost = total = touched = 0
        for seed in range(12):
            rng = np.random.default_rng(seed)
            t = arrivals(rng)
            r = check(cells(rng, len(t), bright), t)
            assert r.checked and r.removed < 0.02
            lost, total = lost + int((~r.keep).sum()), total + len(t)
            touched += bool(r.intervals)
        # measured: nothing removed from any of the 48 files here. Over 50
        # seeds at 40% and 40 at each of 10, 20, 30 and 48%, with the two
        # fluorescence parameters alone: one file lost 450 events, at 10%.
        assert lost / total < 0.001 and touched <= 1


def test_each_population_of_a_parameter_is_followed():
    _, t, x = acquisition(0)
    r = check(x, t)
    assert r.bins.channel.tolist() == [0, 1, 2, 2, 3, 3]
    assert r.populations[:2] == ["FSC-A", "SSC-A"]
    assert r.populations[2:4] == ["FITC-A, around 150", "FITC-A, around 3,000"]
    # each on its own scale, and where the population is: 150 and 3,000
    inverse = Biex().inverse
    assert np.median(r.bins.medians[:, 0]) == pytest.approx(90000, rel=0.01)
    assert inverse(np.median(r.bins.medians[:, 2])) == pytest.approx(150, rel=0.1)
    assert inverse(np.median(r.bins.medians[:, 3])) == pytest.approx(3000, rel=0.05)
    # a population too small to give every bin a median is left alone
    rng = np.random.default_rng(0)
    few = check(cells(rng, len(t), 0.04), t)
    assert few.bins.channel.tolist() == [0, 1, 2, 3]


# --- one defect each ------------------------------------------------------------------

def test_clog():
    """The flow rate falls to 3% of the usual for ten seconds."""
    rng, t, x = acquisition(1)
    stay = ~((t >= 20) & (t <= 30)) | (rng.random(len(t)) < 0.03)
    t, x = t[stay], x[stay]
    r = check(x, t)
    assert len(r.intervals) == 1 and r.counts[qc.RATE] > 0
    lo, hi = r.intervals[0]
    # measured: 19.92 to 30.09 — the bins that straddle the ends go too, and
    # they are the 0.34% lost outside; right to 0.3 s in 29 of 30 seeds
    assert abs(lo - 20) < 2 * BIN and abs(hi - 30) < 2 * BIN
    assert removed_within(r, t, 20, 30) == 1.0
    assert removed_outside(r, t, 20, 30) < 0.005


def test_burst():
    """Five times the usual flow rate for two seconds."""
    rng, t, x = acquisition(2)
    extra = rng.uniform(40, 42, rng.poisson(4 * RATE * 2))
    t = np.sort(np.concatenate([t, extra]))
    x = cells(rng, len(t))
    r = check(x, t)
    assert len(r.intervals) == 1 and r.counts[qc.RATE] > 0
    lo, hi = r.intervals[0]
    # measured: 40.01 to 42.03, 99.7% of the burst; right to 0.3 s in 30 of 30 seeds
    assert abs(lo - 40) < 2 * BIN and abs(hi - 42) < 2 * BIN
    assert removed_within(r, t, 40, 42) > 0.97
    assert removed_outside(r, t, 40, 42) < 0.005


def test_step_in_one_channel():
    """PE is 1.8 times brighter from 45 s to the end."""
    _, t, x = acquisition(3)
    x[t >= 45, 3] *= 1.8
    r = check(x, t)
    assert len(r.intervals) == 1 and r.blamed == ["PE-A"] and r.counts[qc.RATE] == 0
    lo, hi = r.intervals[0]
    # measured: from 44.93; within 0.45 s in 30 of 30 seeds
    assert abs(lo - 45) < 3 * BIN and hi == t.max()
    assert removed_within(r, t, 45, 60) > 0.98
    assert removed_outside(r, t, 45, 60) < 0.005


def test_slow_drift():
    """FITC climbs steadily over the last third, to three times its value.
    A drift is found from where it has gone further than the limit, not
    from where it began."""
    _, t, x = acquisition(4)
    x[:, 2] *= 1 + 2 * np.clip((t - 40) / 20, 0, None)
    r = check(x, t)
    assert len(r.intervals) == 1 and r.blamed == ["FITC-A"]
    lo, hi = r.intervals[0]
    # measured: from 43.8 s, where it is 38% brighter — 81% of the drift;
    # over 30 seeds from between 42.4 and 44.3 s
    assert 40 <= lo < 45 and hi == t.max()
    assert removed_within(r, t, 40, 60) > 0.75
    assert removed_outside(r, t, 40, 60) == 0.0


def test_air_bubble():
    """One second of low scatter."""
    _, t, x = acquisition(5)
    x[(t >= 10) & (t <= 11), :2] *= 0.3
    r = check(x, t)
    assert len(r.intervals) == 1 and set(r.blamed) == {"FSC-A", "SSC-A"}
    lo, hi = r.intervals[0]
    # measured: 9.97 to 11.10; within 0.45 s at both ends in 30 of 30 seeds
    assert abs(lo - 10) < 3 * BIN and abs(hi - 11) < 3 * BIN
    assert removed_within(r, t, 10, 11) > 0.9
    assert removed_outside(r, t, 10, 11) < 0.005


# --- the rules about neighbours ------------------------------------------------------

def test_one_noisy_bin_among_good_ones_stays():
    _, t, x = acquisition(6)
    x[1500:1650, :2] *= 0.3     # exactly the eleventh bin
    r = check(x, t)
    assert r.bins.medians[10, 0] < 40000 and not r.intervals
    x[1650:1800, :2] *= 0.3     # and two in a row
    assert not check(x, t).intervals
    x[1800:1950, :2] *= 0.3     # three are a defect
    assert len(check(x, t).intervals) == 1


def test_a_short_good_island_between_bad_stretches_goes_too():
    _, t, x = acquisition(7)
    x[((t >= 20) & (t < 23)) | ((t >= 23.45) & (t < 26)), 0] *= 0.3
    r = check(x, t)                      # three good bins between: fewer than five
    assert len(r.intervals) == 1 and r.counts[qc.ISLANDS] > 0
    assert removed_within(r, t, 20, 26) > 0.95
    assert len(check(x, t, consecutive=2).intervals) == 2
    _, t, x = acquisition(7)
    x[((t >= 20) & (t < 23)) | ((t >= 26) & (t < 29)), 0] *= 0.3
    r = check(x, t)                      # twenty good bins between: kept
    assert len(r.intervals) == 2 and r.counts[qc.ISLANDS] == 0
    assert removed_within(r, t, 23.5, 25.5) == 0.0


# --- the result ------------------------------------------------------------------------

def test_what_is_reported_is_what_the_gate_removes():
    _, t, x = acquisition(5)
    x[(t >= 10) & (t <= 11), :2] *= 0.3
    r = check(x, t)
    gate = r.gate("Time")
    assert np.array_equal(gate.contains({"Time": t}.__getitem__), r.keep)
    assert r.removed == pytest.approx(1 - r.keep.mean())
    b = r.bins
    assert b.events.sum() == len(t) and len(b.rate) == len(b.start) == len(t) // 150
    assert np.median(b.rate) == pytest.approx(RATE, rel=0.05)
    assert b.medians.shape == b.bad_signal.shape == (len(b.rate), 6)
    assert b.limits.shape == (6, 2) and b.rate_limits[0] < RATE < b.rate_limits[1]
    assert r.summary == ("To remove: 2.0% of 60,237 events, in 1 stretch(es) of time — "
                         "signal 1,200 (FSC-A, SSC-A).")
    assert b.bad[b.bad_signal.any(axis=1)].all()


def test_events_out_of_order_are_put_in_order():
    _, t, x = acquisition(8)
    x[(t >= 10) & (t <= 11), :2] *= 0.3
    r = check(x, t)
    # written a few events out of order, as happens
    swap = np.arange(len(t))
    swap[100:102] = swap[100:102][::-1]
    jittered = check(x[swap], t[swap])
    assert jittered.intervals == r.intervals and np.array_equal(jittered.keep, r.keep[swap])


def test_files_that_cannot_be_checked_say_so():
    _, t, x = acquisition(9)
    none = check(x, None)
    assert not none.checked and "no Time" in none.summary and none.keep.all()
    assert none.intervals == [] and none.removed == 0.0
    constant = check(x, np.zeros(len(t)))
    assert not constant.checked and "same for every event" in constant.reason
    # two acquisitions in one file: the clock starts again half way
    reset = check(x, np.concatenate([t[::2], t[1::2]]))
    assert not reset.checked and "backwards" in reset.reason
    few = check(x[:2000], t[:2000])
    assert not few.checked and "too few" in few.reason
    broken = t.copy()
    broken[5] = np.nan
    assert not check(x, broken).checked


def test_time_channel_is_found_by_name():
    assert qc.time_channel(["FSC-A", "Time"]) == "Time"
    assert qc.time_channel(["TIME", "FSC-A"]) == "TIME"
    assert qc.time_channel(["FSC-A", "Lifetime"]) is None


# --- the gate ---------------------------------------------------------------------------

def test_time_gate():
    t = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    column = {"Time": t}.__getitem__
    gate = TimeGate(["Time"], [[1.0, 2.0], [4.5, 9.0]])
    assert gate.contains(column).tolist() == [False, True, True, False, False, True]
    gate.inside = False
    assert gate.contains(column).tolist() == [True, False, False, True, True, False]
    assert TimeGate(["Time"], [], inside=False).contains(column).all()

    tree = GateNode("s")
    tree.add(GateNode("Cells", RectangleGate(["FSC-A"], [1.0], [None])))
    before = tree_to_dict(tree)
    node = wrap(tree, "Clean", gate)
    assert [n.path for n in tree.walk()] == ["", "Clean", "Clean/Cells"]
    assert time_node(tree) is node and tree.find("Clean/Cells").parent is node

    back = tree_from_dict(json.loads(json.dumps(tree_to_dict(tree))))
    assert tree_to_dict(back) == tree_to_dict(tree)
    g = back.find("Clean").gate
    assert isinstance(g, TimeGate) and g.intervals == [[1.0, 2.0], [4.5, 9.0]] and not g.inside

    copy = copy_tree(tree)
    copy.find("Clean").gate.intervals[0][0] = 0.0
    assert gate.intervals[0][0] == 1.0

    # another sample's tree, given to this one, keeps this one's intervals
    other = copy_tree(tree)
    other.find("Clean").name = "Good"
    keep_times(other, copy)
    assert other.find("Good").gate.intervals == [[0.0, 2.0], [4.5, 9.0]]
    keep_times(other, GateNode("never checked"))
    assert other.find("Good").gate.intervals == []

    masks = tree.apply({"Time": t, "FSC-A": t}.__getitem__, len(t))
    assert masks["Clean"].sum() == 3 and masks["Clean/Cells"].tolist() == [
        False, False, False, True, True, False]

    unwrap(tree, node)
    assert tree_to_dict(tree) == before and time_node(tree) is None
    assert tree.find("Cells").parent is tree


# --- the session --------------------------------------------------------------------------

TIMESTEP = 0.01


@pytest.fixture
def folder(tmp_path):
    """Three samples of 8 s; the second has an air bubble from 3 to 4 s and
    the third a clog from 5 to 7 s. Time is stored in hundredths of a second."""
    for k, name in enumerate(("A1", "B1", "C1")):
        rng = np.random.default_rng(20 + k)
        t = arrivals(rng, 8.0)
        if name == "C1":
            t = t[~((t >= 5) & (t <= 7)) | (rng.random(len(t)) < 0.03)]
        x = cells(rng, len(t))
        if name == "B1":
            x[(t >= 3) & (t <= 4), :2] *= 0.3
        write_fcs(str(tmp_path / f"{name}.fcs"),
                  np.column_stack([x, t / TIMESTEP]).astype(np.float32), CHANNELS + ["Time"],
                  markers=["", "", "CD8", "CD4", ""],
                  keywords={"SPILL": "2,FITC-A,PE-A,1,0.02,0.01,1", "$TIMESTEP": str(TIMESTEP)})
    return str(tmp_path)


def gated_session(folder) -> Session:
    s = Session()
    s.open_folder(folder)
    a = s.entries[0]
    cells_ = a.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [50000.0], [None])))
    cells_.add(GateNode("CD8+", RectangleGate(["Comp-FITC-A"], [1000.0], [None])))
    s.gates_changed(a, everyone=True)
    s.population = "Cells/CD8+"
    return s


def test_session_finds_each_samples_own_defect(folder):
    s = gated_session(folder)
    a, b, c = (s.quality(e) for e in s.entries)
    assert a.channels == ["FSC-A", "SSC-A", "CD8", "CD4"]
    assert a.scales[0] is None and isinstance(a.scales[2], Biex)
    assert not a.intervals
    (lo, hi), = b.intervals
    assert abs(lo * TIMESTEP - 3) < 3 * BIN and abs(hi * TIMESTEP - 4) < 3 * BIN
    (lo, hi), = c.intervals
    assert abs(lo * TIMESTEP - 5) < 2 * BIN and abs(hi * TIMESTEP - 7) < 2 * BIN
    assert c.counts[qc.RATE] > 0 and b.counts[qc.RATE] == 0


def test_session_applies_keeps_apart_saves_and_removes(folder):
    s = gated_session(folder)
    a, b, c = s.entries
    before = [tree_to_dict(e.gates) for e in s.entries]
    counts = [s.load(e).count("Cells") for e in s.entries]
    found = [(e, s.quality(e).intervals) for e in s.entries]
    s.apply_quality(found)

    # every tree alike, every gate under Clean, and the chosen population followed
    assert s.populations() == ["Clean", "Clean/Cells", "Clean/Cells/CD8+"]
    assert s.population == "Clean/Cells/CD8+"
    assert s.load(a).count("Clean") == a.events
    kept = s.quality(b).keep
    assert s.load(b).count("Clean") == int(kept.sum()) < b.events
    assert np.array_equal(s.load(b).mask("Clean"), kept)
    assert s.load(b).count("Clean/Cells") < counts[1]
    assert s.count(b) == s.load(b).count("Clean/Cells/CD8+")

    # a gate drawn on one sample goes to all; its time intervals do not
    time_node(a.gates).add(GateNode("Singlets", RectangleGate(["SSC-A"], [None], [60000.0])))
    s.gates_changed(a, everyone=True)
    assert "Clean/Singlets" in s.populations()
    assert [time_node(e.gates).gate.intervals for e in s.entries] == [
        [], [list(i) for i in found[1][1]], [list(i) for i in found[2][1]]]
    # and the other way: from a sample that has them to one that has none
    time_node(b.gates).name = "Good"
    s.gates_changed(b, everyone=True)
    assert time_node(a.gates).name == "Good" and time_node(a.gates).gate.intervals == []
    assert time_node(c.gates).gate.intervals == [list(i) for i in found[2][1]]
    assert s.load(a).count("Good") == a.events
    time_node(b.gates).name = "Clean"
    time_node(b.gates).children.pop()        # Singlets
    s.gates_changed(b, everyone=True)

    # it is in the project file, and comes back
    with open(os.path.join(folder, "MIL-C.json"), encoding="utf-8") as f:
        saved = json.load(f)["gates"]["B1.fcs"]["children"][0]["gate"]
    assert saved["type"] == "time" and saved["dims"] == ["Time"] and not saved["inside"]
    again = Session()
    again.open_folder(folder)
    assert [tree_to_dict(e.gates) for e in again.entries] == [tree_to_dict(e.gates)
                                                              for e in s.entries]
    assert again.load(again.entries[1]).count("Clean") == int(kept.sum())

    # applying again to one sample leaves the others' as they are
    s.apply_quality([(b, [])])
    assert s.load(b).count("Clean") == b.events
    assert time_node(c.gates).gate.intervals == [list(i) for i in found[2][1]]

    s.remove_quality()
    assert [tree_to_dict(e.gates) for e in s.entries] == before
    assert s.population == "Cells/CD8+" and s.populations() == ["Cells", "Cells/CD8+"]
    assert [s.load(e).count("Cells") for e in s.entries] == counts


def test_all_events_becomes_the_clean_ones(folder):
    s = Session()
    s.open_folder(folder)
    assert s.population == ""
    s.apply_quality([(e, s.quality(e).intervals) for e in s.entries])
    assert s.population == "Clean" and s.count(s.entries[1]) < s.entries[1].events
    s.remove_quality()
    assert s.population == "" and not any(e.gated for e in s.entries)


def test_no_time_parameter_cannot_be_applied(tmp_path):
    rng = np.random.default_rng(0)
    write_fcs(str(tmp_path / "a.fcs"), cells(rng, 4000).astype(np.float32), CHANNELS)
    s = Session()
    s.open_folder(str(tmp_path))
    r = s.quality(s.entries[0])
    assert not r.checked and "no Time" in r.reason
    with pytest.raises(ValueError, match="no Time"):
        s.apply_quality([(s.entries[0], [])])


# --- the window ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def app(tmp_path_factory):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    QtWidgets = pytest.importorskip("PyQt6.QtWidgets")
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


def test_window(app, folder):
    from milc.ui import worker
    from milc.ui.shell import MainShell

    w = MainShell()
    w.open_folder(folder)
    g, s = w.gates, w.session
    assert g.open_quality() is None               # no sample chosen yet
    g.samples.setCurrentRow(1)
    g.add_rectangle(50000, 0, 262144, 262144, "Cells")
    assert s.populations() == ["Cells"]

    d = g.open_quality()
    assert d.result.checked and len(d.result.intervals) == 1
    assert d.summary.text() == d.result.summary
    assert d.channel.currentText() in ("FSC-A", "SSC-A")
    assert not d.every_button.isEnabled() and not d.remove_button.isEnabled()
    shown = d.series()
    assert shown["time"][-1] == pytest.approx(8.0, abs=0.2)          # seconds
    assert np.median(shown["rate"]) == pytest.approx(RATE, rel=0.1)  # events per second
    (lo, hi), = shown["shaded"]
    assert abs(lo - 3) < 3 * BIN and abs(hi - 4) < 3 * BIN
    assert shown["bad_signal"].any() and not shown["bad_rate"].any()
    d.set_channel(d.result.populations[-1])        # the brighter CD4 population
    in_file_units = d.series()
    assert 4000 < np.median(in_file_units["median"]) < 6000
    assert in_file_units["limits"][0] < np.median(in_file_units["median"]) < in_file_units["limits"][1]
    assert not in_file_units["bad_signal"].any()

    # this sample only: every tree gets the population, only this one loses events
    d.apply_here()
    assert s.populations() == ["Clean", "Clean/Cells"]
    assert g.tree.topLevelItem(0).child(0).text(0) == "Clean"
    assert g.tree.topLevelItem(0).child(0).child(0).text(0) == "Cells"
    kept = int(d.result.keep.sum())
    assert g.tree.topLevelItem(0).child(0).text(1) == f"{kept:,d}"
    assert [bool(time_node(e.gates).gate.intervals) for e in s.entries] == [False, True, False]
    assert d.remove_button.isEnabled() and "1 of 3 samples" in d.state.text()

    # every sample, from the table
    d.survey()
    for _ in range(3):
        worker.wait_all()
        app.processEvents()
    assert d.table.rowCount() == 3 and d.every_button.isEnabled()
    assert [d.table.item(r, 0).text() for r in range(3)] == ["A1", "B1", "C1"]
    assert d.table.item(0, 2).text() == "0.0%" and d.table.item(1, 2).text() != "0.0%"
    assert d.table.item(1, 5).text() == ", ".join(d.result.blamed)
    assert d.ticked() == ["A1", "B1", "C1"]
    d.set_ticked("B1", False)
    d.apply_everyone()
    assert [bool(time_node(e.gates).gate.intervals) for e in s.entries] == [False, False, True]
    assert g.tree.topLevelItem(0).child(0).text(1) == f"{s.entries[1].events:,d}"
    d.set_ticked("B1", True)
    d.apply_everyone()
    assert [bool(time_node(e.gates).gate.intervals) for e in s.entries] == [False, True, True]
    assert "2 of 3 samples" in d.state.text()

    # a gate drawn afterwards, for every sample, leaves each its own intervals
    g.select("Clean/Cells")
    g.add_rectangle(0, 0, 100000, 100000, "Small")
    assert [bool(time_node(e.gates).gate.intervals) for e in s.entries] == [False, True, True]
    assert "Clean/Cells/Small" in s.populations()

    # looking at all events shows time, with what was removed
    g.select("")
    assert g.x_axis.currentData() == "Time" and g.y_axis.currentData() == "#histogram"

    d.remove()
    assert s.populations() == ["Cells", "Cells/Small"] and not d.remove_button.isEnabled()
    assert g.tree.topLevelItem(0).child(0).text(0) == "Cells"

    # deleting the population in the tree is the same as Remove: the gates stay
    d.apply_here()
    g.select("Clean")
    g.delete_selected()
    assert s.populations() == ["Cells", "Cells/Small"]
    d.close()


def test_window_without_time(app, tmp_path):
    from milc.ui.shell import MainShell
    rng = np.random.default_rng(0)
    write_fcs(str(tmp_path / "a.fcs"), cells(rng, 4000).astype(np.float32), CHANNELS)
    w = MainShell()
    w.open_folder(str(tmp_path))
    w.gates.samples.setCurrentRow(0)
    d = w.gates.open_quality()
    assert not d.result.checked and d.summary.text().startswith("Cannot be checked")
    assert not d.here_button.isEnabled()
    d.apply_here()
    assert w.session.populations() == []


# --- the command line, and the report -------------------------------------------------------

def _rows(path) -> dict:
    import csv
    with open(path, encoding="utf-8-sig", newline="") as f:
        return {row["Sample"]: row for row in csv.DictReader(f)}


def test_cli_checks_every_sample(folder, tmp_path, capsys):
    from milc.cli import main
    rng = np.random.default_rng(3)
    few = np.column_stack([cells(rng, 500), np.sort(rng.uniform(0, 800, 500))])
    write_fcs(str(tmp_path / "D1.fcs"), few.astype(np.float32), CHANNELS + ["Time"],
              markers=["", "", "CD8", "CD4", ""], keywords={"$TIMESTEP": str(TIMESTEP)})
    write_fcs(str(tmp_path / "other.fcs"), few[:, :2].astype(np.float32), CHANNELS[:2])
    out, pdf, plots = (str(tmp_path / "o" / n) for n in ("qc.csv", "qc.pdf", "plots"))
    assert main(["qc", folder, "--out", out, "--plots", pdf]) == 0
    said = capsys.readouterr().out
    assert "3 of 5 file(s) checked" in said and "not checked: 500 events are too few" in said

    s = Session()
    s.open_folder(folder)
    rows = _rows(out)
    with open(out, "rb") as f:
        assert f.read(3) == b"\xef\xbb\xbf"
    assert list(rows) == ["A1", "B1", "C1", "D1", "other"]
    for e in s.entries[:3]:            # the table is what the session finds
        r, row = s.quality(e), rows[e.name]
        assert row["Checked"] == "yes" and int(row["Events"]) == e.events
        assert int(row["Removed"]) == int((~r.keep).sum())
        assert float(row["Removed %"]) == pytest.approx(100 * r.removed, abs=0.005)
        assert int(row["Flow rate"]) == r.counts[qc.RATE]
    assert rows["A1"]["Removed"] == "0" and rows["A1"]["Intervals"] == ""
    lo, hi = (float(v) for v in rows["C1"]["Intervals"].split("-"))      # the clog, in seconds
    assert abs(lo - 5) < 2 * BIN and abs(hi - 7) < 2 * BIN and rows["C1"]["Time unit"] == "s"
    assert sorted(rows["B1"]["Parameters"].split("; ")) == ["FSC-A", "SSC-A"]
    assert rows["D1"]["Checked"] == "no" and "too few" in rows["D1"]["Note"]
    assert rows["D1"]["Removed"] == "" and rows["other"]["Note"].startswith("different channels")

    # a page for each sample that could be read, or a file for each
    from test_plate import _pdf_pages
    assert _pdf_pages(pdf) == 4
    assert main(["qc", folder, "--plots", plots, "--plot-signals", "1"]) == 0
    assert sorted(os.listdir(plots)) == [f"{n}_qc.png" for n in ("A1", "B1", "C1", "D1")]
    assert all(os.path.getsize(os.path.join(plots, f)) > 2000 for f in os.listdir(plots))

    # the same options as `check`: the bubble is in the scatter, not in CD8
    assert main(["qc", folder, "--channels", "cd8", "--out", out]) == 0
    assert _rows(out)["B1"]["Removed"] == "0" and int(_rows(out)["C1"]["Removed"]) > 0
    assert main(["qc", folder, "--mad", "3", "--consecutive", "1", "--events-per-bin", "200",
                 "--out", out]) == 0
    r = s.quality(s.entries[1], mad=3, consecutive=1, events_per_bin=200)
    assert int(_rows(out)["B1"]["Removed"]) == int((~r.keep).sum())
    assert main(["qc", folder, "--channels", "nope"]) == 2
    assert main(["qc", str(tmp_path / "o")]) == 2
    assert not os.path.exists(os.path.join(folder, "MIL-C.json"))        # only looked at


def test_cli_cleans_before_gating_as_the_session_does(folder, tmp_path, capsys):
    from milc.cli import main
    from milc.dataset import Dataset
    from milc.fcs import read_fcs
    s = gated_session(folder)                       # Cells/CD8+, saved beside the files
    with open(os.path.join(folder, "MIL-C.json"), encoding="utf-8") as f:
        project = f.read()
    wsp = tmp_path / "none.wsp"
    wsp.write_text("<Workspace/>")
    os.makedirs(tmp_path / "o")
    out, stats_csv = str(tmp_path / "o" / "c"), str(tmp_path / "o" / "s.csv")
    assert main(["concat", folder, "--wsp", str(wsp), "--population", "Cells/CD8+", "--qc",
                 "--out", out]) == 0
    said = capsys.readouterr().out
    assert "population: Clean/Cells/CD8+" in said and "removed in all" in said
    assert main(["stats", folder, "--qc", "--population", "Cells;Cells/CD8+",
                 "--out", stats_csv]) == 0
    assert main(["export", folder, "--qc", "--format", "csv",
                 "--out", str(tmp_path / "o" / "x")]) == 0
    assert "removed in all" in capsys.readouterr().out
    # other limits remove otherwise, and the dataset says which they were
    assert main(["concat", folder, "--wsp", str(wsp), "--qc", "--mad", "3", "--qc-channels",
                 "cd8,FSC-A", "--out", out + "3"]) == 0
    other = Dataset.load(out + "3.npz").provenance["quality"]
    assert other["settings"]["mad"] == 3 and other["settings"]["channels"] == ["FSC-A", "CD8"]
    # and none of it touched the project
    with open(os.path.join(folder, "MIL-C.json"), encoding="utf-8") as f:
        assert f.read() == project

    # the window's route: check every sample, apply, concatenate
    s.apply_quality([(e, s.quality(e).intervals) for e in s.entries])
    assert s.population == "Clean/Cells/CD8+"
    by_session = s.build()
    ds = Dataset.load(out + ".npz")
    assert ds.provenance["population"] == "Clean/Cells/CD8+"
    assert ds.provenance["available"] == by_session.provenance["available"]
    assert ds.provenance["quality"] == by_session.provenance["quality"]
    assert np.array_equal(ds.event_index, by_session.event_index)
    assert np.array_equal(ds.X, by_session.X)
    done = ds.provenance["quality"]
    assert done["settings"] == {"mad": 6.0, "consecutive": 5, "events_per_bin": None,
                                "channels": None}
    b = done["samples"]["B1"]
    assert b["removed"] == int((~s.quality(s.entries[1]).keep).sum()) > 0 and b["checked"]
    assert done["samples"]["A1"]["removed"] == 0 and b["timestep"] == TIMESTEP
    assert done != other
    for line in (f"{b['removed'] / b['events']:6.1%} removed ({b['removed']:,d} events",
                 "  0.0% removed"):
        assert line in said
    # beside the dataset: a table, and a keyword of the concatenated file
    rows = _rows(out + "_qc.csv")
    assert [int(rows[n]["Removed"]) for n in ("A1", "B1", "C1")] == [
        done["samples"][n]["removed"] for n in ("A1", "B1", "C1")]
    assert read_fcs(out + ".fcs").keywords["MILC QC"].split(";")[1] == \
        f"{b['removed']}/{b['events']}"
    # statistics and export were of the clean events, by the workspace's paths
    with open(stats_csv, encoding="utf-8-sig") as f:
        lines = [line.split(",") for line in f.read().splitlines()]
    path, count = lines[0].index("Path"), lines[0].index("Count")
    assert [(row[0], row[path], int(row[count])) for row in lines[1:]] == [
        (e.name, p, s.load(e).count(p)) for e in s.entries
        for p in ("Clean/Cells", "Clean/Cells/CD8+")]
    with open(tmp_path / "o" / "x" / "B1_Clean.csv", encoding="utf-8-sig") as f:
        assert len(f.read().splitlines()) - 1 == s.load(s.entries[1]).count("Clean")

    # the window saved it with the project: that is what the command line
    # uses from now on, and it is not applied twice
    capsys.readouterr()
    assert main(["concat", folder, "--wsp", str(wsp), "--population", "Cells/CD8+", "--qc",
                 "--mad", "2", "--out", out]) == 0
    assert "already applied" in capsys.readouterr().out
    again = Dataset.load(out + ".npz")
    assert again.provenance["quality"] == done and np.array_equal(again.X, ds.X)
    assert main(["concat", folder, "--wsp", str(wsp), "--out", out]) == 0   # and without --qc
    assert Dataset.load(out + ".npz").provenance["quality"] == done
    reopened = Session()
    reopened.open_folder(folder)
    assert reopened.quality_settings == s.quality_settings
    assert reopened.quality_settings["checked"] == ["A1", "B1", "C1"]
    reopened.remove_quality()
    assert reopened.quality_settings is None
    assert "quality" not in reopened.build().provenance


def test_cli_cannot_clean_without_time(tmp_path, capsys):
    from milc.cli import main
    rng = np.random.default_rng(0)
    write_fcs(str(tmp_path / "a.fcs"), cells(rng, 4000).astype(np.float32), CHANNELS)
    assert main(["qc", str(tmp_path)]) == 0
    assert "not checked: the file has no Time parameter" in capsys.readouterr().out
    assert main(["export", str(tmp_path), "--qc", "--out", str(tmp_path / "x")]) == 2
    assert "no Time" in capsys.readouterr().err


def test_report_says_what_was_removed():
    from milc import plate
    from milc.dataset import Dataset
    quality = {"population": "Clean",
               "settings": {"mad": 6.0, "consecutive": 5, "events_per_bin": None,
                            "channels": None},
               "samples": {
                   "a": {"events": 1000, "removed": 0, "intervals": [], "timestep": 0.01,
                         "checked": True},
                   "b": {"events": 2000, "removed": 166, "intervals": [[0.0, 2300.0]],
                         "timestep": 0.01, "checked": True},
                   "c": {"events": 4000, "removed": 40, "intervals": [[5.0, 9.0], [50.0, 51.5]],
                         "timestep": 0.0, "checked": True},
                   "d": {"events": 300, "removed": 0, "intervals": [], "timestep": 0.01,
                         "checked": False},
                   "e": {"events": 1000, "removed": 0, "intervals": [], "timestep": 0.01,
                         "checked": True}}}
    ds = Dataset(np.zeros((4, 1), np.float32), ["A"], [""], [False], np.zeros(4, np.int32),
                 np.arange(4), list("abcde"), [{}] * 5,
                 provenance={"population": "Clean/Cells", "quality": quality,
                             "mode": "equal_per_sample", "seed": 0, "take": {}})
    lines = plate.summary(ds)
    at = lines.index("Quality of the acquisition")
    assert lines[at - 1] == "" and at < lines.index("Sampling")
    assert lines[at + 2] == ("Limits: 6 MADs from the acquisition's median; fewer than 5 good "
                             "bins between two removed stretches go too; bins of 1/500 of the "
                             "file, at least 150 events.")
    assert lines[at + 3] == "Parameters judged: every one but Time."
    assert lines[at + 4:at + 9] == [
        "Samples checked: 4 of 5. Removed in all: 206 of 8,300 events (2.5%).",
        "  b: 8.3% (166 of 2,000 events; 0–23 s)",         # the worst first; in seconds,
        "  c: 1.0% (40 of 4,000 events; 5–9, 50–51.5)",    # or as the file stores Time
        "Nothing removed from 2 sample(s): a, e.",
        "Could not be checked, and kept whole: d."]
    quality["settings"].update(mad=4.0, events_per_bin=300, channels=["CD8", "CD4"])
    said = "\n".join(plate.summary(ds))
    assert "Limits: 4 MADs" in said and "bins of 300 events." in said
    assert "Parameters judged: CD8, CD4." in said
    # applied before the limits were kept: only what the gates say
    quality["settings"] = None
    for v in quality["samples"].values():
        v["checked"] = None
    said = "\n".join(plate.summary(ds))
    assert "The limits it was applied with are not recorded." in said
    assert "Samples checked: 5. Removed" in said and "Could not be checked" not in said
    del ds.provenance["quality"]
    assert "Quality of the acquisition" not in plate.summary(ds)
