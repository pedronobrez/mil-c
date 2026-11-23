"""Batch normalisation beyond the concatenated events: one batch as the
target, and the samples themselves normalised for Gates and Statistics."""

import numpy as np
import pytest

from milc.session import Session
from test_ui import app, batches  # noqa: F401  (fixtures)

PATTERN = r"(?P<Day>\d+)(?P<Group>[A-Za-z]+)(?P<Animal>\d+)"


def _session(folder):
    s = Session()
    s.open_folder(folder)
    s.annotate(PATTERN)
    return s


def _median(s, name, column="Comp-FITC-A"):
    e = next(e for e in s.entries if e.name == name)
    return float(np.median(s.load(e).column(column)))


def _dataset_medians(ds, channel="FITC-A"):
    day, group = ds.obs("Day"), ds.obs("Group")
    x = ds.column(channel)
    return {(d, g): float(np.median(x[(day == d) & (group == g)]))
            for d in ("3", "7") for g in ("Ctrl", "Case")}


def test_one_batch_as_the_target_is_left_as_acquired(batches):
    s = _session(batches)
    raw = _dataset_medians(s.build())
    s.set_normalization("Day", ("Group", "Ctrl"), goal="3")
    assert s.normalization == {"batch": "Day", "reference": ["Group", "Ctrl"],
                               "clusters": 1, "goal": "3"}
    fixed = _dataset_medians(s.build())
    # day 3 is untouched, to the float; day 7, 1.6 times as bright, comes to it
    assert fixed[("3", "Ctrl")] == raw[("3", "Ctrl")]
    assert fixed[("3", "Case")] == raw[("3", "Case")]
    assert raw[("7", "Ctrl")] / raw[("3", "Ctrl")] == pytest.approx(1.6, rel=0.05)
    assert fixed[("7", "Ctrl")] / raw[("3", "Ctrl")] == pytest.approx(1.0, rel=0.03)
    # with the mean as the goal both move, and meet between
    s.set_normalization("Day", ("Group", "Ctrl"))
    assert "goal" not in s.normalization
    mean = _dataset_medians(s.build())
    assert raw[("3", "Ctrl")] < mean[("3", "Ctrl")] < raw[("7", "Ctrl")]
    assert mean[("7", "Ctrl")] / mean[("3", "Ctrl")] == pytest.approx(1.0, rel=0.03)


def test_a_target_without_reference_samples_is_refused(batches):
    s = _session(batches)
    s.set_normalization("Day", ("Group", "Ctrl"), goal="14")
    with pytest.raises(ValueError, match="no reference"):
        s.build()


def test_samples_normalised_for_gates_and_statistics(batches):
    s = _session(batches)
    raw3, raw7 = _median(s, "3Ctrl1"), _median(s, "7Ctrl1")
    assert raw7 / raw3 == pytest.approx(1.6, rel=0.06)
    acquired = _median(s, "7Ctrl1", "FITC-A")

    s.set_normalization("Day", ("Group", "Ctrl"), goal="3", everywhere=True)
    assert s.normalised_everywhere()
    e3 = next(e for e in s.entries if e.name == "3Ctrl1")
    e7 = next(e for e in s.entries if e.name == "7Ctrl1")
    assert s.load(e7).normalised and s.load(e3).normalised
    assert _median(s, "3Ctrl1") == raw3                      # the target batch: as acquired
    assert _median(s, "7Ctrl1") / raw3 == pytest.approx(1.0, rel=0.04)
    assert _median(s, "7Ctrl1", "FITC-A") == acquired        # the uncompensated column is the file's
    assert _median(s, "7Ctrl1", "FSC-A") == pytest.approx(90000, rel=0.02)   # scatter untouched

    # the statistics table reads the same normalised values
    table = s.statistics([""], [("CD8", "Comp-FITC-A")], ["Median"])
    by_name = dict(zip(table.column("Sample"), table.column("CD8 · Median")))
    assert by_name["7Ctrl1"] == pytest.approx(_median(s, "7Ctrl1"))
    assert by_name["7Case1"] / by_name["3Case1"] == pytest.approx(1.0, rel=0.06)

    # concatenating does not normalise a second time: the dataset holds
    # exactly what the samples hold
    ds = s.build()
    done = ds.provenance["normalization"]
    assert done["everywhere"] is True and done["goal"] == "3"
    rows = ds.sample_index == ds.samples.index("7Ctrl1")
    again = np.sort(s.load(e7).column("Comp-FITC-A"))
    assert np.isin(ds.column("FITC-A")[rows], again).all()
    assert np.mean(done["distance_after"]) < np.mean(done["distance_before"]) / 3

    # a setting that changes the curves drops what was normalised under the old ones
    before = _median(s, "7Ctrl1")
    s.set_normalization("Day", ("Group", "Ctrl"), everywhere=True)
    assert _median(s, "7Ctrl1") != before
    # turned off, the samples are the files again; and it survives reopening
    s.set_normalization("Day", ("Group", "Ctrl"), goal="3", everywhere=True)
    again = Session()
    again.open_folder(batches)
    assert again.normalization["everywhere"] is True and again.normalization["goal"] == "3"
    assert _median(again, "7Ctrl1") == pytest.approx(before)
    s.set_normalization(None)
    assert _median(s, "7Ctrl1") == raw7 and not s.load(e7).normalised if hasattr(
        s.load(e7), "normalised") else True


def test_without_a_usable_reference_samples_stay_as_acquired(batches):
    s = _session(batches)
    raw = _median(s, "7Ctrl1")
    s.set_normalization("Day", ("Group", "Nobody"), everywhere=True)
    e = next(e for e in s.entries if e.name == "7Ctrl1")
    assert not s.load(e).normalised and _median(s, "7Ctrl1") == raw


def test_the_window_offers_the_target_and_everywhere(app, batches):
    from milc.ui import worker
    from milc.ui.shell import MainShell

    w = MainShell()
    w.open_folder(batches)
    t = w.samples
    t.fill_from_names(PATTERN)
    assert not t.norm_goal.isEnabled() and not t.norm_everywhere.isEnabled()
    t.set_normalisation("Day", "Group", "Ctrl", goal="3", everywhere=True)
    assert [t.norm_goal.itemData(i) for i in range(t.norm_goal.count())] == ["mean", "3", "7"]
    assert t.norm_goal.currentData() == "3" and t.norm_everywhere.isChecked()
    assert "already show normalised" in t.norm_result.text()
    # the Gates tab says what it is showing
    w.gates.refresh()
    w.gates.select_samples(["7Ctrl1"])
    w.gates.set_axes("FITC-A", "PE-A")
    for _ in range(3):
        worker.wait_all()
        app.processEvents()
    assert "batch-normalised" in w.gates.caption.text()
    # through the controls themselves
    t.norm_everywhere.setChecked(False)
    t.norm_goal.setCurrentIndex(0)
    assert w.session.normalization == {"batch": "Day", "reference": ["Group", "Ctrl"],
                                       "clusters": 1}
    w.gates.select_samples(["7Ctrl1"])
    assert "7Ctrl1" in w.gates.caption.text()
    assert "batch-normalised" not in w.gates.caption.text()
    w.close()
