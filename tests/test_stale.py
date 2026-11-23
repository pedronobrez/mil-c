"""What is said when the values under the gates and under the concatenated
events change after the fact — the compensation, or normalising the
samples themselves — and the window fitting a laptop's screen."""

import numpy as np

from milc.compensation import Compensation
from milc.gating import GateNode, RectangleGate
from milc.session import Session
from test_ui import app, batches  # noqa: F401  (fixtures)

PATTERN = r"(?P<Day>\d+)(?P<Group>[A-Za-z]+)(?P<Animal>\d+)"


def _session(folder, gate=True):
    s = Session()
    s.open_folder(folder)
    s.annotate(PATTERN)
    if gate:
        for e in s.usable():
            e.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [1000.0], [None])))
    return s


def _other(s):
    m, _ = s.matrix_in_use(s.usable()[0])
    changed = np.array(m.matrix, dtype=float)
    changed[0, 1] += 0.2
    return Compensation(list(m.channels), changed, "edited")


def test_a_new_matrix_makes_the_dataset_stale_and_the_gates_to_review(batches):
    s = _session(batches)
    s.set_matrix(_other(s))
    assert s.dataset_stale == "" and s.review == "the compensation changed"   # nothing concatenated yet
    s.clear_review()
    before = s.build().column("PE-A").copy()
    assert s.dataset_stale == "" and s.review == ""

    s.set_matrix(None)
    assert s.dataset_stale == "the compensation changed" == s.review
    # the review is the project's: it is there when the folder is opened again
    again = Session()
    again.open_folder(batches)
    assert again.review == "the compensation changed" and again.dataset_stale == ""
    # concatenating again answers the first; only looking at the gates answers the second
    after = s.build().column("PE-A")
    assert s.dataset_stale == "" and s.review and not np.array_equal(before, after)
    s.clear_review()
    again.open_folder(batches)
    assert again.review == ""


def test_without_gates_there_is_nothing_to_review(batches):
    s = _session(batches, gate=False)
    s.build()
    s.set_matrix(_other(s))
    assert s.review == "" and s.dataset_stale == "the compensation changed"


def test_a_new_matrix_refits_the_normalisation(batches):
    s = _session(batches)
    s.set_normalization("Day", ("Group", "Ctrl"), goal="3", everywhere=True)
    assert s.review == "the values were batch-normalised"
    s.clear_review()
    e = next(e for e in s.entries if e.name == "7Ctrl1")
    s.load(e)
    assert s._fit is not None
    s.set_matrix(_other(s))
    assert s._fit is None           # fitted on values compensated the old way


def test_normalising_the_samples_asks_for_the_gates_to_be_reviewed(batches):
    s = _session(batches)
    s.build()
    s.set_normalization("Day", ("Group", "Ctrl"))            # the concatenated events only
    assert s.review == "" and s.dataset_stale == ""
    s.set_normalization("Day", ("Group", "Ctrl"), everywhere=True)
    assert s.review == "the values were batch-normalised" == s.dataset_stale
    s.clear_review()
    s.set_normalization("Day", ("Group", "Ctrl"), everywhere=True)   # the same again: no news
    assert s.review == ""
    s.set_normalization("Day", ("Group", "Ctrl"), goal="3", everywhere=True)
    assert s.review == "the batch normalisation changed"
    s.set_normalization("Day", ("Group", "Ctrl"))
    assert s.review == "the batch normalisation was turned off"


def test_the_window_says_so_and_fits_a_laptop(app, batches):
    from milc.ui import worker
    from milc.ui.shell import MainShell

    def settle():
        for _ in range(4):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(batches)
    w.show()
    s, g = w.session, w.gates
    w.samples.fill_from_names(PATTERN)
    for e in s.usable():
        e.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [1000.0], [None])))
    w.tabs.setCurrentWidget(g)
    settle()
    assert not w.stale.isVisible() and not g.review_bar.isVisible()
    w.samples.concatenate()
    settle()

    g.set_compensation(_other(s))
    settle()
    assert w.stale.isVisible() and "before the compensation changed" in w.stale.text()
    assert "Concatenate again" in w.stale.text()
    w.tabs.setCurrentWidget(g)
    settle()
    assert g.review_bar.isVisible() and "before the compensation changed" in g.review.text()
    g.reviewed_button.click()
    assert not g.review_bar.isVisible() and s.review == ""
    w.samples.concatenate()
    settle()
    assert not w.stale.isVisible()

    # every tab on a 1280-point screen
    for tab in (w.samples, w.gates, w.statistics, w.explore, w.clusters):
        w.tabs.setCurrentWidget(tab)
        settle()
        assert tab.minimumSizeHint().width() <= 1280, type(tab).__name__
    assert w.minimumSizeHint().width() <= 1280
    w.close()
