"""The compensation dialog's way to a matrix from single-stain controls:
choose them, see their detectors, compute off the window's thread, put
the result in the table and apply it like a matrix typed by hand."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

QtWidgets = pytest.importorskip("PyQt6.QtWidgets")

from milc import autocomp  # noqa: E402
from milc.fcs import write_fcs  # noqa: E402

DETECTORS = ["FITC-A", "PE-A", "APC-A"]
TRUTH = np.array([[1.0, 0.20, 0.01], [0.03, 1.0, 0.08], [0.0, 0.02, 1.0]])


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


@pytest.fixture
def folders(tmp_path):
    samples, controls = tmp_path / "samples", tmp_path / "controls"
    samples.mkdir()
    controls.mkdir()
    rng = np.random.default_rng(0)
    for name in ("3Ctrl1", "3Ctrl2", "3Case1"):
        ev = np.column_stack([rng.normal(90000, 9000, 300),
                              rng.lognormal(6, 1, (300, 3))]).astype(np.float32)
        write_fcs(str(samples / f"{name}.fcs"), ev, ["FSC-A"] + DETECTORS,
                  keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
    paths = []
    for name, ev, ch in autocomp.simulate(TRUTH, DETECTORS, n=6000, seed=3):
        p = str(controls / f"{name}.fcs")
        write_fcs(p, ev, ch)
        paths.append(p)
    return str(samples), paths


def settle(app):
    from milc.ui import worker
    for _ in range(3):
        worker.wait_all()
        app.processEvents()


def test_a_matrix_from_controls_goes_through_the_table(app, folders):
    from milc.ui.gates_tab import CompensationDialog
    from milc.ui.shell import MainShell

    samples, controls = folders
    w = MainShell()
    w.open_folder(samples)
    w.show()
    g, s = w.gates, w.session
    w.tabs.setCurrentWidget(g)
    settle(app)
    d = CompensationDialog(s, g.entries, g.entry, g.overlaid(), {}, g.set_compensation, g)
    assert d.controls.isEnabled() and d.controls.text() == "From controls…"
    c = d.controls_dialog()
    assert not c.run_button.isEnabled() and not c.use.isEnabled()

    # the controls: each its detector, the unstained one by its name
    c.add_files(controls)
    settle(app)
    names = [c.files.item(r, 0).text() for r in range(c.files.rowCount())]
    assert len(names) == 4 and c.error == ""
    given = {n: c.primary_of(p) for n, p in zip(names, c._listed())}
    assert given == {"FITC stained control": "FITC-A", "PE stained control": "PE-A",
                     "APC stained control": "APC-A", "Unstained control": None}
    assert c.carried.item(0, 1).text() == "5.00"         # the file's matrix beside
    assert c.carried.item(0, 0).text() == "100" and c.carried.item(2, 0).text() == "0"

    # two controls on one detector cannot be computed
    c.set_primary("APC stained control", "PE-A")
    assert not c.run_button.isEnabled() and "PE-A" in c.status.text()
    c.set_primary("APC stained control", "APC-A")
    assert c.run_button.isEnabled()

    c.run()
    settle(app)
    assert c.error == "" and c.result is not None and c.use.isEnabled()
    found = c.matrix()
    assert found.channels == DETECTORS and np.abs(found.matrix - TRUTH).max() < 0.005
    assert c.found.item(0, 1).text() == f"{100 * found.matrix[0, 1]:.1f}"
    assert "iterations" in c.status.text()

    # into the table, not applied; then to every sample, by the usual road
    d.put(found)
    assert d.changed() and "not applied" in d.state.text()
    assert s.current_matrix()[1] == "file"
    assert d.apply_to(None)
    assert [d.source(e.name) for e in g.entries] == ["shared"] * 3
    # four figures of a percentage: within 0.01 of a point
    assert np.abs(s.matrix_of(g.entries[0]).matrix - found.matrix).max() < 1e-4
    # and to one sample alone
    d.reset(None)
    d.put(found)
    d.select([g.entries[1].name])
    assert d.apply_to(d.selected())
    assert [d.source(e.name) for e in g.entries] == ["file", "edited", "file"]
    c.close()
    d.close()
    w.close()


def test_a_control_that_cannot_be_read_is_said(app, folders, tmp_path):
    from milc.ui.autocomp_dialog import ControlsDialog
    _, controls = folders
    bad = tmp_path / "other.fcs"
    write_fcs(str(bad), np.zeros((10, 2), np.float32), ["FSC-A", "X-A"])
    c = ControlsDialog(DETECTORS)
    c.add_files(controls[:2] + [str(bad)])
    settle(app)
    assert "other.fcs has no FITC-A, PE-A, APC-A" in c.status.text()
    assert c.files.rowCount() == 2 and len(c.paths) == 2
    c.close()
