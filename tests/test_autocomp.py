"""A spillover matrix from single-stain controls (`milc.autocomp`): on
controls simulated with a known matrix, both methods find it, the
AutoSpill refinement leaves less spill than its first step, and every
control is given its own detector."""

from __future__ import annotations

import os

import numpy as np
import pytest

from milc import autocomp as ac
from milc.compensation import Compensation
from milc.fcs import write_fcs

DETECTORS = ["FITC-A", "PE-A", "PerCP-Cy5.5-A", "PE-Cy7-A", "APC-A", "BV421-A"]


def truth() -> np.ndarray:
    S = np.eye(6)
    for (i, j), v in {(0, 1): 0.22, (0, 2): 0.05, (1, 0): 0.02, (1, 2): 0.30, (1, 3): 0.04,
                      (2, 3): 0.12, (2, 4): 0.08, (3, 2): 0.03, (4, 2): 0.15, (4, 3): 0.02,
                      (5, 0): 0.03}.items():
        S[i, j] = v
    return S


def write(folder, sims) -> list[str]:
    paths = []
    for name, events, channels in sims:
        p = os.path.join(folder, name + ".fcs")
        write_fcs(p, events, channels)
        paths.append(p)
    return paths


@pytest.fixture(scope="module")
def clean(tmp_path_factory):
    folder = tmp_path_factory.mktemp("controls")
    return write(str(folder), ac.simulate(truth(), DETECTORS, n=10000, seed=1))


def test_both_methods_find_the_known_matrix(clean):
    S = truth()
    # measured on seeds 0-7 of this simulation: classic within 0.0026 of
    # every coefficient, AutoSpill within 0.0004
    for method, tolerance in ((ac.CLASSIC, 0.005), (ac.AUTOSPILL, 0.001)):
        r = ac.compute(clean, DETECTORS, method)
        assert r.channels == DETECTORS and r.converged
        assert np.allclose(np.diag(r.matrix), 1.0)
        assert np.abs(r.matrix - S).max() < tolerance, method
        assert sorted(r.controls) == sorted(DETECTORS)
    # with the unstained tube as the negative of every control
    r = ac.compute(clean, DETECTORS, ac.CLASSIC, universal_negative=True)
    assert np.abs(r.matrix - S).max() < 0.005


def test_the_refinement_lowers_the_spill_left(clean):
    r = ac.compute(clean, DETECTORS, ac.AUTOSPILL)
    start, final = np.abs(r.residual_start).max(), r.largest_residual()
    # measured: 3.1e-4 after the first step, 3.3e-5 after the refinement
    assert final < 0.7 * start and final < 2e-4
    assert r.history[0].scale == "linear" and r.history[-1].scale == "biexponential"
    assert r.history[-1].largest < 1e-4                  # AutoSpill's own tolerance
    # compensated with what was found, the controls show as little spill
    # as compensated with the truth
    controls = ac.load_controls(clean, DETECTORS)
    ac.assign(controls, DETECTORS)
    left_truth, _ = ac.residual_spill(controls, DETECTORS, truth())
    assert final < 2 * np.abs(left_truth).max() + 1e-4
    # and compensating with nothing leaves the spill itself
    none, rows = ac.residual_spill(controls, DETECTORS, np.eye(6))
    assert rows == DETECTORS and abs(none[1, 2] - 0.30) < 0.01


def test_each_control_gets_its_brightest_detector(clean):
    controls = ac.load_controls(clean, DETECTORS)
    # file names say nothing: the order is shuffled and the names hidden
    order = [3, 0, 5, 6, 1, 4, 2]
    shuffled = [controls[i] for i in order]
    for c in shuffled:
        c.path = os.path.join(os.path.dirname(c.path), "tube.fcs")
    shuffled[3].unstained = True                       # the unstained one, as the user says
    got = ac.assign(shuffled, DETECTORS)
    expected = [DETECTORS[i] if i < 6 else None for i in order]
    assert got == expected
    # by name: the unstained control is recognised
    assert ac.looks_unstained("Unstained control.fcs") and ac.looks_unstained("blank_01.fcs")
    assert ac.looks_unstained("Compensation Controls_Unstained Control_001.fcs")
    assert not ac.looks_unstained("PE stained control.fcs")
    assert not ac.looks_unstained("BUV395 Stained Control_017.fcs")


def test_the_user_overrides_the_assignment(clean):
    pe = next(p for p in clean if os.path.basename(p).startswith("PE stained"))
    apc = next(p for p in clean if os.path.basename(p).startswith("APC stained"))
    # PE's control said to be APC's and the other way round: the matrix
    # follows what was said
    r = ac.compute(clean, DETECTORS, ac.AUTOSPILL, primaries={pe: "APC-A", apc: "PE-A"})
    assert r.controls["APC-A"] == pe and r.controls["PE-A"] == apc
    # leaving a control out (said to be unstained) keeps that row the base's
    base = Compensation(DETECTORS, truth() + 0.001 * (1 - np.eye(6)))
    r = ac.compute(clean, DETECTORS, ac.AUTOSPILL, primaries={pe: None}, base=base)
    assert "PE-A" not in r.controls
    assert np.array_equal(r.matrix[1], base.matrix[1])
    assert any("no control for PE-A" in n for n in r.notes)
    with pytest.raises(ValueError, match="two controls"):
        ac.compute(clean, DETECTORS, ac.CLASSIC, primaries={pe: "APC-A"})


def test_a_scatter_gate_keeps_debris_out(tmp_path):
    paths = write(str(tmp_path), ac.simulate(truth(), DETECTORS, n=10000, seed=2, debris=0.2))
    loose = ac.compute(paths, DETECTORS, ac.AUTOSPILL)
    gated = ac.compute(paths, DETECTORS, ac.AUTOSPILL, gate=True)
    controls = ac.load_controls(paths, DETECTORS, gate=True)
    # 20% debris, and the 95% ellipse of the cells: measured 69% to 71% kept
    assert all(0.65 < c.gated < 0.8 for c in controls)
    # measured: 0.0041 off without the gate, 0.0003 with it
    assert np.abs(loose.matrix - truth()).max() > 0.002
    assert np.abs(gated.matrix - truth()).max() < 0.001


def test_huber_is_robust_and_least_squares_is_not():
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 1000, 2000)
    y = 0.25 * x + 10 + rng.normal(0, 5, 2000)
    y[:100] += 3000                                        # 5% far outliers
    a, b, ok = ac.huber(x, y)
    assert ok and abs(b - 0.25) < 0.01
    assert abs(np.polyfit(x, y, 1)[0] - 0.25) > 0.05


def test_files_without_a_detector_and_the_csv(clean, tmp_path):
    with pytest.raises(ValueError, match="has no BV605-A"):
        ac.load_controls(clean[:1], DETECTORS + ["BV605-A"])
    r = ac.compute(clean, DETECTORS, ac.CLASSIC)
    out = tmp_path / "m.csv"
    ac.write_csv(str(out), r)
    raw = out.read_bytes()
    assert raw.startswith(b"\xef\xbb\xbf")
    lines = raw.decode("utf-8-sig").splitlines()
    assert lines[0] == "," + ",".join(DETECTORS) and lines[2].startswith("PE-A,")
    assert ac.fluorescence(["FSC-A", "FSC-H", "SSC-A", "FITC-A", "FITC-H", "PE-A", "Time"]) \
        == ["FITC-A", "PE-A"]


def test_the_command_line(clean, tmp_path, capsys):
    from milc.cli import main
    folder = os.path.dirname(clean[0])
    out = tmp_path / "spill.csv"
    assert main(["autocomp", folder, "--method", "classic", "--out", str(out)]) == 0
    said = capsys.readouterr().out
    assert "Unstained control" in said and "unstained" in said and out.exists()
    assert any(line.startswith("PE stained control") and line.endswith("PE-A")
               for line in said.splitlines())
    # a control said to be another's detector, and a wrong one
    pe = "PE stained control.fcs"
    assert main(["autocomp", folder, "--assign", f"{pe}=BV999-A"]) == 2
    assert "FILE=DETECTOR" in capsys.readouterr().err
