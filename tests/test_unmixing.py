"""Spectral unmixing on simulated spectral data whose spectra are known:
the spectra read back from simulated controls, OLS and WLS recovering the
abundances, the similarity matrix and complexity index, the CSV, and the
session: unmixed channels as ordinary ones, saved with the project, and
the concatenated events marked stale when the unmixing changes."""

import os

import numpy as np
import pytest

from milc import unmixing as um
from milc.fcs import read_fcs, write_fcs
from milc.gating import GateNode, RectangleGate
from milc.session import Session
from test_ui import app  # noqa: F401  (fixture)

D = 24
DETECTORS = [f"V{i}-A" for i in range(1, 9)] + [f"B{i}-A" for i in range(1, 9)] \
    + [f"R{i}-A" for i in range(1, 9)]
FLUORS = ["BV421", "FITC", "PE", "APC"]
CENTRES = [2.0, 9.0, 13.0, 19.0]


def spectra() -> np.ndarray:
    x = np.arange(D)
    m = np.array([np.exp(-0.5 * ((x - c) / 2.2) ** 2) for c in CENTRES])
    m[2] += 0.3 * np.exp(-0.5 * ((x - 20) / 1.5) ** 2)     # PE seen in the red too
    return m / m.max(axis=1, keepdims=True)


def noisy(clean: np.ndarray, rng, gain: float = 1.0, background: float = 30.0) -> np.ndarray:
    """Poisson-like: variance the signal (times a gain) plus a background."""
    return clean + rng.normal(size=clean.shape) * np.sqrt(gain * np.maximum(clean, 0)
                                                          + background ** 2)


AF_SPECTRUM = np.linspace(1.0, 0.2, D)


def write_spectral(path: str, abundances: np.ndarray, rng, af: float = 0.0) -> None:
    """A spectral file: FSC-A, SSC-A, Time, then the 24 detectors."""
    n = abundances.shape[0]
    clean = abundances @ spectra() + af * AF_SPECTRUM
    det = noisy(clean, rng)
    ev = np.column_stack([rng.normal(80000, 8000, n), rng.normal(30000, 3000, n),
                          np.arange(n, dtype=float), det]).astype(np.float32)
    write_fcs(path, ev, ["FSC-A", "SSC-A", "Time", *DETECTORS],
              ranges=[262144.0] * 3 + [4194304.0] * D)


def single_stain(j: int, n: int, rng) -> np.ndarray:
    a = np.zeros((n, len(FLUORS)))
    positive = rng.random(n) < 0.5
    a[positive, j] = rng.lognormal(10, 0.3, positive.sum())
    return a


# --- spectra and numbers ----------------------------------------------------------


def test_detectors_of_a_spectral_file():
    channels = ["Time", "UV1-A", "SSC-H", "SSC-A", "V1-A", "FSC-H", "FSC-A", "SSC-B-H",
                "SSC-B-A", "B1-A"]
    assert um.detector_channels(channels) == ["UV1-A", "V1-A", "B1-A"]


def test_noise_free_events_are_unmixed_exactly_by_both():
    rng = np.random.default_rng(0)
    a = rng.lognormal(8, 1, size=(500, 4))
    u = um.Unmixing(DETECTORS, FLUORS, spectra())
    y = a @ spectra()
    for method in (um.OLS, um.WLS):
        assert np.allclose(u.unmix(y, method), a, rtol=1e-9, atol=1e-6)
    assert np.allclose(u.residual(y), 0, atol=1e-9)


def test_ols_and_wls_recover_the_abundances():
    rng = np.random.default_rng(1)
    n = 20000
    a = rng.lognormal(8, 1, size=(n, 4)) * (rng.random((n, 4)) < 0.5)
    y = noisy(a @ spectra(), rng)
    u = um.Unmixing(DETECTORS, FLUORS, spectra())
    for method in (um.OLS, um.WLS):
        got = u.unmix(y, method)
        # measured: mean error under 2 units for abundances in the thousands
        assert np.abs((got - a).mean(axis=0)).max() < 5
        for j in range(4):
            on = a[:, j] > 0
            assert np.corrcoef(got[on, j], a[on, j])[0, 1] > 0.99


def test_wls_spreads_a_bright_fluorochrome_less_into_a_neighbour():
    """Cells bright for FITC and without PE: under Poisson-like noise the
    weighted fit leaves PE narrower than ordinary least squares does."""
    rng = np.random.default_rng(2)
    n = 20000
    a = np.zeros((n, 4))
    a[:, 1] = 2e5
    y = noisy(a @ spectra(), rng, gain=5.0)
    u = um.Unmixing(DETECTORS, FLUORS, spectra())
    sd_ols = u.unmix(y, um.OLS)[:, 2].std()
    sd_wls = u.unmix(y, um.WLS)[:, 2].std()
    assert sd_wls < 0.6 * sd_ols                # measured: 0.49
    assert abs(u.unmix(y, um.WLS)[:, 1].mean() - 2e5) < 0.01 * 2e5


def test_spectrum_of_a_single_stain_control():
    rng = np.random.default_rng(3)
    unstained = noisy(np.zeros((5000, D)), rng)
    for j in range(4):
        y = noisy(single_stain(j, 5000, rng) @ spectra(), rng)
        for negative in (None, unstained):
            found = um.spectrum_of(y, negative)
            cos = um.similarity_matrix(np.array([found["spectrum"], spectra()[j]]))[0, 1]
            assert cos > 0.999
            assert found["spectrum"].max() == pytest.approx(1.0)
            assert found["peak"] == int(np.argmax(spectra()[j]))


def test_similarity_and_complexity_index():
    # two spectra [1, 0] and [1, 1]: cosine 1/sqrt(2); singular values
    # sqrt((3 ± sqrt(5)) / 2), whose ratio is the golden ratio squared
    s = np.array([[1.0, 0.0], [1.0, 1.0]])
    assert um.similarity_matrix(s)[0, 1] == pytest.approx(1 / np.sqrt(2))
    assert um.complexity_index(s) == pytest.approx(((1 + np.sqrt(5)) / 2) ** 2)
    assert um.complexity_index(np.eye(3)) == pytest.approx(1.0)
    # scaling a spectrum does not change either: each is taken with its peak at 1
    assert um.complexity_index(s * [[3.0], [0.5]]) == pytest.approx(um.complexity_index(s))
    u = um.Unmixing(DETECTORS, FLUORS, spectra())
    sim = u.similarity()
    assert np.allclose(np.diag(sim), 1) and np.allclose(sim, sim.T)
    # measured for these spectra (FITC and PE overlap most); kept so that a
    # change in either formula shows
    assert sim[1, 2] == pytest.approx(0.4218, abs=1e-4)
    assert u.complexity() == pytest.approx(1.8838, abs=1e-4)


def test_spectra_that_cannot_be_told_apart_are_refused():
    s = spectra()
    with pytest.raises(ValueError):
        um.Unmixing(DETECTORS, ["a", "b", "c", "d", "e"], np.vstack([s, s[0] + s[1]]))
    with pytest.raises(ValueError):
        um.Unmixing(DETECTORS[:3], FLUORS, s[:, :3])


def test_csv_both_ways_round(tmp_path):
    u = um.Unmixing(DETECTORS, FLUORS, spectra())
    path = str(tmp_path / "spectra.csv")
    um.write_csv(path, u)
    with open(path, "rb") as f:
        assert f.read(3) == b"\xef\xbb\xbf"
    dets, names, s = um.read_csv(path)
    assert dets == DETECTORS and names == FLUORS and np.allclose(s, spectra(), atol=1e-5)
    # detectors down the first column, fluorochromes across, out of order
    other = str(tmp_path / "transposed.csv")
    with open(other, "w", encoding="utf-8") as f:
        f.write("Detector," + ",".join(FLUORS) + "\n")
        for d in reversed(range(D)):
            f.write(DETECTORS[d] + "," + ",".join(f"{2 * v:.8g}" for v in spectra()[:, d]) + "\n")
    dets, names, s = um.read_csv(other, DETECTORS)
    assert dets == DETECTORS and names == FLUORS and np.allclose(s, spectra(), atol=1e-6)


# --- from files, into the session ---------------------------------------------------


@pytest.fixture
def spectral(tmp_path):
    """Three samples of a four-colour panel, its four single-stain controls
    and an unstained one, in two folders."""
    rng = np.random.default_rng(7)
    data = tmp_path / "samples"
    controls = tmp_path / "controls"
    data.mkdir()
    controls.mkdir()
    for i in range(3):
        n = 3000
        a = rng.lognormal(7.5, 1.0, size=(n, 4)) * (rng.random((n, 4)) < 0.4)
        write_spectral(str(data / f"S{i + 1}.fcs"), a, rng, af=300.0)
    for j, name in enumerate(FLUORS):
        write_spectral(str(controls / f"{name}.fcs"), single_stain(j, 4000, rng), rng, af=300.0)
    write_spectral(str(controls / "Unstained.fcs"), np.zeros((4000, 4)), rng, af=300.0)
    return str(data), str(controls)


def _unmixing(controls: str, method=um.OLS, af=True) -> um.Unmixing:
    refs = [{"path": os.path.join(controls, f"{n}.fcs"), "fluor": n, "marker": f"CD{j + 1}",
             "negative": "unstained"} for j, n in enumerate(FLUORS)]
    dets, found, noise = um.read_references(refs, unstained=os.path.join(controls, "Unstained.fcs"),
                                            autofluorescence=af)
    assert dets == DETECTORS
    return um.Unmixing.from_references(dets, found, method, noise)


def test_references_from_files(spectral):
    _, controls = spectral
    u = _unmixing(controls)
    assert u.fluors == FLUORS + [um.AF]
    assert u.markers[:4] == ["CD1", "CD2", "CD3", "CD4"]
    for j in range(4):
        assert um.similarity_matrix(np.array([u.spectra[j], spectra()[j]]))[0, 1] > 0.999
    assert um.similarity_matrix(np.array([u.spectra[4], AF_SPECTRUM]))[0, 1] > 0.999


def test_unmixed_channels_are_ordinary_channels(spectral):
    data, controls = spectral
    s = Session()
    s.open_folder(data)
    e = s.usable()[0]
    raw = s.load(e)
    assert raw.fcs.channels[3:] == DETECTORS
    u = _unmixing(controls)
    s.set_unmixing(u)
    sample = s.load(e)
    assert sample.fcs.channels == ["FSC-A", "SSC-A", "Time"] + u.fluors
    assert sample.fcs.markers[3:7] == ["CD1", "CD2", "CD3", "CD4"]
    f = read_fcs(e.path)
    expected = u.unmix(f.events[:, 3:])
    assert np.allclose(sample.column("Comp-PE"), expected[:, 2], rtol=1e-5, atol=1e-2)
    assert np.allclose(sample.column("FSC-A"), f.events[:, 0])
    assert s.fluorescence() == set(u.fluors)
    assert s.scale_of("Comp-APC") is not None and s.scale_of("FSC-A") is None
    labels = dict(s.stat_channels())
    assert labels["CD3"] == "Comp-PE" and labels["FSC-A"] == "FSC-A" and "V1-A" not in labels.values()
    # a gate on an unmixed channel, then the concatenated events
    for x in s.usable():
        x.gates.add(GateNode("PE+", RectangleGate(["Comp-PE"], [1000.0], [None])))
    s.population = "PE+"
    ds = s.build()
    assert ds.channels == ["FSC-A", "SSC-A", "Time"] + u.fluors
    assert ds.compensated == [False] * 3 + [True] * 5
    assert (ds.column("PE") > 1000).all()
    assert s.default_channels() == ["BV421", "FITC", "PE", "APC", "AF"]


def test_unmixing_is_saved_with_the_project(spectral):
    data, controls = spectral
    s = Session()
    s.open_folder(data)
    u = _unmixing(controls, um.WLS)
    s.set_unmixing(u)
    again = Session()
    again.open_folder(data)
    got = again.unmixing
    assert got is not None and got.method == um.WLS and got.fluors == u.fluors
    assert np.allclose(got.spectra, u.spectra) and np.allclose(got.noise, u.noise)
    assert got.markers == u.markers
    e = again.usable()[0]
    assert np.allclose(again.load(e).column("Comp-FITC"), s.load(s.usable()[0]).column("Comp-FITC"))
    again.set_unmixing(None)
    third = Session()
    third.open_folder(data)
    assert third.unmixing is None
    assert third.load(third.usable()[0]).fcs.channels[3:] == DETECTORS


def test_changing_the_unmixing_makes_the_dataset_stale(spectral):
    data, controls = spectral
    s = Session()
    s.open_folder(data)
    for e in s.usable():
        e.gates.add(GateNode("Cells", RectangleGate(["FSC-A"], [1000.0], [None])))
    s.set_unmixing(_unmixing(controls))
    s.population = "Cells"
    s.build()
    assert s.dataset_stale == "" and s.review
    s.clear_review()
    s.set_unmixing(_unmixing(controls, um.WLS))
    assert s.dataset_stale == "the unmixing changed"
    assert s.review == "the unmixing changed"
    s.build()
    assert s.dataset_stale == ""
    s.set_unmixing(None)
    assert s.dataset_stale == "the unmixing changed"


def test_spectra_for_other_detectors_are_refused(spectral, tmp_path):
    data, _ = spectral
    s = Session()
    s.open_folder(data)
    other = um.Unmixing(["X1-A", "X2-A"], ["a"], np.array([[1.0, 0.5]]))
    with pytest.raises(ValueError):
        s.set_unmixing(other)
    assert s.unmixing is None


def test_conventional_files_are_not_touched(tmp_path):
    """Without unmixing nothing changes for a compensated file: its
    channels, its matrix, its values."""
    rng = np.random.default_rng(4)
    ev = np.column_stack([rng.normal(9e4, 9e3, 500), rng.lognormal(6, 1, 500),
                          rng.lognormal(6, 1, 500)]).astype(np.float32)
    write_fcs(str(tmp_path / "a.fcs"), ev, ["FSC-A", "FITC-A", "PE-A"],
              keywords={"SPILL": "2,FITC-A,PE-A,1,0.1,0,1"})
    s = Session()
    s.open_folder(str(tmp_path))
    sample = s.load(s.usable()[0])
    assert s.unmixing is None and not s.unmixes(s.usable()[0])
    assert sample.fcs.channels == ["FSC-A", "FITC-A", "PE-A"]
    assert sample.compensation.name == "acquisition"
    assert s.fluorescence() == {"FITC-A", "PE-A"}
    s.save_project()
    with open(s.project_file(), encoding="utf-8") as f:
        assert '"unmixing": null' in f.read()


# --- the dialog ------------------------------------------------------------------------


def test_the_dialog_computes_shows_and_applies(app, spectral):  # noqa: F811
    from PyQt6 import QtCore

    from milc.ui import worker
    from milc.ui.shell import MainShell
    from milc.ui.unmixing_dialog import UnmixingDialog, apply_to_tab, guess_fluor

    def settle():
        for _ in range(4):
            worker.wait_all()
            app.processEvents()

    assert guess_fluor("x/A3 BUV395 (Beads)_Bead controls_1.fcs") == "BUV395"
    data, controls = spectral
    w = MainShell()
    w.open_folder(data)
    w.show()
    g = w.gates
    w.tabs.setCurrentWidget(g)
    settle()
    assert any(b.text() == "Unmixing…" for b in g.findChildren(type(g.compensation_button)))
    d = UnmixingDialog(w.session, lambda u: apply_to_tab(g, u), g)
    # nothing in the folder of samples looks like a control
    assert d.chosen() == []
    for j, name in enumerate(FLUORS):
        d._add_row(os.path.join(controls, f"{name}.fcs"))
        r = d.controls.rowCount() - 1
        assert d.controls.item(r, 2).text() == name
        d.controls.item(r, 0).setCheckState(QtCore.Qt.CheckState.Checked)
        d.controls.item(r, 3).setText(f"CD{j + 1}")
    d.unstained.addItem("Unstained.fcs", os.path.join(controls, "Unstained.fcs"))
    d.unstained.setCurrentIndex(d.unstained.count() - 1)
    d.autofluorescence.setChecked(True)
    d.compute()
    settle()
    assert d.unmixing is not None and d.unmixing.fluors == FLUORS + ["AF"]
    assert d.similarity.rowCount() == 5 and "Complexity index" in d.complexity.text()
    assert d.check.rowCount() == 4
    # each control is explained by the set: little left over, little spilled
    assert max(c["residual"] for c in d.checks) < 0.05
    assert max(c["leak"] for c in d.checks) < 0.05
    shot = os.environ.get("MILC_SHOT")
    if shot:
        d.show()
        settle()
        for i, tab in enumerate(("spectra", "similarity", "check")):
            d.tabs.setCurrentIndex(i)
            settle()
            d.grab().save(os.path.join(shot, f"unmixing_{tab}.png"))
    d.apply()
    settle()
    assert w.session.unmixing is not None
    labels = [g.x_axis.itemText(i) for i in range(g.x_axis.count())]
    assert "CD3" in labels and "V1-A" not in labels
    assert g.x_axis.findData("Comp-PE") >= 0
    assert d.stop_button.isEnabled()
    if shot:
        w.grab().save(os.path.join(shot, "gates_unmixed.png"))
    assert w.minimumSizeHint().width() <= 1052
    d.stop()
    settle()
    assert w.session.unmixing is None
    assert g.x_axis.findData("V1-A") >= 0
    d.close()
    w.close()
    d.deleteLater()
    w.deleteLater()
    settle()
