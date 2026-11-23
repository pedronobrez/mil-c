"""Proliferation, the cell cycle and kinetics: each model against data
simulated with known answers — the generations and FlowJo's indices of a
dye dilution, the phases of a DNA histogram, the marks of a calcium flux —
and their dialogs driven from the Gates tab without a screen."""

import os

import numpy as np
import pytest

from milc import cellcycle, figures, kinetics, modelfit, proliferation

STYLE = figures.FigureStyle(font="DejaVu Sans")
SHARES = [10, 15, 25, 25, 15, 8, 2]


def _svg(content, tmp_path) -> str:
    path = str(tmp_path / "figure.svg")
    figures.save(content, STYLE, path)
    with open(path, encoding="utf-8") as f:
        return f.read()


# --- proliferation ---------------------------------------------------------------------

def test_indices_are_flowjos_definitions():
    # 100 cells stayed, 100 divided once (200 now), 100 divided twice (400 now)
    found = proliferation.indices([100, 200, 400])
    assert found == pytest.approx({
        "% Divided": 100 * 200 / 300,            # two of the three hundred at the start
        "Division Index": (0 * 100 + 1 * 100 + 2 * 100) / 300,
        "Proliferation Index": (1 * 100 + 2 * 100) / 200,
        "Expansion Index": 700 / 300,
        "Replication Index": 600 / 200})
    none = proliferation.indices([500])
    assert none["% Divided"] == 0 and none["Division Index"] == 0
    assert none["Expansion Index"] == 1 and np.isnan(none["Proliferation Index"])
    assert list(found) == proliferation.INDICES


@pytest.mark.parametrize("sd", [0.05, 0.10, 0.13])       # resolved, touching, run together
@pytest.mark.parametrize("anchored", [True, False])
def test_generations_and_indices_are_recovered(sd, anchored):
    """Within 2.5 points for each generation and for % Divided, 0.05 for
    the Division and Proliferation Index and 4% for the other two — at a
    peak width where the peaks are apart, and where they are not."""
    true = proliferation.indices(SHARES)
    for seed in range(2):
        rng = np.random.default_rng(seed)
        dye = proliferation.simulate(20000, SHARES, 4.5, sd, rng)
        # a control that is a little brighter than the cells that were stimulated
        control = (proliferation.undivided(proliferation.simulate(5000, [1], 4.55, sd, rng))
                   if anchored else None)
        f = proliferation.fit(dye, control)
        assert f.generations == 6 and f.anchored == anchored and f.events == 20000
        assert f.position == pytest.approx(4.5, abs=0.02) and f.sd == pytest.approx(sd, abs=0.01)
        assert np.abs(f.shares[:7] - np.array(SHARES)).max() < 2.5 and f.shares[7:].max() < 0.5
        assert f.cells.sum() == pytest.approx(20000, rel=0.01) and f.modelled > 99
        for name, tolerance in (("% Divided", 2.5), ("Division Index", 0.05),
                                ("Proliferation Index", 0.05),
                                ("Expansion Index", 0.04 * true["Expansion Index"]),
                                ("Replication Index", 0.04 * true["Replication Index"])):
            value, se = f.indices[name]
            assert value == pytest.approx(true[name], abs=tolerance), name
            # the error said is the size of the error made
            assert 0 < se < tolerance and abs(value - true[name]) < 5 * se, name
        assert f.rms_percent < 6 and 0.6 < f.chi2 < 1.6


def test_standard_errors_are_the_scatter_of_repeated_acquisitions():
    values, said = [], []
    for seed in range(12):
        dye = proliferation.simulate(8000, SHARES, 4.5, 0.08, np.random.default_rng(seed))
        f = proliferation.fit(dye, proliferation.Peak(4.5, 0.08))
        values.append([f.indices[n][0] for n in proliferation.INDICES])
        said.append([f.indices[n][1] for n in proliferation.INDICES])
    scatter, said = np.std(values, axis=0, ddof=1), np.mean(said, axis=0)
    assert ((said > 0.5 * scatter) & (said < 2.0 * scatter)).all(), (scatter, said)


def test_a_control_says_which_peak_is_the_undivided_one():
    rng = np.random.default_rng(3)
    shares = [0, 10, 30, 40, 20]                        # every cell divided
    dye = proliferation.simulate(15000, shares, 4.5, 0.07, rng)
    control = proliferation.undivided(proliferation.simulate(6000, [1], 4.5, 0.07, rng))
    assert control.position == pytest.approx(4.5, abs=0.01)
    assert control.sd == pytest.approx(0.07, abs=0.005) and control.events == 6000
    with_control = proliferation.fit(dye, control)
    assert with_control.generations == 4
    assert with_control.indices["% Divided"][0] > 99
    assert with_control.indices["Division Index"][0] == pytest.approx(
        proliferation.indices(shares)["Division Index"], abs=0.05)
    # without one the brightest peak with cells is taken for it: one division is lost
    alone = proliferation.fit(dye)
    assert alone.position == pytest.approx(4.5 - proliferation.HALVING, abs=0.02)
    assert alone.generations == 3 and not alone.anchored


def test_wider_peaks_autofluorescence_and_what_cannot_be_fitted():
    rng = np.random.default_rng(5)
    wider = proliferation.simulate(20000, SHARES, 4.5, 0.05, rng, growth=0.12)
    f = proliferation.fit(wider, grow=True)
    assert f.growth == pytest.approx(0.12, abs=0.04) and f.sd == pytest.approx(0.05, abs=0.006)
    assert np.abs(f.shares[:7] - np.array(SHARES)).max() < 2.5
    assert proliferation.fit(wider).growth == 0.0

    # a dye that halves above a floor: the dim generations crowd against it
    floored = proliferation.simulate(20000, SHARES, 4.0, 0.05, rng, background=300.0)
    f = proliferation.fit(floored, background=300.0)
    assert np.abs(f.shares[:7] - np.array(SHARES)).max() < 2.5
    assert f.indices["Division Index"][0] == pytest.approx(
        proliferation.indices(SHARES)["Division Index"], abs=0.05)

    # zero and negative values have no logarithm: they are not events of the fit
    dye = proliferation.simulate(5000, [50, 50], rng=rng)
    assert proliferation.fit(np.concatenate([dye, [0.0, -5.0, np.nan]])).events == 5000
    with pytest.raises(ValueError, match="too few"):
        proliferation.fit(dye[:20])
    with pytest.raises(ValueError, match="control"):
        proliferation.undivided(-np.ones(100))


def test_proliferation_table_and_figure(tmp_path):
    f = proliferation.fit(proliferation.simulate(8000, SHARES, rng=np.random.default_rng(0)),
                          max_generations=7)
    names, row = proliferation.columns(7), proliferation.row(f, 7)
    assert len(names) == len(row) and names[-1] == "% in generation 7"
    cell = dict(zip(names, row))
    assert cell["Generations"] == 6 and cell["Division Index SE"] == f.indices["Division Index"][1]
    assert cell["Undivided peak"] == pytest.approx(10 ** 4.5, rel=0.05)
    assert sum(row[names.index("% in generation 0"):]) == pytest.approx(100)
    c = proliferation.content(f, "CFSE", "3Case1")
    assert [p.name for p in c.components] == [f"Generation {g}" for g in range(7)]
    assert [n for n, _ in c.category_names()] == [p.name for p in c.components]
    assert c.xticks == [(3.0, "10³"), (4.0, "10⁴")]
    svg = _svg(c, tmp_path)
    assert ">Generation 0 · " in svg and ">CFSE<" in svg and "Division Index" in svg
    assert "Division Index" in proliferation.summary(f)


# --- the cell cycle --------------------------------------------------------------------

#: how far the mean of five simulated acquisitions may be from the truth,
#: in points of a phase, where S is even: Dean-Jett-Fox has no bias to
#: speak of; Watson's is under 1.5 and its mean of five still moves by 1
BIAS = {cellcycle.DJF: 0.8, cellcycle.WATSON: 3.0}
#: and how far any one of them may be
SCATTER = {cellcycle.DJF: 4.0, cellcycle.WATSON: 9.0}


@pytest.mark.parametrize("kind", [cellcycle.WATSON, cellcycle.DJF])
@pytest.mark.parametrize("phases", [(60, 25, 15), (80, 8, 12), (25, 15, 60)])
def test_phases_are_recovered(kind, phases):
    for cv in (0.03, 0.05):
        errors = []
        for seed in range(5):
            dna = cellcycle.simulate(20000, *phases, mean=50000, cv=cv,
                                     rng=np.random.default_rng(seed))
            f = cellcycle.fit(dna, kind)
            errors.append(f.percent - np.array(phases, dtype=float))
            assert f.g1_mean == pytest.approx(50000, rel=0.01)   # also where G2 is the taller
            assert f.g1_cv == pytest.approx(100 * cv, abs=0.5)
            assert f.ratio == pytest.approx(2.0, abs=0.03) and f.events == 20000
            assert f.sub_g1 < 0.5 and f.above_g2 < 0.5
            assert (f.percent_se > 0).all() and (f.percent_se < 5).all()
        errors = np.array(errors)
        assert np.abs(errors.mean(axis=0)).max() < BIAS[kind], (cv, errors.mean(axis=0))
        assert np.abs(errors).max() < SCATTER[kind], (cv, errors)


def test_the_ratio_the_cvs_and_what_lies_outside():
    rng = np.random.default_rng(2)
    dna = cellcycle.simulate(30000, 60, 25, 15, cv=0.04, rng=rng, ratio=1.94, sub_g1=0.10)
    for kind in cellcycle.MODELS:
        f = cellcycle.fit(dna, kind)
        assert f.ratio == pytest.approx(1.94, abs=0.02) and f.sub_g1 == pytest.approx(10, abs=1)
        assert f.fitted[0] < f.g1_mean < f.g2_mean < f.fitted[1]
        assert f.percent.sum() == pytest.approx(100)
        assert np.abs(f.percent - [60, 25, 15]).max() < SCATTER[kind]
        tied = cellcycle.fit(dna, kind, same_cv=True)
        assert tied.g2_cv == pytest.approx(tied.g1_cv, rel=1e-6)
        fixed = cellcycle.fit(dna, kind, ratio=(2.0, 2.0))
        assert fixed.ratio == pytest.approx(2.0, abs=1e-3)
    # told where G1 is, it is taken there even if another peak is taller
    arrested = cellcycle.simulate(20000, 10, 10, 80, mean=40000, cv=0.04, rng=rng)
    assert cellcycle.fit(arrested, g1=41000).g1_mean == pytest.approx(40000, rel=0.02)
    with pytest.raises(ValueError, match="too few"):
        cellcycle.fit(dna[:50])
    with pytest.raises(ValueError, match="unknown model"):
        cellcycle.fit(dna, "fox")


def test_where_each_model_goes_wrong():
    """A fifth of the cells held a quarter of the way through S: a shape
    no parabola has. Dean-Jett-Fox gives several points of G1 to S and
    its chi-square says the fit is poor; Watson follows the histogram."""
    errors = {kind: [] for kind in cellcycle.MODELS}
    chi2 = {kind: [] for kind in cellcycle.MODELS}
    for seed in range(4):
        rng = np.random.default_rng(seed)
        n = 20000
        phase = rng.choice(4, size=n, p=[0.45, 0.20, 0.20, 0.15])
        dna = np.select([phase == 0, phase == 1, phase == 2],
                        [1.0, 1 + rng.random(n), rng.normal(1.25, 0.04, n)], 2.0)
        signal = 50000 * dna * (1 + 0.03 * rng.standard_normal(n))
        for kind in cellcycle.MODELS:
            f = cellcycle.fit(signal, kind)
            errors[kind].append(f.percent - [45, 40, 15])
            chi2[kind].append(f.chi2)
    djf, watson = (np.mean(errors[k], axis=0) for k in (cellcycle.DJF, cellcycle.WATSON))
    assert -9 < djf[0] < -3 and np.mean(chi2[cellcycle.DJF]) > 8
    assert np.abs(watson).max() < 3 and np.mean(chi2[cellcycle.WATSON]) < 2


def test_cell_cycle_table_and_figure(tmp_path):
    dna = cellcycle.simulate(10000, 60, 25, 15, rng=np.random.default_rng(0))
    for kind, name in cellcycle.MODELS.items():
        f = cellcycle.fit(dna, kind)
        row = dict(zip(cellcycle.COLUMNS, cellcycle.row(f), strict=True))
        assert row["% G1"] + row["% S"] + row["% G2"] == pytest.approx(100)
        assert row["G2/G1"] == pytest.approx(row["G2 mean"] / row["G1 mean"])
        c = cellcycle.content(f, "PI", "3Case1")
        assert [p.name for p in c.components] == ["G1", "S", "G2"] and c.fitted == f.fitted
        assert np.allclose(c.model, sum(p.y for p in c.components))
        svg = _svg(c, tmp_path)
        assert ">G1 · " in svg and ">PI<" in svg and name in svg and name in cellcycle.summary(f)


def test_the_uncertainty_of_a_fit_is_that_of_counting():
    # a straight line through counts: the errors of its two numbers, known in closed form
    x = np.arange(50.0)
    def line(t):
        return t[0] + t[1] * x
    theta = np.array([200.0, 3.0])
    cov = modelfit.covariance(line, theta)
    design = np.column_stack([np.ones(50), x])
    bread = np.linalg.inv(design.T @ design)
    exact = bread @ (design.T * line(theta)) @ design @ bread
    assert np.allclose(cov, exact, rtol=1e-5)
    se = modelfit.propagate(lambda t: [t[0] + 10 * t[1], t[0] * t[1]], theta, cov)
    g = np.array([[1.0, 10.0], [3.0, 200.0]])
    assert np.allclose(se, np.sqrt(np.diag(g @ exact @ g.T)), rtol=1e-4)
    assert modelfit.gaussian(np.arange(-50, 51.0), 0.0, 4.0).sum() == pytest.approx(1.0)


# --- kinetics ----------------------------------------------------------------------------

def seconds_of(ticks):
    return ticks * 0.01


def test_slices_statistics_and_time_in_seconds():
    rng = np.random.default_rng(0)
    ticks = np.sort(rng.uniform(0, 10000, 6000))                 # the file's units
    signal = rng.lognormal(5, 0.4, 6000) * (1 + (ticks > 4000))
    blue = rng.lognormal(4, 0.2, 6000)
    c = kinetics.curve(ticks, signal, "s1", step=0.01, width=10.0, smoothing=kinetics.NONE)
    # ten slices of 10 s, from the first event
    assert np.allclose(c.t, seconds_of(ticks).min() + np.arange(5, 100, 10))
    assert c.n.sum() == 6000 and c.group == "s1"
    seconds = ticks * 0.01
    for k in (0, 4, 9):
        inside = (seconds >= ticks.min() * 0.01 + 10 * k) & (
            seconds < ticks.min() * 0.01 + 10 * (k + 1)) if k < 9 else (
            seconds >= ticks.min() * 0.01 + 90)
        assert c.n[k] == inside.sum() and c.y[k] == pytest.approx(np.median(signal[inside]))
    assert np.array_equal(c.smooth, c.y) and (c.se > 0).all()
    # without a $TIMESTEP the axis is the file's own
    raw = kinetics.curve(ticks, signal, width=1000.0, smoothing=kinetics.NONE)
    assert np.allclose(raw.t, c.t * 100, rtol=1e-3) and np.allclose(raw.y, c.y)

    mean = kinetics.curve(ticks, signal, step=0.01, width=10.0, statistic=kinetics.MEAN,
                          span=(0.0, 100.0), smoothing=kinetics.NONE)
    first = seconds < 10
    assert mean.y[0] == pytest.approx(signal[first].mean())
    assert mean.se[0] == pytest.approx(signal[first].std(ddof=1) / np.sqrt(first.sum()))
    above = kinetics.curve(ticks, signal, step=0.01, width=10.0, statistic=kinetics.ABOVE,
                           threshold=200.0, span=(0.0, 100.0), smoothing=kinetics.NONE)
    share = (signal[first] > 200).mean()
    assert above.y[0] == pytest.approx(100 * share)
    assert above.se[0] == pytest.approx(100 * np.sqrt(share * (1 - share) / first.sum()))
    assert above.y[9] > above.y[0] + 30
    ratio = kinetics.curve(ticks, signal, step=0.01, width=10.0, denominator=blue,
                           span=(0.0, 100.0), smoothing=kinetics.NONE)
    assert ratio.y[0] == pytest.approx(np.median(signal[first] / blue[first]))
    # an event with nothing to divide by is not an event of the ratio
    blue[:100] = 0.0
    assert kinetics.curve(ticks, signal, denominator=blue).n.sum() == 5900

    # a slice of a few events has no number, and a default of a hundred slices
    sparse = kinetics.curve(ticks[:150], signal[:150], span=(0.0, 10000.0))
    assert sparse.t.size == kinetics.SLICES and np.isnan(sparse.y[sparse.n < 5]).all()
    with pytest.raises(ValueError, match="same time"):
        kinetics.curve(np.ones(50), signal[:50])
    with pytest.raises(ValueError, match="unknown statistic"):
        kinetics.curve(ticks, signal, statistic="Mode")

    assert kinetics.timestep({"$TIMESTEP": "0.01"}) == 0.01
    assert kinetics.timestep({}) == 0 and kinetics.timestep({"$TIMESTEP": "x"}) == 0
    assert kinetics.time_channel(["FSC-A", " time "]) == " time "
    assert kinetics.time_channel(["FSC-A"]) is None


def test_smoothing_skips_what_is_missing_and_carries_the_error():
    y, se = np.array([1.0, 2.0, 3.0, 4.0, 5.0]), np.ones(5)
    mean, error = kinetics.smoothed(y, se, kinetics.MOVING, 3)
    assert np.allclose(mean, [1.5, 2, 3, 4, 4.5])
    assert np.allclose(error, [np.sqrt(2) / 2, np.sqrt(3) / 3, np.sqrt(3) / 3, np.sqrt(3) / 3,
                               np.sqrt(2) / 2])
    holed, _ = kinetics.smoothed(np.array([1.0, np.nan, 3.0, 5.0]), np.ones(4), kinetics.MOVING, 3)
    assert np.allclose(holed, [1, np.nan, 4, 4], equal_nan=True)
    wide, error = kinetics.smoothed(y, se, kinetics.GAUSSIAN, 1.0)
    w = np.exp(-0.5 * np.arange(-2, 3) ** 2.0)
    assert wide[2] == pytest.approx((w * y).sum() / w.sum())
    assert error[2] == pytest.approx(np.sqrt((w ** 2).sum()) / w.sum())
    assert np.array_equal(kinetics.smoothed(y, se, kinetics.NONE)[0], y)
    assert kinetics.smoothed(y, se, kinetics.MOVING, 4)[0][2] == pytest.approx(3.0)   # 4 → 5 wide


def _flux() -> kinetics.Curve:
    """100 at rest; up by 50 a second from 30.5 to 600 at 40.5; down by
    10 a second to 100 at 90.5."""
    t = np.arange(100) + 0.5
    k = np.arange(100)
    y = np.where(k < 30, 100.0, np.where(k <= 40, 100.0 + 50 * (k - 30),
                                         np.maximum(600.0 - 10 * (k - 40), 100.0)))
    return kinetics.Curve("s", "s", t, np.full(100, 50), y, np.full(100, 2.0), y.copy(),
                          np.full(100, 2.0))


def test_the_marks_of_a_response():
    m = kinetics.measure(_flux(), (0.0, 30.0), (30.0, 100.0))
    assert list(m) == kinetics.COLUMNS
    assert m["Events"] == 5000 and m["Baseline"] == 100 and m["Baseline SD"] == 0
    assert (m["Peak"], m["Peak SE"], m["Peak − baseline"], m["Peak / baseline"]) == (
        600, 2, 500, 6)
    assert (m["Time of peak"], m["Time to peak"]) == (40.5, 10.5)
    assert m["Slope"] == pytest.approx(50) and m["Slope SE"] == pytest.approx(0, abs=1e-9)
    assert m["AUC"] == pytest.approx(0.5 * 10 * 500 + 0.5 * 50 * 500)     # the two triangles
    # the marks decide: a response that ends before the peak has a lower one
    early = kinetics.measure(_flux(), (0.0, 30.0), (30.0, 36.0))
    assert early["Peak"] == 350 and early["Time to peak"] == 5.5
    late = kinetics.measure(_flux(), (0.0, 30.0), (60.0, 100.0))
    assert late["Peak"] == 400 and np.isnan(late["Slope"])               # it only falls
    nothing = kinetics.measure(_flux(), (0.0, 30.0), (200.0, 300.0))
    assert nothing["Baseline"] == 100 and np.isnan(nothing["Peak"]) and np.isnan(nothing["AUC"])
    assert kinetics.default_marks(0.0, 300.0) == ((0.0, 50.0), (50.0, 300.0))


def test_a_simulated_flux_is_measured_and_drawn(tmp_path):
    rng = np.random.default_rng(4)
    seconds = np.sort(rng.uniform(0, 300, 30000))
    rise = np.clip((seconds - 60) / 25, 0, None)
    signal = 400 * (1 + 2.5 * rise * np.exp(1 - rise)) * rng.lognormal(0, 0.35, 30000)
    c = kinetics.curve(seconds, signal, "3Case1", "Case", width=3.0, span=(0.0, 300.0))
    m = kinetics.measure(c, (0.0, 55.0), (60.0, 300.0))
    # the median of the signal is 400 at rest and 3.5 times that 25 s after the stimulus
    assert m["Baseline"] == pytest.approx(400, rel=0.03)
    assert m["Peak / baseline"] == pytest.approx(3.5, rel=0.07)
    assert m["Time to peak"] == pytest.approx(25, abs=6) and m["Slope"] > 20
    # ∫ 400 · 2.5 · u · e^(1 − u) du · 25 s, u from 0 to (300 − 60) / 25
    u = 240 / 25
    area = 400 * 2.5 * 25 * np.e * (1 - (1 + u) * np.exp(-u))
    assert m["AUC"] == pytest.approx(area, rel=0.06)
    content = kinetics.KineticsContent([c, kinetics.curve(seconds, signal / 2, "3Ctrl1", "Ctrl",
                                                         width=3.0)],
                                       "Time (s)", "Median of Indo-1", (0, 55), (60, 300))
    assert content.groups() == ["Case", "Ctrl"]
    svg = _svg(content, tmp_path)
    assert ">Time (s)<" in svg and ">Median of Indo-1<" in svg and ">Ctrl<" in svg


# --- the window -----------------------------------------------------------------------------

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

TRUE = {"3Ctrl1": ([100], (80, 8, 12), 0.2), "3Ctrl2": ([96, 3, 1], (78, 9, 13), 0.2),
        "3Case1": (SHARES, (55, 30, 15), 2.5), "3Case2": ([6, 12, 22, 28, 20, 10, 2],
                                                         (50, 33, 17), 2.0)}


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    QtWidgets = pytest.importorskip("PyQt6.QtWidgets")
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


@pytest.fixture
def gates(app, tmp_path):
    """Four samples, each with a dye that was diluted, a DNA stain and a
    calcium signal that rises 60 s into 300 s of acquisition."""
    from milc.fcs import write_fcs
    from milc.ui import worker
    from milc.ui.shell import MainShell
    rng = np.random.default_rng(1)
    n = 8000
    for name, (shares, phases, flux) in TRUE.items():
        ticks = np.sort(rng.uniform(0, 30000, n))
        rise = np.clip((ticks * 0.01 - 60) / 25, 0, None)
        calcium = 400 * (1 + flux * rise * np.exp(1 - rise)) * rng.lognormal(0, 0.35, n)
        events = np.column_stack([rng.normal(90000, 9000, n),
                                  proliferation.simulate(n, shares, 4.5, 0.07, rng),
                                  cellcycle.simulate(n, *phases, cv=0.04, rng=rng),
                                  calcium, rng.lognormal(5, 0.2, n), ticks]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), events,
                  ["FSC-A", "BL1-A", "YL2-A", "VL1-A", "VL2-A", "Time"],
                  markers=["", "CFSE", "PI", "Indo violet", "Indo blue", ""],
                  keywords={"$TIMESTEP": "0.01"})
    w = MainShell()
    w.open_folder(str(tmp_path))
    w.show()
    w.tabs.setCurrentWidget(w.gates)
    for _ in range(3):
        worker.wait_all()
        app.processEvents()
    w.gates.select_samples(list(TRUE))
    yield w.gates
    w.close()


def _column(dialog, name: str) -> dict:
    return dict(zip(dialog.result.column("Sample"), dialog.result.column(name)))


def test_one_button_opens_the_three_models(gates):
    g = gates
    assert g.model_button.text() == "Model…"
    assert [a.text() for a in g.model_button.menu().actions()] == [
        "Proliferation…", "Cell cycle…", "Kinetics…"]
    assert g.window().minimumSizeHint().width() <= 1280
    from milc.ui import model_dialogs
    for kind, cls in model_dialogs.DIALOGS.items():
        d = g.open_model(kind, show=False)
        assert type(d) is cls and kind in d.windowTitle() and "All events" in d.windowTitle()
        assert d.width() <= 1280 and d.minimumSizeHint().width() <= 1280
        assert d.result.column("Sample") == sorted(TRUE) and not d.problems
        assert d.result.columns[:4] == ["Sample", "Group", "Population", "Path"]
        assert d.compare_button.isEnabled() and d.figure_button.isEnabled()
        assert d.COMPARE in d.value_columns() and "Note" not in d.value_columns()
        d.close()


def test_proliferation_from_gates(gates, tmp_path):
    from milc.ui.figure_dialog import FigureDialog
    g = gates
    g.set_axes("BL1-A", "FSC-A")
    d = g.open_model("Proliferation", show=False)
    assert d.parameter.currentText() == "CFSE" and d.control.currentData() == ""
    assert d.shown == g.entry.name and d.shown in d.said.text()
    assert d.table.rowCount() == 4 and d.table.item(0, 0).text() == "3Case1"
    d.set_control("3Ctrl1")
    assert d.peak.position == pytest.approx(4.5, abs=0.02) and "control" in d.said.text()
    index = _column(d, "Division Index")
    for name, (shares, _, _) in TRUE.items():
        assert index[name] == pytest.approx(proliferation.indices(shares)["Division Index"],
                                            abs=0.06), name
    assert _column(d, "Generations") == {"3Ctrl1": 0, "3Ctrl2": 2, "3Case1": 6, "3Case2": 6}
    assert _column(d, "Group") == {"3Ctrl1": "Ctrl", "3Ctrl2": "Ctrl", "3Case1": "Case",
                                   "3Case2": "Case"}

    # a row of the table is a sample on the plot, and the figure is of that sample
    d.show_sample("3Case2")
    c = d.content()
    assert isinstance(c, proliferation.ProliferationContent) and c.xlabel == "CFSE"
    assert c.summary.startswith("3Case2") and "3Case2" in d.said.text()
    assert len(c.components) == 7 and len(d.plot.listDataItems()) >= 9
    f = FigureDialog(c, d, "proliferation")
    assert f.categories.rowCount() == 7 and not f.colormap.isVisibleTo(f)
    f.close()
    figures.save(c, STYLE, str(tmp_path / "proliferation.pdf"))

    d.generations.setValue(4)
    d.compute()
    assert d.result.columns[-2] == "% in generation 4" and d.fits["3Case1"].generations == 4
    table, found = d.comparison("Division Index", ["Group"])
    assert table.rows[0][1] == "Case vs Ctrl" and 0 < found[0].p < 0.2 and found[0].n == [2, 2]
    path = str(tmp_path / "proliferation.csv")
    d.export(path)
    with open(path, "rb") as f:
        raw = f.read()
    assert raw.startswith(b"\xef\xbb\xbf") and b"Division Index SE" in raw.splitlines()[0]
    assert len(raw.splitlines()) == 5
    d.copy()
    from PyQt6 import QtWidgets
    assert QtWidgets.QApplication.clipboard().text().splitlines()[1].startswith(
        "3Case1\tCase\tAll events\t")
    d.close()


def test_cell_cycle_from_gates(gates, tmp_path):
    g = gates
    g.set_axes("YL2-A", "FSC-A")
    d = g.open_model("Cell cycle", show=False)
    assert d.parameter.currentText() == "PI" and d.model.currentData() == cellcycle.WATSON
    assert "Watson pragmatic" in d.said.text()
    watson = _column(d, "% S")
    d.set_model(cellcycle.DJF)
    assert "Dean-Jett-Fox" in d.said.text()
    djf = _column(d, "% S")
    for name, (_, phases, _) in TRUE.items():
        assert watson[name] == pytest.approx(phases[1], abs=8), name
        assert djf[name] == pytest.approx(phases[1], abs=4), name
        assert _column(d, "G1 mean")[name] == pytest.approx(50000, rel=0.01)
    c = d.content()
    assert isinstance(c, cellcycle.CellCycleContent) and [p.name for p in c.components] == [
        "G1", "S", "G2"]
    figures.save(c, STYLE, str(tmp_path / "cycle.svg"))
    d.ratio_lo.setValue(2.0)
    d.ratio_hi.setValue(2.0)
    d.same_cv.setChecked(True)
    d.compute()
    assert list(_column(d, "G2/G1").values()) == [pytest.approx(2.0, abs=1e-3)] * 4
    table, found = d.comparison("% S", ["Group"])
    assert found[0].groups == ["Case", "Ctrl"] and found[0].mean[0] > found[0].mean[1] + 15
    assert found[0].p < 0.1 and "Case" in table.rows[0][2]

    # a population that is not gated in one of the samples is a row that says so
    g.everyone.setChecked(False)
    g.add_rectangle(60000, 20000, 120000, 130000, "Cells")
    g.select("Cells")
    d2 = g.open_model("Cell cycle", show=False)
    assert "Cells" in d2.windowTitle() and list(d2.fits) == [g.entry.name]
    assert set(d2.problems) == set(TRUE) - {g.entry.name}
    assert "not gated" in d2.problems["3Ctrl1"] and "Note" in d2.result.columns
    assert d2.table.horizontalHeaderItem(d2.table.columnCount() - 1).text() == "Note"
    d2.show_sample("3Ctrl1")
    assert d2.content() is None and "not gated" in d2.said.text() and d2.make_figure() is None
    d.close()
    d2.close()


def test_kinetics_from_gates(gates, tmp_path):
    g = gates
    g.set_axes("Time", "VL1-A")
    d = g.open_model("Kinetics", show=False)
    # not the time against itself: the first parameter that is not Time
    assert d.parameter.currentData() == "FSC-A"
    d.set_parameter("VL1-A")
    assert d._unit == "s" and d._span == pytest.approx((0, 300), abs=1)
    assert d.baseline() == pytest.approx((0, 50), abs=1)
    assert d.response() == pytest.approx((50, 300), abs=1)
    d.set_marks((0, 55), (60, 300))
    fold = _column(d, "Peak / baseline")
    for name, (_, _, flux) in TRUE.items():
        assert fold[name] == pytest.approx(1 + flux, rel=0.12), name
        assert _column(d, "Time to peak")[name] == pytest.approx(25, abs=10)
        assert _column(d, "Baseline")[name] == pytest.approx(400, rel=0.05)
    c = d.content()
    assert isinstance(c, kinetics.KineticsContent) and len(c.series) == 4
    assert c.xlabel == "Time (s)" and c.ylabel == "Median of Indo violet"
    assert (c.baseline, c.response) == ((0, 55), (60, 300))
    figures.save(c, STYLE, str(tmp_path / "kinetics.svg"))

    # dragging a band is typing its ends
    d.regions[1].setRegion((60, 70))
    assert d.response() == pytest.approx((60, 70))
    assert _column(d, "Peak / baseline")["3Case1"] < fold["3Case1"]
    d.set_marks((0, 55), (60, 300))

    # the ratio of two parameters, and the share of events above a threshold
    d.denominator.setCurrentIndex(d.denominator.findData("VL2-A"))
    d.compute()
    assert d.ylabel() == "Median of Indo violet / Indo blue"
    assert _column(d, "Baseline")["3Case1"] == pytest.approx(400 / np.exp(5), rel=0.06)
    d.denominator.setCurrentIndex(0)
    d.statistic.setCurrentIndex(d.statistic.findData(kinetics.ABOVE))
    assert d.threshold.isEnabled()
    d.threshold.setValue(800)
    d.smoothing.setCurrentIndex(d.smoothing.findData(kinetics.GAUSSIAN))
    d.slice.setValue(5.0)
    d.compute()
    assert "above 800" in d.ylabel() and d.fits["3Case1"].t.size == 60
    assert _column(d, "Baseline")["3Case1"] < 5 < 40 < _column(d, "Peak")["3Case1"]
    table, found = d.comparison("Peak − baseline", ["Group"])
    assert found[0].mean[0] > 10 * found[0].mean[1] and found[0].n == [2, 2]
    assert found[0].p < 0.1 and len(table.rows) == 1

    # coloured as the Gates tab colours: by a column, a group is one colour
    g.colour_by.setCurrentIndex(g.colour_by.findData("Group"))
    d2 = g.open_model("Kinetics", show=False)
    assert d2.content().groups() == ["Case", "Ctrl"]
    d.close()
    d2.close()
