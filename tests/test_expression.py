"""Differential expression within the clusters: the moderated t against
limma's own numbers, what it does with missing medians, how often it
cries wolf, and the dialog that asks it."""

import re

import numpy as np
import pytest

from milc import cluster_stats, clusters, compare, expression as ex, figures
from milc.abundance import design
from milc.expression_figures import DOWN, REST, UP, DifferenceHeatmap, ExpressionVolcano
from test_ui import app, folder  # noqa: F401  (fixtures)

# --- against limma ---------------------------------------------------------------------
# R 4.6.1, limma 3.68.5: eBayes(lmFit(y, X, weights)) and its column of the
# group; `legacy=TRUE` where there are NAs. Four samples of "a", four of "b".

Y = np.array([[2.4, 2.43, 2.373, 2.311, 1.955, 1.901, 2.006, 2.134],
              [1.902, 1.876, 2.098, 2.071, 2.021, 1.814, 1.994, 2.139],
              [1.597, 1.863, 1.43, 1.613, 1.447, 1.929, 1.62, 2.081],
              [1.424, 1.372, 1.022, 1.319, 1.993, 2.017, 1.77, 1.928],
              [1.511, 1.596, 2.53, 1.596, 1.984, 2.442, 1.708, 1.944],
              [2.028, 2.016, 1.694, 2.019, 2.34, 1.613, 2.215, 2.03],
              [1.923, 2.24, 2.091, 1.856, 2.009, 2.069, 1.977, 2.082],
              [1.973, 2.267, 2.575, 1.73, 2.081, 1.815, 2.051, 1.525]])
W = np.array([[0.95, 1.39, 0.56, 2.91, 0.88, 2.11, 1.11, 2.66],
              [2.09, 0.66, 2.58, 2.85, 2.74, 1.84, 0.69, 0.82],
              [2.81, 1.79, 0.79, 2.69, 2.03, 1.84, 1.32, 1.41],
              [0.95, 0.4, 2.67, 1.56, 1.78, 1.17, 2.33, 0.37],
              [1.3, 0.38, 0.63, 2.91, 2.08, 1.46, 1.71, 2.66],
              [1.23, 1.89, 2.15, 1.26, 1.7, 2.37, 2.75, 0.71],
              [2.82, 0.31, 2.33, 2.49, 0.67, 1.43, 2.5, 0.34],
              [2.0, 2.44, 1.69, 2.26, 0.91, 0.84, 1.28, 0.78]])
GROUPS = ["a"] * 4 + ["b"] * 4
COEF = [0.37950, -0.00525, -0.14350, -0.64275, -0.21125, -0.11025, -0.00675, 0.26825]


def _missing() -> np.ndarray:
    y = Y.copy()
    y[1, 0] = y[4, 5] = y[4, 6] = y[6, 3] = np.nan
    return y


LIMMA = {
    "plain": dict(
        d0=3.30964214359, coef=COEF,
        t=[4.44353498679, -0.05178701667, -0.93061474123, -5.75051300681, -0.87706413129,
           -0.68575551034, -0.06592123577, 1.38456910137],
        p=[0.0014865012807, 0.9597926289776, 0.3755546279678, 0.0002429709593, 0.4025172292367,
           0.5095825308553, 0.9488345421368, 0.1984665896615]),
    "weights": dict(
        d0=3.91035306186, weights=W,
        coef=[0.34304413682, 0.04931533316, -0.10649423192, -0.70318771846, -0.30786290209,
              -0.12348196189, -0.05158454409, 0.21585749815],
        t=[3.6411106825, 0.5028680423, -0.7407755211, -5.6671142300, -1.6054047449,
           -0.7219147656, -0.5215022905, 1.0234210168],
        p=[0.0045958047177, 0.6260432855066, 0.4760146169693, 0.0002145906794, 0.1397615220426,
           0.4870306856197, 0.6134710402876, 0.3304409755849]),
    "batch": dict(
        d0=2.3989377206, covariate=list("11221122"), coef=COEF,
        t=[4.72573615935, -0.06586554498, -0.86884698774, -7.10573971679, -0.79890487284,
           -0.63714852009, -0.06697318741, 1.27720128811],
        p=[0.0018444040436, 0.9492294809280, 0.4122254706595, 0.0001481157587, 0.4492270744803,
           0.5432333476424, 0.9483771229152, 0.2401485487215]),
    "paired": dict(
        d0=2.53714796599, covariate=list("12341234"), coef=COEF,
        t=[4.23681378166, -0.06528456924, -1.21391427026, -8.46952944044, -0.76987559984,
           -0.67402085286, -0.06925130739, 2.11993176089],
        p=[0.0065340084258, 0.9502399112165, 0.2740037554491, 0.0002255881856, 0.4729257916278,
           0.5274138662832, 0.9472219875863, 0.0821421706584]),
    "missing": dict(
        d0=3.09726780756, y=_missing(),
        coef=[0.3795, 0.023, -0.1435, -0.64275, -0.15575, -0.11025, 0.05041666667, 0.26825],
        t=[4.5430485943, 0.2060647595, -0.9292721261, -5.7958467884, -0.5391779700,
           -0.6842076044, 0.4896661142, 1.3774326827],
        p=[0.0013618533626, 0.8418271982901, 0.3767468987414, 0.0002503361634, 0.6062566298390,
           0.5108958939510, 0.6373600424787, 0.2013175010816],
        df=[9.097267808, 8.097267808, 9.097267808, 9.097267808, 7.097267808, 9.097267808,
            8.097267808, 9.097267808]),
    "missing, weights": dict(
        d0=3.723846646, y=_missing(), weights=W,
        coef=[0.343044136818, 0.090018062397, -0.106494231923, -0.703187718464,
              -0.273997187060, -0.123481961892, -0.002507904376, 0.215857498146],
        t=[3.76791407265, 0.89251758638, -0.74539291572, -5.77072741467, -1.44994387010,
           -0.72334686821, -0.02485090424, 1.02399921057],
        p=[0.0038602676344, 0.3960669314009, 0.4736664046495, 0.0002004291047, 0.1864349073763,
           0.4865022699994, 0.9807330192285, 0.3306350982863]),
}


@pytest.mark.parametrize("case", list(LIMMA))
def test_the_moderated_t_is_limmas(case):
    ref = LIMMA[case]
    X = design(GROUPS, ("a", "b"), ref.get("covariate"))
    f = ex.fit(ref.get("y", Y), X, ref.get("weights"))
    t, df, p, d0 = ex.moderated(f)
    assert d0 == pytest.approx(ref["d0"], rel=1e-8)
    assert f.coefficient == pytest.approx(ref["coef"], rel=1e-8)
    assert t == pytest.approx(ref["t"], rel=1e-8)
    assert p == pytest.approx(ref["p"], rel=1e-7)
    if "df" in ref:
        assert df == pytest.approx(ref["df"], rel=1e-8)


def test_unmoderated_is_the_t_test_and_paired_is_the_paired_t():
    from scipy import stats as st
    t, df, p, d0 = ex.moderated(ex.fit(Y, design(GROUPS, ("a", "b"))), moderate=False)
    ref = st.ttest_ind(Y[:, :4], Y[:, 4:], axis=1)
    assert d0 == 0 and (df == 6).all()
    assert t == pytest.approx(ref.statistic) and p == pytest.approx(ref.pvalue)
    paired = ex.fit(Y, design(GROUPS, ("a", "b"), list("12341234")))
    t, df, p, _ = ex.moderated(paired, moderate=False)
    ref = st.ttest_rel(Y[:, :4], Y[:, 4:], axis=1)
    assert (df == 3).all()
    assert t == pytest.approx(ref.statistic) and p == pytest.approx(ref.pvalue)


# --- medians ---------------------------------------------------------------------------

def test_medians_per_sample_and_cluster():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(500, 3))
    labels = rng.choice([2, 5, 9], 500)
    sample = rng.integers(0, 4, 500)
    labels[(labels == 9) & (sample == 3)] = 2     # a cluster a sample has no cell of
    ids, values, counts = ex.medians(X, labels, sample, 4)
    assert ids == [2, 5, 9] and values.shape == (3, 4, 3) and counts.sum() == 500
    for i, c in enumerate(ids):
        for s in range(4):
            mask = (labels == c) & (sample == s)
            assert counts[i, s] == mask.sum()
            if mask.any():
                assert values[i, s] == pytest.approx(np.median(X[mask], axis=0))
    assert counts[2, 3] == 0 and np.isnan(values[2, 3]).all()


# --- the comparison --------------------------------------------------------------------

def simulate(rng, n=(5, 5), clusters=12, markers=6, shift=0.0, shifted=6):
    """(values clusters x samples x markers, counts, groups, the rows —
    cluster x marker — that were shifted). Each row varies between
    animals with a standard deviation of its own; a median also has the
    sampling error of the cells it is of."""
    samples = sum(n)
    base = np.exp(rng.uniform(np.log(20), np.log(5000), clusters))
    depth = np.exp(rng.uniform(-0.8, 0.8, samples))
    counts = rng.poisson(base[:, None] * depth[None, :]).astype(float)
    sd = np.exp(rng.uniform(np.log(0.05), np.log(0.4), (clusters, 1, markers)))
    values = 2.0 + rng.normal(size=(clusters, samples, markers)) * sd
    with np.errstate(divide="ignore"):
        values += rng.normal(size=values.shape) * (1.0 / np.sqrt(counts))[:, :, None]
    truth = np.zeros((clusters, markers), dtype=bool)
    if shift:
        truth.flat[rng.choice(truth.size, shifted, replace=False)] = True
        values[:, :n[0], :] += shift * truth[:, None, :]
    return values, counts, ["a"] * n[0] + ["b"] * n[1], truth.ravel()


def tally(reps, n=(5, 5), shift=0.0, seed=0, **kw):
    """Over `reps` simulated data sets: the share of true nulls with
    p < 0.05 and p < 0.01, the share of the shifted rows found at
    FDR < 0.05, the false discovery rate had (the mean, over the data
    sets, of the share of what was found that was false), the share of
    rows tested."""
    rng = np.random.default_rng(seed)
    p05 = p01 = nulls = hits = truths = tested = rows = 0
    fdp = []
    for _ in range(reps):
        values, counts, groups, truth = simulate(rng, n, shift=shift)
        res = ex.differential(values, counts, [str(i) for i in range(values.shape[0])],
                              [str(i) for i in range(values.shape[2])], groups, ("a", "b"),
                              **kw)
        p = np.array([r.p for r in res])
        q = np.array([r.p_adjusted for r in res])
        ok = np.isfinite(p)
        tested, rows = tested + ok.sum(), rows + ok.size
        nulls += (ok & ~truth).sum()
        p05 += (p[ok & ~truth] < 0.05).sum()
        p01 += (p[ok & ~truth] < 0.01).sum()
        truths += truth.sum()
        hits += (q[truth] < 0.05).sum()
        fdp.append((q[~truth] < 0.05).sum() / max((q < 0.05).sum(), 1))
    return dict(p05=p05 / nulls, p01=p01 / nulls, power=hits / max(truths, 1),
                fdr=float(np.mean(fdp)), tested=tested / rows)


def test_under_no_difference_p_is_small_as_often_as_it_says():
    null = tally(150, (5, 5), seed=1)                # about 10,000 tests
    assert 0.04 < null["p05"] < 0.06, null
    assert 0.005 < null["p01"] < 0.016, null
    assert 0.035 < tally(150, (3, 3), seed=3)["p05"] < 0.065
    # weighting by cells where animals differ more than cells do: too many
    weighted = tally(150, (5, 5), seed=1, weights=True)
    assert null["p05"] < weighted["p05"] < 0.085, weighted


def test_a_shift_is_found_and_moderation_is_what_finds_it():
    found = tally(60, (5, 5), shift=0.5, seed=2)
    plain = tally(60, (5, 5), shift=0.5, seed=2, moderate=False)
    assert found["power"] > 0.55 and found["fdr"] < 0.10
    assert found["power"] > plain["power"] + 0.05


def _small():
    rng = np.random.default_rng(5)
    values, counts, groups, _ = simulate(rng, (4, 4), clusters=3, markers=2)
    counts[:] = 500
    return values, counts, groups


def test_too_few_cells_is_a_missing_median_and_too_few_samples_is_said():
    values, counts, groups = _small()
    values[0, :4, 0] += 2.0
    counts[1, 0] = 4                      # one sample short: still tested
    counts[2, :3] = 4                     # three of group a short: not tested
    res = ex.differential(values, counts, ["A", "B", "C"], ["m1", "m2"], groups, ("a", "b"))
    assert [(r.cluster, r.marker) for r in res] == [
        ("A", "m1"), ("A", "m2"), ("B", "m1"), ("B", "m2"), ("C", "m1"), ("C", "m2")]
    a, b, c = res[0], res[2], res[4]
    assert a.n == (4, 4) and a.difference > 1.5 and a.p_adjusted < 0.05 and a.note == ""
    assert a.difference == pytest.approx(a.medians[0] - a.medians[1])
    assert b.n == (3, 4) and np.isfinite(b.p) and b.df < a.df
    assert c.n == (1, 4) and np.isnan(c.p) and np.isnan(c.p_adjusted)
    assert "fewer than 2 samples of a with 10 cells" in c.note
    # reported, with what there is to say: the plain difference of the means
    assert c.difference == pytest.approx(c.medians[0] - c.medians[1])
    # the correction counts the tests made, not the rows
    p = [r.p for r in res[:4]]
    assert [r.p_adjusted for r in res[:4]] == pytest.approx(list(compare.adjust(p)))
    # a lower bar lets them in
    low = ex.differential(values, counts, ["A", "B", "C"], ["m1", "m2"], groups, ("a", "b"),
                          min_cells=3)
    assert low[4].n == (4, 4) and np.isfinite(low[4].p)
    none = ex.differential(values, counts, ["A", "B", "C"], ["m1", "m2"], groups, ("a", "b"),
                           correction=compare.NONE)
    assert [r.p_adjusted for r in none[:4]] == pytest.approx(p)


def test_other_groups_are_left_out_and_one_group_is_no_comparison():
    values, counts, groups = _small()
    three = groups[:4] + ["b", "b", "c", "c"]
    res = ex.differential(values, counts, ["A", "B", "C"], ["m1", "m2"], three, ("a", "b"))
    assert all(r.n == (4, 2) for r in res)
    same = ex.differential(values[:, [0, 1, 2, 3, 4, 5]], counts[:, :6], ["A", "B", "C"],
                           ["m1", "m2"], three[:6], ("a", "b"))
    assert [r.p for r in res] == pytest.approx([r.p for r in same])
    alone = ex.differential(values, counts, ["A", "B", "C"], ["m1", "m2"], groups, ("a", "z"))
    assert all(np.isnan(r.p) and "one group" in r.note for r in alone)


def test_a_covariate_takes_its_effect_out_and_cannot_be_the_groups():
    values, counts, groups = _small()
    day = list("12121212")
    values[:, [0, 2, 4, 6], :] += 3.0                 # a day far larger than any group
    values[0, :4, 0] += 0.8
    plain = ex.differential(values, counts, ["A", "B", "C"], ["m1", "m2"], groups, ("a", "b"))
    adjusted = ex.differential(values, counts, ["A", "B", "C"], ["m1", "m2"], groups,
                               ("a", "b"), covariate=day)
    assert adjusted[0].p < plain[0].p and adjusted[0].p_adjusted < 0.05
    assert adjusted[0].difference == pytest.approx(plain[0].difference)
    aliased = ex.differential(values, counts, ["A", "B", "C"], ["m1", "m2"], groups,
                              ("a", "b"), covariate=groups)
    assert all(np.isnan(r.p) and "cannot be told apart" in r.note for r in aliased)
    each = ex.differential(values, counts, ["A", "B", "C"], ["m1", "m2"], groups,
                           ("a", "b"), covariate=list("12345678"))
    assert all(np.isnan(r.p) and "cannot be told apart" in r.note for r in each)
    few = ex.differential(values[:, [0, 1, 4, 5]], counts[:, [0, 1, 4, 5]], ["A", "B", "C"],
                          ["m1", "m2"], list("aabb"), ("a", "b"), covariate=list("1212"))
    assert all(np.isfinite(r.p) for r in few)       # one degree of freedom, and the others'
    two = ex.differential(values[:, [0, 4]], counts[:, [0, 4]], ["A", "B", "C"],
                          ["m1", "m2"], list("ab"), ("a", "b"), min_samples=1)
    assert all(np.isnan(r.p) and "too few samples" in r.note for r in two)


def test_within_panels_and_the_table(tmp_path):
    rng = np.random.default_rng(8)
    values, counts, groups, _ = simulate(rng, (6, 6), clusters=3, markers=2)
    counts[:] = 300
    groups = ["a", "a", "a", "b", "b", "b"] * 2
    days = ["3"] * 6 + ["14"] * 6
    values[1, 6:9, 1] += 3.0                           # only on day 14
    res = ex.within(values, counts, ["A", "B", "C"], ["m1", "m2"], groups, ("a", "b"), days)
    assert [r.panel for r in res] == ["3"] * 6 + ["14"] * 6
    assert res[9].p_adjusted < 0.05 and res[3].p_adjusted > 0.05
    assert [r.p_adjusted for r in res] == pytest.approx(
        list(compare.adjust([r.p for r in res])))
    table = ex.to_table(res, panel="Day")
    assert table.columns[:3] == ["Day", "Cluster", "Marker"]
    assert table.columns[3:5] == ["a (median)", "b (median)"]
    assert "p adjusted (Benjamini-Hochberg)" in table.columns and len(table.rows) == 12
    assert table.column("n")[0] == "3 vs 3"
    path = tmp_path / "ds.csv"
    table.to_csv(str(path))
    assert path.read_bytes().startswith(b"\xef\xbb\xbf")
    assert ex.to_table([]).rows == []


# --- the figures -----------------------------------------------------------------------

def results(panels=False):
    rng = np.random.default_rng(11)
    values, counts, groups, _ = simulate(rng, (6, 6), clusters=4, markers=3)
    counts[:] = 400
    counts[3, :5] = 2                                  # a cluster Case has almost none of
    values[0, :6, 1] += 2.0
    values[2, :6, 0] -= 2.0
    groups = ["Case" if g == "a" else "Ctrl" for g in groups]
    args = (values, counts, ["T", "B", "NK", "rare"], ["CD8 FITC-A", "CD4 PE-A", "TIGIT"],
            groups, ("Case", "Ctrl"))
    if panels:
        return ex.within(*args, ["3", "14"] * 6, min_samples=1)
    return ex.differential(*args)


def test_the_volcano_and_the_heatmap_of_differences(tmp_path):
    found = results()
    volcano = ExpressionVolcano(found)
    pts = volcano.points()
    assert len(pts) == 9                               # the untested cluster is no point
    kinds = {name: kind for _, _, name, kind in pts}
    assert kinds["T · CD4"] == UP and kinds["NK · CD8"] == DOWN and kinds["B · TIGIT"] == REST
    assert volcano.xlabel == "Difference in median (Case − Ctrl)"
    assert volcano.names()[UP] == "Higher in Case" and volcano.fold == 0
    assert ExpressionVolcano([]).limits()[0] > 0

    heat = DifferenceHeatmap(found)
    assert heat.panels() == [None] and heat.clusters() == ["T", "B", "NK", "rare"]
    difference, p = heat.cells()
    assert difference.shape == (4, 3) and np.isnan(difference[3]).all() and np.isnan(p[3]).all()
    assert difference[0, 1] == found[1].difference and p[0, 1] == found[1].p_adjusted
    assert heat.limit() == pytest.approx(max(abs(r.difference) for r in found[:9]))
    assert isinstance(heat, figures.HeatmapContent) and not figures.uses_colormap(heat)
    assert heat != DifferenceHeatmap(found) and DifferenceHeatmap([]).limit() == 1.0

    style = figures.FigureStyle(font="DejaVu Sans")
    for content, said in ((volcano, ["Higher in Case", "T · CD4"]),
                          (heat, ["Difference in median (Case − Ctrl)", "rare", "TIGIT"])):
        assert figures.preview_png(content, style, dpi=40)[:4] == b"\x89PNG"
        path = str(tmp_path / "figure.svg")
        figures.save(content, style, path)
        with open(path, encoding="utf-8") as f:
            svg = f.read()
        assert all(text in svg for text in said)
    starred = re.compile(r">\*+</text>")               # the stars, and none without `annotate`
    assert starred.search(svg)
    figures.save(heat, figures.FigureStyle(font="DejaVu Sans", annotate=False), path)
    with open(path, encoding="utf-8") as f:
        assert not starred.search(f.read())

    panelled = DifferenceHeatmap(results(panels=True))
    assert panelled.panels() == ["3", "14"] and panelled.cells("14")[0].shape == (4, 3)
    assert figures.preview_png(panelled, style, dpi=40)[:4] == b"\x89PNG"
    assert len(ExpressionVolcano(results(panels=True)).panels()) == 2

    # both go on a plate with the other figures
    from milc import plate
    page = plate.Plate(columns=2, width_mm=180, style=style)
    page.add(panelled, title="Differences")
    page.add(volcano, title="Volcano")
    plate.save(page, path)
    with open(path, encoding="utf-8") as f:
        svg = f.read()
    assert "Difference in median (Case − Ctrl)" in svg and "T · CD4" in svg


# --- the window ------------------------------------------------------------------------

def _medians(seed=12):
    rng = np.random.default_rng(seed)
    values, counts, groups, _ = simulate(rng, (6, 6), clusters=4, markers=3)
    counts[:] = 400
    counts[3, :5] = 2
    values[0, :6, 1] += 2.0
    groups = ["Case" if g == "a" else "Ctrl" for g in groups]
    groups[5], groups[11] = "UV", "UV"
    return ex.Medians(["T", "B", "NK", "rare"], [f"s{i}" for i in range(12)],
                      ["CD8 FITC-A", "CD4 PE-A", "TIGIT"], values, counts, "biex"), groups


def test_dialog(app, tmp_path):
    from PyQt6 import QtWidgets
    from milc.ui.expression_dialog import ExpressionDialog
    from milc.ui.figure_dialog import FigureDialog
    m, groups = _medians()
    day = ["3", "7"] * 6
    d = ExpressionDialog(m, groups, "Group", {"Day": day})
    assert [d.first.itemData(i) for i in range(d.first.count())] == ["Case", "Ctrl", "UV"]
    assert (d.first.currentData(), d.second.currentData()) == ("Case", "Ctrl")
    assert d.chosen_markers() == m.markers and d.markers_button.text() == "3 of 3"
    assert not d.weights.isChecked() and d.min_cells.value() == 10
    assert d.correction.currentData() == compare.BH and not d.scale.isVisibleTo(d)
    expected = ex.differential(m.values, m.counts, m.clusters, m.markers, groups,
                               ("Case", "Ctrl"))
    assert len(d.results) == 12 and d.results[0].n == (5, 5)
    assert [r.p for r in d.results[:9]] == pytest.approx([r.p for r in expected[:9]])
    assert d.table.rowCount() == 12 and d.table.columnCount() == len(d.result.columns)
    assert [d.table.horizontalHeaderItem(i).text() for i in range(2)] == ["Cluster", "Marker"]
    # the cluster too rare in Case is a row that says so, and a grey line of the heatmap
    rare = d.results[9]
    assert np.isnan(rare.p) and "fewer than 2 samples of Case" in rare.note
    assert d.result.column("Note")[9] == rare.note
    assert len(d.heatmap.panels) == 1 and np.isnan(d.heatmap.content().cells()[1][3]).all()
    assert len(d.volcano.panels) == 1 and len(d.volcano.content().points()) == 9
    text = d.statement.text()
    assert "moderated t" in text and "9 test(s)" in text and "fewer than 10 cells" in text
    assert "Benjamini-Hochberg" in text and "biexponential" in text
    assert "higher in Case" in d.legend.text() and "FDR < 0.05" in d.legend.text()

    # fewer cells asked for: the rare cluster is tested
    d.min_cells.setValue(2)
    assert np.isfinite(d.results[9].p) and "12 test(s)" in d.statement.text()
    d.min_cells.setValue(10)
    # fewer markers, the day adjusted for, weights, another correction
    d.set_markers(["CD4 PE-A", "TIGIT"])
    assert len(d.results) == 8 and d.markers_button.text() == "2 of 3"
    assert {r.marker for r in d.results} == {"CD4 PE-A", "TIGIT"}
    d.adjust.setCurrentIndex(d.adjust.findData("Day"))
    d.weights.setChecked(True)
    d.correction.setCurrentIndex(d.correction.findData(compare.HOLM))
    again = ex.differential(m.values[:, :, 1:], m.counts, m.clusters, m.markers[1:], groups,
                            ("Case", "Ctrl"), day, True, correction=compare.HOLM)
    assert [r.p_adjusted for r in d.results[:6]] == pytest.approx(
        [r.p_adjusted for r in again[:6]])
    text = d.statement.text()
    assert "Day as a second factor" in text and "weighing" in text and "Holm" in text
    assert "p adjusted (Holm)" in d.result.columns and "adjusted p < 0.05" in d.legend.text()
    assert d.volcano.content().ylabel == "−log10 p (Holm)"

    path = str(tmp_path / "ds.csv")
    d.save_csv(path)
    with open(path, encoding=cluster_stats.CSV_ENCODING) as f:
        lines = f.read().splitlines()
    assert lines[0].startswith("Cluster,Marker,Case (median),Ctrl (median),Difference,p,")
    assert len(lines) == 9
    with open(path, "rb") as f:
        assert f.read(3) == b"\xef\xbb\xbf"
    d.copy()
    assert QtWidgets.QApplication.clipboard().text().startswith("Cluster\tMarker\tCase (median)")
    for view in (d.heatmap, d.volcano):
        figure = FigureDialog(view.content(), d, "expression")
        assert not figure.categories.isVisibleTo(figure)
        assert figure.annotate.isVisibleTo(figure) == (view is d.heatmap)
        figure.close()

    # nothing to compare, or nothing to compare it on, is said, not raised
    d.set_markers([])
    assert d.results == [] and "at least one marker" in d.statement.text()
    assert d.heatmap.content() is None and d.volcano.content() is None
    d.set_markers(m.markers)
    d.first.setCurrentIndex(d.first.findData("Ctrl"))
    assert d.results == [] and "Choose two" in d.statement.text()
    lone = ExpressionDialog(m, ["x"] * 12)
    assert lone.results == [] and lone.table.rowCount() == 0

    # within each day: a panel and a set of rows for each
    w = ExpressionDialog(m, groups, "Group", {}, day, "Day", chosen=["TIGIT"])
    assert w.chosen_markers() == ["TIGIT"]
    assert len(w.results) == 8 and len(w.heatmap.panels) == 2 and len(w.volcano.panels) == 2
    assert w.result.columns[0] == "Day" and "separately within each Day" in w.statement.text()

    # another scale is asked of whoever made the medians
    asked = []

    def reload(scale):
        asked.append(scale)
        return ex.Medians(m.clusters, m.samples, m.markers, m.values / 100, m.counts, scale)

    r = ExpressionDialog(m, groups, reload=reload)
    assert r.scale.isVisibleTo(r) and r.scale.currentData() == "biex"
    before = r.results[1].difference
    r.scale.setCurrentIndex(r.scale.findData("arcsinh"))
    assert asked == ["arcsinh"] and r.results[1].difference == pytest.approx(before / 100)
    assert "arcsinh(x / 150)" in r.statement.text()


def test_from_the_clusters_tab(app, folder):
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
    w.explore.set_markers(["FITC-A", "PE-A"])            # Foxp3 is left to be a state
    w.explore.clusters.setValue(3)
    w.explore.grid.setValue(5)
    w.explore.run(["FlowSOM", "PhenoGraph"])
    settle()
    c, s, ds = w.clusters, w.session, w.session.dataset
    assert c.expression_button.text() == "Differential expression…"
    ids = s.clusters("FlowSOM").ids
    clusters.rename(ds, {ids[0]: "Tregs"})
    clusters.merge(ds, [ids[1], ids[2]])
    c.clusters_changed()
    now = s.clusters("FlowSOM")
    d = c.expression_dialog()
    m = d.medians
    # the clusters as named and merged; every stained marker, the unused one ticked
    assert m.clusters == now.names and "Tregs" in m.clusters and len(m.clusters) == 2
    assert m.markers == ["CD8", "CD4", "Foxp3"] and d.chosen_markers() == ["Foxp3"]
    assert m.samples == ds.samples and m.counts.sum() == len(ds) and m.scale == "biex"
    X = s.features(["FITC-A", "PE-A", "APC-A"])
    mask = (now.labels == now.ids[1]) & (ds.sample_index == 2)
    assert m.counts[1, 2] == mask.sum()
    assert m.values[1, 2] == pytest.approx(np.median(X[mask], axis=0))
    assert d.groups == c._groups() and d.first.currentData() == "Case"
    assert len(d.results) == 2 and d.results[0].cluster == "Tregs"
    assert all(r.n == (2, 2) and np.isfinite(r.p) for r in d.results)
    # on another scale the medians are made again
    d.scale.setCurrentIndex(d.scale.findData("arcsinh"))
    assert d.medians.scale == "arcsinh" and d.medians.values.max() < 20
    assert "arcsinh" in d.statement.text()
    # and of the other clustering, where that is the one shown
    c.set_clustering("PhenoGraph")
    p = c.expression_dialog()
    assert p.medians.clusters == s.clusters("PhenoGraph").names
    assert p.medians.counts.sum() == len(ds)
    assert s.cluster_medians(key="nothing") is None
    w.close()
