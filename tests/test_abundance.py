"""Differential abundance: the model against calculations made another
way, the test against simulations where the truth is known, and the
dialog driven without a screen."""

import os

import numpy as np
import pytest

from milc import abundance, cluster_stats, compare, figures
from milc.volcano import DOWN, REST, UP, VolcanoContent

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def simulate(rng, n=(5, 5), clusters=15, fold=1.0, changed=3):
    """Counts of `clusters` clusters, from 0.1% to 30% of the sample, in
    n[0] samples of group A and n[1] of B. Each cluster's share varies
    between samples as a gamma with a dispersion of its own (0.02 to 0.3),
    the samples' totals vary five-fold, and the counts are multinomial.
    With `fold`, `changed` of the clusters under 10% are that many times
    larger in A; the largest cluster gives up the difference, so that
    every other share is truly the same — it is itself left unjudged.
    Returns the counts, the groups and which clusters truly changed."""
    p = np.geomspace(0.001, 0.3, clusters)
    p /= p.sum()
    phi = np.exp(rng.uniform(np.log(0.02), np.log(0.3), clusters))
    truth = np.zeros(clusters, bool)
    pa = p.copy()
    if fold != 1.0:
        truth[rng.choice(np.flatnonzero(p < 0.1), changed, replace=False)] = True
        pa[truth] *= fold
        pa[-1] -= pa.sum() - 1
    groups = ["A"] * n[0] + ["B"] * n[1]
    counts = np.empty((clusters, sum(n)))
    for i, g in enumerate(groups):
        q = (pa if g == "A" else p) * rng.gamma(1 / phi, phi)
        total = int(np.exp(rng.uniform(np.log(10000), np.log(50000))))
        counts[:, i] = rng.multinomial(total, q / q.sum())
    return counts, groups, truth


def t_test(counts, groups):
    """What the Clusters tab does already: Welch's t on the percentages."""
    pct = 100 * counts / counts.sum(axis=0)
    values = {(k,): {g: [pct[k, i] for i, x in enumerate(groups) if x == g] for g in "AB"}
              for k in range(len(pct))}
    found = compare.compare(values)
    return np.array([c.p for c in found]), np.array([c.p_adjusted for c in found])


def tally(reps, fold, seed, n=(5, 5), clusters=15):
    """Over `reps` simulated data sets, for each method: the share of
    truly unchanged clusters with p < 0.05, the share of truly changed
    ones found at FDR < 0.05, and the false discovery rate that was had:
    the share of a data set's discoveries that are false, averaged."""
    rng = np.random.default_rng(seed)
    names = [str(i) for i in range(clusters)]
    kept = {m: {"p": [], "fdp": [], "hit": 0, "true": 0}
            for m in (abundance.QL, abundance.LRT, "t")}
    for _ in range(reps):
        counts, groups, truth = simulate(rng, n, clusters, fold)
        judged = np.ones(clusters, bool)
        judged[-1] = fold == 1.0
        for m, d in kept.items():
            if m == "t":
                p, q = t_test(counts, groups)
            else:
                found = abundance.differential(counts, names, groups, ("A", "B"), test=m)
                p = np.array([r.p for r in found])
                q = np.array([r.p_adjusted for r in found])
            significant = (q < 0.05) & judged
            d["p"] += list(p[judged & ~truth])
            d["fdp"].append((significant & ~truth).sum() / max(significant.sum(), 1))
            d["hit"] += (significant & truth).sum()
            d["true"] += truth.sum()
    return {m: (float(np.mean(np.array(d["p"]) < 0.05)), d["hit"] / max(d["true"], 1),
                float(np.mean(d["fdp"]))) for m, d in kept.items()}


# --- the model ---------------------------------------------------------------------

def test_log_likelihood_is_the_negative_binomials():
    from scipy import stats as st
    y = np.array([0.0, 1, 7, 150, 30000])
    mu = np.array([0.3, 2.5, 4.0, 180.0, 25000.0])
    for phi in (1e-6, 0.05, 0.4, 3.0):
        expected = st.nbinom.logpmf(y, 1 / phi, 1 / (1 + phi * mu))
        assert np.allclose(abundance.log_likelihood(y, mu, phi), expected, rtol=1e-9, atol=1e-7)
    assert np.allclose(abundance.log_likelihood(y, mu, 0.0), st.poisson.logpmf(y, mu))
    # the deviance is twice the likelihood lost against a perfect fit
    for phi in (0.0, 0.2):
        lost = 2 * (abundance.log_likelihood(y, np.maximum(y, 1e-300), phi)
                    - abundance.log_likelihood(y, mu, phi)).sum()
        assert np.isclose(abundance.deviance(y, mu, phi), lost)


def test_poisson_fit_is_the_ratio_of_pooled_rates():
    rng = np.random.default_rng(1)
    totals = rng.integers(5000, 40000, 8).astype(float)
    groups = ["A"] * 4 + ["B"] * 4
    y = rng.poisson(totals * np.array([[0.01] * 4 + [0.03] * 4, [0.2] * 8]))
    X = abundance.design(groups, ("A", "B"))
    beta, mu = abundance.glm_fit(y, X, np.log(totals), 0.0)
    rate_a, rate_b = y[:, :4].sum(axis=1) / totals[:4].sum(), y[:, 4:].sum(axis=1) / totals[4:].sum()
    assert np.allclose(beta[:, 0], np.log(rate_b)) and np.allclose(beta[:, 1], np.log(rate_a / rate_b))
    assert np.allclose(mu[:, :4], rate_a[:, None] * totals[:4])


def test_negative_binomial_fit_is_the_maximum_of_scipys_likelihood():
    from scipy import optimize, stats as st
    rng = np.random.default_rng(2)
    groups = ["A"] * 5 + ["B"] * 5
    batch = ["x", "y", "z", "x", "y", "z", "x", "y", "z", "x"]
    X = abundance.design(groups, ("A", "B"), batch)
    assert X.shape == (10, 4) and X[:, 1].tolist() == [1] * 5 + [0] * 5
    totals = rng.integers(8000, 40000, 10).astype(float)
    y = rng.negative_binomial(4, 4 / (4 + totals * 0.01)).astype(float)
    phi = 0.3

    def minus(b):
        return -st.nbinom.logpmf(y, 1 / phi, 1 / (1 + phi * np.exp(X @ b + np.log(totals)))).sum()

    best = optimize.minimize(minus, [np.log(0.01), 0, 0, 0], method="BFGS", options={"gtol": 1e-9})
    beta, mu = abundance.glm_fit(y, X, np.log(totals), phi)
    assert np.allclose(beta[0], best.x, atol=2e-5)
    assert np.isclose(-abundance.log_likelihood(y, mu[0], phi).sum(), best.fun)


def test_dispersion_against_a_direct_maximisation():
    from scipy import optimize, stats as st
    rng = np.random.default_rng(3)
    groups = ["A"] * 4 + ["B"] * 4
    X = abundance.design(groups, ("A", "B"))
    totals = rng.integers(8000, 40000, 8).astype(float)
    share = np.array([0.002, 0.01, 0.05, 0.2, 0.004, 0.03])[:, None]
    y = rng.negative_binomial(5, 5 / (5 + totals * share)).astype(float)

    def profile(log_phi):       # Cox-Reid, written out with scipy alone
        phi, total = np.exp(log_phi), 0.0
        for row in y:
            def minus(b):
                mu = np.exp(X @ b + np.log(totals))
                return -st.nbinom.logpmf(row, 1 / phi, 1 / (1 + phi * mu)).sum()
            best = optimize.minimize(minus, [np.log(row.mean() / totals.mean()), 0], method="BFGS")
            mu = np.exp(X @ best.x + np.log(totals))
            info = X.T @ (X * (mu / (1 + phi * mu))[:, None])
            total += -best.fun - 0.5 * np.linalg.slogdet(info)[1]
        return -total

    direct = np.exp(optimize.minimize_scalar(profile, bounds=(np.log(1e-3), np.log(4)),
                                             method="bounded").x)
    common, own = abundance.dispersions(y, X, np.log(totals))
    assert np.isclose(common, direct, rtol=0.03)
    assert np.allclose(abundance.adjusted_profile_likelihood(y, X, np.log(totals), direct).sum(),
                       -profile(np.log(direct)), atol=1e-4)
    # each cluster's own is pulled towards the common one, not onto it
    alone = abundance.dispersions(y, X, np.log(totals), prior_df=0.0)[1]
    assert np.all(np.abs(np.log(own / common)) <= np.abs(np.log(alone / common)) + 1e-9)
    assert np.ptp(own) > 0

    # and with many clusters of one true dispersion, that is what is found
    many = rng.negative_binomial(10, 10 / (10 + totals * np.full((300, 1), 0.02))).astype(float)
    assert np.isclose(abundance.dispersions(many, X, np.log(totals))[0], 0.1, rtol=0.12)


def test_squeeze_finds_the_prior_it_was_simulated_with():
    from scipy.special import polygamma
    assert np.isclose(float(polygamma(1, abundance._trigamma_inverse(0.37))), 0.37)
    rng = np.random.default_rng(4)
    df, prior_df, prior = 8.0, 6.0, 1.7
    true = prior_df * prior / rng.chisquare(prior_df, 4000)
    s2 = true * rng.chisquare(df, 4000) / df
    squeezed, found = abundance.squeeze(s2, np.full(4000, df))
    assert np.isclose(found, prior_df, rtol=0.15)
    centre = np.median((squeezed * (df + found) - df * s2) / found)
    assert np.isclose(centre, prior, rtol=0.1)
    assert np.all(np.abs(np.log(squeezed / prior)) <= np.abs(np.log(s2 / prior)) + 0.1)
    # variances that differ no more than chance allows are one variance
    same, found = abundance.squeeze(np.full(30, 2.0), np.full(30, 50.0))
    assert found == float("inf") and np.ptp(same) == 0 and np.isclose(same[0], 2.0, rtol=0.03)
    # one alone has nobody to learn from
    alone, found = abundance.squeeze(np.array([1.3]), np.array([8.0]))
    assert found == 0.0 and alone[0] == 1.3


# --- calibration and power ---------------------------------------------------------

def test_p_values_under_the_null_and_power_against_the_t_test():
    null = tally(150, 1.0, seed=0)
    # 2,250 p values each: 5% are expected below 0.05, give or take one
    assert 0.03 < null[abundance.QL][0] < 0.07
    assert 0.03 < null[abundance.LRT][0] < 0.08
    assert 0.02 < null["t"][0] < 0.07
    changed = tally(100, 2.0, seed=1)
    ql, lrt, t = changed[abundance.QL], changed[abundance.LRT], changed["t"]
    assert ql[1] > 0.5 and lrt[1] > 0.6
    assert ql[1] > t[1] + 0.15           # the counts find what the percentages miss
    assert ql[2] < 0.15 and lrt[2] < 0.18  # and not by calling everything


# --- what can go wrong ---------------------------------------------------------------

def test_a_cluster_absent_from_a_group_and_other_edges():
    rng = np.random.default_rng(5)
    counts, groups, _ = simulate(rng, clusters=8)
    names = [str(i + 1) for i in range(8)]
    counts[0, :5] = 0                                # none in A
    counts[1, 5:] = 0                                # none in B
    counts[2] = 0                                    # none anywhere
    for test in (abundance.QL, abundance.LRT):
        found = abundance.differential(counts, names, groups, ("A", "B"), test=test)
        assert [r.cluster for r in found] == names and found[0].n == (5, 5)
        assert np.isfinite(found[0].log2_fold) and found[0].log2_fold < -3
        assert np.isfinite(found[1].log2_fold) and found[1].log2_fold > 3
        assert "no cells in A" in found[0].note and "no cells in B" in found[1].note
        assert found[0].p < 0.05 and found[1].p < 0.05
        assert np.isnan(found[2].p) and np.isnan(found[2].log2_fold)
        assert found[2].note == "no cells in either group" and found[2].percent == (0.0, 0.0)
        ok = [r for r in found if r.p == r.p]
        assert len(ok) == 7 and np.allclose([r.p_adjusted for r in ok],
                                            compare.adjust([r.p for r in ok]))
        assert all(r.test == abundance.TESTS[test] and r.dispersion > 0 for r in ok)
    assert np.isclose(found[3].percent[0], (100 * counts[3, :5] / counts[:, :5].sum(axis=0)).mean())

    # the direction is the contrast's
    back = abundance.differential(counts, names, groups, ("B", "A"))
    assert np.allclose([r.log2_fold for r in back[3:]], [-r.log2_fold for r in
                       abundance.differential(counts, names, groups, ("A", "B"))[3:]], atol=1e-6)

    # samples of other groups, and samples with nothing in them, are not there
    wider = np.column_stack([counts, counts[:, :2] * 3, np.zeros(8)])
    more = abundance.differential(wider, names, groups + ["C", "C", "A"], ("A", "B"))
    assert more[4].n == (5, 5)
    assert np.allclose([r.p for r in more[3:]],
                       [r.p for r in abundance.differential(counts, names, groups, ("A", "B"))[3:]])

    # one sample in a group: a fold change, no test, and the reason
    one = abundance.differential(counts[:, 4:], names, groups[4:], ("A", "B"))
    assert one[4].n == (1, 5) and np.isnan(one[4].p) and np.isnan(one[4].p_adjusted)
    assert one[4].note == "a group has fewer than 2 samples" and np.isfinite(one[4].log2_fold)
    none = abundance.differential(counts, names, groups, ("A", "Z"))
    assert none[4].note == "one group: nothing to compare" and np.isnan(none[4].log2_fold)

    empty = abundance.differential(np.zeros((2, 10)) + [[0.0], [5.0]], ["1", "2"], groups,
                                   ("A", "B"))
    assert empty[0].note == "no cells in either group" and np.isfinite(empty[1].p)
    assert abundance.differential(np.zeros((2, 10)), ["1", "2"], groups, ("A", "B"))[0].n == (0, 0)

    # every sample the same: nothing to find, and no error on the way
    flat = abundance.differential(np.tile([[10.0], [200.0], [3000.0]], (1, 6)), ["1", "2", "3"],
                                  ["A"] * 3 + ["B"] * 3, ("A", "B"))
    assert [r.p for r in flat] == [1.0, 1.0, 1.0]
    assert np.allclose([r.log2_fold for r in flat], 0.0, atol=1e-9)


def test_adjusting_for_a_second_factor():
    rng = np.random.default_rng(6)
    groups = ["A", "A", "A", "B", "B", "B"] * 2
    day = ["1"] * 6 + ["2"] * 6
    share = np.array([0.01, 0.05, 0.1, 0.3, 0.02, 0.2])[:, None] * np.ones(12)
    share[1, 6:] *= 4                                # day 2 has four times cluster 2 ...
    share[1, [0, 1, 2, 6, 7, 8]] *= 1.6              # ... and A 1.6 times, on either day
    totals = rng.integers(15000, 40000, 12)
    counts = rng.poisson(totals * share * rng.gamma(40, 1 / 40, share.shape)).astype(float)
    names = [str(i + 1) for i in range(6)]
    plain = abundance.differential(counts, names, groups, ("A", "B"), totals=totals)
    adjusted = abundance.differential(counts, names, groups, ("A", "B"), day, totals=totals)
    assert adjusted[1].p < 0.01 < plain[1].p
    assert np.isclose(adjusted[1].log2_fold, np.log2(1.6), atol=0.35)

    # a factor that is the groups under another name cannot be adjusted for
    same = abundance.differential(counts, names, groups, ("A", "B"), ["x" if g == "A" else "y"
                                                                     for g in groups])
    assert np.isnan(same[1].p) and "cannot be told apart" in same[1].note
    few = abundance.differential(counts[:, :4], names, ["A", "A", "B", "B"], ("A", "B"),
                                 ["1", "2", "3", "4"])
    assert np.isnan(few[1].p) and few[1].note

    # within each day: each day's own model, one adjustment over both
    both = abundance.within(counts, names, groups, ("A", "B"), day)
    first = abundance.differential(counts[:, :6], names, groups[:6], ("A", "B"))
    assert [r.panel for r in both] == ["1"] * 6 + ["2"] * 6
    assert np.allclose([r.p for r in both[:6]], [r.p for r in first])
    assert np.allclose([r.p_adjusted for r in both], compare.adjust([r.p for r in both]))

    table = abundance.to_table(both, panel="Day")
    assert table.columns[:3] == ["Day", "Cluster", "A (mean %)"] and len(table.rows) == 12
    assert all(len(row) == len(table.columns) for row in table.rows)
    assert abundance.to_table(plain).columns[0] == "Cluster"
    assert table.column("p adjusted (Benjamini-Hochberg)")[1] == both[1].p_adjusted


# --- the volcano ---------------------------------------------------------------------

def results():
    rng = np.random.default_rng(7)
    counts, groups, _ = simulate(rng, fold=3.0, changed=2)
    counts[0, 5:] = 0
    return abundance.differential(counts, [str(i + 1) for i in range(15)], groups, ("A", "B"))


def test_volcano_content_and_figure(tmp_path):
    found = results()
    content = VolcanoContent(found)
    pts = content.points()
    assert len(pts) == 15 and content.panels() == [None]
    kinds = {cluster: kind for _, _, cluster, kind in pts}
    for r in found:
        expected = REST if r.p_adjusted >= 0.05 else UP if r.log2_fold > 0 else DOWN
        assert kinds[r.cluster] == expected
    assert UP in kinds.values() and REST in kinds.values()
    x, y, cluster, _ = pts[3]
    assert x == found[3].log2_fold and np.isclose(y, -np.log10(found[3].p_adjusted))
    half, top = content.limits()
    assert half > max(abs(r.log2_fold) for r in found) and top > -np.log10(0.05)
    assert content.xlabel == "log2 fold change (A / B)" and content.ylabel == "−log10 FDR"
    assert content.names()[UP] == "More in A" and content.names()[DOWN] == "More in B"

    style = figures.FigureStyle(font="DejaVu Sans")
    assert figures.preview_png(content, style, dpi=40)[:4] == b"\x89PNG"
    path = str(tmp_path / "volcano.svg")
    figures.save(content, style, path)
    with open(path, encoding="utf-8") as f:
        svg = f.read()
    named = [c for c, k in kinds.items() if k != REST]
    assert "More in A" in svg and "log2 fold change (A / B)" in svg
    assert all(f">{c}</text>" in svg or f">{c}<" in svg for c in named)

    # one panel for each value of "within each"; untested clusters are not points
    for r in found[:7]:
        r.panel = "3"
    for r in found[7:]:
        r.panel = "7"
    found[1].p_adjusted = float("nan")
    panelled = VolcanoContent(found)
    assert panelled.panels() == ["3", "7"]
    assert len(panelled.points("3")) == 6 and len(panelled.points("7")) == 8
    assert figures.preview_png(panelled, style, dpi=40)[:4] == b"\x89PNG"
    assert VolcanoContent([]).panels() == [None] and VolcanoContent([]).limits()[0] > 1


# --- the window ----------------------------------------------------------------------

@pytest.fixture(scope="module")
def app(tmp_path_factory):
    QtWidgets = pytest.importorskip("PyQt6.QtWidgets")
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


def test_dialog(app, tmp_path):
    from PyQt6 import QtWidgets
    from milc.ui.abundance_dialog import AbundanceDialog
    from milc.ui.figure_dialog import FigureDialog
    rng = np.random.default_rng(8)
    counts, groups, _ = simulate(rng, n=(6, 6), clusters=8, fold=3.0, changed=2)
    table = cluster_stats.Table([str(i + 1) for i in range(8)], [f"s{i}" for i in range(12)], counts)
    groups = ["Case" if g == "A" else "Ctrl" for g in groups]
    groups[5], groups[11] = "UV", "UV"
    day = ["3", "7"] * 6
    d = AbundanceDialog(table, groups, "Group", {"Day": day})
    assert [d.first.itemData(i) for i in range(d.first.count())] == ["Case", "Ctrl", "UV"]
    assert (d.first.currentData(), d.second.currentData()) == ("Case", "Ctrl")
    assert len(d.results) == 8 and d.results[0].groups == ("Case", "Ctrl") and d.results[0].n == (5, 5)
    expected = abundance.differential(counts, table.rows, groups, ("Case", "Ctrl"))
    assert np.allclose([r.p for r in d.results], [r.p for r in expected])
    assert d.table.rowCount() == 8 and d.table.columnCount() == len(d.result.columns)
    assert d.table.horizontalHeaderItem(0).text() == "Cluster"
    assert len(d.volcano.panels) == 1 and d.volcano.content().results is d.results
    text = d.statement.text()
    assert "negative-binomial" in text and "quasi-likelihood F test" in text
    assert "not checked against edgeR" in text and "8 test(s)" in text

    # the other two, the other test, adjusted for the day
    d.first.setCurrentIndex(d.first.findData("Ctrl"))
    assert d.results == [] and "Choose two" in d.statement.text() and d.volcano.content() is None
    d.second.setCurrentIndex(d.second.findData("Case"))
    assert np.allclose([r.log2_fold for r in d.results], [-r.log2_fold for r in expected], atol=1e-6)
    d.test.setCurrentIndex(d.test.findData(abundance.LRT))
    d.adjust.setCurrentIndex(d.adjust.findData("Day"))
    with_day = abundance.differential(counts, table.rows, groups, ("Ctrl", "Case"), day,
                                      abundance.LRT)
    assert np.allclose([r.p for r in d.results], [r.p for r in with_day])
    assert "adjusted for Day" in d.statement.text() and "likelihood-ratio" in d.statement.text()

    path = str(tmp_path / "da.csv")
    d.save_csv(path)
    with open(path, encoding=cluster_stats.CSV_ENCODING) as f:
        lines = f.read().splitlines()
    assert lines[0].startswith("Cluster,Ctrl (mean %),Case (mean %),log2 fold change,p,")
    assert len(lines) == 9
    d.copy()
    assert QtWidgets.QApplication.clipboard().text().startswith("Cluster\tCtrl (mean %)")
    figure = FigureDialog(d.volcano.content(), d, "volcano")
    assert not figure.categories.isVisibleTo(figure)
    figure.close()

    # within each day: a panel and a set of rows for each
    w = AbundanceDialog(table, groups, "Group", {}, day, "Day")
    assert len(w.results) == 16 and len(w.volcano.panels) == 2
    assert w.result.columns[0] == "Day" and w.table.rowCount() == 16
    assert "separately within each Day" in w.statement.text()
    assert w.adjust.count() == 1

    # nothing to compare is said, not raised
    lone = AbundanceDialog(table, ["x"] * 12)
    assert lone.results == [] and lone.table.rowCount() == 0


def test_from_the_clusters_tab(app, tmp_path):
    from milc.fcs import write_fcs
    from milc.ui import worker
    from milc.ui.shell import MainShell
    rng = np.random.default_rng(0)
    centres = np.array([[300, 3000, 100], [3000, 200, 100], [200, 300, 4000.0]])
    for name, n in (("3Ctrl1", 300), ("3Ctrl2", 400), ("3Case1", 500), ("3Case2", 600)):
        pop = rng.choice(3, n, p=[0.6, 0.3, 0.1] if "Case" in name else [0.2, 0.3, 0.5])
        fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (n, 3)))
        ev = np.column_stack([rng.normal(90000, 9000, n), fluo]).astype(np.float32)
        write_fcs(str(tmp_path / f"{name}.fcs"), ev, ["FSC-A", "FITC-A", "PE-A", "APC-A"],
                  markers=["", "CD8", "CD4", "Foxp3"],
                  keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(str(tmp_path))
    w.samples.concatenate()
    settle()
    w.explore.clusters.setValue(3)
    w.explore.grid.setValue(5)
    w.explore.run(["FlowSOM"])
    settle()
    c = w.clusters
    assert c.abundance_button.text() == "Differential abundance…"
    d = c.abundance_dialog()
    ds = w.session.dataset
    assert d.counts.values.sum() == len(ds) and d.counts.columns == ds.samples
    assert len(d.results) == 3 and d.results[0].groups == ("Case", "Ctrl")
    assert d.results[0].n == (2, 2) and all(np.isfinite(r.p) for r in d.results)
    assert min(r.p_adjusted for r in d.results) < 0.05     # 60% against 20% of 300 cells
    assert d.facets is None and d.adjust.count() == 1
    d.close()
    w.close()
