"""The designs beyond "different samples, equal variances": the same
subjects in every group (paired t, Wilcoxon, repeated-measures ANOVA,
Friedman) and Welch's ANOVA. Each is set beside scipy where scipy has the
test, and beside a published example and a second calculation by another
route where it does not."""

import inspect
import os

import numpy as np
import pytest
from scipy import stats as st

from milc import compare as c
from milc.stats import StatsTable

# Field, "Discovering Statistics Using SPSS", repeated-measures chapter:
# eight celebrities, the time to retch after each of four foods. The text
# reports F(3, 21) = 3.79, p = .026; Greenhouse-Geisser ε = .533, p = .063.
# The figures were typed from memory of the text, not from the page: that
# they give the published F, ε and both p values is what vouches for them.
BUSHTUCKER = np.array([[8, 7, 1, 6], [9, 5, 2, 5], [6, 2, 3, 8], [5, 3, 1, 9],
                       [8, 4, 5, 8], [7, 5, 6, 7], [10, 2, 7, 2], [12, 6, 8, 1.0]])

# R's PlantGrowth; `oneway.test(weight ~ group, PlantGrowth)` prints
# F = 5.181, num df = 2.000, denom df = 17.128, p-value = 0.01739
PLANTS = [[4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14],
          [4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69],
          [6.31, 5.12, 5.54, 5.50, 5.37, 5.29, 4.92, 6.15, 5.80, 5.26]]


def test_two_conditions_on_the_same_subjects():
    before = [12.0, 15, 9, 20, 14, 11, 17]
    after = [14.0, 19, 9.5, 26, 13, 15, 22]
    name, t, p, note = c.test_paired([before, after])
    assert name == "Paired t" and note == ""
    r = st.ttest_rel(before, after)
    assert t == pytest.approx(r.statistic) and p == pytest.approx(r.pvalue)
    # and by hand: the mean change over its standard error, on n - 1
    d = np.subtract(before, after)
    by_hand = d.mean() / (d.std(ddof=1) / np.sqrt(d.size))
    assert t == pytest.approx(by_hand)
    assert p == pytest.approx(2 * st.t.sf(abs(by_hand), d.size - 1))
    # pairing is the point: the same numbers as two separate groups say less
    assert p < c.test_groups([before, after])[2]
    assert c.test_paired([before, after], c.WELCH)[:3] == (name, t, p)

    name, w, p, note = c.test_paired([before, after], c.RANK)
    r = st.wilcoxon(before, after)
    assert name == "Wilcoxon signed-rank" and note == ""
    assert w == pytest.approx(r.statistic) and p == pytest.approx(r.pvalue)


def test_ties_and_subjects_that_do_not_change():
    a = [5.0, 7, 7, 9, 4, 6, 8, 3]
    b = [6.0, 8, 7, 11, 5, 5, 10, 3]       # two do not change, four change by 1
    name, w, p, note = c.test_paired([a, b], c.RANK)
    r = st.wilcoxon(a, b)
    assert w == pytest.approx(r.statistic) and p == pytest.approx(r.pvalue) and note == ""
    # those that do not change are out of the ranks: five that do cannot reach 0.05
    _, _, p, note = c.test_paired([a[:6], b[:6]], c.RANK)
    assert p == pytest.approx(st.wilcoxon(a[:6], b[:6]).pvalue)
    assert "smallest p possible with these n: 0.062" in note


def test_too_few_to_say():
    # five pairs, all in one direction: the smallest p a signed-rank test has
    _, _, p, note = c.test_paired([[1.0, 2, 3, 4, 5], [2.0, 4, 6, 8, 10]], c.RANK)
    assert p == pytest.approx(2 / 2 ** 5) and "0.062" in note
    assert c.test_paired([[1.0, 2, 3, 4, 5, 6], [2.0, 4, 6, 8, 10, 12]], c.RANK)[3] == ""
    # Friedman, two subjects in three groups, perfectly ordered: p is still 0.135
    _, chi, p, note = c.test_paired([[1.0, 2], [3.0, 4], [5.0, 6]], c.RANK)
    assert chi == pytest.approx(4) and p == pytest.approx(np.exp(-2)) and "0.14" in note
    # what cannot be tested is no p value, with the reason
    for columns, kind, why in (
            ([[1.0, 2, 3]], c.PARAMETRIC, "one group"),
            ([[1.0], [2.0]], c.PARAMETRIC, "fewer than 2 subjects"),
            ([[2.0, 2], [2.0, 2]], c.PARAMETRIC, "same value"),
            ([[1.0, 2, 3], [2.0, 3, 4]], c.PARAMETRIC, "same amount"),
            ([[1.0, 2, 3], [1.0, 2, 3]], c.RANK, "no subject differs"),
            ([[1.0, 2, 3], [2.0, 3, 4], [4.0, 5, 6]], c.PARAMETRIC, "no variation left")):
        _, statistic, p, note = c.test_paired(columns, kind)
        assert np.isnan(p) and np.isnan(statistic) and why in note


def _rm_by_another_route(data):
    """F as the comparison of two least-squares fits — subjects and
    groups against subjects alone — and ε from the eigenvalues of the
    covariance of orthonormal contrasts (Box's form), not from traces of
    the centred covariance as the code has it."""
    from scipy.linalg import null_space
    n, k = data.shape
    y = data.ravel()
    subject = np.kron(np.eye(n), np.ones((k, 1)))
    group = np.kron(np.ones((n, 1)), np.eye(k))

    def rss(x):
        return float(((y - x @ np.linalg.lstsq(x, y, rcond=None)[0]) ** 2).sum())

    full, reduced = rss(np.hstack([subject, group])), rss(subject)
    f = ((reduced - full) / (k - 1)) / (full / ((n - 1) * (k - 1)))
    contrasts = null_space(np.ones((1, k))).T
    roots = np.linalg.eigvalsh(contrasts @ np.cov(data, rowvar=False) @ contrasts.T)
    return f, roots.sum() ** 2 / ((k - 1) * (roots ** 2).sum())


def test_repeated_measures_anova():
    f, df1, df2, epsilon, p, corrected = c.rm_anova(BUSHTUCKER)
    # the published figures, to the digits they are printed with
    assert (df1, df2) == (3, 21)
    assert round(f, 2) == 3.79 and round(p, 3) == 0.026
    assert round(epsilon, 3) == 0.533 and round(corrected, 3) == 0.063
    # the same by another route, to full precision
    f2, epsilon2 = _rm_by_another_route(BUSHTUCKER)
    assert f == pytest.approx(f2) and epsilon == pytest.approx(epsilon2)
    assert corrected == pytest.approx(st.f.sf(f2, epsilon2 * 3, epsilon2 * 21))
    rng = np.random.default_rng(3)
    for n, k in ((5, 3), (4, 6), (12, 4)):   # fewer subjects than groups too
        data = rng.normal(0, 1, (n, k)) * rng.uniform(0.5, 3, k) + rng.normal(0, 2, (n, 1))
        f, _, _, epsilon, _, corrected = c.rm_anova(data)
        f2, epsilon2 = _rm_by_another_route(data)
        assert f == pytest.approx(f2) and epsilon == pytest.approx(epsilon2)
        assert 1 / (k - 1) - 1e-12 <= epsilon <= 1 + 1e-12
    # two groups: it is the paired t test, and there is nothing to correct
    f, _, _, epsilon, p, corrected = c.rm_anova(BUSHTUCKER[:, :2])
    r = st.ttest_rel(BUSHTUCKER[:, 0], BUSHTUCKER[:, 1])
    assert f == pytest.approx(r.statistic ** 2) and p == pytest.approx(r.pvalue)
    assert epsilon == pytest.approx(1) and corrected == pytest.approx(p)

    # what a comparison reports is the corrected p, and it says so
    name, f, p, note = c.test_paired(list(BUSHTUCKER.T))
    assert name == "RM ANOVA (Greenhouse-Geisser)" and round(p, 3) == 0.063
    assert note == "ε = 0.533; p without the correction: 0.026"
    name, chi, p, note = c.test_paired(list(BUSHTUCKER.T), c.RANK)
    r = st.friedmanchisquare(*BUSHTUCKER.T)
    assert name == "Friedman" and note == ""
    assert chi == pytest.approx(r.statistic) and p == pytest.approx(r.pvalue)


def test_welch_anova():
    f, df1, df2, p = c.welch_anova([np.array(g) for g in PLANTS])
    assert df1 == 2 and round(f, 3) == 5.181 and round(df2, 3) == 17.128
    assert float(f"{p:.4g}") == 0.01739
    if "equal_var" in inspect.signature(st.f_oneway).parameters:   # scipy 1.16 on
        r = st.f_oneway(*PLANTS, equal_var=False)
        assert f == pytest.approx(r.statistic) and p == pytest.approx(r.pvalue)
        rng = np.random.default_rng(1)
        groups = [rng.normal(m, s, n) for m, s, n in ((0, 1, 4), (1, 5, 9), (0.5, 0.2, 6), (2, 2, 5))]
        r = st.f_oneway(*groups, equal_var=False)
        assert c.welch_anova(groups)[0] == pytest.approx(r.statistic)
        assert c.welch_anova(groups)[3] == pytest.approx(r.pvalue)
    # two groups: Welch's t test, F = t² on Satterthwaite's degrees of freedom
    a, b = np.array(PLANTS[0]), np.array(PLANTS[1][:6]) * 3
    f, _, df2, p = c.welch_anova([a, b])
    r = st.ttest_ind(a, b, equal_var=False)
    assert f == pytest.approx(r.statistic ** 2) and p == pytest.approx(r.pvalue)
    assert df2 == pytest.approx(r.df)

    # as a kind of test: its own name, Welch's t for two and for the pairs
    assert c.KINDS[c.WELCH] == "Welch t test / Welch ANOVA"
    name, f, p, note = c.test_groups(PLANTS, c.WELCH)
    assert name == "Welch ANOVA" and round(f, 3) == 5.181 and note == ""
    assert c.test_groups(PLANTS[:2], c.WELCH) == c.test_groups(PLANTS[:2], c.PARAMETRIC)
    found = c.compare({("w",): dict(zip("abc", PLANTS))}, c.WELCH, pairwise=True)
    assert [x.test for x in found] == ["Welch ANOVA", "Welch t", "Welch t", "Welch t"]
    # a group that does not vary would weigh infinitely
    _, _, p, note = c.test_groups([[1.0, 1, 1], [2.0, 3, 4], [5.0, 7, 6]], c.WELCH)
    assert np.isnan(p) and "no variation inside a group" in note
    assert "fewer than 2" in c.test_groups([[1.0], [2.0, 3], [5.0, 7]], c.WELCH)[3]


def test_only_whole_subjects_are_compared():
    values = {("x",): {"a": [1.0, 2, 3, 4, 9], "b": [7.0, 2.5, 4.5, 5, 2, np.nan],
                       "c": [3.0, 3, 9, 8, 1]}}
    # b is in another order, lacks a value for E and has an F nobody else has;
    # the fifth sample of c was never given its animal
    subjects = {("x",): {"a": list("ABCDE"), "b": list("FABCDE"), "c": ["A", "B", "C", "D", " "]}}
    found = c.compare(values, pairwise=True, subjects=subjects)
    assert [f.test for f in found] == ["RM ANOVA (Greenhouse-Geisser)"] + ["Paired t"] * 3
    assert [f.groups for f in found] == [["a", "b", "c"], ["a", "b"], ["a", "c"], ["b", "c"]]
    a, b, d = [1.0, 2, 3, 4], [2.5, 4.5, 5, 2], [3.0, 3, 9, 8]
    whole = c.rm_anova(np.column_stack([a, b, d]))
    assert found[0].n == [4, 4, 4] and found[0].mean == [2.5, 3.5, 5.75]
    assert found[0].statistic == pytest.approx(whole[0]) and found[0].p == pytest.approx(whole[5])
    # 5 + 5 + 4 samples had a value and a subject; 3 x 4 were used
    assert found[0].dropped == 2 and "2 sample(s) left out" in found[0].note
    assert found[0].note.startswith("ε = ")
    for f, (x, y), dropped in zip(found[1:], ((a, b), (a, d), (b, d)), (2, 1, 1)):
        assert f.p == pytest.approx(st.ttest_rel(x, y).pvalue) and f.dropped == dropped
    by_rank = c.compare(values, c.RANK, pairwise=True, subjects=subjects)
    assert [f.test for f in by_rank] == ["Friedman"] + ["Wilcoxon signed-rank"] * 3
    assert by_rank[0].p == pytest.approx(st.friedmanchisquare(a, b, d).pvalue)
    assert by_rank[2].p == pytest.approx(st.wilcoxon(a, d).pvalue)
    assert "smallest p possible" in by_rank[2].note and "1 sample(s) left out" in by_rank[2].note
    # nobody left out: nothing said
    whole = c.compare({("x",): {"a": a, "b": b}}, subjects={("x",): {"a": "ABCD", "b": "ABCD"}})
    assert whole[0].dropped == 0 and whole[0].note == ""


def test_adjustment_spans_every_test():
    rng = np.random.default_rng(5)
    values, subjects = {}, {}
    for cluster in "123":
        values[(cluster,)] = {g: list(rng.normal(i * int(cluster), 1, 6))
                              for i, g in enumerate(("d0", "d3", "d7"))}
        subjects[(cluster,)] = {g: list("ABCDEF") for g in ("d0", "d3", "d7")}
    values[("3",)]["d7"][0] = float("nan")
    for correction in c.CORRECTIONS:
        found = c.compare(values, c.PARAMETRIC, correction, pairwise=True, subjects=subjects)
        assert len(found) == 12 and all(np.isfinite(f.p) for f in found)
        assert np.allclose([f.p_adjusted for f in found], c.adjust([f.p for f in found], correction))
    assert c.adjust([f.p for f in found], c.BONFERRONI)[0] == pytest.approx(min(found[0].p * 12, 1))
    # a test that could not be done is not one of the tests adjusted for
    subjects[("2",)]["d3"] = list("AACDEF")
    found = c.compare(values, c.RANK, c.BONFERRONI, pairwise=True, subjects=subjects)
    undone = [f for f in found if np.isnan(f.p)]
    assert [f.groups for f in undone] == [["d0", "d3", "d7"], ["d0", "d3"], ["d3", "d7"]]
    assert all(f.note == "subject 'A' is more than once in group 'd3'" and f.n[0] == 6
               and np.isnan(f.statistic) and np.isnan(f.p_adjusted) for f in undone)
    done = [f for f in found if np.isfinite(f.p)]
    assert len(done) == 9 and all(f.p_adjusted == pytest.approx(min(f.p * 9, 1)) for f in done)


def test_a_pairing_that_is_not_one():
    values = {("x",): {"a": [1.0, 2, 3], "b": [2.0, 3, 5]}, ("y",): {"a": [1.0, 2, 3], "b": [4.0, 4, 6]}}
    ids = {"a": ["A", "B", "C"], "b": ["A", "B", "C"]}
    assert c.pairing_problem(values, {("x",): ids, ("y",): ids}) == ""
    twice = {"a": ["A", "B", "C"], "b": ["A", "B", "B"]}
    assert c.pairing_problem(values, {("x",): ids, ("y",): twice}) == \
        "subject 'B' is more than once in group 'b'"
    apart = {"a": ["A", "B", "C"], "b": ["D", "E", "F"]}
    assert c.pairing_problem(values, {("x",): apart, ("y",): apart}) == "no subject is in every group"
    # shared in one of the things compared is enough to go on: the other says why not
    assert c.pairing_problem(values, {("x",): ids, ("y",): apart}) == ""
    found = c.compare(values, subjects={("x",): ids, ("y",): apart})
    assert np.isfinite(found[0].p) and np.isnan(found[1].p) and found[1].n == [0, 0]
    assert found[1].note.startswith("fewer than 2 subjects in every group") and found[1].dropped == 6
    blank = {"a": ["", "", ""], "b": [None, " ", ""]}
    assert c.pairing_problem(values, {("x",): blank, ("y",): blank}) == "no subject is in every group"


def test_without_pairing_nothing_changed():
    """The figures are what `compare` gave before it knew of pairing."""
    a, b, e = [1.0, 2, 3, 4, 5], [4.0, 5, 6, 7, 9], [9.0, 10, 11]
    values = {("1",): {"7": a, "14": b, "3": e},
              ("2",): {"7": a, "14": [2.0, 2, 3, 9, 9], "3": [1.0, 5, np.nan, 2]}}
    was = {c.PARAMETRIC: [
        ("ANOVA", 17.33065442020668, 0.0005627853302036703, 0.003939497311425692),
        ("Welch t", 7.668115805072326, 0.000282916152097784, 0.002263329216782272),
        ("Welch t", 3.6678898769003854, 0.010523855229344707, 0.06314313137606824),
        ("Welch t", -2.8736848324283986, 0.021524004181540508, 0.10762002090770253),
        ("ANOVA", 0.9703599153140438, 0.411952798733171, 1.0),
        ("Welch t", -0.2390457218668788, 0.8246966511032668, 1.0),
        ("Welch t", -1.1461557146618477, 0.29541964369683094, 1.0),
        ("Welch t", -1.118033988749895, 0.3105420432944341, 1.0)], c.RANK: [
        ("Kruskal-Wallis", 9.13019390581717, 0.010408869987616112, 0.08327095990092889),
        ("Mann-Whitney", 15.0, 0.03571428571428571, 0.2490518326771587),
        ("Mann-Whitney", 14.5, 0.051228003006160656, 0.2561400150308033),
        ("Mann-Whitney", 2.0, 0.0355788332395941, 0.2490518326771587),
        ("Kruskal-Wallis", 1.0971428571428605, 0.5777746134000564, 1.0),
        ("Mann-Whitney", 6.5, 0.879339414891054, 1.0),
        ("Mann-Whitney", 4.0, 0.3563733518637908, 1.0),
        ("Mann-Whitney", 9.5, 0.5947293361935915, 1.0)]}
    for kind, expected in was.items():
        found = c.compare(values, kind, c.HOLM, pairwise=True)
        assert [f.groups for f in found] == [["3", "7", "14"], ["3", "7"], ["3", "14"], ["7", "14"]] * 2
        assert [f.n for f in found] == [[3, 5, 5], [3, 5], [3, 5], [5, 5]] * 2
        assert all(f.note == "" and f.dropped == 0 for f in found)
        for f, (test, statistic, p, adjusted) in zip(found, expected):
            # to the last digits: scipy's own arithmetic may differ by a release
            assert f.test == test and f.statistic == pytest.approx(statistic, rel=1e-12)
            assert f.p == pytest.approx(p, rel=1e-12) and f.p_adjusted == pytest.approx(adjusted, rel=1e-12)
    assert c.compare(values, subjects=None)[0].p == pytest.approx(0.0005627853302036703, rel=1e-12)
    table = c.to_table(c.compare(values), ["Cluster"])
    assert table.columns == ["Cluster", "Groups", "Mean ± SD (n)", "Test", "Statistic", "p",
                             "p adjusted (Benjamini-Hochberg)", "", "Note"]
    assert table.rows[0][:4] == ["1", "3 vs 7 vs 14",
                                 "3: 10 ± 1 (3); 7: 3 ± 1.58 (5); 14: 6.2 ± 1.92 (5)", "ANOVA"]
    assert table.rows[0][-2:] == ["**", ""]
    assert list(c.KINDS)[:2] == [c.PARAMETRIC, c.RANK]


def _tissues() -> StatsTable:
    """Five animals, three tissues each; animal 5 has no lung."""
    t = StatsTable(["Sample", "Group", "Animal", "Day", "Population", "Path", "Count", "% of parent"])
    rng = np.random.default_rng(2)
    for animal in "12345":
        own = rng.normal(0, 6)
        for tissue, centre in (("Blood", 20), ("Spleen", 24), ("Lung", 30)):
            if (animal, tissue) != ("5", "Lung"):
                t.rows.append([f"{tissue}{animal}", tissue, animal, "3", "T", "Cells/T", 100,
                               float(centre + own + rng.normal(0, 1))])
    return t


def test_a_statistics_table_paired():
    t = _tissues()
    assert c.from_table(t, "% of parent", ["Group"]) == c.from_table(t, "% of parent", ["Group"])
    assert len(c.from_table(t, "% of parent", ["Group"])) == 2          # as it always was
    values, names, subjects = c.from_table(t, "% of parent", ["Group"], paired_by="Animal")
    assert names == ["Population"] and subjects[("T",)]["Lung"] == list("1234")
    assert subjects[("T",)]["Blood"] == list("12345") and len(values[("T",)]["Blood"]) == 5
    values, names, subjects = c.from_table(t, "% of parent", ["Group"], "Day", "Animal")
    assert names == ["Day", "Population"] and set(subjects) == set(values) == {("3", "T")}
    found = c.compare(values, pairwise=True, subjects=subjects)
    assert found[0].n == [4, 4, 4] and found[0].dropped == 2
    blood, spleen = values[("3", "T")]["Blood"], values[("3", "T")]["Spleen"]
    pair = next(f for f in found if f.groups == ["Blood", "Spleen"])
    assert pair.n == [5, 5] and pair.dropped == 0
    assert pair.p == pytest.approx(st.ttest_rel(blood, spleen).pvalue)
    # the animals differ more than the tissues: unpaired, the same data say less
    assert pair.p < c.compare(values, pairwise=True)[1].p
    table = c.to_table(found, names, paired_by="Animal")
    assert table.columns[:6] == ["Day", "Population", "Groups", "Mean ± SD (n)", "Test", "Paired by"]
    assert all(len(row) == len(table.columns) and row[5] == "Animal" for row in table.rows)
    assert table.rows[0][4] == "RM ANOVA (Greenhouse-Geisser)" and table.rows[0][6] == found[0].statistic


# --- in the window -----------------------------------------------------------------

@pytest.fixture(scope="module")
def app(tmp_path_factory):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6.QtWidgets")
    from PyQt6 import QtWidgets
    os.environ["MILC_SETTINGS"] = str(tmp_path_factory.mktemp("settings") / "s.ini")
    from milc.ui import prepare
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    prepare(app)
    return app


@pytest.fixture
def folder(tmp_path):
    """Four animals, blood, spleen and lung of each."""
    from milc.fcs import write_fcs
    rng = np.random.default_rng(0)
    centres = np.array([[300, 3000, 100], [3000, 200, 100], [200, 300, 4000.0]])
    for animal in range(1, 5):
        for tissue, share in (("Blood", 0.2), ("Spleen", 0.35), ("Lung", 0.5)):
            p = share + 0.03 * animal
            pop = rng.choice(3, 300, p=[p, (1 - p) / 2, (1 - p) / 2])
            fluo = np.abs(centres[pop] * rng.lognormal(0, 0.25, (300, 3)))
            ev = np.column_stack([rng.normal(90000, 9000, 300), fluo]).astype(np.float32)
            write_fcs(str(tmp_path / f"{animal}{tissue}.fcs"), ev,
                      ["FSC-A", "FITC-A", "PE-A", "APC-A"], markers=["", "CD8", "CD4", "Foxp3"],
                      keywords={"SPILL": "3,FITC-A,PE-A,APC-A,1,0.05,0,0.02,1,0,0,0,1"})
    return str(tmp_path)


def test_paired_in_the_window(app, folder):
    from milc.ui import worker
    from milc.ui.shell import MainShell
    from milc.ui.statistics_tab import CompareDialog

    def settle():
        for _ in range(3):
            worker.wait_all()
            app.processEvents()

    w = MainShell()
    w.open_folder(folder)
    assert w.samples.fill_from_names(r"(?P<Animal>\d+)(?P<Group>[A-Za-z]+)") == ["Animal", "Group"]
    w.samples.concatenate()
    settle()
    e = w.explore
    e.clusters.setValue(3)
    e.grid.setValue(5)
    e.run(["FlowSOM"])
    settle()

    # the clusters: not paired until asked, and then the tests of the same subjects
    k = w.clusters
    assert [k.paired_by.itemData(i) for i in range(k.paired_by.count())] == ["", "Group", "Animal"]
    assert k.paired_by.currentData() == "" and k.paired_by.itemText(2) == "Animal"
    assert k.test.findData(c.WELCH) > 0
    k.test.setCurrentIndex(k.test.findData(c.PARAMETRIC))
    assert k.paired_by.isEnabled() and {f.test for f in k.comparisons()} == {"ANOVA"}
    unpaired = k.comparisons(pairwise=True)
    k.paired_by.setCurrentIndex(k.paired_by.findData("Animal"))
    assert k.pairing_problem() == "" and k.pairing.text() == ""
    found = k.comparisons(pairwise=True)
    assert len(found) == len(unpaired) == 12
    assert {f.test for f in found} == {"RM ANOVA (Greenhouse-Geisser)", "Paired t"}
    assert all(f.n == [4] * len(f.groups) and f.dropped == 0 for f in found)
    # the numbers are those of the frequency table, animal against the same animal
    freq, ds = k._frequency, w.session.dataset
    share = {(m["Animal"], m["Group"]): freq.values[0, j] for j, m in enumerate(ds.meta)}
    blood, lung = ([share[a, g] for a in "1234"] for g in ("Blood", "Lung"))
    pair = next(f for f in found if f.key == (None, freq.rows[0]) and f.groups == ["Blood", "Lung"])
    assert pair.p == pytest.approx(st.ttest_rel(blood, lung).pvalue)
    assert np.allclose([f.p_adjusted for f in found],
                       c.adjust([f.p for f in found], k.correction.currentData()), equal_nan=True)
    table, _ = k.results_table()
    assert table.columns[:5] == ["Cluster", "Groups", "Mean ± SD (n)", "Test", "Paired by"]
    assert table.rows[0][4] == "Animal" and len(table.rows[0]) == len(table.columns)
    assert set(k.freq.content().notes) <= {(None, r) for r in freq.rows}
    k.test.setCurrentIndex(k.test.findData(c.RANK))
    assert {f.test for f in k.comparisons(pairwise=True)} == {"Friedman", "Wilcoxon signed-rank"}
    assert all("smallest p possible" in f.note for f in k.comparisons(pairwise=True)
               if f.test.startswith("Wilcoxon") and f.p == f.p)

    # a column that does not pair the samples: refused, with the reason, and no p on the chart
    k.paired_by.setCurrentIndex(k.paired_by.findData("Group"))
    assert "cannot be paired by Group" in k.pairing.text()
    assert "'Blood' is more than once in group 'Blood'" in k.pairing_problem()
    assert k.comparisons() == [] and k.results_table()[1] == [] and not k.freq.content().notes
    k.paired_by.setCurrentIndex(k.paired_by.findData("Animal"))
    kept = [m["Animal"] for m in ds.meta]
    for m in ds.meta:                      # each tissue from animals of its own
        m["Animal"] += m["Group"]
    k._draw_tables()
    assert "no subject is in every group" in k.pairing.text() and k.comparisons() == []
    for m, animal in zip(ds.meta, kept):
        m["Animal"] = animal
    k.paired_by.setCurrentIndex(0)
    assert k.pairing.text() == "" and {f.test for f in k.comparisons()} == {"Kruskal-Wallis"}
    k.test.setCurrentIndex(0)
    assert not k.paired_by.isEnabled()

    # the statistics table
    s = w.statistics
    s.tick(parameters=["CD4"], statistics=["Median"])
    s.compute()
    settle()
    d = CompareDialog(s.value_columns(), w.session.columns, w.session.group_by, s)
    assert d.paired_by.currentData() == "" and d.paired_by.findData("Animal") > 0
    assert d.test.findData(c.WELCH) > 0
    d.close()
    d = CompareDialog(s.value_columns(), w.session.columns, w.session.group_by, s, for_prism=True)
    assert d.layout().labelForField(d.paired_by) is None      # Prism does its own testing
    d.close()
    plain = s.comparison("CD4 · Median", ["Group"], pairwise=True)[1]
    assert [f.test for f in plain] == ["ANOVA", "Welch t", "Welch t", "Welch t"]
    table, found = s.comparison("CD4 · Median", ["Group"], pairwise=True, paired_by="Animal")
    assert [f.test for f in found] == ["RM ANOVA (Greenhouse-Geisser)"] + ["Paired t"] * 3
    assert table.columns[:5] == ["Population", "Groups", "Mean ± SD (n)", "Test", "Paired by"]
    rows = s.result
    median = {(r[rows.columns.index("Animal")], r[rows.columns.index("Group")]): r[-1]
              for r in rows.rows}
    blood, lung = ([median[a, g] for a in "1234"] for g in ("Blood", "Lung"))
    assert found[1].groups == ["Blood", "Lung"]
    assert found[1].p ==pytest.approx(st.ttest_rel(blood, lung).pvalue)
    assert s.comparison("CD4 · Median", ["Group"], kind=c.WELCH)[1][0].test == "Welch ANOVA"
    with pytest.raises(ValueError, match="cannot be paired by Sample: no subject is in every group"):
        s.comparison("CD4 · Median", ["Group"], paired_by="Sample")
    with pytest.raises(ValueError, match="cannot be paired by Group: subject 'Blood' is more than once"):
        s.comparison("CD4 · Median", ["Group"], paired_by="Group")
    # the other way about is a pairing too: the animals compared, tissue by tissue
    assert s.comparison("CD4 · Median", ["Animal"], paired_by="Group")[1][0].n == [3, 3, 3, 3]
    w.close()
