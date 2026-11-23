"""
Is a difference between groups more than chance? One value per sample —
the share of a cluster, a median, a count — compared across the groups the
samples belong to.

The unit is the sample, never the cell: five animals are five
observations however many cells each gave. Testing cells against cells
would find every difference "significant" and say nothing about animals.

Two groups: Welch's t test, which does not assume equal variances, or the
Mann-Whitney U test. Three or more: one-way ANOVA or Kruskal-Wallis, and
optionally every pair. When many things are tested at once — every
cluster, every population — the p values are adjusted for it.

Where the variances of three or more groups differ, Welch's ANOVA does
for them what Welch's t does for two.

Where the same subject — an animal, a donor — was measured in every
group, the groups are conditions and the test is of each subject's
change: the paired t test or Wilcoxon's signed-rank test for two,
repeated-measures ANOVA or Friedman's test for more. Only subjects with
a value in every group can be used; the others are left out and counted.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb

import numpy as np

from .stats import StatsTable

PARAMETRIC = "parametric"     # Welch t / ANOVA
RANK = "rank"                 # Mann-Whitney / Kruskal-Wallis
WELCH = "welch"               # Welch t / Welch ANOVA
KINDS = {PARAMETRIC: "Welch t test / ANOVA", RANK: "Mann-Whitney / Kruskal-Wallis",
         WELCH: "Welch t test / Welch ANOVA"}

#: the test of each kind, for two groups and for more: between different
#: samples, and on the same subjects. Unequal variances between groups are
#: not a question within subjects, so there Welch's kind is the parametric.
INDEPENDENT = {PARAMETRIC: ("Welch t", "ANOVA"), WELCH: ("Welch t", "Welch ANOVA"),
               RANK: ("Mann-Whitney", "Kruskal-Wallis")}
REPEATED = {PARAMETRIC: ("Paired t", "RM ANOVA (Greenhouse-Geisser)"),
            WELCH: ("Paired t", "RM ANOVA (Greenhouse-Geisser)"),
            RANK: ("Wilcoxon signed-rank", "Friedman")}

BH = "Benjamini-Hochberg"
HOLM = "Holm"
BONFERRONI = "Bonferroni"
NONE = "None"
CORRECTIONS = [BH, HOLM, BONFERRONI, NONE]


@dataclass
class Comparison:
    #: what was tested: a cluster, a population, or (panel, cluster)
    key: tuple
    groups: list[str]
    n: list[int]
    mean: list[float]
    sd: list[float]
    test: str
    statistic: float
    p: float
    p_adjusted: float = float("nan")
    note: str = ""
    #: paired designs: the samples left out for want of a complete subject
    dropped: int = 0


def test_groups(samples: list[np.ndarray], kind: str = PARAMETRIC) -> tuple[str, float, float, str]:
    """(name of the test, statistic, p, note). NaN where it cannot be
    done, with the reason in the note: a p value is never invented."""
    from scipy import stats as st
    samples = [np.asarray(s, dtype=np.float64) for s in samples]
    samples = [s[np.isfinite(s)] for s in samples]
    two = len(samples) == 2
    name = INDEPENDENT.get(kind, INDEPENDENT[RANK])[0 if two else 1]
    nan = float("nan")
    if len(samples) < 2:
        return name, nan, nan, "one group: nothing to compare"
    if min(s.size for s in samples) < 2:
        return name, nan, nan, "a group has fewer than 2 samples"
    if np.ptp(np.concatenate(samples)) == 0:
        return name, nan, nan, "every sample has the same value"
    note = ""
    try:
        if kind in (PARAMETRIC, WELCH) and two:
            if all(np.ptp(s) == 0 for s in samples):
                return name, nan, nan, "no variation inside either group"
            r = st.ttest_ind(samples[0], samples[1], equal_var=False)
        elif kind == PARAMETRIC:
            r = st.f_oneway(*samples)
        elif kind == WELCH:
            if any(np.ptp(s) == 0 for s in samples):
                return name, nan, nan, "no variation inside a group"
            statistic, _, _, p = welch_anova(samples)
            return name, statistic, p, note
        elif two:
            r = st.mannwhitneyu(samples[0], samples[1], alternative="two-sided")
            # with few samples a rank test cannot give a small p at all:
            # 3 against 3 never goes below 0.1, whatever the data
            floor = 2.0 / comb(samples[0].size + samples[1].size, samples[0].size)
            if floor > 0.05:
                note = f"smallest p possible with these n: {floor:.2g}"
        else:
            r = st.kruskal(*samples)
    except ValueError as exc:
        return name, nan, nan, str(exc)
    return name, float(r.statistic), float(r.pvalue), note


def welch_anova(samples: list[np.ndarray]) -> tuple[float, float, float, float]:
    """(F, df between, df within, p) of Welch's one-way ANOVA (Welch 1951):
    each group weighs n / variance, so one that scatters more counts for
    less, and the denominator's degrees of freedom shrink to match. With
    two groups it is Welch's t test, F = t²."""
    from scipy import stats as st
    k = len(samples)
    n = np.array([s.size for s in samples], dtype=np.float64)
    mean = np.array([s.mean() for s in samples])
    w = n / np.array([s.var(ddof=1) for s in samples])
    grand = (w * mean).sum() / w.sum()
    lack = (((1 - w / w.sum()) ** 2) / (n - 1)).sum()
    f = ((w * (mean - grand) ** 2).sum() / (k - 1)) / (1 + 2 * (k - 2) / (k * k - 1) * lack)
    df = (k * k - 1) / (3 * lack)
    return float(f), float(k - 1), float(df), float(st.f.sf(f, k - 1, df))


def rm_anova(data: np.ndarray) -> tuple[float, float, float, float, float, float]:
    """(F, df groups, df error, ε, p, p corrected) of a one-way
    repeated-measures ANOVA; `data[i, j]` is subject i in group j.

    What differs between subjects is taken out before the groups are
    compared. The F test then assumes sphericity — that the differences
    between any two groups scatter alike — and gives too small a p where
    they do not. Greenhouse and Geisser's ε measures the departure from
    the covariance of the groups, from 1 where it holds down to
    1 / (k - 1), and both degrees of freedom are multiplied by it."""
    from scipy import stats as st
    data = np.asarray(data, dtype=np.float64)
    n, k = data.shape
    grand = data.mean()
    groups = n * ((data.mean(axis=0) - grand) ** 2).sum()
    subjects = k * ((data.mean(axis=1) - grand) ** 2).sum()
    error = ((data - grand) ** 2).sum() - groups - subjects
    df1, df2 = k - 1, (k - 1) * (n - 1)
    f = (groups / df1) / (error / df2)
    centre = np.eye(k) - 1.0 / k
    s = centre @ np.cov(data, rowvar=False) @ centre
    epsilon = float(np.trace(s) ** 2 / (df1 * (s * s).sum()))
    return (float(f), float(df1), float(df2), epsilon, float(st.f.sf(f, df1, df2)),
            float(st.f.sf(f, epsilon * df1, epsilon * df2)))


def test_paired(columns: list[np.ndarray], kind: str = PARAMETRIC) -> tuple[str, float, float, str]:
    """As `test_groups`, for the same subjects in every group:
    `columns[j][i]` is subject i in group j, and no value is missing."""
    from scipy import stats as st
    columns = [np.asarray(c, dtype=np.float64) for c in columns]
    two = len(columns) == 2
    name = REPEATED.get(kind, REPEATED[RANK])[0 if two else 1]
    nan = float("nan")
    if len(columns) < 2:
        return name, nan, nan, "one group: nothing to compare"
    n = columns[0].size
    if n < 2:
        return name, nan, nan, "fewer than 2 subjects in every group"
    if np.ptp(np.concatenate(columns)) == 0:
        return name, nan, nan, "every sample has the same value"
    note = ""
    try:
        if two and kind in (PARAMETRIC, WELCH):
            if np.ptp(columns[0] - columns[1]) == 0:
                return name, nan, nan, "every subject changes by the same amount"
            r = st.ttest_rel(columns[0], columns[1])
        elif two:
            # a subject that does not change says nothing of the direction
            # and is left out of the ranks, as Wilcoxon had it
            changed = int(np.count_nonzero(columns[0] - columns[1]))
            if changed == 0:
                return name, nan, nan, "no subject differs between the two"
            r = st.wilcoxon(columns[0], columns[1])
            floor = 2.0 / 2 ** changed
            if floor > 0.05:
                note = f"smallest p possible with these n: {floor:.2g}"
        elif kind in (PARAMETRIC, WELCH):
            data = np.column_stack(columns)
            rest = data - data.mean(axis=0) - data.mean(axis=1, keepdims=True) + data.mean()
            if np.allclose(rest, 0, atol=1e-12 * np.abs(data).max()):
                return name, nan, nan, "no variation left after subjects and groups"
            f, _, _, epsilon, uncorrected, p = rm_anova(data)
            return name, f, p, f"ε = {epsilon:.3g}; p without the correction: {uncorrected:.2g}"
        else:
            r = st.friedmanchisquare(*columns)
            # scipy reads the statistic, n (k - 1) at most, off a chi-square
            floor = float(st.chi2.sf(n * (len(columns) - 1), len(columns) - 1))
            if floor > 0.05:
                note = f"smallest p possible with these n: {floor:.2g}"
    except ValueError as exc:
        return name, nan, nan, str(exc)
    return name, float(r.statistic), float(r.pvalue), note


def _by_subject(values, ids) -> tuple[dict[str, float], str]:
    """One group as subject → value, and the subject found twice, if one
    is. A sample with no value or no subject cannot be paired."""
    out: dict[str, float] = {}
    for value, subject in zip(values, ids):
        subject = "" if subject is None else str(subject).strip()
        if not subject or not np.isfinite(value):
            continue
        if subject in out:
            return out, subject
        out[subject] = float(value)
    return out, ""


def _complete(by_group: dict, ids: dict, chosen: list[str]) -> tuple[list[np.ndarray], int, str]:
    """(the values of the subjects found in every group of `chosen`, in
    the same order in each; how many samples that left out; what is wrong
    with the pairing, when it cannot be done at all)."""
    found = []
    for g in chosen:
        mapped, twice = _by_subject(by_group[g], ids[g])
        if twice:
            return [np.empty(0)] * len(chosen), 0, f"subject {twice!r} is more than once in group {g or '—'!r}"
        found.append(mapped)
    whole = [s for s in found[0] if all(s in m for m in found[1:])] if found else []
    arrays = [np.array([m[s] for s in whole], dtype=np.float64) for m in found]
    return arrays, sum(len(m) for m in found) - len(whole) * len(found), ""


def pairing_problem(values: dict, subjects: dict) -> str:
    """Why `values` cannot be compared as paired by `subjects` — a subject
    twice in a group, or none found in every group — or "" where it can."""
    shared = False
    for key, by_group in values.items():
        names = [g for g in by_group if len(by_group[g])]
        arrays, _, problem = _complete(by_group, subjects[key], names)
        if problem:
            return problem
        shared = shared or bool(len(names) > 1 and arrays[0].size)
    return "" if shared else "no subject is in every group"


def adjust(p, method: str = BH) -> np.ndarray:
    """The p values adjusted for having been many. Those that are NaN —
    tests that could not be done — stay NaN and do not count as tests."""
    p = np.asarray(p, dtype=np.float64)
    out = np.full(p.shape, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    m = ok.size
    if m == 0 or method == NONE:
        out[ok] = p[ok]
        return out
    values = p[ok]
    order = np.argsort(values)
    ranked = values[order]
    if method == BONFERRONI:
        adjusted = np.minimum(ranked * m, 1.0)
    elif method == HOLM:
        adjusted = np.minimum(np.maximum.accumulate(ranked * (m - np.arange(m))), 1.0)
    elif method == BH:
        adjusted = np.minimum(np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1], 1.0)
    else:
        raise ValueError(f"unknown correction {method!r}")
    out[ok[order]] = adjusted
    return out


def ordered(names) -> list[str]:
    names = list(dict.fromkeys(names))
    try:
        return sorted(names, key=float)
    except (TypeError, ValueError):
        return sorted(names)


def compare(values: dict[tuple, dict[str, list[float]]], kind: str = PARAMETRIC,
            correction: str = BH, pairwise: bool = False,
            subjects: dict[tuple, dict[str, list]] | None = None) -> list[Comparison]:
    """`values[key][group]` is one number per sample. One test per key —
    or, with `pairwise` and more than two groups, one per pair of groups
    as well — and one adjustment over all of them.

    With `subjects` — `subjects[key][group]` says whose each value of
    `values[key][group]` is — the design is paired: the tests are those of
    `test_paired`, on the subjects found in every group being compared."""
    out: list[Comparison] = []
    for key, by_group in values.items():
        names = [g for g in ordered(by_group) if len(by_group[g])]
        sets = [names] + ([list(pair) for pair in combinations(names, 2)]
                          if pairwise and len(names) > 2 else [])
        for chosen in sets:
            arrays = [np.asarray(by_group[g], dtype=np.float64) for g in chosen]
            arrays = [a[np.isfinite(a)] for a in arrays]
            dropped = 0
            if subjects is None:
                name, statistic, p, note = test_groups(arrays, kind)
            else:
                whole, dropped, wrong = _complete(by_group, subjects[key], chosen)
                name, statistic, p, note = test_paired(whole, kind)
                if wrong:   # no test is done, and the groups are described as they are
                    statistic, p, note = float("nan"), float("nan"), wrong
                else:
                    arrays = whole
                if dropped:
                    note = "; ".join(t for t in (note, f"{dropped} sample(s) left out: their "
                                                 "subject is not in every group") if t)
            out.append(Comparison(
                key if isinstance(key, tuple) else (key,), list(chosen),
                [int(a.size) for a in arrays],
                [float(a.mean()) if a.size else float("nan") for a in arrays],
                [float(a.std(ddof=1)) if a.size > 1 else float("nan") for a in arrays],
                name, statistic, p, note=note, dropped=dropped))
    for c, adjusted in zip(out, adjust([c.p for c in out], correction)):
        c.p_adjusted = float(adjusted)
    return out


def stars(p: float) -> str:
    if not np.isfinite(p):
        return ""
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 0.01 else \
        "*" if p < 0.05 else "ns"


def label(p: float) -> str:
    """A p value as it is written on a plot."""
    if not np.isfinite(p):
        return ""
    return "p < 0.0001" if p < 1e-4 else f"p = {p:.2g}" if p < 0.01 else f"p = {p:.2f}"


def to_table(comparisons: list[Comparison], key_names: list[str],
             correction: str = BH, paired_by: str | None = None) -> StatsTable:
    """Every comparison as a row: what was compared, each group's mean,
    SD and n, the test, and p before and after adjustment. Of a paired
    design, `paired_by` — what the subjects are — is a column as well, so
    the table says what it is wherever it is taken."""
    adjusted = "p" if correction == NONE else f"p adjusted ({correction})"
    paired = [paired_by] if paired_by else []
    table = StatsTable([*key_names, "Groups", "Mean ± SD (n)", "Test",
                        *(["Paired by"] if paired else []), "Statistic", "p", adjusted, "", "Note"])
    for c in comparisons:
        summary = "; ".join(
            f"{g or '—'}: {m:.3g}" + (f" ± {s:.3g}" if np.isfinite(s) else "") + f" ({n})"
            for g, m, s, n in zip(c.groups, c.mean, c.sd, c.n))
        table.rows.append([*c.key, " vs ".join(g or "—" for g in c.groups), summary, c.test,
                           *paired, c.statistic, c.p, c.p_adjusted, stars(c.p_adjusted), c.note])
    return table


def from_table(table: StatsTable, value: str, group: list[str],
               within: str | None = None, paired_by: str | None = None) -> tuple:
    """The statistics table as values to compare: for each population —
    and each value of `within`, when given — the `value` of every sample,
    by group. `group` may be several columns: their combination is the
    group. Returns the values and the names of the key's parts — and,
    with `paired_by`, what that column holds for each value: the
    `subjects` of `compare`."""
    cols = table.columns
    v, path, name = cols.index(value), cols.index("Path"), cols.index("Population")
    g = [cols.index(c) for c in group]
    w = cols.index(within) if within else None
    s = cols.index(paired_by) if paired_by else None
    # two populations may share a name; where they do, the path tells them apart
    names = {(r[name], r[path]) for r in table.rows}
    ambiguous = len({n for n, _ in names}) < len(names)
    out: dict = {}
    subjects: dict = {}
    for row in table.rows:
        number = row[v]
        if number is None or not isinstance(number, (int, float, np.integer, np.floating)):
            continue
        who = (row[path] or row[name]) if ambiguous else row[name]
        key = (row[w], who) if w is not None else (who,)
        label = " · ".join(str(row[i]) or "—" for i in g)
        out.setdefault(key, {}).setdefault(label, []).append(float(number))
        if s is not None:
            subjects.setdefault(key, {}).setdefault(label, []).append(row[s])
    names = [within, "Population"] if within else ["Population"]
    return (out, names, subjects) if paired_by else (out, names)
