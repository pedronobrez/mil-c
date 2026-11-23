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
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb

import numpy as np

from .stats import StatsTable

PARAMETRIC = "parametric"     # Welch t / ANOVA
RANK = "rank"                 # Mann-Whitney / Kruskal-Wallis
KINDS = {PARAMETRIC: "Welch t test / ANOVA", RANK: "Mann-Whitney / Kruskal-Wallis"}

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


def test_groups(samples: list[np.ndarray], kind: str = PARAMETRIC) -> tuple[str, float, float, str]:
    """(name of the test, statistic, p, note). NaN where it cannot be
    done, with the reason in the note: a p value is never invented."""
    from scipy import stats as st
    samples = [np.asarray(s, dtype=np.float64) for s in samples]
    samples = [s[np.isfinite(s)] for s in samples]
    two = len(samples) == 2
    name = (("Welch t" if two else "ANOVA") if kind == PARAMETRIC
            else ("Mann-Whitney" if two else "Kruskal-Wallis"))
    nan = float("nan")
    if len(samples) < 2:
        return name, nan, nan, "one group: nothing to compare"
    if min(s.size for s in samples) < 2:
        return name, nan, nan, "a group has fewer than 2 samples"
    if np.ptp(np.concatenate(samples)) == 0:
        return name, nan, nan, "every sample has the same value"
    note = ""
    try:
        if kind == PARAMETRIC and two:
            if all(np.ptp(s) == 0 for s in samples):
                return name, nan, nan, "no variation inside either group"
            r = st.ttest_ind(samples[0], samples[1], equal_var=False)
        elif kind == PARAMETRIC:
            r = st.f_oneway(*samples)
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
            correction: str = BH, pairwise: bool = False) -> list[Comparison]:
    """`values[key][group]` is one number per sample. One test per key —
    or, with `pairwise` and more than two groups, one per pair of groups
    as well — and one adjustment over all of them."""
    out: list[Comparison] = []
    for key, by_group in values.items():
        names = [g for g in ordered(by_group) if len(by_group[g])]
        sets = [names] + ([list(pair) for pair in combinations(names, 2)]
                          if pairwise and len(names) > 2 else [])
        for chosen in sets:
            arrays = [np.asarray(by_group[g], dtype=np.float64) for g in chosen]
            arrays = [a[np.isfinite(a)] for a in arrays]
            name, statistic, p, note = test_groups(arrays, kind)
            out.append(Comparison(
                key if isinstance(key, tuple) else (key,), list(chosen),
                [int(a.size) for a in arrays],
                [float(a.mean()) if a.size else float("nan") for a in arrays],
                [float(a.std(ddof=1)) if a.size > 1 else float("nan") for a in arrays],
                name, statistic, p, note=note))
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
             correction: str = BH) -> StatsTable:
    """Every comparison as a row: what was compared, each group's mean,
    SD and n, the test, and p before and after adjustment."""
    adjusted = "p" if correction == NONE else f"p adjusted ({correction})"
    table = StatsTable([*key_names, "Groups", "Mean ± SD (n)", "Test", "Statistic", "p",
                        adjusted, "", "Note"])
    for c in comparisons:
        summary = "; ".join(
            f"{g or '—'}: {m:.3g}" + (f" ± {s:.3g}" if np.isfinite(s) else "") + f" ({n})"
            for g, m, s, n in zip(c.groups, c.mean, c.sd, c.n))
        table.rows.append([*c.key, " vs ".join(g or "—" for g in c.groups), summary, c.test,
                           c.statistic, c.p, c.p_adjusted, stars(c.p_adjusted), c.note])
    return table


def from_table(table: StatsTable, value: str, group: list[str],
               within: str | None = None) -> tuple[dict, list[str]]:
    """The statistics table as values to compare: for each population —
    and each value of `within`, when given — the `value` of every sample,
    by group. `group` may be several columns: their combination is the
    group. Returns the values and the names of the key's parts."""
    cols = table.columns
    v, path, name = cols.index(value), cols.index("Path"), cols.index("Population")
    g = [cols.index(c) for c in group]
    w = cols.index(within) if within else None
    # two populations may share a name; where they do, the path tells them apart
    names = {(r[name], r[path]) for r in table.rows}
    ambiguous = len({n for n, _ in names}) < len(names)
    out: dict = {}
    for row in table.rows:
        number = row[v]
        if number is None or not isinstance(number, (int, float, np.integer, np.floating)):
            continue
        who = (row[path] or row[name]) if ambiguous else row[name]
        key = (row[w], who) if w is not None else (who,)
        label = " · ".join(str(row[i]) or "—" for i in g)
        out.setdefault(key, {}).setdefault(label, []).append(float(number))
    return out, ([within, "Population"] if within else ["Population"])
