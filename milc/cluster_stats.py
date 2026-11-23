"""
What the clusters are and where they are found: the numbers behind the
heatmap, the bars and the comparison between groups.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Table:
    rows: list[str]
    columns: list[str]
    values: np.ndarray

    def to_csv(self, path: str, corner: str = "") -> None:
        with open(path, "w", encoding="utf-8") as f:
            f.write(",".join([corner] + [_q(c) for c in self.columns]) + "\n")
            for name, row in zip(self.rows, self.values):
                f.write(",".join([_q(name)] + [f"{v:.6g}" for v in row]) + "\n")


def _q(text: str) -> str:
    text = str(text)
    return '"' + text.replace('"', '""') + '"' if any(c in text for c in ',"\n') else text


def ordered(names) -> list[str]:
    """Distinct names in the order a person expects: day 3, 7, 14 — not
    14, 3, 7."""
    names = list(dict.fromkeys(names))
    try:
        return sorted(names, key=float)
    except (TypeError, ValueError):
        return sorted(names)


def clusters_of(labels: np.ndarray) -> list[int]:
    return [int(c) for c in np.unique(labels)]


def expression(X: np.ndarray, labels: np.ndarray, names: list[str],
               stat: str = "median") -> Table:
    """Cluster x marker: the median (or mean) of each marker in each
    cluster, in whatever scale `X` is in."""
    fn = {"median": np.median, "mean": np.mean}[stat]
    ids = clusters_of(labels)
    values = np.array([fn(X[labels == c], axis=0) for c in ids])
    return Table([str(c) for c in ids], list(names), values)


def counts(labels: np.ndarray, sample_index: np.ndarray,
           samples: list[str]) -> Table:
    """Cluster x sample: events."""
    ids = clusters_of(labels)
    pos = {c: i for i, c in enumerate(ids)}
    out = np.zeros((len(ids), len(samples)))
    np.add.at(out, (np.array([pos[int(c)] for c in labels]), sample_index), 1)
    return Table([str(c) for c in ids], list(samples), out)


def frequencies(labels: np.ndarray, sample_index: np.ndarray,
                samples: list[str]) -> Table:
    """Cluster x sample: percent of the sample's events in the cluster.
    Each column sums to 100; a sample with no events is NaN, not zero."""
    t = counts(labels, sample_index, samples)
    totals = t.values.sum(axis=0)
    pct = np.full_like(t.values, np.nan)
    np.divide(100 * t.values, totals, out=pct, where=totals > 0)
    return Table(t.rows, t.columns, pct)


def by_group(freq: Table, groups: list[str]) -> tuple[Table, Table, Table]:
    """Cluster x group: the mean, the standard deviation (n-1) and the
    number of samples behind them. `groups` names each column of `freq`.
    A group of one sample has no standard deviation: NaN."""
    names = ordered(groups)
    shape = (len(freq.rows), len(names))
    mean, sd, n = np.full(shape, np.nan), np.full(shape, np.nan), np.zeros(shape)
    for j, g in enumerate(names):
        cols = freq.values[:, [i for i, x in enumerate(groups) if x == g]]
        cols = cols[:, ~np.isnan(cols).all(axis=0)]
        n[:, j] = cols.shape[1]
        if cols.shape[1]:
            mean[:, j] = cols.mean(axis=1)
        if cols.shape[1] > 1:
            sd[:, j] = cols.std(axis=1, ddof=1)
    return (Table(freq.rows, names, mean), Table(freq.rows, names, sd),
            Table(freq.rows, names, n))


def scaled(table: Table, per_column: bool = True) -> np.ndarray:
    """The table's values from 0 to 1 for colouring: per marker, so a dim
    marker shows its own differences, or over the whole table."""
    v = table.values.astype(np.float64)
    if per_column:
        lo, hi = np.nanmin(v, axis=0), np.nanmax(v, axis=0)
    else:
        lo, hi = np.nanmin(v), np.nanmax(v)
    span = np.where(hi > lo, hi - lo, 1.0)
    return (v - lo) / span


def selection(mask: np.ndarray, sample_index: np.ndarray, samples: list[str],
              groups: list[str]) -> dict:
    """Any set of events — a cluster, a lasso — as numbers: how many, what
    share of everything, and what share of each sample, by group."""
    n = int(mask.sum())
    per_sample = np.bincount(sample_index[mask], minlength=len(samples))
    totals = np.bincount(sample_index, minlength=len(samples))
    pct = np.full(len(samples), np.nan)
    np.divide(100.0 * per_sample, totals, out=pct, where=totals > 0)
    by_group = {}
    for g in ordered(groups):
        vals = np.array([pct[i] for i, x in enumerate(groups)
                         if x == g and not np.isnan(pct[i])])
        if vals.size:
            by_group[g] = (float(vals.mean()),
                           float(vals.std(ddof=1)) if vals.size > 1 else float("nan"),
                           int(vals.size))
    return {"n": n, "share": n / mask.size if mask.size else 0.0,
            "percent_by_sample": pct, "by_group": by_group}
