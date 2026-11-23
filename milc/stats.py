"""
The statistics table: for every population of every sample, how many
events and what share of its parent, and — for the markers asked for — the
figures a cytometrist reports.

Every statistic is computed on the values as the file stores them, after
compensation and before any display scale: a median fluorescence is a
fluorescence, not a position on a biexponential axis. The definitions are
FlowJo's, so a number here can be set beside one from there:

  Median, Mean, SD       the usual ones (SD over n - 1)
  Geometric mean         exp(mean(ln x)) over the events with x > 0; the
                         others have no logarithm and are left out
  CV                     100 x SD / Mean
  Robust SD              half the distance between the 15.87th and the
                         84.13th percentile: the SD of a normal
                         distribution with that spread, deaf to outliers
  Robust CV              100 x Robust SD / Median
  MAD                    median of the absolute deviations from the median
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .gating import GateNode


def _geometric_mean(x: np.ndarray) -> float:
    positive = x[x > 0]
    return float(np.exp(np.log(positive).mean())) if positive.size else float("nan")


def _sd(x: np.ndarray) -> float:
    return float(x.std(ddof=1)) if x.size > 1 else float("nan")


def _robust_sd(x: np.ndarray) -> float:
    lo, hi = np.percentile(x, [15.87, 84.13])
    return float((hi - lo) / 2)


def _ratio(a: float, b: float) -> float:
    return 100.0 * a / b if b not in (0.0,) and np.isfinite(b) else float("nan")


#: name → function of the events' values. In the order they are offered.
STATISTICS = {
    "Median": lambda x: float(np.median(x)),
    "Mean": lambda x: float(x.mean()),
    "Geometric mean": _geometric_mean,
    "SD": _sd,
    "CV": lambda x: _ratio(_sd(x), float(x.mean())),
    "Robust SD": _robust_sd,
    "Robust CV": lambda x: _ratio(_robust_sd(x), float(np.median(x))),
    "MAD": lambda x: float(np.median(np.abs(x - np.median(x)))),
}

#: the columns every row has, whatever markers are asked for
COUNTS = ["Count", "% of parent", "% of grandparent", "% of total"]


def compute(name: str, values: np.ndarray) -> float:
    """One statistic of one set of values; NaN where there are no events."""
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return float("nan")
    return STATISTICS[name](values)


@dataclass
class StatsTable:
    """Rows of mixed text and numbers under named columns."""
    columns: list[str]
    rows: list[list] = field(default_factory=list)

    def column(self, name: str) -> list:
        i = self.columns.index(name)
        return [r[i] for r in self.rows]

    def to_csv(self, path: str) -> None:
        def cell(v) -> str:
            if v is None or (isinstance(v, float) and np.isnan(v)):
                return ""
            if isinstance(v, (float, np.floating)):
                return f"{v:.6g}"
            text = str(v)
            return '"' + text.replace('"', '""') + '"' if any(c in text for c in ',"\n') else text
        with open(path, "w", encoding="utf-8") as f:
            f.write(",".join(cell(c) for c in self.columns) + "\n")
            for row in self.rows:
                f.write(",".join(cell(v) for v in row) + "\n")

    def to_text(self) -> str:
        """Tab-separated, for the clipboard: pastes into Excel and Prism."""
        def cell(v) -> str:
            if v is None or (isinstance(v, float) and np.isnan(v)):
                return ""
            return f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v)
        return "\n".join("\t".join(cell(v) for v in row)
                         for row in [self.columns] + self.rows)


def population_rows(sample, populations: list[str] | None = None,
                    channels: list[tuple[str, str]] | None = None,
                    statistics: list[str] | None = None) -> list[list]:
    """One row per population of `sample`: its path, name, the four counts,
    then one value per (channel, statistic) in that order.

    `populations` are gate paths ("" is every event); None is all of the
    sample's. `channels` are (label, column name) pairs, the column name
    being what the sample answers to — `Comp-PE-A` for a compensated one."""
    masks = sample.masks()
    tree: GateNode = sample.gates
    total = len(sample)
    paths = populations if populations is not None else [n.path for n in tree.walk()]
    rows = []
    for path in paths:
        if path not in masks:
            continue
        node = tree.find(path) if path else tree
        mask = masks[path]
        n = int(mask.sum())
        parent = node.parent
        grand = parent.parent if parent is not None else None
        row = [path, node.name if path else "All events", n,
               _ratio(n, int(masks[parent.path].sum())) if parent is not None else float("nan"),
               _ratio(n, int(masks[grand.path].sum())) if grand is not None else float("nan"),
               _ratio(n, total)]
        for _label, column in channels or []:
            values = sample.column(column)[mask]
            row += [compute(s, values) for s in statistics or []]
        rows.append(row)
    return rows


def headers(annotations: list[str], channels: list[tuple[str, str]] | None,
            statistics: list[str] | None) -> list[str]:
    return (["Sample", *annotations, "Population", "Path", *COUNTS]
            + [f"{label} · {s}" for label, _ in channels or [] for s in statistics or []])
