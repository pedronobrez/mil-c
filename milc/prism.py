"""
Per-sample values in the shape GraphPad Prism wants, so that nobody has to
rearrange a table by hand before plotting it.

Prism does not take a long table (one row per sample and population). It
takes one of its own kinds:

  Column    one column per group, each sample a row in its group's column
  Grouped   one row per cluster or population, one column per group, and
            inside each column one subcolumn per sample (the replicates)

Groups rarely have the same number of samples. A grouped table has one
number of subcolumns for all its columns, so the smaller groups are padded
with blank cells; a blank is also what a missing value (NaN) becomes.

Two ways in: a `.pzfx` file, which Prism opens with every table already
made, and tab-separated text to paste into a table created in Prism.

The `.pzfx` written here follows the files of the R package `pzfx`
(github.com/Yue-Jiang/pzfx, `write_pzfx`), which Prism opens: the same
elements in the same order, the same attributes. One thing goes beyond it:
that package marks every table `TableType="OneWay"`, replicates or not,
while here a grouped table is `"TwoWay"`, which is what Prism itself
writes for one.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime, timezone

import numpy as np

from .cluster_stats import Table, ordered
from .stats import StatsTable

COLUMN = "column"
GROUPED = "grouped"

_TABLE_TYPE = {COLUMN: "OneWay", GROUPED: "TwoWay"}
_NAN = float("nan")


@dataclass
class PrismTable:
    """One Prism data table.

    `columns` are the groups. For a column table `values[group]` is one
    number per sample and `rows` is empty; for a grouped table
    `values[row][group]` is one number per sample, the same sample at the
    same position in every row. Lists may differ in length; NaN is a
    blank cell."""
    title: str
    kind: str
    columns: list[str]
    values: list
    rows: list[str] = field(default_factory=list)

    def replicates(self) -> int:
        """Samples in the largest group: the subcolumns of a grouped
        table, the data rows of a column table."""
        cells = self.values if self.kind == COLUMN else [c for row in self.values for c in row]
        return max((len(c) for c in cells), default=0)


def _blank(v) -> bool:
    return v is None or not math.isfinite(v)


def _padded(cell, n: int) -> list[float]:
    cell = [_NAN if _blank(v) else float(v) for v in cell]
    return cell + [_NAN] * (n - len(cell))


def _trimmed(cell: list[float]) -> list[float]:
    while cell and math.isnan(cell[-1]):
        cell = cell[:-1]
    return cell


def _subcolumns(table: PrismTable) -> list[list[list[float]]]:
    """`[group][subcolumn][row]`: the table as Prism stores it."""
    n = table.replicates()
    if table.kind == COLUMN:
        return [[_padded(cell, n)] for cell in table.values]
    rows = [[_padded(cell, n) for cell in row] for row in table.values]
    return [[[row[g][k] for row in rows] for k in range(n)]
            for g in range(len(table.columns))]


# --- .pzfx ---------------------------------------------------------------

def _decimals(table: PrismTable) -> int:
    """Decimal places Prism shows at first: enough for three significant
    figures of the smallest value, so that a cluster holding 0.012 % of a
    sample is not displayed as 0.00. The file keeps every digit either way."""
    flat = np.array([v for group in _subcolumns(table) for sub in group for v in sub])
    flat = np.abs(flat[np.isfinite(flat) & (flat != 0)])
    if flat.size == 0:
        return 2
    return int(np.clip(2 - np.floor(np.log10(flat.min())), 0, 6))


def _fill(parent: ET.Element, cells: list) -> None:
    sub = ET.SubElement(parent, "Subcolumn")
    for v in cells:
        d = ET.SubElement(sub, "d")
        if isinstance(v, str):
            d.text = v
        elif not math.isnan(v):
            d.text = repr(v)


def write_pzfx(path: str, tables: list[PrismTable]) -> None:
    """A Prism file holding `tables`, in that order."""
    root = ET.Element("GraphPadPrismFile", PrismXMLVersion="5.00")
    created = ET.SubElement(root, "Created")
    ET.SubElement(created, "OriginalVersion", CreatedByProgram="GraphPad Prism",
                  CreatedByVersion="6.0f.254", Login="",
                  DateTime=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00"))
    ET.SubElement(ET.SubElement(root, "InfoSequence"), "Ref", ID="Info0", Selected="1")
    info = ET.SubElement(root, "Info", ID="Info0")
    ET.SubElement(info, "Title").text = "Project info 1"
    ET.SubElement(ET.SubElement(info, "Notes"), "Font", Color="#000000", Face="Helvetica")
    sequence = ET.SubElement(root, "TableSequence")
    for i, table in enumerate(tables):
        ET.SubElement(sequence, "Ref", ID=f"Table{i}", **({"Selected": "1"} if i == 0 else {}))
    for i, table in enumerate(tables):
        n = 1 if table.kind == COLUMN else max(table.replicates(), 1)
        node = ET.SubElement(root, "Table", ID=f"Table{i}", XFormat="none",
                             YFormat="replicates", Replicates=str(n),
                             TableType=_TABLE_TYPE[table.kind],
                             EVFormat="AsteriskAfterNumber")
        ET.SubElement(node, "Title").text = table.title
        if table.kind == GROUPED:
            _fill(ET.SubElement(node, "RowTitlesColumn", Width="120"),
                  [str(r) for r in table.rows])
        decimals = str(_decimals(table))
        for name, group in zip(table.columns, _subcolumns(table)):
            column = ET.SubElement(node, "YColumn", Width=str(89 * n), Decimals=decimals,
                                   Subcolumns=str(n))
            ET.SubElement(column, "Title").text = str(name)
            # a grouped table without samples still needs its one subcolumn
            for cells in group or [[_NAN] * len(table.rows)]:
                _fill(column, cells)
    ET.indent(root)
    ET.ElementTree(root).write(path, encoding="UTF-8", xml_declaration=True)


def read_pzfx(path: str) -> list[PrismTable]:
    """The tables of a file written by `write_pzfx`. The blanks that pad a
    short group are not values and do not come back: every list ends at
    its last number."""
    root = ET.parse(path).getroot()
    for node in root.iter():
        node.tag = node.tag.rpartition("}")[2]    # Prism's own files are namespaced

    def text(node) -> str:
        return "".join(node.itertext()) if node is not None else ""

    def number(d) -> float:
        return float(text(d)) if text(d).strip() else _NAN

    out = []
    for node in root.iter("Table"):
        columns = node.findall("YColumn")
        subs = [[[number(d) for d in sub.findall("d")] for sub in col.findall("Subcolumn")]
                for col in columns]
        names = [text(col.find("Title")) for col in columns]
        title = text(node.find("Title"))
        if node.get("TableType") == "TwoWay":
            rows = [text(d) for d in node.findall("RowTitlesColumn/Subcolumn/d")]
            values = [[_trimmed([sub[r] for sub in group]) for group in subs]
                      for r in range(len(rows))]
            out.append(PrismTable(title, GROUPED, names, values, rows))
        else:
            out.append(PrismTable(title, COLUMN, names, [_trimmed(group[0]) for group in subs]))
    return out


# --- clipboard -----------------------------------------------------------

def to_clipboard(table: PrismTable) -> str:
    """Tab-separated text to paste into a Prism table of the same kind.

    Column: one line per data row, one cell per group. Paste with the
    cursor in the first row of column A.

    Grouped: one line per row — its title, then every subcolumn of the
    first group, of the second, and so on, each group padded to the
    largest. Create the table with `table.replicates()` replicate
    subcolumns and paste with the cursor in the row-title cell of row 1.

    There is no line of column titles. Prism pastes the block cell by cell
    from the cursor rightwards and, in its own words, "does not pay
    attention to the difference between subcolumns and data sets", so a
    title line could only be right for a table without subcolumns, and
    only if Prism were known to lift it into the title row — which its
    documentation does not say. The group names are in `table.columns`,
    in order, to be typed or pasted into the title row separately."""
    def cell(v: float) -> str:
        return "" if math.isnan(v) else f"{v:.6g}"

    groups = _subcolumns(table)
    lines = []
    for r in range(len(table.rows) if table.kind == GROUPED else table.replicates()):
        title = [str(table.rows[r])] if table.kind == GROUPED else []
        lines.append("\t".join(title + [cell(sub[r]) for group in groups for sub in group]))
    return "\n".join(lines)


# --- from the program's own tables ---------------------------------------

def from_frequencies(freq: Table, groups: list[str], facets: list[str] | None = None,
                     title: str = "Cluster frequencies") -> list[PrismTable]:
    """The cluster x sample frequencies as a grouped table: a row per
    cluster, a column per group, a subcolumn per sample. `groups[j]` — and
    `facets[j]`, when given — describe column j of `freq`. With facets
    (the day, the organ) there is one table per facet value, named after
    it, and every table has the same columns even where a group has no
    sample, so that they line up in Prism."""
    names = ordered(groups)

    def table(name: str, chosen: list[int]) -> PrismTable:
        members = [[j for j in chosen if groups[j] == g] for g in names]
        return PrismTable(name, GROUPED, [str(g) for g in names],
                          [[[float(row[j]) for j in m] for m in members] for row in freq.values],
                          list(freq.rows))

    if facets is None:
        return [table(title, list(range(len(groups))))]
    return [table(f"{title} · {f}", [j for j, x in enumerate(facets) if x == f])
            for f in ordered(facets)]


def from_stats_table(table: StatsTable, value: str, group_columns: list[str],
                     within: str | None = None, title: str | None = None) -> list[PrismTable]:
    """The statistics table (one row per sample and population) as Prism
    tables of `value`, a numeric column. A group is the combination of the
    `group_columns`, joined with " · " as in `compare.from_table`. With
    `within`, one table per value of that column, named after it.

    Several populations give grouped tables: a row per population, in the
    table's own order, a column per group, a subcolumn per sample — a
    sample keeps its subcolumn down the rows, blank where it has no value.
    A single population gives column tables instead, named after it and
    without the blanks: that is the kind Prism runs a t test or a one-way
    ANOVA on. `title` defaults to the name of the value."""
    cols = table.columns
    v, path, name = cols.index(value), cols.index("Path"), cols.index("Population")
    g = [cols.index(c) for c in group_columns]
    w = cols.index(within) if within else None
    # two populations may share a name; where they do, the path tells them apart
    pairs = {(r[name], r[path]) for r in table.rows}
    ambiguous = len({n for n, _ in pairs}) < len(pairs)

    populations: list[str] = []
    # facet → group → sample → population → value. A sample is everything
    # to the left of "Population": its name alone may not be unique.
    found: dict = {}
    for row in table.rows:
        who = str((row[path] or row[name]) if ambiguous else row[name])
        if who not in populations:
            populations.append(who)
        number = row[v]
        if isinstance(number, (int, float, np.integer, np.floating)):
            label = " · ".join(str(row[i]) or "—" for i in g)
            facet = str(row[w]) if w is not None else ""
            found.setdefault(facet, {}).setdefault(label, {}) \
                 .setdefault(tuple(row[:name]), {})[who] = float(number)

    base = title or value
    names = ordered(label for by_group in found.values() for label in by_group)
    out = []
    for facet in ordered(found):
        head = f"{base} · {facet}" if within else base
        samples = [list(found[facet].get(label, {}).values()) for label in names]
        if len(populations) == 1:
            cells = [[s[populations[0]] for s in group] for group in samples]
            out.append(PrismTable(f"{head} · {populations[0]}", COLUMN, names,
                                  [[x for x in cell if not _blank(x)] for cell in cells]))
        else:
            out.append(PrismTable(head, GROUPED, names,
                                  [[[s.get(p, _NAN) for s in group] for group in samples]
                                   for p in populations], list(populations)))
    return out
