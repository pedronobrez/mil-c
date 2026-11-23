import math
import xml.etree.ElementTree as ET

import numpy as np

from milc import prism
from milc.cluster_stats import Table
from milc.prism import COLUMN, GROUPED, PrismTable
from milc.stats import StatsTable

NAN = float("nan")


def same(a, b) -> bool:
    """Nested lists equal, NaN equal to NaN."""
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    return a == b or (math.isnan(a) and math.isnan(b))


def column_table():
    return PrismTable("CD8 · Median", COLUMN, ["Control", "Treated", "KO"],
                      [[1.5, 2.25, 3.0], [4.0, NAN, 6.125, 7.0], [8.0]])


def grouped_table():
    # Control has 2 samples, Treated 3; one of Treated's has no value in cluster 2
    return PrismTable("Frequencies", GROUPED, ["Control", "Treated"],
                      [[[10.0, 12.5], [20.0, 21.0, 22.0]],
                       [[0.012, 0.5], [1.0, NAN, 3.0]]],
                      ["1", "2"])


# --- .pzfx ---------------------------------------------------------------

def test_round_trip_column(tmp_path):
    path = str(tmp_path / "c.pzfx")
    prism.write_pzfx(path, [column_table()])
    (back,) = prism.read_pzfx(path)
    assert (back.title, back.kind, back.columns, back.rows) == \
        ("CD8 · Median", COLUMN, ["Control", "Treated", "KO"], [])
    assert same(back.values, column_table().values)


def test_round_trip_grouped(tmp_path):
    path = str(tmp_path / "g.pzfx")
    prism.write_pzfx(path, [grouped_table()])
    (back,) = prism.read_pzfx(path)
    assert (back.title, back.kind, back.columns, back.rows) == \
        ("Frequencies", GROUPED, ["Control", "Treated"], ["1", "2"])
    assert same(back.values, grouped_table().values)


def test_trailing_blanks_are_padding(tmp_path):
    path = str(tmp_path / "p.pzfx")
    prism.write_pzfx(path, [PrismTable("t", COLUMN, ["A", "B"], [[1.0, NAN], [2.0, 3.0]])])
    assert same(prism.read_pzfx(path)[0].values, [[1.0], [2.0, 3.0]])


def test_structure_follows_the_pzfx_package(tmp_path):
    path = str(tmp_path / "s.pzfx")
    prism.write_pzfx(path, [column_table(), grouped_table()])
    root = ET.parse(path).getroot()
    assert root.tag == "GraphPadPrismFile" and root.get("PrismXMLVersion") == "5.00"
    assert [c.tag for c in root] == ["Created", "InfoSequence", "Info", "TableSequence",
                                     "Table", "Table"]
    assert root.find("Created/OriginalVersion").get("CreatedByProgram") == "GraphPad Prism"
    refs = root.findall("TableSequence/Ref")
    assert [r.get("ID") for r in refs] == ["Table0", "Table1"]
    assert [r.get("Selected") for r in refs] == ["1", None]

    one, two = root.findall("Table")
    assert one.attrib == {"ID": "Table0", "XFormat": "none", "YFormat": "replicates",
                          "Replicates": "1", "TableType": "OneWay",
                          "EVFormat": "AsteriskAfterNumber"}
    assert one.find("RowTitlesColumn") is None
    assert [c.tag for c in one] == ["Title", "YColumn", "YColumn", "YColumn"]
    treated = one.findall("YColumn")[1]
    assert treated.get("Subcolumns") == "1" and treated.find("Title").text == "Treated"
    # every column as long as the longest; a blank is an empty <d/>
    assert [d.text for d in treated.findall("Subcolumn/d")] == ["4.0", None, "6.125", "7.0"]
    assert [d.text for d in one.findall("YColumn")[2].findall("Subcolumn/d")] == \
        ["8.0", None, None, None]

    assert two.get("TableType") == "TwoWay" and two.get("Replicates") == "3"
    assert [c.tag for c in two] == ["Title", "RowTitlesColumn", "YColumn", "YColumn"]
    assert [d.text for d in two.findall("RowTitlesColumn/Subcolumn/d")] == ["1", "2"]
    control = two.find("YColumn")
    assert control.get("Subcolumns") == "3" and control.get("Width") == str(3 * 89)
    # one subcolumn per sample, one <d> per row; the third pads the smaller group
    assert [[d.text for d in s.findall("d")] for s in control.findall("Subcolumn")] == \
        [["10.0", "0.012"], ["12.5", "0.5"], [None, None]]
    # 0.012 needs four decimals to show three significant figures
    assert control.get("Decimals") == "4"


def test_titles_are_escaped_and_utf8(tmp_path):
    path = str(tmp_path / "e.pzfx")
    table = PrismTable("CD4 & CD8 <live> · Média", GROUPED, ["WT & KO", "<ctrl>"],
                       [[[1.0], [2.0]]], ["CD8 · Median"])
    prism.write_pzfx(path, [table])
    raw = open(path, "rb").read()
    assert raw.startswith(b"<?xml version='1.0' encoding='UTF-8'?>")
    text = raw.decode("utf-8")
    assert "CD4 &amp; CD8 &lt;live&gt; · Média" in text
    assert "<Title>WT &amp; KO</Title>" in text and "<d>CD8 · Median</d>" in text
    (back,) = prism.read_pzfx(path)
    assert (back.title, back.columns, back.rows) == (table.title, table.columns, table.rows)


def test_several_tables(tmp_path):
    path = str(tmp_path / "m.pzfx")
    tables = [column_table(), grouped_table(), PrismTable("empty", COLUMN, ["A"], [[]])]
    prism.write_pzfx(path, tables)
    back = prism.read_pzfx(path)
    assert [(t.title, t.kind) for t in back] == \
        [("CD8 · Median", COLUMN), ("Frequencies", GROUPED), ("empty", COLUMN)]
    assert all(same(b.values, t.values) for b, t in zip(back, tables))


# --- clipboard -----------------------------------------------------------

def test_clipboard_column():
    assert prism.to_clipboard(column_table()) == "\n".join([
        "1.5\t4\t8",
        "2.25\t\t",
        "3\t6.125\t",
        "\t7\t",
    ])


def test_clipboard_grouped():
    table = grouped_table()
    assert table.replicates() == 3
    assert prism.to_clipboard(table) == "\n".join([
        "1\t10\t12.5\t\t20\t21\t22",
        "2\t0.012\t0.5\t\t1\t\t3",
    ])


# --- builders ------------------------------------------------------------

def frequencies():
    #                 s1   s2   s3   s4   s5
    return Table(["1", "2"], ["s1", "s2", "s3", "s4", "s5"],
                 np.array([[60.0, 70.0, 30.0, 40.0, NAN],
                           [40.0, 30.0, 70.0, 60.0, NAN]]))


def test_from_frequencies():
    (t,) = prism.from_frequencies(frequencies(), ["KO", "WT", "KO", "WT", "KO"])
    assert (t.title, t.kind, t.rows, t.columns) == \
        ("Cluster frequencies", GROUPED, ["1", "2"], ["KO", "WT"])
    # s5 had no events: a blank subcolumn, not a zero
    assert same(t.values, [[[60.0, 30.0, NAN], [70.0, 40.0]],
                           [[40.0, 70.0, NAN], [30.0, 60.0]]])


def test_from_frequencies_by_facet():
    tables = prism.from_frequencies(frequencies(), ["KO", "WT", "KO", "WT", "KO"],
                                    facets=["14", "14", "3", "3", "3"], title="Freq")
    assert [t.title for t in tables] == ["Freq · 3", "Freq · 14"]      # numeric order
    assert all(t.columns == ["KO", "WT"] and t.rows == ["1", "2"] for t in tables)
    assert same(tables[0].values, [[[30.0, NAN], [40.0]], [[70.0, NAN], [60.0]]])
    assert same(tables[1].values, [[[60.0], [70.0]], [[40.0], [30.0]]])


def test_from_frequencies_keeps_a_group_absent_from_a_facet():
    tables = prism.from_frequencies(frequencies(), ["KO", "WT", "KO", "KO", "KO"],
                                    facets=["a", "a", "b", "b", "b"])
    assert tables[1].columns == ["KO", "WT"]
    assert same(tables[1].values[0], [[30.0, 40.0, NAN], []])


def stats_table():
    t = StatsTable(["Sample", "Group", "Day", "Population", "Path", "Count", "CD8 · Median"])
    for sample, group, day, t_cells, cd8 in [
            ("a.fcs", "WT", "7", 100.0, 10.0), ("b.fcs", "WT", "7", 110.0, None),
            ("c.fcs", "KO", "7", 200.0, 20.0), ("d.fcs", "KO", "14", 210.0, 21.0),
            ("e.fcs", "KO", "7", 220.0, 22.0)]:
        t.rows.append([sample, group, day, "T cells", "T cells", 1, t_cells])
        t.rows.append([sample, group, day, "CD8", "T cells/CD8", 1, cd8])
    return t


def test_from_stats_table():
    (t,) = prism.from_stats_table(stats_table(), "CD8 · Median", ["Group"])
    assert (t.title, t.kind, t.rows, t.columns) == \
        ("CD8 · Median", GROUPED, ["T cells", "CD8"], ["KO", "WT"])
    # b.fcs has no CD8 value: its subcolumn stays, blank, so samples line up
    assert same(t.values, [[[200.0, 210.0, 220.0], [100.0, 110.0]],
                           [[20.0, 21.0, 22.0], [10.0, NAN]]])


def test_from_stats_table_combined_groups_and_within():
    tables = prism.from_stats_table(stats_table(), "CD8 · Median", ["Group", "Day"])
    assert tables[0].columns == ["KO · 7", "KO · 14", "WT · 7"]      # day 7 before day 14
    assert same(tables[0].values[0], [[200.0, 220.0], [210.0], [100.0, 110.0]])

    tables = prism.from_stats_table(stats_table(), "CD8 · Median", ["Group"], within="Day",
                                    title="MFI")
    assert [t.title for t in tables] == ["MFI · 7", "MFI · 14"]
    assert all(t.columns == ["KO", "WT"] for t in tables)
    assert same(tables[0].values, [[[200.0, 220.0], [100.0, 110.0]], [[20.0, 22.0], [10.0, NAN]]])
    assert same(tables[1].values, [[[210.0], []], [[21.0], []]])


def test_from_stats_table_single_population_is_a_column_table():
    table = stats_table()
    table.rows = [r for r in table.rows if r[3] == "CD8"]
    (t,) = prism.from_stats_table(table, "CD8 · Median", ["Group"])
    assert (t.title, t.kind, t.rows, t.columns) == \
        ("CD8 · Median · CD8", COLUMN, [], ["KO", "WT"])
    assert same(t.values, [[20.0, 21.0, 22.0], [10.0]])


def test_from_stats_table_same_name_twice_uses_the_path():
    t = StatsTable(["Sample", "Group", "Population", "Path", "Count"])
    t.rows = [["a", "WT", "Live", "A/Live", 5], ["a", "WT", "Live", "B/Live", 7],
              ["a", "KO", "Live", "A/Live", 6]]
    (out,) = prism.from_stats_table(t, "Count", ["Group"])
    assert out.rows == ["A/Live", "B/Live"]
    # two samples named alike in different groups are still two samples
    assert same(out.values, [[[6.0], [5.0]], [[NAN], [7.0]]])


def test_builders_feed_the_writer(tmp_path):
    path = str(tmp_path / "b.pzfx")
    tables = (prism.from_frequencies(frequencies(), ["KO", "WT", "KO", "WT", "KO"])
              + prism.from_stats_table(stats_table(), "CD8 · Median", ["Group"], within="Day"))
    prism.write_pzfx(path, tables)
    back = prism.read_pzfx(path)
    assert [t.title for t in back] == ["Cluster frequencies", "CD8 · Median · 7",
                                       "CD8 · Median · 14"]
    assert same(back[0].values, [[[60.0, 30.0], [70.0, 40.0]], [[40.0, 70.0], [30.0, 60.0]]])
