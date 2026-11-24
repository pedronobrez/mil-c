---
title: Tables, exports and Prism
---

MIL-C writes its numbers out in four ways: the statistics table as CSV or on the clipboard, the cluster tables as CSV, any column as a GraphPad Prism file, and the events of a population as FCS or CSV files. This page says what each one holds.

Every CSV the program writes is UTF-8 with a byte-order mark. Without the mark, Excel on Windows reads a CSV in the system's own code page and `CD8 · Median` arrives garbled; with it, Excel, Prism, R and pandas all read it correctly. Numbers are written with six significant figures; a missing value is an empty cell.

## The statistics table

From [[statistics-tab|the Statistics tab]]:

- **Export CSV…** writes the table as shown, every column, to `statistics.csv` or a name you choose. The folder is remembered.
- **Copy** puts the same table on the clipboard, tab-separated with its header line, to paste into Excel or Prism.

The table is long: one row per sample and population, with the sample columns beside the name. That is the shape R, Python and Excel's pivot tables want. Prism wants another shape, which **Prism…** makes.

The [[command-line]] gives the same file with `milc stats`, and a Prism file with `--prism`.

## Prism files

Prism does not take a long table. It takes one of its own kinds:

- **Column table** — one column per group, each sample a row in its group's column. This is the kind Prism runs a t test or a one-way ANOVA on.
- **Grouped table** — one row per population or cluster, one column per group, and inside each column one subcolumn per sample (the replicates).

Groups rarely have the same number of samples. A grouped table has one number of subcolumns for all its columns, so the smaller groups are padded with blank cells; a missing value is blank too. Each table opens with enough decimal places shown for three significant figures of its smallest value, so a cluster holding 0.012% of a sample is not displayed as 0.00; the file keeps every digit either way.

### From the Statistics tab

**Prism…** opens a dialog titled "A table for Prism":

| Field | What it does |
|---|---|
| **Values** | The column of the table to put in the file. Starts on % of parent |
| **Columns** | The sample columns that make Prism's columns. Several ticked: each combination is a column |
| **A table for each** | Optional: one table per value of this column (the day, say), named after it |

**Save…** asks for the file (`statistics.pzfx` by default). With several populations in the table, each table is grouped: a row per population, a column per group, a subcolumn per sample; a sample keeps its subcolumn down the rows, blank where it has no value. With a single population, each table is a column table named after the value and the population, without blanks. The statistics table must have been computed first. Prism does the testing itself, so the dialog has no test.

### From the Clusters tab

**Prism…** in [[clusters-tab|the Clusters tab]] writes the share of each cluster in each sample as grouped tables titled "Cluster frequency (% of sample)": clusters down, groups across, a subcolumn per sample. With a **within each** column there is one table per value, named after it, and every table has the same columns even where a group has no sample, so they line up in Prism. The default name is `cluster_frequency.pzfx`.

### What was checked

The `.pzfx` follows the files of the R package `pzfx` (the same elements in the same order), except that a grouped table is marked as a two-way table, as Prism itself writes one. The files round-trip through MIL-C's own reader. They have been opened in GraphPad Prism by hand before the 1.0 release.

The order of groups, here and in every other table, is the one set in [[sample-order]].

## The cluster tables

**File ▸ Export tables…** asks for a folder and writes two files about the FlowSOM clusters:

| File | Rows | Columns | Values |
|---|---|---|---|
| `cluster_expression.csv` | Clusters, by name | The markers FlowSOM ran on | The median of each marker in each cluster, over all the concatenated events, on the scale the clustering ran on (0 to 1023 for the biexponential) |
| `cluster_frequency.csv` | Clusters, by name | Samples | Percent of the sample's concatenated events in the cluster; each column sums to 100 |

The clusters are as the Clusters tab shows them: merged and renamed. These tables are FlowSOM's; when only PhenoGraph has been run, the menu answers "Run FlowSOM first". The same two tables come from `milc analyze` as `*_mfi.csv` and `*_freq.csv`.

## The events of a population

**Export…** in [[gates-tab|the Gates tab]] writes the events of the population shown, one file per sample. Its dialog, "Export *population*":

| Field | What it does |
|---|---|
| **Format** | **FCS 3.1** or **CSV** |
| **Every sample ticked in Samples** | Ticked: every sample ticked in the Samples tab. Unticked: only the samples selected in Gates |
| **Compensated values** | Ticked (the default): fluorescence as this program works with it, in columns named `Comp-…`, so nothing downstream compensates them again. Unticked: the events as acquired; an FCS file then carries the matrix in use as its `$SPILLOVER`, for whatever reads it next to apply |
| **Into** | The folder: by default a folder named `exported` inside the data folder; **Choose…** to pick another |

Each file is named `sample_population.fcs` (or `.csv`), for example `sample_01_CD4.fcs`. A CSV has one header line with the stain names (the channel name where there is no stain) and one line per event. An FCS file keeps the original's `$CYT`, `$DATE`, `$BTIM`, `$ETIM` and `TUBE NAME`, and adds `MILC POPULATION` and `MILC SOURCE` (the file it came from). A sample without that population is skipped.

The gates themselves, rather than their events, go out by **File ▸ Export gates…**: see [[workspaces-and-gating-ml]].
