---
title: The Statistics tab
---

The Statistics tab is a table: one row for each population of each sample, with the counts and percentages every row has and, for the parameters you ask for, the statistics a cytometrist reports. The samples' own columns from [[samples-tab|the Samples tab]] (Group, Day, Donor…) sit beside the sample's name, so the table can be sorted here and compared by group here or anywhere else.

It reads the individual files with their gates, not the concatenated events: every event of the population is counted, not only the ones drawn for [[concatenation]].

## Choosing what to compute

The panel on the left has three lists and the buttons.

- **Populations** — the gating tree, with a tick box on each population. Only populations that every gated sample has are listed, because each needs a row in every sample. A gate drawn since the last visit comes in ticked; one you unticked stays unticked. "All events" is the root, the ungated file.
- **Parameters** — every parameter of the samples, by its stain name where there is one. Fluorescence is taken compensated and in the file's own units, as FlowJo reports it: the column read is the `Comp-` one wherever the parameter is compensated (see [[compensation]]). With [[spectral-unmixing]] applied, the unmixed fluorochromes are listed instead of the detectors they replace.
- **Statistics** — tick any of the eight below. What you tick is remembered for next time; Median is ticked the first time.
- **Only the samples ticked in Samples** — ticked by default. Untick it to have a row for every usable sample.
- **Compute** — builds the table, in the background; the status bar shows which sample is being read.

Count and the three percentages are always computed, even with no parameter ticked.

## The statistics

Every statistic is computed on the events' values as the file stores them, after compensation and before any display scale: a median fluorescence is a fluorescence, not a position on a biexponential axis. The definitions are FlowJo's, so a number here can be set beside one from there.

| Statistic | Definition |
|---|---|
| Median | The middle value of the population's events |
| Mean | The arithmetic mean |
| Geometric mean | exp(mean(ln x)) over the events with x > 0. Events at or below zero have no logarithm and are left out |
| SD | Standard deviation, over n − 1 |
| CV | 100 × SD / Mean |
| Robust SD | Half the distance between the 15.87th and the 84.13th percentile: the SD of a normal distribution with that spread, deaf to outliers |
| Robust CV | 100 × Robust SD / Median |
| MAD | Median of the absolute deviations from the median |

A population with no events in a sample gives an empty cell, not a zero. A CV whose mean (or median) is zero is empty too.

## Columns

| Column | Meaning |
|---|---|
| Sample | The sample's name |
| Your columns | One for each column of the Samples table, as typed there |
| Population | The population's name; "All events" for the root |
| Path | The population's full gate path. Shown only when two populations share a name |
| Count | Events in the population |
| % of parent | Count over the parent population's count, × 100 |
| % of grandparent | Count over the parent's parent, × 100 (empty for the first level of gates) |
| % of total | Count over every event of the file, × 100 |
| *Marker* · *Statistic* | One column for each parameter and statistic ticked, e.g. `CD4 · Median` |

Click a header to sort by it; numbers sort as numbers. Counts are written with thousands separators, percentages and CVs with two decimals, other numbers with three significant figures (one decimal at 100 and above). The exported and copied values keep six significant figures.

The caption above the table says how many rows, samples and populations it holds, and whether the samples were batch-normalised as they were read (see [[normalisation]]). When the gates, the samples or their columns change after Compute, a warning appears beside it: **Computed before the last change — Compute again**. The table still shows what it was computed from until you do.

## Taking the table elsewhere

- **Export CSV…** — the whole table as `statistics.csv` (or a name you choose). The file is UTF-8 with a byte-order mark, so Excel on Windows reads `CD4 · Median` correctly.
- **Copy** — the whole table, tab-separated with its header, to paste into Excel or Prism.
- **Prism…** — one column of the table as a GraphPad Prism file: see [[tables-and-prism]].
- **Compare groups…** — test a column between groups of samples: below.

The same table can be made without the window: `milc stats` on the [[command-line]].

## Compare groups

**Compare groups…** opens a dialog that asks what to compare and how.

| Field | What it does |
|---|---|
| **Compare** | The column of the table to test: Count, any percentage, or any statistic column. Starts on % of parent |
| **Between** | The sample columns that make the groups. Tick one (Group) or several: then each combination is a group, written `Treated · 7` |
| **Within each** | Optional. A separate comparison for each value of this column — the groups within each day, say |
| **Paired by** | Optional. The column that says which animal or donor each sample is of, where the same ones were measured in every group |
| **Test** | **Welch t test / ANOVA**, **Mann-Whitney / Kruskal-Wallis**, or **Welch t test / Welch ANOVA** |
| **Adjust p by** | **Benjamini-Hochberg**, **Holm**, **Bonferroni** or **None** |
| **With more than two groups, every pair as well** | Ticked by default: besides the test of all groups together, one test for each pair |

One test is made for each population (and each value of *Within each*), and the p values are adjusted over all the tests together, pairs included. The unit is the sample, never the cell: five animals are five observations however many cells each gave. The tests, the pairing and the notes are explained in [[cluster-statistics]], which uses the same machinery for clusters.

The result opens in its own window, a table titled after the value compared:

| Column | Meaning |
|---|---|
| *Within each* | The value of that column, when one was chosen |
| Population | What was tested |
| Groups | The groups in this test, `Control vs Treated` |
| Mean ± SD (n) | Each group's mean, SD and number of samples |
| Test | The test that was made: Welch t, ANOVA, Welch ANOVA, Mann-Whitney, Kruskal-Wallis, Paired t, Wilcoxon signed-rank, RM ANOVA (Greenhouse-Geisser) or Friedman |
| Paired by | The pairing column, in a paired design |
| Statistic | t, F, U, H, W or χ² of the test |
| p | The p value of that test alone |
| p adjusted (*method*) | After adjustment for the number of tests |
| (unnamed) | Stars from the adjusted p: `****` < 0.0001, `***` < 0.001, `**` < 0.01, `*` < 0.05, `ns` otherwise |
| Note | Why a test could not be made, or what to know about it |

A blank p is a test that could not be done; the Note says why. p values below 0.0001 are shown as `< 0.0001`. **Export CSV…** (as `group_comparison.csv`), **Copy** and **Close** are under the table.

## Limits

- Statistics follow FlowJo's definitions but were not set beside FlowJo's numbers on the same files (see [[measured-facts]]).
- The tests are one value per sample, by group, paired or not. There are no mixed models.
- A population that some gated sample lacks is not offered at all; give that sample the gate first.
