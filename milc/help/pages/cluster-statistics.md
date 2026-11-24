---
title: Comparing clusters between groups
---

Is a cluster a larger share of the samples in one group than in another, more than chance would make it? [[clusters-tab|The Clusters tab]] answers this three ways:

- the **Test** in the tab's second row, on each sample's percentage in each cluster — described here;
- **Differential abundance…**, on the counts, by a negative-binomial model: see [[abundance]];
- **Differential expression…**, on what the cells of each cluster express: see [[differential-expression]].

The same tests, on any column of the statistics table, are under **Compare groups…** in [[statistics-tab|the Statistics tab]].

## The unit is the sample

Every test here takes one value per sample — the share of the sample's events in the cluster — and compares those values across groups. Five animals are five observations however many cells each gave. Testing cells against cells would find every difference "significant" and say nothing about animals. The n of a group is its number of samples.

The share is of the sample's **concatenated** events: what [[concatenation]] drew from it, after the sampling.

## The tests

Choose under **Test**:

| Choice | Two groups | Three or more |
|---|---|---|
| **Welch t test / ANOVA** | Welch's t test, which does not assume equal variances | One-way ANOVA |
| **Mann-Whitney / Kruskal-Wallis** | Mann-Whitney U test (two-sided) | Kruskal-Wallis |
| **Welch t test / Welch ANOVA** | Welch's t test | Welch's ANOVA (Welch, 1951): each group weighs n / variance, for groups whose variances differ |

The first is parametric and the most powerful when the values are roughly normal within each group; the rank tests make no assumption about the shape but cannot give a small p with few samples. With three samples against three, Mann-Whitney never goes below p = 0.1, whatever the data; the Note says so ("smallest p possible with these n").

With more than two groups, **Results…** (and **Compare groups…** in Statistics, when its box is ticked) also tests every pair of groups, beside the test of all of them together. The chart shows only the test of all groups.

### Paired designs

Where the same animals or donors were measured in every group — before and after, or several conditions of the same donor — choose under **Paired by** the column that says which one each sample is of. The groups are then conditions and the test is of each subject's change:

| Choice | Two groups | Three or more |
|---|---|---|
| Welch t (either kind) | Paired t test | Repeated-measures ANOVA, Greenhouse-Geisser corrected |
| Mann-Whitney / Kruskal-Wallis | Wilcoxon signed-rank test | Friedman test |

Only subjects with a value in every group being compared are used; the others are left out and counted in the Note. If a subject appears twice in one group, or no subject is in every group, there is no test, and a warning beside the controls says why: "No test: the samples cannot be paired by Donor — …".

The repeated-measures ANOVA removes what differs between subjects before comparing the groups. Its F test assumes that the differences between any two groups scatter alike (sphericity) and gives too small a p where they do not; the Greenhouse-Geisser ε (from 1 when that holds down to 1/(k − 1)) multiplies both degrees of freedom. The Note gives ε and the p without the correction.

## Adjusting for many tests

Every cluster is a test, and with twenty clusters one will be "significant" at 0.05 by chance alone. The p values are therefore adjusted over all the tests made together — every cluster, every panel, every pair:

| Choice | What it controls |
|---|---|
| **Benjamini-Hochberg** (default) | The false discovery rate: among the clusters called significant, the expected share that are not. The adjusted p is often called the FDR or q value |
| **Holm** | The family-wise error: the chance of any false call at all. Stricter than Benjamini-Hochberg, never stricter than Bonferroni |
| **Bonferroni** | The same, more simply: p × the number of tests, at most 1 |
| **None** | No adjustment |

Tests that could not be done do not count as tests.

## On the chart

With a test chosen, each cluster's bars in the **Frequency** panel carry the adjusted p of the test of all groups — `p = 0.003`, `p < 0.0001` — where it is below 0.05, and `ns` where it is not. Every number is in **Results…**.

## Results…

**Results…** opens the table "Clusters between groups": a row for every comparison, pair by pair where there are more than two groups.

| Column | Meaning |
|---|---|
| *within each* | The panel, when the tab has one |
| Cluster | The cluster's name |
| Groups | The groups compared, `Control vs Treated` |
| Mean ± SD (n) | Each group's mean share (%), SD and number of samples |
| Test | The test made |
| Paired by | The pairing column, in a paired design |
| Statistic | t, F, U, H, W or χ² |
| p | Of that test alone |
| p adjusted (*method*) | After the correction chosen |
| (unnamed) | `****` < 0.0001, `***` < 0.001, `**` < 0.01, `*` < 0.05, `ns` |
| Note | Why a test could not be made, or what to know about it |

A blank p is a test that could not be done; the Note says why: one group only; a group with fewer than 2 samples; every sample with the same value; no variation inside the groups. A p value is never invented. **Export CSV…** writes `cluster_comparison.csv`; **Copy** puts the table on the clipboard.

## Choosing between the three analyses

| Question | Use |
|---|---|
| Does the share of a cluster differ, with clear effects and enough samples? | **Test** in the tab: quick, shown on the chart |
| Does the share differ, with rare clusters or samples of very different size? | Differential abundance, [[abundance]]: it models the counts and found about twice as many true changes in simulation |
| Do the cells of a cluster express more of a marker? | [[differential-expression]] |

None of them is a mixed model: one value per sample, compared by group, paired or not, optionally adjusted for one more factor in the two differential analyses.
