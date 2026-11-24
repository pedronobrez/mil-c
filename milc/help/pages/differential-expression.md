---
title: Differential expression
---

Differential expression asks whether a cluster's cells carry more or less of a marker in one group than in another. It is the companion of [[abundance|differential abundance]], which asks whether there are more or fewer of the cells. The method is diffcyt's DS-limma (Weber et al., *Commun Biol* 2019): each sample's median of each marker in each cluster, compared between groups by a linear model and limma's moderated t (Smyth, 2004).

Open it with **Differential expression…** in [[clusters-tab|the Clusters tab]]. It uses the clusters as the tab shows them and the grouping chosen there under **Compare** and **within each**.

## The controls

The first row:

| Control | What it does |
|---|---|
| **Compare** *column* | The group the difference is of. A positive difference is higher in this group |
| **against** | The reference group |
| **Adjust for** | Optional: another column of the samples — a batch, a day. Its effect is estimated beside the group's and taken out. The column that says which animal or donor each sample is of makes the comparison a paired one |
| **Markers** | A menu of tick boxes, labelled "3 of 8" and so on. Every stained marker is offered; those the clustering did not use are ticked to begin with (all of them if it used every one) |
| **Scale** | **Biexponential** or **arcsinh (cofactor 150)**: the scale the medians are taken on, and the difference is in. Starts on the scale the clustering ran on |

The second row:

| Control | Default | What it does |
|---|---|---|
| **Cells in a sample, at least** | 10 | A sample with fewer cells of a cluster than this gives no median for it: a median of three cells is noise |
| **Samples in a group, at least** | 2 | A cluster with fewer such samples than this in either group is not tested, and the note says so |
| **Weigh by cells** | off | Each median weighs as the number of cells it is of, as in diffcyt |
| **Correction** | Benjamini-Hochberg | Over every cluster and marker tested: Benjamini-Hochberg, Holm, Bonferroni or None |

Every change recomputes at once. With a **within each** column in the Clusters tab, the analysis — prior and all — is made separately within each of its values, and p is adjusted over every panel's tests together.

## The method

1. For each sample, cluster and marker: the median of the marker among the sample's cells of the cluster, on the chosen scale.
2. For each cluster and marker, a linear model `median ~ intercept + group (+ the factor adjusted for)` over the samples of the two groups. The group's coefficient is the **difference** in median between the groups, in the units of the scale.
3. Each row's residual variance is moved towards what all the rows suggest (limma's empirical Bayes: `lmFit` + `eBayes`, no trend, not robust). This is worth extra degrees of freedom where there are few samples and many rows. The moderated t then gives the p.
4. The p values are adjusted over every cluster and marker tested together.

The factor adjusted for is a fixed effect: a column for each of its levels. With the animal as the factor and every animal in both groups, the unmoderated test is the paired t test.

### Reading the difference

- **Biexponential**: the medians are on the plots' scale as numbers from 0 to 1023; where the scale is logarithmic, 100 is about half a decade.
- **arcsinh**: asinh(x / 150); where the marker is bright, 0.69 is a doubling.

## The result

The table on the left, one row per cluster and marker:

| Column | Meaning |
|---|---|
| *within each* | The panel, when there is one |
| Cluster, Marker | What was tested |
| *Group A* (median), *Group B* (median) | Each group's mean of its samples' medians |
| Difference | The group's coefficient: first less second, adjusted for the rest |
| p | The moderated t test's p |
| p adjusted (*method*) | After the correction chosen |
| (unnamed) | Stars from the adjusted p |
| n | Samples with a median in each group, `5 vs 4` |
| Moderated t | The statistic |
| df | Residual degrees of freedom plus those the other rows were worth |
| Note | Why a row was not tested |

On the right, two tabs:

- **Heatmap** — clusters down, markers across, each cell the difference: red higher in the first group, blue in the second, white for none, grey not tested. Stars where the adjusted p is below 0.05 (`*` < 0.05, `**` < 0.01, `***` < 0.001, `****` < 0.0001). The line under it says what the deepest colour is.
- **Volcano** — every cluster and marker a point, the difference across and −log10 of the adjusted p up, with one guide line at no difference; significant points are coloured and named `cluster · marker`.

**Figure…**, at the corner of the tabs, makes a figure of the one shown (see [[figures]]). A paragraph under them says what was done and what it cannot show. **Export CSV…** (as `differential_expression.csv`), **Copy** and **Close** are at the bottom.

## What it cannot say

- A median is of the cluster as it was drawn. A marker the clusters were built on is held within a narrow range inside each cluster by construction: what a group does to it shows as cells moving between clusters (abundance), and only its remainder here. Prefer markers the clustering did not use.
- A minority of the cluster switching a marker on — a bimodal shift — moves the median little or not at all.
- The cells are those that were concatenated, not all that were acquired.

## How well it is calibrated

The moderated t agrees with limma 3.68 to about seven significant figures (see [[measured-facts]]). In simulations with five or three samples a group, unweighted, p < 0.05 came about as often as it should under no difference (5.1% and 5.4%); with weights about 6.5%, too often, which is why **Weigh by cells** is off. Two samples a group found next to nothing. diffcyt itself was not run, and limma's newer estimator for rows with missing samples is not implemented (this follows `eBayes(legacy = TRUE)`).
