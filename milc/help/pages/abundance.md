---
title: Differential abundance
---

Differential abundance asks whether a cluster is a larger or smaller share of the sample in one group than in another. It is asked of the **counts** — how many of each sample's cells fall in the cluster — not of the percentages, by a negative-binomial model in the manner of diffcyt's DA-edgeR (Weber et al., *Commun Biol* 2019).

Open it with **Differential abundance…** in [[clusters-tab|the Clusters tab]]. It works on the clusters as the tab shows them — merged, renamed, with any drawn by hand — and on the grouping chosen there under **Compare** and **within each**.

## Why counts and not percentages

A percentage forgets what it was made of: 0.2% from 109 cells and 0.2% from 50,000 are the same number to a t test, and the percentages of a rare cluster are far from normal. A count model knows that a sample with few cells gives a less certain share. In simulations with five samples a group it found about twice as many true two-fold changes as Welch's t test on the percentages, which the Clusters tab offers too (see [[cluster-statistics]]).

## The controls

| Control | What it does |
|---|---|
| **Compare** *column* | The group in the numerator. A positive fold change means more in this group |
| **against** | The reference group |
| **Adjust for** | Optional: another column of the samples — a batch, a day, the animals' sex. Its effect is estimated beside the group's and taken out, as an additive term of the model. Only columns whose value differs between samples, and that are not the grouping or the panel column, are offered |
| **Test** | **Quasi-likelihood F test** (the default) or **Likelihood-ratio test** |

Every change recomputes at once. Only the samples of the two groups compared are used, for the dispersion too. Where the Clusters tab has a **within each** column, the whole analysis — model, dispersion and all — is made separately within each of its values, and the p values are then adjusted over every panel's tests together.

## The method

For each cluster, the count in each sample is modelled as

`count ~ NB(mean = total × exp(X b), variance = mean + dispersion × mean²)`

where *total* is the sample's number of cells in all clusters (the offset) and *X* holds an intercept, the group and, if chosen, the factor adjusted for. The group's coefficient is the log fold change, and it is tested.

- **Fit** — iteratively reweighted least squares with step halving.
- **Dispersion** — one common dispersion for all clusters, maximising the sum of their Cox-Reid adjusted profile likelihoods (edgeR's `estimateGLMCommonDisp`), searched on a grid from 2⁻²⁰ to 2³.
- **Quasi-likelihood F test** — the model is fitted at the common dispersion; each cluster's quasi-dispersion (deviance over residual degrees of freedom) is squeezed towards the others' by limma's empirical Bayes, and F = likelihood ratio / squeezed quasi-dispersion. This is edgeR's `glmQLFit`/`glmQLFTest` as before version 4. It allows for the dispersion being an estimate from few samples, and is the better calibrated.
- **Likelihood-ratio test** — what `diffcyt::testDA_edgeR` calls: a dispersion for each cluster, shrunk towards the common one with 10 prior degrees of freedom, then a chi-square test. It finds more, and is a little too generous with few samples.
- **Fold change** — from a second fit with a prior count of 0.125 per sample (scaled by its total), as `glmFit(prior.count = 0.125)` does. A cluster absent from a whole group then gets a large but finite fold change, with a note. The test does not use the prior count.
- **Adjustment** — Benjamini-Hochberg over every cluster tested (and every panel).
- No TMM normalisation factors and no filter on rare clusters. A cluster with no cells in either group is not tested.

## The result

The table on the left has one row per cluster (and per panel):

| Column | Meaning |
|---|---|
| *within each* | The panel, when the tab has one |
| Cluster | The cluster's name |
| *Group A* (mean %) | The mean, over the group's samples, of the cluster's percent of the sample |
| *Group B* (mean %) | The same for the reference |
| log2 fold change | The group's coefficient in log2: 1 is twice as large a share, −1 half |
| p | The test's p value |
| p adjusted (Benjamini-Hochberg) | The false discovery rate |
| (unnamed) | Stars from the adjusted p, as in [[statistics-tab]] |
| n | Samples in each group, `5 vs 5` |
| Test | The test used |
| Statistic | F or the likelihood ratio |
| Dispersion | Of the model the test was made with |
| Note | Why no test was made, or that a group had no cells |

The **Volcano** on the right draws every tested cluster: log2 fold change across, −log10 of the adjusted p up, with a dashed line at FDR 0.05 and guide lines at ±1 (two-fold). Clusters below FDR 0.05 are coloured — "More in *A*" or "More in *B*" — and named; the rest are grey. With panels, one volcano per panel on the same axes. **Figure…** makes it a figure (see [[figures]]).

Under them, a paragraph says what was done and what it does not show; it is part of the result. **Export CSV…** (as `differential_abundance.csv`), **Copy** and **Close** are at the bottom.

## Limits

- The cells counted are those that were concatenated, after the sampling, not all that were acquired (see [[concatenation]]).
- What is tested is a cluster's **share** of the sample: one cluster that grows makes every other a little smaller.
- In simulations, five samples a group gave p < 0.05 about as often as it should under no difference with the quasi-likelihood test; with two or three a group, or with the likelihood-ratio test, somewhat too often, and very small p values more so. Two samples a group is not enough for either test.
- Set beside edgeR 4.10 called the same way, it gives the same numbers; against edgeR's and diffcyt's present defaults, 10 to 20% of discoveries differ (see [[measured-facts]]).

For whether the cells of a cluster express more or less of a marker, see [[differential-expression]].
