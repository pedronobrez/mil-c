---
title: Glossary
---

The terms you meet in MIL-C, in alphabetical order, each with the page where it is used.

## A – C

- **Adjusted p** — a p value corrected for the number of tests made at once, so that testing twenty clusters does not produce one "significant" result by chance alone. See [[cluster-statistics]].
- **ANOVA** — analysis of variance: a test of whether the means of three or more groups differ. Welch's ANOVA allows unequal variances; repeated-measures ANOVA is for the same subjects in every group. See [[cluster-statistics]].
- **Arcsinh** — the inverse hyperbolic sine, asinh(x / cofactor): linear near zero, logarithmic for large values, defined for negative values. One of the display and analysis scales. See [[scales-and-transforms]].
- **Autocompensation** — computing the spillover matrix from single-stain controls rather than typing it. See [[autocompensation]].
- **Autofluorescence** — the light cells give off unstained; in spectral unmixing it can be one more spectrum. See [[spectral-unmixing]].
- **Backgating** — showing where a population's events fall on each gate above it. See [[gates-tab]].
- **Batch normalisation** — bringing samples acquired on different days (batches) onto one another, learnt from reference samples that should look alike in every batch. See [[normalisation]].
- **Benjamini-Hochberg** — the adjustment that controls the false discovery rate; the default in MIL-C. See [[cluster-statistics]].
- **Biexponential** — FlowJo's display scale, linear around zero and logarithmic above, with a width and a number of decades. MIL-C's default scale for markers, and the one the algorithms are fed, as numbers from 0 to 1023. See [[scales-and-transforms]].
- **Bonferroni, Holm** — adjustments that control the chance of any false call at all; stricter than Benjamini-Hochberg. See [[cluster-statistics]].
- **Boolean gate** — a population made of others: the events in all of them, in any of them, or in none. See [[drawing-gates]].
- **Cell cycle** — G1, S and G2 fitted to a DNA-content histogram. See [[cell-cycle]].
- **Cluster** — a group of events alike in the markers, found by an algorithm (FlowSOM, PhenoGraph) or drawn by hand. See [[clustering-methods]] and [[clusters-tab]].
- **Cofactor** — the divisor in arcsinh(x / cofactor), which sets where the scale turns from linear to logarithmic: 150 is usual for fluorescence, 5 for mass cytometry. See [[scales-and-transforms]].
- **Compensation** — removing from each fluorescence channel the light that spilt into it from other fluorochromes, by the spillover matrix. Compensated channels are named `Comp-…`. See [[compensation]].
- **Concatenation** — putting the events of many samples into one set, each event remembering its sample, for the embeddings and clusterings. See [[concatenation]].
- **Consensus clustering** — FlowSOM's way of merging its map's nodes into metaclusters: many resampled hierarchical clusterings, combined. See [[clustering-methods]].
- **Count** — the number of events in a population. See [[statistics-tab]].
- **CV** — coefficient of variation: 100 × SD / mean. The robust CV is 100 × robust SD / median. See [[statistics-tab]].

## D – H

- **Differential abundance** — whether a cluster is a larger or smaller share of the sample in one group than another, tested on counts. See [[abundance]].
- **Differential expression** — whether a cluster's cells express more or less of a marker in one group than another, tested on each sample's median. See [[differential-expression]].
- **Dispersion** — in a negative-binomial model, how much more counts vary between samples than chance sampling of cells would make them. See [[abundance]].
- **Downsampling** — drawing fewer events than a sample has. In MIL-C it is done by the sampling plan, which balances the samples. See [[concatenation]].
- **Embedding** — a two-dimensional map of the events, such as UMAP or tSNE. See [[embeddings]].
- **Event** — one particle that passed the laser: usually a cell. One row of an FCS file. See [[opening-data]].
- **FCS** — Flow Cytometry Standard, the file format cytometers write. MIL-C reads versions 2.0, 3.0 and 3.1; your files are only read, and exports are new files. See [[opening-data]].
- **FDR** — false discovery rate: among the results called significant, the expected share that are not. The Benjamini-Hochberg adjusted p is often called the FDR or q value. See [[cluster-statistics]].
- **FlowSOM** — clustering by a self-organising map whose nodes are merged into a chosen number of metaclusters. See [[clustering-methods]].
- **FMO** — fluorescence minus one: a control stained with every fluorochrome but one, to see where that one's negatives end. MIL-C has no FMO-specific tool; overlay the control on a sample in Gates to place a gate. See [[gates-tab]].
- **Fold change** — the ratio between groups; in differential abundance it is given as log2, so 1 is twice the share and −1 half. See [[abundance]].
- **Gate** — a region on one or two parameters that selects events; gates nested in one another form the gating tree. See [[drawing-gates]] and [[auto-gates]].
- **Gating-ML** — the standard file format for gates, read by FlowKit, flowCore and Cytobank. See [[workspaces-and-gating-ml]].
- **Geometric mean** — exp of the mean of the logarithms; events at or below zero are left out. See [[statistics-tab]].
- **Greenhouse-Geisser** — the correction of a repeated-measures ANOVA for groups whose differences do not scatter alike. See [[cluster-statistics]].
- **Group** — the samples that share a value of a column (or of several) in the Samples table. See [[samples-tab]] and [[sample-order]].

## J – M

- **Jaccard** — the share of two sets that they have in common. PhenoGraph weighs each join by the Jaccard index of the two events' neighbourhoods. See [[clustering-methods]].
- **Kinetics** — a parameter followed against time, as in calcium flux. See [[kinetics]].
- **Kruskal-Wallis, Mann-Whitney** — rank tests, for three or more groups and for two, that make no assumption about the shape of the values. See [[cluster-statistics]].
- **Leiden, Louvain** — algorithms that cut a graph into communities of high modularity; PhenoGraph's last step. Leiden is the default. See [[clustering-methods]].
- **Logicle** — a biexponential-type scale with standard parameters (W, M, A, T). See [[scales-and-transforms]].
- **MAD** — median absolute deviation: the median of the distances from the median. A statistic in the table, and the unit of the limits in quality control. See [[statistics-tab]] and [[quality-control]].
- **Median** — the middle value of the events: the usual figure for a marker's level, robust to outliers. See [[statistics-tab]].
- **Metacluster** — a FlowSOM cluster: a group of the map's nodes. In the Clusters tab, simply a cluster. See [[clustering-methods]].
- **MFI** — mean (or median) fluorescence intensity. MIL-C reports the median, mean or geometric mean by name rather than "MFI". See [[statistics-tab]].
- **Moderated t** — limma's t statistic, whose variance is moved towards what all the rows suggest; worth extra degrees of freedom with few samples. See [[differential-expression]].
- **Modularity** — how much more weight a partition of a graph keeps inside its communities than chance would. PhenoGraph maximises it. See [[clustering-methods]].

## N – R

- **Negative binomial** — the count distribution used for differential abundance: Poisson-like, with extra variation between samples. See [[abundance]].
- **Neighbour preservation** — the share of each event's 15 nearest neighbours that are still its neighbours on a map: how MIL-C judges an embedding. See [[embeddings]].
- **Paired design** — the same animals or donors measured in every group; the tests then look at each one's change. See [[cluster-statistics]].
- **PCA** — principal component analysis: the directions along which the events vary most, each a readable mix of markers. See [[embeddings]].
- **Perplexity** — tSNE's setting, roughly the number of neighbours each event weighs. See [[embeddings]].
- **PhenoGraph** — clustering without choosing the number: a nearest-neighbour graph, Jaccard weights and communities. See [[clustering-methods]].
- **Plate** — several figure panels on one page. See [[the-plate]].
- **Population** — the events a gate selects, with its count and percentages. See [[gates-tab]].
- **Prism, .pzfx** — GraphPad Prism and its file format, which MIL-C writes. See [[tables-and-prism]].
- **Proliferation** — generations of a dye diluted by each division, fitted to a histogram. See [[proliferation]].
- **Quality control** — finding stretches of an acquisition where the flow or the signal went wrong, and leaving them out. See [[quality-control]].
- **Resolution** — PhenoGraph's setting: above 1, more and smaller clusters; below, fewer and larger. See [[clustering-methods]].
- **Robust SD** — half the distance between the 15.87th and 84.13th percentiles: the SD of a normal distribution with that spread, deaf to outliers. See [[statistics-tab]].

## S – W

- **Sampling plan** — how many events each sample gives to the concatenation, and why. See [[concatenation]].
- **Seed** — the number that fixes an algorithm's random choices, so the same seed gives the same result. See [[explore-tab]].
- **Self-organising map** — a grid of nodes trained to cover the events, the first step of FlowSOM. See [[clustering-methods]].
- **Single-stain control** — a sample stained with one fluorochrome only, from which spillover or a spectrum is measured. See [[autocompensation]] and [[spectral-unmixing]].
- **Spillover** — the share of a fluorochrome's light read in another fluorochrome's detector. The spillover matrix holds it for every pair. See [[compensation]].
- **Spectral unmixing** — on a spectral cytometer, recovering each fluorochrome's signal from many detectors and their reference spectra. See [[spectral-unmixing]].
- **SD** — standard deviation, over n − 1. See [[statistics-tab]].
- **tSNE** — t-distributed stochastic neighbour embedding, a map that keeps local neighbourhoods. See [[embeddings]].
- **UMAP** — uniform manifold approximation and projection, a map built from a nearest-neighbour graph. See [[embeddings]].
- **Volcano plot** — fold change or difference across, −log10 of the adjusted p up; the clusters that matter end in the top corners. See [[abundance]] and [[differential-expression]].
- **Welch's t test** — a t test that does not assume the two groups have equal variances; MIL-C's default for two groups. See [[cluster-statistics]].
- **Workspace** — a FlowJo `.wsp` file: samples, compensation, scales and gates. MIL-C reads and writes them. See [[workspaces-and-gating-ml]].
