---
title: The Explore tab
---

The Explore tab runs the algorithms on the concatenated events and shows the result: any two axes, coloured by group, sample, cluster or the expression of a marker, in one panel or one per group. Until the samples have been concatenated in [[samples-tab|the Samples tab]] (see [[concatenation]]) it says so and its buttons are off.

The tab has three parts: the settings on the left, the plot in the middle, and how to look at it on the right.

## Left: what to compute

### Markers

The list of stained fluorescence channels — scatter and time are not markers. The ones with a stain name are ticked to begin with. These are the markers the algorithms see: untick what was used to gate (viability, lineage), which would otherwise dominate the map. At least two must be ticked.

### Scale and seed

| Control | Choices | Meaning |
|---|---|---|
| **Scale** | **Biexponential** (default), **Arcsinh, cofactor 150**, **Logicle, W 0.5, M 4.5** | How the values are transformed before the algorithms see them. Biexponential is each marker on the scale its axis has in Gates — FlowJo's biexponential unless another was chosen — as numbers from 0 to 1023, which is what FlowJo feeds its plugins. Logicle is the standard Logicle on 0 to 1023. See [[scales-and-transforms]] |
| **Seed** | 0 to 999,999, default 0 | The random seed. The same seed, markers and settings give the same result |

Nothing is scaled further: the markers arrive on one scale already.

### The algorithms

Each has its own settings and **Run** button. They run in the background, and several can run at once; the button says "Running…" meanwhile.

| Algorithm | Settings (default) | Gives |
|---|---|---|
| UMAP | **Neighbours** (15), **Min. distance** (0.5) | Two axes, UMAP1 and UMAP2 |
| tSNE | **Perplexity** (30) | Two axes, tSNE1 and tSNE2 |
| FlowSOM | **Clusters** (5), **Grid** (10: a 10 × 10 map) | A cluster for every event |
| PhenoGraph | **Neighbours** (30), **Resolution** (1.0), **Method** (Leiden or Louvain), **Min. cluster size** (10) | A cluster for every event, without being told how many |
| PCA | none | Two axes, PCA1 and PCA2, each with its share of the variance |

**Run all three** runs UMAP, tSNE and FlowSOM together. What each algorithm does, and what its settings change, is in [[embeddings]] and [[clustering-methods]]. Running one again replaces its previous result; the seed, markers, scale and settings are kept with the result and saved with the concatenated file.

Under the button, a few lines describe what has been computed: the number of events and samples; for UMAP and tSNE, the share of each event's 15 nearest neighbours the map keeps; for PhenoGraph, the number of clusters and the modularity; for PCA, each component's share of the variance and the three markers that weigh most in it, e.g. `PCA1, 41%: +0.72 CD4, −0.55 CD8, +0.30 CD25`.

## Middle: the plot

Every concatenated event as a dot. A new embedding is shown as soon as it is computed; a new clustering colours the plot.

## Right: how to look at it

| Control | What it does |
|---|---|
| **X**, **Y** | The axes: UMAP1/2, tSNE1/2, PCA1/2 (with the share of variance in the label), or any parameter. A marker is drawn on its scale from Gates and ticked in decades |
| **Colour by** | **One colour**, the grouping chosen in Samples, any column of the samples that differs between them, **Sample**, **FlowSOM cluster**, **PhenoGraph cluster**, or — after a separator — the expression of any marker |
| **Split by** | **One panel**, or one panel per value of any of the same columns, samples or clusters. In each panel the other events are drawn faint behind, so where a group sits is read against where everyone is |
| **and by** | A second subdivision: one panel for every combination, `Treated · 7`. Rows are the values of the first, columns of the second |
| **Ramp** | The colour ramp for expression |
| **Point size** | 1 to 14, default 4 |
| **Opacity** | 5% to 100%, default 80% |

The **Legend** lists each colour with its number of events; for expression it gives the range drawn. The colours of expression span the 1st to the 99th percentile, so a few bright events do not wash out the rest.

**Figure…** turns the plot as shown into a figure — size, type, colours, legend — in SVG, PDF, PNG or TIFF: see [[figures]].

The order in which groups appear in the legend and panels is set in [[sample-order]].

## Where the results go

The clusters appear in [[clusters-tab|the Clusters tab]], with their heatmap, frequencies and tests. The embeddings are offered there as the map. The [[command-line]] does the same with `milc analyze`.
