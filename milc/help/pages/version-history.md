---
title: Version history
---

## 1.0.0 — 24 November 2026

The first release. It covers the route it was written to replace — gate in FlowJo, DownSample,
concatenate, run the UMAP, FlowSOM and ClusterExplorer plugins — and a good deal past it. What
1.0 contains, by area:

### Files and samples

- FCS 2.0, 3.0 and 3.1 read and written with numpy alone; FlowJo 10 workspaces read — gates,
  compensation and FlowJo's counts to compare — including quadrants, boolean combinations,
  ellipses and curly quadrants ([[opening-data]]).
- The Samples table is yours to type in: free columns, values typed into many rows at once,
  columns read from file names with a pattern, a group made of one column or several
  ([[samples-tab]]).
- One order of the groups for every chart, legend, panel, table and Prism file, kept with the
  project ([[sample-order]]).
- Everything done to a folder kept beside the files in `MIL-C.json` ([[projects-and-saving]]).
- Re-acquired tubes flagged before they are counted twice.

### Concatenation, normalisation, quality

- Balanced concatenation with the plan on the table before anything is drawn: the same number
  from every sample, the same total from every group, or a ceiling; a minimum; a seeded draw that
  does not depend on the other samples ([[concatenation]]).
- Batch normalisation in the manner of CytoNorm, on each marker whole or within FlowSOM clusters,
  onto the mean of the batches or onto one batch, optionally applied to the samples themselves;
  a before-and-after view on held-out events ([[normalisation]]).
- Quality control of each acquisition over time — flow rate, signal, lasting drifts and short
  shifts — as a Clean population at the top of every tree, with its limits in the dialog
  ([[quality-control]]).

### Gates

- Rectangle, polygon, lasso, ellipse, range and quadrant gates, quadrants whose arms move one by
  one, boolean combinations, curly quadrants drawn and saved; a gate for every sample or for one
  ([[drawing-gates]]).
- Automatic gates — a valley, a density contour, a robust ellipse, a singlet band — fitted again
  by themselves when the population they are in changes ([[auto-gates]]).
- Biexponential, Logicle, arcsinh, logarithmic and linear scales per parameter
  ([[scales-and-transforms]]).
- A compensation matrix per sample or shared, edited in a table ([[compensation]]); a matrix
  from single-stain controls by AutoSpill's method or the classic one ([[autocompensation]]);
  spectral unmixing by OLS or WLS ([[spectral-unmixing]]).
- Overlays, violins and boxes, a matrix of every pair of parameters, backgating
  ([[gates-tab]]).
- Models of one population: cell cycle ([[cell-cycle]]), proliferation ([[proliferation]]) and
  kinetics ([[kinetics]]).
- Gates written out as a FlowJo workspace or as Gating-ML 2.0
  ([[workspaces-and-gating-ml]]).

### Statistics

- A statistics table of every population of every sample: counts, percentages, medians, means,
  CVs and robust measures by FlowJo's definitions ([[statistics-tab]]).
- Tests between groups, paired and repeated-measures designs, Welch's ANOVA, with correction for
  many tests.
- Tables sent to GraphPad Prism as `.pzfx` files ([[tables-and-prism]]).
- Differential abundance of clusters in the manner of diffcyt's DA-edgeR ([[abundance]]) and
  differential expression within clusters in the manner of DS-limma
  ([[differential-expression]]).

### Explore and Clusters

- UMAP, tSNE, FlowSOM, PhenoGraph-style clustering with Leiden or Louvain and a minimum cluster
  size, and PCA, run in the background ([[embeddings]], [[clustering-methods]]).
- Maps coloured and split by any column, or by two ([[explore-tab]]).
- The map, heatmap, tree and frequencies of the clusters, linked; clusters named, annotated,
  merged and turned into populations; tests per cluster ([[clusters-tab]],
  [[cluster-statistics]]).

### Figures and reports

- Every plot as a figure in millimetres, with its type, colours and labels, as SVG, PDF, PNG or
  TIFF; figures of the Gates tab carry their gates ([[figures]]).
- Plates of several panels, saved with the project, and a PDF report of how the analysis was
  done ([[the-plate]]).

### Platform

- Installers for macOS (Apple Silicon), Windows (`.msi` and a portable `.zip`) and Linux (`.deb`
  and `.tar.gz`), built and tested by continuous integration on all three systems, with
  checksums ([[installation]]).
- A command line that does what the window does, and a self-test ([[command-line]]).
- A light and a dark appearance; a layout that fits a laptop's screen and Windows' wider frames.
- This manual, in English and Portuguese, opened on the page for where you are, and printable as
  one PDF ([[welcome]]).

### Checked against

FlowJo on real acquisitions, and the R and Python packages each method follows — edgeR,
PhenoGraph, PeacoQC, flowAI, CytoNorm, AutoSpill, CytoML, SpectroFlo's exported values and
AutoSpectral. The results, with where each was measured, are in [[measured-facts]].
