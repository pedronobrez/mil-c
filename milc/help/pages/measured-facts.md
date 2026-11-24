---
title: Measured against other software
---

Most methods in MIL-C were written from the papers, not copied from the programs that made them known. Each was then set beside the program it takes after, on the same input, and what is written here is what was measured. Where a number was not measured, it does not appear; where something was not checked, this page says so. The reasoning behind this is in [[design-principles]].

The comparisons with R and Python packages were run on 23 November 2026 with R 4.6.1, Bioconductor 3.23 and Python 3.14. The scripts are in `tests/reference/` of the source code; they need R and are not part of the test suite. "A real experiment" below means acquisitions that are not distributed with MIL-C; the public datasets are named.

## Against FlowJo

On a real experiment of 12 gated samples with five markers, and two concatenated files of about 2,000 events each. The test suite asserts these with a margin.

| What | MIL-C against FlowJo |
|---|---|
| Gate counts, 72 gates | Within 0.77% on every gate; within 3 events on the final population. FlowJo evaluates polygons on a 256-bin grid by a rule that was not reproduced; MIL-C uses the exact geometry |
| Biexponential transform | Medians per cluster within 0.69 of a channel (of 1,023) |
| FlowSOM clusters (adjusted Rand index) | 0.86 and 0.81 against FlowJo's; two runs of MIL-C with different seeds agree at 0.81 and 0.79. FlowJo's run is as close to MIL-C's as two of MIL-C's are to each other |
| UMAP: share of each event's 15 nearest neighbours kept | 47.0% and 47.9%; FlowJo 47.2% and 47.7% |
| tSNE: the same measure | 58.8%; FlowJo 58.0% |

UMAP is `umap-learn`, not the `uwot` the FlowJo plugin calls: the pictures differ, the neighbourhoods they keep do not. See [[embeddings]] and [[clustering-methods]].

**Not set beside FlowJo:** the [[statistics-tab|statistics]] (they follow FlowJo's definitions); curly quadrants (FlowJo writes no curve); proliferation and cell cycle (nor beside ModFit). Ellipses read from a workspace were checked once, on a 100-event published file with linear axes.

## Differential abundance against edgeR and diffcyt

[[abundance|Differential abundance]] beside edgeR 4.10.7 and diffcyt 1.32.1, on 600 simulated count tables (7,600 clusters: 5 v 5 and 3 v 3 samples, 8 or 15 clusters, half with three clusters changed two-fold) and on four tables from a real two-group experiment (46 clusters).

| edgeR called as | Fold change (log2) | p value | Calls at FDR 0.05 |
|---|---|---|---|
| MIL-C imitates (common dispersion, QL test as before edgeR 4) | Within 0.0004 | Within 0.06% in 95% of clusters | 616 shared, 0 differ |
| MIL-C imitates (cluster-wise dispersion, 10 prior df, LRT) | Within 0.0003 | Within 2% in 95% | 724 shared, 1 differs |
| edgeR 4 defaults (QL) | Median 0.0005 | A factor of 2 in 5% of clusters | About one discovery in five differs |
| `diffcyt::testDA_edgeR` (LRT) | Median 0.00008 | A factor of 2 in 5% | About one discovery in ten differs |

The common dispersion is within 0.4% of edgeR's. The difference from the defaults is one of method: they estimate the prior degrees of freedom of the dispersions, MIL-C fixes them at 10. The comparison found two defects in MIL-C's F test (degrees of freedom not capped at the pooled ones; the mean variance taken geometrically), and both were corrected. Calibration on the same tables: edgeR has the same slight excess of small p values under no difference as MIL-C, and at its defaults a little more.

## Differential expression against limma

[[differential-expression|Differential expression]]: coefficients, moderated t and p agree with limma 3.68.5 to about seven significant figures on an 8 × 8 matrix in six designs (with and without a covariate, weights and missing values). With missing medians it follows `eBayes(legacy = TRUE)`; limma's newer estimator gives slightly different prior degrees of freedom (3.52 against 3.10 on the test matrix). diffcyt itself was not run.

## PhenoGraph against the authors' package

[[clustering-methods|PhenoGraph]] beside PhenoGraph 1.5.7 (with its Louvain and with Leiden), k = 30, resolution 1, three runs each:

- On mixtures of six Gaussians in 8 dimensions (populations from 40% to 1%): MIL-C's Leiden and the package's Leiden gave the same partitions on well-separated mixtures (ARI 1.00 and 0.997); on overlapping ones, 0.75 to 0.78, where the package agrees with itself at 0.75 to 0.77.
- On 20,000 real events with six markers (measured while MIL-C still used Louvain): the same number of clusters give or take two, and MIL-C agrees with the package (ARI about 0.7) as well as each agrees with itself from run to run.
- Leiden's modularity beside leidenalg's on the same graph: 0.8065 against 0.8036.
- MIL-C took 1 to 7 seconds a run, the package 27 to 183 (not a careful benchmark).

## Quality control against PeacoQC and flowAI

[[quality-control|Quality control]] beside PeacoQC 1.22.0 and flowAI 1.42.0:

- **Known defects** (synthetic and planted in real files): a step or a sudden drop is removed as PeacoQC removes it; a clog is removed entirely, where PeacoQC, which has no flow-rate check, removes nothing; of a slow drift, 76 to 97% is removed against PeacoQC's 74 to 98%. On defect-free stretches MIL-C removed 0 to 0.7% of the rest, PeacoQC 0 to 2%. flowAI, on raw values at its defaults, found few of the planted defects.
- **32 real files with no known defect:** MIL-C removes 5.5% of the events, PeacoQC 5.1%, flowAI 12.0% — not the same events: 56% of what PeacoQC removes MIL-C removes too. Three files of which PeacoQC removes 12 to 31% are left alone, no channel's median moving there by more than 0.05 of its spread. With no truth, this does not say who is right.
- **Copies of real files with a defect of known size:** a shift of 0.3 of the spread over 5% of the file is found in 24 of 32; PeacoQC found it in 8 of 8 smaller files but removed 6.3% of the same files with no defect, where MIL-C removes 0.04%. The real files lose the same with the drift limit anywhere from 0.20 to 0.33 and the shift limit from 0.20 to 0.30: the default 0.25 is the middle of both.

## Batch normalisation against CytoNorm

[[normalisation|Normalisation]] beside CytoNorm 2.0.13 on simulated batches with a known effect and on a real three-batch experiment. Without clusters it is CytoNorm's curve inside the range the references cover (mean difference about 1e-5 of the axis; a thousandth of the correction made). With clusters, the two programs are as far from each other as each is from itself with another seed. Beyond the outermost quantiles they differ on purpose: CytoNorm continues the spline's slope, MIL-C continues as a shift.

## Compensation from controls against AutoSpill

[[autocompensation|Autocompensation]] beside the AutoSpill R package (0.2.0) and `flowStats::spillover`: on the same events, the AutoSpill matrix to 2e-5 on simulated controls and 3e-4 on the real single-stain controls shipped with flowCore, in as many iterations. On simulated controls with a known matrix, MIL-C is within 0.04 percentage points (the classic method within 0.3). Nothing was measured on a modern multi-laser instrument, on beads or on a spectral one, and no matrix was set beside FlowJo's or FACSDiva's.

## Spectral unmixing against SpectroFlo and AutoSpectral

[[spectral-unmixing|Unmixing]], on a public Aurora bead dataset (42 controls, 64 detectors) with every file raw and as SpectroFlo 3.3 unmixed it: with SpectroFlo's spectra MIL-C reproduces its unmixed values to 1e-7; the spectra MIL-C takes from the same raw controls have a median cosine of 0.99994 with SpectroFlo's (lowest 0.966); the ratios of the beads' unmixed medians lie between 0.991 and 1.023. OLS, per-event WLS, the similarity matrix and the condition number agree with AutoSpectral's R functions to 1e-9 or better. Cytek's own WLS and complexity index are not published and were not compared.

## Cell cycle, proliferation and kinetics

- [[cell-cycle|Cell cycle]], on 223 public BrdU/PI files: fitted on PI alone, the aggregates found on ungated events match what a PI-A/PI-H singlet gate removes (24.8% against 25.6%, r = 0.90). On the 39 controls, Dean-Jett-Fox gives 6.6 points less S than the BrdU+ cells (SD 5.8) — BrdU also counts cells that left S during the pulse — and Watson 3.9 points more (SD 13.5).
- [[proliferation|Proliferation]] agrees with R's flowFit on its own CFSE and CellTrace Violet data within 0.5 points per generation; on its Cell Proliferation Dye sample they differ by up to 6.5 points.
- [[kinetics|Kinetics]] agrees with the same definitions written again in R to 1e-13. That checks the arithmetic, not a real calcium flux.

## Gates, scales and files

- Logicle matches FlowUtils to 15 decimal places (see [[scales-and-transforms]]).
- Curly quadrants, on a workspace MIL-C wrote with 40,000 simulated events, put the same events in every corner as CytoML 2.24.0 except up to 65 of a corner, all from cytolib's coarser table of the biexponential; with that table in its place, at most 1.
- A workspace written by MIL-C was read back by MIL-C and by FlowKit and gated the same events (32 samples, 192 populations), and was opened in FlowJo by hand before the 1.0 release. Gating-ML files validate against the standard's schema (see [[workspaces-and-gating-ml]]).
- [[auto-gates|Automatic gates]]: on test data, valley gates agree with hand-drawn ones at Jaccard 0.95 to 0.97; scatter gates come out looser.
- Prism files follow the R package `pzfx` and round-trip through MIL-C's reader; they were opened in Prism by hand before the 1.0 release (see [[tables-and-prism]]).
- Welch's and repeated-measures ANOVA reproduce textbook and R examples and independent calculations on random data (see [[cluster-statistics]]).
