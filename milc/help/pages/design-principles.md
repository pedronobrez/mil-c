---
title: Design principles
---

MIL-C was written to replace one route through FlowJo — gate, downsample, concatenate, run the UMAP, FlowSOM and ClusterExplorer plugins — and to do three things that route did poorly: figures that can be published as they are, samples that weigh the same in the analysis, and a record of how every result was made. The choices below follow from that. They explain why the program behaves as it does where it differs from what you may expect.

## Balance is decided before anything is drawn

A treated sample often carries several times the cells of its control, and a map built on everything is a map of the treated group. So the sampling plan is on the table before anything is concatenated: what each sample has, what it will give, and which sample set the limit. There are three ways to balance — the same number from every sample, the same total from every group, or a ceiling per sample — and a minimum below which a sample is left out rather than allowed to drag every other down. The cost of balancing is read before it is paid. See [[concatenation]].

## The gates are the one definition of a population

A population is whatever its gates select, and every part of the program asks the gating tree: the [[statistics-tab|statistics]], the export of events, the concatenation, the [[quality-control|quality]] cleaning (which becomes a first population, Clean, that everything else starts from), and even a cluster, which becomes a population by a gate that places each event on the trained FlowSOM map (see [[clusters-tab]]). A gate drawn on one sample is every sample's unless you say otherwise. Gates come from the FlowJo workspace or are drawn here, are kept beside the files, and go out again as a workspace or Gating-ML (see [[workspaces-and-gating-ml]]).

## Nothing is changed behind your back

- The `.fcs` files are only read. Exports, concatenated files and projects are written as new files.
- An algorithm's result is never written over. Cluster names, notes, colours and merges are kept on top of it, so a merge can be undone and nothing the algorithm said is lost.
- When the values under something change after the fact — another compensation matrix, batch-normalised samples, a gate moved — the program says so instead of quietly recomputing or quietly keeping the old: the Gates tab marks the gates drawn on the old values until you press **Reviewed**, a line above the tabs says the concatenated events are out of date, and the statistics table says **Computed before the last change — Compute again**.
- A gate the program proposed is fitted again when its population changes; a gate you drew, or a proposed one you then moved, is not overwritten — a bar offers **Refit** or **Keep**.
- A p value is never invented. Where a test cannot be made, the cell is blank and the note says why.
- The order of the groups is an order of display. Changing it redraws; nothing is computed again (see [[sample-order]]).

## Numbers in the units they were measured in

Statistics are of compensated values in the file's own units, by FlowJo's definitions: a median fluorescence is a fluorescence, not a position on a display axis. Display scales are chosen where the data is shown or fed to an algorithm, and the concatenated events are stored untransformed. See [[scales-and-transforms]].

## The sample is the unit

Tests between groups take one value per sample. Five animals are five observations however many cells each gave; testing cells against cells would find every difference significant and say nothing about animals. Where many things are tested at once — every cluster, every population — the p values are adjusted for it. See [[cluster-statistics]].

## FlowJo's defaults, so a first result is comparable

The algorithms are fed what FlowJo feeds its plugins — markers on the biexponential scale, 0 to 1023 — and run with the plugins' defaults: UMAP with 15 neighbours and minimum distance 0.5, tSNE with perplexity 30, FlowSOM on a 10 × 10 grid. A first map can then be set beside the one it replaces. See [[embeddings]] and [[clustering-methods]].

## Measured, not claimed

Methods were written from the papers and then set beside the programs that made them known, on the same input. What was measured is written down, including where MIL-C differs and where nothing was checked; the comparison twice found defects, which were corrected. See [[measured-facts]].

## A record

Seeds, settings, markers, the sampling plan, the normalisation and the quality limits are saved with the result, and a report of the plates states them on a page of their own. The same seed and settings give the same result. See [[projects-and-saving]].

## The core does not need the window

The analysis — reading files, gating, sampling, the algorithms, the statistics, the figures — is a Python package that never imports Qt. The window is one way to drive it; the [[command-line]] is another, and does everything the window does without it. The same functions can be called from a Python script. This keeps the methods testable on their own and the results reproducible without clicking.

## No R, no licence, no network

UMAP, FlowSOM, PhenoGraph, differential abundance and expression, quality control and normalisation are written in Python, so no R installation is needed, unlike FlowJo's plugins. The program is free, runs offline and is native on Apple Silicon. Collaboration in the cloud was deliberately left out: it does not fit a desktop program that runs without the internet.

## Figures are the product

Every plot becomes a figure in millimetres, in the type and colours a journal asks for, with editable text in SVG and PDF, and the preview is drawn by the same code that writes the file. Several panels go on one page in [[the-plate]]. See [[figures]].

## One interface language, a translated manual

The interface is in English and stays in English, so a control is called the same thing in every lab, screenshot and paper's methods section. This manual is also in Portuguese; there, control names stay in English and are explained.
