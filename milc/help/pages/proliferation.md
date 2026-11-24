---
title: Proliferation
---

Proliferation is read from the dilution of a dye (CFSE, CellTrace Violet and the like). A cell that divides gives each daughter half its dye, so the generations are peaks a factor of two apart: evenly spaced on a logarithmic axis. A model of the histogram says how many cells are in each generation, and from those come FlowJo's indices.

## Where it is

In the [[gates-tab]], choose the population (live single cells, usually), select the samples in the list, then **Model…** ▸ **Proliferation…**. The dialog fits every sample selected and shows one at a time. The **Parameter** starts on the parameter on X. The values are those the files store, compensated, on no display scale; the model works on their logarithm, so events at or below zero are left out.

## The settings

| Control | Default | What it does |
|---|---|---|
| **Parameter** | the one on X | The dye. |
| **Undivided peak** | Fitted in each sample | A sample that was not stimulated: its peak says which peak of the others is the undivided one. |
| **Generations, at most** | 8 | Peaks in the model after the undivided one, 1 to 12. |
| **Autofluorescence** | 0 | The signal of cells without the dye, in the units of the file: only the dye above it halves. 0: none. |
| **Peaks widen with each generation** | off | Unticked, every peak has one width on the logarithmic axis. Ticked, the width may grow by the same fraction, up to 30%, with each generation. |

## The model

The histogram of log10 of the dye (256 bins) is fitted as a sum of Gaussians, one for each generation *g* = 0, 1, 2…, at

`position(g) = position(0) − g · log10(2)`

all of one width (or one that grows with *g*, when ticked). With **Autofluorescence** above zero, the dye above it is what halves: `position(g) = log10(background + (10^position(0) − background) / 2^g)`. The number of events in each generation is found by non-negative least squares; only the position of generation 0, the width and its growth are searched for.

Which peak is generation 0 is the hard part:

- **With a control** chosen in **Undivided peak**, a single Gaussian is fitted to the control's histogram, and generation 0 of every sample is looked for within a third of a division of it (log10(2) / 3). Cells that were stimulated and did not divide are seldom exactly as bright as cells left alone, but never a whole division away. The width is always the sample's own. The control's peak is drawn on the plot as a dashed line.
- **Without one**, every position is tried and, of those that fit as well as the best (within 5%), the dimmest is kept: the brightest peak that holds cells is taken to be the undivided one. That is wrong where every cell has divided. Use a control whenever there is one.

## The table

One row per sample, with your columns from [[samples-tab]]. *N(g)* is the number of events in generation *g*; the cells that were there at the start and whose descendants those are number *N(g) / 2^g*.

| Column | Meaning | Formula |
|---|---|---|
| **Generations** | The highest generation holding at least 0.5% of the modelled events | |
| **% Divided** | Share of the cells at the start that divided at least once (precursor frequency) | 100 · Σ(g≥1) N/2^g ÷ Σ(g≥0) N/2^g |
| **Division Index** | Divisions per cell at the start, those that never divided included | Σ(g≥1) g·N/2^g ÷ Σ(g≥0) N/2^g |
| **Proliferation Index** | Divisions per cell that divided | Σ(g≥1) g·N/2^g ÷ Σ(g≥1) N/2^g |
| **Expansion Index** | How many times the culture grew | Σ(g≥0) N ÷ Σ(g≥0) N/2^g |
| **Replication Index** | How many times the cells that divided grew | Σ(g≥1) N ÷ Σ(g≥1) N/2^g |
| **Undivided peak** | Where generation 0 is, in the units of the file | |
| **Peak SD (decades)** | The width of a peak | |
| **RMS**, **RMS %** | The residual, in events per bin and as % of the tallest bin | |
| **% modelled** | Share of the events the generations account for | |
| **% in generation 0**, **1**… | Share of the modelled events in each generation | |

These are FlowJo's definitions. R's flowFit calls the Expansion Index its "proliferation index". Each index comes with a standard error (**SE**) from counting: how far the same cells acquired again would move the number. It does not say whether the peaks are Gaussian; the RMS does. Cells that died are not seen, so every index is of the cells that are there to be counted.

## The plot and the buttons

The plot shows the histogram on a logarithmic axis, each generation and the model, for the sample chosen in the table; the line above it gives the generations, % divided, the two main indices, the RMS and where the undivided peak is. **Figure…** makes it a figure ([[figures]]); **Compare groups…** tests one column (the **Division Index** at first) between groups of samples ([[statistics-tab]]); **Export CSV…** and **Copy** take the table out ([[tables-and-prism]]).

## Limits

The model assumes Gaussian peaks of one width on the log axis, which holds for a well-stained, compensated dye. A dye that bleeds into a neighbouring channel, or a population with a broad undivided peak, gives peaks that are not Gaussian: look at the RMS and at the plot. Without a control, a culture in which every cell divided is read one generation off. No comparison with FlowJo's own proliferation platform has been made.
