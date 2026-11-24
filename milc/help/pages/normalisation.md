---
title: Batch normalisation
---

Samples acquired on different days differ for reasons that are not biology: a laser, a lot of
antibody, a voltage, a cleaning. If every batch carries reference samples that should look alike
— aliquots of one control, stained and acquired with each batch — then whatever separates the
references between batches is technical, and can be learnt and taken out of every sample of the
batch. MIL-C does this in the manner of CytoNorm (Van Gassen et al., Cytometry A 2020), before
the samples are concatenated.

## The settings

In the **Batch normalisation** section of the [[samples-tab]]. The columns offered are those of
the Samples table, so first add a column that says which batch each sample was acquired in.

| Control | Meaning |
|---|---|
| **Batches** | the column that names each sample's batch — a day, a run. **None** turns normalisation off |
| **Reference** | a column and one of its values: the samples that should look alike in every batch, e.g. Group = Control |
| **Within** | `each marker whole` (1), or a number of FlowSOM clusters (2 to 30) to fit the correction within |
| **Bring to** | `The mean of the batches`, or one batch, `Day 1, left as acquired`, that the others are brought onto |
| **Also in Gates and Statistics** | normalise the samples as they are read, not only the concatenated events |

The other controls are enabled once **Batches** names a column.

Choose the reference with care. Whatever separates the reference samples between batches is
removed from every sample of the batch — including any real difference between them. Aliquots
of one control sample are the safe choice; a control group of different animals is an assumption
that the group does not change between batches.

## Before and after concatenating

Under the settings, a line says what will be done:

> At the next Concatenate: 1: 2 reference sample(s); 2: 2 reference sample(s). Whatever
> separates the Control samples between batches will be removed from every sample — assuming
> they should not differ.

or what stands in the way:

- `Needs samples with Group = Control in at least two batches; there are such samples in 1.`
- `Batch 3, which the others are to be brought onto, has no sample with Group = Control.`

A reference sample counts only if it is ticked and has the chosen population.

After **Concatenate**, the line reports the result, for example:

> Normalised. References apart by 0.0220 before and 0.0135 after, on average over the markers
> (share of the axis, on reference events held out). Most moved: CD4, 0.0410 → 0.0120.

It adds how many batch-and-cluster pairs had too few reference events and were left as they
were, and which samples were not normalised because their batch had no reference.

## Looking at it: Before and after…

**Before and after…** appears once a normalisation has been done. It opens **Batch
normalisation: before and after**: for the **Marker** chosen, the reference events of every
batch drawn over one another as smoothed histograms, before (top) and after (bottom), each curve
scaled to its own peak. Where the curves lie on one another the batches agree.

The events shown are held out: half of each batch's reference events are set aside, the curves
are fitted on the other half, and "after" is measured on the half they never saw. Measured on the
events it was fitted on, the agreement would be perfect by construction and say nothing.
**Figure: before…** and **Figure: after…** make either panel into a figure ([[figures]]).

## What is done

For each batch and each marker:

1. The reference events of the batch are pooled — every event of the chosen population of every
   reference sample of that batch, read from the files, not only the ones drawn.
2. Their 99 quantiles, at 1% to 99%, are taken on the marker's own display scale (the one chosen
   in [[scales-and-transforms]], otherwise FlowJo's biexponential).
3. The goal is the mean of the batches' quantiles, or the chosen batch's own.
4. A monotone curve (a Fritsch-Carlson spline) through the pairs (batch quantile, goal quantile)
   is applied to every event of every sample of the batch. It never reverses the order of two
   events.

Only the compensated channels — the fluorescence markers — are normalised; scatter and time are
left as they are.

With **Within** above 1, a FlowSOM map (10 × 10, trained on the same number of reference events
from every batch) is cut into that many clusters, the curves are fitted within each cluster, and
each event is corrected with the curve of the cluster it falls in. This is for a shift that
differs between cell types. It needs many reference events: a batch with fewer than 50 in a
cluster is not corrected there and is listed in the result.

Beyond the outermost quantiles nothing is known about the batch effect, so the correction
continues as the shift it was at the last quantile, rather than along the spline's end slope.

## Also in Gates and Statistics

By default only the concatenated events are normalised, and the Gates and Statistics tabs show
the values as acquired. Ticking **Also in Gates and Statistics** normalises each sample as it is
read: plots, gates, counts and the statistics table then use normalised values. Gates drawn
before were drawn on the values as acquired, so the Gates tab asks for them to be reviewed
([[gates-tab]]), and the concatenated events are marked out of date until you concatenate again.
The curves are always learnt from the files as acquired, and are refitted at each
**Concatenate**, not while you edit gates: a plot does not shift under a gate as it is drawn.

## How it compares with CytoNorm

Set beside the R package (CytoNorm 2.0.13) on 2026-11-23, on the same events:

- Without clusters, the two corrections differ inside the references' 1%–99% by 1e-5 to 2e-5 of
  the axis on average, where the mean correction is 0.02 to 0.04: the same curve.
- With ten clusters the two programs' maps are not the same map, and their results are as far
  apart as two seeds of either; both leave the references equally close (0.013 to 0.014 of the
  axis, from 0.022).
- Beyond the references' range they part on purpose. CytoNorm's end slope moved single values by
  up to 1.1 of an axis running from 0 to 1; MIL-C's by 0.17 at most.

The details are in [[measured-facts]].

## From the command line

```
milc concat DATA --wsp analysis.wsp --normalize-by Day --reference Group=Control \
    --normalize-clusters 1 --normalize-to 1 --out concat
```

`--normalize-to` is optional (default: the mean of the batches). The command prints how far
apart the references were before and after. See [[command-line]].
