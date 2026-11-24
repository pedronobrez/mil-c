---
title: Kinetics
---

Kinetics follows a parameter through the time of the acquisition, as in a calcium flux: cells run for a while, a stimulus is added, and the signal — or the ratio of two signals, as with Indo-1 — rises and falls. The acquisition is cut into slices of time, one number is taken from the events of each slice, and the series is measured between a baseline and a response.

## Where it is

In the [[gates-tab]], choose the population, select the samples in the list, then **Model…** ▸ **Kinetics…**. The files need a Time parameter. The **Parameter** starts on the parameter on X, unless that is Time itself. The values are those the files store, compensated, on no display scale.

## The settings

| Control | Default | What it does |
|---|---|---|
| **Parameter** | the one on X | What is followed. |
| **Divided by** | — | A second parameter: the ratio of the two is followed, event by event (violet over blue of Indo-1, say). Events whose denominator is not above zero are left out. |
| **Statistic** | Median | **Median**, **Mean** or **% above threshold** of each slice. |
| **Threshold** | 0 | For **% above threshold**: events above this, in the units of the file (or of the ratio), are the responding ones. |
| **Slice** | 100 slices | How long a slice of time is, in the units of the axis. Left at its lowest, the acquisition is cut into 100 slices. |
| **Smoothing** | Moving average, 5 | **None**, **Moving average** or **Gaussian**; the number is in slices: how many a moving average spans (made odd, centred), or the SD of the Gaussian. |
| **Baseline** | the first sixth | From and to: the slices before the stimulus. |
| **Response** | the rest | From and to: the slices after it. |

The **Baseline** and **Response** bands can also be dragged on the plot. All the samples selected are cut into the same slices, over the span of all of them.

Time is in seconds where every file says how long its unit of Time is (`$TIMESTEP`); otherwise in the files' own units, and the axis says so.

## The series

In each slice the statistic is taken over the events of the population. A slice with fewer than 5 events has no number. The series is then smoothed, each slice weighted alike and an empty slice skipped. The points on the plot are the slices, the line is the smoothed series, and each sample has its colour (or each value of the **Colour by** column of the Gates tab).

Each slice has a standard error: SD/√n for a mean, 1.2533 × robust SD/√n for a median (exact for a normal distribution, a guide otherwise), √(p(1−p)/n) for a percentage. The error of a smoothed slice follows from those of the slices in it.

## The table

| Column | Meaning |
|---|---|
| **Events** | Events of the population followed. |
| **Baseline**, **Baseline SD** | Mean and SD of the slices of the baseline, not smoothed: how much a slice moves by itself. |
| **Peak**, **Peak SE** | The highest value of the smoothed series in the response, with its standard error. |
| **Peak − baseline**, **Peak / baseline** | The rise. |
| **Time of peak** | On the axis of the plot. |
| **Time to peak** | Counted from the start of the response, where the stimulus is. |
| **Slope**, **Slope SE** | Of a straight line through the slices (not smoothed) from the start of the response to the peak, by least squares, per unit of time. Needs three slices. |
| **AUC** | Area between the smoothed series and the baseline over the response (trapezia): signal × time. |

Your columns from [[samples-tab]] come first, so the table can go straight to **Compare groups…** (it starts on **Peak − baseline**; see [[statistics-tab]]), **Export CSV…** or **Copy** ([[tables-and-prism]]). **Figure…** makes the plot a figure, with the bands ([[figures]]).

## How it was checked

On 23 public files with a Time parameter, the median of a parameter and of a ratio of two in 100 slices, a moving average of 5, baseline, peak, time of peak and AUC were computed again in R from the same definitions, with the files read by flowCore; the largest relative difference was 3.5 × 10⁻¹⁴ (`docs/comparacao.md`, section 6). That checks the arithmetic, not a real calcium flux: no public calcium-flux data could be obtained, and the R code repeats definitions written here.
