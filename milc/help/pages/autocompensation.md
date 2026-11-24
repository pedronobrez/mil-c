---
title: Autocompensation from single-stain controls
---

A spillover matrix can be computed from single-stain controls: one FCS file per fluorochrome, each stained with that fluorochrome alone, and optionally an unstained one. In a control stained with fluorochrome *i* alone, every detector is a straight line against the control's own detector, its primary; the slope of that line is the spillover of *i* into that detector, which is row *i* of the matrix (see [[compensation]]).

## Where it is

In the [[gates-tab]], **Compensation…**, then **From controls…**. The dialog **Compensation from single-stain controls** computes a matrix for the detectors of the matrix shown in the Compensation dialog. Nothing is applied by computing it: the result goes into the compensation table, from where it is given to the samples chosen, like a matrix typed by hand.

The same is on the command line: `milc autocomp CONTROLS --out spillover.csv` (see [[command-line]]).

## The controls

| Control | What it does |
|---|---|
| **Add files…** | Choose the control FCS files. They are read off the window's thread. |
| **Remove** | Takes the selected controls out of the list. |
| Method | **AutoSpill (robust, refined)** or **Classic (positive vs negative)**, described below. |
| **Gate on scatter** | Keeps the main population on forward against side scatter, as the robust ellipse holding 95% of it ([[auto-gates]]): debris and clumps out. |
| **Unstained as negative** | Classic only: the unstained control's events are every control's negative, in place of each control's own negative events. |
| **Compute** | Computes the matrix from the controls in the list. |

The table lists each **Control**, its **Events**, its **Detector** and the detector it is **Brightest in**. Each control is given the detector it is brightest in, and no detector is given twice; change it in the **Detector** list where that is wrong, or choose **Unstained**. A file whose name says it is unstained (`unstained`, `no stain`, `blank`, `negative control`, or `US` or `neg` as a word) starts as the unstained one. Up to 100,000 events of each control are used, chosen at random.

Brightness is how far the brightest events of a control (its 99.5th percentile) stand above the unstained control's median in each detector, or above the control's own 20th percentile where there is no unstained one. The detectors are then matched to the controls so that each control gets the one where it is relatively brightest.

## The result

The status line says how many detectors were computed, and for AutoSpill how many iterations ran, whether they converged, and the largest spill left in the controls once compensated with the matrix found. Notes below list what went differently (a robust fit that did not converge, a control not gated on scatter and why, detectors with no control).

Two matrices are shown side by side, in percent: **From the controls** and the matrix the samples carry (titled after where it is: the sample's file or the FlowJo workspace). A cell in colour differs by more than 1 percentage point between the two; its tooltip gives the samples' value. Rows for detectors without a control are copied from the samples' matrix.

| Button | What it does |
|---|---|
| **Put in the compensation table** | Closes this dialog and types the matrix found into the compensation table, to four significant figures. Then select the samples and apply it there. |
| **Save as CSV…** | The matrix as proportions (not percentages), detectors in the first row and column, as a `$SPILLOVER` keyword holds them. |
| **Cancel** | Closes without changing the table. |

## The two methods

### Classic

What cytometer software has always done. In each control the valley between the negative and the positive events of the primary detector is found on a biexponential scale (two means, then the lowest point of a smoothed histogram between them; a valley higher than 80% of the lower peak is no valley). The positive events are those above it, the negative ones those below — the same cells, the same tube, the same autofluorescence. The top and bottom 1% are left out. The row is the median of the positives minus the median of the negatives, in every detector, divided by the same difference in the primary detector.

The unstained control's events are used as the negative where **Unstained as negative** is ticked, or where a control has fewer than 50 negative events of its own (beads often have none). A control with one population only and an unstained control is split at the unstained control's 99.5th percentile. At least 50 positive and 50 negative events are needed.

### AutoSpill

The method of Roca et al. (2021, *Nature Communications* 12: 2890), written from the authors' R package (`autospill` 0.2.0).

1. A robust linear regression (Huber's M-estimator, as `MASS::rlm`) of every detector on the primary one, over all the control's events with the top and bottom 1% left out, gives a first matrix. Where a robust fit does not converge, ordinary least squares is used.
2. The controls are compensated with it, the same regression on the compensated events measures the spill that is left, and the matrix is corrected by it — again and again. First on the linear scale until no slope left is above 0.01, then on the biexponential scale until none is above 0.0001, or until nothing improves any more. The step taken is the whole correction at first and a tenth of it once progress stalls; at most 100 iterations.

The package's thresholds, steps and stopping rules are kept as they are.

## How it was checked

Measured on 23 November 2026 against AutoSpill 0.2.0 and flowStats 4.24.0 (`docs/comparacao.md`, section 5):

- On the same events, the AutoSpill here gives the package's matrix to within 2 × 10⁻⁵ on simulated controls and 3 × 10⁻⁴ on the real controls that flowCore ships (a FACSCalibur, four detectors), with the same number of iterations.
- On simulated controls with a known matrix it finds the true matrix to within 0.04 percentage points; the classic method to within 0.3.
- With 15% debris in the controls, gating on scatter brought the AutoSpill error from 0.26 to 0.02 points.

Not checked: controls from a modern multi-laser instrument, beads, spectral instruments, the controls of the AutoSpill paper (they could not be downloaded), and no matrix computed here has been put beside FlowJo's or FACSDiva's for the same controls. Before trusting a matrix for an experiment, compare it with the one the acquisition software computed, and check the detector each control was given.
