---
title: Spectral unmixing
---

A spectral cytometer (Cytek Aurora and Northern Lights, Sony ID7000, BD FACSymphony A5 SE) has more detectors than fluorochromes. Each fluorochrome has a spectrum: how much of it each detector sees, scaled so that its highest detector is 1. An event's detector values are taken to be a sum of the spectra, and unmixing finds how much of each fluorochrome makes that sum. It takes the place of [[compensation]]: while unmixing is on, no compensation matrix is used.

## Where it is

**Unmixing…** in the [[gates-tab]], beside **Compensation…**. The dialog **Spectral unmixing** has the reference controls on the left and three tabs on the right.

## The reference controls

The list holds every file of the folder; **Add files…** adds controls kept elsewhere. Files whose name contains `(Beads)` or `(Cells)`, as SpectroFlo names them, start ticked, and the fluorochrome is guessed from the name (`A3 BUV395 (Beads)_Bead controls_1.fcs` → `BUV395`).

| Column | Meaning |
|---|---|
| **Use** | Tick each single-stain control to use. |
| **File** | The file. |
| **Fluorochrome** | Its name, which becomes the name of the unmixed channel. Required. |
| **Marker** | The stain on it, optional: shown as the channel's marker everywhere. |
| **Negative** | **Own dim**: the control's own dimmest events. **Unstained**: the unstained control's events. |

**Tick all** and **Untick all** do what they say. Below the list:

| Control | What it does |
|---|---|
| **Unstained control** | Unstained cells or negative beads. It is the negative of the controls set to **Unstained**, the floor of the weights for WLS, and the source of the autofluorescence spectrum. A file with `unstained` or `negative` in its name is chosen at the start, one of cells first. |
| **Autofluorescence as one more spectrum** | The unstained control's own spectrum is unmixed as a channel of its own, named `AF`. |
| **Method** | **Ordinary least squares (OLS)** or **Weighted least squares (WLS)**. |
| **Compute spectra** | Reads the controls and computes the spectra. |
| **Read CSV…** | Spectra from a table: one fluorochrome per row and one detector per column, or the other way round (told from where the detector names are). Each is scaled to a peak of 1. |
| **Save CSV…** | The spectra, one fluorochrome per row. |

The detectors are every channel of the first control but scatter, time and the like; of those, only the areas (`-A`) where the file has areas. The dialog does not let the detectors be chosen by hand; a CSV of spectra does.

Controls in the folder of samples are also samples. Untick them in [[samples-tab]] to keep them out of the concatenated events.

## How a spectrum is measured

For each control: the peak detector is the one where its brightest events (99.5th percentile) stand furthest above the negative. The positive events are the brightest 5% on that detector; the negative events are the dimmest 30% of the same control, or the unstained control's events. Events at 99% of a detector's range or above are left out, as saturation bends a spectrum. The spectrum is the positive median minus the negative median, in every detector, divided by its largest value. Fewer than 20 usable events, or positives no brighter than the negatives, is an error that names the control.

The autofluorescence spectrum is the unstained control's median in every detector, divided by its largest value.

## The two methods

- **OLS**: `a = y · pinv(M)`, one matrix for every event. Fast.
- **WLS**: each detector of each event is weighted by 1 / its variance, taken as the value itself (Poisson-like), never below a floor: the unstained control's robust spread (MAD) in that detector, and never below 125. A bright detector, which is noisy, then counts for less than a dim one. One small system per event, so it is slow: about 13 s for 222,000 events and 43 spectra, against 0.02 s for OLS. This is the per-event WLS of the AutoSpectral package, with its floor of 125.

## The tabs

| Tab | What it shows |
|---|---|
| **Spectra** | Every spectrum across the detectors. **Show** draws one over the others, which fade. |
| **Similarity** | The cosine similarity of every two spectra: 1 is the same shape, 0 nothing in common. Pairs above 0.98 are in bold, as Cytek warns above that. Under the table, the **complexity index** — the condition number of the matrix of spectra, largest singular value over smallest — and the most alike pair. The larger the index, the more each fluorochrome's noise spreads into the others. |
| **Check** | Each control unmixed with the whole set, its brightest 5% of events: **Peak** detector, **Positives**, **Unexplained** (the median share of the signal the spectra do not explain) and **Spill** with **Into** (the largest median that lands in another channel, as a share of its own, and which channel). Both near 0 means the set fits its controls. Click a row to show that spectrum. |

## Applying

**Apply** unmixes every sample that has the detectors of the spectra. The unmixed channels, named after the fluorochromes, take the place of the detectors (with their heights and widths); scatter, time and anything the spectra do not cover stay as they were. To the rest of the program they are ordinary compensated fluorescence: they get the biexponential axis, `Comp-` in gates, normalisation, statistics, clusters and figures.

Gates drawn on detectors that unmixing replaces would have nothing to cut, so **Apply** refuses while any exist and names them: remove those gates first. **Stop unmixing** compensates the samples again as before.

Like a new matrix, applying or stopping changes the values under the gates (a bar in the Gates tab asks for them to be reviewed) and makes concatenated events out of date. The unmixing is saved with the project.

## How it was checked

Measured on 23 November 2026 on a public Aurora data set (42 single-stain bead controls, 64 detectors) that holds each file both raw and as SpectroFlo 3.3 unmixed it (`docs/comparacao.md`, section 8):

- With SpectroFlo's own spectra, the OLS here reproduces the values SpectroFlo wrote to within 1.3 × 10⁻⁷ of the largest value.
- The spectra measured here from the same controls have a median cosine similarity of 0.99994 with SpectroFlo's; the worst was 0.966 (ViaKrome 808).
- OLS, WLS, the similarity and the condition number agree with the AutoSpectral package's functions to 10⁻⁹ or better.

Not checked: Cytek's own WLS and complexity index (its formula is not published, and the data set used OLS), files from the Sony ID7000 or the FACSymphony A5 SE, and tissue samples. Not available: a spillover-spreading matrix, per-cell or several autofluorescence spectra, and a different unmixing for each sample.
