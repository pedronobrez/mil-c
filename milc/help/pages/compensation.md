---
title: Compensation
---

A fluorochrome is seen not only in its own detector but, a little, in others. The spillover matrix says how much: row *i* is fluorochrome *i*, and the value in column *j* is the fraction of it that detector *j* sees. What the cytometer measures is the true signal times the matrix; compensation multiplies by the inverse of the matrix to get the true signal back. A compensated parameter is named `Comp-` followed by the detector in the gates, as FlowJo writes it.

## Where a sample's matrix comes from

Each sample is compensated with one matrix. In order of precedence:

| Source | Shown as | Meaning |
|---|---|---|
| Edited here for this sample | `edited` | A matrix changed in the dialog below and given to this sample alone. |
| Edited here for every sample | `shared` | A matrix changed in the dialog and given to all. |
| The FlowJo workspace | `workspace` | The matrix the workspace gives this sample, or the one every sample of the workspace shares (see [[workspaces-and-gating-ml]]). |
| The FCS file | `file` | The `$SPILLOVER` (or `SPILL`) keyword the cytometer wrote. |
| None | `none` | The file carries no matrix: its values are used as acquired. |

The caption under the bar of the [[gates-tab]] says which one the sample shown uses. With spectral unmixing on, no matrix is used at all (see [[spectral-unmixing]]).

## The Compensation dialog

**Compensation…** in the Gates tab opens it. If no sample has a matrix, the status bar says `These files carry no compensation matrix` and nothing opens.

On the left is the list of samples, with the source of each one's matrix in the **Matrix** column (bold when edited here). Under it, a line says how many different matrices the samples have, for example `3 different matrices among the 12 samples.`

On the right is the matrix of the sample selected, in percent: each cell is the percentage of the row's fluorochrome seen in the column's detector. The rows are labelled with the detector and its marker. The diagonal is 100 and cannot be edited. Type over any other number to change it; a comma is accepted as the decimal sign.

Above the table a line says what is shown: the sample's own matrix and where it comes from, or `… matrix with your changes — not applied yet`.

| Button | What it does |
|---|---|
| **Apply to** *sample* / **Apply to the N selected samples** | Gives the matrix in the table to the samples selected in the list; the others keep theirs. Select several with Shift or Cmd/Ctrl. |
| **Apply to every sample** | One matrix for every sample, in place of their own. |
| **Back to the original** | The samples selected go back to the matrix of their file, or of the workspace. Enabled when one of them has an edited matrix. |
| **Discard changes** | Forgets what was typed and shows the sample's matrix again. |
| **From controls…** | Computes a matrix from single-stain controls and puts it in the table, to apply like one typed by hand. See [[autocompensation]]. |
| **Close** | Closes the dialog. What was applied stays applied. |

Choosing another sample in the list while there are changes not yet applied keeps the changes in the table, so the same edit can be given to other samples.

A matrix is refused when a cell is not a number (`Every cell has to be a number`), when it cannot be inverted (`two fluorochromes would be indistinguishable`), or when a sample it is given to lacks one of its detectors.

## What happens when a matrix changes

- Every value compensated with it changes, so every count is computed again.
- The gates were drawn on the old values: a bar above the plot says so until **Reviewed** is pressed. Rectangles and polygons stay where they are in the new values and may no longer sit on their population. Gates fitted by **Auto…** are fitted again (see [[auto-gates]]).
- The concatenated events, if any, were taken from the old values; the window says they are out of date until they are made again (see [[concatenation]]).
- The edited matrices are saved with the project (see [[projects-and-saving]]) and are written into a FlowJo workspace exported with **File ▸ Export gates…**.

An FCS file exported with **Export…** and **Compensated values** unticked carries the matrix in use in its keywords; ticked, the values are written already compensated, in columns named `Comp-…` (see [[tables-and-prism]]).

## When to edit

Edit a matrix when a plot of two fluorescence parameters shows a population leaning along a diagonal that should be upright (under-compensated) or bending away from it (over-compensated). **Matrix…** in the Gates tab shows every pair at once, which is where spillover shows best. A typed change of a percentage point or two is usual; a large change usually means the controls should be acquired again. Check the scale too: on a biexponential axis a narrow **Linear width** can make a correctly compensated negative population look skewed (see [[scales-and-transforms]]).

The matrix files of different acquisition days are often different. The summary line of the dialog shows that at once; deciding which one is right is a matter of the controls of each day.
