---
title: The plate
---

A plate is one page of several panels, as a figure in a paper is: a map, a heatmap, two gating plots and a table of frequencies, lettered A, B, C… and written as one file. It can also be written as a report: a PDF whose first page says what the analysis was, followed by the plates, a page each. A project can hold as many plates as are needed, each with a name.

## Putting panels on a plate

Panels come from the Figure dialog (see [[figures]]). Make a figure look as it should, then press **Add to plate**, with the plate chosen in the list beside it — an existing plate, or **New plate**. The status bar says which letter it became. Open the plate with **View ▸ Plate…**.

A panel keeps what it showed and how it was styled: its labels, legend, colours and ramp. The page decides the type, the size and the resolution of every panel, which is what makes six plots one figure.

## The Plate window

The preview on the left is drawn by the code that writes the file. On the right:

| Part | What it does |
|---|---|
| **Plate** | The list of the project's plates, and **New**, **Rename…**, **Duplicate** and **Delete**. A new plate has the type and size of the current one and no panels. **Duplicate** makes a copy to change apart from the original. **Delete** asks first. |
| **Panels, in order** | The panels of the plate, lettered. A panel whose plot had nothing to show at the last refresh is marked stale. |
| **Move up**, **Move down**, **Remove** | Reorder the panels, or take one off. |
| **Refresh from the current analysis** | Draws each panel again from the plot it came from, as that plot is now. A panel whose plot has nothing to show keeps what it had and is marked stale. |

**The panel selected**:

| Control | What it does |
|---|---|
| **Title** | A title over the panel; none by default. |
| **Takes** | How many columns and rows of the grid it takes, 1 to 6 each. |

**Page**:

| Control | Default | What it does |
|---|---|---|
| **Panels per row** | 2 | Columns of the grid, 1 to 6. |
| **Fit** | | **Single column, 85 mm**, **1.5 columns, 120 mm**, **Double column, 180 mm**, or **Custom**. |
| **Width × height** | 180 mm × Automatic | **Automatic** height is three quarters of a cell's width for each row. |
| **Letters: A, B, C…** | on | A bold letter at the top left of each panel, 4 points larger than the type. After Z come AA, AB… |
| **Font**, **Size**, **Resolution** | | The type of every panel, and the resolution of a PNG or TIFF. |

Panels are placed in order, each in the first free place that holds it, so the letters read across and then down whatever the panels span.

## Writing it

- **Export…** writes the current plate as **SVG**, **PDF**, **PNG** or **TIFF**, by the file's extension.
- **Report…** writes a PDF. If several plates have panels, it first asks which go in, a page each. When events have been concatenated, the first page or pages say what the analysis was, from what the concatenated data records: how many events from how many samples, the population, quality control, the sampling (mode, seed, what was taken from each sample and what was left out), batch normalisation and the runs (UMAP, tSNE, FlowSOM and others) with their settings. Nothing there is recomputed or interpreted. The PDF records the version of MIL-C that wrote it.

## Saved with the project

The plates are written with the project as they are, in the folder's `MIL-C.json` and, for their arrays, `MIL-C plates.npz` beside it (see [[projects-and-saving]]). A panel read back from a project is what it was when it was saved: it no longer follows a plot, and **Refresh from the current analysis** leaves it alone. That is deliberate: the page of a paper must not change because the analysis behind it was run once more. To bring such a panel up to date, make the figure again and add it.

A panel that cannot be read back (a kind of plot this version does not know) is shown as missing on the page, and what was stored is kept, so saving again loses nothing.
