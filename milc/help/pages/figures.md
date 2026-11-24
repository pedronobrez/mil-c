---
title: Figures
---

Every plot that matters can be made into a figure for a paper: its size in millimetres, its type, its words, its colours and its legend, written as SVG, PDF, PNG or TIFF. The window's own plots are for looking; a figure is drawn again for printing, from the same data, so what is exported is what was on the screen, restyled. In SVG and PDF text stays text, so a label can still be edited in Illustrator or Inkscape.

## Where the Figure… buttons are

| Where | What becomes a figure |
|---|---|
| [[gates-tab]] | The plot shown: pseudocolour, contours, histograms, overlays, violins and boxes, with the gates drawn on it. |
| **Matrix…** and **Backgate…** in the Gates tab | The grid of every pair; the population on each gate above it. |
| **Model…** dialogs | The fit of [[proliferation]], [[cell-cycle]] or [[kinetics]]. |
| [[explore-tab]] | The map, coloured as it is shown. |
| [[clusters-tab]] | The heatmap, the frequencies by group, the tree. |
| [[abundance]] and [[differential-expression]] | The volcano plots and the heatmap of differences. |
| [[normalisation]] | **Figure: before…** and **Figure: after…**. |

## The Figure dialog

The preview on the left is drawn by the code that writes the file: what is seen is what is written. Under it a line gives the size, and the size in pixels a PNG or TIFF will have. Controls the figure has no use for (a colour ramp on a histogram, say) are not shown.

| Control | What it does |
|---|---|
| **Preset** | **Last used**, **Defaults**, or a look saved under a name. **Save as…** keeps the current look under a name; **Delete** forgets one. A preset carries the look, not this figure's own words (title, axis labels, legend title). |
| **Fit** | A size: **Single column, 85 × 65 mm**, **1.5 columns, 120 × 90 mm**, **Double column, 180 × 110 mm**, **Square, 85 × 85 mm**, or **Custom**. |
| **Width × height** | In millimetres, 20 to 600. |
| **Resolution** | For PNG and TIFF, 72 to 1200 dpi (300 by default). A vector file has none. |
| **Font**, **Size** | Any font the system has (Arial by default), in points (8 by default). |
| **Title**, **X label**, **Y label** | Empty: the plot's own labels. |
| **Frame** | **Left and bottom**, **Box** or **None**; with **Ticks** and **Grid**. |
| **Point size**, **Opacity** | Of the dots of a map; the opacity is that of the bars on a frequency chart. |
| **Points as a bitmap inside vector files** | For maps of very many events: text and axes stay vector, the cloud of points becomes an image, so the file can still be opened. |
| **Legend** | **Right**, **Below**, **Inside** or **Hidden**; **Counts** adds the number of events of each category. |
| **Legend title** | Empty: the plot's own. |
| **Colour ramp** | For anything coloured by a value. It starts on the ramp the window shows. |
| **Numbers in the cells** | On a heatmap. |
| **Gates, with their names and shares** | On a plot of the Gates tab, described below. |
| **Transparent background** | No white behind the figure. |
| Categories | Each category with its colour and label: double-click a colour to change it, type over a label to rename it in the figure only. |

**Export…** asks for a file; the format is taken from its extension or from the type chosen in the file dialog: **SVG**, **PDF**, **PNG** or **TIFF**. **Close** keeps the look as the last used. **Add to plate**, with the list beside it, puts the figure as it looks now on a plate: an existing one, or **New plate** (see [[the-plate]]).

The look is remembered from one figure to the next: whoever exports one figure for a manuscript is about to export six that must match. A new figure starts from the last look, with its own words empty and the ramp the window shows.

## Gates on a figure

A figure of the Gates tab carries the gates of the population shown, drawn as the plot draws them: each in black, with its name and its share of the population (`CD4+  41.2%`) on a pale ground, in the same place as on the screen. A quadrant's four corners each have their name and share. On a histogram the ranges are drawn and labelled at the top. The gates never move the ends of the axes. Untick **Gates, with their names and shares** for the plot alone.

Only the gates drawn on the plot are drawn on the figure: a gate on other parameters, or a combination, has no place on it. With several samples overlaid the gates are those of the first sample, as on the screen.

## Colours

Categories (samples, groups, clusters) take Paul Tol's qualitative colours, which stay distinguishable with the common colour-vision deficiencies. The order of the categories is the order chosen for the samples and groups (see [[sample-order]]). Ramps for values are the matplotlib ones (viridis, magma, cividis, turbo…), or a ramp of your own colours.
