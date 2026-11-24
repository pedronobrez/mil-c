---
title: Welcome to MIL-C
---

MIL-C is a desktop program for flow cytometry. It takes a folder of `.fcs` files and, when there
is one, the FlowJo workspace that gates them; draws the same number of cells from every sample;
runs UMAP, tSNE, FlowSOM and PhenoGraph-style clustering on what it drew; compares groups; and
turns any plot into a figure at the size, in the type and in the colours a journal asks for.

It was written to replace a FlowJo route — gate, DownSample, concatenate, then the UMAP, FlowSOM
and ClusterExplorer plugins — that worked but left three things to be desired: figures that had
to be fought, a concatenation that silently favoured the sample with the most cells, and no
record of the seeds and settings behind a result. MIL-C is open source (GPL-3.0-or-later), free,
and needs no R installation. It reads FlowJo workspaces and was checked against FlowJo's
numbers, but contains none of FlowJo's code and is not affiliated with FlowJo or BD.

## The window

The main window has five tabs, read from left to right as the analysis goes:

| Tab | What it is for |
|---|---|
| **Samples** | which files take part, how they are grouped, how many events each gives, batch normalisation |
| **Gates** | one sample at a time: plots, the gating tree, compensation, scales, quality, models |
| **Statistics** | a row for each population of each sample, tests between groups, Prism files |
| **Explore** | UMAP, tSNE, FlowSOM, PhenoGraph and PCA on the concatenated events |
| **Clusters** | the clusters: map, heatmap, tree and frequencies, linked; abundance and expression |

## The manual, section by section

| Section | Pages |
|---|---|
| Getting started | [[installation]] · [[opening-data]] · [[projects-and-saving]] |
| Samples | [[samples-tab]] · [[concatenation]] · [[sample-order]] · [[normalisation]] · [[quality-control]] |
| Gates | [[gates-tab]] · [[drawing-gates]] · [[auto-gates]] · [[scales-and-transforms]] · [[compensation]] · [[autocompensation]] · [[spectral-unmixing]] · [[cell-cycle]] · [[proliferation]] · [[kinetics]] |
| Statistics | [[statistics-tab]] · [[abundance]] · [[differential-expression]] · [[tables-and-prism]] |
| Explore | [[explore-tab]] · [[embeddings]] |
| Clusters | [[clusters-tab]] · [[clustering-methods]] · [[cluster-statistics]] |
| Reference | [[figures]] · [[the-plate]] · [[workspaces-and-gating-ml]] · [[command-line]] · [[keyboard-shortcuts]] · [[measured-facts]] · [[design-principles]] · [[troubleshooting]] · [[glossary]] · [[version-history]] |

A first session usually reads [[installation]], [[opening-data]], [[samples-tab]] and
[[concatenation]], then [[explore-tab]] and [[clusters-tab]].

## Getting to the manual

- **Help ▸ Manual** opens this window. Its shortcut is `F1` on Windows and Linux and `⌘?` on
  macOS.
- The manual opens on the page for where you are: the tab shown, or the dialog in front. Press
  `F1` in the Quality dialog and you land on [[quality-control]]; in the Samples tab, on
  [[samples-tab]].
- **Help ▸ Export manual as PDF…** writes every page, in the order of the contents, as one A4
  document. Links between pages stay links inside the PDF. The **Export as PDF…** button in the
  manual's own toolbar does the same.

## Reading it

- The contents are on the left; the page on the right. **◀ Back** and **Forward ▶** move through
  the pages you have opened, and **Contents** comes back here.
- A link in the text opens another page. At the foot of every page, **Linked from** lists the
  pages that link to it, so you can go back the way a reader would have come.
- The search box (`Ctrl+F`, `⌘F` on macOS) looks through every page as you type. Each word
  matches by its beginning — `compens` finds compensation and compensated — and a page is listed
  only if it holds every word you typed, best match first; a word in a page's title counts for
  more. `Return` opens the first hit.
- The **English / Português** switch in the toolbar changes the language of the manual and is
  remembered. The program itself stays in English, so the Portuguese pages keep the names of
  buttons and menus in English.

## Conventions

- **Bold** is a control on the screen, written exactly as it is labelled: **Concatenate**,
  **Order…**, **Apply to every sample**.
- `Code` is something typed, or a file name: `MIL-C.json`, `milc concat`.
- A path through menus is written `File ▸ Open folder…`.
- Examples use invented names — files `sample_01.fcs`, groups Control and Treated.
- Where a number was measured (against FlowJo, an R package, simulated data), the page says so
  and where; the measurements together are in [[measured-facts]].

## What MIL-C is not

- Not all of FlowJo. The gates are rectangles, polygons, lassos, ellipses, ranges, quadrants
  and boolean combinations; FlowJo's layout editor, 96-well plate layouts and keyword editing
  are not here.
- Not an acquisition program. It reads FCS 2.0, 3.0 and 3.1 files written by a cytometer; it
  does not drive one.
- Not a statistics package. The tests compare one value per sample between groups, paired or
  not; mixed models are left to other software, and [[tables-and-prism]] says how to take the
  numbers there.
- Not FlowJo point for point. Its maps and clusters are close to FlowJo's in the ways that can be
  measured, not identical; [[design-principles]] explains the choices.
