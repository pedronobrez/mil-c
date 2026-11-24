---
title: Projects and saving
---

There is no **Save project** command in MIL-C. What you do to a folder is written, as you do
it, into a file beside the `.fcs` files; opening the folder again brings it all back. The
concatenated events, which are the expensive part to rebuild, are saved by hand.

## The project file: MIL-C.json

`MIL-C.json` sits in the folder of `.fcs` files. It is written whenever you change something —
type into the Samples table, tick a sample, draw a gate, choose a scale — and never just because
a folder was opened: a folder that was only looked at is left as it was. It holds:

| What | Set in |
|---|---|
| the columns of annotation, their values for every file, and which columns make the groups | [[samples-tab]] |
| which files are unticked | [[samples-tab]] |
| the population chosen for concatenation | [[samples-tab]] |
| the order of each column's values | [[sample-order]] |
| batch normalisation settings | [[normalisation]] |
| quality-control settings, and which samples were checked | [[quality-control]] |
| the gates of every sample whose gates were drawn or changed here | [[drawing-gates]] |
| the scale chosen for each parameter | [[scales-and-transforms]] |
| an edited compensation matrix, shared or per sample | [[compensation]] |
| spectral unmixing | [[spectral-unmixing]] |
| whether the gates still need a look after their values changed | [[gates-tab]] |
| the plates of panels | [[the-plate]] |

It is plain JSON: readable, and small enough to keep under version control or send with the data.
Files are named in it by their file name, so the folder can be moved or copied as a whole.

The plates keep their arrays in a second file, `MIL-C plates.npz`, beside the first. Keep the two
together. A folder analysed before the project file existed may hold `MIL-C gates.json`; it is
read when there is no `MIL-C.json`, and the next change writes a `MIL-C.json`.

What is not in the project: the FlowJo workspace (read, never written to), the concatenated
events and what was computed on them (save those below), and the settings of the Downsampling
section — **Balance**, **N**, **Minimum** and **Seed** start at their defaults in a new session.
The seed actually used is recorded in each saved concatenation.

The folder has to be writable. A folder on a read-only share opens, but changes cannot be kept.

## The concatenated events

**Concatenate** in the Samples tab builds a new set of events in memory. To keep it, with the
embeddings and clusterings run on it:

- `File ▸ Save concatenated file…` (`Ctrl+S`, `⌘S` on macOS). The dialog proposes `concat.npz`
  in the folder, and offers two formats:

| Format | Holds | For |
|---|---|---|
| MIL-C dataset (`.npz`) | the events, which sample and which event of the file each came from, every result (UMAP, tSNE, FlowSOM, PhenoGraph, PCA), how each was run, and how the set was built | reopening in MIL-C |
| FCS (`.fcs`) | the events, compensated channels named `Comp-…`, and a last column `SampleID` (1, 2, 3… in the order of the samples) | FlowJo or any other program |

- `File ▸ Open concatenated file…` reads either back. A `.npz` opens with its results, so the
  maps and clusters are there without running anything again.

With nothing concatenated yet, Save says `Nothing concatenated yet` in the status bar.

The `.npz` records how the events were drawn: the population, the balance mode, what each sample
had and gave, which were excluded and why, the sample that set the limit, the seed, and —
where they were applied — the quality-control and batch-normalisation records. The FCS file
carries the sample names (keyword `MILC SAMPLES`) and, where quality control was applied, what it
removed from each sample (keyword `MILC QC`).

## Other files MIL-C writes

| Command | Writes |
|---|---|
| `File ▸ Export tables…` | `cluster_expression.csv` and `cluster_frequency.csv`, after FlowSOM ([[cluster-statistics]]) |
| `File ▸ Export gates…` | a FlowJo workspace or Gating-ML 2.0 ([[workspaces-and-gating-ml]]) |
| **Export…** in Gates | the events of a population, one FCS or CSV file per sample ([[gates-tab]]) |
| **Figure…** on any plot | SVG, PDF, PNG or TIFF ([[figures]]) |

CSV files are written with a byte-order mark, so that Excel opens accented names correctly.

## When saved events go out of date

The concatenated events are a copy of values taken at one moment. If the values under them
change afterwards — a new compensation matrix, other unmixing, batch normalisation applied to
the samples themselves (or turned off again) — a line appears above the tabs:

> The concatenated events were taken before … : Explore and Clusters still show the old values.
> Concatenate again in Samples to bring them up to date.

Nothing is recomputed behind your back; the line stays until you press **Concatenate** again.

## Settings of the program

What belongs to you rather than to a folder — the appearance (`View ▸ Appearance`), the folder
opened last, the language of the manual — is kept in the system's usual place for application
settings, under the name MIL-C.
