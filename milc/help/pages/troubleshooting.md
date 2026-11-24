---
title: Troubleshooting
---

What goes wrong most often, what MIL-C says when it does, and what to do. The messages are
quoted as the program writes them; most appear in the status bar at the bottom of the window, or
in a line of warning colour inside the tab or dialog concerned. An error inside work done in the
background (a concatenation, a UMAP) is shown in the status bar as the kind of error followed by
its message, e.g. `ValueError: …`.

## The first start

**macOS says the application cannot be opened, or is from an unidentified developer.** The
installers are not signed. Right-click (or Control-click) **MIL-C** in Applications, choose
**Open**, and **Open** again in the dialog. If macOS offers only **Move to Bin**, close it, open
**System Settings ▸ Privacy & Security**, scroll to the line about MIL-C and choose **Open
Anyway**. It is asked once. From a terminal, `xattr -dr com.apple.quarantine
/Applications/MIL-C.app` does the same.

**Windows says it protected your PC (SmartScreen).** Choose **More info**, then **Run anyway**.
Where an installer is not allowed at all, use the `.zip`: unpack it and run `MIL-C.exe`.

**It does not start on an Intel Mac.** The macOS build is for Apple Silicon only; run it from
source instead ([[installation]]).

**It does not start on an older Linux.** The bundle needs a C library as recent as Ubuntu
22.04's. Install with the `.deb`, which brings the system libraries Qt needs, or run from
source.

**The first UMAP or FlowSOM takes a long time.** The numerical code is compiled the first time
and kept: measured in the macOS bundle, 16 s for the first UMAP, 5 s for each after it.

**Is this installation sound?** Run `MIL-C selftest` ([[command-line]]). On Windows add
`--log selftest.log` and read the file: a windowed program has no console.

## Opening files

| Message | What to do |
|---|---|
| `No .fcs files in` the folder | the files must be directly inside the folder chosen, not in subfolders |
| Note `cannot be read: not an FCS file` | the file is damaged or is not FCS; it is left out and the rest open |
| Note `cannot be read: $MODE …: only list mode is supported` | histogram-mode FCS files are not read |
| Note `cannot be read: $DATATYPE … is not supported` | an unusual number format; export the file again from the acquisition software |
| Note `different channels from the others` | the file's channels differ from most files' — often another program's concatenated export. Move it out of the folder |
| `Open a folder first` | a workspace is read for a folder: open the folder, then **Workspace…** |
| Note `no gates in the workspace` | the workspace does not list this file. Tick **Gate the files the workspace does not list**, or choose **All events** |
| Status `N population(s) not read (…)` after a workspace | those gates are of a kind MIL-C does not read; what is under them is missing too |

Populations shown in italics in the Gates tab were read from the workspace by rules not checked
against FlowJo; FlowJo's own count is in the tooltip to compare ([[workspaces-and-gating-ml]]).

## Concatenating

**Concatenate is greyed out.** The plan takes no events: no sample is ticked, the population is
not in any ticked sample, or the counts are still running (`counting…` in the Note column). Read
the Note column and the summary under the settings ([[samples-tab]]).

**One sample sets a very low number for everyone.** Its Note reads `sets the limit`. Raise
**Minimum** above its size to leave it out, untick it, or balance by group
([[concatenation]]).

**`… look like the same animal acquired twice: untick one, or it counts double.`** Two ticked
files have names that differ only by a mark at the end, the usual naming of a re-acquired tube.
Untick the one you do not want.

**A line above the tabs says the concatenated events were taken before something changed.** The
compensation, unmixing or batch normalisation changed after you concatenated. Explore and
Clusters still show the old values until you press **Concatenate** again.

## Batch normalisation

| Message | What to do |
|---|---|
| `Needs samples with Group = Control in at least two batches; there are such samples in 1.` | the reference value is missing from all but one batch, or those samples are unticked or lack the population |
| `Batch 3, which the others are to be brought onto, has no sample with Group = Control.` | choose another **Bring to**, or the mean of the batches |
| `… batch-and-cluster pair(s) had too few reference events and were left as they were` | lower **Within**, or use more reference events: it takes 50 per batch in each cluster |
| `Not normalised, their batch having no reference: …` | those samples' batch has no reference sample |

See [[normalisation]].

## Quality control

`Cannot be checked:` followed by the reason — `the file has no Time parameter`, `… events are too
few (it takes …)`, `Time is not a number for some events`, `Time is the same for every event`,
`Time runs backwards: several acquisitions in one file, or a clock that was reset`. Such a file is
kept whole. `That is more than a fifth of the file: look at the plots before applying.` is a
warning, not an error: something large went wrong, or the limits are too tight for this file.
See [[quality-control]].

## Gates and compensation

- `That gate needs two parameters: choose a Y axis` — only a range can be drawn on a histogram.
- `Every cell has to be a number` / `This matrix cannot be inverted: …` — in the Compensation
  dialog, a cell was left empty, or two rows are the same ([[compensation]]).
- `These files carry no compensation matrix` — after **Compensation…**: the files have no
  `$SPILL` keyword and no workspace gives them a matrix, so there is nothing to edit. Read a
  workspace that has one ([[compensation]]); `milc autocomp` can compute a matrix from
  single-stain controls ([[autocompensation]]).
- A bar in the Gates tab says the gates were drawn on values that have since changed. Look at
  them and press **Reviewed** ([[gates-tab]]).
- `No gate proposed: …` — the automatic gate found nothing to propose; the reason follows
  ([[auto-gates]]).

## Saving and exporting

| Message | Meaning |
|---|---|
| `Nothing concatenated yet` | Save needs a concatenation first |
| `Run FlowSOM first: the tables are per cluster` | `File ▸ Export tables…` needs FlowSOM |
| `No sample has gates to export` | draw gates or read a workspace first |
| `Nothing on the plot yet` | **Figure…** needs a plot |
| `cannot write .xyz: use .svg, .pdf, .png, .tiff` | a figure is written in one of those formats |

Changes are written to `MIL-C.json` in the data folder. If the folder is read-only (a shared
drive, a CD), copy it somewhere writable first ([[projects-and-saving]]).

## Still stuck

Report the problem at `https://github.com/pedronobrez/mil-c/issues`, with the version (`milc
--version`), the system, the message as written, and what you did before it. Do not attach data
you are not free to share; an invented file that shows the same problem is better. The terms
used in the messages are in the [[glossary]].
