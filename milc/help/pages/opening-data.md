---
title: Opening data
---

MIL-C works on a folder: every `.fcs` file in it is a sample, a FlowJo workspace beside them
gates them, and what you do to them is kept in the same folder (see [[projects-and-saving]]).
There is no import step and no copy of the files: they are read where they are, and never
written to.

## Opening a folder

- `File ▸ Open folder…` (`Ctrl+O`, `⌘O` on macOS), or the **Open folder…** button at the top
  right of the Samples tab. The dialog starts in the folder opened last time.
- From the command line, `milc gui path/to/folder` opens the window with the folder already open.

Every file ending in `.fcs` directly inside the folder is listed; subfolders are not looked into.
The files are sorted as a person counts — `3Ctrl1` before `14Ctrl1`, `sample_2` before
`sample_10` — and only their headers are read at first: the number of events, the channels and
the stains. The events themselves are read when a sample is first needed, by a count, a plot or
a concatenation. The status bar says how many files were found, or `No .fcs files in` the folder.

MIL-C reads FCS 2.0, 3.0 and 3.1 files in list mode, the format every cytometer writes.

## Files that are left out

Two kinds of file are listed in the Samples table but cannot take part. Their box is unticked
and greyed, and the **Note** column says why:

| Note | Why |
|---|---|
| `cannot be read: …` | the header or the data could not be read; the rest of the message is the reason, such as `not an FCS file` or `$DATATYPE … is not supported` |
| `different channels from the others` | the files of an experiment share their channels; one that does not — another program's concatenated export, say — is set aside. The set of channels most files have is the one kept |

A broken file never stops the folder from opening.

## The first guess at groups

A new folder starts with one column of annotation, **Group**, filled by a guess from the file
names: the run of letters that takes few different values across the files. `Control_1`,
`Control_2`, `Treated_1`, `Treated_2` give Control and Treated. Where no run of letters
qualifies the column starts empty. The guess is only a starting point: type over it, add other
columns, or read them out of the names with a pattern, as [[samples-tab]] describes.

## The FlowJo workspace

If the folder holds a `.wsp` file, the first one in alphabetical order is read with the folder.
To use another, or to add one later, use **Workspace…** in the Samples tab or `File ▸ Open
workspace…`; a folder has to be open first, otherwise the status bar says `Open a folder first`.

From the workspace MIL-C takes, for every file it names (matched by file name):

- the gating tree, with FlowJo's count for every population kept to compare;
- the compensation matrix that sample was given.

The axes are drawn on MIL-C's own scales — FlowJo's biexponential for fluorescence, the values
as stored for scatter and time — until you choose others ([[scales-and-transforms]]).

Rectangles, polygons, ranges, quadrants, boolean combinations, ellipses and FlowJo's curly
quadrants are read. The last four are read by rules taken from FlowKit and CytoML, not checked
against FlowJo on real files; populations read that way are shown in italics in the Gates tab,
with FlowJo's own count in the tooltip. After a workspace is read, the status bar reports
anything worth knowing: populations of a kind that could not be read, and how many populations
were read by unverified rules. [[workspaces-and-gating-ml]] goes into the details.

Under the folder's name, the Samples tab says what is in use: `Gates and compensation from`
the workspace, or, without one, `No workspace: every event of every file, compensated with the
file's own matrix`.

## Files the workspace does not list

A workspace often lists only some of the files in the folder — a tube acquired after the
analysis was set up, for instance. When every sample of the workspace has the same gates and the
same compensation matrix, the checkbox **Gate the files the workspace does not list** appears in
the Samples tab. Ticked (the default), the other files are given those gates and that matrix,
which is what FlowJo does when a file is dropped into a group; their Note says `gates from the
group`. Nobody has looked at those gates on those files, so look at them in [[gates-tab]].
Unticked, they have no gates, and can only take part when **Population** is **All events**; with
a gated population their Note says `no gates in the workspace`.

## What is concatenated by default

The **Population** of the Samples tab starts at the deepest gate that every gated sample has —
the end of the longest chain, such as `Cells/Single/Live/CD45`, rather than a gate drawn beside
it. The one chosen last time is remembered with the project, as long as the gates still have
it.

## Opening a concatenated file

`File ▸ Open concatenated file…` opens events that were concatenated before, skipping the
Samples tab:

- a `.npz` saved by MIL-C (`File ▸ Save concatenated file…`), with everything computed on it;
- a concatenated `.fcs`, MIL-C's own or FlowJo's. Channels named `Comp-…` or `FJComp-…` are
  taken as compensated, and a `SampleID` column, where there is one, says which sample each event
  came from.

The window then moves to the Explore tab. See [[concatenation]] and [[embeddings]].

## Next

With the folder open, the table of [[samples-tab]] is where the files are grouped and chosen;
[[quality-control]] can clean each acquisition before anything is drawn from it.
