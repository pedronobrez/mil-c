---
title: The Samples tab
---

The Samples tab is the first of the window and the one the others depend on. It holds the table
of files — one row per `.fcs` file of the folder — with whatever you know about each sample that
the file does not record, and the panel on the right that decides which events go into the
concatenated set: from which population, from which samples, how many from each, and whether
batches are aligned first. Nothing is drawn from the files until **Concatenate** is pressed; until
then every number in the table is a plan you can read and change.

## Above the table

On the left, the folder open and, under it, where the gates and compensation come from: the
workspace's name, or `No workspace: every event of every file, compensated with the file's own
matrix`. On the right, five buttons in a row that wraps when the window is narrow:

| Button | What it does |
|---|---|
| **Add column…** | a new column of annotation: Day, Tissue, Donor, Treatment… |
| **Order…** | which value of a column comes first on every chart and table ([[sample-order]]) |
| **Fill from file names…** | reads columns out of the file names with a pattern (below) |
| **Workspace…** | reads a FlowJo `.wsp` and applies its gates and compensation ([[opening-data]]) |
| **Open folder…** | opens another folder |

## The table

| Column | Meaning | Editable |
|---|---|---|
| (box) | ticked: the sample takes part in the concatenation and the statistics | yes |
| **Sample** | the file name, without `.fcs` | no |
| **Group**, and every column you add | your annotation: free text | yes |
| **Events** | events in the file, from its header | no |
| **In gate** | events of the chosen **Population** in this file | no |
| **Takes** | events the plan draws from this sample | no |
| **Note** | why a sample gives nothing, or something to know about it | no |

To type into a column of yours, double-click a cell or select it and start typing. With several
rows selected (Shift-click, or `Ctrl`/`⌘`-click), what you type into one goes into all of them:
ten samples of one group are one entry, not ten. Values are trimmed of spaces at either end.

Double-click the header of a column of yours to rename it; right-click it for **Rename…** and
**Remove**. The table's own column names — Sample, Events, In gate, Takes, Note — cannot be given
to a column of yours, nor can two columns share a name; the status bar says so if you try.

The **Note** column can say:

| Note | Meaning |
|---|---|
| `sets the limit` | the sample whose size decided everyone's number ([[concatenation]]) |
| `gates from the group` | the workspace does not list this file; it was given the gates the listed ones share |
| `your gates` | the gates of this sample were drawn or changed in MIL-C |
| `counting…` | the population is being counted in this file |
| `N events, below the minimum of M` | fewer events than **Minimum**: left out |
| `no gates in the workspace` | the file has no gates, so it has no events of a gated population |
| `cannot be read: …`, `different channels from the others` | the file cannot take part ([[opening-data]]) |

## Fill from file names

**Fill from file names…** opens a dialog with a regular expression proposed from the first file
name, in which every named group `(?P<Name>…)` becomes a column called Name: one for each run
of digits and each run of letters, the first run of letters called Group. For a file
`7Treated3` (day 7, Treated, animal 3) the proposal is

```
(?P<Number1>\d+)(?P<Group>[A-Za-z]+)(?P<Number2>\d+)
```

Rename the groups to what they are — `(?P<Day>\d+)(?P<Group>[A-Za-z]+)(?P<Animal>\d+)` — and
delete the ones you do not want. `\d+` matches a number, `[A-Za-z]+` a run of letters. The preview shows what the pattern
finds in the first eight names; under it, how many of all the names it matches. A file the
pattern does not match keeps what it had. **OK** is enabled once the pattern has at least one
named group and matches at least one name. Columns that do not exist yet are added.

## The panel on the right

The settings scroll when the window is short; **Concatenate** stays in sight under them.

**Samples**

- **Population**: which events of each sample are concatenated — **All events**, or one gate of
  the tree, indented by depth. It has to exist in every sample that takes part; the list holds the
  gates every gated sample shares. Choosing one counts it in every file (the status bar shows
  `Counting …`), and fills **In gate**.
- **Include**: how many samples are ticked, as `8 of 12 samples`. Its menu ticks **Every
  sample** or **None**, and has a submenu per column (with between 2 and 40 values) to tick or
  untick every sample of one value at once — to concatenate two groups only, untick the others.
- **Gate the files the workspace does not list**: shown only when it applies; see
  [[opening-data]].
- **Group by**: which column, or columns together, make a group — for balancing the same total
  from every group, and for every comparison later. Tick two and each combination is a group:
  `Treated · 7`. With none ticked the button reads `Nothing: one group`.

**Downsampling** — **Balance**, **N**, **Minimum** and **Seed**, described in
[[concatenation]].

**Batch normalisation** — **Batches**, **Reference**, **Within**, **Bring to** and **Also in
Gates and Statistics**, described in [[normalisation]].

## The summary and the warnings

Under the settings, the plan in one paragraph: the total, the number of samples, and the share
of what they have, e.g. **72,000 events** from 6 samples — 31% of what they have; then the total
of each group, and `Limited by` the sample or group that set the number. In warning colour
below it:

- the plan's own notes (a sample short of **N**, samples within a group giving unequal numbers);
- files without gates that are left out of a gated population;
- two ticked files whose names differ only by a mark at the end — `sample_04` and `sample_04'`
  — which is how a re-acquisition of one tube is usually named: `look like the same animal
  acquired twice: untick one, or it counts double`.

## Concatenate

**Concatenate** draws the plan, builds the concatenated set, and switches to the Explore tab. It
is disabled while the plan takes nothing or a concatenation is running; the status bar shows
`Reading sample (3/12)` as it goes, and `Concatenated N events from M samples` at the end. See
[[concatenation]] for what is built, and [[explore-tab]] for what comes next.

Every change made in this tab is saved at once into the folder's project file
([[projects-and-saving]]), except the four settings of Downsampling.
