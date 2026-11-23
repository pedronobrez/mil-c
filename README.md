<img src="packaging/icons/png/MIL-C-128.png" width="96" align="right" alt="">

# MIL-C

**Cytometry: balanced concatenation, UMAP, tSNE and FlowSOM, and figures that can be published.**

MIL-C takes a folder of `.fcs` files and the FlowJo workspace that gates them, draws the same
number of cells from every sample, runs the dimensionality reduction and the clustering, and
lets you look at the result from every side — then turns any plot into a figure at the size, in
the type and in the colours a journal asks for.

It is the cytometry member of the MIL suite, beside MIL-X (untargeted LC-MS) and MIL-Q (targeted
LC-MS), and looks and behaves like them.

> **Not affiliated with FlowJo or BD.** MIL-C reads FlowJo workspaces and was checked against
> FlowJo's results; it contains none of FlowJo's code.

## Why

The usual route — gate in FlowJo, DownSample, concatenate, run the UMAP, FlowSOM and
ClusterExplorer plugins — works, and leaves three things to be desired:

- **Figures.** Little control over legend, colours and type, and exports whose resolution has to
  be fought. Here every plot becomes a figure in millimetres, with editable text, as SVG, PDF,
  PNG or TIFF.
- **Balance.** An treated sample usually carries several times the cells of its control, and a map
  built on everything is a map of the treated group. Here the sampling plan is on the table
  before anything is concatenated: what each sample has, what it will give, and which one set
  the limit.
- **A record.** Seeds, parameters and the sampling plan are saved with the result.

## What it does

**Samples** — open a folder; a workspace beside the files is picked up.
- Gates and the compensation matrix come from the FlowJo workspace (`.wsp`). Files the workspace
  does not list are given the gates and the matrix the listed ones share, when they do share them.
- The table is yours to type in. It starts with one column, Group, guessed from the file names;
  add your own — day, tissue, donor — type a value into several rows at once, or fill columns
  from the names with a pattern. Any column, or several together, makes the groups.
- Three ways to balance: the same number from every sample, the same total from every group, or
  a ceiling per sample; and a minimum below which a sample is left out rather than allowed to
  drag every other down.

**Gates** — open one sample and look at it.
- Any two parameters as a pseudocolour or contour plot, or one as a histogram; markers on
  FlowJo's biexponential scale.
- Select several samples to overlay them: histograms, on one baseline or staggered, or contours
  — a colour for each sample or for each value of one of your columns.
- The gating tree with events and percent of parent; choose a population to see it.
- Draw a rectangle, a polygon, an ellipse or a range; drag its handles to move, stretch or turn
  it; rename; delete.
- A quadrant gate: one click where two lines cross makes four populations, and moving a line
  moves them all.
- Combine populations: the events in all of them, in any of them, or in none of them.
- The scale of each axis is yours: biexponential with its width and decades, arcsinh,
  logarithmic or linear — for one parameter or for every marker.
- See the compensation matrix in use and type over any cell; every sample is compensated again.
- Export the population shown, of the samples selected or of all, as FCS or CSV.
- A gate drawn on one sample is every sample's, as on a FlowJo group — or only this one's.
- What is drawn is kept beside the files (`MIL-C gates.json`) and is there next time.

**Statistics** — a row for each population of each sample.
- Count and percent of parent, of grandparent and of all events.
- For any parameter: median, mean, geometric mean, SD, CV, robust SD, robust CV, MAD — on
  compensated values in the file's own units, by FlowJo's definitions.
- The samples' own columns sit beside the name; sort by any header; export as CSV or copy
  into Excel or Prism.

**Explore** — UMAP, tSNE and FlowSOM, run in the background with FlowJo's defaults.
- Any two axes: an embedding, or two markers on a biexponential scale ticked in decades.
- Colour by group, day, animal, sample, cluster, or the expression of a marker.
- Split into one panel per group — or per group and day — with everyone else faint behind.

**Clusters** — what ClusterExplorer shows, linked.
- The map, the heatmap of expression, the FlowSOM tree and the frequency of each cluster by group
  (mean, standard deviation and every sample as a dot), in one panel or in one per value of a
  second column: the groups compared within each day.
- Click a cluster in any of the four and it is marked in the other three.
- The colour ramp is yours to choose: a named one, or two or three colours of your own.
- Shift+drag on the map outlines any set of events and reports the same numbers for it.

**Figure…** — on every plot.
- Width and height in millimetres, resolution, font and size, title and axis labels, frame,
  point size and opacity, legend position, colour ramp.
- Recolour or rename any category: double-click its colour, type over its label.
- Save the look under a name and apply it to the next figure.
- The preview is drawn by the code that writes the file.

## What it does not do

- **Gating is not all of FlowJo's.** Rectangles, polygons, ellipses, quadrants, ranges and
  boolean combinations are drawn here; from a FlowJo workspace only rectangles and polygons are
  read. No offset or curly quadrants, no automatic gates, no logicle scale, and one
  compensation matrix for all the samples rather than one each.
- **No statistics between groups.** Means and standard deviations, no test.
- **One figure at a time.** No multi-panel plates.
- **Its maps are not FlowJo's maps point for point.** See below.

## How close it is to FlowJo

Measured on one real experiment (BD FACSVerse, 12 gated samples, two concatenated files of 2,040
and 2,100 events, five markers). The tests in `tests/test_real.py` assert these with a margin.

| | MIL-C | FlowJo |
| --- | --- | --- |
| Gate counts, 72 gates | within 0.77% of FlowJo on every gate; within 3 events on the final population | — |
| Biexponential transform | medians per cluster within 0.69 of a channel (of 1,023) | — |
| FlowSOM, agreement with FlowJo's clusters (adjusted Rand) | 0.86 and 0.81 on average | — |
| FlowSOM, two runs of MIL-C with different seeds | 0.81 and 0.79 on average | — |
| UMAP, share of each event's 15 nearest neighbours kept | 47.0% and 47.9% | 47.2% and 47.7% |
| tSNE, same measure | 58.8% | 58.0% |

- FlowJo evaluates polygon gates on a 256-bin grid by a rule that was not reproduced; MIL-C uses
  the exact geometry. That is the 0.77%.
- FlowJo hands its plugins values on its biexponential scale, 0 to 1,023, not raw values. MIL-C
  feeds its algorithms the same by default.
- FlowSOM is stochastic. FlowJo's run is as close to MIL-C's as two of MIL-C's own are to each
  other.
- UMAP is `umap-learn`, not the `uwot` the FlowJo plugin calls: the pictures differ, the
  neighbourhoods they keep do not.

## Install and run

**Installers** are on the releases page: a disk image for macOS (Apple Silicon) and an installer
for Windows, with a portable `.zip` beside it for machines where an installer cannot run.
Nothing is signed. On macOS, right-click ▸ Open the first time; on Windows, choose "More info",
then "Run anyway". The first UMAP after installing takes a while longer than the ones after it:
the numerical code is compiled once and kept.

**From source**, Python 3.10 or newer.

```bash
pip install -e .
```

```bash
python -m milc
```

On macOS, double-clicking `MIL-C.command` does the same. To open a folder directly:

```bash
python -m milc gui path/to/folder
```

## The command line

Everything the window does, without the window.

```bash
python -m milc inspect DATA                      # what a folder holds
python -m milc counts  DATA --wsp analysis.wsp   # every gate, against FlowJo's count
python -m milc stats   DATA --channels CD8,CD4 --stats "median,robust cv" --out stats.csv
python -m milc export  DATA --population "Cells/Live" --format csv --out exported
python -m milc concat  DATA --wsp analysis.wsp --out concat \
    --pattern '(?P<Day>\d+)(?P<Group>[A-Za-z]+)' --group-by Day,Group
python -m milc analyze concat.npz --channels CD8,CD4,Foxp3 --umap --tsne --flowsom --out run
```

`concat --dry-run` prints the sampling plan and writes nothing. `analyze` also takes any
concatenated `.fcs`, FlowJo's included.

## Layout

```
milc/
  fcs.py            FCS 2.0 / 3.0 / 3.1, read and written with numpy alone
  compensation.py   spillover matrices, from the file or from the workspace
  transforms.py     linear, arcsinh, FlowJo's biexponential; decade ticks
  gating.py         rectangle and polygon gates, and the tree they form
  wsp.py            a FlowJo 10 workspace: files, matrices, transforms, gates
  sample.py         one file with its matrix and its gates
  sampling.py       the plan: how many events each sample gives, and why
  dataset.py        the concatenated events, and where each came from
  embed.py          UMAP and tSNE; how well a map keeps neighbours
  flowsom.py        the self-organising map, consensus metaclustering, the tree
  density.py        smoothed histograms and probability contours, for overlays
  cluster_stats.py  expression and frequency tables
  stats.py          the statistics table: counts, percentages, medians, CVs
  figures.py        a plot and a style into SVG, PDF, PNG or TIFF (matplotlib)
  session.py        the open analysis; what the window reads and writes
  cli.py            the command line
  ui/               the window (PyQt6, pyqtgraph). Nothing outside imports Qt.
run.py              what a packaged MIL-C starts from
packaging/
  make_icon.py      the mark and the application icon
  milc.spec         the PyInstaller build, for macOS and Windows
  make_dmg.sh       the macOS disk image
  milc.wxs          the Windows installer (WiX)
tests/              pytest; test_real.py runs only where the acquisitions are
```

## Building the installers

```bash
bash packaging/make_dmg.sh      # macOS: builds MIL-C.app, tests it, wraps it in dist/*.dmg
```

The Windows installer cannot be built on a Mac: a bundle is the interpreter and every compiled
extension for the machine that built it. `.github/workflows/ci.yml` builds both on GitHub's
machines — on a tag `v*`, which also creates the release, or on request from the Actions tab.
Every build is asked to run the whole route on invented data (`MIL-C selftest`) and to open its
window before it is wrapped.

## Tests

```bash
python -m pytest
```

`tests/test_real.py` is asserted against acquisitions that are not in the repository and never
will be; it skips itself where the `Data/` folder is missing.

## The mark

The suite has one animal — the cow's head, face on. MIL-X is the head alone, MIL-Q the head
inside a ring, MIL-C the head inside a **gate**: the polygon drawn around the population that
matters, with a few events left outside. The violet is the 405 nm laser. `packaging/make_icon.py`
draws it from the suite's shared shapes.

## Licence

**GPL-3.0-or-later** — see `LICENSE`. Copyright (C) 2026 Pedro Nobre.

MIL-C is free software: you may use, change and redistribute it, and anything distributed that
is built on it has to stay open under the same licence. It comes with no warranty.

`milc/ui/style.py` and `milc/ui/theme.py` come from MIL-Q, by the same author, where they are
under the MIT licence; here they are part of this program and carry its licence. The window is
built on PyQt6, itself GPL-3.0; the other dependencies are under permissive licences.

## Credit

FlowSOM is by Sofie Van Gassen and colleagues (Cytometry A, 2015); UMAP by Leland McInnes, John
Healy and James Melville; tSNE by Laurens van der Maaten and Geoffrey Hinton, here through
openTSNE. If MIL-C is useful to your work, cite their papers — the methods are theirs.
