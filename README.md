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
- Batch normalisation in the manner of CytoNorm: name the column that says which batch a
  sample was acquired in and the samples that should look alike in every batch; each batch is
  brought onto the others before the samples are concatenated, marker by marker or within
  FlowSOM clusters. Afterwards, look at every batch's references over one another, marker by
  marker, before and after. The batches meet at their mean, or one batch is the target and is
  left exactly as acquired. Optionally the samples themselves are normalised as they are read,
  so Gates and Statistics show normalised values too.
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
- The scale of each axis is yours: biexponential with its width and decades, Logicle (W, M, A
  and T; the same numbers as FlowUtils to 15 decimal places), arcsinh, logarithmic or linear —
  for one parameter or for every marker.
- Each sample is compensated with its own matrix: its file's, or the workspace's.
  **Compensation…** lists every sample with where its matrix comes from (file, workspace,
  edited, shared); type over any cell and apply it to the samples selected or to all, or go
  back to the original. The caption under the plot says which matrix is in use.
- When the values under the gates change after the fact — another matrix, or the samples
  themselves batch-normalised — the tab says the gates were drawn on the old values until you
  press Reviewed, and a line above the tabs says the concatenated events are out of date until
  they are concatenated again.
- Export the population shown, of the samples selected or of all, as FCS or CSV.
- A gate drawn on one sample is every sample's, as on a FlowJo group — or only this one's.
- What is drawn is kept beside the files (`MIL-C gates.json`) and is there next time.
- **Quality…** checks each acquisition in the order it was acquired for four things: a flow
  rate far from the usual one (a clog), an abrupt change of signal (a step, a bubble), a
  moderate shift that lasts (a drift: a population's level, followed over a long window, away
  from where it stays for longest), and a moderate shift that does not (a shift: a short
  stretch whose level stands apart from the bins on both sides). The bad stretches are removed
  as a first population, Clean, that every gate, statistic and concatenation then starts from.
  The limits are in the dialog — MAD limit (6), consecutive bins (5), events per bin, drift
  limit and shift limit (both 0.25 of the spread of the parameter's events; off does not look)
  — with Check again, and what is applied is recorded with the limits used. Where more than a
  fifth of a file would go, look at the plots before applying. No R needed; set beside PeacoQC
  and flowAI in `docs/comparacao.md`, with the limits measured on copies of real files with
  defects of known size.
- A **quadrant's arms can be moved one by one**: drag the small square on an arm to tilt that
  arm alone; the round handle at the centre moves all four. Until an arm is moved it is the
  plain quadrant, and old projects open unchanged.
- **Auto…** proposes a gate from the data, drawn dashed until you name it, and then an ordinary
  gate: the valley between the two peaks of a histogram (it proposes nothing where there is one
  peak only — use a control); one population on two parameters, the one you click or the
  largest, as a 90% density contour or the ellipse of a robust Gaussian fit; a band along the
  diagonal of area against height, for single cells. *Fit on each sample separately* fits the
  gate again on each sample's own events instead of copying the coordinates. On the test data
  the valleys agree with the hand-drawn Live and CD45 gates (Jaccard 0.95 to 0.97); the scatter
  gates are looser than the hand-drawn ones and should be looked at.
- **Model…** fits a model to the population shown, in every sample selected, and gives a table
  with a row for each sample (CSV, tests between groups) and a figure. *Proliferation* (CFSE,
  CellTrace Violet): the generations as Gaussians a halving apart, anchored on an unstimulated
  control if one is chosen; % in each generation and FlowJo's Division, Proliferation,
  Expansion and Replication Index and % Divided, each with a standard error. *Cell cycle* (PI,
  DAPI, Hoechst, linear): G1 and G2 as Gaussians and S by Watson's pragmatic method or
  Dean-Jett-Fox. *Kinetics* (calcium flux): a statistic of a parameter, or of the ratio of two,
  in slices of time; baseline, peak, time to peak, slope and area under the curve.
- **File ▸ Export gates…** writes the gates out: a FlowJo workspace (every sample with its
  file, compensation, scales and gating tree) or Gating-ML 2.0, the standard FlowKit,
  flowCore/CytoML and Cytobank read. What a format has no word for — a population made of a
  FlowSOM cluster — is left out with what is under it and listed in the report shown after.
- With several samples selected on a histogram, show each sample — or each group pooled — as
  a violin (median and quartiles marked) or a box (Tukey whiskers, at most 200 outliers a box).
- **Matrix…**: every pair of the chosen parameters for the population shown, as density with
  histograms on the diagonal — where spillover shows.
- **Backgate…**: the population shown, highlighted over its parent's events on each gate above
  it, with the outlines drawn.
- From a FlowJo workspace, quadrants, boolean combinations and ellipses are read as well as
  rectangles and polygons. Those three are read by rules taken from FlowKit and CytoML, not
  checked against FlowJo on our own files: such populations are in italics, with FlowJo's own
  count in the tooltip to compare.

**Statistics** — a row for each population of each sample.
- Count and percent of parent, of grandparent and of all events.
- For any parameter: median, mean, geometric mean, SD, CV, robust SD, robust CV, MAD — on
  compensated values in the file's own units, by FlowJo's definitions.
- The samples' own columns sit beside the name; sort by any header; export as CSV or copy
  into Excel or Prism.
- Compare groups on any column of the table, with the same tests and corrections.
- Where the same animals or donors were measured in every group, choose the column that names
  them under *Paired by*: the tests become the paired t or Wilcoxon signed-rank test for two
  groups and repeated-measures ANOVA (Greenhouse-Geisser corrected) or Friedman for more, on
  the subjects present in every group. Welch's ANOVA is offered for three or more independent
  groups with unequal variances.
- Send any column to GraphPad Prism as a `.pzfx` file: populations down, groups across, a
  subcolumn for each sample. The cluster frequencies go the same way from the Clusters tab.

**Explore** — UMAP, tSNE and FlowSOM, run in the background with FlowJo's defaults.
- PhenoGraph-style clustering finds the clusters without being told how many: nearest
  neighbours, Jaccard weights, Leiden communities (Louvain can be chosen), a resolution to make
  them coarser or finer, and a minimum size below which a cluster joins the one it is most
  connected to. At resolution 1 it splits large populations whether or not they have parts —
  as the original package does (15 clusters from both on the 3,488 test events): lower the
  resolution (0.5 is a good start, at the risk of losing a rare population in a neighbour),
  raise the neighbours, or merge the pieces in Clusters.
- PCA gives two axes that can be read, each with its share of the variance and the markers that
  weigh in it.
- Any two axes: an embedding, or two markers on a biexponential scale ticked in decades.
- Colour by group, day, animal, sample, cluster, or the expression of a marker.
- Split into one panel per group — or per group and day — with everyone else faint behind.

**Clusters** — what ClusterExplorer shows, linked.
- The map, the heatmap of expression, the FlowSOM tree and the frequency of each cluster by group
  (mean, standard deviation and every sample as a dot), in one panel or in one per value of a
  second column: the groups compared within each day.
- Click a cluster in any of the four and it is marked in the other three.
- The colour ramp is yours to choose: a named one, or two or three colours of your own.
- **Order…** (here and in Samples) says which group comes first — Ctrl before Case, or the
  other way — on every chart, legend, panel and table and in the Prism file; it is kept with
  the project. Left alone, values come as a person counts: day 3, 7, 14.
- Test the groups in each cluster: Welch's t or Mann-Whitney for two, ANOVA or Kruskal-Wallis
  for more, with p adjusted for the number of clusters and written on the chart.
- Shift+drag on the map outlines any set of events and reports the same numbers for it.
- Clusters can be renamed, annotated and merged (the algorithm's own clusters are kept
  underneath, so merges can be undone), names can be proposed from the markers each cluster is
  positive for, and a lasso can be kept as a cluster of its own. A FlowSOM cluster can be turned
  into a population of the full files, available in Gates and Statistics, by assigning every
  event to its nearest node of the trained map.
- Where both FlowSOM and PhenoGraph were run, choose which clustering the tab shows.
- **Differential abundance…**: two groups compared on each cluster's cell counts by a
  negative-binomial model in the manner of diffcyt's DA-edgeR (quasi-likelihood F test, optional
  adjustment for a batch or day), with a table and a volcano. It tests shares of the
  concatenated events, gives edgeR's numbers when edgeR is called the same way, and is not
  to be trusted with two samples a group.
- **Differential expression…**: for each cluster and marker, the median expression of each
  sample compared between two groups (diffcyt's DS-limma: a linear model and limma's moderated
  t, the same numbers as limma 3.68). Adjust for a batch or day, or pair by animal; choose the
  markers and the scale, the fewest cells a sample needs in a cluster (10), and the correction
  over all cluster–marker pairs. A table, a heatmap of the differences with stars, and a
  volcano. Prefer markers the clustering did not use; a median does not see a minority of
  cells switching a marker on.

**Figure…** — on every plot.
- Width and height in millimetres, resolution, font and size, title and axis labels, frame,
  point size and opacity, legend position, colour ramp.
- Recolour or rename any category: double-click its colour, type over its label.
- Save the look under a name and apply it to the next figure.
- The preview is drawn by the code that writes the file.
- **Add to plate**: several panels on one page (View ▸ Plate…), lettered A, B, C…, in one type,
  at a journal's column width, refreshed from the plots they came from when the analysis
  changes. There can be several plates, each with a name, and they are saved with the project
  (`MIL-C.json` and `MIL-C plates.npz`, to be kept together): a saved panel is a snapshot, and
  no longer follows the analysis once the project is reopened. A report writes the plates as a
  PDF behind a page stating the quality control, population, sampling plan, seed, normalisation
  and the settings of each run.

## What it does not do

- **Gating is not all of FlowJo's.** Rectangles, polygons, ellipses, quadrants, ranges and
  boolean combinations are drawn here and read from a FlowJo workspace. FlowJo's *curly*
  quadrants are not read: the workspace does not store their curve. The automatic gates
  propose; where a population has no valley beside it (Tcrb, in the test data) they propose
  nothing, and that gate still needs judgement or a control.
- **The models are simple ones.** Proliferation takes the generations to be exactly a halving
  apart; the cell cycle has no aneuploid population, no debris or doublet model and no Fox
  peak for a synchronised S.
- **Tests between groups stop short of mixed models.** One value per sample, compared by group,
  paired or not; counts by a negative-binomial model written here, not by edgeR itself.
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

**Installers** are on the releases page: a disk image for macOS (Apple Silicon), an installer
for Windows with a portable `.zip` beside it for machines where an installer cannot run, and
for Linux (x86-64, Ubuntu 22.04 or anything as recent) a `.deb` and a `.tar.gz`:

```bash
sudo apt install ./MIL-C-0.1.0-linux-x86_64.deb      # Debian, Ubuntu, Mint: menu entry and `milc`
tar -xzf MIL-C-0.1.0-linux-x86_64.tar.gz && MIL-C-0.1.0/install.sh   # any other, into ~/.local
```

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
python -m milc qc      DATA --out qc.csv --plots qc.pdf   # each acquisition over time
python -m milc counts  DATA --wsp analysis.wsp   # every gate, against FlowJo's count
python -m milc stats   DATA --channels CD8,CD4 --stats "median,robust cv" --out stats.csv
python -m milc export  DATA --population "Cells/Live" --format csv --out exported
python -m milc gates   DATA --wsp analysis.wsp --out gates.wsp   # or gates.xml: Gating-ML 2.0
python -m milc concat  DATA --wsp analysis.wsp --out concat \
    --pattern '(?P<Day>\d+)(?P<Group>[A-Za-z]+)' --group-by Day,Group
python -m milc analyze concat.npz --channels CD8,CD4,Foxp3 --umap --tsne --flowsom --out run
```

`concat --dry-run` prints the sampling plan and writes nothing. `analyze` also takes any
concatenated `.fcs`, FlowJo's included; `--phenograph` takes `--k`, `--resolution`, `--method
leiden|louvain` and `--min-cluster-size`, and `--transform` takes `biex`, `arcsinh` or `logicle`.

`qc` only looks. To leave the bad stretches out, give `--qc` (with `--mad`, `--consecutive`,
`--events-per-bin`, `--drift`, `--shift`, `--qc-channels` if the defaults do not suit) to `concat`, `stats` or
`export`: Clean becomes the first population for that run, gate paths are written as before,
and `concat` also writes `NAME_qc.csv`. A folder where it was applied in the window is used
as it is.

## Layout

```
milc/
  fcs.py            FCS 2.0 / 3.0 / 3.1, read and written with numpy alone
  compensation.py   spillover matrices, from the file or from the workspace
  transforms.py     linear, arcsinh, FlowJo's biexponential, Logicle; decade ticks
  gating.py         the gates — drawn, by time, by cluster — and the tree they form
  wsp.py            a FlowJo 10 workspace: files, matrices, transforms, gates
  gate_export.py    the gates made ready for another program, and what was left out
  wsp_write.py      write a FlowJo 10 workspace
  gatingml.py       write Gating-ML 2.0
  autogate.py       gates proposed from the data: valley, contour, ellipse, band
  sample.py         one file with its matrix and its gates
  sampling.py       the plan: how many events each sample gives, and why
  dataset.py        the concatenated events, and where each came from
  order.py          which group comes first: counting order, or the user's
  embed.py          UMAP, tSNE and PCA; how well a map keeps neighbours
  graphcluster.py   PhenoGraph-style clustering: neighbours, Jaccard, Leiden / Louvain
  clusters.py       clusters as the user has them: names, notes, merges
  qc.py             quality of an acquisition over time
  qc_plots.py       the flow rate and the signal against time, as figures
  abundance.py      differential abundance of clusters, on counts
  volcano.py        the volcano plot of it
  expression.py     differential expression within clusters, on medians
  expression_figures.py   its heatmap of differences and its volcano
  views.py          violins and boxes, the matrix of pairs, backgating
  modelfit.py       what the three models share: histogram, fit, uncertainty
  proliferation.py  generations of a dye diluted by division
  cellcycle.py      G1, S and G2 from DNA content
  kinetics.py       a parameter against time
  plate.py          several panels on a page, and the PDF report
  plate_store.py    plates kept with the project
  flowsom.py        the self-organising map, consensus metaclustering, the tree
  density.py        smoothed histograms and probability contours, for overlays
  cluster_stats.py  expression and frequency tables
  stats.py          the statistics table: counts, percentages, medians, CVs
  compare.py        tests between groups, and adjustment for many of them
  prism.py          tables in the shape GraphPad Prism takes, as .pzfx
  normalize.py      batch normalisation: quantile curves per batch, marker and cluster
  figures.py        a plot and a style into SVG, PDF, PNG or TIFF (matplotlib)
  session.py        the open analysis; what the window reads and writes
  cli.py            the command line
  ui/               the window (PyQt6, pyqtgraph). Nothing outside imports Qt.
run.py              what a packaged MIL-C starts from
packaging/
  make_icon.py      the mark and the application icon
  milc.spec         the PyInstaller build, for macOS, Windows and Linux
  make_dmg.sh       the macOS disk image
  milc.wxs          the Windows installer (WiX)
  make_linux.sh     the Linux package (.deb) and archive (.tar.gz)
tests/              pytest; test_real.py runs only where the acquisitions are
```

## Building the installers

```bash
bash packaging/make_dmg.sh      # macOS: builds MIL-C.app, tests it, wraps it in dist/*.dmg
```

The Windows and Linux ones cannot be built on a Mac: a bundle is the interpreter and every
compiled extension for the machine that built it. `.github/workflows/ci.yml` builds all three
on GitHub's machines (`packaging/make_linux.sh` wraps the Linux one, and the build installs the
package and the archive before it lets them out) — on a tag `v*`, which also creates the release, or on request from the Actions tab.
Every build is asked to run the whole route on invented data (`MIL-C selftest`) and to open its
window before it is wrapped.

## What was set beside the original, and what was not

Four methods here were written from the papers. On 23 November 2026 each was run beside the
software it takes after, on the same input; `docs/comparacao.md` has the tables, and
`tests/reference/` the scripts (they need R with Bioconductor and are not part of the tests).

- **Differential abundance** gives edgeR 4.10's numbers when edgeR is called the way this
  imitates: log fold changes within 0.0004, p values within 0.06% for 95% of clusters, the same
  calls at FDR 0.05 (0 of 616 not shared). The comparison found and fixed two defects in the F
  test. Against edgeR's present defaults and `diffcyt::testDA_edgeR`, 10 to 20% of discoveries
  are not shared: they estimate the prior degrees of freedom, this fixes them at 10.
- **Differential expression** gives limma 3.68's coefficients, t and p to 7 figures on an
  8 × 8 matrix in six designs. With missing medians it follows `eBayes(legacy=TRUE)`; limma's
  newer estimator for unequal degrees of freedom is not implemented, and diffcyt itself was
  not run.
- **PhenoGraph**: on 20,000 real events, this and the authors' package agree with each other
  as well as each agrees with itself from one run to the next (ARI about 0.7), with as many
  clusters give or take two. Leiden's modularity was set beside leidenalg's (0.8065 against
  0.8036 on the test events).
- **Quality control**: on a planted step it removes what PeacoQC removes; a clog it finds and
  PeacoQC, which has no flow-rate check, does not; of a slow drift it removes as much (76 to
  97% against 74 to 98%; before the slow pass, 53 to 81%). On the 32 real files it removes
  5.5% of the events (5.3% before the short pass), PeacoQC 5.1% and flowAI 12.0%, with no truth to say who is right, and
  not the same events: 56% of what PeacoQC removes this removes too. It leaves alone three
  files of which PeacoQC removes 12 to 31% with no channel moving, and takes more than a fifth
  of two (7Ctrl5, 14Case1) where several channels shift by a quarter to more than one spread.
  On copies of the real files with a defect of known size: a shift of 0.3 of the spread over
  5% of the file is found in 24 of 32 (20 without the short pass; PeacoQC found it in 8 of 8
  small ones, and removed 6.3% of the same files with no defect, where this removes 0.04%); a
  lasting shift of 0.3 in 80% of the cases, of 0.25 in 56%. The real files lose the same
  (5.5%) with the drift limit anywhere from 0.20 to 0.33 and the shift limit from 0.20 to
  0.30: 0.25 is the middle of both. Where a shift becomes a defect the data cannot say.
- **A workspace written here has not been opened in FlowJo.** It was read back by this program
  and by FlowKit and gated the same events: 32 samples, 192 populations, no event of
  difference. Its elements are those FlowJo 10.10 writes for samples, matrices, scales and
  gates; window layout, tables and cytometer defaults are not written. Open it in FlowJo and
  compare the counts before relying on it. Time is written in seconds, as FlowJo's own
  workspaces imply. Gating-ML files validate against the standard's schema, and FlowKit gates
  the same events from them.
- **Proliferation** agrees with R's flowFit on its own CFSE and CellTrace Violet data (each
  generation within 0.5 points); on its Cell Proliferation Dye sample the two differ by up to
  6.5 points. **Cell cycle** and **kinetics** were checked on simulations only, and neither
  was set beside FlowJo or ModFit. The test data has no proliferation dye and no DNA stain.
- **Batch normalisation**: without clusters it is CytoNorm 2.0's curve inside the range the
  references cover (apart by a thousandth of the correction on the real data). With clusters
  the two programs are as far from each other as each is from itself with another seed, and
  bring the reference samples equally close. Beyond the outermost quantiles they differ on
  purpose — CytoNorm continues the spline's end slope, this continues as a shift.
- **Logicle** was set beside FlowUtils, not beside FlowJo's own display. No workspace here
  declares a logicle scale: reading one was tested on XML written by hand.
- **Prism files** follow the layout of the R package `pzfx`, round-trip through this program's
  own reader, and have not been opened in Prism by the author.
- **Statistics** follow FlowJo's definitions and were not set beside FlowJo's numbers.
- **Quadrants, booleans and ellipses from a workspace.** Quadrant and boolean counts rest on
  arithmetic already verified for rectangles. Ellipses were checked against FlowJo once, on a
  100-event published file with linear axes; on biexponential axes nothing was measured, and an
  ellipse near zero may sit up to about one display channel off. FlowJo's count beside each is
  the check.
- **Repeated-measures and Welch's ANOVA** reproduce textbook and R examples typed from memory,
  and independent calculations on random data.

## Tests

```bash
python -m pytest
```

`tests/test_real.py` is asserted against acquisitions that are not in the repository and never
will be; it skips itself where the `Data/` folder is missing.

## Licence

**GPL-3.0-or-later** — see `LICENSE`. Copyright (C) 2026 Pedro Nobre.

MIL-C is free software: you may use, change and redistribute it, and anything distributed that
is built on it has to stay open under the same licence. It comes with no warranty.

`milc/ui/style.py` and `milc/ui/theme.py` come from MIL-Q, by the same author, where they are
under the MIT licence; here they are part of this program and carry its licence. The window is
built on PyQt6, itself GPL-3.0; the other dependencies are under permissive licences.
