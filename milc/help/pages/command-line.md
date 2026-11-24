---
title: The command line
---

Everything the window does to a folder can be done without it: for a script, a server, or a run
that has to be repeated exactly. The command line uses the same session as the window — the
same gates, matrices, project file and quality control — so its numbers are the window's.

## Starting it

| Installed as | Command |
|---|---|
| from source, or `pip install` | `milc …` or `python -m milc …` |
| Linux `.deb` or `install.sh` | `milc …` |
| Windows | `MIL-C.exe …`, in the program's folder |
| macOS | `/Applications/MIL-C.app/Contents/MacOS/MIL-C …` |

With no command at all, the window opens. `milc --version` prints the version; `milc --help`
lists the commands, and `milc COMMAND --help` the options of one. Below, `DATA` is a folder of
`.fcs` files.

A command that cannot do what was asked prints why on the error stream and exits with status 2.

## inspect: what a folder holds

```
milc inspect DATA [-v]
```

One line per file: name, events, parameters, cytometer and date, from the headers alone.
`-v` / `--verbose` also lists every parameter with its stain.

## counts: every gate against FlowJo

```
milc counts DATA --wsp analysis.wsp [--listed-only] [--pattern REGEX] [--qc …]
```

Each sample's gating tree with MIL-C's count of every population beside FlowJo's, the
difference, and at the end the largest difference over all gates. Files the workspace does not
list are named at the end.

Options shared by `counts` and `concat`:

| Option | Meaning |
|---|---|
| `--wsp FILE` | the FlowJo workspace (required) |
| `--listed-only` | only the files the workspace lists; by default the others are given the gates and matrix those share |
| `--pattern REGEX` | columns read from each file name; every named group is a column: `'(?P<Day>\d+)(?P<Group>[A-Za-z]+)'`. Without it, one column, Group, guessed from the names or as typed in the window |

## concat: balanced downsampling and concatenation

```
milc concat DATA --wsp analysis.wsp [--population PATH] [--mode MODE] [--n N]
     [--min-events N] [--group-by COLUMNS] [--seed N] [--dry-run] [--out NAME]
     [--normalize-by COLUMN --reference COLUMN=VALUE [--normalize-clusters N]
      [--normalize-to BATCH]] [--qc …]
```

| Option | Default | Meaning |
|---|---|---|
| `--population` | the only end gate, if there is one | the gate's path, e.g. `Cells/Live/CD45` |
| `--mode` | `equal_per_sample` | `equal_per_sample`, `equal_per_group` or `ceiling` |
| `--n` | the data decides | per sample, or per group for `equal_per_group`; required by `ceiling` |
| `--min-events` | 1 | samples with fewer are excluded |
| `--group-by` | `Group` | column(s) that define a group, comma-separated |
| `--seed` | 0 | of the random draw |
| `--dry-run` | | print the plan only |
| `--normalize-by` | | the column that names each sample's batch |
| `--reference` | | the samples that should look alike in every batch, e.g. `Group=Control` |
| `--normalize-clusters` | 1 | fit within this many FlowSOM clusters; 1 aligns each marker whole |
| `--normalize-to` | the mean | a batch the others are brought onto, left as acquired |
| `--out` | `concat` | writes `NAME.npz` and `NAME.fcs` (and `NAME_qc.csv` with `--qc`) |

It prints the plan — each sample with its group, what it has and what it gives — the group
totals, and the share kept. See [[concatenation]] and [[normalisation]].

## analyze: UMAP, tSNE, FlowSOM, PhenoGraph, PCA

```
milc analyze concat.npz --channels CD8,CD4,CD19 [--umap] [--tsne] [--flowsom]
     [--phenograph] [--pca] [--out analysis]
```

Takes a `.npz` from `concat`, or any concatenated `.fcs` (FlowJo's too; `--sample-column`, default
`SampleID`, names the column that says which sample). Without `--channels` it lists the channels
and stops. Channels may be named by detector or by stain.

| Option | Default | Meaning |
|---|---|---|
| `--transform` | `biex` | `biex` (what FlowJo feeds its plugins), `arcsinh` or `logicle` |
| `--cofactor` | 150 | for `arcsinh` |
| `--neighbors`, `--min-dist` | 15, 0.5 | UMAP |
| `--perplexity` | 30 | tSNE |
| `--clusters`, `--grid` | 5, 10 | FlowSOM metaclusters and grid side |
| `--k`, `--resolution` | 30, 1.0 | PhenoGraph neighbours; above 1, more and smaller clusters |
| `--method` | `leiden` | or `louvain` |
| `--min-cluster-size` | 10 | smaller PhenoGraph clusters join their most connected neighbour |
| `--seed`, `--group-by` | 0, `Group` | |

It writes `NAME.npz` (the set with its results) and `NAME_events.csv` (every result per event),
and for each clustering the medians per cluster (`_mfi.csv`), frequencies per sample
(`_freq.csv`) and, where there are groups, their mean and SD (`_freq_mean.csv`, `_freq_sd.csv`);
PCA writes its loadings. See [[embeddings]] and [[clustering-methods]].

## stats: the statistics table

```
milc stats DATA [--wsp FILE] [--population "A;B"] [--channels CD8,CD4]
     [--stats Median,CV] [--out statistics.csv] [--prism FILE.pzfx]
```

A row per population of each sample: counts, percentages, and the statistics chosen —
`Median`, `Mean`, `Geometric mean`, `SD`, `CV`, `Robust SD`, `Robust CV`, `MAD` — of each
channel. `--prism` also writes a Prism file of one column (`--prism-value`, default `% of
parent`), grouped by `--group-by`, with a table for each value of `--within`. See
[[statistics-tab]] and [[tables-and-prism]].

## export: the events of one population

```
milc export DATA [--population PATH] [--format fcs|csv] [--raw] [--out exported]
```

A file per sample. `--raw` writes the values as acquired, not compensated; an FCS file then
carries the matrix.

## gates: the gates for another program

```
milc gates DATA [--wsp FILE] --out gates.wsp      # or gates.xml
```

A FlowJo workspace if the name ends in `.wsp`, otherwise Gating-ML 2.0 (`--format wsp|gatingml`
overrides the name). `--no-counts` skips opening the files to write each population's count.
Without `--wsp`, the gates come from the workspace in the folder or those saved beside the
files. See [[workspaces-and-gating-ml]].

## qc: the quality of every acquisition

```
milc qc DATA [--out qc.csv] [--plots qc.pdf] [--plot-format png|pdf] [--plot-signals 3]
```

For every file, what would be removed and what found it; nothing is removed. `--plots FILE.pdf`
writes one PDF with a page per sample; any other path is a folder that gets a file per sample.
See [[quality-control]].

The limits, on `qc` and wherever `--qc` is accepted: `--mad` (6), `--consecutive` (5),
`--events-per-bin` (auto), `--drift` (0.25), `--shift` (0.25), and `--qc-channels` (on `qc`
also `--channels`) for the parameters to judge. `--qc` on `counts`, `concat`, `stats`, `export`
and `gates` removes the bad stretches for that run only: Clean becomes the first population,
gate paths are written as before, and the project beside the files is not changed. A folder where
quality control was applied in the window is used as it is.

## autocomp: a matrix from single-stain controls

```
milc autocomp CONTROLS… [--method autospill|classic] [--gate] [--out spillover.csv]
```

The control files, or a folder of them: one per fluorochrome and optionally an unstained one.
`--assign FILE=DETECTOR` (repeatable) sets a control's detector, `FILE=unstained` marks the
unstained; `--detectors`, `--universal-negative` and `--max-events` (100,000; 0 for all) refine
it. See [[autocompensation]].

## selftest and gui

`milc selftest` runs the whole route on invented data and says whether this build is sound
(`--no-window` skips the window, `--log FILE` writes to a file and never raises); see
[[installation]]. `milc gui [FOLDER]` opens the window, with a folder if one is given.
