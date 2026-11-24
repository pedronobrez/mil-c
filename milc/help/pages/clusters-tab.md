---
title: The Clusters tab
---

The Clusters tab shows the clusters from four sides at once: the **Map** says where a cluster is, the **Expression** heatmap what it expresses, the **Frequency** bars in which samples and groups it is found, and the **FlowSOM tree** how the map's nodes relate. A click on a cluster in any of them marks it in all four. It is empty until FlowSOM or PhenoGraph has been run in [[explore-tab|the Explore tab]]; how those work is in [[clustering-methods]].

## The first row: how it looks

| Control | What it does |
|---|---|
| **Clusters** | Which clustering is shown, FlowSOM or PhenoGraph. Appears only when both have been run |
| **Map** | What the map is drawn on: UMAP or tSNE where they have been run, or the first two markers the clustering ran on (e.g. `CD3 × CD4`) |
| **Heatmap** | **Median** or **Mean** of each marker in each cluster |
| **Colours** | The colour ramp of the heatmap, of the tree coloured by a marker and of expression on the maps: a named ramp or two or three colours of your own |
| **Colour per marker** | Ticked (default): each column of the heatmap runs over the whole ramp by itself, so a dim marker shows its own differences. Unticked: one ramp for the whole table |
| **Numbers** | The value of each heatmap cell written in it, here and in its figure |
| **Other markers too (n)** | Also show the stained markers the clustering was not run on. Says "(none)" and is off when it ran on all of them |
| **Tree** | Colour the FlowSOM tree by **Metacluster** or by any marker |
| **Clear selection** | Let go of the cluster or lasso selected |

## The second row: what is compared

| Control | What it does |
|---|---|
| **Compare** | The sample column whose groups the bars compare by colour. Starts on the grouping chosen in Samples. With several columns chosen there, their combination is offered too. "All samples" when the samples have no column that differs |
| **within each** | A second subdivision: one panel of bars for each of its values, so the groups are compared within each — the groups within each day, say |
| **Order…** | Which group comes first, on every chart, legend, panel and table: see [[sample-order]] |
| **Test** | **No test**, or a test of the groups' share of each cluster: **Welch t test / ANOVA**, **Mann-Whitney / Kruskal-Wallis**, **Welch t test / Welch ANOVA** |
| (correction) | Beside the test: **Benjamini-Hochberg**, **Holm**, **Bonferroni** or **None** |
| **Paired by** | The column that says which animal or donor a sample is of, where the same ones were measured in every group |
| **Results…** | Every comparison as a table, pair by pair where there are more than two groups |
| **Prism…** | The share of each cluster in each sample as a GraphPad Prism file: see [[tables-and-prism]] |
| **Differential abundance…** | Two groups compared on the clusters' cell counts: see [[abundance]] |
| **Differential expression…** | Two groups compared on what each cluster's cells express: see [[differential-expression]] |

What the tests do and how to read them is in [[cluster-statistics]].

## The four panels

Each panel has its title and a **Figure…** button that makes it a figure (see [[figures]]). The splitters between them can be dragged.

- **Map** — every concatenated event, coloured by its cluster. With a cluster selected, the others turn grey and stay where they are, so the cluster is read against the whole map.
- **Expression** — clusters down, markers across; each cell is the median (or mean) of the marker over all the cluster's concatenated events, on the scale the clustering ran on (0 to 1023 for the biexponential). Colour runs per marker unless **Colour per marker** is unticked.
- **Frequency** — for each cluster, the share of each sample's events in it: a bar per group at the mean, a whisker at one standard deviation, and every sample as a dot — five animals are five dots, not just a bar. The legend gives each group's n. With **within each**, one panel for each of its values, on one scale. With a test chosen, the adjusted p is written above each cluster where it is below 0.05, and "ns" where it is not.
- **FlowSOM tree** — a node for each map unit that holds events, as large as it is populated, joined by the minimal spanning tree of the nodes; coloured by metacluster or by a marker. Shown only for FlowSOM.

## Selecting

- **Click** a cluster in the heatmap, the bars or the tree to select it; click it again to let it go.
- **Shift+drag** on the map outlines any set of events (a lasso).

The line under the panels describes the selection: its number of events and share of all, the cluster's note, and its mean percentage of the sample in each group (± SD, n). For a lasso it also names the three clusters most of its events belong to.

## Working the clusters

The row of buttons at the bottom acts on the selected cluster. The algorithm's own result is never written over: names, notes, colours and merges are kept on top of it, saved with the run, so a merge can always be undone. A new run of the algorithm starts clean.

| Button | What it does |
|---|---|
| **Rename…** | Call the selected cluster something. The name is shown wherever its number was. Empty: back to its number. Two clusters cannot share a name |
| **Note…** | A note kept with the cluster, shown in the line under the panels |
| **Colour…** | A colour of your choosing for the cluster on the map, the tree and their figures, or the palette's own again |
| **Merge with…** | Make one cluster of the selected one and others ticked in a list. It takes the lowest number and keeps that one's name and note |
| **Undo merges** | Back to the clusters the algorithm gave; names and notes of merged clusters go too |
| **Suggest names…** | Propose a name for each cluster from the markers it is positive for — `CD4+ Foxp3+` — to accept, edit or refuse |
| **Make population** | The selected FlowSOM cluster as a population of the files, in Gates and Statistics |
| **Save selection…** | Keep the lassoed events as a cluster of their own; the other events stay in their clusters |
| **Remove** | Remove the selected hand-drawn cluster: its events go back to the clusters they were in |

### How names are suggested

For each marker the clustering ran on, the clusters' values are sorted and cut in two at the widest gap between neighbours: those above are positive, those below negative. A marker whose widest gap is narrower than a tenth of the 0–1023 biexponential scale (about half a decade), or 1 on an arcsinh scale, separates nothing and is left out. A cluster positive for none is named by the markers it is negative for (`CD4- CD8-`); two that the markers do not tell apart get ` (2)`, ` (3)` after the name. The dialog shows each cluster's present name beside the suggestion; type over any, or clear one to leave that cluster as it is, then **Apply**.

### A cluster as a population

**Make population** creates, in every file, a population named after the cluster (`FlowSOM 3` for a numbered one) under the population the concatenation was drawn from. Every event of every file — not only the ones drawn for the concatenation — goes to its nearest node of the trained FlowSOM map, and belongs to the population if that node is one of the cluster's. The population then appears in [[gates-tab|Gates]] and [[statistics-tab|Statistics]] like any gate. It is refused for a hand-drawn cluster (the map it was outlined on has coordinates only for the events drawn), for PhenoGraph clusters, when the trained map is not in memory (a concatenated file opened from disk does not carry it: run FlowSOM again), and when the map was trained on batch-normalised events but the files are not normalised.
