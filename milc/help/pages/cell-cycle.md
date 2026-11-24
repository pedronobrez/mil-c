---
title: Cell cycle
---

The cell cycle is read from the DNA content of each cell (propidium iodide, DAPI, Hoechst…) on a linear scale. Cells with one complement of DNA make a peak (G0/G1), cells with two make a peak at about twice the signal (G2/M), and the cells copying their DNA (S) lie between. A model of the histogram says how many cells are in each phase.

## Where it is

In the [[gates-tab]], choose the population (single cells, unless aggregates are modelled), select the samples in the list, then **Model…** ▸ **Cell cycle…**. The dialog fits every sample selected and shows one at a time.

The values are those the files store, compensated, on no display scale. The **Parameter** starts on the parameter on X if it is linear, otherwise on the first linear one.

## The settings

| Control | Default | What it does |
|---|---|---|
| **Parameter** | | The DNA signal. |
| **Model** | Watson pragmatic | **Watson pragmatic** or **Dean-Jett-Fox**, below. |
| **G2 / G1** | 1.90 to 2.10 | The mean of G2 is held between these multiples of the mean of G1: a dye is seldom exactly proportional to DNA. The same number twice fixes it. |
| **G1 near** | Found | About where the G1 peak is, in the units of the file, when the one found is not it. With a second population: where the diploid G1 is. |
| **G2 has the CV of G1** | off | Unticked, each peak has its own CV. Ticked, one CV for both: steadier where G2 is small. |
| **Second (aneuploid) population** | off | A second G1, S and G2 (Dean-Jett-Fox only). |
| **DNA index near** | Found | About where the aneuploid G1 is, as a multiple of the diploid G1 (1.5: half as much DNA again). |
| **Debris** | None | **None**, **Exponential** or **Single cut** (Dean-Jett-Fox only). |
| **Aggregates (doublets, triplets)** | off | Two and three cells stuck together (Dean-Jett-Fox only). |

The fit is run again a moment after any setting changes.

## The models

G0/G1 and G2/M are Gaussians in both. The histogram has 256 bins and is fitted between 4 standard deviations below G1 and 4 above G2. The number of cells in each component is found by non-negative least squares; only the positions and widths are searched for. The fit is by least squares on the counts, unweighted, as FlowJo's and ModFit's are.

- **Watson pragmatic** (after Watson, Chambers and Smith, 1987) assumes nothing of the shape of S. G1 is fitted on the side of its peak that S does not reach (the left), G2 on its right side, and S is what is left between the two means once both are taken away. The S lying under the sides fitted is allowed for, and the peaks are fitted again four times with it. It is less precise, since it uses half of each peak, and biased where S piles up against a peak.
- **Dean-Jett-Fox** (Dean and Jett, 1974; Fox, 1980) gives S a shape: a second-degree polynomial of DNA content between the two means, never negative, broadened by the same measurement error as the peaks. All three are fitted at once to the whole histogram. It is the more precise where S is smooth, and wrong where it is not.

With Dean-Jett-Fox, three more components can be fitted with the rest:

- **A second population**: its own G1, S and G2, its G1 at the DNA index times the first's (searched within 10% of where it is found or given), its G2 at the same multiple of its G1 as the first's, its own CV. The lower G1 is taken as diploid unless **G1 near** says otherwise. A second G1 near 2 cannot be told from the diploid G2: give **DNA index near** in that case.
- **Debris**: **Exponential**, a falling curve whose rate is fitted, or **Single cut**, nuclei cut once at a random point, worked out from the intact cells found (after Bagwell et al., 1991). The fit then reaches down to a fifth of G1.
- **Aggregates**: doublets and triplets, the distribution of single cells convolved with itself, as ModFit's aggregate model does. G1+G1 lies on G2, G1+G2 and G2+G2 above it. With this, single cells need not be gated.

Debris and aggregates take their shape from the cells found, so the fit is repeated four times, each time with the shapes the one before gave.

## The table

One row per sample, with your columns from [[samples-tab]]. The main columns:

| Column | Meaning |
|---|---|
| **% G1**, **% S**, **% G2** | Of the cells of the first (diploid) population, each with its **SE**. |
| **G1 mean**, **G1 CV**, **G2 mean**, **G2 CV** | Position and CV of each peak. |
| **G2/G1** | The ratio of the two means, with its SE. |
| **% sub-G1**, **% above G2** | Of all events, outside the fitted range (the dotted lines on the plot). |
| **RMS**, **RMS %** | The residual, in events per bin and as a percentage of the tallest bin. |
| **Chi²** | Reduced chi-square: near 1 when nothing is left but counting. |

With the extras: **% aneuploid** (of the cells of both populations), **DNA index**, the aneuploid phases and CV; **% debris** (of all events modelled); **% aggregates** (of events: an aggregate is one event).

The standard errors are those of counting: how far the same cells acquired again would move the number. They do not say whether the model is right; RMS and chi-square do. A sample that cannot be fitted has its reason in a **Note** column.

## The plot and the buttons

The plot shows the histogram, each component and the model of the sample chosen in the table, with a line above it summarising the fit. **Figure…** makes it a figure ([[figures]]); **Compare groups…** tests one column between groups of samples, each sample one observation ([[statistics-tab]]); **Export CSV…** writes the table; **Copy** puts it on the clipboard for Excel or Prism ([[tables-and-prism]]).

## How far to trust it

On histograms of 20,000 simulated cells with known phases: Dean-Jett-Fox is within 0.3 points at CVs of 3 to 5% where S is smooth, and varies by 0.5 to 1.5 points from one acquisition to the next; Watson is within about 1 point on average, with 1.5 to 3 points of scatter. Where S has a shape no parabola has, Dean-Jett-Fox can be off by 6 to 12 points, and its chi-square (9 to 30 instead of 1) says so; Watson's chi-square does not.

On 223 public BrdU/PI files, where BrdU gives each cell's phase without a model (`docs/comparacao.md`, section 6): in the 39 controls Dean-Jett-Fox gave 6.6 points less S than BrdU (SD 5.8), the expected direction, since BrdU also marks cells that passed through S during the pulse; Watson was 3.9 points over, with an SD of 13.5. Fitted without a singlet gate, the aggregates found (24.8% on average) matched what a singlet gate removes (25.6%; r = 0.90). Nothing was compared with ModFit or FlowJo, which were not available.
