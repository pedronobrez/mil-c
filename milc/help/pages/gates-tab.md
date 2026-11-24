---
title: The Gates tab
---

The Gates tab shows one sample at a time: any two of its parameters, or the histogram of one, for one population of its gating tree. It is where gates are drawn, moved, renamed and deleted, where the compensation and the scale of each axis are set, and where the models of one parameter (proliferation, cell cycle, kinetics) are fitted.

A gate drawn here belongs to the population being looked at, as in FlowJo. While **Apply to every sample** is ticked it also belongs to every other sample, the way a gate on a FlowJo group does. What is drawn is kept with the project beside the files (see [[projects-and-saving]]) and takes the place of the workspace's gates from then on.

## The left side

| Part | What it is |
|---|---|
| **Samples** | The usable files of the folder, each with its number of events. Click one to look at it. Select several (Shift-click, or Cmd/Ctrl-click) to overlay them. |
| **Populations** | The gating tree of the sample shown, with three columns: **Population**, **Events** and **% parent**. Click a population to look at it. |
| **New gate on the population shown** | The nine buttons that draw and edit gates, described in [[drawing-gates]]. |
| **Apply to every sample** | Ticked (the default), a gate drawn, moved, renamed or deleted on this sample is drawn, moved, renamed or deleted on all of them. Unticked, only here. |
| **Back to the workspace's gates** | Forgets every gate drawn or moved here and goes back to the gates of the FlowJo workspace, or to no gates at all. Enabled once something was changed. |

A population shown in *italics* was read from a FlowJo workspace by a rule that has not been checked against FlowJo itself; its tooltip gives FlowJo's own count beside the count here. See [[workspaces-and-gating-ml]].

When a population is chosen, the axes jump to the two parameters of its first gate, so a population is first seen the way its gates were drawn. After a gate is drawn, the axes stay where they are.

## The bar above the plot

| Control | What it does |
|---|---|
| **X**, **Y** | The parameters on each axis, named by their marker where the file gives one. **Y** ends with **Histogram**, which shows the parameter on X alone. |
| **Scale…** (beside X and beside Y) | The scale of that axis: biexponential, Logicle, arcsinh, logarithmic or linear. See [[scales-and-transforms]]. |
| **Plot** | **Pseudocolour** or **Contour**, for two parameters and one sample. |
| **Colour by** | With several samples overlaid: a colour for each **Sample**, or for each value of one of your columns (see [[samples-tab]]). |
| Overlay choice | With several samples on a histogram: **Overlaid**, **Staggered**, **Violins** or **Boxes**, and **One for each group** to pool the samples of each value of **Colour by** into one violin or box. |
| **Quality…** | Looks for stretches of the acquisition that went wrong. See [[quality-control]]. |
| **Compensation…** | The spillover matrix of each sample. See [[compensation]]. |
| **Unmixing…** | Spectral unmixing in place of compensation. See [[spectral-unmixing]]. |
| **Export…** | The events of the population shown, as FCS 3.1 or CSV files. See [[tables-and-prism]]. |
| **Figure…** | The plot as a figure for a paper. See [[figures]]. |

Under the bar is a caption that says what is drawn, for example `sample_01 · CD4+: 12,034 events, 41.2% of Lymphocytes — 40,000 drawn · compensated with its file's matrix`. It also says when samples are overlaid (and whose gates are shown: the first sample's), when values are batch-normalised, and which matrix compensates the sample.

## What the plot shows

- **Pseudocolour**: each event is a dot coloured by how crowded its neighbourhood is, on a logarithmic ramp (a 160 × 160 grid over the plot). The densest dots are drawn on top. Above 40,000 events a random 40,000 are drawn (always the same ones); every count and percentage still uses all the events.
- **Contour**: lines that enclose 90, 75, 50, 25 and 10% of the events; the events outside the outermost line are drawn as dots.
- **Histogram**: 256 bins across the axis, the count of events on the left. With several samples, each is scaled to its own tallest bin (`% of the tallest bin`); **Staggered** puts them one above the other.
- **Violins** and **Boxes**: each sample's distribution of the parameter on X, one row each, the axis across as on the histogram, so a range gate stays where it is.

Each gate of the population shown is drawn on the plot with its name and its share of the population, when the plot's two parameters are the gate's. Gates on other parameters are not drawn, but they are still in the tree.

## Beside the caption

| Button | What it opens |
|---|---|
| **Matrix…** | Every pair of parameters of the population shown, in one grid: where spillover shows. Tick the parameters in the list (at most 12; the first ticked are drawn). It starts with the stained fluorescence parameters. Each cell is coloured by the density of events on its own logarithmic ramp. |
| **Backgate…** | The population shown on each gate above it: grey is every event the gate was applied to, colour the events of the population, wherever they fall. Not available on **All events**, which has no gate above it. |
| **Model…** | A menu: **Proliferation…**, **Cell cycle…** and **Kinetics…**, fitted to the population shown in every sample selected. See [[proliferation]], [[cell-cycle]] and [[kinetics]]. |

Both **Matrix…** and **Backgate…** have their own **Figure…** button.

## Warnings above the plot

Two bars appear over the plot when the gates may no longer sit where they should.

- **These gates were drawn before …** appears after something changed the values every gate cuts: a new compensation matrix, unmixing, or batch normalisation applied to the gates. Look at the gates, move what needs moving, then press **Reviewed**.
- **… was fitted automatically and then moved by hand …** appears when a gate proposed by **Auto…** was moved by hand and the population it is in has changed since. **Refit** fits it again on the population as it is now and gives up the hand's work; **Keep** leaves it where it was moved to and stops saying so. See [[auto-gates]].

## Keys

**Esc** cancels the gate being drawn. The other shortcuts of the window are in [[keyboard-shortcuts]].
