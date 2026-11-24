---
title: Drawing and editing gates
---

Every gate is drawn on the population shown in the [[gates-tab]], against the two parameters on the plot (or the one parameter of a histogram). It becomes a child of that population in the tree. Choose the population first, then the axes, then the tool.

## The tools

The nine buttons under **New gate on the population shown**:

| Button | What it does |
|---|---|
| **Rectangle** | Drag across the events. On a histogram it draws a range of X. |
| **Polygon** | Click each corner; double-click to close (three corners at least). |
| **Lasso** | Hold the mouse button down and draw round the events freehand; let go to close. It becomes a polygon of a few dozen corners, moved like any other. |
| **Quadrant** | Click where the two lines should cross: four populations at once. |
| **Ellipse** | Drag the box the ellipse fits in; then turn and stretch it by its handles. |
| **Auto…** | A menu of gates proposed from the events. See [[auto-gates]]. |
| **Combine…** | A population made of others of the population shown: in all of them, any of them, or none of them. |
| **Rename** | A new name for the population selected in the tree. |
| **Delete** | Deletes the population selected, and everything gated inside it. |

The drawing tools stay pressed until the gate is drawn; pressing another tool, pressing the same one again, or **Esc** cancels. Every tool except **Rectangle** needs two parameters: on a histogram the status bar says `That gate needs two parameters: choose a Y axis`.

When a shape is closed, a box asks for its name. The name offered is `Gate`, or `Gate 2`, `Gate 3`… where that is taken: two populations of one parent cannot share a name.

## What each gate keeps

- A **rectangle** is kept in the units the file stores, so it does not depend on the scale of the axes. An event on its lower edge is inside, one on its upper edge is outside. A range on a histogram is a rectangle on one parameter.
- A **polygon** (and a lasso) is tested on the scale it was drawn on: a side that is straight on a biexponential axis is the line the events are compared with, as FlowJo does. If the scale of that axis is changed later, the polygon follows: its corners stay at the same values and its sides are straight on the new scale. See [[scales-and-transforms]].
- An **ellipse** keeps its centre, half-axes and angle in the coordinates it was drawn in. An ellipse cannot be carried to another scale without becoming another shape, so after its axis is given another scale it is drawn as the curve it has become, with its name and share, but it can no longer be dragged.
- A **quadrant** cuts the plot at one point, in the file's units. An event exactly on a line goes to the high side.

## The lasso

What the hand draws has hundreds of points. It is thinned to a polygon by the Ramer–Douglas–Peucker method: corners are kept only where the drawn path would otherwise stray more than 0.4% of the plot's width or height from the polygon. Where that still leaves more than 60 corners, the tolerance is loosened by half again, until there are 60 or fewer. A mouse resting in one place adds nothing. Where the hand went back over its own track, the polygon follows it.

## Quadrants

A quadrant makes four populations, named after the first word of each axis label and the side of each line they are on: `CD4+ CD8-`, `CD4+ CD8+`, `CD4- CD8+`, `CD4- CD8-`. They are numbered internally as FlowJo numbers them: Q1 upper left, Q2 upper right, Q3 lower right, Q4 lower left. Each corner's name and share are written in its corner of the plot.

- Drag either line to move the crossing; the four populations move together.
- Each arm has a small square handle. Drag one to tilt that arm alone. From then on the four arms are half-lines of their own (FlowJo's quadrant with offset arms), and a round handle at the centre moves all four.
- An arm cannot be turned onto or past the arm next to it; if it is let go there it goes back, and the status bar says so. Every event is always in exactly one corner.

A FlowJo **curly quadrant**, whose upper and right arms bend away from the centre, is read from a workspace but not drawn here; it can be moved by its centre. See [[workspaces-and-gating-ml]].

## Moving, renaming, deleting

- Drag a gate's outline to move it, or its handles to reshape it. The counts in the tree and the percentages on the plot follow when the mouse is let go.
- With **Apply to every sample** ticked, a move, a rename or a delete is made on every sample. Unticked, only on the one shown.
- **Rename** follows the population to every sample (when ticked), and a combination made of it is renamed with it.
- **Delete** removes the population, every population under it, and every combination beside it that was made of it (it would otherwise be left empty). Deleting the population made by quality control takes quality control off instead (see [[quality-control]]).
- **Back to the workspace's gates** (under the tools) forgets everything drawn or moved and returns to the workspace's gates.

## Combinations

**Combine…** opens **Combine populations**. Tick the populations of the one shown that are to be combined, and choose under **Events**:

| Choice | Events kept |
|---|---|
| **In all of them (AND)** | inside every population ticked |
| **In any of them (OR)** | inside at least one |
| **In none of them (NOT)** | inside none |

The **Name** box offers a name made from the choice, such as `CD4+ and CD69+` or `not CD4+ or CD8+`. A combination has no shape and is not drawn on the plot. If a population it names is renamed away or deleted, the combination is empty rather than a guess.

## Drawing on several samples

Selecting several samples in **Samples** overlays them on the plot, but the gates shown and drawn are those of the first sample (the caption says so). With **Apply to every sample** ticked that makes no difference: every sample gets the gate.

To place a gate differently in each sample, untick **Apply to every sample** and draw or move it on each sample in turn, or use **Fit on each sample separately** in the [[auto-gates]] menu.
