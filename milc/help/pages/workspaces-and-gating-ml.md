---
title: FlowJo workspaces and Gating-ML
---

MIL-C reads a FlowJo 10 workspace (`.wsp`) for what the analysis needs from it — which files, the compensation matrix each was given, and each sample's gating tree with the counts FlowJo reported — and writes the gates back out for another program, as a FlowJo workspace or as Gating-ML 2.0, the ISAC standard. Gating-ML is written, not read.

## Opening a workspace

A workspace in the folder is picked up when the folder is opened (the first `.wsp`, by name, if there are several). Another can be chosen with **Workspace…** in the [[samples-tab]] or **File ▸ Open workspace…**; a folder must be open first. The line under the folder's name says which workspace the gates and compensation come from.

The workspace's samples are matched to the files of the folder by file name: where FlowJo found the file (often a path on another computer) does not matter. A file the workspace lists gets its gates and its matrix. A file it does not list gets no gates — unless every sample of the workspace has the same gates and the same matrix, in which case **Gate the files the workspace does not list** appears in the Samples tab. Ticked, those files are given the shared gates and matrix, as FlowJo does when a file is dropped into the group. Nobody has looked at those gates on those files: check them.

Gates drawn or moved in MIL-C are kept with the project (see [[projects-and-saving]]) and take the place of the workspace's from then on. **Back to the workspace's gates** in the [[gates-tab]] forgets them.

## What is read, and how well

| Gate in FlowJo | Read as | Checked |
|---|---|---|
| Rectangle, range | rectangle | against FlowJo's counts |
| Polygon | polygon, tested on the axis' biexponential or Logicle scale, as FlowJo does | against FlowJo's counts |
| Quadrant | four quadrant populations | by the rules of other readers |
| And, Or, Not | combination | by the rules of other readers |
| Ellipse | ellipse, converted from FlowJo's 256-channel grid | by the rules of other readers |
| Curly quadrant | curly quadrant | against CytoML (below) |

Rectangles and polygons are read by rules checked against FlowJo's own counts on a real experiment, 72 populations. FlowJo counts on a grid by a rule not reproduced here, so exact equality is not expected: a count is flagged as different when it is more than 1.2% (or 3 events) from FlowJo's, a margin set from what was measured on those 72. The other kinds are read the way two other readers of the format do it, FlowKit (Python) and CytoML with cytolib (R), where the two agree. No workspace at hand had them, so:

- a population read that way is shown in *italics* in the Gates tab, and its tooltip gives FlowJo's count beside the count here: comparing the two is the check;
- what could not be read with confidence is not read. The status bar lists how many populations were left out and of which kinds; everything under such a population goes with it;
- plug-in nodes (UMAP, FlowSOM, cluster explorers) are skipped: their results are computed again here, not imported.

Gates on Time are converted from FlowJo's seconds to the file's units with the file's `$TIMESTEP`.

### The curly quadrant

FlowJo's curly quadrant has its upper and right arms bend away from the centre, following the spread of a compensated negative population. The workspace holds only the centre and which corner, not the curve. MIL-C draws the curve as CytoML does: on 255 channels of each axis, the right arm rises by 0.001·d², d the channels from the centre, and the upper arm leans right by as much; the other two arms are straight. Against CytoML on 40,000 simulated events the corners hold the same events but for up to 65, which comes from cytolib's own biexponential table; with that table in place of this one, at most 1 event differs. **The curve has not been compared with FlowJo itself**: no workspace with curly quadrants and its files could be obtained. Another reader fitted a different curve to FlowJo's counts, and the two can differ by a quarter of a small corner, so compare FlowJo's count with the count here (`docs/comparacao.md`, section 7).

## Exporting the gates

**File ▸ Export gates…** opens **Export gates**:

| Control | What it does |
|---|---|
| **Format** | **FlowJo workspace (.wsp)** or **Gating-ML 2.0 (.xml)**. The extension of the file follows. |
| **Only the samples ticked in Samples** | Unticked, every sample with gates. |
| **Write each population's count** | Opens every file to write the count of every population here, so the other program has a count to compare with its own. |
| **Into** | The file; **Choose…** to browse. |

When it is done, **Gates exported** says what was written and lists what was not, and why.

### As a FlowJo workspace

One workspace with every sample: its file, its compensation matrix, the scale of each parameter and its gates, claiming FlowJo 10.8.1.

- A rectangle is written as it is; a quadrant as FlowJo writes it, one rectangle open on two sides for each corner.
- A polygon is written in the file's units on the scale it was drawn on; drawn on a scale FlowJo is not known to test on (arcsinh, logarithmic), corners are added along its sides until they are the same line on any scale.
- An ellipse is written on FlowJo's 256 channels, by its foci and four points of its edge, on linear and biexponential axes; on others it becomes a polygon.
- A curly quadrant is written as FlowJo writes it, centre and corner; a quadrant with tilted arms becomes the polygon of each corner.
- A combination becomes an AndNode, OrNode or NotNode; Time is written in seconds.

It has been read back by MIL-C and by FlowKit, event for event. It was also opened in FlowJo by hand before the 1.0 release. Compare the counts there when a workspace matters: a gate FlowJo has no word for is left out (see the report shown after export).

### As Gating-ML 2.0

Gating-ML has no samples in it, so there is one file where every sample has the same gates and the same matrix, and otherwise one file per sample, named after it. The matrix is a `spectrumMatrix`; rectangles have no transformation; polygons are written on the scale they were drawn on (`logicle`, `fasinh`, `flog`, or none for linear). FlowJo's biexponential has no Gating-ML form: it is written as the nearest Logicle, with corners added where a side would bend, so the polygon holds the same events. An ellipse is an `EllipsoidGate`, except on a biexponential axis, where it becomes a polygon. Each gate carries its name, path and count in a `custom_info`. The file is valid against the standard's schema, and FlowKit reads it and gates the same events. flowCore and other programs read it too.

### What cannot be exported

A cluster turned into a population (see [[clusters-tab]]) is the events nearest some nodes of a map that no other program has: it is left out with everything under it. The population made by [[quality-control]] is written as ranges of Time and their union. A population that is none of several is written as the negation of their union.

## On the command line

`milc gates DATA --out FILE` does the same without the window: `FILE.wsp` for a workspace, anything else for Gating-ML. `--wsp` takes the gates from a given workspace, `--format` overrides the extension, `--no-counts` skips opening the files. See [[command-line]].
