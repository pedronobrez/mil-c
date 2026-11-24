---
title: Scales and transforms
---

A scale is how a measured value becomes a position on an axis. Fluorescence that has been compensated has negative values and spans five decades, so a linear axis crushes every dim population against zero and a logarithmic one cannot show the negatives at all. The scales below are the usual answers. In MIL-C the scale belongs to a parameter, not to a plot: every plot of that parameter, in every tab, uses it, and so do the algorithms of [[embeddings]] and [[clustering-methods]], which are fed each marker on its own scale.

## Which scale a parameter starts on

- Fluorescence (every channel the compensation matrix covers, or every unmixed channel) starts on FlowJo's biexponential with its usual settings.
- Scatter, time and every other parameter start linear, as stored.

## Changing it

Press **Scale…** beside **X** or **Y** in the [[gates-tab]]. The dialog **Scale of** *parameter* has:

| Control | What it does |
|---|---|
| **Scale** | **Biexponential**, **Logicle**, **Arcsinh**, **Logarithmic** or **Linear**. The settings below change with it. |
| **Every fluorescence parameter, not just this one** | Gives the same scale to every fluorescence parameter at once. |
| **The usual one** | Goes back to the starting scale above, and closes. |
| **OK** / **Cancel** | Apply the scale, or leave it as it was. |

The scales are saved with the project (see [[projects-and-saving]]).

## The scales and their settings

### Biexponential

FlowJo's biexponential, computed the way cytolib translated FlowJo's routine: a table of 4,096 points read by interpolation, from 0 at the bottom of the axis to 1 at the top, with the top at 262,144.

| Setting | Range | Default | Meaning |
|---|---|---|---|
| **Linear width** | 1 to 5,000 | 100 | How wide the nearly linear part around zero is. Larger pulls a smeared negative population together. |
| **Negative decades** | 0 to 3 | 0 | More room below zero, in decades. |
| **Positive decades** | 2 to 7 | 4.4185 | The decades of the positive side, FlowJo's default. |

Its shape (monotonic, zero where expected, the full range at the top) is tested. It has **not** been compared with values taken out of FlowJo itself.

### Logicle

The published Logicle (Parks, Roederer and Moore, 2006), in the parameters of Moore and Parks (2012). Checked against FlowUtils' C translation of the authors' own code.

| Setting | Default | Meaning |
|---|---|---|
| **Linear width (W), decades** | 0.5 | The decades given to the nearly linear part around zero. Between 0 and M / 2. |
| **Decades (M)** | 4.5 | The decades the whole scale spans. Between 1 and 8. |
| **Negative decades (A)** | 0 | More room below zero. Between −W and M − 2W. |
| **Top of scale (T)** | 262,144 | The value at the top of the axis; 262,144 for an 18-bit cytometer. |

The dialog keeps W and A inside their allowed ranges as M and W change.

### Arcsinh

`asinh(x / cofactor)`. **Cofactor** from 1 to 10,000, default 150, the usual value for conventional fluorescence; 5 is usual for mass cytometry. The axis runs from the value of −10 × cofactor, to leave room for a compensated negative population, up to the value of 262,144.

### Logarithmic

A plain logarithm from **Lowest value** (0.01 to 10,000, default 10) to 262,144. Zero and every negative value are drawn at the bottom edge — which is the reason biexponential scales exist.

### Linear

The values as stored, from zero to the parameter's range in the file.

On every scale but linear the axis is marked at 0, 10², 10³… where there is room for the label.

## What changes when a scale changes

- **Polygons**, lassos and polygons proposed by **Auto…** follow the axis: their corners keep their values and their sides become straight on the new scale, which is where they are then tested. Counts can move a little as a result.
- **Quadrant arms** that were tilted follow the same way; a curly quadrant bends in the channels of the new scale.
- **Rectangles**, ranges and plain quadrants are kept in the file's units and do not move.
- **Ellipses** keep the scale they were drawn on. They are drawn as the curve they have become on the new axis and can no longer be dragged; go back to the old scale to move one.
- Gates fitted by **Auto…** are fitted again on the new values (see [[auto-gates]]).
- Every count elsewhere is computed again, and the other tabs redraw their plots of that parameter.

The note at the bottom of the dialog says the same: polygons drawn on this parameter are tested on the new scale, and the algorithms in Explore are fed each marker on its own scale.

## Which to use

For most conventional and spectral panels the default biexponential is the right start. Widen **Linear width** (or W) when a negative population is smeared into a broad band, and narrow it when two dim populations run together. Arcsinh with a cofactor of 150 is what many R pipelines use, and gives similar pictures. Use the logarithmic scale only for parameters that are never negative. Scatter is best left linear. Compensation itself is described in [[compensation]].
