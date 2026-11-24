---
title: Gates proposed from the events
---

**Auto…** in the [[gates-tab]] proposes a gate from the events of the population shown: the valley between two peaks of a histogram, the outline or the ellipse of one population on two parameters, or the band of singlets along a diagonal. The proposal is drawn dashed on the plot first. Accepted, it is an ordinary range, polygon or ellipse: it can be moved, renamed and deleted like one drawn by hand (see [[drawing-gates]]). When the events do not show what a method looks for, nothing is proposed and a message says why.

Everything is measured on the plot's own coordinates, after the scale of each axis (see [[scales-and-transforms]]): a population is a peak where it is looked at, not in the units the file stores.

## The menu

| Item | What is proposed |
|---|---|
| **Valley on X: the positive side** | A range of X from the valley up. |
| **Valley on X: the negative side** | A range of X from the valley down. |
| **Population at a click: contour** | Click on a population: the contour holding 90% of its events, as a polygon. |
| **Population at a click: ellipse** | Click on a population: the ellipse of a Gaussian fitted to its core. |
| **Largest population: contour** | The same contour, for the population with the most events. |
| **Largest population: ellipse** | The same ellipse, for the population with the most events. |
| **Singlets: band along the diagonal** | A band three standard deviations either side of a straight line of Y against X. |
| **Singlets: a tighter band** | The same, with the spread measured on the cleaner side of the line only. |
| **Fit on each sample separately** | A setting, ticked or not, for the next gate accepted: see below. |

After a proposal a box shows what was found (for example `valley at 18% of the lower peak; 32.4% of the events above it`) and asks for a name. The name offered is the first word of the X label with `+` or `-` for a valley, `Population` for a contour or an ellipse, and `Singlets` for a band. **OK** accepts; **Cancel** throws the proposal away.

The valley works on X alone and can be used on a histogram. The others need two parameters.

## How each is found

### Valley between two peaks

1. The histogram of X is taken over where the events are: from their 0.1% to their 99.9% point, widened by 5%, in 256 bins. Events piled at the very ends of the axis are left out.
2. It is smoothed by a Gaussian whose width depends on the number of events: 3 bins for 5,000 events or more, up to 10 bins for few events.
3. A peak counts only if it stands out from its surroundings by five times what counting noise alone does to a smoothed bin (and by at least 0.1% of the tallest). A small population is a low, wide peak and still counts.
4. The two peaks that stand out most are the two populations; the threshold is the lowest point between them.

Nothing is proposed when there are fewer than 100 events, when there is only one peak, when the valley is still higher than 70% of the lower peak, or when one side would hold less than 1% of the events. A threshold is then a matter of judgement, or of a control.

### Populations on two parameters

The plot is cut into populations first. The events are counted on a 128 × 128 grid (the cells at the edges, where events pile up, are emptied) and smoothed. Every cell is joined to the peak it reaches by walking uphill. A peak whose pass to a neighbour is higher than 85% of its own height, or whose ground holds less than 1% of the events, is a ripple of the neighbour and is merged into it. Each remaining peak and its ground is one population.

- **Contour**: the cells of the chosen population that hold 90% of its events, densest first, and the line round them. Because the plot was cut into populations first, the line stops at the valley to a neighbour instead of running round both. It is thinned to at most 32 corners, so it can be moved by hand. The message says how many events are inside and how many other populations are on the plot.
- **Ellipse**: a two-dimensional Gaussian fitted to the population's core. The fit starts from the densest half of the population, then repeats: keep the events within the distance that holds 90% of a Gaussian, take their mean and covariance (enlarged for the tails cut off), until nothing moves. Stragglers and other populations are outside the cut and do not pull on it. The ellipse drawn holds 95% of that Gaussian. This is what `flowCore::norm2Filter` or a one-cluster `flowClust` are used for in R.

A click that falls between populations chooses the nearest peak.

### Singlets band

For singlets on a pulse's area against its height: a doublet has more area for its height and falls off the line. The line is fitted by Tukey's biweight regression, which gives events far from it no say. The spread about it is the median absolute distance of the events the fit listens to, scaled to a standard deviation. The band reaches three of those either side, from the lowest to the highest X of the events it holds. This is the method of `flowStats::singletGate` (a robust linear model and a prediction band).

**A tighter band** is for a line with a fan of doublets on one side, which widens the spread: the line is moved onto the ridge, where the events are densest across it, and the spread is measured on its cleaner side only.

Nothing is proposed when the events do not follow a line: when the spread about the best line is more than half of the spread of Y itself.

## One gate for all, or one for each sample

With **Apply to every sample** ticked and **Fit on each sample separately** unticked, the gate is fitted on the sample shown and the same gate is given to every sample. It is remembered that it was fitted there: when the population it is in changes (a gate above it moved, a new compensation, a new scale) it is fitted again on that sample, and given to all again.

With both ticked, each sample gets the gate fitted on its own events, so it sits where that sample's population is. Each sample then keeps its own: moving it on one sample moves it there only. Where nothing can be found in a sample, that sample keeps the gate of the sample shown, and the status bar and the population's tooltip say so (`not fitted to this sample: …`).

A fitted gate that was moved by hand is not fitted again by itself when its population changes, so the hand's work is not lost unasked. The bar above the plot says which gates are in that state; **Refit** fits them again, **Keep** keeps them where they are.

## Limits

These are proposals, not a substitute for looking. A marker with no valley gets no valley gate; a population that sits on the shoulder of another may be merged into it; a contour is a contour of density, not of biology. Look at each proposal before accepting it, and at the gates of every sample after **Fit on each sample separately**. The same robust ellipse is the scatter gate of [[autocompensation]].
