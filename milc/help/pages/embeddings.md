---
title: UMAP, tSNE and PCA
---

An embedding puts every event on a plane so that events alike in many markers lie close together. [[explore-tab|The Explore tab]] offers three: UMAP and tSNE, which draw the map most people expect, and PCA, whose two axes can be read. All three run on the concatenated events (see [[concatenation]]), on the markers ticked and the scale chosen in Explore.

The defaults are the ones the FlowJo plugins ran with, so a first map is comparable to the one it replaces.

## What the algorithms are fed

| Scale | Values |
|---|---|
| Biexponential (default) | Each marker on the scale its axis has in Gates — FlowJo's biexponential unless another was chosen there — as numbers from 0 to 1023. This is what FlowJo hands its plugins. A marker shown linear in Gates is put on 0 to 1023 linearly; one shown as arcsinh in Gates is fed as its arcsinh |
| Arcsinh, cofactor 150 | asinh(x / 150) of every marker |
| Logicle, W 0.5, M 4.5 | The standard Logicle (top 262,144) of every marker, on 0 to 1023 |

Nothing else is scaled or centred, except that PCA centres the values. Fluorescence is compensated (see [[compensation]]). The scale matters: on a scale that squashes the negatives, events near zero merge; on a linear one, bright markers dominate. See [[scales-and-transforms]].

## UMAP

UMAP (McInnes et al., 2018) builds a graph of each event's nearest neighbours in the marker space and lays it out on the plane so that the graph is kept. MIL-C calls the `umap-learn` package, Euclidean distance, two dimensions, its original ("compatibility") layout.

| Setting | Default | Effect |
|---|---|---|
| **Neighbours** | 15 | How many neighbours each event is joined to. Fewer: finer local structure, more scattered islands. More: a broader picture where large populations keep their relative places. Never more than the number of events less one |
| **Min. distance** | 0.5 | How close events may be packed on the map. Lower: tight, clumped islands. Higher: evenly spread blobs |
| **Seed** | 0 | With a seed UMAP runs on one thread, so it is reproducible and slower |

## tSNE

tSNE (van der Maaten and Hinton, 2008) matches, on the plane, the probability that two events are neighbours in the marker space. MIL-C calls the `openTSNE` package, Euclidean distance, 1,000 iterations (250 of early exaggeration, then 750), on all cores.

| Setting | Default | Effect |
|---|---|---|
| **Perplexity** | 30 | Roughly the number of neighbours each event weighs. Lower: small, tight groups; higher: smoother, larger structure. Capped at a third of the number of events |
| **Seed** | 0 | The random starting layout |

Distances between islands on a tSNE map mean little; sizes of islands mean little too.

## PCA

The principal components are the directions along which the events vary most, found by singular value decomposition of the centred values. Unlike UMAP and tSNE they are the same every time, and each axis can be read: it is so much of this marker and so much of that.

- The axes are labelled with their share of the total variance, `PCA1 (41%)`.
- Under the run buttons, each component is described by its three heaviest markers with their weights (loadings): `PCA1, 41%: +0.72 CD4, −0.55 CD8, +0.30 CD25`. A marker with a large positive weight is high on the right of that axis.
- The sign of a component is arbitrary; MIL-C makes its heaviest marker positive, so a rerun points the same way.
- There is no setting and no seed. Two components are computed.

## How good is a map

Neither UMAP nor tSNE gives the same picture twice without the same seed, and never the same picture as another implementation. What can be compared is how well a map keeps the neighbours an event had in the marker space. After each UMAP and tSNE run, Explore reports:

> UMAP keeps 47% of each event's 15 nearest neighbours

That is: for each event, of its 15 nearest neighbours in the markers, the share still among its 15 nearest on the map, averaged over the events (over a random 20,000 where there are more). 1 would be a perfect map; chance is about 15 divided by the number of events. FlowJo's own maps of the same events keep about the same share (see [[measured-facts]]).

## Reading a map

- Islands are groups of events alike in the markers fed; their placement relative to each other is only loosely meaningful.
- A map built on everything is dominated by the groups with the most events. The balanced sampling in [[concatenation]] is there so that every sample contributes the same.
- Colour by a marker to see what separates the islands, and by sample or group to see whether an island is one sample's.
- Clustering is not done on the map: FlowSOM and PhenoGraph run on the markers themselves. See [[clustering-methods]].
