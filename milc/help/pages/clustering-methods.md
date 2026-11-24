---
title: FlowSOM and PhenoGraph
---

MIL-C has two ways of sorting the concatenated events into clusters, both run from [[explore-tab|the Explore tab]] and both shown in [[clusters-tab|the Clusters tab]]. **FlowSOM** is told how many clusters to make; **PhenoGraph** finds the number itself. Both run on the markers ticked in Explore, on the scale chosen there (see [[embeddings]] for what the algorithms are fed), not on a UMAP or tSNE map. Nothing is scaled further: as in FlowJo's plugins, the markers arrive on one scale already.

Every event ends in exactly one cluster, numbered from 1, the largest first, so a rerun with another seed does not shuffle the names for nothing.

## FlowSOM

FlowSOM (Van Gassen et al., *Cytometry A* 2015) is written after the R package's own steps.

1. **The self-organising map.** A grid of nodes (10 × 10 by default) is placed in the marker space, starting at randomly chosen events. The events are presented one at a time, ten passes over all of them; each pulls its nearest node, and the nodes near it on the grid, towards itself. The neighbourhood is a bubble that shrinks from the 0.67 quantile of the grid distances to zero, and the learning rate falls from 0.05 to 0.01.
2. **Each event to its node.** Every event goes to its nearest node.
3. **Metaclusters.** The nodes are merged into the number of clusters asked for by consensus hierarchical clustering, as ConsensusClusterPlus does: a hundred resamples of 90% of the nodes, each cut by average linkage; how often two nodes fall together becomes their similarity, and that is cut once more, by average linkage, into the final clusters.

Each event's cluster is its node's metacluster.

| Setting | Default | Effect |
|---|---|---|
| **Clusters** | 5 | The number of metaclusters, 2 to 60. Ask for more than you expect: merging two in the Clusters tab is easy, and splitting one is not |
| **Grid** | 10 | The map is this many nodes on each side (100 nodes at 10). Larger grids see finer structure and take longer |
| **Seed** | 0 | The starting nodes, the order of the events and the resamples |

**The tree.** The Clusters tab draws the minimal spanning tree of the nodes that hold events, each node as large as it is populated. The layout puts two nodes as far apart as they are steps apart along the tree — what Kamada-Kawai does in the R package — and is the same every time.

R's random numbers are not reproduced, so the clusters match FlowJo's in substance and not event for event: on a real experiment, the adjusted Rand index between MIL-C's and FlowJo's clusters was 0.86 and 0.81, against 0.81 and 0.79 between two runs of MIL-C with different seeds (see [[measured-facts]]).

Only a FlowSOM cluster can become a population of the full files (**Make population** in the Clusters tab), because the trained map can place any event.

## PhenoGraph

PhenoGraph (Levine et al., *Cell* 2015) finds clusters as communities of a graph. It is written from the papers, not from the PhenoGraph package: the steps are its steps, the labels are not promised to be its labels.

1. **Neighbours.** Each event is joined to its *k* nearest neighbours by Euclidean distance: exactly (a k-d tree) up to 50,000 events or 10 markers, by nearest-neighbour descent beyond, on one thread so that the same seed gives the same neighbours.
2. **Jaccard weights.** Each join weighs the share of the two events' neighbourhoods they have in common. An event counts in its own neighbourhood, and a pair is joined when either lists the other, at the full weight.
3. **Communities.** The graph is cut into the communities that maximise modularity, by the **Leiden** algorithm (Traag et al., 2019) or the **Louvain** method (Blondel et al., 2008), once, from the seed.
4. **Small clusters.** A cluster of fewer events than the minimum is taken into the cluster it is most strongly joined to, smallest first. The original package leaves such events unassigned; here every event must be in a cluster, for the tables of frequencies.

| Setting | Default | Effect |
|---|---|---|
| **Neighbours** (*k*) | 30 | Each event is joined to this many. Fewer finds smaller clusters; more joins them |
| **Resolution** | 1.0 | Above 1, more and smaller clusters; below 1, fewer and larger |
| **Method** | Leiden | Leiden finds a slightly better cut and never leaves a cluster in disconnected pieces; Louvain is what the original PhenoGraph used |
| **Min. cluster size** | 10 | A cluster of fewer events is merged into its neighbour; 1 keeps every cluster |
| **Seed** | 0 | The neighbour search beyond 50,000 events and the order the communities are searched in |

After a run, Explore reports the number of clusters and the modularity, with how it was made, e.g. `PhenoGraph: 14 clusters, modularity 0.85 (Leiden, k 30, resolution 1, 2 of fewer than 10 events merged)`.

### Why it may give more clusters than you expect

Modularity has a scale of its own and no way of saying "one population". At resolution 1, a large population comes out in several pieces whether or not it has parts: one Gaussian of 3,500 events on 8 markers comes out in 13 clusters at k = 30, here and in the original package alike. The better the optimiser, the more surely the large population is cut. The hint under the settings says what to do:

- lower **Resolution** (0.5 is a good start), at the risk of losing a rare population in a neighbour;
- or raise **Neighbours** (100 to 150), at a higher cost in time;
- or keep the pieces and merge them in the Clusters tab.

A population much smaller than *k* is taken into its neighbour at any resolution.

## Which to use

| | FlowSOM | PhenoGraph |
|---|---|---|
| Number of clusters | You choose | Found from the data and the resolution |
| Time | Grows with events × nodes × passes | Dominated by the neighbour search, which grows with events and markers |
| Reproducibility | Same seed, same result | Same seed, same result |
| As FlowJo's plugin | Yes, in substance | No FlowJo equivalent |
| Can become a population of the files | Yes | No |

Both can be run on the same events; the Clusters tab then offers a choice between them. Neither is the truth: two runs of either, with different seeds, agree about as well as either agrees with the original package (see [[measured-facts]]). Look at the [[cluster-statistics|clusters' statistics]] across seeds before trusting a small cluster.
