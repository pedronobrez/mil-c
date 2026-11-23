"""
Batch normalisation in the manner of CytoNorm (Van Gassen et al.,
Cytometry A 2020).

Samples acquired on different days differ for reasons that are not
biology: a laser, a lot of antibody, a voltage. If every batch carries
reference samples that should look alike — aliquots of one control — then
whatever separates them is technical, and can be learnt and taken out:
per batch and per marker, a monotone curve from that batch's quantiles to
the quantiles of a common goal, applied afterwards to every sample of the
batch. A shift need not be the same for every cell type, so the curves
are fitted within FlowSOM metaclusters of the reference events and each
event is corrected with the curve of the cluster it falls in.

Written after the R package's steps — 99 quantiles from 0.01 to 0.99, the
goal the mean of the batches' quantiles or one batch's own, a
Fritsch-Carlson monotone spline, the identity where a batch has under 50
events in a cluster — and NOT checked against its output. It knowingly
differs in three places: beyond the outermost quantiles the correction
continues as a shift (R continues the spline's end slope); the map is
trained on the same number of events from every batch; and the map is
this project's FlowSOM, 10 x 10 by default where CytoNorm's is 15 x 15.

Everything here is in transformed units: the caller applies the display
transform before and inverts it after.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator

from .flowsom import FlowSOM, flowsom

#: what the map is trained on at most, all batches together (CytoNorm's nCells)
_SOM_EVENTS = 1_000_000


@dataclass
class Model:
    batches: list[str]
    n_markers: int
    #: what the map gave, which can be fewer than was asked for
    n_clusters: int
    #: the probabilities the quantiles are taken at
    probs: np.ndarray
    #: batches x clusters x quantiles x markers; NaN where there were too
    #: few events to estimate them
    quantiles: np.ndarray
    #: "mean" or the name of a batch
    goal: str
    #: clusters x quantiles x markers; NaN where no batch could give them
    goal_quantiles: np.ndarray
    #: None when n_clusters is 1
    som: FlowSOM | None
    #: (batch, cluster from 1, reference events there): the pairs left as
    #: they were
    fallbacks: list[tuple[str, int, int]]


def _probs(n_quantiles: int) -> np.ndarray:
    # 0 and 1 are left out: the minimum and the maximum of a sample are
    # its outliers, and a curve pinned on them follows the outliers
    return np.arange(1, n_quantiles + 1) / (n_quantiles + 1)


def _goal(quantiles: np.ndarray, batches: list[str], goal: str) -> np.ndarray:
    """`quantiles` is batches x ... : the goal's, without the first axis."""
    if goal == "mean":
        some = ~np.isnan(quantiles).all(axis=0)
        return np.where(some, np.nansum(quantiles, axis=0)
                        / np.maximum((~np.isnan(quantiles)).sum(axis=0), 1), np.nan)
    if goal not in batches:
        raise KeyError(f"goal {goal!r} is neither 'mean' nor a batch: {batches}")
    return quantiles[batches.index(goal)]


def fit(reference: dict[str, np.ndarray], n_clusters: int = 1,
        n_quantiles: int = 99, goal: str = "mean", seed: int = 0,
        grid: int = 10, min_events: int = 50) -> Model:
    """`reference` gives, per batch, the events x markers of its reference
    samples pooled, already transformed.

    With `n_clusters` 1 nothing is clustered and each marker is aligned as
    a whole. Above 1, a `grid` x `grid` FlowSOM is trained on the same
    number of events from every batch and the curves are fitted within
    its metaclusters.

    `goal` is "mean", the mean of the batches' quantiles, or the name of a
    batch every other is brought onto; that batch is then left unchanged.

    A batch with fewer than `min_events` reference events in a cluster
    cannot give `n_quantiles` quantiles worth the name: it is left out of
    the mean, its events in that cluster are not corrected, and the pair
    is listed in `Model.fallbacks`. When the goal is a batch and it is
    that batch which lacks the events, no batch is corrected in that
    cluster and all of them are listed.
    """
    batches = list(reference)
    data = [np.asarray(reference[b], dtype=np.float64) for b in batches]
    n_markers = data[0].shape[1]
    for b, d in zip(batches, data):
        if d.ndim != 2 or d.shape[1] != n_markers:
            raise ValueError(f"batch {b!r}: expected events x {n_markers} markers, "
                             f"got {d.shape}")
        if len(d) == 0:
            raise ValueError(f"batch {b!r} has no reference events")
    if goal != "mean" and goal not in batches:
        raise KeyError(f"goal {goal!r} is neither 'mean' nor a batch: {batches}")

    som = None
    labels = [np.ones(len(d), dtype=np.int64) for d in data]
    if n_clusters > 1:
        # the same from every batch, so the map is not the map of
        # whichever batch had the most reference events
        rng = np.random.default_rng(seed)
        each = min(min(len(d) for d in data), _SOM_EVENTS // len(data))
        pooled = np.vstack([d[rng.choice(len(d), size=each, replace=False)]
                            for d in data])
        som = flowsom(pooled, n_clusters, grid, grid, seed=seed)
        labels = [som.map(d) for d in data]
    n_clusters = 1 if som is None else som.n_clusters

    probs = _probs(n_quantiles)
    quantiles = np.full((len(batches), n_clusters, n_quantiles, n_markers), np.nan)
    counts = np.zeros((len(batches), n_clusters), dtype=np.int64)
    for i, (d, lab) in enumerate(zip(data, labels)):
        for c in range(n_clusters):
            rows = d[lab == c + 1]
            counts[i, c] = len(rows)
            if len(rows) >= min_events:
                quantiles[i, c] = np.quantile(rows, probs, axis=0)
    goal_quantiles = _goal(quantiles, batches, goal)

    identity = (counts < min_events) | np.isnan(goal_quantiles).all(axis=(1, 2))
    fallbacks = [(batches[i], int(c) + 1, int(counts[i, c]))
                 for i, c in zip(*np.nonzero(identity))]
    return Model(batches, n_markers, n_clusters, probs, quantiles, goal,
                 goal_quantiles, som, fallbacks)


def _curve(x: np.ndarray, source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """`x` carried by the monotone curve through (source, target), both
    sorted."""
    # values piled on one number (a floor, a discrete channel) repeat a
    # quantile; the spline wants each knot once, with the mean of its goals
    knots, where = np.unique(source, return_inverse=True)
    goals = np.bincount(where, weights=target) / np.bincount(where)
    # Outside what the reference covered nothing is known about the batch
    # effect, so the correction continues as the shift it was at the last
    # knot — slope 1 — and not along the spline's end slope, which a
    # sparse tail can make steep or flat and would then be applied to the
    # brightest events of every sample.
    out = np.where(x < knots[0], x + (goals[0] - knots[0]),
                   x + (goals[-1] - knots[-1]))
    inside = (x >= knots[0]) & (x <= knots[-1])
    if knots.size > 1 and inside.any():
        out[inside] = PchipInterpolator(knots, goals)(x[inside])
    return out


def apply(model: Model, X: np.ndarray, batch: str) -> np.ndarray:
    """One sample of `batch`, events x markers and already transformed,
    brought onto the goal: the same shape back, and the same dtype when
    that is a floating one. Each event is corrected within the cluster
    the model's map puts it in; within a cluster, each marker's
    correction never reverses the order of two events."""
    if batch not in model.batches:
        raise KeyError(f"no batch {batch!r} in the model: {model.batches}")
    X = np.asarray(X)
    if X.ndim != 2 or X.shape[1] != model.n_markers:
        raise ValueError(f"expected events x {model.n_markers} markers, got {X.shape}")
    out = X.astype(np.float64)
    labels = (np.ones(len(out), dtype=np.int64) if model.som is None
              else model.som.map(out))
    b = model.batches.index(batch)
    for c in range(model.n_clusters):
        source, target = model.quantiles[b, c], model.goal_quantiles[c]
        rows = labels == c + 1
        if np.isnan(source).any() or np.isnan(target).any() or not rows.any():
            continue
        for j in range(model.n_markers):
            out[rows, j] = _curve(out[rows, j], source[:, j], target[:, j])
    return out.astype(X.dtype if np.issubdtype(X.dtype, np.floating) else np.float64)


def distance(samples: dict[str, np.ndarray], n_quantiles: int = 99,
             goal: str = "mean") -> np.ndarray:
    """Per marker, how far apart the batches are: the mean, over batches
    and over `n_quantiles` quantiles from 1/(n+1) to n/(n+1), of the
    absolute difference between a batch's quantile and the goal's — the
    mean of the batches' quantiles, or one batch's.

    It is in the units the events are given in (transformed ones), and is
    close to the mean earth mover's distance from each batch to the goal,
    the tails beyond the outermost quantiles left out. Each marker is
    taken whole, not within clusters: zero means the batches' marginal
    distributions coincide. Samples that truly differ in composition
    keep a distance that no normalisation should remove."""
    batches = list(samples)
    q = np.stack([np.quantile(np.asarray(samples[b], dtype=np.float64),
                              _probs(n_quantiles), axis=0) for b in batches])
    return np.abs(q - _goal(q, batches, goal)).mean(axis=(0, 1))
