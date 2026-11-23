"""
FlowSOM (Van Gassen et al., 2015): a self-organising map of the events,
then the map's nodes merged into metaclusters by consensus hierarchical
clustering.

Written after the R package's own steps — the online SOM with a bubble
neighbourhood shrinking from the 0.67 quantile of the grid distances to
zero, the learning rate falling from 0.05 to 0.01, ten passes; then
ConsensusClusterPlus's scheme, a hundred 90% resamples of the nodes under
average linkage. R's random numbers are not reproduced, so the clusters
match FlowJo's in substance and not event for event: see tests/test_real.py
for how closely.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit


@dataclass
class FlowSOM:
    xdim: int
    ydim: int
    #: nodes x features, in the space the events were given in
    codes: np.ndarray
    #: per event: its node, from 0
    node: np.ndarray
    #: per node: its metacluster, from 1
    node_cluster: np.ndarray
    #: (i, j) pairs of nodes: the minimal spanning tree of the codes
    mst: np.ndarray
    seed: int = 0

    @property
    def labels(self) -> np.ndarray:
        """Per event: its metacluster, from 1."""
        return self.node_cluster[self.node]

    @property
    def n_clusters(self) -> int:
        return int(self.node_cluster.max())

    def node_sizes(self) -> np.ndarray:
        return np.bincount(self.node, minlength=self.codes.shape[0])

    def grid(self) -> np.ndarray:
        """nodes x 2: where each node sits on the map."""
        return _grid(self.xdim, self.ydim)

    def map(self, X: np.ndarray) -> np.ndarray:
        """New events on this map: their metaclusters."""
        return self.node_cluster[_nearest(np.ascontiguousarray(X, dtype=np.float64),
                                          self.codes)]


def _grid(xdim: int, ydim: int) -> np.ndarray:
    gx, gy = np.meshgrid(np.arange(xdim), np.arange(ydim))
    return np.column_stack([gx.ravel(), gy.ravel()]).astype(np.float64)


@njit(cache=True)
def _train(data, codes, nhbr, alpha0, alpha1, radius0, radius1, rlen, seed):
    np.random.seed(seed)
    n, p = data.shape
    ncodes = codes.shape[0]
    niter = rlen * n
    for k in range(niter):
        i = np.random.randint(n)
        best, best_d = 0, np.inf
        for c in range(ncodes):
            d = 0.0
            for j in range(p):
                t = data[i, j] - codes[c, j]
                d += t * t
            if d < best_d:
                best, best_d = c, d
        frac = k / niter
        threshold = radius0 - (radius0 - radius1) * frac
        alpha = alpha0 - (alpha0 - alpha1) * frac
        for c in range(ncodes):
            if nhbr[c, best] > threshold:
                continue
            for j in range(p):
                codes[c, j] += alpha * (data[i, j] - codes[c, j])
    return codes


@njit(cache=True)
def _nearest(data, codes):
    n, p = data.shape
    out = np.empty(n, dtype=np.int64)
    for i in range(n):
        best, best_d = 0, np.inf
        for c in range(codes.shape[0]):
            d = 0.0
            for j in range(p):
                t = data[i, j] - codes[c, j]
                d += t * t
            if d < best_d:
                best, best_d = c, d
        out[i] = best
    return out


def consensus_clusters(codes: np.ndarray, k: int, reps: int = 100,
                       p_item: float = 0.9, seed: int = 0) -> np.ndarray:
    """Per row of `codes` a cluster from 1 to k, numbered from the
    largest group of rows to the smallest."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist, squareform

    n = codes.shape[0]
    k = max(1, min(k, n))
    if k == 1:
        return np.ones(n, dtype=np.int64)
    rng = np.random.default_rng(seed)
    together = np.zeros((n, n))
    sampled = np.zeros((n, n))
    m = max(k, int(np.floor(p_item * n)))
    for _ in range(reps):
        rows = np.sort(rng.choice(n, size=min(m, n), replace=False))
        cut = fcluster(linkage(codes[rows], "average"), k, "maxclust")
        sampled[np.ix_(rows, rows)] += 1
        together[np.ix_(rows, rows)] += cut[:, None] == cut[None, :]
    consensus = np.divide(together, sampled, out=np.zeros_like(together),
                          where=sampled > 0)
    np.fill_diagonal(consensus, 1.0)
    dist = squareform(1.0 - consensus, checks=False)
    labels = fcluster(linkage(dist, "average"), k, "maxclust")
    return labels


def flowsom(X: np.ndarray, n_clusters: int = 5, xdim: int = 10, ydim: int = 10,
            rlen: int = 10, seed: int = 0) -> FlowSOM:
    """`X` is events x features, already transformed: nothing is scaled
    here, as FlowJo's plugin does not scale either."""
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.spatial.distance import cdist, pdist, squareform

    X = np.ascontiguousarray(X, dtype=np.float64)
    n = X.shape[0]
    ncodes = xdim * ydim
    rng = np.random.default_rng(seed)
    codes = X[rng.choice(n, size=ncodes, replace=n < ncodes)].copy()
    grid = _grid(xdim, ydim)
    nhbr = cdist(grid, grid, "chebyshev")
    radius = float(np.quantile(nhbr, 0.67))
    codes = _train(X, codes, nhbr, 0.05, 0.01, radius, 0.0, int(rlen),
                   int(seed) % (2 ** 31))
    node = _nearest(X, codes)

    node_cluster = consensus_clusters(codes, n_clusters, seed=seed)
    # numbered by how many events they hold, so 1 is always the largest
    # and a rerun with another seed does not shuffle the names for nothing
    sizes = np.bincount(node_cluster[node], minlength=node_cluster.max() + 1)
    order = np.argsort(-sizes[1:], kind="stable") + 1
    rename = np.zeros(sizes.size, dtype=np.int64)
    rename[order] = np.arange(1, order.size + 1)
    node_cluster = rename[node_cluster]

    tree = minimum_spanning_tree(squareform(pdist(codes))).tocoo()
    mst = np.column_stack([tree.row, tree.col]).astype(np.int64)
    return FlowSOM(xdim, ydim, codes, node, node_cluster, mst, int(seed))


def adjusted_rand(a: np.ndarray, b: np.ndarray) -> float:
    """How alike two labellings of the same events are, whatever their
    cluster names: 1 identical, about 0 unrelated."""
    from sklearn.metrics import adjusted_rand_score
    return float(adjusted_rand_score(a, b))


# --- the tree ----------------------------------------------------------------

def tree(codes: np.ndarray, sizes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The minimal spanning tree of the nodes that hold events, laid out
    flat: (positions, nodes x 2 with NaN for an empty node; edges, pairs
    of node indices). The layout puts two nodes as far apart as they are
    steps apart along the tree — what Kamada-Kawai does in the R package —
    and is the same every time."""
    from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path
    from scipy.spatial.distance import pdist, squareform

    used = np.flatnonzero(np.asarray(sizes) > 0)
    pos = np.full((codes.shape[0], 2), np.nan)
    if used.size == 0:
        return pos, np.zeros((0, 2), dtype=np.int64)
    if used.size == 1:
        pos[used[0]] = 0.0
        return pos, np.zeros((0, 2), dtype=np.int64)
    mst = minimum_spanning_tree(squareform(pdist(codes[used]))).tocoo()
    edges = np.column_stack([used[mst.row], used[mst.col]]).astype(np.int64)
    hops = shortest_path(mst, directed=False, unweighted=True)
    pos[used] = _stress_layout(hops)
    return pos, edges


def _stress_layout(D: np.ndarray, iterations: int = 300) -> np.ndarray:
    """Points in the plane whose distances are as close to `D` as they
    can be (SMACOF), started from classical scaling so there is nothing
    random in it."""
    n = D.shape[0]
    J = np.eye(n) - 1.0 / n
    w, v = np.linalg.eigh(-0.5 * J @ (D ** 2) @ J)
    top = np.argsort(w)[::-1][:2]
    X = v[:, top] * np.sqrt(np.maximum(w[top], 1e-12))
    for _ in range(iterations):
        diff = X[:, None, :] - X[None, :, :]
        dist = np.sqrt((diff ** 2).sum(axis=2))
        np.fill_diagonal(dist, 1.0)
        B = -D / np.maximum(dist, 1e-12)
        np.fill_diagonal(B, 0.0)
        np.fill_diagonal(B, -B.sum(axis=1))
        X = B @ X / n
    return X
