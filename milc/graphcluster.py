"""
Clusters without being told how many: the method of PhenoGraph (Levine et
al., Cell 2015). Each event is joined to its k nearest neighbours, the
join weighs what share of their two neighbourhoods the pair has in common
(Jaccard), and the graph is cut into the communities that maximise
modularity, by the Louvain method (Blondel et al., 2008).

Written from the two papers, not from the PhenoGraph package, which was
not at hand to compare with: the steps are its steps, the labels of a
given file are not promised to be its labels. Where this is known to
differ: an event counts in its own neighbourhood, so that no join ever
weighs nothing and no event is left alone; a pair is joined when either
lists the other; Louvain is run once, from the seed, where the package
keeps the best of several runs; and no cluster is discarded for being
small. What is checked is in tests/test_graphcluster.py: mixtures whose
truth is known, a graph whose best partition is known, and the modularity
against its definition written out in full.

Modularity has a scale of its own: among many events a large uniform
population scores higher cut in two than whole (measured: 8,000 events
of one Gaussian among 20,000, at resolution 1), and a population much
smaller than k is taken into its neighbour. `resolution` below 1 joins,
above 1 divides.

Louvain can leave a community in two pieces that only touch through
events that have since moved away; Leiden's refinement mends that and is
not implemented here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit


@dataclass
class Clustering:
    #: per event: its cluster, from 1, the largest first
    labels: np.ndarray
    #: of the partition on the graph, at the resolution asked for
    modularity: float
    k: int
    resolution: float = 1.0
    seed: int = 0

    @property
    def n_clusters(self) -> int:
        return int(self.labels.max())


# --- the graph ---------------------------------------------------------------

#: below either, the exact search is as fast as the approximate one
EXACT_EVENTS, EXACT_FEATURES = 50_000, 10


def neighbours(X: np.ndarray, k: int = 30, seed: int = 0) -> np.ndarray:
    """events x k: each event's nearest other events, by Euclidean
    distance. Exact (a k-d tree) for small or low-dimensional data and
    by nearest-neighbour descent beyond. Measured at 200,000 events and
    k = 30 on ten cores: the tree takes 2 s on 5 markers, 12 s on 10 and
    30 s on 15; the descent 11 s on 15, finding 99.8% of the true
    neighbours, after some 25 s compiling the first time in a session.
    It is kept to one thread, which is what makes the same seed give the
    same neighbours."""
    n, p = X.shape
    k = max(1, min(int(k), n - 1))
    if n <= EXACT_EVENTS or p <= EXACT_FEATURES:
        from sklearn.neighbors import NearestNeighbors
        found = NearestNeighbors(n_neighbors=k + 1, algorithm="kd_tree", n_jobs=-1) \
            .fit(X).kneighbors(X, return_distance=False)
    else:
        from pynndescent import NNDescent
        found = NNDescent(X, n_neighbors=k + 1, random_state=int(seed),
                          n_jobs=1).neighbor_graph[0]
    own = found == np.arange(n)[:, None]
    # identical events can crowd one out of its own list
    own[~own.any(axis=1), -1] = True
    return found[~own].reshape(n, k).astype(np.int64)


@njit(cache=True)
def _jaccard(hood):
    """`hood` is events x m, each row sorted: the event and its
    neighbours. Per event and member of its row, the share of the two
    rows that is common to both."""
    n, m = hood.shape
    out = np.empty((n, m))
    for i in range(n):
        for a in range(m):
            j = hood[i, a]
            common, p, q = 0, 0, 0
            while p < m and q < m:
                if hood[i, p] == hood[j, q]:
                    common += 1
                    p += 1
                    q += 1
                elif hood[i, p] < hood[j, q]:
                    p += 1
                else:
                    q += 1
            out[i, a] = common / (2 * m - common)
    return out


def jaccard_graph(nearest: np.ndarray):
    """The symmetric sparse matrix of the joins: for a pair where either
    lists the other, the Jaccard similarity of their neighbourhoods,
    each event counted in its own."""
    from scipy import sparse
    n = nearest.shape[0]
    hood = np.sort(np.column_stack([np.arange(n), nearest]), axis=1)
    weight = _jaccard(hood)
    rows = np.repeat(np.arange(n), hood.shape[1])
    cols = hood.ravel()
    apart = rows != cols
    graph = sparse.csr_matrix((weight.ravel()[apart], (rows[apart], cols[apart])),
                              shape=(n, n))
    return graph.maximum(graph.T).tocsr()


# --- the communities ---------------------------------------------------------

@njit(cache=True)
def _move(indptr, indices, weights, degree, total, resolution, seed):
    """Louvain's first phase: every node starts alone and, in a random
    order and until none moves, goes to the neighbouring community that
    raises modularity the most. Returns the community of each node."""
    n = degree.size
    community = np.arange(n)
    inside = degree.copy()          # per community: the degrees of its nodes
    link = np.zeros(n)              # from the node in hand to each community
    seen = np.empty(n, dtype=np.int64)
    np.random.seed(seed)
    order = np.random.permutation(n)
    tolerance = 1e-12 * total       # a gain smaller than rounding is no gain
    moved = True
    while moved:
        moved = False
        for i in order:
            own = community[i]
            n_seen = 0
            for e in range(indptr[i], indptr[i + 1]):
                j = indices[e]
                if j == i:
                    continue
                c = community[j]
                if link[c] == 0.0:
                    seen[n_seen] = c
                    n_seen += 1
                link[c] += weights[e]
            inside[own] -= degree[i]
            cost = resolution * degree[i] / total
            best, best_gain = own, link[own] - cost * inside[own]
            for s in range(n_seen):
                c = seen[s]
                gain = link[c] - cost * inside[c]
                if gain > best_gain + tolerance:
                    best, best_gain = c, gain
                link[c] = 0.0
            inside[best] += degree[i]
            if best != own:
                community[i] = best
                moved = True
    return community


def louvain(graph, resolution: float = 1.0, seed: int = 0) -> tuple[np.ndarray, float]:
    """(community of each node, from 0; modularity). `graph` is a
    symmetric sparse matrix of positive weights; a weight on the
    diagonal counts twice the loop it stands for, which is what merging
    nodes produces. A higher `resolution` gives more, smaller
    communities."""
    from scipy import sparse
    A = sparse.csr_matrix(graph, dtype=np.float64)
    A.sum_duplicates()
    member = np.arange(A.shape[0])
    level = 0
    while A.nnz:
        degree = np.asarray(A.sum(axis=1)).ravel()
        community = _move(A.indptr, A.indices, A.data, degree, degree.sum(),
                          float(resolution), (int(seed) + level) % (2 ** 31))
        kept, community = np.unique(community, return_inverse=True)
        if kept.size == A.shape[0]:
            break
        member = community[member]
        # the second phase: each community becomes one node
        n = A.shape[0]
        fold = sparse.csr_matrix((np.ones(n), (np.arange(n), community)),
                                 shape=(n, kept.size))
        A = (fold.T @ A @ fold).tocsr()
        level += 1
    return member, modularity(graph, member, resolution)


def modularity(graph, labels: np.ndarray, resolution: float = 1.0) -> float:
    """Newman's modularity of a partition of a symmetric weighted graph:
    the share of the weight that falls inside the communities, less what
    would fall there if the same degrees were joined at random."""
    from scipy import sparse
    A = sparse.csr_matrix(graph, dtype=np.float64)
    total = A.sum()
    if total == 0:
        return 0.0
    _, c = np.unique(labels, return_inverse=True)
    n = A.shape[0]
    fold = sparse.csr_matrix((np.ones(n), (np.arange(n), c)), shape=(n, c.max() + 1))
    within = (fold.T @ A @ fold).diagonal()
    degree = np.bincount(c, weights=np.asarray(A.sum(axis=1)).ravel())
    return float(np.sum(within / total - resolution * (degree / total) ** 2))


def by_size(labels: np.ndarray) -> np.ndarray:
    """The same partition named from 1, the largest cluster first, so a
    rerun with another seed does not shuffle the names for nothing."""
    _, c = np.unique(labels, return_inverse=True)
    order = np.argsort(-np.bincount(c), kind="stable")
    rename = np.empty(order.size, dtype=np.int64)
    rename[order] = np.arange(1, order.size + 1)
    return rename[c]


def phenograph(X: np.ndarray, k: int = 30, resolution: float = 1.0,
               seed: int = 0) -> Clustering:
    """`X` is events x features, already transformed; nothing is scaled
    here. `k` is how many neighbours each event is joined to: fewer
    finds smaller clusters."""
    X = np.ascontiguousarray(X, dtype=np.float64)
    nearest = neighbours(X, k, seed)
    graph = jaccard_graph(nearest)
    member, q = louvain(graph, resolution, seed)
    return Clustering(by_size(member), q, nearest.shape[1], float(resolution), int(seed))
