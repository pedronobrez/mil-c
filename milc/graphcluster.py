"""
Clusters without being told how many: the method of PhenoGraph (Levine et
al., Cell 2015). Each event is joined to its k nearest neighbours, the
join weighs what share of their two neighbourhoods the pair has in common
(Jaccard), and the graph is cut into the communities that maximise
modularity: by the Leiden algorithm (Traag, Waltman & van Eck, Sci. Rep.
2019) or, as the package did first, by the Louvain method (Blondel et
al., 2008).

Written from the papers, not from the PhenoGraph package: the steps are
its steps, the labels of a given file are not promised to be its labels.
Where this is known to differ: an event counts in its own neighbourhood,
so that no join ever weighs nothing and no event is left alone; a pair is
joined when either lists the other, at the full weight, where the
package halves the weight of a join only one of the two lists; the
communities are found once, from the seed, where the package keeps the
best of several runs of Louvain; and a cluster below the minimum size is
taken into the cluster it is most joined to, where the package leaves
its events unassigned (-1) — here every event has to be in a cluster,
for the tables of frequencies. What is checked is in
tests/test_graphcluster.py: mixtures whose truth is known, a graph whose
best partition is known, the modularity against its definition written
out in full, and that no community of Leiden's is in pieces.

Louvain can leave a community in two pieces that only touch through
events that have since moved away. Leiden refines each community into
connected sub-communities before merging, which mends that, and finds a
slightly higher modularity: on the 3,488 real events this was first run
on (k = 30, three seeds) 0.8068 against Louvain's 0.7989, and 0.8050
from the leidenalg package on the same graph.

Modularity has a scale of its own, and no way of saying "one
population": one Gaussian of 3,500 events on 8 markers comes out in 13
clusters at k = 30 and resolution 1, on this graph and on the package's
alike. A population is cut when its share of the events and its size
against k are large enough, so the same mixture found whole at 3,500
events is in pieces at 10,000, and the better the optimiser the more
surely. The real events above give 15 clusters, here and from the
package (1.5.7, by Louvain or by Leiden) alike, 93% of the events in a
cluster that lies inside one of the 5 that k = 150 gives; the same 5
populations redrawn as plain Gaussians give 7. `resolution` below 1 joins, above 1 divides: 0.5 kept
the large populations whole in every mixture tried from 3,500 to 30,000
events, but with fewer events, or populations that overlap, it took a
rare one (1%) into its neighbour, which 1 had kept; that is why 1 is
still what is offered first. A larger k joins too, at a higher cost
(k = 100 to 150 on those 3,488). A population much smaller than k is
taken into its neighbour at any resolution.
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
    method: str = "leiden"
    min_cluster_size: int = 1
    #: clusters below the minimum that were taken into another
    absorbed: int = 0

    @property
    def n_clusters(self) -> int:
        return int(self.labels.max())

    def settings(self) -> dict:
        """What a dataset keeps of the run, beside the labels."""
        return {"k": self.k, "resolution": self.resolution, "method": self.method,
                "min_cluster_size": self.min_cluster_size, "absorbed": self.absorbed,
                "n_clusters": self.n_clusters, "modularity": self.modularity}


def describe(run: dict) -> str:
    """How a kept run was made, in a few words. Runs kept before Leiden
    was written name no method: they were Louvain's, with no minimum."""
    said = (f"{run.get('method', 'louvain').capitalize()}, k {run['k']}, "
            f"resolution {run['resolution']:g}")
    if run.get("absorbed"):
        said += (f", {run['absorbed']} of fewer than {run['min_cluster_size']} "
                 "events merged")
    return said


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


# --- Leiden -------------------------------------------------------------------

#: passes of Leiden at most; the second and third still gain a little
#: (measured in the module's notes), later ones next to nothing
PASSES = 3

#: how random the refinement's choice among the merges that do not lower
#: modularity is, in units of weight; the value of Traag et al.
THETA = 0.01


@njit(cache=True)
def _move_queue(indptr, indices, weights, degree, community, total, resolution, seed):
    """Leiden's local moving, on the partition given in `community` (ids
    below the number of nodes), changed in place. As Louvain's, but only
    nodes whose neighbourhood has changed are visited again, and a node
    may leave for a community of its own."""
    n = degree.size
    inside = np.zeros(n)
    members = np.zeros(n, dtype=np.int64)
    for i in range(n):
        inside[community[i]] += degree[i]
        members[community[i]] += 1
    empty = np.empty(n, dtype=np.int64)
    n_empty = 0
    for c in range(n):
        if members[c] == 0:
            empty[n_empty] = c
            n_empty += 1
    link = np.zeros(n)
    seen = np.empty(n, dtype=np.int64)
    np.random.seed(seed)
    queue = np.random.permutation(n)
    waiting = np.ones(n, dtype=np.bool_)
    head, count = 0, n
    tolerance = 1e-12 * total
    while count:
        i = queue[head]
        head = (head + 1) % n
        count -= 1
        waiting[i] = False
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
        if best_gain < -tolerance and members[own] > 1 and n_empty:
            n_empty -= 1
            best = empty[n_empty]
        inside[best] += degree[i]
        if best != own:
            community[i] = best
            members[best] += 1
            members[own] -= 1
            if members[own] == 0:
                empty[n_empty] = own
                n_empty += 1
            for e in range(indptr[i], indptr[i + 1]):
                j = indices[e]
                if community[j] != best and not waiting[j]:
                    waiting[j] = True
                    queue[(head + count) % n] = j
                    count += 1


@njit(cache=True)
def _refine(indptr, indices, weights, degree, community, total, resolution, theta, seed):
    """Leiden's refinement: inside each community every node starts
    alone again and, once, may merge into a sub-community it is joined
    to, if both are well enough tied to the rest of the community and
    modularity does not fall. A merge needs a join, so every
    sub-community is connected. Returns the sub-community of each node."""
    n = degree.size
    refined = np.arange(n)
    size = np.ones(n, dtype=np.int64)
    sub = degree.copy()             # per sub-community: the degrees of its nodes
    whole = np.zeros(n)             # the same per community
    out = np.zeros(n)               # per sub-community: weight to the rest of its community
    for i in range(n):
        whole[community[i]] += degree[i]
        for e in range(indptr[i], indptr[i + 1]):
            j = indices[e]
            if j != i and community[j] == community[i]:
                out[i] += weights[e]
    link = np.zeros(n)
    seen = np.empty(n, dtype=np.int64)
    gains = np.empty(n)
    np.random.seed(seed)
    scale = resolution / total
    for i in np.random.permutation(n):
        if refined[i] != i or size[i] != 1:
            continue
        rest = whole[community[i]]
        if out[i] < scale * degree[i] * (rest - degree[i]):
            continue
        n_seen = 0
        for e in range(indptr[i], indptr[i + 1]):
            j = indices[e]
            if j == i or community[j] != community[i]:
                continue
            t = refined[j]
            if link[t] == 0.0:
                seen[n_seen] = t
                n_seen += 1
            link[t] += weights[e]
        n_fit, best_gain = 0, 0.0
        for s in range(n_seen):
            t = seen[s]
            gain = link[t] - scale * degree[i] * sub[t]
            if gain >= 0.0 and out[t] >= scale * sub[t] * (rest - sub[t]):
                seen[n_fit] = t
                gains[n_fit] = gain
                n_fit += 1
                best_gain = max(best_gain, gain)
            else:
                link[t] = 0.0
        if n_fit:
            # staying alone is one of the choices, at a gain of nothing
            chance = np.exp(-best_gain / theta)
            for s in range(n_fit):
                gains[s] = np.exp((gains[s] - best_gain) / theta)
                chance += gains[s]
            draw = np.random.random() * chance
            to = -1
            for s in range(n_fit):
                draw -= gains[s]
                if draw < 0.0:
                    to = seen[s]
                    break
            if to >= 0:
                refined[i] = to
                size[to] += 1
                sub[to] += degree[i]
                out[to] += out[i] - 2.0 * link[to]
            for s in range(n_fit):
                link[seen[s]] = 0.0
    return refined


def _leiden_pass(A0, member, resolution, seed):
    """One pass from the partition `member` of the nodes: move, refine,
    merge the sub-communities into nodes, and again until every
    community is one node."""
    from scipy import sparse
    A = A0
    where = np.arange(A.shape[0])   # of each node of the graph: its node in A
    _, start = np.unique(member, return_inverse=True)
    level = 0
    while True:
        n = A.shape[0]
        degree = np.asarray(A.sum(axis=1)).ravel()
        total = float(degree.sum())
        community = start.astype(np.int64)
        _move_queue(A.indptr, A.indices, A.data, degree, community, total,
                    float(resolution), (int(seed) + 2 * level) % (2 ** 31))
        kept, community = np.unique(community, return_inverse=True)
        if kept.size == n:
            return community[where]
        refined = _refine(A.indptr, A.indices, A.data, degree, community, total,
                          float(resolution), THETA, (int(seed) + 2 * level + 1) % (2 ** 31))
        kept, refined = np.unique(refined, return_inverse=True)
        if kept.size == n:   # nothing to merge; not met, and Louvain's answer if it were
            return _components_apart(A0, community[where])
        start = np.empty(kept.size, dtype=np.int64)
        start[refined] = community
        fold = sparse.csr_matrix((np.ones(n), (np.arange(n), refined)),
                                 shape=(n, kept.size))
        A = (fold.T @ A @ fold).tocsr()
        where = refined[where]
        level += 1


def _components_apart(graph, labels):
    """The partition with every community that is in several pieces cut
    into them."""
    pieces = components(graph, labels)
    return np.unique(np.column_stack([labels, pieces]), axis=0, return_inverse=True)[1].ravel()


def components(graph, labels) -> np.ndarray:
    """Of each node: the connected piece of its community it is in,
    numbered over the whole graph. A partition whose communities are all
    connected has as many pieces as communities."""
    from scipy import sparse
    A = sparse.coo_matrix(graph)
    labels = np.asarray(labels)
    same = labels[A.row] == labels[A.col]
    inner = sparse.csr_matrix((np.ones(same.sum()), (A.row[same], A.col[same])),
                              shape=A.shape)
    return sparse.csgraph.connected_components(inner, directed=False)[1]


def leiden(graph, resolution: float = 1.0, seed: int = 0,
           passes: int = PASSES) -> tuple[np.ndarray, float]:
    """As `louvain`, by the Leiden algorithm (Traag, Waltman & van Eck,
    Sci. Rep. 2019): between Louvain's two phases each community is
    refined into connected sub-communities, and it is those that become
    the nodes, so that no community is left in pieces. Each pass starts
    from the partition of the one before; they stop early when one
    changes nothing."""
    from scipy import sparse
    A = sparse.csr_matrix(graph, dtype=np.float64)
    A.sum_duplicates()
    member = np.arange(A.shape[0])
    if not A.nnz:
        return member, 0.0
    quality = modularity(A, member, resolution)
    for turn in range(max(1, int(passes))):
        found = _leiden_pass(A, member, resolution, int(seed) + 7919 * turn)
        q = modularity(A, found, resolution)
        if q <= quality + 1e-12:
            break
        member, quality = found, q
    return np.unique(member, return_inverse=True)[1], quality


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


def absorb_small(graph, labels: np.ndarray, min_size: int) -> tuple[np.ndarray, int]:
    """(the partition, how many clusters were taken into another). A
    cluster of fewer than `min_size` events joins the cluster it shares
    the most weight of joins with, the smallest first, until none is
    left below; one joined to no other stays as it is."""
    from scipy import sparse
    _, c = np.unique(labels, return_inverse=True)
    sizes = np.bincount(c).astype(np.int64)
    if min_size <= 1 or sizes.min() >= min_size:
        return c, 0
    n = c.size
    fold = sparse.csr_matrix((np.ones(n), (np.arange(n), c)), shape=(n, sizes.size))
    between = (fold.T @ sparse.csr_matrix(graph) @ fold).toarray()
    np.fill_diagonal(between, 0.0)
    into = np.arange(sizes.size)
    alive = np.ones(sizes.size, dtype=bool)
    taken = 0
    while True:
        small = np.flatnonzero(alive & (sizes < min_size) & (between.max(axis=1) > 0))
        if not small.size:
            break
        a = small[np.argmin(sizes[small])]
        b = int(np.argmax(between[a]))
        between[b] += between[a]
        between[:, b] += between[:, a]
        between[a] = between[:, a] = 0.0
        between[b, b] = 0.0
        sizes[b] += sizes[a]
        alive[a] = False
        into[into == a] = b
        taken += 1
    return np.unique(into[c], return_inverse=True)[1], taken


METHODS = {"leiden": leiden, "louvain": louvain}


def phenograph(X: np.ndarray, k: int = 30, resolution: float = 1.0,
               seed: int = 0, *, method: str = "leiden",
               min_cluster_size: int = 10) -> Clustering:
    """`X` is events x features, already transformed; nothing is scaled
    here. `k` is how many neighbours each event is joined to: fewer
    finds smaller clusters. Clusters of fewer than `min_cluster_size`
    events are taken into the cluster they are most joined to (1 keeps
    them all)."""
    X = np.ascontiguousarray(X, dtype=np.float64)
    nearest = neighbours(X, k, seed)
    graph = jaccard_graph(nearest)
    member, q = METHODS[method](graph, resolution, seed)
    member, absorbed = absorb_small(graph, member, int(min_cluster_size))
    if absorbed:
        q = modularity(graph, member, resolution)
    return Clustering(by_size(member), q, nearest.shape[1], float(resolution), int(seed),
                      method, int(min_cluster_size), absorbed)
