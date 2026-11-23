"""
Two-dimensional maps of the events: UMAP and tSNE.

The defaults are the ones the FlowJo plugins ran with, so a first map is
comparable to the one it replaces. Neither algorithm gives the same
picture twice without the same seed, and never the same picture as
another implementation: what can be compared is how well each keeps the
neighbours an event has in the measured space (`neighbour_preservation`).

The heavy imports are made inside the functions: reading files and
planning a concatenation should not wait for numba to compile.
"""

from __future__ import annotations

import numpy as np


def umap(X: np.ndarray, n_neighbors: int = 15, min_dist: float = 0.5,
         metric: str = "euclidean", seed: int = 0,
         n_components: int = 2) -> np.ndarray:
    import warnings

    import umap as _umap
    with warnings.catch_warnings():
        # a seed makes it single-threaded, and it says so every time
        warnings.simplefilter("ignore", UserWarning)
        warnings.simplefilter("ignore", FutureWarning)
        kw = dict(n_neighbors=min(n_neighbors, max(2, X.shape[0] - 1)),
                  min_dist=min_dist, metric=metric,
                  n_components=n_components, random_state=int(seed))
        try:
            # 0.6 announces a new default layout; stay on the one measured
            model = _umap.UMAP(compatibility_layout=True, **kw)
        except TypeError:
            model = _umap.UMAP(**kw)
        return np.asarray(model.fit_transform(np.asarray(X, dtype=np.float32)),
                          dtype=np.float32)


def tsne(X: np.ndarray, perplexity: float = 30.0, n_iter: int = 1000,
         metric: str = "euclidean", seed: int = 0) -> np.ndarray:
    from openTSNE import TSNE
    n = X.shape[0]
    model = TSNE(perplexity=min(perplexity, max(1.0, (n - 1) / 3)),
                 n_iter=max(int(n_iter) - 250, 250),  # 250 more of early exaggeration
                 metric=metric, random_state=int(seed), n_jobs=-1)
    return np.asarray(model.fit(np.asarray(X, dtype=np.float64)), dtype=np.float32)


def neighbour_preservation(X: np.ndarray, Y: np.ndarray, k: int = 15,
                           max_events: int = 20000, seed: int = 0) -> float:
    """Of each event's `k` nearest neighbours in the measured space `X`,
    the share that are still among its `k` nearest on the map `Y`; averaged
    over the events. 1 is a perfect map, k/n is chance."""
    from sklearn.neighbors import NearestNeighbors
    n = X.shape[0]
    if n > max_events:
        rows = np.random.default_rng(seed).choice(n, max_events, replace=False)
        X, Y, n = X[rows], Y[rows], max_events
    k = min(k, n - 1)
    a = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X, return_distance=False)[:, 1:]
    b = NearestNeighbors(n_neighbors=k + 1).fit(Y).kneighbors(Y, return_distance=False)[:, 1:]
    kept = sum(len(np.intersect1d(a[i], b[i], assume_unique=True)) for i in range(n))
    return kept / (n * k)
