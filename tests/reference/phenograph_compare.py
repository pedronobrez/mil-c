"""milc.graphcluster.phenograph beside the authors' PhenoGraph package
(PyPI `PhenoGraph`), on the same events.

    <venv>/bin/python tests/reference/phenograph_compare.py

Only the public call `graphcluster.phenograph(X, k=, resolution=, seed=)`
is used, so this can be run again whatever the function does inside
(Louvain or Leiden). The package is run both ways it offers: its own
Louvain (the C++ program of Blondel et al. it ships, best of many random
starts, no seed) and Leiden (leidenalg).
"""

from __future__ import annotations

import time

import numpy as np
from sklearn.metrics import adjusted_rand_score as ari
from sklearn.metrics import normalized_mutual_info_score as nmi

import common
from milc import graphcluster

K = 30
RUNS = 3
SHARES = (0.40, 0.25, 0.18, 0.10, 0.06, 0.01)


def mixture(shares, apart, n=6000, p=8, seed=0):
    """As tests/test_graphcluster.py: Gaussian populations of unit spread,
    their centres `apart` from one another along different axes."""
    rng = np.random.default_rng(seed)
    truth = rng.choice(len(shares), n, p=shares)
    centres = np.zeros((len(shares), p))
    for i in range(len(shares)):
        centres[i, i % p] = apart / np.sqrt(2) * (1 + i // p)
    return centres[truth] + rng.normal(0, 1, (n, p)), truth


def modularity(graph, labels) -> float:
    """Newman's, of a partition of a symmetric weighted graph, written out
    here so that both programs' partitions are scored by the same code on
    the same graph."""
    from scipy import sparse
    A = sparse.csr_matrix(graph, dtype=np.float64)
    A = (A + A.T) / 2
    total = A.sum()
    _, c = np.unique(labels, return_inverse=True)
    n = A.shape[0]
    fold = sparse.csr_matrix((np.ones(n), (np.arange(n), c)), shape=(n, c.max() + 1))
    within = (fold.T @ A @ fold).diagonal()
    degree = np.bincount(c, weights=np.asarray(A.sum(axis=1)).ravel())
    return float(np.sum(within / total - (degree / total) ** 2))


def pairs(a: list, b: list | None = None, score=ari) -> float:
    """Mean score between every run of `a` and every run of `b` (or every
    two different runs of `a`)."""
    if b is None:
        found = [score(x, y) for i, x in enumerate(a) for y in a[i + 1:]]
    else:
        found = [score(x, y) for x in a for y in b]
    return float(np.mean(found)) if found else float("nan")


def case(name: str, X: np.ndarray, truth, rep: common.Report) -> None:
    import phenograph
    runs: dict[str, dict] = {}

    t0 = time.time()
    found = [graphcluster.phenograph(X, k=K, resolution=1.0, seed=s) for s in range(RUNS)]
    method = getattr(found[0], "method", "louvain")
    runs[f"ours ({method})"] = {"labels": [r.labels for r in found],
                                "q": [r.modularity for r in found],
                                "seconds": (time.time() - t0) / RUNS}
    graph = None
    for algo in ("louvain", "leiden"):
        labels, q = [], []
        t0 = time.time()
        for s in range(RUNS):
            lab, graph, quality = phenograph.cluster(X, k=K, clustering_algo=algo, seed=s,
                                                     n_jobs=1)
            labels.append(lab)
            q.append(float(quality))
        runs[f"package ({algo})"] = {"labels": labels, "q": q,
                                     "seconds": (time.time() - t0) / RUNS}

    rows = []
    for who, r in runs.items():
        n = [len(set(lab[lab >= 0])) for lab in r["labels"]] if "package" in who \
            else [int(lab.max()) for lab in r["labels"]]
        unassigned = max(int((lab < 0).sum()) for lab in r["labels"]) if "package" in who else 0
        rows.append([name, who, len(X), "/".join(map(str, n)), unassigned,
                     float(np.mean(r["q"])),
                     float(np.mean([modularity(graph, lab) for lab in r["labels"]])),
                     pairs(r["labels"]),
                     pairs(r["labels"], [truth]) if truth is not None else float("nan"),
                     pairs(r["labels"], [truth], nmi) if truth is not None else float("nan"),
                     r["seconds"]])
    rep.table(f"{name}: each program on its own ({RUNS} runs)",
              ["data", "program", "events", "clusters (each run)", "events left unassigned",
               "modularity it reports", "modularity on the package's graph",
               "ARI between its own runs", "ARI with the truth", "NMI with the truth",
               "seconds per run"], rows)
    mine = next(iter(runs))
    rows = [[name, f"{mine} vs {who}", pairs(runs[mine]["labels"], r["labels"]),
             pairs(runs[mine]["labels"], r["labels"], nmi)]
            for who, r in runs.items() if who != mine]
    rows.append([name, "package (louvain) vs package (leiden)",
                 pairs(runs["package (louvain)"]["labels"], runs["package (leiden)"]["labels"]),
                 pairs(runs["package (louvain)"]["labels"], runs["package (leiden)"]["labels"],
                       nmi)])
    rep.table(f"{name}: one program against the other (mean over {RUNS} x {RUNS} runs)",
              ["data", "pair", "ARI", "NMI"], rows)


def main() -> None:
    import phenograph
    rep = common.Report("phenograph")
    rep.say("# graphcluster.phenograph beside the PhenoGraph package\n")
    rep.say(f"PhenoGraph {phenograph.__version__}, k = {K}, resolution 1, {RUNS} runs each")
    for apart, n in ((10, 6000), (6, 6000), (4.5, 6000), (10, 20000)):
        X, truth = mixture(SHARES, apart, n)
        case(f"mixture of 6, centres {apart:g} SD apart", X, truth, rep)
    if common.has_real():
        X, _, _ = common.real_events()
        rng = np.random.default_rng(0)
        take = np.sort(rng.choice(len(X), min(20000, len(X)), replace=False))
        case("real CD45+ events, 6 stains (biex)", X[take], None, rep)
    else:
        rep.say("\nno real data (MILC_REF_DATA not set)")
    rep.close()


if __name__ == "__main__":      # the package starts processes: this guard is needed
    main()
