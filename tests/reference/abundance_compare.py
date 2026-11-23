"""milc.abundance.differential beside edgeR and diffcyt, on the same count
tables: simulated ones whose truth is known and the real experiment's.

    python3 tests/reference/abundance_compare.py [tables per setting]

See README.md here for what has to be installed.
"""

from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

import common
from milc import abundance

#: which of ours each edgeR variant is set beside (see abundance_edger.R)
VARIANTS = {"ql_matched": abundance.QL, "ql_default": abundance.QL,
            "lrt_matched": abundance.LRT, "lrt_disp": abundance.LRT,
            "diffcyt": abundance.LRT}
SETTINGS = [((5, 5), 15), ((3, 3), 15), ((5, 5), 8)]


def simulate(rng, n=(5, 5), clusters=15, fold=1.0, changed=3):
    """The simulation of tests/test_abundance.py, copied: shares from 0.1%
    to 30%, a dispersion of its own for each cluster (0.02 to 0.3), totals
    varying five-fold, counts multinomial; with `fold`, `changed` clusters
    are that many times larger in A and the largest gives up the
    difference."""
    p = np.geomspace(0.001, 0.3, clusters)
    p /= p.sum()
    phi = np.exp(rng.uniform(np.log(0.02), np.log(0.3), clusters))
    truth = np.zeros(clusters, bool)
    pa = p.copy()
    if fold != 1.0:
        truth[rng.choice(np.flatnonzero(p < 0.1), changed, replace=False)] = True
        pa[truth] *= fold
        pa[-1] -= pa.sum() - 1
    groups = ["A"] * n[0] + ["B"] * n[1]
    counts = np.empty((clusters, sum(n)))
    for i, g in enumerate(groups):
        q = (pa if g == "A" else p) * rng.gamma(1 / phi, phi)
        total = int(np.exp(rng.uniform(np.log(10000), np.log(50000))))
        counts[:, i] = rng.multinomial(total, q / q.sum())
    return counts, groups, truth


def simulated(reps: int) -> list[dict]:
    out = []
    rng = np.random.default_rng(2026)
    for n, clusters in SETTINGS:
        for fold in (1.0, 2.0):
            for r in range(reps):
                counts, groups, truth = simulate(rng, n, clusters, fold)
                judged = np.ones(clusters, bool)
                judged[-1] = fold == 1.0
                out.append({"id": f"sim {n[0]}v{n[1]} k{clusters} fold{fold:g} #{r}",
                            "kind": f"{n[0]} v {n[1]}, {clusters} clusters",
                            "fold": fold, "counts": counts.tolist(), "groups": groups,
                            "contrast": ["A", "B"], "covariate": None,
                            "truth": truth.tolist(), "judged": judged.tolist()})
    return out


def real() -> list[dict]:
    """CD45+ events of the 30 samples, 12 FlowSOM metaclusters on six
    stains, counted per sample; Case against Ctrl within each day, and over
    the three days with the day adjusted for."""
    from milc.flowsom import flowsom
    X, which, entries = common.real_events()
    labels = flowsom(X, 12, 10, 10, seed=0).labels
    counts = np.array([[np.sum((labels == c) & (which == i)) for i in range(len(entries))]
                       for c in range(1, labels.max() + 1)], dtype=float)
    group = [e.meta["Group"] for e in entries]
    day = [e.meta["Day"] for e in entries]
    out = []
    for d in ("3", "7", "14"):
        cols = [i for i, x in enumerate(day) if x == d]
        out.append({"id": f"real day {d}", "kind": "real", "fold": None,
                    "counts": counts[:, cols].tolist(), "groups": [group[i] for i in cols],
                    "contrast": ["Case", "Ctrl"], "covariate": None})
    out.append({"id": "real all days, day adjusted", "kind": "real", "fold": None,
                "counts": counts.tolist(), "groups": group, "contrast": ["Case", "Ctrl"],
                "covariate": day})
    return out


def ours(t: dict, test: str) -> dict[str, np.ndarray]:
    counts = np.array(t["counts"])
    found = abundance.differential(counts, [str(i) for i in range(len(counts))], t["groups"],
                                   tuple(t["contrast"]), t["covariate"], test)
    return {k: np.array([getattr(r, k) for r in found])
            for k in ("log2_fold", "dispersion", "statistic", "p", "p_adjusted")}


def read_edger(path: str) -> dict:
    out: dict = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            d = out.setdefault((row["table"], row["variant"]), {k: [] for k in
                               ("log2fc", "dispersion", "statistic", "p", "fdr", "common",
                                "prior_df")})
            for k in d:
                d[k].append(float(row[k]) if row[k] not in ("NA", "") else float("nan"))
    return {k: {a: np.array(b) for a, b in v.items()} for k, v in out.items()}


def spearman(a, b) -> float:
    from scipy.stats import spearmanr
    return float(spearmanr(a, b)[0])


def agreement(tables, mine, theirs, variant) -> list:
    """One row: how far ours and one edgeR variant are over every cluster
    of `tables`."""
    test = VARIANTS[variant]
    fc, disp, p, calls_o, calls_e, rho = [], [], [], [], [], []
    for t in tables:
        o, e = mine[(t["id"], test)], theirs[(t["id"], variant)]
        ok = np.isfinite(o["p"]) & np.isfinite(e["p"])
        fc += list(np.abs(o["log2_fold"] - e["log2fc"])[ok])
        with np.errstate(divide="ignore", invalid="ignore"):
            disp += list(np.abs(np.log2(o["dispersion"] / e["dispersion"]))[ok])
            p += list(np.abs(np.log10(np.maximum(o["p"], 1e-300))
                             - np.log10(np.maximum(e["p"], 1e-300)))[ok])
        calls_o += list(o["p_adjusted"][ok] < 0.05)
        calls_e += list(e["fdr"][ok] < 0.05)
        if ok.sum() > 2:
            rho.append(spearman(o["p"][ok], e["p"][ok]))
    fc, disp, p = np.array(fc), np.array(disp), np.array(p)
    disp = disp[np.isfinite(disp)]
    co, ce = np.array(calls_o), np.array(calls_e)
    q = lambda v, x: float(np.quantile(v, x)) if v.size else float("nan")
    return [variant, test.upper(), len(fc), q(fc, .5), q(fc, .95), q(fc, 1),
            q(disp, .5), q(disp, .95), q(p, .5), q(p, .95), q(p, 1), int((p > 0.1).sum()),
            float(np.median(rho)), int((co & ce).sum()), int((co & ~ce).sum()),
            int((~co & ce).sum())]


AGREEMENT = ["edgeR variant", "ours", "clusters", "|Δlog2FC| median", "p95", "max",
             "|log2 disp. ratio| median", "p95", "|Δlog10 p| median", "p95", "max",
             "clusters with |Δlog10 p| > 0.1",
             "Spearman of p (median table)", "FDR<0.05 both", "only ours", "only edgeR"]


def calibration(tables, p_of, q_of) -> tuple[float, float, float]:
    """(share of unchanged clusters with p < 0.05 under the null; share of
    changed clusters found at FDR < 0.05; false discovery rate had)."""
    null_p, fdp, hit, true = [], [], 0, 0
    for t in tables:
        truth, judged = np.array(t["truth"]), np.array(t["judged"])
        p, q = p_of(t), q_of(t)
        if t["fold"] == 1.0:
            null_p += list(p[judged])
        else:
            sig = (q < 0.05) & judged
            fdp.append((sig & ~truth).sum() / max(sig.sum(), 1))
            hit += (sig & truth).sum()
            true += truth.sum()
    return (float(np.mean(np.array(null_p) < 0.05)), hit / max(true, 1), float(np.mean(fdp)))


def main(reps: int = 100) -> None:
    folder = common.work("abundance")
    rep = common.Report("abundance")
    tables = simulated(reps) + (real() if common.has_real() else [])
    with open(os.path.join(folder, "tables.json"), "w") as f:
        json.dump(tables, f)
    rep.say("# abundance.differential beside edgeR\n")
    rep.say("versions: " + common.rscript("abundance_edger.R", folder).strip().splitlines()[-1])
    theirs = read_edger(os.path.join(folder, "edger.csv"))
    mine = {(t["id"], test): ours(t, test) for t in tables
            for test in (abundance.QL, abundance.LRT)}
    variants = [v for v in VARIANTS if (tables[0]["id"], v) in theirs]

    sim = [t for t in tables if t["kind"] != "real"]
    rep.say(f"\n{len(sim)} simulated tables ({reps} for each of {len(SETTINGS)} settings, "
            "with no change and with a two-fold one)"
            + ("" if common.has_real() else "; no real data (MILC_REF_DATA not set)"))
    rep.table("Simulated tables: ours against each edgeR variant", AGREEMENT,
              [agreement(sim, mine, theirs, v) for v in variants])

    # the common dispersion and the prior degrees of freedom, table by table
    rows = []
    for v, test, label in (("ql_matched", abundance.QL, "common dispersion (log2 ratio)"),
                           ("lrt_disp", abundance.LRT, None)):
        if label:
            d = np.array([np.log2(mine[(t["id"], test)]["dispersion"][0]
                                  / theirs[(t["id"], v)]["common"][0]) for t in sim])
            rows.append([label, float(np.median(d)), float(np.quantile(np.abs(d), .95)),
                         float(np.abs(d).max())])
    rep.table("Simulated tables: the common dispersion, ours / estimateGLMCommonDisp",
              ["quantity", "median", "p95 of |.|", "max of |.|"], rows)
    pd = np.array([theirs[(t["id"], "lrt_disp")]["prior_df"][0] for t in sim])
    pq = np.array([theirs[(t["id"], "ql_matched")]["prior_df"][0] for t in sim])
    at = lambda v, x: float(np.sort(v)[int(x * (len(v) - 1))])     # quantiles, infinity allowed
    rep.say(f"\nedgeR's own prior degrees of freedom on these tables: estimateDisp "
            f"(for the dispersions) median {at(pd, .5):.3g}, 5%–95% "
            f"{at(pd, .05):.3g}–{at(pd, .95):.3g}, infinite in "
            f"{100 * np.mean(np.isinf(pd)):.0f}% of tables (ours: fixed at "
            f"{abundance.PRIOR_DF:g}); glmQLFit legacy (for the quasi-dispersions) median "
            f"{at(pq, .5):.3g}, infinite in {100 * np.mean(np.isinf(pq)):.0f}% of tables")

    rows = []
    for n, clusters in SETTINGS:
        kind = f"{n[0]} v {n[1]}, {clusters} clusters"
        these = [t for t in sim if t["kind"] == kind]
        for name, p_of, q_of in (
                [(f"ours {m.upper()}", lambda t, m=m: mine[(t["id"], m)]["p"],
                  lambda t, m=m: mine[(t["id"], m)]["p_adjusted"])
                 for m in (abundance.QL, abundance.LRT)]
                + [(f"edgeR {v}", lambda t, v=v: theirs[(t["id"], v)]["p"],
                    lambda t, v=v: theirs[(t["id"], v)]["fdr"]) for v in variants]):
            a, b, c = calibration(these, p_of, q_of)
            rows.append([kind, name, 100 * a, 100 * b, 100 * c])
    rep.table(f"The same simulated tables, each method on its own ({reps} null and "
              f"{reps} two-fold tables per setting)",
              ["setting", "method", "null p<0.05 (%)", "power at FDR<0.05 (%)",
               "FDR had (%)"], rows)

    if common.has_real():
        realt = [t for t in tables if t["kind"] == "real"]
        rep.table("Real tables (12 FlowSOM clusters of CD45+, Case vs Ctrl): ours against "
                  "each edgeR variant", AGREEMENT,
                  [agreement(realt, mine, theirs, v) for v in variants])
        rows = []
        for t in realt:
            o, e, d = (mine[(t["id"], abundance.QL)], theirs[(t["id"], "ql_matched")],
                       theirs[(t["id"], "ql_default")])
            ol, el = mine[(t["id"], abundance.LRT)], theirs[(t["id"], "diffcyt")]
            for k in range(len(o["p"])):
                rows.append([t["id"], k + 1, o["log2_fold"][k], e["log2fc"][k], o["p"][k],
                             e["p"][k], d["p"][k], ol["p"][k], el["p"][k]])
        rep.table("Real tables, cluster by cluster",
                  ["table", "cluster", "log2FC ours", "log2FC edgeR", "p ours QL",
                   "p edgeR ql_matched", "p edgeR ql_default", "p ours LRT", "p diffcyt"], rows)
    rep.close()


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 100)
