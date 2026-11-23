"""milc.normalize beside the CytoNorm package, on the same events.

    python3 tests/reference/normalize_compare.py

Two data sets. A simulated one whose truth is known: three batches, three
populations, a batch effect that differs between populations, batch "a"
without any — so that, brought onto "a", every event should come back to
the value it had before the effect. And the real experiment: the day as
the batch, the Ctrl animals as the reference, live cells, eight stains on
the biex scale, the goal the mean of the days.

Each is normalised without clusters and with ten metaclusters, by both
programs, each from two seeds: what separates two runs of one program is
the yardstick for what separates the two programs.
"""

from __future__ import annotations

import json
import os

import numpy as np

import common
from milc import normalize

CLUSTERS = 10
FLUORESCENCE = ["FITC-A", "PE-A", "PerCP-Cy5.5-A", "PE-Cy7-A", "APC-A", "APC-Cy7-A",
                "BV421-A", "BV510-A"]


def simulated() -> dict:
    """Events on a scale like an arcsinh's (0 to 5), four markers. The
    populations sit at different places and each batch moves and stretches
    each population differently — the case clusters are for."""
    rng = np.random.default_rng(11)
    centres = np.array([[1.0, 1.0, 3.5, 0.8], [3.2, 1.2, 1.0, 2.8], [1.2, 3.4, 2.6, 3.0]])
    spread = 0.25
    # shifts of up to 0.6 of a population's spread: enough to matter, not
    # enough for a batch's population to become a cluster of its own
    shift = {"a": np.zeros((3, 4)),
             "b": 0.5 * np.array([[0.30, 0.10, -0.20, 0.00], [-0.15, 0.25, 0.10, 0.20],
                                  [0.05, -0.25, 0.30, -0.10]]),
             "c": 0.5 * np.array([[-0.20, 0.00, 0.25, 0.15], [0.20, -0.10, -0.15, -0.25],
                                  [-0.10, 0.30, 0.00, 0.20]])}
    scale = {"a": 1.0, "b": 1.15, "c": 0.9}

    def draw(n, shares):
        which = rng.choice(3, n, p=shares)
        return centres[which] + rng.normal(0, spread, (n, 4)), which

    def effect(clean, which, b):
        return centres[which] + (clean - centres[which]) * scale[b] + shift[b][which]

    ref, x, truth = {}, {}, {}
    for b in "abc":
        clean, which = draw(20000, [0.5, 0.3, 0.2])
        ref[b] = effect(clean, which, b)
        # the samples to normalise have another composition than the reference
        clean, which = draw(10000, [0.2, 0.3, 0.5])
        x[b], truth[b] = effect(clean, which, b), clean
    return {"name": "simulated", "channels": ["M1", "M2", "M3", "M4"], "batches": list("abc"),
            "ref": ref, "x": x, "truth": truth, "goal": "a"}


def real() -> dict:
    from milc.transforms import Biex
    s, entries = common.real_session()
    t = Biex()
    ref, x = {}, {}
    for e in entries:
        sm = s.load(e)
        cols = [sm.fcs.index(c) for c in FLUORESCENCE]
        ev = t(sm.compensated[sm.mask(common.LIVE)][:, cols].astype(np.float64))
        (ref if e.meta["Group"] == "Ctrl" else x).setdefault(e.meta["Day"], []).append(ev)
    batches = ["3", "7", "14"]
    return {"name": "real", "channels": [c.replace("-", "_") for c in FLUORESCENCE],
            "batches": batches, "ref": {b: np.vstack(ref[b]) for b in batches},
            "x": {b: np.vstack(x[b]) for b in batches}, "truth": None, "goal": "mean"}


def runs(goal: str) -> list[dict]:
    return [{"name": "noclusters", "clusters": 1, "grid": 0, "seed": 1, "goal": goal},
            {"name": "c10_15x15_s1", "clusters": CLUSTERS, "grid": 15, "seed": 1, "goal": goal},
            {"name": "c10_15x15_s2", "clusters": CLUSTERS, "grid": 15, "seed": 2, "goal": goal},
            {"name": "c10_10x10_s1", "clusters": CLUSTERS, "grid": 10, "seed": 1, "goal": goal}]


def read(folder: str, run: str, data: str, kind: str, batch: str, p: int) -> np.ndarray:
    path = os.path.join(folder, f"{run}__{data}_{kind}_{batch}.bin")
    return np.fromfile(path, dtype="<f4").reshape(-1, p).astype(np.float64)


def uncorrected(model, reference: dict) -> float:
    """The share of the reference events that fall in a (batch, cluster)
    pair the model leaves as it was."""
    left = {(b, c) for b, c, _ in model.fallbacks}
    n = sum(int((model.som.map(reference[b]) == c).sum()) for b, c in left)
    return n / sum(len(v) for v in reference.values())


def gap(a: dict, b: dict, before: dict) -> list[float]:
    """Between two normalisations of the same events: the mean and the
    99th percentile of the absolute difference, and the mean as a share of
    the mean correction made (by the first)."""
    d = np.abs(np.vstack([a[k] - b[k] for k in a]))
    moved = np.abs(np.vstack([a[k] - before[k] for k in a])).mean()
    return [float(d.mean()), float(np.quantile(d, 0.99)), float(d.max()),
            float(100 * d.mean() / moved)]


def main() -> None:
    folder = common.work("normalize")
    rep = common.Report("normalize")
    rep.say("# normalize.fit / apply beside the CytoNorm package\n")
    datasets = [simulated()] + ([real()] if common.has_real() else [])
    manifest = {"datasets": []}
    for d in datasets:
        for b in d["batches"]:
            # 32-bit, as an FCS file holds them: both programs get these very numbers
            for kind in ("ref", "x"):
                d[kind][b] = d[kind][b].astype("<f4")
                d[kind][b].tofile(os.path.join(folder, f"{d['name']}_{kind}_{b}.bin"))
                d[kind][b] = d[kind][b].astype(np.float64)
        manifest["datasets"].append({"name": d["name"], "channels": d["channels"],
                                     "batches": d["batches"], "runs": runs(d["goal"])})
    with open(os.path.join(folder, "manifest.json"), "w") as f:
        json.dump(manifest, f)
    out = common.rscript("normalize_cytonorm.R", folder)
    rep.say("versions: " + out.strip().splitlines()[-1])
    found = {tuple(line.split()[:2]): line.split()[-1] for line in out.splitlines()
             if " metaclusters " in line}

    for d in datasets:
        name, batches, p = d["name"], d["batches"], len(d["channels"])
        n_ref = {b: len(d["ref"][b]) for b in batches}
        rep.say(f"\n## {name}: batches {', '.join(batches)}; reference events "
                f"{', '.join(f'{n:,d}' for n in n_ref.values())}; events to normalise "
                f"{', '.join(f'{len(d['x'][b]):,d}' for b in batches)}; {p} markers; "
                f"goal {d['goal']!r}")
        theirs = {r["name"]: {k: {b: read(folder, r["name"], name, k, b, p) for b in batches}
                              for k in ("ref", "x")} for r in runs(d["goal"])}
        ours, models = {}, {}
        for label, kw in (("noclusters", dict(n_clusters=1)),
                          ("c10_s0", dict(n_clusters=CLUSTERS, seed=0)),
                          ("c10_s1", dict(n_clusters=CLUSTERS, seed=1)),
                          ("c10_15x15_s0", dict(n_clusters=CLUSTERS, seed=0, grid=15))):
            m = normalize.fit(d["ref"], goal=d["goal"], **kw)
            models[label] = m
            ours[label] = {k: {b: normalize.apply(m, d[k][b], b) for b in batches}
                           for k in ("ref", "x")}

        # --- without clusters: the same curve? --------------------------------
        rows = []
        for b in batches:
            x, o, t = d["x"][b], ours["noclusters"]["x"][b], theirs["noclusters"]["x"][b]
            lo, hi = np.quantile(d["ref"][b], [0.01, 0.99], axis=0)
            inside = (x >= lo) & (x <= hi)
            moved = np.abs(o - x)[inside].mean()
            g = np.abs(o - t)
            rows.append([b, float(moved), float(g[inside].mean()),
                         float(np.quantile(g[inside], 0.99)), float(g[inside].max()),
                         100 * float(1 - inside.mean()),
                         float(g[~inside].mean()) if (~inside).any() else float("nan"),
                         float(g[~inside].max()) if (~inside).any() else float("nan")])
        rep.table(f"{name}, no clusters: ours against QuantileNorm, the events to normalise "
                  "(scale units)",
                  ["batch", "mean correction (ours), inside", "mean |ours − CytoNorm| inside "
                   "the reference's 1–99%", "p99 inside", "max inside", "values outside %",
                   "mean |ours − CytoNorm| outside", "max outside"], rows)

        # --- with clusters -------------------------------------------------------
        before = d["x"]
        rows = [["ours seed 0 vs ours seed 1 (10 x 10)"]
                + gap(ours["c10_s0"]["x"], ours["c10_s1"]["x"], before),
                ["CytoNorm seed 1 vs CytoNorm seed 2 (15 x 15)"]
                + gap(theirs["c10_15x15_s1"]["x"], theirs["c10_15x15_s2"]["x"], before),
                ["ours (10 x 10, seed 0) vs CytoNorm (15 x 15, seed 1)"]
                + gap(ours["c10_s0"]["x"], theirs["c10_15x15_s1"]["x"], before),
                ["ours (10 x 10, seed 0) vs CytoNorm (10 x 10, seed 1)"]
                + gap(ours["c10_s0"]["x"], theirs["c10_10x10_s1"]["x"], before),
                ["ours (15 x 15, seed 0) vs CytoNorm (15 x 15, seed 1)"]
                + gap(ours["c10_15x15_s0"]["x"], theirs["c10_15x15_s1"]["x"], before),
                ["ours with clusters vs ours without"]
                + gap(ours["c10_s0"]["x"], ours["noclusters"]["x"], before),
                ["CytoNorm with clusters vs CytoNorm without"]
                + gap(theirs["c10_15x15_s1"]["x"], theirs["noclusters"]["x"], before)]
        rep.table(f"{name}, {CLUSTERS} metaclusters: how far apart two normalisations of the "
                  "same events are (scale units)",
                  ["pair", "mean |difference|", "p99", "max", "mean as % of the mean correction"],
                  rows)
        rep.say(f"\nmetaclusters found — ours: "
                + ", ".join(f"{k} {m.n_clusters}" for k, m in models.items() if k != "noclusters")
                + "; CytoNorm: "
                + ", ".join(f"{r} {n}" for (dn, r), n in found.items() if dn == name and r != "noclusters")
                + f"; (batch, cluster) pairs ours left uncorrected for want of "
                f"{50} reference events, and the share of the reference events in them: "
                + ", ".join(f"{k} {len(m.fallbacks)} ({100 * uncorrected(m, d['ref']):.2f}%)"
                            for k, m in models.items() if k != "noclusters"))

        # --- what it is for: are the batches closer? -------------------------------
        every = [("before", {k: d[k] for k in ("ref", "x")}),
                 ("ours, no clusters", ours["noclusters"]),
                 ("CytoNorm, no clusters", theirs["noclusters"]),
                 ("ours, 10 clusters (10 x 10)", ours["c10_s0"]),
                 ("ours, 10 clusters (15 x 15)", ours["c10_15x15_s0"]),
                 ("CytoNorm, 10 clusters (15 x 15)", theirs["c10_15x15_s1"]),
                 ("CytoNorm, 10 clusters (10 x 10)", theirs["c10_10x10_s1"])]
        rows = []
        for label, got in every:
            row = [label, float(normalize.distance(got["ref"], goal=d["goal"]).mean()),
                   float(normalize.distance(got["x"], goal=d["goal"]).mean())]
            if d["truth"] is not None:
                others = [b for b in batches if b != d["goal"]]
                err = np.abs(np.vstack([got["x"][b] - d["truth"][b] for b in others]))
                row += [float(err.mean()), float(np.quantile(err, 0.99))]
            # how far a value was carried, and whether it left the scale: beyond
            # the outermost quantiles CytoNorm follows the spline's end slope
            x = np.vstack([d["x"][b] for b in batches])
            y = np.vstack([got["x"][b] for b in batches])
            width = x.max() - x.min()
            row += [float(np.abs(y - x).max()), float(y.min()), float(y.max()),
                    int(((y < x.min() - 0.1 * width) | (y > x.max() + 0.1 * width)).sum())]
            rows.append(row)
        header = ["normalisation", "references: distance between batches (mean over markers)",
                  "events to normalise: distance between batches"]
        if d["truth"] is not None:
            header += ["mean |normalised − true value| (batches b, c)", "p99"]
        header += ["largest correction of any value", "lowest value after", "highest after",
                   "values carried more than 10% of the scale beyond the input's range"]
        rep.say(f"\nthe events to normalise span {min(d['x'][b].min() for b in batches):.3g} "
                f"to {max(d['x'][b].max() for b in batches):.3g}")
        rep.table(f"{name}: what is left of the batch effect (normalize.distance, scale units)",
                  header, rows)
    rep.close()


if __name__ == "__main__":
    main()
