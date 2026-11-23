"""milc.qc.check beside PeacoQC (and flowAI): which events each removes.

    python3 tests/reference/qc_compare.py

Three kinds of acquisition: synthetic ones with one known defect each (as
tests/test_qc.py makes them), every real file as it is, and one clean real
file with a defect put into it. `qc.check` is called as the program calls
it (`Session.quality`): every parameter but Time, fluorescence
compensated, each on its display scale. PeacoQC gets those same values.
"""

from __future__ import annotations

import json
import os

import numpy as np

import common
from milc import qc
from milc.transforms import Biex


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    """Of the two sets of removed events; 1 when neither removes any."""
    union = (a | b).sum()
    return float((a & b).sum() / union) if union else 1.0


# --- the acquisitions --------------------------------------------------------------

def synthetic() -> list[dict]:
    """60 s at 1,000 events/s, two scatters and two stains with two
    populations each; one defect each, its events known."""
    channels = ["FSC-A", "SSC-A", "FITC-A", "PE-A"]
    scales = [None, None, Biex(), Biex()]

    def cells(rng, n):
        bright = rng.random(n) < 0.4
        return np.column_stack([
            rng.normal(90000, 9000, n), rng.normal(40000, 6000, n),
            np.where(bright, 3000, 150) * rng.lognormal(0, 0.3, n) + rng.normal(0, 40, n),
            np.where(bright, 200, 5000) * rng.lognormal(0, 0.3, n) + rng.normal(0, 40, n)])

    def make(seed):
        rng = np.random.default_rng(seed)
        t = np.sort(rng.uniform(0, 60, rng.poisson(60000)))
        return rng, t, cells(rng, len(t))

    out = []

    def add(name, t, x, truth):
        # Time in ticks of 10 ms, as a cytometer writes it ($TIMESTEP 0.01)
        out.append({"name": name, "kind": "synthetic", "raw": x, "time": t * 100.0,
                    "channels": channels, "scales": scales, "compensate": None,
                    "timestep": 0.01, "truth": truth})

    _, t, x = make(0)
    add("synthetic clean", t, x, None)
    rng, t, x = make(1)
    stay = ~((t >= 20) & (t <= 30)) | (rng.random(len(t)) < 0.03)
    add("synthetic clog (rate 3% for 10 s)", t[stay], x[stay], ((t >= 20) & (t <= 30))[stay])
    _, t, x = make(3)
    x[t >= 45, 3] *= 1.8
    add("synthetic step (PE x1.8, last 15 s)", t, x, t >= 45)
    _, t, x = make(4)
    x[:, 2] *= 1 + 2 * np.clip((t - 40) / 20, 0, None)
    add("synthetic drift (FITC to x3, last 20 s)", t, x, t >= 40)
    _, t, x = make(5)
    x[(t >= 10) & (t <= 12)] *= 0.5
    add("synthetic dip (all x0.5 for 2 s)", t, x, (t >= 10) & (t <= 12))
    return out


def real() -> list[dict]:
    s, _ = common.real_session()
    out = []
    for e in s.entries:
        if e.problem or not s.available(e):
            continue
        sample = s.load(e)
        f = sample.fcs
        time = qc.time_channel(f.channels)
        if time is None:
            continue
        names = [c for c in f.channels if c != time]
        cols = [f.index(c) for c in names]
        comp = sample.compensation
        out.append({"name": e.name, "kind": "real",
                    "raw": np.asarray(f.events[:, cols], dtype=np.float64),
                    "time": np.asarray(sample.column(time), dtype=np.float64),
                    "channels": names, "scales": [s.scale_of(c) for c in names],
                    "compensate": (lambda x, comp=comp, names=names: comp.apply(x, names))
                    if comp is not None else None,
                    "fluorescence": [i for i, c in enumerate(names)
                                     if comp is not None and c in comp.channels],
                    "timestep": float(f.keywords.get("$TIMESTEP", 1.0)), "truth": None})
    return out


def injected(base: dict) -> list[dict]:
    """The same file three times, each with one defect written into the
    values as acquired, over a stretch of its events in time order."""
    order = np.argsort(base["time"], kind="stable")
    rank = np.empty(len(order), dtype=np.int64)
    rank[order] = np.arange(len(order))
    share = rank / len(rank)
    fl = base["fluorescence"]
    out = []

    def add(what, raw, time, truth):
        out.append({**base, "name": f"{base['name']} + {what}", "kind": "injected",
                    "raw": raw, "time": time, "truth": truth})

    raw = base["raw"].copy()
    hit = (share >= 0.40) & (share < 0.50)
    raw[np.ix_(hit, fl)] *= 0.5
    add("step (fluorescence x0.5, 40-50% of the events)", raw, base["time"], hit)

    raw = base["raw"].copy()
    hit = share >= 0.70
    raw[np.ix_(hit, fl)] *= (1 - 0.5 * (share[hit] - 0.70) / 0.30)[:, None]
    add("drift (fluorescence down to x0.5 over the last 30%)", raw, base["time"], hit)

    # a clog: the events of 60-63% arrive ten times slower; nothing else changes
    t = base["time"][order]
    gaps = np.diff(t, prepend=t[0])
    slow = (np.arange(len(t)) >= 0.60 * len(t)) & (np.arange(len(t)) < 0.63 * len(t))
    gaps[slow] *= 10
    time = np.empty_like(t)
    time[order] = t[0] + np.cumsum(gaps)
    add("clog (60-63% of the events ten times slower)", base["raw"], time, slow[rank])
    return out


# --- ours, and the files for R ---------------------------------------------------

def run_ours(a: dict) -> qc.Result:
    x = a["compensate"](a["raw"]) if a["compensate"] else a["raw"]
    a["scaled"] = np.column_stack([sc(x[:, j]) if sc is not None else x[:, j]
                                   for j, sc in enumerate(a["scales"])])
    return qc.check(x, a["time"], a["channels"], a["scales"])


def export(folder: str, acquisitions: list[dict]) -> None:
    manifest = []
    for i, a in enumerate(acquisitions):
        stem = f"a{i:03d}"
        a["scaled"].astype("<f4").tofile(os.path.join(folder, stem + ".scaled.bin"))
        a["raw"].astype("<f4").tofile(os.path.join(folder, stem + ".raw.bin"))
        a["time"].astype("<f8").tofile(os.path.join(folder, stem + ".time.bin"))
        a["stem"] = stem
        manifest.append({"name": stem, "n": len(a["time"]), "channels": a["channels"],
                         "scaled": stem + ".scaled.bin", "raw": stem + ".raw.bin",
                         "time": stem + ".time.bin", "timestep": a["timestep"]})
    with open(os.path.join(folder, "manifest.json"), "w") as f:
        json.dump(manifest, f)


def flags(folder: str, stem: str, which: str, n: int) -> np.ndarray | None:
    """True where the program removes the event; None where it did not run."""
    path = os.path.join(folder, f"{stem}.{which}.bin")
    if not os.path.exists(path):
        return None
    good = np.fromfile(path, dtype=np.uint8)
    return good[:n] == 0


PROGRAMS = [("PeacoQC", "peacoqc_all"), ("PeacoQC MAD only", "peacoqc_MAD"), ("flowAI", "flowai")]


def main() -> None:
    folder = common.work("qc")
    for old in os.listdir(folder):
        if old.endswith(".bin"):
            os.remove(os.path.join(folder, old))
    rep = common.Report("qc")
    rep.say("# qc.check beside PeacoQC and flowAI\n")
    acquisitions = synthetic()
    if common.has_real():
        files = real()
        results = [run_ours(a) for a in files]
        # the largest file ours leaves whole is the one the defects are put into
        clean = [a for a, r in zip(files, results) if r.checked and r.removed < 0.005]
        acquisitions += files
        if clean:
            acquisitions += injected(max(clean, key=lambda a: len(a["time"])))
    else:
        rep.say("no real data (MILC_REF_DATA not set): synthetic acquisitions only\n")
    found = [run_ours(a) for a in acquisitions]
    export(folder, acquisitions)
    rep.say("versions: " + common.rscript("qc_peacoqc.R", folder).strip().splitlines()[-1])
    with open(os.path.join(folder, "notes.json")) as f:
        notes = json.load(f)

    rows, truth_rows, summary = [], [], {name: [] for name, _ in PROGRAMS}
    for a, r in zip(acquisitions, found):
        n = len(a["time"])
        ours = ~r.keep
        row = [a["name"], n, 100 * ours.mean() if r.checked else float("nan"),
               100 * r.counts.get(qc.RATE, 0) / n, 100 * r.counts.get(qc.SIGNAL, 0) / n]
        theirs = {}
        for name, which in PROGRAMS:
            theirs[name] = flags(folder, a["stem"], which, n)
            row.append(100 * theirs[name].mean() if theirs[name] is not None else float("nan"))
        for name in ("PeacoQC", "flowAI"):
            row.append(jaccard(ours, theirs[name]) if theirs[name] is not None else float("nan"))
        row.append(notes.get(a["stem"], ""))
        rows.append(row)
        if a["kind"] == "real":
            for name, _ in PROGRAMS:
                if theirs[name] is not None:
                    summary[name].append((ours.mean(), theirs[name].mean(),
                                          jaccard(ours, theirs[name]), (ours & theirs[name]).sum(),
                                          ours.sum(), theirs[name].sum(), n))
        if a["truth"] is not None:
            truth = a["truth"]
            line = [a["name"], 100 * truth.mean()]
            for who, removed in [("ours", ours)] + [(k, v) for k, v in theirs.items()]:
                if removed is None:
                    line += [float("nan"), float("nan")]
                else:
                    line += [100 * removed[truth].mean(), 100 * removed[~truth].mean()]
            truth_rows.append(line)

    header = ["acquisition", "events", "ours removed %", "ours: flow rate %", "ours: signal %",
              "PeacoQC %", "PeacoQC MAD only %", "flowAI (rate+signal) %",
              "Jaccard ours/PeacoQC", "Jaccard ours/flowAI", "PeacoQC's own account"]
    for kind, title in (("synthetic", "Synthetic acquisitions"), ("real", "Real files"),
                        ("injected", "A clean real file with a defect put in")):
        some = [row for a, row in zip(acquisitions, rows) if a["kind"] == kind]
        if some:
            rep.table(title, header, some)
    if truth_rows:
        rep.table("Where the defect is known: share of the defect's events removed, and of "
                  "the others",
                  ["acquisition", "defect (% of events)", "ours: of the defect %",
                   "ours: of the rest %", "PeacoQC: defect %", "PeacoQC: rest %",
                   "PeacoQC MAD only: defect %", "MAD only: rest %", "flowAI: defect %",
                   "flowAI: rest %"], truth_rows)
    out = []
    for name, got in summary.items():
        if not got:
            continue
        g = np.array(got, dtype=float)
        out.append([name, len(g), 100 * g[:, 4].sum() / g[:, 6].sum(),
                    100 * g[:, 5].sum() / g[:, 6].sum(), 100 * np.median(g[:, 0]),
                    100 * np.median(g[:, 1]),
                    g[:, 3].sum() / max(g[:, 4].sum() + g[:, 5].sum() - g[:, 3].sum(), 1),
                    100 * g[:, 3].sum() / max(g[:, 4].sum(), 1),
                    100 * g[:, 3].sum() / max(g[:, 5].sum(), 1),
                    int(((g[:, 0] < 0.01) & (g[:, 1] < 0.01)).sum()),
                    int(((g[:, 0] >= 0.01) & (g[:, 1] >= 0.01)).sum()),
                    int(((g[:, 0] >= 0.01) & (g[:, 1] < 0.01)).sum()),
                    int(((g[:, 0] < 0.01) & (g[:, 1] >= 0.01)).sum())])
    if out:
        rep.table("Real files together: ours against each program",
                  ["program", "files", "ours removed % (all events)", "it removed %",
                   "ours, median file %", "it, median file %", "Jaccard over all events",
                   "of ours' removed, it removes too %", "of its removed, ours removes too %",
                   "files: both under 1%", "both 1% or more", "only ours 1% or more",
                   "only it 1% or more"], out)
    rep.close()


if __name__ == "__main__":
    main()
