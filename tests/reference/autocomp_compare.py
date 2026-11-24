"""milc.autocomp beside AutoSpill (carlosproca/autospill) and
flowStats::spillover, on the same single-stain controls.

    python3 tests/reference/autocomp_compare.py

Three sets of controls:

- `sim6`: six detectors, simulated with a known spillover matrix
  (`autocomp.simulate`), with 15% debris and 1% stray bright events, and
  an unstained control;
- `sim12`: twelve detectors, a denser matrix, clean;
- `compdata`: the real single-stain controls flowCore ships
  (`extdata/compdata`, a FACSCalibur, four detectors and an unstained
  tube), read by MIL-C — which applies their log amplification — and
  written again as linear values so that both programs read the same
  numbers.

On each, MIL-C's two methods and AutoSpill's two steps are run on every
event (the algorithms alone, the same events), on the events AutoSpill's
own scatter gate keeps (its gate written out by the R side and given to
MIL-C), and MIL-C on its own scatter gate. Every matrix is judged by its
distance to the truth where there is one, and everywhere by the spill it
leaves in the controls once they are compensated with it (the robust
slope of every detector on the control's own, linear scale).
"""

from __future__ import annotations

import json
import os
import shutil
import time

import numpy as np

import common
from milc import autocomp
from milc.fcs import read_fcs, write_fcs

R_LIB = os.environ.get("R_LIBS_USER", "")
COMPDATA = os.path.join(R_LIB, "flowCore", "extdata", "compdata", "data")


def sim6():
    det = ["FITC-A", "PE-A", "PerCP-Cy5.5-A", "PE-Cy7-A", "APC-A", "BV421-A"]
    S = np.eye(6)
    for (i, j), v in {(0, 1): 0.22, (0, 2): 0.05, (1, 0): 0.02, (1, 2): 0.30, (1, 3): 0.04,
                      (2, 3): 0.12, (2, 4): 0.08, (3, 2): 0.03, (4, 2): 0.15, (4, 3): 0.02,
                      (5, 0): 0.03, (1, 4): 0.01, (2, 1): 0.01}.items():
        S[i, j] = v
    return det, S, autocomp.simulate(S, det, n=20000, seed=6, debris=0.15, outliers=0.01)


def sim12():
    det = ["BUV395-A", "BUV737-A", "BV421-A", "BV510-A", "BV605-A", "BV786-A", "FITC-A",
           "PerCP-Cy5.5-A", "PE-A", "PE-Cy7-A", "APC-A", "APC-Cy7-A"]
    rng = np.random.default_rng(12)
    S = np.eye(12)
    for i in range(12):
        for j in range(12):
            if i != j and abs(i - j) <= 3:
                S[i, j] = rng.uniform(0.0, 0.35) / abs(i - j)
    return det, S, autocomp.simulate(S, det, n=15000, seed=12)


def compdata():
    names = {"060909.001": None, "060909.002": "FL1-H", "060909.003": "FL2-H",
             "060909.004": "FL4-H", "060909.005": "FL3-H"}          # flowCore's comp_match
    keep = ["FSC-H", "SSC-H", "FL1-H", "FL2-H", "FL3-H", "FL4-H"]
    out = []
    for file, dye in names.items():
        f = read_fcs(os.path.join(COMPDATA, file))
        ev = f.events[:, [f.index(c) for c in keep]].astype(np.float64)
        out.append((f"{file.replace('.', '_')}_{dye or 'unstained'}", ev, keep, dye))
    return keep[2:], None, out


def write_set(name, detectors, controls, scatter, folder):
    fcs = os.path.join(folder, "fcs")
    shutil.rmtree(fcs, ignore_errors=True)
    os.makedirs(fcs)
    files = []
    for item in controls:
        if len(item) == 3:
            cname, ev, ch = item
            dye = None if cname.startswith("Unstained") else detectors[
                [d.rsplit("-", 1)[0] for d in detectors].index(cname.replace(" stained control", ""))]
        else:
            cname, ev, ch, dye = item
        path = os.path.join(fcs, cname.replace(" ", "_") + ".fcs")
        write_fcs(path, ev, ch, ranges=[262144.0] * len(ch))
        files.append({"file": os.path.basename(path), "dye": dye})
    return {"name": name, "files": files, "scatter": scatter, "detectors": detectors}


def read_r(path, detectors):
    import csv
    with open(path, encoding="utf-8") as f:
        rows = list(csv.reader(f))
    cols = rows[0][1:]
    m = {r[0]: [float(v) for v in r[1:]] for r in rows[1:]}
    out = np.eye(len(detectors))
    for i, a in enumerate(detectors):
        if a in m:
            for j, b in enumerate(detectors):
                out[i, j] = m[a][cols.index(b)]
    return out


def ours(folder, d, method, gate=False, rows=None, **kw):
    paths, primaries = [], {}
    for f in d["files"]:
        p = os.path.join(folder, "fcs", f["file"])
        paths.append(p)
        primaries[p] = f["dye"]
    controls = autocomp.load_controls(paths, d["detectors"], gate=gate, max_events=None)
    if rows is not None:
        for c, f in zip(controls, d["files"]):
            if f["dye"] is not None:
                c.events = c.events[rows[f["dye"]]]
    t = time.time()
    r = autocomp.compute(paths, d["detectors"], method, primaries=primaries, controls=controls,
                         **kw)
    return r, time.time() - t, controls


def offdiag(a, b):
    m = ~np.eye(len(a), dtype=bool)
    return float(np.abs(a - b)[m].max())


def main() -> None:
    folder = common.work("autocomp")
    rep = common.Report("autocomp")
    sets = []
    for name, make, scatter in (("sim6", sim6, ["FSC-A", "SSC-A"]),
                                ("sim12", sim12, ["FSC-A", "SSC-A"]),
                                ("compdata", compdata, ["FSC-H", "SSC-H"])):
        if name == "compdata" and not os.path.isdir(COMPDATA):
            rep.say("flowCore's compdata not found (R_LIBS_USER): left out")
            continue
        det, truth, controls = make()
        os.makedirs(os.path.join(folder, name), exist_ok=True)
        d = write_set(name, det, controls, scatter, os.path.join(folder, name))
        d["truth"] = truth
        sets.append(d)
    with open(os.path.join(folder, "manifest.json"), "w") as f:
        json.dump({"datasets": [{k: v for k, v in d.items() if k != "truth"} for d in sets]}, f)
    t = time.time()
    printed = common.rscript("autocomp_autospill.R", folder)
    rep.say(f"R side: {time.time() - t:.0f} s\n```\n{printed.strip()}\n```")

    for d in sets:
        here = os.path.join(folder, d["name"])
        det, truth = d["detectors"], d["truth"]
        rep.say(f"\n### {d['name']}: {len(det)} detectors, "
                f"{sum(f['dye'] is not None for f in d['files'])} stained controls"
                + (" and an unstained one" if any(f["dye"] is None for f in d["files"]) else ""))
        mats, seconds, iters = {}, {}, {}
        for method in (autocomp.CLASSIC, autocomp.AUTOSPILL):
            r, s, controls_all = ours(here, d, method)
            mats[f"MIL-C {method}, all events"] = r.matrix
            seconds[f"MIL-C {method}, all events"] = s
            iters[f"MIL-C {method}, all events"] = len(r.history) or None
            r, s, _ = ours(here, d, method, gate=True)
            mats[f"MIL-C {method}, MIL-C's scatter gate"] = r.matrix
            seconds[f"MIL-C {method}, MIL-C's scatter gate"] = s
            iters[f"MIL-C {method}, MIL-C's scatter gate"] = len(r.history) or None
        if any(f["dye"] is None for f in d["files"]):
            r, s, _ = ours(here, d, autocomp.CLASSIC, universal_negative=True)
            key = "MIL-C classic, unstained as negative, all events"
            mats[key], seconds[key] = r.matrix, s
        for run in ("all", "own"):
            label = "all events" if run == "all" else "AutoSpill's scatter gate"
            for step, file in (("first step", "initial.csv"), ("refined", "final.csv")):
                p = os.path.join(here, run, file)
                if os.path.exists(p):
                    mats[f"AutoSpill {step}, {label}"] = read_r(p, det)
            p = os.path.join(here, run, "convergence.csv")
            if os.path.exists(p):
                with open(p) as f:
                    iters[f"AutoSpill refined, {label}"] = sum(1 for _ in f) - 1
                with open(os.path.join(here, run, "seconds.txt")) as f:
                    seconds[f"AutoSpill refined, {label}"] = float(f.read())
        if os.path.exists(os.path.join(here, "posnegpop_all.csv")):
            mats["AutoSpill posnegpop, all events"] = read_r(
                os.path.join(here, "posnegpop_all.csv"), det)
        if os.path.exists(os.path.join(here, "flowstats.csv")):
            mats["flowStats::spillover, all events"] = read_r(
                os.path.join(here, "flowstats.csv"), det)
        # MIL-C on AutoSpill's own gate: the same events on both sides
        gates = {}
        for f in d["files"]:
            p = os.path.join(here, f"gate_{f['dye']}.csv")
            if f["dye"] is not None and os.path.exists(p):
                gates[f["dye"]] = np.loadtxt(p, skiprows=1, dtype=int).reshape(-1) - 1
        if len(gates) == sum(f["dye"] is not None for f in d["files"]):
            for method in (autocomp.CLASSIC, autocomp.AUTOSPILL):
                r, s, _ = ours(here, d, method, rows=gates)
                key = f"MIL-C {method}, AutoSpill's scatter gate"
                mats[key], seconds[key], iters[key] = r.matrix, s, len(r.history) or None
        # every matrix judged on the same events: every event of every control
        judge = controls_all
        for c, f in zip(judge, d["files"]):
            c.primary, c.unstained = f["dye"], f["dye"] is None
        rows = []
        for k, m in mats.items():
            left, _ = autocomp.residual_spill(judge, det, m)
            rows.append([k, offdiag(m, truth) if truth is not None else float("nan"),
                         float(np.abs(left).max()), float(np.sqrt((left ** 2).mean())),
                         iters.get(k) or "", f"{seconds[k]:.1f}" if k in seconds else ""])
        rep.table(f"{d['name']}: each matrix", ["matrix", "largest error against the truth",
                                                "largest spill left", "rms spill left",
                                                "iterations", "seconds"], rows)
        pairs = [("MIL-C autospill, all events", "AutoSpill refined, all events"),
                 ("MIL-C autospill, AutoSpill's scatter gate", "AutoSpill refined, AutoSpill's scatter gate"),
                 ("MIL-C classic, all events", "AutoSpill posnegpop, all events"),
                 ("MIL-C classic, unstained as negative, all events",
                  "flowStats::spillover, all events"),
                 ("AutoSpill first step, all events", "AutoSpill refined, all events")]
        rows = [[a, b, offdiag(mats[a], mats[b])] for a, b in pairs if a in mats and b in mats]
        rep.table(f"{d['name']}: between programs", ["", "", "largest difference"], rows)
        if d["name"] == "compdata":
            names = list(mats)
            rep.say("\nThe matrices (rows: fluorochrome; columns: detector):\n")
            for k in names:
                rep.say(f"{k}:\n```\n{np.array2string(mats[k], precision=4, suppress_small=True)}\n```")
    rep.close()


if __name__ == "__main__":
    main()
