"""milc.qc.check on known-bad copies of the real files: how much of a defect
of known size it finds, at which limits, and what it takes from the files
as they are.

    python3 tests/reference/qc_calibrate.py

No R. Every defect is written into a copy of a real file, in memory, on the
values `qc.check` judges (compensated, each parameter on its display
scale), as a move of `size` times the spread of that parameter's events
(the MAD, scaled to a standard deviation) over a stretch of the events in
the order of Time:

- short: every FSC and SSC moved over 1, 2, 5 or 10% of the events, at a
  place (fixed per file and duration) the unchanged file keeps whole
  under every limit tried — the shift limit's table;
- lasting: every stain over the first 30%; PE-A alone over the last 25%;
  every stain, a ramp from nothing to `size` over the last 30%; every
  FSC and SSC over 40-60% — the drift limit's table, on the files the
  default check leaves almost whole (under 1%).

Found: half of the defect's events or more removed. "Rest": what is
removed of the other events and not of the unchanged file.
"""

from __future__ import annotations

import os
import zlib
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import common
import qc_compare
from milc import qc

SHORT_SIZES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.75]
SHORT_DURATIONS = [0.01, 0.02, 0.05, 0.10]
SHIFT_LIMITS = [0.0, 0.2, 0.25, 0.3]
DRIFT_SIZES = [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
DRIFT_LIMITS = [0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.0]
KINDS = ["start/stains", "end/PE-A", "ramp/stains", "middle/scatter"]


def _judged(a: dict) -> np.ndarray:
    return a["compensate"](a["raw"]) if a["compensate"] else a["raw"]


def _check(a: dict, x: np.ndarray, **kw) -> qc.Result:
    return qc.check(x, a["time"], a["channels"], a["scales"], **kw)


def move(a: dict, x: np.ndarray, rows: np.ndarray, size: float, which: str,
         amount: np.ndarray | None = None) -> np.ndarray:
    """A copy of `x` with the parameters `which` names ("scatter", "stains"
    or one parameter) moved over `rows` by `size` spreads, times `amount`."""
    x = x.copy()
    for j, (c, sc) in enumerate(zip(a["channels"], a["scales"])):
        wanted = {"scatter": c.startswith("FSC") or c.startswith("SSC"),
                  "stains": j in a["fluorescence"]}.get(which, c == which)
        if not wanted:
            continue
        v = sc(x[:, j]) if sc is not None else x[:, j]
        by = size * qc._spread(v) * (1.0 if amount is None else amount)
        moved = v[rows] + by
        x[rows, j] = sc.inverse(moved) if sc is not None else moved
    return x


def short_power(a: dict) -> list[dict]:
    x = _judged(a)
    n = len(x)
    order = np.argsort(a["time"], kind="stable")
    rng = np.random.default_rng(zlib.crc32(a["name"].encode()))
    base = {lim: _check(a, x, shift=lim).keep for lim in SHIFT_LIMITS}
    whole = np.logical_and.reduce(list(base.values()))[order]
    out = []
    for duration in SHORT_DURATIONS:
        length = int(round(duration * n))
        free = np.convolve(~whole, np.ones(length), "valid") == 0
        lo, hi = int(0.02 * n), int(0.98 * n) - length
        starts = np.flatnonzero(free[lo:hi]) + lo
        if not len(starts):
            continue
        rows = order[int(rng.choice(starts)):][:length]
        hit = np.zeros(n, bool)
        hit[rows] = True
        for size in SHORT_SIZES:
            moved = move(a, x, rows, size, "scatter")
            for lim in SHIFT_LIMITS:
                r = _check(a, moved, shift=lim)
                gone = ~r.keep
                out.append({"file": a["name"], "duration": duration, "size": size, "limit": lim,
                            "defect": float(gone[hit].mean()),
                            "rest": float((gone & base[lim])[~hit].mean())})
    return out


def lasting(a: dict, x: np.ndarray, kind: str, size: float):
    n = len(x)
    order = np.argsort(a["time"], kind="stable")
    share = np.empty(n)
    share[order] = np.arange(n) / n
    where, which = kind.split("/")
    amount = None
    if where == "start":
        hit = share < 0.30
    elif where == "end":
        hit = share >= 0.75
    elif where == "ramp":
        hit = share >= 0.70
        amount = np.clip((share[hit] - 0.70) / 0.30, 0, 1)
    else:
        hit = (share >= 0.40) & (share < 0.60)
    return move(a, x, np.flatnonzero(hit), size, which, amount), hit


def drift_power(a: dict) -> list[dict]:
    x = _judged(a)
    runs = [(lim, {"drift": lim}) for lim in DRIFT_LIMITS] + [("drift only", {"shift": 0.0})]
    base = {lim: _check(a, x, **kw).keep for lim, kw in runs}
    out = []
    for kind in KINDS:
        for size in DRIFT_SIZES:
            moved, hit = lasting(a, x, kind, size)
            for lim, kw in runs:
                gone = ~_check(a, moved, **kw).keep
                out.append({"file": a["name"], "kind": kind, "size": size, "limit": lim,
                            "defect": float(gone[hit].mean()),
                            "rest": float((gone & base[lim])[~hit].mean())})
    return out


def as_is(a: dict) -> list[dict]:
    x = _judged(a)
    out = []
    for drift in DRIFT_LIMITS:
        for shift in SHIFT_LIMITS:
            if drift != qc.DRIFT_LIMIT and shift not in (0.0, qc.SHIFT_LIMIT):
                continue
            r = _check(a, x, drift=drift, shift=shift)
            out.append({"file": a["name"], "n": len(x), "drift": drift, "shift": shift,
                        "removed": int((~r.keep).sum())})
    return out


def _all(a: dict) -> tuple:
    r = as_is(a)
    base = [v for v in r if v["drift"] == qc.DRIFT_LIMIT and v["shift"] == qc.SHIFT_LIMIT][0]
    lasting_rows = drift_power(a) if base["removed"] < 0.01 * base["n"] else []
    return r, short_power(a), lasting_rows


def _cell(values: list[float]) -> str:
    v = np.array(values)
    return f"{100 * v.mean():.0f} ({int((v >= 0.5).sum())}/{len(v)})"


def main() -> None:
    rep = common.Report("qc_calibrate")
    rep.say("# qc.check on known-bad copies of the real files\n")
    if not common.has_real():
        rep.say("no real data (MILC_REF_DATA not set): nothing to do\n")
        rep.close()
        return
    files = qc_compare.real()
    for a in files:     # compensated once, here: a compensation does not go to another process
        a["raw"], a["compensate"] = _judged(a), None
    with ProcessPoolExecutor(int(os.environ.get("MILC_REF_WORKERS", "4"))) as ex:
        done = list(ex.map(_all, files))
    plain = [v for d in done for v in d[0]]
    short = [v for d in done for v in d[1]]
    drift = [v for d in done for v in d[2]]
    total = sum(v["n"] for v in plain if v["drift"] == qc.DRIFT_LIMIT and v["shift"] == 0.0)

    rows = []
    for d, s in sorted({(v["drift"], v["shift"]) for v in plain}):
        these = [v for v in plain if (v["drift"], v["shift"]) == (d, s)]
        sel = {v["file"]: v["removed"] / v["n"] for v in these}
        removed = sum(v["removed"] for v in these)
        rows.append([d, s, 100 * removed / total, sum(x >= 0.01 for x in sel.values()),
                     sum(x >= 0.05 for x in sel.values()), sum(x >= 0.2 for x in sel.values())]
                    + [100 * sel.get(f, float("nan"))
                       for f in ("7Ctrl5", "14Case1", "14Case3", "14Case4`", "7Case5")])
    rep.table(f"The {len(files)} real files as they are",
              ["drift limit", "shift limit", "removed % (all events)", "files with 1% or more",
               "5% or more", "20% or more", "7Ctrl5 %", "14Case1 %", "14Case3 %", "14Case4` %",
               "7Case5 %"], rows)

    for lim in SHIFT_LIMITS:
        rows = []
        for duration in SHORT_DURATIONS:
            sel = [v for v in short if v["limit"] == lim and v["duration"] == duration]
            rows.append([f"{100 * duration:g}%"]
                        + [_cell([v["defect"] for v in sel if v["size"] == s]) for s in SHORT_SIZES]
                        + [100 * np.mean([v["rest"] for v in sel])])
        rep.table(f"A short shift of every FSC and SSC, shift limit {lim:g}: % of the defect "
                  "removed (files with half of it or more)",
                  ["duration"] + [f"size {s:g}" for s in SHORT_SIZES] + ["rest %"], rows)

    bases = sorted({v["file"] for v in drift})
    rows = []
    for lim in DRIFT_LIMITS + ["drift only"]:
        sel = [v for v in drift if v["limit"] == lim]
        rows.append([lim] + [_cell([v["defect"] for v in sel if v["size"] == s])
                             for s in DRIFT_SIZES]
                    + [100 * np.mean([v["rest"] for v in sel])])
    rep.table(f"A lasting shift, four kinds pooled, on {len(bases)} files ({', '.join(bases)}): "
              "% of the defect removed (cases with half of it or more)",
              ["drift limit"] + [f"size {s:g}" for s in DRIFT_SIZES] + ["rest %"], rows)
    for kind in KINDS:
        rows = [[lim] + [100 * np.mean([v["defect"] for v in drift if v["limit"] == lim
                                        and v["kind"] == kind and v["size"] == s])
                         for s in DRIFT_SIZES]
                for lim in DRIFT_LIMITS + ["drift only"]]
        rep.table(f"Lasting, {kind}: % of the defect removed",
                  ["drift limit"] + [f"size {s:g}" for s in DRIFT_SIZES], rows)
    rep.close()


if __name__ == "__main__":
    main()
