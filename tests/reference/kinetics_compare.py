"""milc.kinetics beside the same definitions written again in R
(`kinetics_flowcore.R`, the files read by flowCore).

    python3 -I tests/reference/kinetics_compare.py /path/to/folder/of/fcs [every]

No public calcium flux could be reached when this was written (FlowRepository
did not answer; Zenodo had none), so this checks the arithmetic on real
files with a Time parameter — the BrdU/PI files of `cellcycle_compare.py`,
every `every`-th of them (10 by default) — following FITC-A and the ratio
FITC-A / PI-A. It says whether the slices, medians, smoothing and the
measures are computed as defined, not whether the definitions suit a flux.
"""

from __future__ import annotations

import csv
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # run with -I too
import common

from milc import kinetics
from milc.fcs import read_fcs
from milc.qc_plots import timestep


def main(folder: str, every: int = 10) -> None:
    out = common.work("kinetics")
    files = sorted(glob.glob(os.path.join(folder, "**", "*.fcs"), recursive=True))[::every]
    listing = os.path.join(out, "files.txt")
    with open(listing, "w", encoding="utf-8") as f:
        f.write("\n".join(files) + "\n")
    theirs_path = os.path.join(out, "flowcore.csv")
    print(common.rscript("kinetics_flowcore.R", listing, theirs_path).strip())
    with open(theirs_path, encoding="utf-8") as f:
        theirs = {(r["file"], r["what"]): r for r in csv.DictReader(f)}
    worst: dict = {}
    for path in files:
        fcs = read_fcs(path)
        col = {c: fcs.events[:, fcs.index(c)].astype(np.float64) for c in fcs.channels}
        for what in ("FITC-A", "ratio"):
            c = kinetics.curve(col["Time"], col["FITC-A"], step=timestep(fcs.keywords),
                               denominator=col["PI-A"] if what == "ratio" else None)
            half = (c.t[1] - c.t[0]) / 2
            m = kinetics.measure(c, *kinetics.default_marks(c.t[0] - half, c.t[-1] + half))
            ref = theirs[(path, what)]
            pairs = [(m["Baseline"], ref["baseline"]), (m["Peak"], ref["peak"]),
                     (m["Time of peak"], ref["time_of_peak"]), (m["AUC"], ref["auc"]),
                     (c.y[49], ref["median_slice_50"])]
            for name, (a, b) in zip(("Baseline", "Peak", "Time of peak", "AUC",
                                     "median of slice 50"), pairs, strict=True):
                rel = abs(a - float(b)) / max(abs(float(b)), 1e-12)
                worst[(what, name)] = max(worst.get((what, name), 0.0), rel)
    lines = [f"{len(files)} files; largest relative difference, ours against R:"]
    lines += [f"  {what:7s} {name:20s} {v:.1e}" for (what, name), v in worst.items()]
    report = "\n".join(lines)
    with open(os.path.join(out, "report.md"), "w", encoding="utf-8") as f:
        f.write(report + "\n")
    print(report)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 10)
