"""milc.cellcycle.fit on PI alone, beside the phases BrdU gives in the
same cells.

    python3 -I tests/reference/cellcycle_compare.py /path/to/unzipped/folder

The data: Yang & Hill, "The ORFIUS complex regulates ORC2 localization at
replication origins" — 223 FCS files of four ovarian cancer lines (OV8,
FT194, UWB1.289 with and without BRCA1, CaOV3, SKOV3) pulsed with BrdU,
stained with anti-BrdU (FITC) and propidium iodide, on a BD Fortessa;
Dryad doi:10.5061/dryad.sn02v6x89 (CC0), also on Zenodo, record 10499157.
Unzip `Flow_cytometry_data_for_Yang_and_Hill_ORFIUS_complex_manuscript.zip`
into a folder of its own and give that folder.

No program that fits DNA histograms the way ModFit or FlowJo do could be
run here, so the reference is the BrdU in the same file, read the way a
bivariate BrdU/PI plot is gated (and written here, not by the authors):

  single cells   PI-H / PI-A above 0.82 of its median at the G1 peak
  BrdU+          log10 FITC-A more than 0.4 above the BrdU- band, which is
                 drawn through the mode of the G1 peak and that of the G2
                 peak; BrdU+ = S
  G1, G2         BrdU- single cells below and above the midpoint of the
                 two peaks, between 0.8 x G1 and 1.2 x G2

This is not the truth for DNA content: cells that took up BrdU and went on
to G2 (or divided) during the pulse count as S, and early-S cells with
little BrdU as G1. Where a drug stops replication (XL413, AZD6738,
olaparib), cells sit in S without taking up BrdU and the reference breaks;
the summary is therefore also given for the control samples only.
"""

from __future__ import annotations

import csv
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # run with -I too
import common

from milc import cellcycle as cc
from milc.fcs import read_fcs

DELTA = 0.4
E = cc.Extras
FITS = {
    "singlets, Watson": (True, cc.WATSON, None),
    "singlets, DJF": (True, cc.DJF, None),
    "singlets, DJF + debris": (True, cc.DJF, E(debris=cc.SINGLE_CUT)),
    "all events, DJF": (False, cc.DJF, None),
    "all events, DJF + aggregates": (False, cc.DJF, E(aggregates=True)),
    "all events, DJF + debris + aggregates": (False, cc.DJF,
                                              E(debris=cc.SINGLE_CUT, aggregates=True)),
}
CONTROL = re.compile(r"dmso|sicon|controi|shcon", re.IGNORECASE)


def _mode(v, lo, hi):
    h, e = np.histogram(v, bins=60, range=(lo, hi))
    k = np.convolve(h, np.ones(3) / 3, mode="same").argmax()
    return (e[k] + e[k + 1]) / 2


def brdu(pa, ph, br):
    """(single-cell mask, % G1 / S / G2 by BrdU, G1, G2)."""
    ok = (pa > 0) & (ph > 0) & (pa < pa.max())
    m1, _ = cc._find_g1(pa[ok], 256)
    near = ok & (np.abs(pa / m1 - 1) < 0.06)
    ratio = ph / np.maximum(pa, 1)
    single = ok & (ratio > 0.82 * np.median(ratio[near]))
    lb = np.log10(np.clip(br, 1, None))
    g1n = _mode(lb[single & near], 1.5, 5.0)
    m2 = 2 * m1
    for _ in range(2):
        m2 = float(np.median(pa[single & (np.abs(pa / m2 - 1) < 0.08)]))
    g2n = _mode(lb[single & (np.abs(pa / m2 - 1) < 0.05)], 1.5, 5.0)
    positive = lb > g1n + (g2n - g1n) * (pa - m1) / (m2 - m1) + DELTA
    inside = single & (pa > 0.8 * m1) & (pa < 1.2 * m2)
    mid = (m1 + m2) / 2
    n = np.array([(inside & ~positive & (pa < mid)).sum(), (inside & positive).sum(),
                  (inside & ~positive & (pa >= mid)).sum()], dtype=float)
    return single, 100 * n / n.sum(), m1, m2


def main(folder: str) -> None:
    out = common.work("cellcycle")
    rows = []
    for path in sorted(glob.glob(os.path.join(folder, "**", "*.fcs"), recursive=True)):
        f = read_fcs(path)
        col = {c: f.events[:, f.index(c)].astype(float) for c in f.channels}
        pa, ph, br = col["PI-A"], col["PI-H"], col["FITC-A"]
        single, truth, m1, _ = brdu(pa, ph, br)
        row = {"file": os.path.relpath(path, folder),
               "control": bool(CONTROL.search(path)) and "dox" not in path.lower()}
        row.update({f"{p} BrdU": truth[k] for k, p in enumerate(cc.PHASES)})
        for name, (singles, kind, extras) in FITS.items():
            fit = cc.fit(pa[single] if singles else pa, kind, extras=extras)
            row.update({f"{p} {name}": fit.percent[k] for k, p in enumerate(cc.PHASES)})
            row[f"chi2 {name}"] = fit.chi2
            if extras is not None and extras.aggregates:
                # the aggregates found against what the singlet gate drops, over
                # the same range of DNA
                lo, hi = fit.fitted
                lo = max(lo, 0.75 * m1)
                span = (pa >= lo) & (pa <= hi) & (pa < pa.max())
                keep = (fit.x >= lo) & (fit.x <= hi)
                agg = fit.components[fit.names.index("Aggregates")][keep].sum()
                row[f"aggregates {name}"] = 100 * agg / fit.model[keep].sum()
                row["dropped by the singlet gate"] = 100 * (span & ~single).sum() / span.sum()
        rows.append(row)
    with open(os.path.join(out, "cellcycle.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    lines = [f"{len(rows)} files; model minus BrdU, in points: mean ± SD over files", ""]
    for title, chosen in (("all files", rows), ("controls", [r for r in rows if r["control"]])):
        lines.append(f"{title} ({len(chosen)}); BrdU: " + ", ".join(
            f"{p} {np.mean([r[f'{p} BrdU'] for r in chosen]):.1f}" for p in cc.PHASES))
        for name in FITS:
            d = {p: np.array([r[f"{p} {name}"] - r[f"{p} BrdU"] for r in chosen])
                 for p in cc.PHASES}
            chi2 = np.median([r[f"chi2 {name}"] for r in chosen])
            lines.append(f"  {name:40s} " + "  ".join(
                f"{p} {d[p].mean():+5.1f} ± {d[p].std(ddof=1):4.1f}" for p in cc.PHASES)
                + f"  median χ² {chi2:.2f}")
        lines.append("")
    name = "all events, DJF + debris + aggregates"
    a = np.array([(r[f"aggregates {name}"], r["dropped by the singlet gate"]) for r in rows])
    lines.append(f"aggregates found {a[:, 0].mean():.1f}% (SD {a[:, 0].std():.1f}), dropped by "
                 f"the singlet gate {a[:, 1].mean():.1f}% (SD {a[:, 1].std():.1f}); "
                 f"r = {np.corrcoef(a.T)[0, 1]:.2f}")
    for p in ("S", "G2"):
        d = np.array([r[f"{p} {name}"] - r[f"{p} singlets, DJF"] for r in rows])
        plain = np.array([r[f"{p} all events, DJF"] - r[f"{p} singlets, DJF"] for r in rows])
        lines.append(f"{p}: all events with aggregates minus singlets {d.mean():+.2f} ± "
                     f"{d.std():.2f}; all events without them {plain.mean():+.2f} ± "
                     f"{plain.std():.2f}")
    report = "\n".join(lines)
    with open(os.path.join(out, "report.md"), "w", encoding="utf-8") as f:
        f.write(report + "\n")
    print(report)


if __name__ == "__main__":
    main(sys.argv[1])
