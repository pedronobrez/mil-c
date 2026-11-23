"""milc.proliferation.fit beside R's flowFit, on flowFit's own data.

    python3 tests/reference/proliferation_compare.py

The data are Quah & Parish's CD4 T cells labelled with CFSE, CellTrace
Violet (CTV) and Cell Proliferation Dye (CPD), stimulated for three days,
and one sample of all three left unstimulated: the `QuahAndParish` flowSet
of Bioconductor's flowFitExampleData 1.22.0. flowFit itself (1.24.0) left
Bioconductor after release 3.10 and is installed from that release's
archive; with flowCore 2.x its NAMESPACE must stop importing `Data`, which
flowCore no longer has (the one change made to it).

Both fit Gaussians of one width, evenly spaced on the logarithm of the
dye. They differ in two things: flowFit also fits the distance between
generations, where ours is log10(2) by definition (or, with an
autofluorescence, what halving above it gives); and flowFit starts from
the control's peak and lets it go, where ours keeps it within a third of
a division. flowFit's "proliferation index" is FlowJo's Expansion Index.
"""

from __future__ import annotations

import csv
import os

import numpy as np

import common
from milc import proliferation

DYES = {"<FITC-A>": "CFSE", "<APC-A>": "CPD", "<Alexa Fluor 405-A>": "CTV"}
SHOWN = 9          # generations 0 to 8


def events(folder: str, k: int) -> tuple[list[str], np.ndarray]:
    path = os.path.join(folder, f"quah_{k}.csv")
    with open(path, encoding="utf-8") as f:
        names = next(csv.reader(f))
    return names, np.loadtxt(path, delimiter=",", skiprows=1)


def main() -> None:
    folder = common.work("proliferation")
    print(common.rscript("proliferation_flowfit.R", folder).strip().splitlines()[-1])
    with open(os.path.join(folder, "flowfit.csv"), encoding="utf-8") as f:
        theirs = list(csv.DictReader(f))
    names, control = events(folder, 1)
    for row in theirs:
        channel, k = row["channel"], int(row["sample"])
        j = names.index(channel)
        _, data = events(folder, k)
        # the package stores log10 of the fluorescence; 0 is an event with none
        dye = 10.0 ** data[:, j][data[:, j] > 0]
        peak = proliferation.undivided(10.0 ** control[:, j][control[:, j] > 0])
        percent = np.array([float(v) for v in row["percent"].split()])[:SHOWN]
        percent = 100 * percent / percent.sum()
        print(f"\n{DYES[channel]} ({channel}), {dye.size:,d} events; control's peak at "
              f"{peak.position:.3f} ± {peak.sd:.3f} decades")
        print(f"  flowFit: peak {float(row['position']):.3f}, width {float(row['size']):.3f}, "
              f"distance {float(row['distance']):.3f} (log10(2) is 0.301); its index "
              f"{float(row['index']):.3f}")
        print("  flowFit   % in generations 0-8:", np.round(percent, 1))
        for label, kw in (("control", dict(control=peak)), ("no control", {}),
                          ("control, autofluorescence 150", dict(control=peak, background=150.0))):
            f = proliferation.fit(dye, max_generations=10, **kw)
            said = ", ".join(f"{name} {value:.3f} ± {se:.3f}"
                             for name, (value, se) in f.indices.items())
            print(f"  MIL-C, {label}: peak {f.position:.3f}, width {f.sd:.3f}, RMS {f.rms:.1f} "
                  f"({f.rms_percent:.1f}% of the tallest bin)")
            print("            % in generations 0-8:", np.round(f.shares[:SHOWN], 1),
                  f" furthest from flowFit: {np.abs(f.shares[:SHOWN] - percent).max():.1f} points")
            print(f"            {said}")


if __name__ == "__main__":
    main()
