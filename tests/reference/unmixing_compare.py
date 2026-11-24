"""`milc.unmixing` beside SpectroFlo and AutoSpectral.

The data: "AutoSpectral: Bead unmixed" (Burton, Dooley, Liston; Mendeley
Data, doi 10.17632/x7z4bdvv7y.1, CC BY 4.0): a 5-laser Cytek Aurora (64
detectors), 42 single-stain bead controls, negative beads and unstained
spleen cells, each as acquired (`Raw`) and as SpectroFlo 3.3.0 unmixed it
(`Unmixed`), event for event. Only the bead controls are needed (about
260 MB of the 7.4 GB archive); `--fetch` takes them out of the archive
with HTTP range requests, without downloading the rest.

1. SpectroFlo. Its unmixed values are a fixed linear map of the raw ones,
   U (64 x 43), found here by least squares; since its fit is ordinary
   least squares, its reference spectra are pinv(U). They are set beside
   the spectra `read_references` takes from the same raw controls, and
   the unmixed values of both are compared on the same events.
2. AutoSpectral's own R functions (OLS, per-event WLS, cosine similarity,
   condition number) on the same spectra and events as ours.

    export MILC_REF_SPECTRAL=/some/empty/folder       # the controls go here
    export AUTOSPECTRAL_R=/path/to/AutoSpectral/R     # optional: part 2
    python3 tests/reference/unmixing_compare.py --fetch   # once
    python3 tests/reference/unmixing_compare.py
"""

from __future__ import annotations

import glob
import io
import json
import os
import re
import subprocess
import sys
import zipfile

import numpy as np

from common import rscript, work  # noqa: E402  (also puts `milc` on the path)

from milc import unmixing as um  # noqa: E402
from milc.fcs import read_fcs  # noqa: E402

DATA = os.environ.get("MILC_REF_SPECTRAL", "")
AUTOSPECTRAL = os.environ.get("AUTOSPECTRAL_R", "")
ARCHIVE = ("https://data.mendeley.com/public-files/datasets/x7z4bdvv7y/files/"
           "a3057ddd-564d-49ef-86bf-a1077bf961ba/file_downloaded")
#: file names that are not the fluorochrome's channel name
ALIAS = {"RB780": "Real Blue 780"}


class _Remote(io.RawIOBase):
    """A file on a server that answers range requests, read with curl."""

    def __init__(self, url):
        head = subprocess.run(["curl", "-sIL", url], capture_output=True, text=True).stdout
        self.url = [ln.split(" ", 1)[1].strip() for ln in head.splitlines()
                    if ln.lower().startswith("location:")][-1] if "ocation:" in head else url
        self.size = int([ln for ln in head.splitlines()
                         if ln.lower().startswith("content-length")][-1].split(":")[1])
        self.pos = 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        data = subprocess.run(["curl", "-s", "-r", f"{self.pos}-{self.pos + n - 1}", self.url],
                              capture_output=True, check=True).stdout
        b[:len(data)] = data
        self.pos += len(data)
        return len(data)


def fetch(folder: str) -> None:
    z = zipfile.ZipFile(io.BufferedReader(_Remote(ARCHIVE), buffer_size=1 << 20))
    for kind in ("Raw", "Unmixed"):
        os.makedirs(os.path.join(folder, kind), exist_ok=True)
    for info in z.infolist():
        parts = info.filename.replace("\\", "/").split("/")
        if len(parts) < 2 or parts[1] != "Bead controls" or not parts[-1].endswith(".fcs"):
            continue
        dest = os.path.join(folder, parts[0], parts[-1])
        if not os.path.exists(dest):
            with z.open(info) as src, open(dest, "wb") as out:
                out.write(src.read())
        print("have", dest)


def _files(kind: str) -> dict[str, str]:
    """Well → file: `A3` → `.../A3 BUV395 (Beads)_Bead controls_1.fcs`."""
    out = {}
    for p in glob.glob(os.path.join(DATA, kind, "*.fcs")) + glob.glob(
            os.path.join(DATA, f"{kind}*.fcs")):
        name = p.replace("\\", "/").split("/")[-1]
        out[name.split(" ")[0]] = p
    return out


def _title(path: str) -> tuple[str, str]:
    name = path.replace("\\", "/").split("/")[-1]
    m = re.match(r"^\S+ (.*) \((Beads|Cells)\)", name)
    return m.group(1), m.group(2)


def cos(a, b) -> float:
    return float(a @ b / np.linalg.norm(a) / np.linalg.norm(b))


def spectroflo(report: list[str]) -> dict:
    raw, unmixed = _files("Raw"), _files("Unmixed")
    wells = sorted(set(raw) & set(unmixed))
    first = read_fcs(raw[wells[0]])
    det = um.detector_channels(first.channels)
    maps, fluors = {}, None
    for w in wells:
        fr, fu = read_fcs(raw[w]), read_fcs(unmixed[w])
        fluors = [c for c in fu.channels if c not in fr.channels]
        Y = fr.events[:, [fr.index(c) for c in det]].astype(float)
        Z = fu.events[:, [fu.index(c) for c in fluors]].astype(float)
        U = np.linalg.lstsq(Y, Z, rcond=None)[0]
        fit = float(np.abs(Y @ U - Z).max() / np.abs(Z).max())
        maps[w] = (U, fit)
    U0 = maps[wells[0]][0]
    spread = max(float(np.abs(U - U0).max() / np.abs(U0).max()) for U, _ in maps.values())
    worst_fit = max(f for _, f in maps.values())
    M_sf = np.linalg.pinv(U0)
    report += [f"{len(wells)} files, {len(det)} detectors, {len(fluors)} unmixed channels "
               f"({', '.join(fluors[:3])}, …, {fluors[-1]})",
               f"SpectroFlo's unmixed values as a linear map of the raw ones: largest error "
               f"{worst_fit:.2e} of the largest value; the map differs between files by at most "
               f"{spread:.2e} of its largest entry",
               f"its spectra pinv(U): row maxima {M_sf.max(axis=1).min():.6f} to "
               f"{M_sf.max(axis=1).max():.6f}"]
    sf = dict(zip(fluors, M_sf))
    controls, unstained, negative = [], None, None
    for w in wells:
        title, kind = _title(raw[w])
        if title.startswith("Unstained"):
            unstained = raw[w]
        elif title.startswith("Negative"):
            negative = raw[w]
        else:
            controls.append({"path": raw[w], "fluor": ALIAS.get(title, title) + "-A", "well": w})
    order = [c for c in fluors if c != "AF-A"]
    controls.sort(key=lambda c: order.index(c["fluor"]))
    results = {}
    for label, how in (("internal", {}), ("negative beads", {"negative": "unstained"})):
        detectors, refs, _ = um.read_references(
            [dict(c, **how) for c in controls], det,
            unstained=negative if how else None)
        cells = read_fcs(unstained)
        idx = [cells.index(d) for d in detectors]
        af = um.autofluorescence_of(cells.events[:, idx], [cells.ranges[i] for i in idx])
        refs.append(um.Reference("AF-A", af, source="unstained cells"))
        ours = um.Unmixing.from_references(detectors, refs)
        cosines = {r.fluor: cos(r.spectrum, sf[r.fluor]) for r in refs}
        diff = {r.fluor: float(np.abs(r.spectrum - sf[r.fluor]).max()) for r in refs}
        peaks = [r.fluor for r in refs[:-1] if r.peak != detectors[int(np.argmax(sf[r.fluor]))]]
        results[label] = {"unmixing": ours, "cosine": cosines, "maxdiff": diff, "peaks": peaks}
        worst = sorted(cosines.items(), key=lambda kv: kv[1])[:4]
        report += [f"\nspectra from the raw controls, negative = {label}:",
                   f"  cosine with SpectroFlo's: median {np.median(list(cosines.values())):.5f}, "
                   f"lowest {', '.join(f'{k} {v:.4f}' for k, v in worst)}; AF {cosines['AF-A']:.5f}",
                   f"  largest difference in one detector (peak = 1): median "
                   f"{np.median(list(diff.values())):.4f}, largest {max(diff.values()):.3f} "
                   f"({max(diff, key=diff.get)})",
                   f"  peak on another detector than SpectroFlo's: {', '.join(peaks) or 'none'}",
                   f"  complexity index: ours {ours.complexity():.2f}, SpectroFlo's spectra "
                   f"{um.complexity_index(M_sf):.2f}"]
    # unmixed values on the same events: each control's positive events in
    # its own channel, and the unstained cells in every channel
    ours = results["internal"]["unmixing"]
    sfu = um.Unmixing(det, fluors, M_sf)
    ratios, leaks_ours, leaks_sf = [], [], []
    for c in controls:
        f = read_fcs(c["path"])
        Y = f.events[:, [f.index(d) for d in det]].astype(float)
        a, b = ours.unmix(Y), sfu.unmix(Y)
        j = fluors.index(c["fluor"])
        top = b[:, j] >= np.quantile(b[:, j], 0.95)
        ratios.append(float(np.median(a[top, j]) / np.median(b[top, j])))
        others = [k for k in range(len(fluors)) if k != j and fluors[k] != "AF-A"]
        own = np.median(b[top, j])
        leaks_ours.append(float(np.abs(np.median(a[top][:, others], axis=0)).max() / own))
        leaks_sf.append(float(np.abs(np.median(b[top][:, others], axis=0)).max() / own))
    f0, fu = read_fcs(controls[0]["path"]), read_fcs(unmixed[controls[0]["well"]])
    exact = sfu.unmix(f0.events[:, [f0.index(d) for d in det]].astype(float))
    stored = fu.events[:, [fu.index(c) for c in fluors]].astype(float)
    report += ["\nunmixed values, ours (spectra from the controls, internal negative, AF) "
               "against SpectroFlo's, on the same events:",
               f"  positive beads, own channel, median ours / SpectroFlo: "
               f"{np.median(ratios):.4f} (from {min(ratios):.4f} to {max(ratios):.4f})",
               f"  positive beads, largest median in another channel, as a share of the own "
               f"channel: ours median {np.median(leaks_ours):.4f} (largest {max(leaks_ours):.4f}),"
               f" SpectroFlo median {np.median(leaks_sf):.4f} (largest {max(leaks_sf):.4f})",
               "  largest for SpectroFlo: " + ", ".join(
                   f"{controls[i]['fluor']} {leaks_sf[i]:.3f} (ours {leaks_ours[i]:.3f})"
                   for i in np.argsort(leaks_sf)[::-1][:3]),
               f"  with SpectroFlo's spectra, our OLS gives its stored values to "
               f"{np.abs(exact - stored).max() / np.abs(stored).max():.2e} of the largest"]
    # OLS and WLS on the same events: how far each spreads what should be 0
    import time
    cells = read_fcs(unstained)
    Y = cells.events[:, [cells.index(d) for d in det]].astype(float)
    t0 = time.perf_counter()
    a_ols = ours.unmix(Y, "ols")
    t1 = time.perf_counter()
    a_wls = ours.unmix(Y, "wls")
    t2 = time.perf_counter()

    def robust_sd(x):
        return 1.4826 * np.median(np.abs(x - np.median(x, axis=0)), axis=0)

    keep = [k for k, c in enumerate(fluors) if c != "AF-A"]
    ratio = robust_sd(a_wls)[keep] / robust_sd(a_ols)[keep]
    spread_pos = []
    for c in controls[:: max(1, len(controls) // 10)]:
        f = read_fcs(c["path"])
        Yc = f.events[:, [f.index(d) for d in det]].astype(float)
        o, w_ = ours.unmix(Yc, "ols"), ours.unmix(Yc, "wls")
        j = fluors.index(c["fluor"])
        top = o[:, j] >= np.quantile(o[:, j], 0.95)
        others = [k for k in keep if k != j]
        spread_pos.append(float(np.median(robust_sd(w_[top][:, others])
                                          / robust_sd(o[top][:, others]))))
    report += [f"\nWLS against OLS (our spectra): unstained cells ({Y.shape[0]} events), robust "
               f"SD of each channel, WLS / OLS: median {np.median(ratio):.3f} "
               f"(from {ratio.min():.3f} to {ratio.max():.3f}); brightest 5% of "
               f"{len(spread_pos)} bead controls, the other channels: median "
               f"{np.median(spread_pos):.3f} (from {min(spread_pos):.3f} to {max(spread_pos):.3f})",
               f"  time for {Y.shape[0]} events x {len(det)} detectors x {len(fluors)} spectra: "
               f"OLS {t1 - t0:.2f} s, WLS {t2 - t1:.1f} s"]
    return {"detectors": det, "fluors": fluors, "M_sf": M_sf, "unstained": unstained,
            "controls": controls, "ratios": ratios, "results": {
                k: {"cosine": v["cosine"], "maxdiff": v["maxdiff"], "peaks": v["peaks"]}
                for k, v in results.items()}}


def autospectral(found: dict, report: list[str], n: int = 3000) -> None:
    out = work("unmixing")
    det, fluors, M = found["detectors"], found["fluors"], found["M_sf"]
    rows = []
    for path in [found["unstained"]] + [c["path"] for c in found["controls"][:6]]:
        f = read_fcs(path)
        rows.append(f.events[:n // 7, [f.index(d) for d in det]].astype(float))
    Y = np.concatenate(rows)
    with open(os.path.join(out, "spectra.csv"), "w", encoding="utf-8") as f:
        f.write(",".join(["name", *det]) + "\n")
        for name, row in zip(fluors, M):
            f.write(",".join([f'"{name}"', *(repr(float(v)) for v in row)]) + "\n")
    np.savetxt(os.path.join(out, "raw.csv"), Y, delimiter=",", header=",".join(det),
               comments="", fmt="%.17g")
    print(rscript("unmixing_autospectral.R", AUTOSPECTRAL, out), end="")
    u = um.Unmixing(det, fluors, M)

    def load(name):
        return np.loadtxt(os.path.join(out, name), delimiter=",", skiprows=1)

    ols, wls = load("ols.csv"), load("wls_cell.csv")
    sim = np.loadtxt(os.path.join(out, "similarity.csv"), delimiter=",", skiprows=1,
                     usecols=range(1, len(fluors) + 1))
    ci_rounded, ci = (float(x) for x in open(os.path.join(out, "complexity.txt")).read().split())
    scale = np.abs(ols).max()
    report += [f"\nAutoSpectral's R functions on SpectroFlo's spectra and {Y.shape[0]} events:",
               f"  OLS: largest difference {np.abs(u.unmix(Y, 'ols') - ols).max() / scale:.1e} "
               "of the largest value",
               f"  per-event WLS: largest difference "
               f"{np.abs(u.unmix(Y, 'wls') - wls).max() / np.abs(wls).max():.1e} of the largest",
               f"  cosine similarity: largest difference {np.abs(u.similarity() - sim).max():.1e}",
               f"  condition number: ours {u.complexity():.6f}, AutoSpectral {ci:.6f} "
               f"(it prints {ci_rounded:g})"]


def library(report: list[str]) -> None:
    """The Aurora library AutoSpectral ships, read as a CSV of spectra."""
    path = os.path.join(DATA, "Aurora_spectral_reference_library.csv")
    if not os.path.exists(path):
        return
    dets, names, spectra = um.read_csv(path)
    report.append(f"\nAutoSpectral's Aurora library read by `read_csv`: {len(names)} spectra x "
                  f"{len(dets)} detectors")


def main() -> None:
    if "--fetch" in sys.argv:
        fetch(DATA)
        return
    if not DATA or not os.path.isdir(DATA):
        sys.exit("MILC_REF_SPECTRAL: the folder with the bead controls (see --fetch)")
    report: list[str] = []
    found = spectroflo(report)
    if AUTOSPECTRAL and os.path.isdir(AUTOSPECTRAL):
        autospectral(found, report)
    library(report)
    text = "\n".join(report)
    print(text)
    out = work("unmixing")
    with open(os.path.join(out, "report.md"), "w", encoding="utf-8") as f:
        f.write(text + "\n")
    with open(os.path.join(out, "report.json"), "w", encoding="utf-8") as f:
        json.dump({k: v for k, v in found.items() if k in ("ratios", "results")}, f, indent=1)


if __name__ == "__main__":
    main()
