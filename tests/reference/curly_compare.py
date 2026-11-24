"""MIL-C's curly quadrant beside CytoML's, on one workspace.

The workspace and its file are made by `tests/test_curly_quadrant.py`
(`curly_files`): simulated events with the spread of photon counting,
compensated, and three curly crossings, written by MIL-C's own workspace
writer in the form FlowJo writes a curly quadrant. It was not made by
FlowJo. CytoML reads it and gives, event by event, what it puts in each
population (`curly_cytoml.R`); MIL-C reads it with `wsp.read_wsp`.

Printed, for each corner: CytoML's count, ours, the events the two put
differently, and — to show how much the unknown curve matters — the
count with the curve another program fitted to FlowJo's counts (GateLab:
1.5 for the power, 0.012 for the bend, on 256 channels; not checked
here) and with the arms straight.

Usage: python3 tests/reference/curly_compare.py   (R with CytoML on R_LIBS_USER)
"""

from __future__ import annotations

import csv
import importlib.util
import json
import os

import numpy as np

from common import HERE, rscript, work

from scipy.interpolate import CubicSpline

from milc import wsp
from milc.gating import CurlyQuadrantGate, copy_tree
from milc.transforms import Biex

spec = importlib.util.spec_from_file_location(
    "curly_test", os.path.join(os.path.dirname(HERE), "test_curly_quadrant.py"))
curly_test = importlib.util.module_from_spec(spec)
spec.loader.exec_module(curly_test)


def cytolib_log_root(b: float, w: float) -> float:
    """cytolib's biexpTrans::logRoot as it is, with its casts to whole
    numbers: it stops as soon as a Newton step is under one, so the root —
    and the negative side of the scale — is not the exact one."""
    trunc = lambda v: float(abs(int(v)))  # noqa: E731
    lo, hi = 0.0, b
    d = (lo + hi) / 2
    last = trunc(lo - hi)
    fb = -2 * np.log(b) + w * b
    f = 2.0 * np.log(d) + w * b + fb
    df = 2 / d + w
    if w == 0:
        return b
    for _ in range(100):
        if ((d - hi) * df - f) * ((d - lo) * df - f) >= 0 or trunc(2 * f) > trunc(last * df):
            dx = (hi - lo) / 2
            d = lo + dx
            if d == lo:
                return d
        else:
            dx = f / df
            t = d
            d -= dx
            if d == t:
                return d
        if trunc(dx) < 1.0e-12:
            return d
        last = dx
        f = 2 * np.log(d) + w * d + fb
        df = 2 / d + w
        if f < 0:
            lo = d
        else:
            hi = d
    return d


class CytolibBiex(Biex):
    """FlowJo's biexponential as cytolib builds it for a curly quadrant's
    curve (transformation.cpp, biexpTrans::computCalTbl after
    setTransformedScale(255)): its table on 255 channels, not 4096 — the
    channel of zero, a whole number, then falls elsewhere and the decades
    are stretched to meet it — and read through a natural spline. Only to
    tell where the two programs part."""

    def _table(self):
        channels = 255
        ln10 = np.log(10.0)
        decades = self.pos
        width = np.log10(-self.width)
        if width < 0.5 or width > 3:
            width = 0.5
        decades -= width / 2
        extra = max(self.neg, 0.0) + width / 2
        zero = min(int(extra * channels / (extra + decades)), channels // 2)
        if zero > 0:
            decades = extra * channels / zero
        width /= 2 * decades
        pos_range = ln10 * decades
        minimum = self.max_range / np.exp(pos_range)
        neg_range = cytolib_log_root(pos_range, width)
        n = channels + 1
        j = np.arange(n, dtype=np.float64)
        positive = np.exp(j / n * pos_range)
        negative = np.exp(j / n * -neg_range) * np.exp((pos_range + neg_range)
                                                       * (width + extra / decades))
        s = positive[zero] - negative[zero]
        data = np.empty(n)
        data[zero:] = minimum * (positive[zero:] - negative[zero:] - s)
        data[:zero] = -data[2 * zero - np.arange(zero)]
        self._forward = CubicSpline(data, j / channels, bc_type="natural")
        self._back = CubicSpline(j / channels, data, bc_type="natural")
        return data, j / channels

    def __call__(self, x):
        return self._forward(np.asarray(x, dtype=np.float64))

    def inverse(self, y):
        return self._back(np.asarray(y, dtype=np.float64))


def as_cytolib(tree):
    out = copy_tree(tree)
    for node in out.walk():
        g = node.gate
        if isinstance(g, CurlyQuadrantGate):
            g.scales = [CytolibBiex(s.length, s.max_range, s.neg, s.width, s.pos)
                        if isinstance(s, Biex) else s for s in g.scales]
    return out


def with_curve(tree, bend, power, channels):
    out = copy_tree(tree)
    for node in out.walk():
        if isinstance(node.gate, CurlyQuadrantGate):
            node.gate.bend, node.gate.power, node.gate.channels = bend, power, channels
    return out


def main() -> None:
    folder = work("curly")
    fcs, path = curly_test.curly_files(folder)
    printed = rscript("curly_cytoml.R", path, folder, folder)
    with open(os.path.join(folder, "cytoml_members.csv")) as f:
        rows = list(csv.reader(f))
    names = rows[0]
    members = np.array(rows[1:], dtype=np.int8).astype(bool)
    theirs = {n: members[:, k] for k, n in enumerate(names)}

    read = wsp.read_wsp(path).samples[0]
    variants = {"ours": read.gates, "cytolib": as_cytolib(read.gates),
                "gatelab": with_curve(read.gates, 0.012, 1.5, 256.0),
                "straight": with_curve(read.gates, 0.0, 1.0, 255.0)}
    masks = {k: curly_test.opened(fcs, t).masks() for k, t in variants.items()}
    table = []
    for node in read.gates.walk():
        if node.gate is None:
            continue
        p = node.path
        ours, cyto = masks["ours"][p], theirs[p]
        table.append({"population": p, "kind": type(node.gate).__name__,
                      "cytoml": int(cyto.sum()), "ours": int(ours.sum()),
                      "differ": int((ours != cyto).sum()),
                      "differ_cytolib": int((masks["cytolib"][p] != cyto).sum()),
                      "gatelab": int(masks["gatelab"][p].sum()),
                      "straight": int(masks["straight"][p].sum())})
    lines = ["| population | CytoML | MIL-C | differ | differ, cytolib's table "
             "| GateLab's curve | straight |", "|---|---:|---:|---:|---:|---:|---:|"]
    for r in table:
        lines.append(f"| {r['population']} | {r['cytoml']:,d} | {r['ours']:,d} | {r['differ']} "
                     f"| {r['differ_cytolib']} | {r['gatelab']:,d} | {r['straight']:,d} |")
    report = "\n".join(lines) + "\n\n" + printed
    print(report)
    with open(os.path.join(folder, "report.md"), "w") as f:
        f.write(report)
    with open(os.path.join(folder, "report.json"), "w") as f:
        json.dump(table, f, indent=1)


if __name__ == "__main__":
    main()
