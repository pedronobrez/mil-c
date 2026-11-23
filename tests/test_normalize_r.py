"""
The curves, against R.

CytoNorm's QuantileNorm is, for one cluster, three lines of base R:
quantiles at 0.01 to 0.99, the mean of the batches' quantiles as the goal,
and `splinefun(method = "monoH.FC")` from one to the other. Written out
here in R and run on the same numbers as `normalize.fit` / `apply`.

Measured on 2026-11-23 on the real experiment (three batches, 891 to
2,551 reference events, eight markers): inside the range the references
cover, the two differ by 1e-5 of the axis on average, against a mean
correction of 0.014 — the same curve — and by at most 0.005. Beyond the
outermost quantiles they differ by up to 0.21 of the axis, on purpose:
R continues the spline's end slope, this continues as a shift.

`PchipInterpolator` and R's Fritsch-Carlson spline are both monotone and
pass through the same knots; they choose the slopes at the knots a little
differently, which shows where knots are far apart.

The whole CytoNorm package (with its clustering) is not what this compares
with: that is `tests/reference/normalize_compare.py`, run by hand, and
written up in `docs/comparacao.md`. Skipped where there is no Rscript.
"""

import shutil
import subprocess

import numpy as np
import pytest

from milc import normalize

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="no R here")

R = """
args <- commandArgs(trailingOnly = TRUE); dir <- args[1]
batches <- c("a", "b", "c"); probs <- (1:99) / 100
ref <- lapply(batches, function(b) as.matrix(read.csv(file.path(dir, paste0("ref_", b, ".csv")), header = FALSE)))
q <- lapply(ref, function(m) apply(m, 2, quantile, probs = probs, type = 7))
goal <- Reduce(`+`, q) / length(q)
for (i in seq_along(batches)) {
  x <- as.matrix(read.csv(file.path(dir, paste0("x_", batches[i], ".csv")), header = FALSE))
  out <- x
  for (j in seq_len(ncol(x))) out[, j] <- splinefun(q[[i]][, j], goal[, j], method = "monoH.FC")(x[, j])
  write.table(out, file.path(dir, paste0("r_", batches[i], ".csv")), sep = ",", row.names = FALSE, col.names = FALSE)
}
"""


def test_the_curve_is_r_s_curve(tmp_path):
    rng = np.random.default_rng(0)

    def events(n, shift, scale):
        two = rng.random(n) < 0.4
        base = np.where(two[:, None], rng.normal(0.7, 0.05, (n, 3)), rng.normal(0.3, 0.04, (n, 3)))
        return base * scale + shift

    effects = {"a": (0.0, 1.0), "b": (0.05, 1.1), "c": (-0.03, 0.9)}
    reference = {b: events(3000, *e) for b, e in effects.items()}
    others = {b: events(2000, *e) for b, e in effects.items()}
    model = normalize.fit(reference)
    for b in effects:
        np.savetxt(tmp_path / f"ref_{b}.csv", reference[b], delimiter=",", fmt="%.12f")
        np.savetxt(tmp_path / f"x_{b}.csv", others[b], delimiter=",", fmt="%.12f")
    (tmp_path / "core.R").write_text(R)
    subprocess.run(["Rscript", str(tmp_path / "core.R"), str(tmp_path)], check=True,
                   capture_output=True, timeout=120)
    for b in effects:
        ours = normalize.apply(model, others[b], b)
        theirs = np.loadtxt(tmp_path / f"r_{b}.csv", delimiter=",")
        lo, hi = np.quantile(reference[b], [0.01, 0.99], axis=0)
        inside = (others[b] >= lo) & (others[b] <= hi)
        moved = np.abs(ours - others[b])[inside].mean()
        gap = np.abs(ours - theirs)[inside]
        # the same curve: on these numbers the two are apart by about a
        # hundredth of what the curve does (a thousandth on the real data)
        assert gap.mean() < moved / 50 and np.quantile(gap, 0.99) < moved / 5
        # the largest gaps are in the empty stretch between two populations,
        # where the knots are far apart and two monotone splines through the
        # same points need not agree; few events are there
        assert gap.max() < 0.05
        assert inside.mean() > 0.9
