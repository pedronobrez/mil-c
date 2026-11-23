"""
Differential expression within the clusters: does a cluster's cells carry
more or less of a marker in one group than in the other? The companion of
`abundance`, which asks whether there are more or fewer of the cells.

The method is diffcyt's DS-limma (Weber et al., Commun Biol 2019). Each
sample gives one number for each cluster and marker — the median of the
marker, on the scale it was transformed to, among the sample's cells of
the cluster — and those numbers are compared by a linear model

    median ~ intercept + group (+ one more factor, additive)

fitted for each cluster and marker. The group's coefficient is the
difference in median between the groups, in the units of the scale, and
it is tested by limma's moderated t (Smyth, Stat Appl Genet Mol Biol
2004): every row's residual variance is moved towards what all the rows
suggest, which is worth extra degrees of freedom where there are few
samples and many rows.

What is implemented, and where it is not diffcyt
------------------------------------------------
* The fit, the empirical-Bayes prior (d0 and s0², `abundance.squeeze`)
  and the moderated t are limma's `lmFit` + `eBayes`, and were checked
  against limma 3.68.5 on a small matrix (`tests/test_expression.py`:
  the numbers are R's), with and without a covariate, weights and
  missing values. No variance trend (`trend=FALSE`; diffcyt's default is
  TRUE), no robust option.
* Rows with missing samples have fewer residual degrees of freedom than
  the others. The prior is then fitted as `eBayes(legacy=TRUE)` does;
  limma's present default for unequal degrees of freedom is a newer
  estimator, which gives slightly different numbers (on the test
  matrix, 3.52 prior degrees of freedom against 3.10).
* A sample with fewer than `min_cells` cells of a cluster gives no
  median for it: a median of three cells is noise. (diffcyt keeps every
  median and drops whole clusters that are small in too many samples.)
  A cluster left with fewer than `min_samples` samples in either group
  is not tested, and says so; it is not dropped.
* Weights, optional and off unless asked for: each sample's median
  weighs as the number of cells it is the median of (diffcyt's default).
  That is right where medians differ because of how few cells they are
  of, and wrong where the animals themselves differ, which a count says
  nothing about: see the table below. Within a row the weights are
  scaled to a mean of 1, so that the residual variances of a large
  cluster and a small one are on one scale when the prior is fitted.
  (diffcyt gives limma the raw counts; the plain t statistics are the
  same, the moderation is not.)
* The other factor — a batch, a day, or the animal for a paired design —
  is a fixed effect: a column for each of its levels. (diffcyt's
  `block_id` is a random effect, by `duplicateCorrelation`.) With the
  animal as the factor and every animal in both groups the unmoderated
  test is the paired t test.
* Only the samples of the two groups compared are used.
* p values are adjusted over every cluster and marker tested together.

What it cannot say
------------------
A median is of the cluster as it was drawn. A marker the clusters were
built on is held within a narrow range inside each cluster by
construction: what a group does to it shows as cells moving between
clusters (abundance), and only its remainder as a shift here. And a
difference in a bimodal marker — a quarter of the cluster switching it
on — moves the median little or not at all.

How well calibrated it is
-------------------------
Simulated (`simulate` and `tally` in `tests/test_expression.py`): 12
clusters x 6 markers, each row varying between animals with a standard
deviation of its own (0.05 to 0.4, log-uniform), cell counts from 20 to
5,000 between clusters and five-fold between samples, each median with
the sampling error of its count on top (1 / sqrt(cells)). Under the null,
the share of p < 0.05 and of p < 0.01. With a shift of 0.5 in 6 of the 72
rows, the share of those found at FDR < 0.05 and the false discovery rate
had. "Plain t" is the same model without the moderation. 1,000 data sets
each.

                           null p < 0.05   p < 0.01   power   FDR had
    5 v 5   unweighted         5.1%          1.1%      69%      4.9%
            weighted           6.6%          1.5%      71%      7.1%
            plain t            5.0%          1.0%      58%      4.5%
    3 v 3   unweighted         5.4%          1.1%      33%      5.0%
            weighted           6.5%          1.5%      38%      6.3%
            plain t            5.0%          0.9%      16%      3.6%
    2 v 2   unweighted         5.7%          1.1%       3%      3.6%
            weighted           6.7%          1.5%       5%      5.1%
            plain t            4.9%          1.0%       1%      4.4%

So the moderation finds more than the plain t, twice as much with three
a group, at a cost of slightly too many small p values with few samples;
the weights buy little and cost calibration, which is why they are off
by default; and two a group finds next to nothing.

The simulation shares the model's assumptions — medians that vary
between animals as a normal does, variances not too unlike between
rows. Nothing in a simulation says real animals do.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import compare
from .abundance import design, squeeze
from .stats import StatsTable

#: fewer cells of a cluster than this in a sample: no median from it
MIN_CELLS = 10
#: fewer samples than this with a median in either group: no test
MIN_SAMPLES = 2


@dataclass
class Medians:
    """What each sample's cells of each cluster express."""
    clusters: list[str]
    samples: list[str]
    markers: list[str]
    #: clusters x samples x markers; NaN where the sample has no cell of
    #: the cluster
    values: np.ndarray
    #: clusters x samples: the cells each median is of
    counts: np.ndarray
    #: the scale the values are on: `biex` or `arcsinh`
    scale: str = ""


@dataclass
class Result:
    cluster: str
    marker: str
    #: (the group the difference is of, the reference)
    groups: tuple[str, str]
    #: the samples of each group that gave a median
    n: tuple[int, int]
    #: each group's mean of its samples' medians
    medians: tuple[float, float]
    #: the group's coefficient: first less second, adjusted for the rest
    difference: float = float("nan")
    t: float = float("nan")
    #: residual degrees of freedom plus those the other rows were worth
    df: float = float("nan")
    p: float = float("nan")
    p_adjusted: float = float("nan")
    note: str = ""
    #: the value of "within each" this was computed in, if any
    panel: str | None = None


def medians(X: np.ndarray, labels: np.ndarray, sample_index: np.ndarray,
            n_samples: int) -> tuple[list[int], np.ndarray, np.ndarray]:
    """(the clusters, clusters x samples x columns of `X`: the median of
    each column among each sample's events of each cluster, clusters x
    samples: how many events that is). No events, NaN."""
    ids, code = np.unique(labels, return_inverse=True)
    code = code * n_samples + sample_index
    order = np.argsort(code, kind="stable")
    found, starts = np.unique(code[order], return_index=True)
    ends = np.append(starts[1:], code.size)
    values = np.full((ids.size * n_samples, X.shape[1]), np.nan)
    counts = np.zeros(ids.size * n_samples)
    for k, a, b in zip(found, starts, ends):
        values[k] = np.median(X[order[a:b]], axis=0)
        counts[k] = b - a
    return ([int(c) for c in ids], values.reshape(ids.size, n_samples, X.shape[1]),
            counts.reshape(ids.size, n_samples))


# --- the model -------------------------------------------------------------------

@dataclass
class Fit:
    """One linear model per row. NaN where a row could not be fitted;
    `note` says why."""
    coefficient: np.ndarray
    #: the coefficient's standard error over the residual standard deviation
    unscaled: np.ndarray
    #: residual variance and its degrees of freedom
    s2: np.ndarray
    df: np.ndarray
    note: list[str]


def _tested_design(X: np.ndarray):
    """`X` with the column tested first and the others as an orthonormal
    basis of what they span — a level of the covariate that no sample
    left has is no column — or why the coefficient cannot be estimated."""
    g = X[:, 1]
    u, s, _ = np.linalg.svd(np.delete(X, 1, axis=1), full_matrices=False)
    D = np.column_stack([g, u[:, s > 1e-10 * max(s.max(), 1.0)]])
    if np.linalg.matrix_rank(D) < D.shape[1]:
        return None, "what is adjusted for cannot be told apart from the groups"
    if X.shape[0] <= D.shape[1]:
        return None, "too few samples" + (" to adjust for anything" if X.shape[1] > 2 else "")
    return D, ""


def fit(y: np.ndarray, X: np.ndarray, weights: np.ndarray | None = None) -> Fit:
    """Least squares of each row of `y` (rows x samples) on `X` (samples x
    coefficients) for the coefficient of its second column, as limma's
    `lmFit`. A NaN in `y` is a sample the row does not have. `weights`,
    of `y`'s shape, are relative: a sample of weight 2 counts as two."""
    y = np.atleast_2d(np.asarray(y, dtype=np.float64))
    X = np.asarray(X, dtype=np.float64)
    rows = y.shape[0]
    w = np.ones(y.shape) if weights is None else np.asarray(weights, dtype=np.float64)
    have = np.isfinite(y) & np.isfinite(w) & (w > 0)
    out = Fit(*(np.full(rows, np.nan) for _ in range(4)), [""] * rows)
    for pattern in np.unique(have, axis=0):
        which = np.flatnonzero((have == pattern).all(axis=1))
        D, reason = _tested_design(X[pattern]) if pattern.any() else (None, "no samples")
        if D is None:
            for i in which:
                out.note[i] = reason
            continue
        df = D.shape[0] - D.shape[1]
        for i in which:
            wi, yi = w[i, pattern], y[i, pattern]
            inverse = np.linalg.inv(D.T @ (wi[:, None] * D))
            beta = inverse @ (D.T @ (wi * yi))
            residual = yi - D @ beta
            out.coefficient[i], out.unscaled[i] = beta[0], np.sqrt(inverse[0, 0])
            out.s2[i], out.df[i] = float(wi @ residual ** 2) / df, df
    return out


def moderated(f: Fit, moderate: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """(t, its degrees of freedom, two-sided p, the prior degrees of
    freedom) of each row's coefficient, with the residual variances
    squeezed towards each other as limma's `eBayes` does; without
    `moderate`, the ordinary t test of each row by itself."""
    from scipy import stats as st
    ok = np.isfinite(f.s2) & (f.df > 0)
    s2, df, prior_df = f.s2.copy(), f.df.copy(), 0.0
    if moderate and ok.sum() > 1:
        s2[ok], prior_df = squeeze(f.s2[ok], f.df[ok])
        # all the rows together cannot be worth more than they have
        df[ok] = np.minimum(f.df[ok] + prior_df, f.df[ok].sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        t = f.coefficient / (np.sqrt(s2) * f.unscaled)
    t = np.where(ok & (f.coefficient == 0), 0.0, t)   # no difference: whatever the variance
    t[~ok] = np.nan
    p = np.where(ok, 2.0 * st.t.sf(np.abs(t), np.where(ok, df, 1.0)), np.nan)
    return t, np.where(ok, df, np.nan), p, prior_df


# --- the comparison ------------------------------------------------------------------

def differential(values: np.ndarray, counts: np.ndarray, clusters: list[str],
                 markers: list[str], groups: list[str], contrast: tuple[str, str],
                 covariate: list[str] | None = None, weights: bool = False,
                 min_cells: int = MIN_CELLS, min_samples: int = MIN_SAMPLES,
                 correction: str = compare.BH, moderate: bool = True) -> list[Result]:
    """The median of each marker in each cluster in `contrast[0]` against
    `contrast[1]`: one `Result` for each cluster and marker, the clusters
    in order and the markers within each. `values` is clusters x samples
    x markers and `counts` clusters x samples, as `medians` gives them;
    `groups` and `covariate` name each sample's group and, optionally,
    the level of a factor to adjust for — the animal, for a paired
    design.

    Where no test can be made the p is NaN and the note says why; the
    difference is then of the groups' plain means, where both have one."""
    values = np.asarray(values, dtype=np.float64)
    counts = np.asarray(counts, dtype=np.float64)
    keep = [i for i, g in enumerate(groups) if g in contrast]
    kept = [groups[i] for i in keep]
    X = design(kept, contrast, [covariate[i] for i in keep] if covariate is not None else None)
    first = X[:, 1] == 1
    n_markers = len(markers)
    y = np.full((len(clusters) * n_markers, len(keep)), np.nan)
    w = np.ones(y.shape)
    out: list[Result] = []
    for c, cluster in enumerate(clusters):
        enough = counts[c, keep] >= max(min_cells, 1)
        for m, marker in enumerate(markers):
            row = np.where(enough, values[c, keep, m], np.nan)
            have = np.isfinite(row)
            n = (int((have & first).sum()), int((have & ~first).sum()))
            with np.errstate(invalid="ignore"):
                means = tuple(float(row[have & side].mean()) if (have & side).any()
                              else float("nan") for side in (first, ~first))
            r = Result(str(cluster), str(marker), tuple(contrast), n, means)
            r.difference = means[0] - means[1]
            if contrast[0] == contrast[1] or min(kept.count(g) for g in contrast) == 0:
                r.note = "one group: nothing to compare"
            elif min(n) < min_samples:
                short = contrast[0] if n[0] < min_samples else contrast[1]
                r.note = (f"fewer than {min_samples} samples of {short or '—'} with "
                          f"{min_cells} cells or more")
            else:
                y[c * n_markers + m] = row
                if weights:
                    w[c * n_markers + m] = np.where(
                        have, counts[c, keep] / counts[c, keep][have].mean(), 1.0)
            out.append(r)
    f = fit(y, X, w if weights else None)
    t, df, p, _ = moderated(f, moderate)
    for i, r in enumerate(out):
        if np.isfinite(f.coefficient[i]):
            r.difference = float(f.coefficient[i])
        r.t, r.df, r.p = float(t[i]), float(df[i]), float(p[i])
        r.note = r.note or f.note[i]
    return adjusted(out, correction)


def adjusted(results: list[Result], correction: str = compare.BH) -> list[Result]:
    """The same results with p adjusted over all of them."""
    for r, p in zip(results, compare.adjust([r.p for r in results], correction)):
        r.p_adjusted = float(p)
    return results


def within(values: np.ndarray, counts: np.ndarray, clusters: list[str], markers: list[str],
           groups: list[str], contrast: tuple[str, str], panels: list[str],
           covariate: list[str] | None = None, weights: bool = False,
           min_cells: int = MIN_CELLS, min_samples: int = MIN_SAMPLES,
           correction: str = compare.BH) -> list[Result]:
    """The comparison made separately among the samples of each value of
    `panels` — within each day, say — prior and all, with the p values
    then adjusted over every panel's tests together."""
    values, counts = np.asarray(values), np.asarray(counts)
    out: list[Result] = []
    for panel in compare.ordered(panels):
        keep = [i for i, p in enumerate(panels) if p == panel]
        found = differential(values[:, keep], counts[:, keep], clusters, markers,
                             [groups[i] for i in keep], contrast,
                             [covariate[i] for i in keep] if covariate is not None else None,
                             weights, min_cells, min_samples, compare.NONE)
        for r in found:
            r.panel = panel
        out += found
    return adjusted(out, correction)


def to_table(results: list[Result], correction: str = compare.BH,
             panel: str | None = None) -> StatsTable:
    """Every cluster and marker as a row; `panel` names the first column
    where the results are of several panels."""
    adjusted_name = "p" if correction == compare.NONE else f"p adjusted ({correction})"
    a, b = results[0].groups if results else ("", "")
    table = StatsTable(([panel] if panel else [])
                       + ["Cluster", "Marker", f"{a or '—'} (median)", f"{b or '—'} (median)",
                          "Difference", "p", adjusted_name, "", "n", "Moderated t", "df",
                          "Note"])
    for r in results:
        table.rows.append(([r.panel] if panel else [])
                          + [r.cluster, r.marker, r.medians[0], r.medians[1], r.difference,
                             r.p, r.p_adjusted, compare.stars(r.p_adjusted),
                             f"{r.n[0]} vs {r.n[1]}", r.t, r.df, r.note])
    return table
