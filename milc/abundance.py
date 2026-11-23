"""
Differential abundance: is a cluster a larger or smaller share of the
sample in one group than in the other? Asked of the *counts*, not of the
percentages.

A percentage forgets what it was made of: 0.2% from 109 cells and 0.2%
from 50,000 are the same number to a t test, and a rare cluster's
percentages are far from normal. The field's answer (diffcyt's DA-edgeR:
Weber et al., Commun Biol 2019) is a negative-binomial model of each
cluster's count per sample, with the sample's total as offset:

    count ~ NB(mean = total x exp(X b), variance = mean + dispersion x mean^2)

X has an intercept, the group, and — optionally — one more factor to
adjust for (a batch, a day) as additive terms. The group's coefficient is
the log fold change, and it is tested.

What is implemented, and where it is not edgeR
----------------------------------------------
Nothing here was run beside edgeR or diffcyt: R is not a dependency and
the numbers were not compared. The method is "in the manner of", and this
is the manner exactly:

* The GLM is fitted by iteratively reweighted least squares with step
  halving. (edgeR: Levenberg-damped Fisher scoring; same maximum.)
* The offset is the log of the sample's total over *all* clusters. No
  TMM normalisation factors (diffcyt's default has none either), so a
  cluster that grows makes every other one a smaller share: the test is
  of shares, not of absolute numbers.
* Dispersion. One common dispersion for all clusters, maximising the sum
  of the Cox-Reid adjusted profile likelihoods (edgeR's
  `estimateGLMCommonDisp`), found on a grid of 2^-20 ... 2^3 in steps of
  2^0.25 with a parabola through the best point. No trended dispersion
  (diffcyt's default has none: there is no mean-dispersion trend to speak
  of across a dozen clusters).
* Test, by default `QL`: a quasi-likelihood F test as edgeR's
  `glmQLFit`/`glmQLFTest` before version 4. The NB model is fitted at the
  common dispersion; each cluster's quasi-dispersion is its deviance over
  its residual degrees of freedom (samples whose fitted value is zero do
  not count, as in edgeR), squeezed towards the others' by limma's
  empirical Bayes (`squeezeVar`: Smyth 2004) — without its trend on
  abundance, without its robust option, without edgeR 4's bias-adjusted
  deviances and without `poisson.bound`. F = likelihood ratio / squeezed
  quasi-dispersion on (1, residual + prior) degrees of freedom.
* `LRT` is the other test, and is what diffcyt's `testDA_edgeR` calls:
  cluster-wise dispersions maximising each cluster's adjusted profile
  likelihood plus `PRIOR_DF / residual df` times the average of all
  clusters' (edgeR's weighted likelihood, with the prior degrees of
  freedom fixed at 10 as in `estimateGLMTagwiseDisp`; `estimateDisp`
  estimates them from the data), then a chi-square likelihood-ratio test.
* The fold change reported is from a second fit with a prior count of
  0.125 per sample, scaled by its total, as `glmFit(prior.count=0.125)`
  does: a cluster absent from a whole group gets a large finite fold
  change instead of an infinite one, and a note. The test does not use
  the prior count.
* Only the samples of the two groups compared are used, for the
  dispersion too. (edgeR would fit every group and test a contrast.)
* No filter on rare clusters (diffcyt drops those with fewer than 3 cells
  in half the samples); a cluster with no cells in either group is not
  tested.

How well calibrated it is
-------------------------
Simulated (`simulate` and `tally` in `tests/test_abundance.py`): clusters
from 0.1% to 30% of the sample, totals varying five-fold between samples,
each cluster's share varying between animals with a dispersion of its own
(0.02 to 0.3, log-uniform), counts multinomial. Under the null, the share
of p < 0.05 (5% is right). With a true two-fold change in 3 clusters, the
share of those found at FDR < 0.05 and the false discovery rate actually
had. "t" is Welch's test on the percentages (`compare.compare`), which
the Clusters tab already offers. 2,000 data sets for the first line of
each block, 1,000 for the others.

                               null p < 0.05    power    FDR had
    5 v 5, 15 clusters   QL        5.3%          65%       4.9%
                         LRT       5.3%          72%       6.3%
                         t         4.5%          34%       2.9%
    3 v 3, 15 clusters   QL        5.8%          36%       4.7%
                         LRT       6.0%          45%       6.7%
                         t         3.4%           6%       1.7%
    5 v 5, 8 clusters    QL        5.6%          67%       3.7%
                         LRT       5.7%          73%       4.1%
                         t         4.3%          35%       1.8%

So the counts find about twice what the percentages do, and both tests
give somewhat too many small p values: QL the fewer, which is why it is
the default. Further out it is worse: under the null, p < 0.01 came 1.1%
of the time with QL and 1.5% with LRT (5 v 5), and with 2 v 2 (400 data
sets) 1.7% and 2.3%. Two samples a group is not enough for either.

The simulation and the model share their assumption — counts that vary
between animals much as a negative binomial, with dispersions not too
unlike between clusters. Nothing in a simulation says real animals do.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import compare
from .stats import StatsTable

QL = "ql"
LRT = "lrt"
TESTS = {QL: "Quasi-likelihood F test", LRT: "Likelihood-ratio test"}

#: as edgeR's glmFit: added to every count, in proportion to the sample's
#: total, for the fold change that is reported — never for the test
PRIOR_COUNT = 0.125
#: the weight of all clusters' likelihood in one cluster's dispersion, as
#: degrees of freedom (LRT only)
PRIOR_DF = 10.0
#: the dispersions tried
GRID = 2.0 ** np.arange(-20.0, 3.25, 0.25)

#: a fitted value below this, where the count is zero, is a zero: the
#: coefficient went to minus infinity and the sample says nothing more
_ZERO = 1e-5


@dataclass
class Result:
    cluster: str
    #: (the group in the numerator, the reference)
    groups: tuple[str, str]
    n: tuple[int, int]
    #: each group's mean percent of the sample, for reading beside the fold
    percent: tuple[float, float]
    log2_fold: float = float("nan")
    test: str = ""
    statistic: float = float("nan")
    p: float = float("nan")
    p_adjusted: float = float("nan")
    #: of the NB model the test was made with
    dispersion: float = float("nan")
    note: str = ""
    #: the value of "within each" this was computed in, if any
    panel: str | None = None


# --- the model -------------------------------------------------------------------

def log_likelihood(y: np.ndarray, mu: np.ndarray, dispersion) -> np.ndarray:
    """Of each count under a negative binomial of mean `mu` and variance
    mu + dispersion x mu^2; Poisson where the dispersion is 0."""
    from scipy.special import gammaln, xlogy
    y, mu = np.asarray(y, dtype=np.float64), np.asarray(mu, dtype=np.float64)
    phi = np.broadcast_to(np.asarray(dispersion, dtype=np.float64), np.broadcast(y, mu).shape)
    poisson = xlogy(y, mu) - mu - gammaln(y + 1)
    r = 1.0 / np.where(phi > 0, phi, 1.0)
    nb = (gammaln(y + r) - gammaln(r) - gammaln(y + 1)
          - r * np.log1p(phi * mu) + xlogy(y, phi * mu) - y * np.log1p(phi * mu))
    return np.where(phi > 0, nb, poisson)


def deviance(y: np.ndarray, mu: np.ndarray, dispersion) -> np.ndarray:
    """Twice the log-likelihood lost by `mu` against a mean equal to each
    count, summed over the last axis."""
    from scipy.special import xlogy
    phi = np.broadcast_to(np.asarray(dispersion, dtype=np.float64), np.broadcast(y, mu).shape)
    safe = np.where(phi > 0, phi, 1.0)
    nb = xlogy(y, y / mu) - (y + 1 / safe) * (np.log1p(safe * y) - np.log1p(safe * mu))
    poisson = xlogy(y, y / mu) - (y - mu)
    return 2.0 * np.where(phi > 0, nb, poisson).sum(axis=-1)


def glm_fit(y: np.ndarray, X: np.ndarray, offset: np.ndarray, dispersion=0.0,
            max_iter: int = 200, tol: float = 1e-9) -> tuple[np.ndarray, np.ndarray]:
    """Negative-binomial regressions with a log link, one for each row of
    `y` (rows x samples) on the same `X` (samples x coefficients):
    mean = exp(X b + offset). `dispersion` is one number or one per row.
    Returns the coefficients (rows x coefficients) and the fitted means.

    Iteratively reweighted least squares. A step that makes the deviance
    worse is halved; a coefficient on its way to minus infinity (a group
    with no cells) stops when the deviance no longer changes."""
    y = np.atleast_2d(np.asarray(y, dtype=np.float64))
    X = np.asarray(X, dtype=np.float64)
    phi = np.broadcast_to(np.asarray(dispersion, dtype=np.float64), (y.shape[0],))[:, None]

    def mean(beta):
        return np.exp(np.clip(beta @ X.T + offset, -700.0, 700.0))

    beta = np.linalg.lstsq(X, (np.log(y + 0.5) - offset).T, rcond=None)[0].T
    mu = mean(beta)
    dev = deviance(y, mu, phi)
    for _ in range(max_iter):
        w = mu / (1.0 + phi * mu)
        A = np.einsum("kn,ni,nj->kij", w, X, X)
        b = ((y - mu) / (1.0 + phi * mu)) @ X
        # damping changes the path, not where it ends: at the maximum b is 0
        A += (1e-10 * np.trace(A, axis1=1, axis2=2) / X.shape[1] + 1e-300)[:, None, None] \
            * np.eye(X.shape[1])
        step = np.linalg.solve(A, b[:, :, None])[:, :, 0]
        for _ in range(30):
            new = beta + step
            new_mu = mean(new)
            new_dev = deviance(y, new_mu, phi)
            worse = ~(new_dev <= dev + 1e-9 * (np.abs(dev) + 1.0))
            if not worse.any():
                break
            step[worse] *= 0.5
        new[worse], new_mu[worse], new_dev[worse] = beta[worse], mu[worse], dev[worse]
        done = np.abs(dev - new_dev) <= tol * (np.abs(dev) + 1.0)
        beta, mu, dev = new, new_mu, new_dev
        if done.all():
            break
    return beta, mu


def _zeros(y: np.ndarray, mu: np.ndarray) -> np.ndarray:
    return (y == 0) & (mu < _ZERO)


def adjusted_profile_likelihood(y: np.ndarray, X: np.ndarray, offset: np.ndarray,
                                dispersion) -> np.ndarray:
    """Cox and Reid's: each row's log-likelihood at the coefficients that
    maximise it for this dispersion, less half the log determinant of
    their information — the correction for having estimated them, which
    is what keeps a dispersion from few samples from coming out too small.

    A sample whose fitted value is zero has no information; it is given a
    token weight, the same for every dispersion, so that it shifts the
    likelihood by a constant and does not move its maximum."""
    y = np.atleast_2d(np.asarray(y, dtype=np.float64))
    phi = np.broadcast_to(np.asarray(dispersion, dtype=np.float64), (y.shape[0],))[:, None]
    _, mu = glm_fit(y, X, offset, phi[:, 0])
    w = np.where(_zeros(y, mu), 1e-10, mu / (1.0 + phi * mu))
    logdet = np.linalg.slogdet(np.einsum("kn,ni,nj->kij", w, X, X))[1]
    return log_likelihood(y, mu, phi).sum(axis=1) - 0.5 * logdet


def _argmax(curve: np.ndarray) -> float:
    """Where on `GRID` the curve is highest, as a dispersion: the vertex
    of the parabola through the best point and its neighbours."""
    i = int(np.argmax(curve))
    shift = 0.0
    if 0 < i < curve.size - 1:
        a, b, c = curve[i - 1], curve[i], curve[i + 1]
        if a + c - 2 * b < 0:
            shift = 0.5 * (a - c) / (a + c - 2 * b)
    step = np.log(GRID[1] / GRID[0])
    return float(GRID[i] * np.exp(shift * step))


def dispersions(y: np.ndarray, X: np.ndarray, offset: np.ndarray,
                prior_df: float = PRIOR_DF) -> tuple[float, np.ndarray]:
    """(the common dispersion, one for each row shrunk towards it).

    The common one maximises the sum of the rows' adjusted profile
    likelihoods. A row's own maximises its likelihood plus the average of
    all rows' weighted as `prior_df` degrees of freedom against the
    residual ones it has itself: with 5 against 5, its own evidence counts
    for 8 and everybody's for 10."""
    y = np.atleast_2d(np.asarray(y, dtype=np.float64))
    k, n = y.shape
    apl = adjusted_profile_likelihood(np.tile(y, (GRID.size, 1)), X, offset,
                                      np.repeat(GRID, k)).reshape(GRID.size, k)
    shared = apl.mean(axis=1)
    weight = prior_df / max(n - X.shape[1], 1)
    return _argmax(shared), np.array([_argmax(apl[:, i] + weight * shared) for i in range(k)])


def _trigamma_inverse(x: float) -> float:
    from scipy.special import polygamma
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = float(polygamma(1, y))
        step = tri * (1.0 - tri / x) / float(polygamma(2, y))
        y += step
        if -step / y < 1e-8:
            break
    return y


def squeeze(s2: np.ndarray, df: np.ndarray) -> tuple[np.ndarray, float]:
    """Variances from few degrees of freedom, each moved towards what the
    others suggest: limma's `squeezeVar` (Smyth 2004), no trend, not
    robust. Returns the squeezed variances and the prior degrees of
    freedom the others were worth — 0 when there are too few to say,
    infinity when they do not differ more than chance would have them."""
    from scipy.special import digamma, polygamma
    s2, df = np.asarray(s2, dtype=np.float64), np.asarray(df, dtype=np.float64)
    ok = (df > 0) & np.isfinite(s2)
    if ok.sum() < 2 or not (s2[ok] > 0).any():
        return s2.copy(), 0.0
    x = np.maximum(s2[ok], 1e-5 * np.median(s2[ok][s2[ok] > 0]))
    e = np.log(x) - digamma(df[ok] / 2) + np.log(df[ok] / 2)
    spread = e.var(ddof=1) - polygamma(1, df[ok] / 2).mean()
    if spread <= 0:
        return np.full(s2.shape, float(np.exp(e.mean()))), float("inf")
    prior_df = 2.0 * _trigamma_inverse(spread)
    prior = np.exp(e.mean() + digamma(prior_df / 2) - np.log(prior_df / 2))
    return (df * np.where(ok, s2, 0.0) + prior_df * prior) / (df + prior_df), prior_df


# --- the comparison ------------------------------------------------------------------

def design(groups: list[str], contrast: tuple[str, str],
           covariate: list[str] | None = None) -> np.ndarray:
    """samples x coefficients: an intercept, 1 where the sample is of
    `contrast[0]`, and a column for each level of the covariate after its
    first."""
    columns = [np.ones(len(groups)), np.array([g == contrast[0] for g in groups], dtype=float)]
    if covariate is not None:
        for level in compare.ordered(covariate)[1:]:
            columns.append(np.array([c == level for c in covariate], dtype=float))
    return np.column_stack(columns)


def _untestable(n: tuple[int, int], X: np.ndarray, adjusted: bool) -> str:
    if min(n) == 0:
        return "one group: nothing to compare"
    if min(n) < 2:
        return "a group has fewer than 2 samples"
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return "what is adjusted for cannot be told apart from the groups"
    if X.shape[0] <= X.shape[1]:
        return "too few samples" + (" to adjust for anything" if adjusted else "")
    return ""


def _fold(y: np.ndarray, X: np.ndarray, totals: np.ndarray, dispersion) -> np.ndarray:
    """The group's coefficient in log2, from counts with the prior count."""
    prior = PRIOR_COUNT * totals / totals.mean()
    beta, _ = glm_fit(y + prior, X, np.log(totals + 2 * prior), dispersion)
    return beta[:, 1] / np.log(2.0)


def _residual_df(y: np.ndarray, mu: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Samples less coefficients, among the samples that can still say
    something: a group with no cells at all is fitted exactly and cannot
    speak of the variation (edgeR does the same)."""
    zero = _zeros(y, mu)
    full = X.shape[0] - X.shape[1]
    return np.array([full if not z.any() else
                     (~z).sum() - np.linalg.matrix_rank(X[~z]) for z in zero], dtype=float)


def _test(y: np.ndarray, X: np.ndarray, offset: np.ndarray, test: str):
    """(dispersion used for each row, statistic, p)."""
    from scipy import stats as st
    common, own = dispersions(y, X, offset)
    phi = np.full(y.shape[0], common) if test == QL else own
    null = np.delete(X, 1, axis=1)
    _, mu = glm_fit(y, X, offset, phi)
    without = deviance(y, glm_fit(y, null, offset, phi)[1], phi[:, None])
    lr = without - deviance(y, mu, phi[:, None])
    # a difference within what the fits were iterated to is no difference
    lr = np.where(lr > 1e-7 * (without + 1.0), lr, 0.0)
    if test == LRT:
        return phi, lr, st.chi2.sf(lr, 1)
    df = _residual_df(y, mu, X)
    s2 = np.divide(deviance(y, mu, phi[:, None]), df, out=np.full(df.shape, np.nan),
                   where=df > 0)
    squeezed, prior_df = squeeze(s2, df)
    with np.errstate(divide="ignore", invalid="ignore"):
        # no difference at all is p = 1 whatever the variance, zero included
        f = np.where(lr > 0, lr / squeezed, 0.0)
    total = np.minimum(df + prior_df, 1e12)     # infinity is a chi-square; scipy says NaN
    p = np.where(total > 0, st.f.sf(f, 1, np.where(total > 0, total, 1.0)), np.nan)
    return phi, f, p


def differential(counts: np.ndarray, clusters: list[str], groups: list[str],
                 contrast: tuple[str, str], covariate: list[str] | None = None,
                 test: str = QL, correction: str = compare.BH,
                 totals: np.ndarray | None = None) -> list[Result]:
    """Each cluster's share of the sample in `contrast[0]` against
    `contrast[1]`: one `Result` per row of `counts` (clusters x samples).
    `groups` and `covariate` name each sample's group and, optionally, the
    level of a factor to adjust for. `totals` are the samples' sizes where
    they are not the column sums. Samples of other groups, and samples
    with no events, are left out.

    Where no test can be made the p is NaN and the note says why; the fold
    change is then from the pooled counts (a Poisson fit)."""
    counts = np.asarray(counts, dtype=np.float64)
    totals = counts.sum(axis=0) if totals is None else np.asarray(totals, dtype=np.float64)
    keep = [i for i, g in enumerate(groups) if g in contrast and totals[i] > 0]
    y, lib = counts[:, keep], totals[keep]
    kept = [groups[i] for i in keep]
    n = (kept.count(contrast[0]), kept.count(contrast[1]))
    X = design(kept, contrast, [covariate[i] for i in keep] if covariate is not None else None)
    members = (X[:, 1] == 1, X[:, 1] == 0)
    pct = [(100 * y[:, m] / lib[m]).mean(axis=1) if m.any() else np.full(len(clusters), np.nan)
           for m in members]
    out = [Result(str(c), tuple(contrast), n, (float(pct[0][i]), float(pct[1][i])))
           for i, c in enumerate(clusters)]
    reason = _untestable(n, X, covariate is not None)
    if min(n) == 0:
        return _noted(out, reason)
    present = np.flatnonzero(y.sum(axis=1) > 0)
    for i in np.flatnonzero(y.sum(axis=1) == 0):
        out[i].note = "no cells in either group"
    if present.size == 0:
        return _noted(out, reason)
    if reason:
        for i, fold in zip(present, _fold(y[present], X[:, :2], lib, 0.0)):
            out[i].log2_fold = float(fold)
        return _noted(out, reason)
    phi, statistic, p = _test(y[present], X, np.log(lib), test)
    fold = _fold(y[present], X, lib, phi)
    for k, i in enumerate(present):
        r = out[i]
        r.log2_fold, r.dispersion = float(fold[k]), float(phi[k])
        r.test, r.statistic, r.p = TESTS[test], float(statistic[k]), float(p[k])
        empty = [g for g, m in zip(contrast, members) if y[i][m].sum() == 0]
        if empty:
            r.note = (f"no cells in {empty[0] or '—'}: the fold change is held finite by a "
                      f"prior count of {PRIOR_COUNT}")
    return adjusted(out, correction)


def _noted(results: list[Result], note: str) -> list[Result]:
    for r in results:
        r.note = r.note or note
    return results


def adjusted(results: list[Result], correction: str = compare.BH) -> list[Result]:
    """The same results with p adjusted over all of them."""
    for r, p in zip(results, compare.adjust([r.p for r in results], correction)):
        r.p_adjusted = float(p)
    return results


def within(counts: np.ndarray, clusters: list[str], groups: list[str],
           contrast: tuple[str, str], panels: list[str],
           covariate: list[str] | None = None, test: str = QL,
           correction: str = compare.BH) -> list[Result]:
    """The comparison made separately among the samples of each value of
    `panels` — within each day, say — model, dispersion and all, with the
    p values then adjusted over every panel's tests together."""
    counts = np.asarray(counts, dtype=np.float64)
    out: list[Result] = []
    for panel in compare.ordered(panels):
        keep = [i for i, p in enumerate(panels) if p == panel]
        found = differential(counts[:, keep], clusters, [groups[i] for i in keep], contrast,
                             [covariate[i] for i in keep] if covariate is not None else None,
                             test, compare.NONE)
        for r in found:
            r.panel = panel
        out += found
    return adjusted(out, correction)


def to_table(results: list[Result], correction: str = compare.BH,
             panel: str | None = None) -> StatsTable:
    """Every cluster as a row; `panel` names the first column where the
    results are of several panels."""
    adjusted_name = "p" if correction == compare.NONE else f"p adjusted ({correction})"
    a, b = results[0].groups if results else ("", "")
    table = StatsTable(([panel] if panel else [])
                       + ["Cluster", f"{a or '—'} (mean %)", f"{b or '—'} (mean %)",
                          "log2 fold change", "p", adjusted_name, "", "n", "Test", "Statistic",
                          "Dispersion", "Note"])
    for r in results:
        table.rows.append(([r.panel] if panel else [])
                          + [r.cluster, r.percent[0], r.percent[1], r.log2_fold, r.p,
                             r.p_adjusted, compare.stars(r.p_adjusted), f"{r.n[0]} vs {r.n[1]}",
                             r.test, r.statistic, r.dispersion, r.note])
    return table
