import numpy as np
import pytest

from milc import normalize as nz

# two populations, apart on the first two markers and alike on the third
_A = np.array([1.0, 1.0, 2.0])
_B = np.array([3.0, 2.5, 2.0])
# per batch: what the instrument did to each marker that day
_EFFECT = {
    "d1": (np.array([1.0, 1.0, 1.0]), np.array([0.0, 0.0, 0.0])),
    "d2": (np.array([1.15, 0.9, 1.1]), np.array([0.4, -0.3, 0.25])),
    "d3": (np.array([0.9, 1.1, 1.0]), np.array([-0.3, 0.2, -0.4])),
}


def _mix(rng, n, frac_a):
    """events x 3 and, per event, whether it is of population A."""
    is_a = rng.random(n) < frac_a
    X = rng.normal(0, 0.25, (n, 3)) + np.where(is_a[:, None], _A, _B)
    return X, is_a


def _batches(seed, n=4000, frac_a=0.5):
    rng = np.random.default_rng(seed)
    out, truth = {}, {}
    for b, (scale, shift) in _EFFECT.items():
        X, truth[b] = _mix(rng, n, frac_a)
        out[b] = X * scale + shift
    return out, truth


def test_batches_coincide_after_normalisation():
    ref, _ = _batches(1)
    model = nz.fit(ref)
    assert model.batches == ["d1", "d2", "d3"] and model.goal == "mean"
    assert (model.n_markers, model.n_clusters, model.som) == (3, 1, None)
    assert model.quantiles.shape == (3, 1, 99, 3) and model.fallbacks == []
    before = nz.distance(ref)
    after = nz.distance({b: nz.apply(model, X, b) for b, X in ref.items()})
    print("reference: distance before", before.round(4), "after", after.round(5))
    assert np.all(before > 0.15) and np.all(after < before / 100)


def test_other_samples_lose_the_batch_and_keep_the_biology():
    ref, _ = _batches(1)
    model = nz.fit(ref)
    # samples the model never saw, with population A enriched
    new, truth = _batches(2, n=3000, frac_a=0.8)
    norm = {b: nz.apply(model, X, b) for b, X in new.items()}
    before, after = nz.distance(new), nz.distance(norm)
    print("other samples: distance before", before.round(4), "after", after.round(4))
    assert np.all(after < before / 8)
    for pop in (True, False):
        means = np.array([norm[b][truth[b] == pop].mean(axis=0) for b in norm])
        assert np.ptp(means, axis=0).max() < 0.06
    # the proportions are the samples' own, not the reference's: the two
    # populations are still where a threshold between them finds 80%
    cut = (norm["d1"][truth["d1"]][:, 0].mean() + norm["d1"][~truth["d1"]][:, 0].mean()) / 2
    for b in norm:
        assert abs((norm[b][:, 0] < cut).mean() - 0.8) < 0.02


def _population_specific(seed, frac_a=0.5, n=4000):
    """Two batches; in the second, population A alone is 0.8 brighter on
    the third marker, where A and B otherwise lie on top of each other."""
    rng = np.random.default_rng(seed)
    out, truth = {}, {}
    for b, shift in (("d1", 0.0), ("d2", 0.8)):
        X, truth[b] = _mix(rng, n, frac_a)
        X[truth[b], 2] += shift
        out[b] = X
    return out, truth


def _gap(norm, truth, pop):
    """Between the batches, the third marker's mean in one population."""
    return abs(norm["d1"][truth["d1"] == pop][:, 2].mean()
               - norm["d2"][truth["d2"] == pop][:, 2].mean())


def test_an_effect_on_one_population_needs_clusters():
    ref, _ = _population_specific(3)
    new, truth = _population_specific(4, frac_a=0.7, n=3000)
    whole = nz.fit(ref, n_clusters=1)
    split = nz.fit(ref, n_clusters=2, grid=5, seed=1)
    assert split.n_clusters == 2 and split.som is not None and split.fallbacks == []
    assert split.quantiles.shape == (2, 2, 99, 3)
    # the map's two clusters are the two populations
    assert (split.som.map(new["d1"]) == split.som.map(new["d1"][:1])[0]).mean() \
        == pytest.approx(truth["d1"].mean() if truth["d1"][0] else 1 - truth["d1"].mean(),
                         abs=0.01)
    for model, ok in ((whole, False), (split, True)):
        norm = {b: nz.apply(model, X, b) for b, X in new.items()}
        gaps = [_gap(norm, truth, True), _gap(norm, truth, False)]
        print("clusters", model.n_clusters, "gap left in A, B:", np.round(gaps, 3))
        assert (max(gaps) < 0.05) == ok
        if not ok:
            assert min(gaps) > 0.2  # one curve for both: neither is right


def test_goal_batch_is_left_unchanged():
    ref, _ = _batches(5)
    new, _ = _batches(6, n=1000, frac_a=0.3)
    for k in (1, 2):
        model = nz.fit(ref, n_clusters=k, goal="d2", grid=5)
        assert model.goal == "d2"
        assert np.allclose(nz.apply(model, new["d2"], "d2"), new["d2"], atol=1e-9)
        norm = {b: nz.apply(model, X, b) for b, X in ref.items()}
        assert np.all(nz.distance(norm, goal="d2") < nz.distance(ref, goal="d2") / 20)
    with pytest.raises(KeyError, match="d9"):
        nz.fit(ref, goal="d9")


def test_identical_batches_are_left_alone():
    X, _ = _mix(np.random.default_rng(7), 3000, 0.5)
    ref = {"a": X, "b": X.copy(), "c": X.copy()}
    assert np.allclose(nz.distance(ref), 0)
    wide = np.random.default_rng(8).uniform(-2, 6, (500, 3))  # beyond the reference too
    for k in (1, 2):
        model = nz.fit(ref, n_clusters=k, grid=5)
        assert np.allclose(nz.apply(model, wide, "b"), wide, atol=1e-9)


def test_monotone_and_a_shift_beyond_the_reference():
    ref, _ = _batches(9)
    model = nz.fit(ref)
    x = np.linspace(-20, 30, 5001)
    X = np.column_stack([x, x, x])
    for b in ref:
        Y = nz.apply(model, X, b)
        assert np.all(np.diff(Y, axis=0) >= 0)
        # slope 1 past the outermost quantiles, at both ends
        low, high = x < ref[b].min(), x > ref[b].max()
        assert np.allclose(np.diff(Y[low], axis=0), np.diff(X[low], axis=0))
        assert np.allclose(np.diff(Y[high], axis=0), np.diff(X[high], axis=0))


def test_tied_values_do_not_break_the_curve():
    rng = np.random.default_rng(10)
    ref = {b: np.maximum(rng.normal(m, 1, (2000, 1)), 0) for b, m in (("a", 0.5), ("b", 1.0))}
    model = nz.fit(ref)  # a third of the events sit at 0 exactly
    Y = nz.apply(model, np.linspace(-1, 5, 400)[:, None], "a")
    assert np.isfinite(Y).all() and np.all(np.diff(Y[:, 0]) >= 0)


def test_shape_and_dtype():
    ref, _ = _batches(11, n=1000)
    model = nz.fit(ref, n_clusters=2, grid=4)
    X32 = ref["d2"][:137].astype(np.float32)
    Y = nz.apply(model, X32, "d2")
    assert Y.shape == X32.shape and Y.dtype == np.float32
    assert nz.apply(model, ref["d2"][:5], "d2").dtype == np.float64
    assert nz.apply(model, np.zeros((0, 3)), "d2").shape == (0, 3)
    assert not np.shares_memory(Y, X32)
    with pytest.raises(ValueError):
        nz.apply(model, np.zeros((4, 2)), "d2")


def test_too_few_events_fall_back_to_identity():
    ref, _ = _batches(12, n=2000)
    ref["d3"] = ref["d3"][:30]
    model = nz.fit(ref)
    assert model.fallbacks == [("d3", 1, 30)]
    assert np.isnan(model.quantiles[2]).all()
    assert np.array_equal(nz.apply(model, ref["d3"], "d3"), ref["d3"])
    # the others still meet, on a goal the short batch had no say in
    assert np.allclose(model.goal_quantiles, model.quantiles[:2].mean(axis=0))
    rest = {b: ref[b] for b in ("d1", "d2")}
    norm = {b: nz.apply(model, X, b) for b, X in rest.items()}
    assert np.all(nz.distance(norm) < nz.distance(rest) / 50)

    # a threshold of one's own
    assert nz.fit(ref, min_events=20).fallbacks == []

    # the goal batch itself too short: nothing to aim at, nothing moves
    model = nz.fit(ref, goal="d3")
    assert model.fallbacks == [("d1", 1, 2000), ("d2", 1, 2000), ("d3", 1, 30)]
    assert np.array_equal(nz.apply(model, ref["d1"], "d1"), ref["d1"])


def test_a_cluster_one_batch_lacks():
    """Per cluster, not per batch: a batch without one population is
    corrected in the other."""
    ref, truth = _batches(13, n=2000)
    ref["d2"] = np.vstack([ref["d2"][truth["d2"]], ref["d2"][~truth["d2"]][:10]])
    model = nz.fit(ref, n_clusters=2, grid=5)
    assert len(model.fallbacks) == 1
    batch, cluster, n = model.fallbacks[0]
    assert batch == "d2" and 5 <= n <= 15
    Y = nz.apply(model, ref["d2"], "d2")
    lacking = model.som.map(ref["d2"]) == cluster
    assert np.array_equal(Y[lacking], ref["d2"][lacking])
    assert not np.allclose(Y[~lacking], ref["d2"][~lacking])


def test_unknown_batch_and_bad_input():
    ref, _ = _batches(14, n=500)
    model = nz.fit(ref)
    with pytest.raises(KeyError, match="d7"):
        nz.apply(model, ref["d1"], "d7")
    with pytest.raises(ValueError, match="markers"):
        nz.fit({"a": ref["d1"], "b": ref["d2"][:, :2]})
    with pytest.raises(ValueError, match="no reference events"):
        nz.fit({"a": ref["d1"], "b": ref["d2"][:0]})


def test_same_seed_same_model():
    ref, _ = _batches(15, n=1500)
    a = nz.fit(ref, n_clusters=2, grid=5, seed=3)
    b = nz.fit(ref, n_clusters=2, grid=5, seed=3)
    assert np.array_equal(a.quantiles, b.quantiles)
    assert np.array_equal(a.som.codes, b.som.codes)
    assert np.array_equal(nz.apply(a, ref["d2"], "d2"), nz.apply(b, ref["d2"], "d2"))
