"""E-P distance -> P(ON) curve recovery (stage7_validation/ep_curve.py)."""
import numpy as np

from simlive.stage7_validation import ep_curve as ep


def test_kernel_columns_sum_to_one_and_noise_free_is_a_histogram():
    for dims in (2, 3):
        K = ep.kernel_matrix(40.0, dims, 25.0)
        assert np.allclose(K.sum(0), 1.0) and (K >= 0).all()
    K0 = ep.kernel_matrix(9.0, 3, 0.0)
    g = np.argmin(np.abs(ep.D_GRID - 410.0))
    b = np.argmax(K0[:, g])
    assert ep.EDGES[b] <= 410.0 <= ep.EDGES[b + 1]
    # 2D projection shortens distances: mean measured distance of a true 500 nm pair is 500 x pi/4 (isotropic) without noise
    K2 = ep.kernel_matrix(9.0, 2, 0.0)
    g5 = np.argmin(np.abs(ep.D_GRID - 500.0))
    centres = 0.5 * (ep.EDGES[:-1] + ep.EDGES[1:])
    assert abs((K2[:, g5] * centres).sum() - 500.0 * np.pi / 4) < 12.0


def test_measure_distance_matches_the_kernel():
    rng = np.random.default_rng(0)
    d = np.full(200000, 600.0)
    for dims, shift in ((3, 0.0), (2, 0.0), (3, 80.0)):
        m = ep.measure_distance(d, ep.Measurement(sigma_nm=60.0, dims=dims, shift_nm=shift), rng)
        h = np.histogram(np.clip(m, 0, 2999.9), bins=ep.EDGES)[0] / len(m)
        K = ep.kernel_matrix(60.0, dims, shift)
        g = np.argmin(np.abs(ep.D_GRID - 600.0))
        assert np.abs(h - K[:, g]).max() < 0.004, (dims, shift)


def test_distance_process_has_the_requested_contact_fraction():
    rng = np.random.default_rng(1)
    d = ep.ep_distance(400000, 5.0, rng, 400.0, 0.10, 240.0, 30.0)
    assert abs((d < 400.0).mean() - 0.10) < 0.01


def test_fit_recovers_a_known_curve_through_2d_noise_and_beats_the_naive_fit():
    rng = np.random.default_rng(2)
    n = 300000
    d = np.where(rng.random(n) < 0.12, rng.uniform(100, 380, n), 400 * 2 * np.exp(0.35 * rng.normal(size=n)))
    truth = (0.05, 8.0, 350.0)
    y = (rng.random(n) < ep.curve(d, *truth, 4.0)).astype(float)
    cells = rng.integers(0, 300, n)
    meas = ep.Measurement(sigma_nm=70.0, dims=2)
    dm = ep.measure_distance(d, meas, rng)
    agg = ep.Aggregates.build(cells, dm, y)
    K = ep.kernel_matrix(70.0, 2, 0.0)
    prior = ep.population_prior(agg.n.sum(0), K)
    f = ep.fit_curve(agg, K, prior)
    nv = ep.fit_curve(agg, K, prior, naive=True)
    assert abs(f.fold / 8.0 - 1) < 0.2 and abs(f.d_half / 350.0 - 1) < 0.15
    grid = ep.D_GRID
    err = lambda q: np.abs(q.g(grid) - ep.curve(grid, *truth, 4.0))[(grid > 100) & (grid < 1000)].mean()   # noqa: E731
    assert err(f) < err(nv)


def test_null_population_gives_no_coupling_and_the_ci_contains_one():
    s = ep.EPSetup(n_cells=300)
    pop = ep.simulate_population(s, 5, coupled=False, with_trace=False)
    a = ep.analyse(pop, ep.Measurement(sigma_nm=50.0, dims=3), n_boot=60, seed=1)
    lo, hi = a["contrast_ci"]
    assert lo <= 1.0 <= hi, (lo, hi)


def test_attenuation_correction_undoes_a_known_sensitivity_and_false_positive_rate():
    s = ep.EPSetup(n_cells=500)
    pop = ep.simulate_population(s, 6, with_trace=False)
    rng = np.random.default_rng(3)
    se, fp = 0.7, 0.03
    y = np.where(pop.on_true == 1, rng.random(len(pop.on_true)) < se, rng.random(len(pop.on_true)) < fp).astype(float)
    m = ep.Measurement(sigma_nm=30.0, dims=3)
    raw = ep.analyse(pop, m, y=y)["fit"]
    cor = ep.analyse(pop, m, y=y, correct_calls=(se, fp))["fit"]
    ref = ep.analyse(pop, m)["fit"]
    assert abs(np.log(cor.fold / ref.fold)) < abs(np.log(raw.fold / ref.fold))
    assert abs(cor.p_far / ref.p_far - 1) < 0.25


def test_wrong_picks_are_corrected_when_their_share_is_supplied():
    # the same curve as above, measured with 20% wrong picks (distance uniform in 0-2 um whatever the truth)
    rng = np.random.default_rng(4)
    n = 300000
    d = np.where(rng.random(n) < 0.12, rng.uniform(100, 380, n), 400 * 2 * np.exp(0.35 * rng.normal(size=n)))
    truth = (0.05, 8.0, 350.0)
    y = (rng.random(n) < ep.curve(d, *truth, 4.0)).astype(float)
    cells = rng.integers(0, 300, n)
    dm = ep.measure_distance(d, ep.Measurement(sigma_nm=70.0, dims=2, outlier_rate=0.2), rng)
    agg = ep.Aggregates.build(cells, dm, y)
    K = ep.kernel_matrix(70.0, 2, 0.0)
    assert np.allclose(ep.with_outliers(K, 0.2).sum(0), 1.0)               # still a distribution over measured bins
    fits = {}
    for rate in (0.0, 0.2):
        Km = ep.with_outliers(K, rate)
        fits[rate] = ep.fit_curve(agg, Km, ep.population_prior(agg.n.sum(0), Km))
    want = ep.curve(150.0, *truth, 4.0) / ep.curve(1000.0, *truth, 4.0)
    assert fits[0.0].contrast() < 0.85 * want                               # uncorrected: the effect is diluted
    assert abs(fits[0.2].contrast() / want - 1) < 0.1                         # corrected: recovered
