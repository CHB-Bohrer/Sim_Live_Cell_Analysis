"""Tests of the chromatin analysis tools on synthetic data with known answers (no GPU)."""
import numpy as np

from simlive.stage1_chromatin import analysis as A


def random_walk(n, seed=0, step=1.0):
    rng = np.random.default_rng(seed)
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return np.cumsum(d * step, axis=0).astype(np.float32)


def test_contact_map_is_symmetric_and_counts_adjacent_monomers():
    confs = [random_walk(600, s) for s in range(3)]
    m = A.contact_map(confs, bin_size=10, cutoff=3.0)
    assert m.shape == (60, 60) and np.allclose(m, m.T)
    assert np.diag(m).min() > 0                      # neighbours along the chain are always in contact


def test_contact_probability_decays_with_separation_for_a_random_walk():
    confs = [random_walk(3000, s) for s in range(4)]
    s, p = A.contact_probability(A.separation_histogram(confs, cutoff=3.0))
    assert p[0] > p[len(p) // 2] > p[-1]
    slope = A.loglog_slope(s, p, 20, 800)
    assert -2.0 < slope < -0.9, slope                # ideal chain: contact probability ~ s^-1.5


def test_msd_of_brownian_beads_grows_linearly():
    rng = np.random.default_rng(1)
    traj = np.cumsum(rng.normal(size=(2000, 50, 3)), axis=0).astype(np.float32)
    lags = np.array([1, 2, 4, 8, 16, 32, 64])
    m = A.msd(traj, lags)
    assert abs(A.msd_exponent(lags, m, 1, 64) - 1.0) < 0.05
    assert abs(m[0] - 3.0) < 0.2                      # 3 dimensions x variance 1


def test_insulation_ratio_detects_a_synthetic_boundary():
    nb, b = 100, 50
    m = np.ones((nb, nb))
    m[:b, b:] = 0.4
    m[b:, :b] = 0.4                                   # contacts across the boundary are 2.5-fold lower
    r = A.boundary_insulation(m, [b * 10], bin_size=10, window_kb=150)
    assert abs(r[0]["ratio"] - 2.5) < 1e-6


def test_corner_peak_enrichment_detects_a_synthetic_dot():
    oe = np.ones((120, 120))
    oe[30, 80] = oe[80, 30] = 4.0
    e = A.corner_peak_enrichment(oe, [300, 800], bin_size=10)
    assert len(e) == 1 and e[0] > 1.2
