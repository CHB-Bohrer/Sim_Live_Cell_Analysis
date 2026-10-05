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


def test_insulation_is_one_for_a_map_with_only_distance_decay():
    """Regression: on a raw map the within/across ratio is >> 1 even with NO boundary, because 'within' pixels are
    closer to the diagonal. On observed/expected it must be ~1."""
    nb = 200
    i, j = np.indices((nb, nb))
    m = 1.0 / (1.0 + np.abs(i - j)) ** 1.2          # pure power-law decay, no domains
    raw = A.boundary_insulation(m, [1000], bin_size=10, window_kb=150)[0]["ratio"]
    oe = A.boundary_insulation(m, [1000], bin_size=10, window_kb=150, use_oe=True)[0]["ratio"]
    assert raw > 1.5, "raw ratio is biased by distance decay (this is why use_oe exists)"
    assert abs(oe - 1.0) < 0.05


def test_loop_dot_enrichment_finds_designed_loop_and_not_empty_positions():
    oe = np.ones((200, 200))
    oe[49:52, 79:82] = 5.0        # a dot spans a few bins; the metric averages the 3x3 centre
    oe[79:82, 49:52] = 5.0
    on = A.loop_dot_enrichment(oe, [(500, 800)], bin_size=10)
    off = A.loop_dot_enrichment(oe, [(1000, 1300)], bin_size=10)
    assert on[0] > 1.5 and abs(off[0] - 1.0) < 1e-9


def test_across_insulation_is_one_without_boundaries_and_detects_one_with():
    nb = 200
    i, j = np.indices((nb, nb))
    m = 1.0 / (1.0 + np.abs(i - j)) ** 1.2            # distance decay only
    assert abs(A.across_insulation(m, [1000], bin_size=10)[0] - 1.0) < 0.05
    m2 = m.copy()
    m2[85:100, 100:115] *= 0.4                          # contacts across position 1000 reduced 2.5-fold
    m2[100:115, 85:100] *= 0.4
    assert A.across_insulation(m2, [1000], bin_size=10)[0] > 1.5


def test_anchor_contact_enrichment_detects_a_held_loop():
    rng = np.random.default_rng(0)
    T = 1000
    beads = np.array([0, 10, 20])
    traj = rng.normal(0, 20, size=(T, 3, 3)).astype(np.float32)     # far apart most of the time
    held = rng.random(T) < 0.25
    traj[held, 2] = traj[held, 0] + 0.7                              # anchors together 25% of the time
    s = np.array([5.0, 50.0]); p = np.array([0.05, 0.01])           # P(s=20) is ~0.03
    e = A.anchor_contact_enrichment(traj, beads, [(0, 20)], (s, p))
    assert e[0] > 5


def test_pooled_anchor_contact_is_not_dominated_by_an_unformed_loop():
    rng = np.random.default_rng(1)
    T, beads = 2000, np.array([0, 20, 40, 60])
    traj = rng.normal(0, 20, size=(T, 4, 3)).astype(np.float32)
    held = rng.random(T) < 0.3
    traj[held, 1] = traj[held, 0] + 0.7        # loop (0, 20) formed 30% of the time; loop (40, 60) never forms
    s = np.array([5.0, 100.0]); p = np.array([0.02, 0.01])
    obs, exp = A.anchor_contact_pooled(traj, beads, [(0, 20), (40, 60)], (s, p))
    assert abs(obs - 0.15) < 0.03 and obs / exp > 8   # a per-loop median would be ~0.15 vs 0.3: pooled is stable


def test_pooled_insulation_fold_matches_a_synthetic_fold():
    nb = 200
    i, j = np.indices((nb, nb))
    m = 1.0 / (1.0 + np.abs(i - j)) ** 1.2
    m[85:100, 100:115] *= 0.5
    m[100:115, 85:100] *= 0.5
    assert abs(A.insulation_fold_pooled(m, [1000], 10) - 2.0) < 0.15
