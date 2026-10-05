"""Unit tests for the 1D loop-extrusion engine (no GPU needed)."""
import numpy as np

from simlive.stage1_chromatin.lef1d import LEF1D, synthetic_tad_boundaries


def test_lefs_never_overlap_and_stay_ordered():
    lef = LEF1D(2000, 40, lifetime_steps=300, seed=1)
    for _ in range(50):
        lef.step(20)
        assert (lef.left < lef.right).all()
        legs = np.concatenate([lef.left, lef.right])
        assert len(np.unique(legs)) == len(legs), "two LEF legs on the same monomer"
        assert legs.min() >= 0 and legs.max() < 2000


def test_free_lef_processivity_is_about_two_times_speed_times_lifetime():
    # one isolated LEF on a huge chain: mean loop size at unbinding = 2 * lifetime (speed 1 monomer/step/leg)
    lifetime = 100
    n_trials, sizes = 400, []
    for s in range(n_trials):
        lef = LEF1D(20000, 1, lifetime_steps=lifetime, seed=s)
        age = 0
        while True:
            before = lef.loop_sizes()[0]
            lef.step(1)
            age += 1
            if lef.loop_sizes()[0] < before or age > 5000:  # re-placed (loop size reset): the previous loop ended
                sizes.append(before)
                break
    assert abs(np.mean(sizes) - 2 * lifetime) < 0.25 * 2 * lifetime


def test_convergent_ctcf_pair_holds_a_persistent_loop():
    sites = [{"pos": 400, "blocks_left": 1.0, "blocks_right": 0.0},    # left anchor stops left-moving legs
             {"pos": 600, "blocks_left": 0.0, "blocks_right": 1.0}]    # right anchor stops right-moving legs
    lef = LEF1D(1000, 6, lifetime_steps=2000, ctcf=sites, seed=3)
    lef.step(3000)
    b = lef.bonds()
    held = ((b[:, 0] == 400) & (b[:, 1] == 600)).sum()
    assert held >= 1, f"expected a LEF anchored at (400, 600); bonds: {b.tolist()}"


def test_impermeable_boundary_is_never_crossed_by_a_leg_starting_inside():
    sites = [{"pos": 500, "blocks_left": 1.0, "blocks_right": 1.0}]
    lef = LEF1D(1000, 10, lifetime_steps=10 ** 9, ctcf=sites, seed=5)   # no unbinding: legs only drift
    # keep only LEFs that started on the right side of the boundary; their left legs must never get below 500
    inside = lef.left > 500
    lef.step(5000)
    assert (lef.left[inside] >= 500).all()


def test_synthetic_boundaries_are_ordered_and_spaced():
    sites = synthetic_tad_boundaries(46710, 800.0, np.random.default_rng(0), stall=0.9)
    pos = np.array([s["pos"] for s in sites])
    assert (np.diff(pos) >= 150).all() and 30 < len(pos) < 120
