"""Probe binding model: the exact simulation must agree with the exact theory (no GPU needed)."""
import numpy as np
import pandas as pd
import pytest

from simlive.stage1_chromatin import probes as P


def test_simulation_matches_theory():
    p = P.ProbeParams(n_probes=40, k_bind=0.1, k_scan=0.4, k_unbind=0.2, k_off=0.15)
    th = P.summary(p)
    tr = P.simulate_trace(p, 60000.0, seed=3)
    d = P.sample_trace(tr, np.arange(0, 60000.0, 1.0), 1e-3)
    a = d.n_attached.to_numpy()
    assert abs(a.mean() - th["mean_attached"]) < 0.03 * th["mean_attached"]
    assert abs(a.std() / a.mean() - th["cv_attached"]) < 0.1 * th["cv_attached"]
    assert abs(d.n_scanning.mean() - th["mean_scanning"]) < 0.04 * th["mean_scanning"]
    assert (d.n_free + d.n_bound + d.n_scanning == p.n_probes).all()          # probes are conserved


def test_dwell_time_of_a_single_probe():
    """One probe, no abort: attached time per visit = 1/k_scan + 1/k_off (mean), measured from the event list."""
    p = P.ProbeParams(n_probes=1, k_bind=1.0, k_scan=2.0, k_unbind=0.0, k_off=0.5)
    t, F, B, S = P.simulate_trace(p, 40000.0, seed=1)
    att = (B + S) > 0
    starts, ends = t[1:][att[1:] & ~att[:-1]], t[1:][~att[1:] & att[:-1]]
    ends = ends[ends > starts[0]]
    n = min(len(starts), len(ends))
    assert abs((ends[:n] - starts[:n]).mean() - (1 / 2.0 + 1 / 0.5)) < 0.1
    assert abs(P.summary(p)["residence_s"] - 2.5) < 1e-12


def test_capacity_limits_and_saturates():
    p = P.ProbeParams(n_probes=50, k_bind=5.0, k_scan=0.5, k_unbind=0.0, k_off=0.05, capacity=8)
    t, F, B, S = P.simulate_trace(p, 5000.0, seed=2)
    assert (B + S).max() <= 8 and (B + S).mean() > 6                            # strong binding: locus nearly full


def test_exposure_average_and_frame_sampling():
    tr = (np.array([0.0, 1.0, 3.0]), np.array([2, 1, 0]), np.array([0, 1, 2]), np.array([0, 0, 0]))   # attached 0,1,2
    d = P.sample_trace(tr, [0.0, 1.0, 3.5], 2.0)
    assert list(d.n_attached) == [0, 1, 2]
    assert np.allclose(d.mean_attached, [0.5, 1.0, 2.0])      # [0,2]: (0*1+1*1)/2; [1,3]: 1 throughout; [3.5,5.5]: 2 throughout


def test_occupancy_table_reproducible_and_independent_per_locus():
    cells = pd.DataFrame({"t": np.tile(np.arange(5), 2), "cell_id": np.repeat([1, 2], 5)})
    cfg = {"acquisition": {"frame_interval_s": 10, "exposure_s": 0.1}, "loci": {"n_loci": 2}}
    a = P.simulate_probe_occupancy(cells, cfg, np.random.default_rng(1))
    b = P.simulate_probe_occupancy(cells, cfg, np.random.default_rng(1))
    pd.testing.assert_frame_equal(a, b)
    assert len(a) == 2 * 5 * 2 and set(a.locus_name) == {"promoter", "enhancer"}
    assert not np.array_equal(a[a.locus_id == 0].n_attached.to_numpy(), a[a.locus_id == 1].n_attached.to_numpy())


def test_bad_parameters_rejected():
    with pytest.raises(ValueError):
        P.ProbeParams(k_off=0.0)
    with pytest.raises(ValueError):
        P.ProbeParams(capacity=0)


def test_sweep_has_theory_and_simulation():
    df = P.sweep(P.DEFAULTS["promoter"], {"k_off": [0.05, 0.2]}, reps=2, duration_s=4000.0, seed=1)
    assert len(df) == 4 and {"th_mean_attached", "sim_mean_attached", "k_off"} <= set(df.columns)
    m = df.groupby("k_off").sim_mean_attached.mean()
    assert m[0.05] > m[0.2]                                                  # longer scanning -> more probes attached


def test_locus_brightness_follows_attached_probes():
    from simlive.io.runs import REPO_ROOT, load_config
    from simlive.stage3_microscopy.render import render_loci
    cfg = load_config(REPO_ROOT / "configs" / "loci_demo.yaml", ["geometry.fov_um=[10, 10]", "acquisition.n_frames=2",
                                                                "optics.read_noise_e=0", "optics.background_photons=0",
                                                                "loci.bleach_tau_s=null", "loci.probes.enabled=true"])
    loci = pd.DataFrame({"t": [0, 1], "cell_id": [1, 1], "locus_id": [0, 0], "y_um": [5.0, 5.0], "x_um": [5.0, 5.0]})
    occ = pd.DataFrame({"t": [0, 1], "cell_id": [1, 1], "locus_id": [0, 0], "mean_attached": [10.0, 20.0]})
    img = render_loci(loci, cfg, np.random.default_rng(0), 0, occupancy=occ)
    signal = (img - cfg["optics"]["offset_adu"]).reshape(2, -1).sum(1)
    assert abs(signal[1] / signal[0] - 2.0) < 0.1                       # shot noise only (~1% at these counts)
    expect = 300 * cfg["acquisition"]["exposure_s"] * 10 * cfg["optics"]["quantum_efficiency"]
    assert abs(signal[0] - expect) < 0.1 * expect
