"""E-P curve with polymer distances: calibration, acquisition model, history estimator (stage7_validation/ep_polymer.py, ep_acquisition.py, ep_history.py).
Needs the saved validation simulations ($SIMLIVE_DATA/chromatin/validation); skipped without them."""
import numpy as np
import pytest

from simlive.stage7_validation import ep_acquisition as aq
from simlive.stage7_validation import ep_curve as ep
from simlive.stage7_validation import ep_history as eh
from simlive.stage7_validation import ep_polymer as pp

pytestmark = pytest.mark.skipif(not (pp.VAL_DIR / pp.SIMS[0] / "tracked_positions.npy").exists(), reason="validation simulations not on this machine")


@pytest.fixture(scope="module")
def pool():
    return pp.load_pool(sims=pp.SIMS[:4])


def test_calibration_hits_the_contact_fraction_and_the_msd(pool):
    cal = pp.calibrate(pool, working_nm=400.0, contact_fraction=0.10, D_um2_s=0.002, lag_s=60.0)
    d = pool.distances() * cal.nm_per_unit
    assert abs((d < 400.0).mean() - 0.10) < 1e-3
    assert abs(float(cal.msd_single_s(60.0)) / (6 * 0.002 * 60 * 1e6) - 1) < 0.02            # MSD(60 s) = 6 D t
    assert 0.4 < cal.alpha < 0.7 and cal.block_s > 0
    other = pp.calibrate(pool, applies_to="pair")
    assert other.block_s > cal.block_s                                                         # the pair moves faster than one locus: a longer block for the same D


def test_population_has_isotropic_vectors_and_the_right_contact_fraction(pool):
    cal = pp.calibrate(pool)
    pop = pp.simulate_population(ep.EPSetup(n_cells=60), pool, cal, 3, with_trace=False)
    v = pop.vec_nm
    assert abs(np.mean(v[:, 2] ** 2) / np.mean((v ** 2).sum(1)) * 3 - 1) < 0.06                         # no preferred direction
    assert abs(np.mean(pop.d_true < 400.0) - 0.10) < 0.02
    assert pop.blocks.shape[0] == 60 and np.allclose(np.linalg.norm(pop.blocks[0, 0]), np.linalg.norm(pop.vec_nm[0]))


def test_acquisition_error_budget_is_ordered_sensibly(pool):
    cal = pp.calibrate(pool)
    pop = pp.simulate_population(ep.EPSetup(n_cells=80), pool, cal, 4, with_trace=False)
    m = {k: aq.measure(pop, pop.blocks, aq.Acquisition(kind=k), seed=1) for k in ("2d", "zstack", "biplane")}
    assert m["2d"].dims == 2 and m["zstack"].dims == 3
    assert m["zstack"].parts["mean_time_offset_s"] > 1.0 > m["biplane"].parts["mean_time_offset_s"]       # a z-stack takes seconds, the biplane one moment
    assert m["zstack"].parts["exposure_blur_nm"] > m["biplane"].parts["exposure_blur_nm"]
    assert m["zstack"].parts["piezo_z_nm"] > 0 == m["biplane"].parts["piezo_z_nm"]
    rich = aq.measure(pop, pop.blocks, aq.Acquisition(kind="biplane", photon_scale=10.0), seed=1)
    assert rich.parts["localization_z_nm"] < m["biplane"].parts["localization_z_nm"]                   # more photons, better z
    simult = aq.measure(pop, pop.blocks, aq.Acquisition(kind="biplane", colour_delay_s=0.0), seed=1)
    assert simult.parts["colour_delay_nm"] == 0.0 and simult.sigma_xy < m["biplane"].sigma_xy


def test_kernel_aniso_matches_the_analytic_kernel_when_isotropic():
    K = aq.kernel_aniso(60.0, 60.0, 3, 0.0)
    Ka = ep.kernel_matrix(60.0, 3, 0.0)
    cols = ep.D_GRID <= 2500.0                      # beyond the last bin edge the two kernels treat the overflow differently
    assert np.abs(K - Ka)[:, cols].sum(0).max() < 0.35 and np.abs(K - Ka)[:, cols].sum(0).mean() < 0.1         # Monte Carlo (30k draws) against the analytic one: L1 distance per column


def test_history_fit_recovers_a_slow_response_to_fast_fluctuations(pool):
    """With the polymer's fast distance fluctuations the instantaneous curve is nearly flat but the history fit finds the promoter's ~45 s memory and a large effect."""
    cal = pp.calibrate(pool)
    pop = pp.simulate_population(ep.EPSetup(n_cells=500), pool, cal, 6, with_trace=False, keep_blocks=False)
    K0 = aq.kernel_aniso(9.0, 9.0, 3, 0.0)
    h = eh.analyse_history(pop, pop.d_true, K0)
    inst = ep.analyse_measured(pop, pop.d_true, K0)
    assert h["fit"].tau_s in (20.0, 45.0, 90.0) and h["contrast"] > 2.0 * inst["contrast"] > 0
    null = pp.simulate_population(ep.EPSetup(n_cells=400), pool, cal, 7, coupled=False, with_trace=False, keep_blocks=False)
    a = eh.analyse_history(null, null.d_true, K0, n_boot=40, seed=2)
    lo, hi = a["contrast_ci"]
    assert lo <= 1.0 <= hi, (lo, hi)
