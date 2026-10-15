"""The shared polymer -> nm / seconds calibration (simlive.calibration.polymer) and how a library config reads it."""
import json

import numpy as np
import pytest

from simlive.calibration import polymer as P


def fake_pool(n_sims=2, T=400, pairs=30, D_unit2=0.05, seed=0):
    """Random-walk pair vectors and MSDs from a free random walk of per-block variance 2 D_unit2 per axis (MSD(lag) = 6 D_unit2 lag)."""
    rng = np.random.default_rng(seed)
    rel = rng.normal(size=(n_sims, T, pairs, 3)).astype(np.float32) * 3.0
    lags = P.LAGS
    return P.PolymerPool(rel, 6 * D_unit2 * lags, 12 * D_unit2 * lags, ("a", "b"), 100.0)


def test_calibration_reproduces_the_requested_contact_fraction_and_msd():
    pool = fake_pool()
    cal = P.calibrate(pool, working_nm=400.0, contact_fraction=0.1, D_um2_s=0.002, lag_s=60.0)
    d_nm = pool.distances() * cal.nm_per_unit
    assert abs((d_nm < 400.0).mean() - 0.1) < 0.005
    # the single-locus MSD at lag_s is 6 D lag_s
    msd_at_lag = float(cal.msd_single_s(60.0))
    assert np.isclose(msd_at_lag, 6 * 0.002 * 60.0 * 1e6, rtol=0.02)


def test_given_numbers_override_the_derivation():
    cal = P.calibrate(fake_pool(), nm_per_unit=60.0, block_s=2.0)
    assert cal.nm_per_unit == 60.0 and cal.block_s == 2.0


def test_json_roundtrip_keeps_everything(tmp_path):
    cal = P.calibrate(fake_pool())
    f = cal.save(tmp_path / "x.json")
    back = P.Calibration.from_json(f.read_text(encoding="utf-8"))
    assert back == cal
    assert back.sims == ("a", "b") and back.working_nm == 400.0 and back.contact_fraction == 0.1


def test_scale_from_config_precedence(tmp_path):
    cal = P.calibrate(fake_pool())
    f = cal.save(tmp_path / "mine.json")
    nm, blk, src = P.scale_from_config({"calibration": str(f)})
    assert (nm, blk) == (cal.nm_per_unit, cal.block_s) and "calibration" in src
    nm, blk, src = P.scale_from_config({"calibration": str(f), "block_duration_s": 3.0})
    assert nm == cal.nm_per_unit and blk == 3.0 and "overridden" in src
    nm, blk, src = P.scale_from_config({"nm_per_unit": 50, "block_duration_s": 10})
    assert (nm, blk) == (50.0, 10.0)
    with pytest.raises(KeyError):
        P.scale_from_config({"nm_per_unit": 50})
    with pytest.raises(FileNotFoundError):
        P.scale_from_config({"calibration": str(tmp_path / "missing.json")})


def test_the_committed_default_calibration_is_loadable_and_self_consistent():
    cal = P.load_calibration("default")
    assert 20 < cal.nm_per_unit < 200 and 0.05 < cal.block_s < 30
    assert cal.applies_to == "locus" and cal.D_um2_s == 0.002 and cal.lag_s == 60.0
    assert np.isclose(float(cal.msd_single_s(cal.lag_s)), 6 * cal.D_um2_s * cal.lag_s * 1e6, rtol=0.05)
    json.loads(P.calibration_path("default").read_text(encoding="utf-8"))
    assert P.default_scale()[:2] == (cal.nm_per_unit, cal.block_s)
