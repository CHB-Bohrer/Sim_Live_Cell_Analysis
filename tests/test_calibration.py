"""Calibration from simulated dark / flat-field / bead acquisitions: blind estimates must match the simulation's true instrument."""
import numpy as np
import pytest

from simlive.io.runs import REPO_ROOT, load_config
from simlive.stage3_microscopy import errors as ERR
from simlive.stage3_microscopy.calibration import fit as F
from simlive.stage3_microscopy.calibration import sim as S
from simlive.stage3_microscopy.calibration.model import Calibration
from simlive.stage3_microscopy.calibration.run import run_calibration, settings_for
from simlive.stage3_microscopy.render import entropy_for


def _cfg(*extra):
    return load_config(REPO_ROOT / "configs" / "loci_demo.yaml", ["geometry.fov_um=[30, 30]", *extra])


FAST = {"dark": {"n_frames": 40, "exposures_s": [0.01, 1.0, 5.0], "extra_frames": 8},
        "flat": {"levels_photons": [0, 300, 1000, 3000, 8000, 20000], "n_frames": 10},
        "flicker": {"level_photons": 3000, "n_frames": 80},
        "beads": {"n_beads": 30, "photons": 30000, "n_frames": 6, "n_frames_3d": 4, "z_um": 1.0}}


@pytest.fixture(scope="module")
def cal2d():
    cfg = _cfg()
    st = {**settings_for(cfg), **FAST}
    cal, table = run_calibration(cfg, entropy_for(3), st, progress=lambda *_: None)
    return cfg, cal, table


def _row(table, q, ch=None):
    t = table[table.quantity == q]
    return t[t.channel == ch].iloc[0] if ch else t.iloc[0]


def test_camera_properties_are_recovered(cal2d):
    cfg, cal, table = cal2d
    for ch in ("nucleus", "locus0"):
        assert abs(_row(table, "gain", ch).error_pct) < 3
        assert abs(_row(table, "offset", ch).error) < 0.2
        assert abs(_row(table, "read noise", ch).error_pct) < 5
        assert abs(_row(table, "row noise", ch).error_pct) < 15
        assert abs(_row(table, "offset pattern (DSNU)", ch).error_pct) < 10
        assert abs(_row(table, "dark current", ch).error_pct) < 15
        assert abs(_row(table, "pixel gain spread (PRNU)", ch).error_pct) < 10


def test_hot_pixels_are_found_without_false_alarms(cal2d):
    cfg, cal, table = cal2d
    r = _row(table, "hot pixels found", "nucleus")
    assert r.estimate == r.truth and "0 false" in r.note


def test_illumination_pattern_and_flicker_are_recovered(cal2d):
    cfg, cal, table = cal2d
    ch = "locus1"
    assert abs(_row(table, "vignetting (loss at the corners)", ch).error) < 0.01
    assert abs(_row(table, "illumination tilt", ch).error) < 0.01
    assert abs(_row(table, "flicker", ch).error_pct) < 25            # 80 frames of an AR(1) series: a rough estimate


def test_chromatic_shift_and_psf_width_from_beads(cal2d):
    cfg, cal, table = cal2d
    for ch in ("locus0", "locus1"):
        assert abs(_row(table, "chromatic shift y", ch).error) < 0.03
        assert abs(_row(table, "chromatic shift x", ch).error) < 0.03
    for ch in ("nucleus", "locus0", "locus1"):
        assert abs(_row(table, "PSF width (lateral sigma)", ch).error_pct) < 3
    assert _row(table, "stage drift (rms error of the estimated path)").estimate < 0.02


def test_correction_removes_fixed_patterns_and_hot_pixels(cal2d):
    cfg, cal, table = cal2d
    ch = "nucleus"
    E = ERR.ErrorSet(cfg)
    flat = S.acquire_flat(cfg, E, entropy_for(3), ch, 5000.0, 6, tag=77).mean(axis=0)          # a new uniform slide, other session
    corr = cal.correct(flat, ch)
    inner = (slice(20, -20), slice(20, -20))
    assert corr[inner].std() / corr[inner].mean() < 0.5 * flat[inner].std() / flat[inner].mean()
    hot = cal.maps[ch]["hot"]
    assert corr[hot].mean() < flat[hot].mean() - 10                                                # hot pixels were replaced


def test_registration_correction_moves_positions_back(cal2d):
    cfg, cal, table = cal2d
    y, x = cal.correct_positions(np.array([100.25]), np.array([50.0 - 0.15]), "locus0")
    assert abs(y[0] - 100.0) < 0.03 and abs(x[0] - 50.0) < 0.03


def test_calibration_round_trips_through_disk(cal2d, tmp_path):
    cfg, cal, table = cal2d
    cal.save(tmp_path)
    c2 = Calibration.load(tmp_path)
    assert c2.channels == cal.channels and np.array_equal(c2.maps["locus0"]["hot"], cal.maps["locus0"]["hot"])
    assert np.allclose(c2.maps["nucleus"]["offset"], cal.maps["nucleus"]["offset"])
    assert c2.shift_px("locus1") == cal.shift_px("locus1")


def test_calibration_data_has_independent_noise_but_the_same_instrument():
    cfg = _cfg()
    E = ERR.ErrorSet(cfg)
    a = S.acquire_dark(cfg, E, 5, "nucleus", 0.1, 8, tag=0)
    b = S.acquire_dark(cfg, E, 5, "nucleus", 0.1, 8, tag=1)
    assert not np.array_equal(a, b)
    tame = lambda m: np.clip(m, *np.percentile(m, [1, 99.5])).ravel()       # cosmic-ray blobs would dominate a plain correlation
    ma, mb = tame(a.mean(0)), tame(b.mean(0))
    assert np.corrcoef(ma, mb)[0, 1] > 0.5                                  # the fixed pattern (offsets, hot pixels) is shared
    c = S.acquire_dark(cfg, ERR.ErrorSet(cfg), 6, "nucleus", 0.1, 8, tag=0)  # another instrument
    assert np.corrcoef(ma, tame(c.mean(0)))[0, 1] < 0.2


def test_3d_beads_give_axial_psf_and_axial_chromatic_shift():
    cfg = _cfg("geometry.fov_um=[22, 22]", "acquisition.z_stack.enabled=true", "acquisition.z_stack.n_slices=11", "acquisition.z_stack.step_um=0.5",
               "geometry.z.z_center_um=7.0", "imaging_errors.stage_drift.drift_z_um_per_frame=0.02", "imaging_errors.chromatic_shift.shift_z_um={locus0: 0.3, locus1: -0.2}")
    st = {**settings_for(cfg), **FAST}
    st["dark"] = {"n_frames": 12, "exposures_s": [0.01, 1.0], "extra_frames": 4}
    st["flat"] = {"levels_photons": [0, 1000, 5000], "n_frames": 6}
    st["flicker"] = {"level_photons": 3000, "n_frames": 20}
    st["beads"] = {"n_beads": 12, "photons": 40000, "n_frames": 3, "n_frames_3d": 3, "z_um": 1.0}
    cal, table = run_calibration(cfg, entropy_for(2), st, progress=lambda *_: None)
    for ch, z in (("locus0", 0.3), ("locus1", -0.2)):
        r = _row(table, "chromatic shift z", ch)
        assert abs(r.estimate - z) < 0.06
    r = _row(table, "PSF width (axial sigma)", "nucleus")
    assert abs(r.error_pct) < 12
    assert "drift_z_um_series" in cal.beads
