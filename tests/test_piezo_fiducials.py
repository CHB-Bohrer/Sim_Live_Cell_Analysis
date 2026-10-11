"""Piezo plane-position errors, fiducial beads on the slide, and the corrections that use them."""
import numpy as np
import pandas as pd
import tifffile

from simlive.io.runs import REPO_ROOT, load_config
from simlive.stage3_microscopy.calibration.sim import _empty_cells
from simlive.stage3_microscopy.fiducials import make_fiducials
from simlive.stage3_microscopy.render import Imager
from simlive.stage6_analysis.fiducials import correct_with_beads, correct_with_plane_positions, score_variants, track_fiducials

CFG = REPO_ROOT / "configs" / "loci_3d.yaml"


def _cfg(*extra):
    return load_config(CFG, ["geometry.fov_um=[24, 24]", "geometry.cell_radius_um=4.0", "cells.n_cells=1", "acquisition.n_frames=8",
                             "optics.bleach_tau_s=null", "acquisition.z_stack.n_slices=15", "acquisition.z_stack.center_um=7.0",
                             "fiducials.enabled=false", "imaging_errors.enabled=true", *extra])


def _scene(n_t=8):
    cells = pd.DataFrame({"t": range(n_t), "cell_id": 1, "y_um": 12.0, "x_um": 12.0, "radius_um": 4.0, "aspect": 1.0, "angle_rad": 0.0,
                          "bound_radius_um": 4.5, "parent_id": 0, **{f"{a}{k}": 0.0 for k in (2, 3, 4, 5) for a in ("amp", "phase")},
                          "z_um": 7.0, "z_semi_um": 4.0})
    truth = pd.DataFrame([{"t": t, "cell_id": 1, "locus_id": k, "y_um": 12.0, "x_um": 12.0 + k, "z_um": 7.0} for t in range(n_t) for k in (0, 1)])
    return cells, truth


def test_piezo_errors_have_the_modelled_structure():
    cfg = _cfg("imaging_errors.piezo.repeatability_nm=40", "imaging_errors.piezo.bow_um=0.2", "imaging_errors.piezo.hysteresis_um=0.1",
               "imaging_errors.piezo.settling_um=0.0")
    cells, truth = _scene(40)
    im = Imager(cfg, cells, truth, entropy=1, n_frames=40)
    err = np.array([im.piezo_error_um(t) for t in range(40)])
    up, down = err[::2], err[1::2]
    assert [im.scan_dir(t) for t in range(4)] == [1, -1, 1, -1]
    assert abs(up.mean() - down.mean() - 0.1) < 0.02                          # hysteresis: up-scan planes are higher by the set amount
    mid = err.mean(axis=0)
    assert mid[7] > mid[0] + 0.15 and abs(mid[0] - 0.0) < 0.05               # bow: zero at the ends, 0.2 um in the middle (+ the hysteresis mean 0)
    assert abs(np.concatenate([up - up.mean(0), down - down.mean(0)]).std() - 0.04) < 0.008     # random step error ~ 40 nm (per scan direction)
    z = im.plane_z(3)
    assert np.allclose(z - im.zk, im.piezo_error_um(3))
    assert list(im.time_index(1)) == list(range(14, -1, -1))                  # down-scan: the top plane is recorded first
    assert abs((im.piezo_readout(3) - z).std() - 0.005) < 0.003               # the position log: true heights + 5 nm of sensor noise
    off = Imager(_cfg("imaging_errors.piezo.enabled=false"), cells, truth, entropy=1)
    assert np.array_equal(off.plane_z(3), off.zk)


def test_planes_are_where_the_piezo_put_them_not_where_it_was_asked():
    cfg = _cfg("imaging_errors.piezo.repeatability_nm=150", "imaging_errors.piezo.bow_um=0", "imaging_errors.piezo.hysteresis_um=0",
               "imaging_errors.piezo.settling_um=0", "acquisition.z_stack.n_slices=25", "acquisition.z_stack.step_um=0.25")
    cells, truth = _scene()
    im = Imager(cfg, cells, truth, entropy=2)
    E = im.E.with_enabled(["psf_blur", "piezo"])
    vol = im.signal("locus0", 0, E)[0].astype(float)
    prof = vol.sum(axis=(1, 2))
    prof = prof / prof.max()
    sz = im.sigma_z_um("locus0", E)
    near = prof > 0.05
    expected_actual = np.exp(-0.5 * ((im.plane_z(0, E) - 7.0) / sz) ** 2)      # brightness of each plane follows the TRUE plane height ...
    expected_nominal = np.exp(-0.5 * ((im.zk - 7.0) / sz) ** 2)                 # ... not the nominal one
    assert np.abs(prof[near] - expected_actual[near] / expected_actual.max() * prof.max()).max() < 0.02 or         np.abs(prof[near] - expected_actual[near]).max() < 0.03
    assert np.abs(prof[near] - expected_nominal[near]).max() > 0.1


def test_fiducials_avoid_nuclei_and_are_only_seen_if_the_stack_reaches_the_slide():
    cfg = _cfg("fiducials.enabled=true", "fiducials.n_beads=3")
    cells, truth = _scene()
    fid = make_fiducials(cfg, cells, np.random.default_rng(3))
    assert len(fid) == 3 and (np.hypot(fid.y_um - 12, fid.x_um - 12) > 4.5 + 3.0 - 1e-9).all()
    assert (fid.z_um < 0.5).all()
    quiet = ["psf_blur"]
    high = Imager(cfg, cells, truth, entropy=1, fiducials=fid)                               # stack 3.4 to 10.6 um: the slide is out of range
    cfg_low = _cfg("fiducials.enabled=true", "acquisition.z_stack.center_um=5.6", "acquisition.z_stack.n_slices=15")   # 0.0 to 11.2 um
    low = Imager(cfg_low, cells, truth, entropy=1, fiducials=fid)
    y, x = int(fid.y_um[0] / high.px_um), int(fid.x_um[0] / high.px_um)
    a = high.signal("nucleus", 0, high.E.with_enabled(quiet))[0][:, y - 2:y + 3, x - 2:x + 3].max()
    b = low.signal("nucleus", 0, low.E.with_enabled(quiet))[0][:, y - 2:y + 3, x - 2:x + 3].max()
    assert b > 20 * max(a, 1e-3) and b > 1000          # out of range the bead only leaves its faint out-of-focus haze


def test_beads_measure_a_whole_picture_shift(tmp_path):
    cfg = _cfg("fiducials.enabled=true", "fiducials.n_beads=5", "geometry.fov_um=[40, 40]", "imaging_errors.piezo.enabled=false", "acquisition.z_stack.center_um=5.6",
               "imaging_errors.stage_drift.drift_z_um_per_frame=0.3", "imaging_errors.stage_drift.jitter_z_um=0.0",
               "imaging_errors.defocus.enabled=false", "imaging_errors.cosmic_rays.enabled=false", "imaging_errors.hot_pixels.enabled=false")
    cells, truth = _scene()
    empty = _empty_cells(cfg, True)
    fid = make_fiducials(cfg, cells, np.random.default_rng(4))
    im = Imager(cfg, empty, None, None, None, 7, fiducials=fid, channels=["nucleus"])
    stack = np.stack([np.rint(im.frame(t, channels=["nucleus"])["images"]["nucleus"]).astype(np.uint16) for t in range(8)])
    tifffile.imwrite(tmp_path / "nucleus_zstack.tif", stack, photometric="minisblack")
    shifts, fits = track_fiducials(tmp_path / "nucleus_zstack.tif", im.zk.tolist(), None)
    true = np.array([im.drift_z_um(t) for t in range(8)])
    true = true - true[0]
    assert shifts.n_beads.iloc[0] >= 4
    assert np.abs(shifts.dz_um.to_numpy() - true).max() < 0.06 and np.abs(true).max() > 0.3
    ty = np.array([im.drift_px(t)[0] for t in range(8)])
    assert np.abs(shifts.dy_px.to_numpy() - (ty - ty[0])).max() < 0.1


def test_corrections_and_scoring():
    pos = pd.DataFrame({"cell_id": 1, "t": [0, 0], "locus_id": [0, 1], "lab_y_px": [10.0, 12.0], "lab_x_px": [10.0, 10.0], "z_um": [5.0, 6.0]})
    shifts = pd.DataFrame({"t": [0], "dy_px": [1.0], "dx_px": [0.0], "dz_um": [0.5]})
    c = correct_with_beads(pos, shifts)
    assert np.allclose(c.lab_y_px, [9.0, 11.0]) and np.allclose(c.z_um, [4.5, 5.5])
    log = pd.DataFrame({"t": 0, "k": [0, 1, 2], "readout_um": [0.0, 1.1, 2.0]})
    p2 = pd.DataFrame({"cell_id": 1, "t": [0], "locus_id": [0], "lab_y_px": [0.0], "lab_x_px": [0.0], "z_um": [0.5]})
    assert abs(correct_with_plane_positions(p2, log, [0.0, 1.0, 2.0]).z_um.iloc[0] - 0.55) < 1e-9        # halfway between logged 0.0 and 1.1
    truth = pd.DataFrame({"cell_id": 1, "t": [0, 0], "locus_id": [0, 1], "y_um": [1.0, 1.26], "x_um": [1.0, 1.0], "z_um": [5.0, 6.0]})
    same = pos.assign(lab_y_px=truth.y_um / 0.13, lab_x_px=truth.x_um / 0.13)
    t = score_variants({"perfect": same}, truth, 0.13)
    assert t.iloc[0]["distance error, typical (nm)"] < 1e-6
