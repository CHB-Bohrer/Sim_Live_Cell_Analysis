"""3D z-stacks: nuclei heights, 3D loci, the 3D PSF, axial errors, and 3D localization."""
import numpy as np
import pandas as pd

from simlive.io.runs import REPO_ROOT, load_config
from simlive.stage1_chromatin.toy_loci import simulate_toy_loci
from simlive.stage2_cells.motion import add_z, simulate_cells, z_enabled
from simlive.stage3_microscopy.render import Imager, axial_sigma_um, truth_in_image

CFG = REPO_ROOT / "configs" / "loci_3d.yaml"


def _cfg(*extra):
    return load_config(CFG, ["geometry.fov_um=[24, 24]", "geometry.cell_radius_um=4.0", "cells.n_cells=2", "acquisition.n_frames=4", "optics.bleach_tau_s=null",
                             "imaging_errors.enabled=false", *extra])


def _scene(cfg, z_ratio=1.0):
    cells = pd.DataFrame({"t": range(4), "cell_id": 1, "y_um": 12.0, "x_um": 12.0, "radius_um": 5.0, "aspect": 1.0, "angle_rad": 0.0,
                          "bound_radius_um": 5.5, "parent_id": 0, **{f"{a}{k}": 0.0 for k in (2, 3, 4, 5) for a in ("amp", "phase")},
                          "z_um": 7.0, "z_semi_um": 5.0 * z_ratio})
    truth = pd.DataFrame([{"t": t, "cell_id": 1, "locus_id": k, "y_um": 12.0 + k, "x_um": 12.0, "z_um": 7.0 + 0.5 * k} for t in range(4) for k in (0, 1)])
    return cells, truth


def test_add_z_is_a_separate_stream_and_keeps_nuclei_above_the_coverslip():
    cfg = _cfg()
    cells = simulate_cells(cfg, np.random.default_rng(1))
    z1 = add_z(cells, cfg, np.random.default_rng(2))
    z2 = add_z(cells, cfg, np.random.default_rng(2))
    pd.testing.assert_frame_equal(z1, z2)
    assert z_enabled(cfg) and (z1.z_um - z1.z_semi_um >= 0.19).all()
    assert np.allclose(z1.y_um, cells.y_um) and np.allclose(z1.x_um, cells.x_um)           # the lateral simulation is untouched
    # a cell's height changes slowly through time
    g = z1[z1.cell_id == 1].sort_values("t")
    assert g.z_um.diff().abs().max() < 1.0


def test_toy_loci_get_a_height_inside_the_nucleus():
    cfg = _cfg()
    cells = add_z(simulate_cells(cfg, np.random.default_rng(1)), cfg, np.random.default_rng(2))
    cfg["loci"]["n_loci"] = 2
    truth = simulate_toy_loci(cells, cfg, np.random.default_rng(3))
    m = truth.merge(cells[["t", "cell_id", "z_um", "z_semi_um"]], on=["t", "cell_id"], suffixes=("", "_cell"))
    assert ((m.z_um - m.z_um_cell).abs() <= 0.7 * m.z_semi_um + 1e-9).all()


def test_nucleus_is_a_3d_ellipsoid_with_the_right_volume_and_slice_pattern():
    cfg = _cfg()
    cells, truth = _scene(cfg)
    im = Imager(cfg, cells, truth, entropy=1)
    assert im.is3d and im.nz == 15 and abs(im.zk.mean() - 7.0) < 1e-9
    lab = im.frame(0, channels=["nucleus"])["labels"]
    vox = (lab == 1).sum() * im.px_um**2 * im.zstep
    # the stack (15 x 0.8 um = 11.2 um) holds the whole 10 um sphere
    assert abs(vox - 4 / 3 * np.pi * 5.0**3) < 0.1 * 4 / 3 * np.pi * 5.0**3
    areas = [(lab[k] == 1).sum() for k in range(im.nz)]
    assert np.argmax(areas) == im.nz // 2 and areas[0] < areas[im.nz // 2] / 3          # largest in the middle, small near the poles


def test_point_source_axial_profile_matches_the_axial_psf():
    cfg = _cfg("acquisition.z_stack.n_slices=25", "acquisition.z_stack.step_um=0.25")
    cells, truth = _scene(cfg)
    im = Imager(cfg, cells, truth, entropy=1)
    E = im.E.with_enabled(["psf_blur"])
    vol = im.signal("locus0", 0, E)[0].astype(float)
    k0 = int(np.argmax(vol.max(axis=(1, 2))))
    assert abs(im.zk[k0] - 7.0) <= 0.125 + 1e-9                                           # brightest plane = the one nearest the true z
    prof = vol.sum(axis=(1, 2))
    z = im.zk
    sz = np.sqrt((prof * (z - (prof * z).sum() / prof.sum()) ** 2).sum() / prof.sum())
    assert abs(sz - axial_sigma_um(cfg["optics"], 580)) < 0.1 * axial_sigma_um(cfg["optics"], 580)
    # lateral position is unaffected by the stack
    yy, xx = np.indices(vol.shape[1:])
    img = vol.sum(0)
    assert abs((img * yy).sum() / img.sum() - 12.0 / im.px_um) < 0.1


def test_confocal_has_a_narrower_axial_psf_than_widefield():
    wf = axial_sigma_um({"wavelength_nm": 580, "NA": 1.2, "modality": "widefield"})
    cf = axial_sigma_um({"wavelength_nm": 580, "NA": 1.2, "modality": "confocal"})
    assert cf < wf and 0.3 < wf < 0.6


def test_focus_error_and_chromatic_z_shift_displace_objects_along_z():
    cfg = _cfg("imaging_errors.enabled=true", "imaging_errors.chromatic_shift.shift_z_um={locus0: 0.8, locus1: 0.0}",
               "imaging_errors.defocus.rms_um=0.5", "acquisition.z_stack.n_slices=25", "acquisition.z_stack.step_um=0.25")
    cells, truth = _scene(cfg)
    im = Imager(cfg, cells, truth, entropy=7)
    t = int(np.argmax([abs(im.defocus_um(t)) for t in range(4)]))
    E = im.E.with_enabled(["psf_blur", "defocus", "chromatic_shift"])
    vol = im.signal("locus0", t, E)[0].astype(float)
    prof = vol.sum(axis=(1, 2))
    zc = (prof * im.zk).sum() / prof.sum()
    assert abs(zc - (7.0 + 0.8 + im.defocus_um(t, E))) < 0.05


def test_slice_skew_moves_a_locus_between_the_first_and_last_slice():
    cfg = _cfg("acquisition.z_stack.slice_time_s=2.0", "acquisition.z_stack.n_slices=9", "acquisition.frame_interval_s=20")
    cells, truth = _scene(cfg)
    truth.loc[truth.t >= 1, "x_um"] += 2.0                                                 # the locus jumps 2 um between frames 0 and 1
    im = Imager(cfg, cells, truth, entropy=1)
    base = im.E.with_enabled(["psf_blur"])
    skew = im.E.with_enabled(["psf_blur", "slice_skew"])
    P0 = im._points_at(im.loci[0].by_t, 0, base)
    P1 = im._points_at(im.loci[0].by_t, 0, skew)
    assert np.allclose(P0[1][0], P0[1][-1])
    assert (P1[1][-1] - P1[1][0])[0] > 0.0 and abs((P1[1][-1] - P1[1][0])[0] * im.px_um - 2.0 * (8 * 2.0 / 20)) < 1e-6   # linear in slice time


def test_frames_are_reproducible_in_3d_and_images_are_stacks():
    cfg = _cfg("imaging_errors.enabled=true")
    cells, truth = _scene(cfg)
    a = Imager(cfg, cells, truth, entropy=3).frame(2)
    b = Imager(cfg, cells, truth, entropy=3).frame(2)
    for c in a["images"]:
        assert a["images"][c].shape == (15, 96 * 0 + a["images"][c].shape[1], a["images"][c].shape[2]) and np.array_equal(a["images"][c], b["images"][c])
    assert a["labels"].ndim == 3
    # every slice has independent noise
    n = a["images"]["locus0"]
    assert not np.array_equal(n[0], n[1])
    mip = Imager.project(n, "max")
    assert mip.shape == n.shape[1:] and (mip >= n[0]).all()


def test_truth_in_image_includes_axial_drift_and_the_drawn_z():
    cfg = _cfg("imaging_errors.enabled=true", "imaging_errors.stage_drift.drift_z_um_per_frame=0.5")
    cells, truth = _scene(cfg)
    im = Imager(cfg, cells, truth, entropy=5)
    it = im.imaging_truth()
    ti = truth_in_image(truth, it, im.px_um, im.static_truth()["channel_shift_px"], im.static_truth()["z_stack"]["channel_shift_z_um"])
    r = ti[(ti.t == 3) & (ti.locus_id == 0)].iloc[0]
    assert abs(r.z_um - (7.0 + it.drift_z_um[3])) < 1e-9
    assert abs(r.z_drawn_um - (r.z_um + 0.15 + it.defocus_um[3])) < 1e-9


def test_3d_localizer_finds_a_noise_free_spot_to_sub_voxel_accuracy(tmp_path):
    import tifffile
    from simlive.stage5_linking.celldata import CellData
    from simlive.stage6_analysis.localize3d import localize_loci_3d

    cfg = _cfg("acquisition.z_stack.n_slices=15")
    cells, truth = _scene(cfg)
    im = Imager(cfg, cells, truth, entropy=1)
    E = im.E.with_enabled(["psf_blur", "shot_noise"])
    stack = np.stack([im.frame(t, E, ["locus0"])["images"]["locus0"] for t in range(2)])
    tifffile.imwrite(tmp_path / "locus0_zstack.tif", np.rint(stack).astype(np.uint16), photometric="minisblack")
    S = 61
    fr = pd.DataFrame({"t": [0, 1], "crop_y0": [int(12 / im.px_um) - 30] * 2, "crop_x0": [int(12 / im.px_um) - 30] * 2})
    mask = np.ones((2, S, S), bool)
    cell = CellData(1, fr, mask, {"locus0": np.zeros((2, S, S), np.uint16)}, S)
    pos = localize_loci_3d(cell, {"locus0": str(tmp_path / "locus0_zstack.tif")}, im.zk.tolist(), n_loci=1, min_snr=3.0)
    t0 = truth[(truth.t == 0) & (truth.locus_id == 0)].iloc[0]
    p = pos[pos.t == 0].iloc[0]
    assert abs(p.lab_y_px - t0.y_um / im.px_um) < 0.5 and abs(p.lab_x_px - t0.x_um / im.px_um) < 0.5 and abs(p.z_um - t0.z_um) < 0.3
