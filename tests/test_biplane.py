"""Biplane acquisition (Bruker Vutara style): two planes 600 nm apart per stage position, half the light each, registration error of the second half."""
import numpy as np
import pandas as pd

from simlive.io.runs import REPO_ROOT, load_config
from simlive.stage3_microscopy.render import Imager

CFG = REPO_ROOT / "configs" / "loci_3d.yaml"


def _cfg(biplane: bool, *extra):
    base = ["geometry.fov_um=[20, 20]", "geometry.cell_radius_um=4.0", "cells.n_cells=1", "acquisition.n_frames=4", "optics.bleach_tau_s=null",
            "acquisition.z_stack.center_um=7.0", "fiducials.enabled=false", "imaging_errors.enabled=true"]
    if biplane:
        base += ["acquisition.z_stack.n_slices=6", "acquisition.z_stack.biplane={enabled: true, separation_um: 0.6, photon_fraction: 0.5}"]
    else:
        base += ["acquisition.z_stack.n_slices=12", "acquisition.z_stack.step_um=0.6"]
    return load_config(CFG, base + list(extra))


def _scene(n_t=4):
    cells = pd.DataFrame({"t": range(n_t), "cell_id": 1, "y_um": 10.0, "x_um": 10.0, "radius_um": 4.0, "aspect": 1.0, "angle_rad": 0.0,
                          "bound_radius_um": 4.5, "parent_id": 0, **{f"{a}{k}": 0.0 for k in (2, 3, 4, 5) for a in ("amp", "phase")},
                          "z_um": 7.0, "z_semi_um": 4.0})
    truth = pd.DataFrame([{"t": t, "cell_id": 1, "locus_id": k, "y_um": 10.0, "x_um": 10.0 + k, "z_um": 7.0} for t in range(n_t) for k in (0, 1)])
    return cells, truth


def test_planes_come_in_pairs_that_tile_the_sample_evenly():
    cells, truth = _scene()
    im = Imager(_cfg(True), cells, truth, entropy=1)
    assert im.biplane and im.n_stage == 6 and im.nz == 12
    assert np.allclose(np.diff(im.zk), 0.6)                            # pairs 0.6 um apart, pairs 1.2 um apart -> evenly 0.6 um overall
    assert np.allclose(im.zk[1::2] - im.zk[0::2], 0.6)
    ti = im.time_index(0)
    assert list(ti[0::2]) == list(ti[1::2]) == list(range(6))          # both planes of a pair are recorded at the same moment
    assert list(im.time_index(1)[0::2]) == list(range(5, -1, -1))      # down-scan
    plain = Imager(_cfg(False), cells, truth, entropy=1)
    assert np.allclose(plain.zk, im.zk) and not plain.biplane


def test_the_stage_error_is_shared_by_the_two_planes_of_a_pair():
    cells, truth = _scene()
    im = Imager(_cfg(True), cells, truth, entropy=1)
    e = im.piezo_error_um(0)
    assert len(e) == 12 and np.array_equal(e[0::2], e[1::2])


def test_each_plane_gets_half_of_the_light():
    cells, truth = _scene()
    E = ["psf_blur"]
    bi = Imager(_cfg(True), cells, truth, entropy=1)
    pl = Imager(_cfg(False), cells, truth, entropy=1)
    a = bi.signal("locus0", 0, bi.E.with_enabled(E))[0].sum(axis=(1, 2))
    b = pl.signal("locus0", 0, pl.E.with_enabled(E))[0].sum(axis=(1, 2))
    assert np.allclose(a, 0.5 * b, rtol=1e-4, atol=1e-3) and b.max() > 100


def test_registration_error_shifts_the_second_plane_and_its_height():
    cells, truth = _scene()
    cfg = _cfg(True, "imaging_errors.biplane_registration.enabled=true", "imaging_errors.biplane_registration.shift_px=[0.5, 0.0]",
               "imaging_errors.biplane_registration.separation_error_um=0.05", "imaging_errors.piezo.enabled=false")
    im = Imager(cfg, cells, truth, entropy=1)
    E = im.E.with_enabled(["psf_blur", "biplane_registration"])
    vol = im.signal("locus0", 0, E)[0].astype(float)
    k = int(np.argmin(np.abs(im.zk - 7.0)))                           # the two planes around the locus: k (lower of a pair) and k + 1 or k - 1
    ka, kb = (k, k + 1) if k % 2 == 0 else (k - 1, k)
    yy = np.arange(vol.shape[1])[:, None]

    def cy(p):
        return float((vol[p] * yy).sum() / vol[p].sum())

    assert abs((cy(kb) - cy(ka)) - 0.5) < 0.02                         # the second plane's image is displaced by 0.5 px in y
    assert np.allclose(im.plane_z(0, E)[1::2] - im.zk[1::2], 0.05) and np.allclose(im.plane_z(0, E)[0::2], im.zk[0::2])
    off = Imager(_cfg(True, "imaging_errors.biplane_registration.enabled=false"), cells, truth, entropy=1)
    assert not off.plane_offset_px().any()
    assert im.static_truth()["z_stack"]["biplane"] is True


def test_beads_in_the_image_measure_the_registration_offset(tmp_path):
    import tifffile

    from simlive.stage3_microscopy.calibration.sim import _empty_cells
    from simlive.stage3_microscopy.fiducials import make_fiducials
    from simlive.stage6_analysis.biplane import estimate_plane_offset

    def render(shift):
        cfg = _cfg(True, "geometry.fov_um=[40, 40]", "acquisition.z_stack.center_um=3.3", "fiducials.enabled=true", "fiducials.n_beads=5",
                   "imaging_errors.biplane_registration.enabled=true", f"imaging_errors.biplane_registration.shift_px={shift}",
                   "imaging_errors.biplane_registration.separation_error_um=0.0", "imaging_errors.piezo.enabled=false",
                   "imaging_errors.cosmic_rays.enabled=false", "imaging_errors.hot_pixels.enabled=false", "imaging_errors.defocus.enabled=false")
        cells, _ = _scene()
        fid = make_fiducials(cfg, cells, np.random.default_rng(4))
        im = Imager(cfg, _empty_cells(cfg, True), None, None, None, 7, fiducials=fid, channels=["nucleus"])
        stack = np.stack([np.rint(im.frame(t, channels=["nucleus"])["images"]["nucleus"]).astype(np.uint16) for t in range(3)])
        path = tmp_path / f"s{abs(hash(str(shift)))}.tif"
        tifffile.imwrite(path, stack, photometric="minisblack")
        return path

    (dy, dx), tab = estimate_plane_offset(render([0.4, -0.3]))
    assert abs(dy - 0.4) < 0.06 and abs(dx + 0.3) < 0.06 and tab.bead.nunique() >= 3
    (dy0, dx0), _ = estimate_plane_offset(render([0.0, 0.0]))
    assert abs(dy0) < 0.06 and abs(dx0) < 0.06
