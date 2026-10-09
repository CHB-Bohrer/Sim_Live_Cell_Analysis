"""Image-error model (stage 3): every source is switchable, deterministic per frame, and does what it says."""
import numpy as np
import pandas as pd

from simlive.io.runs import REPO_ROOT, load_config
from simlive.stage3_microscopy import errors as ERR
from simlive.stage3_microscopy.render import Imager, cells_in_image, truth_in_image

CFG = REPO_ROOT / "configs" / "loci_demo.yaml"


def _cfg(*extra):
    return load_config(CFG, ["geometry.fov_um=[24, 24]", "acquisition.n_frames=6", "cells.n_cells=1", "optics.bleach_tau_s=null",
                             *extra])


def _scene(n_t=6):
    cells = pd.DataFrame({"t": range(n_t), "cell_id": 1, "y_um": 12.0, "x_um": 12.0, "radius_um": 5.0, "aspect": 1.0,
                          "angle_rad": 0.0, "bound_radius_um": 5.5, "parent_id": 0,
                          **{f"{a}{k}": 0.0 for k in (2, 3, 4, 5) for a in ("amp", "phase")}})
    truth = pd.DataFrame([{"t": t, "cell_id": 1, "locus_id": k, "y_um": 12.0 + 1.5 * (k - 0.5), "x_um": 12.0 + 0.5 * k}
                          for t in range(n_t) for k in (0, 1)])
    return cells, truth


def test_every_source_has_documentation_and_a_place_in_the_order():
    assert set(ERR.ORDER) == set(ERR.SOURCES) and len(ERR.ORDER) >= 20
    for s in ERR.SOURCES.values():
        assert s["group"] in ERR.GROUPS and s["title"] and s["what"] and s["effect"]


def test_master_switch_off_keeps_only_the_original_model():
    E = ERR.ErrorSet(_cfg("imaging_errors.enabled=false"))
    assert E.enabled_ids() == ERR.BASELINE and not E.master
    full = ERR.ErrorSet(_cfg())
    assert full.master and len(full.enabled_ids()) == len(ERR.ORDER)


def test_frames_are_reproducible_and_independent_of_other_frames():
    cfg = _cfg()
    cells, truth = _scene()
    a = Imager(cfg, cells, truth, entropy=7)
    b = Imager(cfg, cells, truth, entropy=7)
    f3 = a.frame(3)["images"]
    a.frame(1)                                                              # rendering another frame first changes nothing
    for c in f3:
        assert np.array_equal(f3[c], b.frame(3)["images"][c])
    assert not np.array_equal(f3["nucleus"], Imager(cfg, cells, truth, entropy=8).frame(3)["images"]["nucleus"])


def test_switching_one_source_off_does_not_change_the_noise_of_another():
    cfg = _cfg()
    cells, truth = _scene()
    im = Imager(cfg, cells, truth, entropy=3)
    only_shot = im.E.with_enabled(["shot_noise", "qe", "psf_blur"])
    plus_read = im.E.with_enabled(["shot_noise", "qe", "psf_blur", "read_noise"])
    a = im.frame(2, only_shot, ["locus0"])["images"]["locus0"]
    b = im.frame(2, plus_read, ["locus0"])["images"]["locus0"]
    assert not np.array_equal(a, b)
    assert abs(float((b - a).mean())) < 0.2                                # only the (zero-mean) read noise differs


def test_stage_drift_moves_nuclei_labels_and_spots_together_and_truth_follows():
    cfg = _cfg("imaging_errors.stage_drift.drift_px_per_frame=1.5")
    cells, truth = _scene()
    im = Imager(cfg, cells, truth, entropy=5)
    px, t = im.px_um, 5
    dy, dx = im.drift_px(t)
    assert abs(dy) + abs(dx) > 0.5
    lab = im.frame(t, channels=["nucleus"])["labels"]
    yy, xx = np.nonzero(lab)
    assert abs(yy.mean() - (12.0 / px + dy)) < 0.6 and abs(xx.mean() - (12.0 / px + dx)) < 0.6
    ti = truth_in_image(truth, im.imaging_truth(), px, im.static_truth()["channel_shift_px"])
    r = ti[(ti.t == t) & (ti.locus_id == 0)].iloc[0]
    assert abs(r.y_um - (truth[(truth.t == t) & (truth.locus_id == 0)].y_um.iloc[0] + dy * px)) < 1e-9
    ci = cells_in_image(cells, im.imaging_truth(), px)
    assert abs(ci[ci.t == t].y_um.iloc[0] - (12.0 + dy * px)) < 1e-9


def test_chromatic_shift_offsets_a_locus_channel_by_the_configured_amount():
    cfg = _cfg("imaging_errors.chromatic_shift.shift_px={locus0: [2.0, -1.0], locus1: [0.0, 0.0]}")
    cells, truth = _scene()
    im = Imager(cfg, cells, truth, entropy=1)
    clean = im.E.with_enabled(["psf_blur", "chromatic_shift"])               # no noise: the centroid is exact
    sig = im.signal("locus0", 0, clean)[0].astype(float)
    yy, xx = np.indices(sig.shape)
    cy, cx = (sig * yy).sum() / sig.sum(), (sig * xx).sum() / sig.sum()
    t0 = truth[(truth.t == 0) & (truth.locus_id == 0)].iloc[0]
    assert abs(cy - (t0.y_um / im.px_um + 2.0)) < 0.05 and abs(cx - (t0.x_um / im.px_um - 1.0)) < 0.05


def test_defocus_widens_and_dims_a_spot_but_conserves_its_photons():
    cfg = _cfg("imaging_errors.defocus.rms_um=0.4")
    cells, truth = _scene()
    im = Imager(cfg, cells, truth, entropy=2)
    ts = np.argsort([abs(im.defocus_um(t)) for t in range(6)])
    t_good, t_bad = int(ts[0]), int(ts[-1])
    assert abs(im.defocus_um(t_bad)) > abs(im.defocus_um(t_good))
    E = im.E.with_enabled(["psf_blur", "defocus"])
    a, b = (im.signal("locus0", t, E)[0].astype(float) for t in (t_good, t_bad))
    assert b.max() < a.max() and abs(a.sum() - b.sum()) < 0.02 * a.sum()


def test_haze_and_crosstalk_conserve_light_and_move_it_where_they_should():
    cfg = _cfg("imaging_errors.crosstalk.pairs=[[nucleus, locus1, 0.1]]")
    cells, truth = _scene()
    im = Imager(cfg, cells, truth, entropy=4)
    E0 = im.E.with_enabled(["psf_blur"])
    Eh = im.E.with_enabled(["psf_blur", "haze"])
    s0, sh = im.signal("nucleus", 0, E0)[0].astype(float), im.signal("nucleus", 0, Eh)[0].astype(float)
    assert sh.max() < s0.max() and 0.9 * s0.sum() < sh.sum() <= s0.sum() * 1.01        # broader; light conserved (to interpolation accuracy) except what leaves the field
    ghost = (im.frame(0, im.E.with_enabled(["psf_blur", "crosstalk"]), ["locus1"])["images"]["locus1"]
             - im.frame(0, E0, ["locus1"])["images"]["locus1"])
    assert ghost.sum() > 0.05 * s0.sum()                                                # the nucleus now shows up in the locus1 channel


def test_camera_noise_matches_theory():
    cfg = _cfg()
    E = ERR.ErrorSet(cfg).with_enabled(["qe", "shot_noise", "read_noise"])
    E.params["read_noise"]["pixel_cv"] = 0.0
    adu = ERR.Camera((200, 200), E, cfg["optics"], 0.1, 1, 0)(np.full((200, 200), 400.0), 0)
    o = cfg["optics"]
    qe, gain, off, rn = o["quantum_efficiency"], o["gain_adu_per_e"], o["offset_adu"], o["read_noise_e"]
    assert abs(adu.mean() - (400 * qe * gain + off)) < 0.5
    assert abs(adu.var() - gain**2 * (400 * qe + rn**2)) < 0.1 * gain**2 * (400 * qe + rn**2)


def test_adc_rounds_to_integers_and_saturates():
    cfg = _cfg()
    E = ERR.ErrorSet(cfg).with_enabled(["adc"])
    cam = ERR.Camera((20, 20), E, cfg["optics"], 0.1, 1, 0)
    a = cam(np.full((20, 20), 1e9), 0)
    assert np.all(a == a.max()) and a.max() == cfg["optics"]["offset_adu"] + 30000 * cfg["optics"]["gain_adu_per_e"]
    b = cam(np.full((20, 20), 3.3), 0)
    assert np.all(b == np.rint(b))


def test_the_imaging_truth_table_and_saved_json_describe_what_happened():
    cfg = _cfg()
    cells, truth = _scene()
    im = Imager(cfg, cells, truth, entropy=9)
    it = im.imaging_truth()
    assert len(it) == 6 and {"drift_y_px", "defocus_um", "flicker_nucleus", "sigma_px_locus1"} <= set(it.columns)
    st = im.static_truth()
    assert st["channel_shift_px"]["locus0"] == (0.25, -0.15) and st["sources"]["stage_drift"]["on"]
