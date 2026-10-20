"""CPU tests for stage 2 (motion), stage 3 (render) and the tracking scorer."""
import numpy as np

from simlive.io.runs import load_config, REPO_ROOT
from simlive.stage2_cells.motion import simulate_cells
from simlive.stage3_microscopy.render import render_nuclei
from simlive.stage7_validation.tracking_metrics import ctc_ground_truth, evaluate, gt_table, write_ctc

CFG = REPO_ROOT / "configs" / "tracking_demo.yaml"


def small_cfg(**over):
    return load_config(CFG, ["acquisition.n_frames=12", "cells.n_cells=8", "imaging_errors.enabled=false", *[f"{k}={v}" for k, v in over.items()]])


def test_motion_is_deterministic_and_collision_free():
    cfg = small_cfg()
    a = simulate_cells(cfg, np.random.default_rng(5))
    b = simulate_cells(cfg, np.random.default_rng(5))
    assert a.equals(b)
    for _, g in a.groupby("t"):
        p, r = g[["y_um", "x_um"]].to_numpy(), g["bound_radius_um"].to_numpy()
        d = np.hypot(*(p[:, None] - p[None]).transpose(2, 0, 1))
        np.fill_diagonal(d, np.inf)
        assert (d >= r[:, None] + r[None] - 1e-6).all(), "cells overlap"


def test_division_creates_daughters_with_parent():
    cells = simulate_cells(small_cfg(**{"cells.p_divide_per_frame": 0.3, "cells.max_cells": 20}),
                           np.random.default_rng(1))
    kids = cells[cells.parent_id > 0]
    assert len(kids) > 0
    t = gt_table(cells).set_index("label")
    for cid, row in t[t.parent > 0].iterrows():
        assert t.loc[row.parent, "t2"] == row.t1 - 1  # parent ends the frame before daughters start


def test_shapes_are_noncircular_and_masks_intact():
    cfg = small_cfg()
    cells = simulate_cells(cfg, np.random.default_rng(4))
    assert (cells.bound_radius_um > cells.radius_um * 1.05).mean() > 0.5  # clearly not circles
    _, lab = render_nuclei(cells, cfg, np.random.default_rng(4))
    for t, g in cells.groupby("t"):
        assert set(np.unique(lab[t])) - {0} == set(g.cell_id)  # every cell rendered, none overwritten away


def test_zero_deformation_ellipse_has_expected_area():
    cfg = small_cfg(**{"geometry.shape.deform_amp": 0.0, "geometry.shape.aspect_std": 0.0})
    cells = simulate_cells(cfg, np.random.default_rng(4))
    _, lab = render_nuclei(cells, cfg, np.random.default_rng(4))
    expected = np.pi * (cfg["geometry"]["cell_radius_um"] / (cfg["optics"]["pixel_size_nm"] / 1000)) ** 2
    areas = [(lab[0] == i).sum() for i in cells[cells.t == 0].cell_id]
    assert np.allclose(areas, expected, rtol=0.1)  # area-preserving ellipse


def test_texture_adds_structure_but_keeps_mean_brightness():
    base = {"acquisition.n_frames": 3, "cells.n_cells": 4, "optics.bleach_tau_s": "null", "optics.read_noise_e": 0.0,
            "optics.background_photons": 0, "optics.photons_per_px_s": 20000, "cells.p_divide_per_frame": 0.0}
    cfg_u = small_cfg(**base, **{"nucleus_texture.contrast": 0.0, "nucleus_texture.n_nucleoli": 0.0})
    cfg_t = small_cfg(**base, **{"nucleus_texture.contrast": 0.4, "nucleus_texture.n_nucleoli": 2.0})
    cells = simulate_cells(cfg_u, np.random.default_rng(1))
    img_u, lab = render_nuclei(cells, cfg_u, np.random.default_rng(2))
    img_t, _ = render_nuclei(cells, cfg_t, np.random.default_rng(2))
    inner = lab[0] > 0
    from scipy.ndimage import binary_erosion

    inner = binary_erosion(inner, iterations=3)  # stay away from the blurred edge
    cv = lambda im: np.std(im[0][inner] / np.mean(im[0][inner]))
    assert cv(img_t) > cv(img_u) * 2, "texture should add intra-nuclear variation"
    assert abs(img_t[0][lab[0] > 0].mean() / img_u[0][lab[0] > 0].mean() - 1) < 0.15  # same overall brightness
    # texture is attached to the nucleus: identical config and seed reproduce the same image
    img_t2, _ = render_nuclei(cells, cfg_t, np.random.default_rng(2))
    assert np.array_equal(img_t, img_t2)


def test_ground_truth_scores_perfectly_against_itself(tmp_path):
    cfg = small_cfg()
    cells = simulate_cells(cfg, np.random.default_rng(2))
    _, lab = render_nuclei(cells, cfg, np.random.default_rng(3))
    gt = write_ctc(lab, gt_table(cells), tmp_path / "gt")
    pred = write_ctc(lab, gt_table(cells), tmp_path / "pred")
    m = evaluate(gt, pred)
    assert m["CTCMetrics"]["TRA"] == 1.0 and m["CTCMetrics"]["DET"] == 1.0


def test_cells_can_walk_out_of_the_field_and_back(tmp_path):
    """geometry.outside_um: cells live in a larger world; the CTC truth splits a cell's stretches in view into linked tracks."""
    cfg = small_cfg(**{"acquisition.n_frames": 60, "geometry.outside_um": 30, "motion.speed_scale": 3.0})
    cells = simulate_cells(cfg, np.random.default_rng(4))
    fov = np.array(cfg["geometry"]["fov_um"])
    assert (cells.y_um < 0).any() or (cells.x_um < 0).any() or (cells.y_um > fov[0]).any() or (cells.x_um > fov[1]).any()
    _, lab = render_nuclei(cells, cfg, np.random.default_rng(3))
    masks, tab = ctc_ground_truth(lab, cells)
    for r in tab.itertuples():           # every track is present in every frame from t1 to t2 (what traccuracy requires)
        assert all((masks[t] == r.label).any() for t in range(r.t1, r.t2 + 1))
    gt = write_ctc(masks, tab, tmp_path / "gt")
    m = evaluate(gt, write_ctc(masks, tab, tmp_path / "pred"))
    assert m["CTCMetrics"]["DET"] == 1.0


def test_cell_layout_is_redrawn_instead_of_failing_and_unchanged_when_it_fit():
    from simlive.stage2_cells.motion import place_cells, placement_problem

    biplane = REPO_ROOT / "configs" / "loci_demo.yaml"
    # 4 nuclei in the 40 um biplane field: one at a time at random used to block the last one (run aborted); a 2 x 2 layout fits
    assert placement_problem(load_config(biplane, ["instrument=vutara_vxl_biplane", "cells.n_cells=4", "seed=71"])) is None
    msg = placement_problem(load_config(biplane, ["instrument=vutara_vxl_biplane", "cells.n_cells=8", "seed=71"]))
    assert msg and "Could not place 8 nuclei" in msg and "40 x 40 um" in msg

    # a layout that fitted on the first try is exactly the old one (same random numbers in the same order)
    bound, lo, hi = np.full(6, 7.0), np.zeros(2), np.array([80.0, 80.0])
    old_rng, new_rng = np.random.default_rng(5), np.random.default_rng(5)
    old = np.zeros((6, 2))
    for i in range(6):
        for _ in range(2000):
            p = old_rng.uniform(lo + bound[i], hi - bound[i])
            if all(np.hypot(*(p - old[j])) >= bound[i] + bound[j] for j in range(i)):
                old[i] = p
                break
    assert np.array_equal(place_cells(bound, lo, hi, new_rng), old)
    assert old_rng.uniform() == new_rng.uniform()          # and the stream continues from the same place
