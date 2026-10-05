"""CPU tests for stage 2 (motion), stage 3 (render) and the tracking scorer."""
import numpy as np

from simlive.io.runs import load_config, REPO_ROOT
from simlive.stage2_cells.motion import simulate_cells
from simlive.stage3_microscopy.render import render_nuclei
from simlive.stage7_validation.tracking_metrics import evaluate, gt_table, write_ctc

CFG = REPO_ROOT / "configs" / "tracking_demo.yaml"


def small_cfg(**over):
    return load_config(CFG, ["acquisition.n_frames=12", "cells.n_cells=8", *[f"{k}={v}" for k, v in over.items()]])


def test_motion_is_deterministic_and_collision_free():
    cfg = small_cfg()
    a = simulate_cells(cfg, np.random.default_rng(5))
    b = simulate_cells(cfg, np.random.default_rng(5))
    assert a.equals(b)
    for _, g in a.groupby("t"):
        p, r = g[["y_um", "x_um"]].to_numpy(), g["radius_um"].to_numpy()
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


def test_ground_truth_scores_perfectly_against_itself(tmp_path):
    cfg = small_cfg()
    cells = simulate_cells(cfg, np.random.default_rng(2))
    _, lab = render_nuclei(cells, cfg, np.random.default_rng(3))
    gt = write_ctc(lab, gt_table(cells), tmp_path / "gt")
    pred = write_ctc(lab, gt_table(cells), tmp_path / "pred")
    m = evaluate(gt, pred)
    assert m["CTCMetrics"]["TRA"] == 1.0 and m["CTCMetrics"]["DET"] == 1.0
