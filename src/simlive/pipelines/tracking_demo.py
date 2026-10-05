"""Cells-only tracking test as a function (used by scripts/run_tracking_demo.py and the dashboard).

stage 2 (motion) -> stage 3 (nuclear channel) -> stage 4 (segment + track, twice) -> stage 7 (score).
Each stage reads the previous stage's files from disk.
"""
from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import tifffile

from simlive.io.runs import run_id_for, stage_dir, write_provenance
from simlive.stage2_cells.motion import simulate_cells
from simlive.stage3_microscopy.render import render_nuclei
from simlive.stage4_segtrack.segment import segment
from simlive.stage4_segtrack.trackers import make_tracker
from simlive.stage7_validation.identity import identity_summary, match_ids
from simlive.stage7_validation.tracking_metrics import evaluate, gt_table, write_ctc

VARIANTS = ("gt_masks", "auto_masks")


def run_tracking_demo(cfg: dict, out_root: Path, progress: Callable[[str], None] = print) -> Path:
    logging.getLogger("traccuracy").setLevel(logging.CRITICAL)
    seed = int(cfg["seed"])
    run = Path(out_root) / run_id_for(cfg)
    s2, s3, s4, s7 = (run / n for n in ("stage2_cells", "stage3_microscopy", "stage4_segtrack", "stage7_validation"))
    ss = np.random.SeedSequence(seed).spawn(2)  # independent streams: motion, rendering
    progress(f"Run folder: {run}")

    progress("Stage 2: simulating cell motion and shape")
    stage_dir(run, "stage2_cells")
    simulate_cells(cfg, np.random.default_rng(ss[0])).to_csv(s2 / "cells.csv", index=False)
    write_provenance(s2, cfg, seed)

    progress("Stage 3: rendering nuclear channel")
    stage_dir(run, "stage3_microscopy")
    cells = pd.read_csv(s2 / "cells.csv")
    img, lab = render_nuclei(cells, cfg, np.random.default_rng(ss[1]))
    tifffile.imwrite(s3 / "nucleus.tif", img)
    tifffile.imwrite(s3 / "labels.tif", lab)
    write_provenance(s3, cfg, seed)

    progress("Stage 4: segmenting")
    stage_dir(run, "stage4_segtrack")
    img, gt_lab = tifffile.imread(s3 / "nucleus.tif"), tifffile.imread(s3 / "labels.tif")
    seg_cfg = dict(cfg["segmentation"]); method = seg_cfg.pop("method")
    auto_lab = segment(img, method, **seg_cfg)
    tifffile.imwrite(s4 / "auto_segmentation.tif", auto_lab)
    tcfg = dict(cfg["tracking"]); tname = tcfg.pop("tracker")
    tracker = make_tracker(tname, **tcfg)

    stage_dir(run, "stage7_validation")
    gt_dir = write_ctc(gt_lab, gt_table(cells), s7 / "ctc_gt")
    metrics = {}
    for variant, masks in (("gt_masks", gt_lab), ("auto_masks", auto_lab)):
        progress(f"Stage 4: tracking on {variant} ({tname})")
        t0 = time.time()
        res = tracker.track(img, masks)
        tifffile.imwrite(s4 / f"tracked_{variant}.tif", res.masks)
        res.table.to_csv(s4 / f"tracks_{variant}.csv", index=False)
        progress(f"Stage 7: scoring {variant} (tracking took {time.time() - t0:.0f}s)")
        pred_dir = write_ctc(res.masks, res.table, s7 / f"ctc_pred_{variant}")
        m = evaluate(gt_dir, pred_dir)
        matches = match_ids(gt_lab, res.masks)
        totals, per, events = identity_summary(matches)
        matches.to_csv(s7 / f"id_matches_{variant}.csv", index=False)
        per.to_csv(s7 / f"id_per_cell_{variant}.csv", index=False)
        events.to_csv(s7 / f"id_switch_events_{variant}.csv", index=False)
        m["Identity"] = totals
        metrics[variant] = m
    write_provenance(s4, cfg, seed)
    write_provenance(s7, cfg, seed)
    (s7 / "metrics.json").write_text(json.dumps(metrics, indent=2, default=str))
    progress("Done")
    return run
