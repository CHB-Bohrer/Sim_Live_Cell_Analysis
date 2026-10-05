"""Cells-only tracking test: simulate moving cells, render nuclei, segment, track, score against ground truth.

    scripts\\run.cmd python scripts\\run_tracking_demo.py
    scripts\\run.cmd python scripts\\run_tracking_demo.py --set motion.D_um2_s=0.05 --set seed=2

Tracking is run twice: on GROUND-TRUTH masks (isolates linking error) and on AUTOMATIC-segmentation masks
(realistic). Both are scored against the true cell IDs with Cell Tracking Challenge metrics.
Outputs in data/runs/<run_id>/ (one folder per stage; each has params.yaml, seed.txt, provenance.json).
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from simlive.io.runs import REPO_ROOT, load_config, run_id_for, stage_dir, write_provenance
from simlive.stage2_cells.motion import simulate_cells
from simlive.stage3_microscopy.render import render_nuclei
from simlive.stage4_segtrack.segment import segment
from simlive.stage4_segtrack.trackers import make_tracker
from simlive.stage7_validation.tracking_metrics import evaluate, gt_table, write_ctc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(REPO_ROOT / "configs" / "tracking_demo.yaml"))
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL", help="override, e.g. motion.D_um2_s=0.05")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "data" / "runs")
    a = ap.parse_args()

    cfg = load_config(a.config, a.set)
    seed = int(cfg["seed"])
    run = a.out / run_id_for(cfg)
    s2, s3, s4, s7 = (run / n for n in ("stage2_cells", "stage3_microscopy", "stage4_segtrack", "stage7_validation"))
    ss = np.random.SeedSequence(seed).spawn(2)  # independent streams: motion, rendering
    print(f"Run: {run}")

    # ---- stage 2: cell motion -> cells.csv ----
    stage_dir(run, "stage2_cells")
    simulate_cells(cfg, np.random.default_rng(ss[0])).to_csv(s2 / "cells.csv", index=False)
    write_provenance(s2, cfg, seed)

    # ---- stage 3: render (reads cells.csv) -> nucleus.tif, labels.tif ----
    stage_dir(run, "stage3_microscopy")
    cells = pd.read_csv(s2 / "cells.csv")
    img, lab = render_nuclei(cells, cfg, np.random.default_rng(ss[1]))
    tifffile.imwrite(s3 / "nucleus.tif", img)
    tifffile.imwrite(s3 / "labels.tif", lab)
    write_provenance(s3, cfg, seed)
    n_div = int((cells.groupby("cell_id").parent_id.first() > 0).sum() // 2)
    print(f"  {cells.cell_id.nunique()} cell IDs, {cfg['acquisition']['n_frames']} frames, "
          f"{img.shape[1]}x{img.shape[2]} px, {n_div} divisions")

    # ---- stage 4: segment + track (reads stage 3 files) ----
    stage_dir(run, "stage4_segtrack")
    img, gt_lab = tifffile.imread(s3 / "nucleus.tif"), tifffile.imread(s3 / "labels.tif")
    seg_cfg = dict(cfg["segmentation"]); method = seg_cfg.pop("method")
    auto_lab = segment(img, method, **seg_cfg)
    tifffile.imwrite(s4 / "auto_segmentation.tif", auto_lab)
    tcfg = dict(cfg["tracking"]); tname = tcfg.pop("tracker")
    tracker = make_tracker(tname, **tcfg)

    # ---- stage 7: score both variants ----
    stage_dir(run, "stage7_validation")
    gt_dir = write_ctc(gt_lab, gt_table(cells), s7 / "ctc_gt")
    metrics = {}
    for variant, masks in (("gt_masks", gt_lab), ("auto_masks", auto_lab)):
        t0 = time.time()
        res = tracker.track(img, masks)
        tifffile.imwrite(s4 / f"tracked_{variant}.tif", res.masks)
        res.table.to_csv(s4 / f"tracks_{variant}.csv", index=False)
        pred_dir = write_ctc(res.masks, res.table, s7 / f"ctc_pred_{variant}")
        metrics[variant] = evaluate(gt_dir, pred_dir)
        print(f"  tracked [{variant}] in {time.time() - t0:.0f}s")
    write_provenance(s4, cfg, seed)
    write_provenance(s7, cfg, seed)
    (s7 / "metrics.json").write_text(json.dumps(metrics, indent=2, default=str))

    print("\n=== Tracking accuracy vs ground truth ===")
    for variant, m in metrics.items():
        print(f"[{variant}]")
        for group, vals in m.items():
            if isinstance(vals, dict):
                print(f"  {group}: " + ", ".join(f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}"
                                                  for k, v in vals.items()))
    print(f"\nView:  scripts\\run.cmd python scripts\\view_run.py {run} --labels stage4_segtrack\\tracked_auto_masks.tif")


if __name__ == "__main__":
    main()
