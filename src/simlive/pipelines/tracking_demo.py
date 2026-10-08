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
from simlive.stage1_chromatin.library_loci import simulate_library_loci
from simlive.stage1_chromatin.toy_loci import simulate_toy_loci
from simlive.stage3_microscopy.render import render_loci, render_nuclei
from simlive.stage4_segtrack.segment import segment
from simlive.stage5_linking.celldata import map_cells
from simlive.stage5_linking.isolate import isolate_cells
from simlive.stage6_analysis.localize import localize_loci, locus_distances, score_localization
from simlive.stage4_segtrack.trackers import make_tracker
from simlive.stage7_validation.identity import identity_summary, match_ids
from simlive.stage7_validation.tracking_metrics import evaluate, gt_table, write_ctc

VARIANTS = ("gt_masks", "auto_masks")


def run_tracking_demo(cfg: dict, out_root: Path, progress: Callable[[str], None] = print) -> Path:
    logging.getLogger("traccuracy").setLevel(logging.CRITICAL)
    seed = int(cfg["seed"])
    run = Path(out_root) / run_id_for(cfg)
    s2, s3, s4, s7 = (run / n for n in ("stage2_cells", "stage3_microscopy", "stage4_segtrack", "stage7_validation"))
    ss = np.random.SeedSequence(seed).spawn(4)  # independent streams: motion, rendering, loci, locus rendering
    progress(f"Run folder: {run}")

    progress("Stage 2: simulating cell motion and shape")
    stage_dir(run, "stage2_cells")
    simulate_cells(cfg, np.random.default_rng(ss[0])).to_csv(s2 / "cells.csv", index=False)
    write_provenance(s2, cfg, seed)

    progress("Stage 3: rendering nuclear channel")
    stage_dir(run, "stage3_microscopy")
    cells = pd.read_csv(s2 / "cells.csv")
    img, lab = render_nuclei(cells, cfg, np.random.default_rng(ss[1]))
    tifffile.imwrite(s3 / "nucleus.tif", _stored(img, cfg))
    tifffile.imwrite(s3 / "labels.tif", lab)
    if cfg.get("loci"):
        progress("Stage 1 (stand-in): simulating loci inside each nucleus; stage 3: rendering locus channel")
        stage_dir(run, "stage1_chromatin")
        simulate = simulate_library_loci if cfg["loci"].get("source") == "library" else simulate_toy_loci
        truth = simulate(cells, cfg, np.random.default_rng(ss[2]))
        truth.to_csv(run / "stage1_chromatin" / "loci_truth.csv", index=False)
        write_provenance(run / "stage1_chromatin", cfg, seed)
        for k, rk in enumerate(ss[3].spawn(int(cfg["loci"]["n_loci"]))):  # one image channel (colour) per locus
            tifffile.imwrite(s3 / f"locus{k}.tif", _stored(render_loci(truth, cfg, np.random.default_rng(rk), locus_id=k), cfg))
    write_provenance(s3, cfg, seed)

    progress("Stage 4: segmenting")
    stage_dir(run, "stage4_segtrack")
    img, gt_lab = tifffile.imread(s3 / "nucleus.tif"), tifffile.imread(s3 / "labels.tif")
    seg_cfg = dict(cfg["segmentation"]); method = seg_cfg.pop("method")
    if seg_cfg.get("diameter_px") == "auto":  # expected nuclear diameter in pixels, from the geometry settings
        seg_cfg["diameter_px"] = 2 * cfg["geometry"]["cell_radius_um"] / (cfg["optics"]["pixel_size_nm"] / 1000.0)
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
    if cfg.get("loci"):
        _isolate_and_analyze(cfg, run, progress)
    progress("Done")
    return run


ID_SOURCES = {"truth": ("stage3_microscopy", "labels.tif"),
              "tracked_gt_masks": ("stage4_segtrack", "tracked_gt_masks.tif"),
              "tracked_auto_masks": ("stage4_segtrack", "tracked_auto_masks.tif")}


def _isolate_and_analyze(cfg: dict, run: Path, progress) -> None:
    """Stage 5 (isolate each cell) + stage 6 (placeholder per-cell locus analysis), for each identity source."""
    iso = cfg.get("isolation", {})
    sources = iso.get("sources", ["truth", "tracked_auto_masks"])
    s3 = run / "stage3_microscopy"
    px_um = cfg["optics"]["pixel_size_nm"] / 1000.0
    images = {"nucleus": s3 / "nucleus.tif",
              **{f"locus{k}": s3 / f"locus{k}.tif" for k in range(int(cfg["loci"]["n_loci"]))}}
    truth = pd.read_csv(run / "stage1_chromatin" / "loci_truth.csv")
    ana = cfg.get("locus_analysis", {})
    scores = {}
    for src in sources:
        d, f = ID_SOURCES[src]
        progress(f"Stage 5: isolating cells (identity source: {src})")
        t0 = time.time()
        out5 = stage_dir(run, "stage5_cells") / src
        summ = isolate_cells(images, run / d / f, out5, margin=iso.get("margin", 0.3), n_workers=iso.get("n_workers"),
                             min_frames=iso.get("min_frames", 1))
        progress(f"Stage 6: locating loci in {len(summ)} cells in parallel ({time.time() - t0:.0f}s so far)")
        out6 = stage_dir(run, "stage6_analysis") / src
        out6.mkdir(parents=True, exist_ok=True)
        pos = map_cells(localize_loci, out5, n_workers=iso.get("n_workers"), px_um=px_um,
                        n_loci=cfg["loci"]["n_loci"], **ana)
        pos.to_csv(out6 / "loci_positions.csv", index=False)
        dist = locus_distances(pos, px_um) if len(pos) else pd.DataFrame(columns=["cell_id", "t", "dist_um"])
        dist.to_csv(out6 / "locus_distances.csv", index=False)
        if src == "truth":  # IDs are the true cell IDs, so detections can be scored against the true loci
            err, sc = score_localization(pos, truth, px_um)
            err.to_csv(out6 / "localization_errors.csv", index=False)
            scores[src] = sc
        scores.setdefault(src, {})["peak_worker_memory_mb"] = summ.attrs.get("peak_worker_mb")
    try:  # peak working set of this (main) process over the whole run, incl. segmentation/tracking input arrays
        import psutil
        scores["pipeline_main_process_peak_mb"] = round(psutil.Process().memory_info().peak_wset / 2**20, 1)
    except Exception:
        pass
    (run / "stage6_analysis" / "scores.json").write_text(json.dumps(scores, indent=2, default=str))
    write_provenance(run / "stage5_cells", cfg, int(cfg["seed"]))
    write_provenance(run / "stage6_analysis", cfg, int(cfg["seed"]))


def _stored(img: np.ndarray, cfg: dict) -> np.ndarray:
    """Cameras output integers: with image_dtype: uint16 store counts as 16-bit (half the memory/disk of float32)."""
    if cfg.get("image_dtype") == "uint16":
        return np.rint(img).clip(0, 65535).astype(np.uint16)
    return img
