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
from simlive.stage2_cells.motion import add_z, simulate_cells, z_enabled
from simlive.stage1_chromatin.library_loci import simulate_library_loci
from simlive.stage1_chromatin.probes import simulate_probe_occupancy
from simlive.stage1_chromatin.toy_loci import simulate_toy_loci
from simlive.stage1_chromatin.transcription import ms2_spots, simulate_transcription, traces_from_truth
from simlive.stage3_microscopy.calibration.run import run_calibration
from simlive.stage3_microscopy.render import Imager, cells_in_image, entropy_for, truth_in_image
from simlive.stage4_segtrack.segment import segment
from simlive.stage5_linking.celldata import map_cells
from simlive.stage5_linking.isolate import isolate_cells
from simlive.stage6_analysis.localize import localize_loci, locus_distances, score_localization
from simlive.stage6_analysis.localize3d import localize_loci_3d, locus_distances_3d, score_localization_3d
from simlive.stage4_segtrack.trackers import make_tracker
from simlive.stage7_validation.identity import identity_summary, match_ids
from simlive.stage7_validation.tracking_metrics import evaluate, gt_table, write_ctc

VARIANTS = ("gt_masks", "auto_masks")


def run_tracking_demo(cfg: dict, out_root: Path, progress: Callable[[str], None] = print) -> Path:
    logging.getLogger("traccuracy").setLevel(logging.CRITICAL)
    seed = int(cfg["seed"])
    run = Path(out_root) / run_id_for(cfg)
    s2, s3, s4, s7 = (run / n for n in ("stage2_cells", "stage3_microscopy", "stage4_segtrack", "stage7_validation"))
    ss = np.random.SeedSequence(seed).spawn(8)  # independent streams: motion, rendering, loci, (unused), probes, transcription, cell heights (z), calibration
    progress(f"Run folder: {run}")

    progress("Stage 2: simulating cell motion and shape")
    stage_dir(run, "stage2_cells")
    cells0 = simulate_cells(cfg, np.random.default_rng(ss[0]))
    if z_enabled(cfg):                                  # z-stacks: give every nucleus a height and a thickness (own random stream)
        cells0 = add_z(cells0, cfg, np.random.default_rng(ss[6]))
    cells0.to_csv(s2 / "cells.csv", index=False)
    write_provenance(s2, cfg, seed)

    cells = pd.read_csv(s2 / "cells.csv")
    truth = occ = spots = None
    if cfg.get("loci"):
        progress("Stage 1: simulating loci inside each nucleus")
        stage_dir(run, "stage1_chromatin")
        traces = None
        if cfg["loci"].get("source") == "library":
            truth, traces = simulate_library_loci(cells, cfg, np.random.default_rng(ss[2]), return_traces=True)
        else:
            truth = simulate_toy_loci(cells, cfg, np.random.default_rng(ss[2]))
        truth.to_csv(run / "stage1_chromatin" / "loci_truth.csv", index=False)
        write_provenance(run / "stage1_chromatin", cfg, seed)
        if (cfg["loci"].get("probes") or {}).get("enabled"):
            progress("Stage 1b: stochastic probe binding at the promoter and enhancer (sets each locus' brightness)")
            occ = simulate_probe_occupancy(cells, cfg, np.random.default_rng(ss[4]))
            occ.to_csv(run / "stage1_chromatin" / "probe_occupancy.csv", index=False)
        if (cfg["loci"].get("transcription") or {}).get("enabled"):
            progress("Stage 1c: stochastic transcription driven by promoter-enhancer distance (the MS2 signal)")
            tx, tx_events = simulate_transcription(cells, traces or traces_from_truth(truth, cfg), cfg,
                                                   np.random.default_rng(ss[5]))
            tx.to_csv(run / "stage1_chromatin" / "transcription.csv", index=False)
            tx_events.to_csv(run / "stage1_chromatin" / "transcription_events.csv", index=False)
            spots = ms2_spots(tx, truth, cfg)

    stage_dir(run, "stage3_microscopy")
    imager = Imager(cfg, cells, truth, spots, occ, entropy_for(seed))
    on = imager.E.enabled_ids()
    progress("Stage 3: rendering " + " + ".join(imager.channels) + (" (MS2 spots are in the nuclear channel)" if spots is not None else "")
             + f" with {len(on)} image-error sources on: {', '.join(on)}")
    mv = imager.render_movie(progress=progress, dtype=np.uint16 if cfg.get("image_dtype") == "uint16" else np.float32)
    for name, arr in mv["images"].items():                                   # one image channel (colour) per locus + the nucleus
        if imager.is3d:        # full stacks (T, Z, Y, X) for 3D analysis; the 2D projection feeds segmentation, tracking and isolation
            tifffile.imwrite(s3 / f"{name}_zstack.tif", arr, photometric="minisblack")
            arr = imager.project(arr, imager.projection).astype(arr.dtype)
        tifffile.imwrite(s3 / f"{name}.tif", arr, photometric="minisblack")
    if imager.is3d:
        tifffile.imwrite(s3 / "labels_3d.tif", mv["labels"], photometric="minisblack")
        mv["labels"] = mv["labels"].max(axis=1)                              # nuclei never overlap laterally, so max over z keeps the IDs
    tifffile.imwrite(s3 / "labels.tif", mv["labels"], photometric="minisblack")
    # ground truth of the imaging itself: what drift / focus / power / blur each frame really had, and the truth in image coordinates
    it = imager.imaging_truth()
    it.to_csv(s3 / "imaging_truth.csv", index=False)
    px_um = cfg["optics"]["pixel_size_nm"] / 1000.0
    st_truth = imager.static_truth()
    cells_in_image(cells, it, px_um).to_csv(s3 / "cells_image.csv", index=False)
    if truth is not None:
        truth_in_image(truth, it, px_um, st_truth["channel_shift_px"],
                       (st_truth.get("z_stack") or {}).get("channel_shift_z_um")).to_csv(s3 / "loci_truth_image.csv", index=False)
    (s3 / "imaging_errors.json").write_text(json.dumps(st_truth, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    del mv
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
        # calibration runs here, after the tracker: torch is already loaded, so numpy's LAPACK calls cannot break its import
        cal_dir = _run_calibration(cfg, run, progress) if (cfg.get("calibration") or {}).get("simulate") else None
        _isolate_and_analyze(cfg, run, progress, cal_dir)
    progress("Done")
    return run


ID_SOURCES = {"truth": ("stage3_microscopy", "labels.tif"),
              "tracked_gt_masks": ("stage4_segtrack", "tracked_gt_masks.tif"),
              "tracked_auto_masks": ("stage4_segtrack", "tracked_auto_masks.tif")}


def _run_calibration(cfg: dict, run: Path, progress) -> Path:
    """Simulate dark / flat-field / bead acquisitions of this run's instrument, estimate the instrument from them, save the calibration."""
    progress("Calibration: simulating dark, flat-field and bead acquisitions of this instrument and estimating it from them")
    t0 = time.time()
    cal, table = run_calibration(cfg, entropy_for(int(cfg["seed"])), progress=progress)
    d = run / "stage3_microscopy" / "calibration"
    cal.save(d)
    table.to_csv(d / "comparison.csv", index=False)
    progress(f"Calibration: done in {time.time() - t0:.0f}s ({len(table)} quantities compared with the truth)")
    return d


def _isolate_and_analyze(cfg: dict, run: Path, progress, cal_dir: Path | None = None) -> None:
    """Stage 5 (isolate each cell) + stage 6 (placeholder per-cell locus analysis), for each identity source."""
    iso = cfg.get("isolation", {})
    sources = iso.get("sources", ["truth", "tracked_auto_masks"])
    s3 = run / "stage3_microscopy"
    px_um = cfg["optics"]["pixel_size_nm"] / 1000.0
    images = {"nucleus": s3 / "nucleus.tif",
              **{f"locus{k}": s3 / f"locus{k}.tif" for k in range(int(cfg["loci"]["n_loci"]))}}
    ti = s3 / "loci_truth_image.csv"        # truth moved by the stage drift (a localizer cannot know it); chromatic shift stays an error
    truth = pd.read_csv(ti if ti.exists() else run / "stage1_chromatin" / "loci_truth.csv")
    ej = s3 / "imaging_errors.json"
    zinfo = json.loads(ej.read_text()).get("z_stack") if ej.exists() else None
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
        apply_cal = cal_dir is not None and (cfg.get("calibration") or {}).get("apply", True)
        if apply_cal:           # the same analysis on dark/flat/hot-pixel-corrected crops with the chromatic shift removed
            progress("Stage 6: locating loci again with the calibration applied")
            posc = map_cells(localize_loci, out5, n_workers=iso.get("n_workers"), px_um=px_um, n_loci=cfg["loci"]["n_loci"],
                             calib=str(cal_dir), **ana)
            posc.to_csv(out6 / "loci_positions_cal.csv", index=False)
            (locus_distances(posc, px_um) if len(posc) else pd.DataFrame(columns=["cell_id", "t", "dist_um"])).to_csv(
                out6 / "locus_distances_cal.csv", index=False)
            if src == "truth":
                errc, scc = score_localization(posc, truth, px_um)
                errc.to_csv(out6 / "localization_errors_cal.csv", index=False)
                scores["truth_calibrated"] = scc
        if zinfo is not None:   # z-stacks: locate the loci in 3D from the full stacks (the 2D result above uses the projections)
            progress(f"Stage 6: locating loci in 3D in {len(summ)} cells")
            stacks = {f"locus{k}": str(s3 / f"locus{k}_zstack.tif") for k in range(int(cfg["loci"]["n_loci"]))}
            pos3 = map_cells(localize_loci_3d, out5, n_workers=iso.get("n_workers"), stacks=stacks, z_um=zinfo["z_positions_um"],
                             n_loci=cfg["loci"]["n_loci"], **ana)
            pos3.to_csv(out6 / "loci_positions_3d.csv", index=False)
            (locus_distances_3d(pos3, px_um) if len(pos3) else pd.DataFrame(columns=["cell_id", "t", "dist_um"])).to_csv(
                out6 / "locus_distances_3d.csv", index=False)
            if src == "truth" and "z_um" in truth.columns:
                err3, sc3 = score_localization_3d(pos3, truth, px_um)
                err3.to_csv(out6 / "localization_errors_3d.csv", index=False)
                scores["truth_3d"] = sc3
            if apply_cal:
                pos3c = map_cells(localize_loci_3d, out5, n_workers=iso.get("n_workers"), stacks=stacks, z_um=zinfo["z_positions_um"],
                                  n_loci=cfg["loci"]["n_loci"], calib=str(cal_dir), **ana)
                pos3c.to_csv(out6 / "loci_positions_3d_cal.csv", index=False)
                (locus_distances_3d(pos3c, px_um) if len(pos3c) else pd.DataFrame(columns=["cell_id", "t", "dist_um"])).to_csv(
                    out6 / "locus_distances_3d_cal.csv", index=False)
                if src == "truth" and "z_um" in truth.columns:
                    err3c, sc3c = score_localization_3d(pos3c, truth, px_um)
                    err3c.to_csv(out6 / "localization_errors_3d_cal.csv", index=False)
                    scores["truth_3d_calibrated"] = sc3c
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
