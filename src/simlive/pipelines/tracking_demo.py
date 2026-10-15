"""Cells-only tracking test as a function (used by scripts/run_tracking_demo.py and the dashboard).

stage 2 (motion) -> stage 3 (nuclear channel) -> stage 4 (segment + track, twice) -> stage 7 (score).
Each stage reads the previous stage's files from disk.
"""
from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass
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
from simlive.stage3_microscopy.fiducials import make_fiducials
from simlive.stage3_microscopy.render import Imager, entropy_for, save_imaging_truth
from simlive.stage4_segtrack.segment import segment
from simlive.stage5_linking.celldata import map_cells
from simlive.stage5_linking.isolate import isolate_cells
from simlive.stage6_analysis.localize import localize_loci, locus_distances, photons_in_frame, score_localization
from simlive.stage6_analysis.fiducials import (correct_with_beads, correct_with_plane_positions, score_variants, track_fiducials)
from simlive.stage6_analysis.localize3d import localize_loci_3d, locus_distances_3d, score_localization_3d
from simlive.stage6_analysis.ms2 import ms2_trace, ms2_truth, score_ms2
from simlive.stage4_segtrack.trackers import make_tracker
from simlive.stage7_validation.identity import identity_summary, match_ids
from simlive.stage7_validation.propagation import VARIANT_OF_SOURCE, propagation_scores
from simlive.stage7_validation.tracking_metrics import ctc_ground_truth, evaluate, write_ctc

VARIANTS = ("gt_masks", "auto_masks")


def run_tracking_demo(cfg: dict, out_root: Path, progress: Callable[[str], None] = print) -> Path:
    logging.getLogger("traccuracy").setLevel(logging.CRITICAL)
    seed = int(cfg["seed"])
    run = Path(out_root) / run_id_for(cfg)
    s2, s3, s4, s7 = (run / n for n in ("stage2_cells", "stage3_microscopy", "stage4_segtrack", "stage7_validation"))
    ss = np.random.SeedSequence(seed).spawn(8)  # independent streams: motion, rendering, loci, (unused), probes, transcription, cell heights (z), fiducial beads
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
    fid = make_fiducials(cfg, cells, np.random.default_rng(ss[7]))          # beads on the slide in every image (None if off)
    if fid is not None:
        fid.to_csv(s3 / "fiducials.csv", index=False)
    imager = Imager(cfg, cells, truth, spots, occ, entropy_for(seed), fiducials=fid)
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
    st_truth = save_imaging_truth(imager, cells, truth, s3)
    if imager.is3d and imager.E.on("piezo"):               # where the planes really were, and what the piezo's position sensor logged
        imager.piezo_truth().to_csv(s3 / "piezo_truth.csv", index=False)
    (s3 / "imaging_errors.json").write_text(json.dumps(st_truth, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    del mv
    write_provenance(s3, cfg, seed)

    if (cfg.get("tracking") or {}).get("skip"):          # quick runs (sweeps of the imaging / spot analysis): true IDs only
        progress("Stage 4: skipped (tracking.skip); the analysis uses the true cell IDs only")
        if cfg.get("loci"):
            cal_dir = _run_calibration(cfg, run, progress) if (cfg.get("calibration") or {}).get("simulate") else None
            cfg = {**cfg, "isolation": {**(cfg.get("isolation") or {}), "sources": ["truth"]}}
            _isolate_and_analyze(cfg, run, progress, cal_dir)
        progress("Done")
        return run

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
    gt_dir = write_ctc(*ctc_ground_truth(gt_lab, cells), s7 / "ctc_gt")      # cells that leave and re-enter: one CTC track per stretch
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


def nominal_axial_sigma(cfg: dict) -> dict:
    """Axial PSF sigma (um) of every locus colour from the microscope's specification (NA, wavelength, modality): what an analysis may know."""
    from simlive.stage3_microscopy.render import axial_sigma_um

    return {f"locus{k}": axial_sigma_um(cfg["optics"], wl) for k, wl in enumerate(cfg["loci"]["wavelengths_nm"][:int(cfg["loci"]["n_loci"])])}


EMPTY_DIST = ["cell_id", "t", "dist_um"]


@dataclass
class _Analysis:
    """What every locus-analysis step of `_isolate_and_analyze` needs: the run, the instrument numbers and the (identity-independent) truth."""
    cfg: dict
    run: Path
    progress: Callable[[str], None]
    cal_dir: Path | None
    px_um: float
    truth: pd.DataFrame
    zinfo: dict | None
    ana: dict
    iso: dict

    @property
    def s3(self) -> Path:
        return self.run / "stage3_microscopy"

    @property
    def n_workers(self):
        return self.iso.get("n_workers")

    @property
    def n_loci(self) -> int:
        return int(self.cfg["loci"]["n_loci"])

    @property
    def apply_cal(self) -> bool:
        return self.cal_dir is not None and bool((self.cfg.get("calibration") or {}).get("apply", True))

    @property
    def has_transcription(self) -> bool:
        return (self.run / "stage1_chromatin" / "transcription.csv").exists()


def _write_locus_results(A: _Analysis, out6: Path, suffix: str, pos: pd.DataFrame, dist_fn, score_fn=None) -> dict | None:
    """Write loci_positions<suffix>.csv and locus_distances<suffix>.csv (and localization_errors<suffix>.csv when `score_fn` is given);
    returns the scores. Suffixes: '' (2D), '_cal' (calibration applied), '_3d', '_3d_cal'."""
    pos.to_csv(out6 / f"loci_positions{suffix}.csv", index=False)
    (dist_fn(pos, A.px_um) if len(pos) else pd.DataFrame(columns=EMPTY_DIST)).to_csv(out6 / f"locus_distances{suffix}.csv", index=False)
    if score_fn is None:
        return None
    err, sc = score_fn(pos, A.truth, A.px_um)
    err.to_csv(out6 / f"localization_errors{suffix}.csv", index=False)
    return sc


def _locate_2d(A: _Analysis, src: str, out5: Path, out6: Path, scores: dict) -> None:
    """2D spot finding on every isolated cell, then (identity = truth) the score against the true loci, or (tracked IDs) stage 7b."""
    pos = map_cells(localize_loci, out5, n_workers=A.n_workers, px_um=A.px_um, n_loci=A.n_loci, **A.ana)
    sc = _write_locus_results(A, out6, "", pos, locus_distances, score_localization if src == "truth" else None)
    if sc is not None:
        scores[src] = sc
    posc = None
    if A.apply_cal:           # the same analysis on dark/flat/hot-pixel-corrected crops with the chromatic shift removed
        A.progress("Stage 6: locating loci again with the calibration applied")
        posc = map_cells(localize_loci, out5, n_workers=A.n_workers, px_um=A.px_um, n_loci=A.n_loci, calib=str(A.cal_dir), **A.ana)
        scc = _write_locus_results(A, out6, "_cal", posc, locus_distances, score_localization if src == "truth" else None)
        if scc is not None:
            scores["truth_calibrated"] = scc
    if src != "truth":
        _propagate(A, src, out6, pos, posc, scores)


def _propagate(A: _Analysis, src: str, out6: Path, pos: pd.DataFrame, posc: pd.DataFrame | None, scores: dict) -> None:
    """Stage 7b: score the locus detections of a tracked identity source through the tracker's IDs (see stage7_validation/propagation.py)."""
    mf = A.run / "stage7_validation" / f"id_matches_{VARIANT_OF_SOURCE.get(src, '')}.csv"
    if not mf.exists():
        return
    matches = pd.read_csv(mf)
    out = scores.setdefault("propagation", {}).setdefault(src, {})
    for key, p, stem in (("uncalibrated", pos, ""), ("calibrated", posc, "_cal")):
        if p is None:
            continue
        err, sc = propagation_scores(p, matches, A.truth, A.px_um)
        err.to_csv(out6 / f"localization_errors_via_tracks{stem}.csv", index=False)
        out[key] = sc


def _trace_ms2(A: _Analysis, src: str, out5: Path, out6: Path, scores: dict) -> None:
    """The MS2 spot in the nuclear channel: Chen et al. 2025 intensity trace per cell (with and without the calibration)."""
    A.progress("Stage 6: finding the MS2 spot and its intensity trace in every cell")
    ms2_t = ms2_truth(A.run) if src == "truth" else None
    kw = {k: v for k, v in A.ana.items() if k != "dilate_px"}
    for tag, cdir in (("", None),) + ((("_cal", str(A.cal_dir)),) if A.apply_cal else ()):
        tr = map_cells(ms2_trace, out5, n_workers=A.n_workers, calib=cdir, **kw)
        tr.to_csv(out6 / f"ms2_traces{tag}.csv", index=False)
        if ms2_t is not None and len(tr):
            scores["ms2" + ("_calibrated" if tag else "")] = score_ms2(tr, ms2_t, A.px_um)[1]


def _locate_3d(A: _Analysis, src: str, out5: Path, out6: Path, n_cells: int, scores: dict) -> pd.DataFrame:
    """z-stacks: locate the loci in 3D from the full stacks (the 2D result uses the projections). Returns the 3D positions to correct with the beads
    (calibrated if available)."""
    A.progress(f"Stage 6: locating loci in 3D in {n_cells} cells")
    stacks = {f"locus{k}": str(A.s3 / f"locus{k}_zstack.tif") for k in range(A.n_loci)}
    ana3 = {**A.ana, "axial_sigma_um": nominal_axial_sigma(A.cfg), "px_um": A.px_um}   # nominal axial PSF (median_z: auto); px for z scaling
    z = A.zinfo["z_positions_um"]
    score3 = score_localization_3d if src == "truth" and "z_um" in A.truth.columns else None
    pos3 = map_cells(localize_loci_3d, out5, n_workers=A.n_workers, stacks=stacks, z_um=z, n_loci=A.n_loci, **ana3)
    sc = _write_locus_results(A, out6, "_3d", pos3, locus_distances_3d, score3)
    if sc is not None:
        scores["truth_3d"] = sc
    if not A.apply_cal:
        return pos3
    pos3c = map_cells(localize_loci_3d, out5, n_workers=A.n_workers, stacks=stacks, z_um=z, n_loci=A.n_loci, calib=str(A.cal_dir), **ana3)
    scc = _write_locus_results(A, out6, "_3d_cal", pos3c, locus_distances_3d, score3)
    if scc is not None:
        scores["truth_3d_calibrated"] = scc
    return pos3c


def _correct_with_fiducials(A: _Analysis, out6: Path, basepos: pd.DataFrame, scores: dict) -> None:
    """Track the fiducial beads (and the piezo log) and score the corrected heights (identity = truth, z-stacks with beads only)."""
    A.progress("Stage 6: tracking the fiducial beads and correcting the heights")
    s3, z = A.s3, A.zinfo["z_positions_um"]
    seg = A.run / "stage4_segtrack" / "auto_segmentation.tif"
    if not seg.exists():                 # tracking.skip: no automatic segmentation was made, the true nuclei say where the beads are NOT
        seg = s3 / "labels.tif"
    shifts, _ = track_fiducials(s3 / "nucleus_zstack.tif", z, seg)
    shifts.to_csv(out6 / "fiducial_shifts.csv", index=False)
    variants = {"uncorrected": basepos, "fiducial beads": correct_with_beads(basepos, shifts)}
    pt = s3 / "piezo_truth.csv"
    if pt.exists():
        log = pd.read_csv(pt)
        variants["piezo position log"] = correct_with_plane_positions(basepos, log, z)
        zlog = {int(t): g.sort_values("k").readout_um.to_numpy() for t, g in log.groupby("t")}
        shifts_l, _ = track_fiducials(s3 / "nucleus_zstack.tif", z, seg, z_per_frame=zlog)
        variants["beads + piezo log"] = correct_with_beads(variants["piezo position log"], shifts_l)
    tab = score_variants(variants, pd.read_csv(A.run / "stage1_chromatin" / "loci_truth.csv"), A.px_um)
    tab.to_csv(out6 / "fiducial_corrections.csv", index=False)
    scores["fiducial_corrections"] = tab.to_dict("records")


def _isolate_and_analyze(cfg: dict, run: Path, progress, cal_dir: Path | None = None) -> None:
    """Stage 5 (isolate each cell) + stage 6 (per-cell locus analysis) + stage 7b (tracked IDs vs the truth), for each identity source."""
    iso = cfg.get("isolation", {})
    sources = iso.get("sources", ["truth", "tracked_auto_masks"])
    s3 = run / "stage3_microscopy"
    px_um = cfg["optics"]["pixel_size_nm"] / 1000.0
    images = {"nucleus": s3 / "nucleus.tif",
              **{f"locus{k}": s3 / f"locus{k}.tif" for k in range(int(cfg["loci"]["n_loci"]))}}
    ti = s3 / "loci_truth_image.csv"        # truth moved by the stage drift (a localizer cannot know it); chromatic shift stays an error
    truth = pd.read_csv(ti if ti.exists() else run / "stage1_chromatin" / "loci_truth.csv")
    itp = s3 / "imaging_truth.csv"
    truth = photons_in_frame(truth, pd.read_csv(itp) if itp.exists() else None)     # true photons per spot and frame, to score intensities
    ej = s3 / "imaging_errors.json"
    zinfo = json.loads(ej.read_text()).get("z_stack") if ej.exists() else None
    A = _Analysis(cfg, run, progress, cal_dir, px_um, truth, zinfo, cfg.get("locus_analysis", {}), iso)
    scores = {}
    for src in sources:
        d, f = ID_SOURCES[src]
        progress(f"Stage 5: isolating cells (identity source: {src})")
        t0 = time.time()
        out5 = stage_dir(run, "stage5_cells") / src
        summ = isolate_cells(images, run / d / f, out5, margin=iso.get("margin", 0.3), n_workers=A.n_workers,
                             min_frames=iso.get("min_frames", 1))
        progress(f"Stage 6: locating loci in {len(summ)} cells in parallel ({time.time() - t0:.0f}s so far)")
        out6 = stage_dir(run, "stage6_analysis") / src
        out6.mkdir(parents=True, exist_ok=True)
        _locate_2d(A, src, out5, out6, scores)
        if A.has_transcription:
            _trace_ms2(A, src, out5, out6, scores)
        if zinfo is not None:
            basepos = _locate_3d(A, src, out5, out6, len(summ), scores)
            if src == "truth" and (s3 / "fiducials.csv").exists() and "z_um" in truth.columns:
                _correct_with_fiducials(A, out6, basepos, scores)
        scores.setdefault(src, {})["peak_worker_memory_mb"] = summ.attrs.get("peak_worker_mb")
    progress("Stage 6: photobleaching correction (background removed first; added as intensity_corr)")
    from simlive.stage6_analysis.bleach import bleach_correct_run

    scores["bleaching"] = bleach_correct_run(run, tuple(sources))
    if A.has_transcription:      # MS2 trace -> Pol II initiation events and bursts (Chen et al. 2025)
        from simlive.stage6_analysis.deconvolve import deconvolve_run

        progress("Stage 6: deconvolving the MS2 traces into Pol II initiation events (Chen et al. 2025)")
        scores["ms2_deconvolution"] = deconvolve_run(run, ("truth",), progress=lambda s: None)
    try:  # peak working set of this (main) process over the whole run, incl. segmentation/tracking input arrays
        import psutil
        scores["pipeline_main_process_peak_mb"] = round(psutil.Process().memory_info().peak_wset / 2**20, 1)
    except Exception:
        pass
    (run / "stage6_analysis" / "scores.json").write_text(json.dumps(scores, indent=2, default=str))
    write_provenance(run / "stage5_cells", cfg, int(cfg["seed"]))
    write_provenance(run / "stage6_analysis", cfg, int(cfg["seed"]))
