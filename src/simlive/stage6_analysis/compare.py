"""Compare locus-finding methods on a finished run: re-run only stage 6 on the saved single-cell crops (true cell IDs) with each variant
and score it against the simulation's truth. Used by scripts/compare_localizers.py and the dashboard's 🔎 Spot finding tab.

Output: <run>/stage6_analysis/method_comparison.csv (one row per variant, dimension and calibration) and
<run>/stage6_analysis/comparison/<variant>[_cal][_3d].csv (per-detection errors, for drill-down).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import yaml

from simlive.stage5_linking.celldata import map_cells
from simlive.stage6_analysis.localize import localize_loci, score_localization
from simlive.stage6_analysis.localize3d import localize_loci_3d, score_localization_3d

# key -> (label, settings that differ from the defaults in spots.DEFAULTS). Every median_peaks variant starts from the defaults (not from the
# run's own locus_analysis block, except dilate_px), so the table means the same thing for every run.
VARIANTS = {
    "original": ("Original placeholder: brightest pixel after smoothing + 5x5 centroid", {"method": "brightest"}),
    "paper_literal": ("As close to Bennabi 2025 as the scale allows: median 3x3, one background per cell, Otsu movie threshold on the "
                      "filtered image, area, smoothing, local maxima, brightest, integer pixel",
                      {"method": "median_peaks", "background": "cell", "threshold": "otsu", "threshold_on": "filtered", "select": "brightest",
                       "refine": "none"}),
    "default": ("Default: Bennabi 2025 steps + local background, noise threshold on the smoothed image, track choice, Gaussian fit",
                {"method": "median_peaks"}),
    "no_median": ("... without the median filter (k = 1)", {"method": "median_peaks", "median_px": 1}),
    "median5": ("... median 5x5", {"method": "median_peaks", "median_px": 5}),
    "cell_background": ("... one background value per cell instead of the local background", {"method": "median_peaks", "background": "cell"}),
    "otsu": ("... the paper's Otsu threshold instead of the noise threshold", {"method": "median_peaks", "threshold": "otsu"}),
    "on_filtered": ("... threshold on the median-filtered image (the paper's order)", {"method": "median_peaks", "threshold_on": "filtered"}),
    "brightest_choice": ("... brightest candidate instead of the track choice", {"method": "median_peaks", "select": "brightest"}),
    "per_frame": ("... each frame's own threshold instead of one per movie", {"method": "median_peaks", "threshold_scope": "frame"}),
    "integer": ("... integer peak position, no sub-pixel step (as the paper)", {"method": "median_peaks", "refine": "none"}),
    "centroid": ("... centroid instead of the Gaussian fit", {"method": "median_peaks", "refine": "centroid"}),
}
VARIANTS_3D = {
    "original": ("Original placeholder (3D): brightest voxel + centroid", {"method": "brightest"}),
    "default": ("Default in 3D (axial median: auto = only when the planes are closer than the axial PSF sigma)", {"method": "median_peaks"}),
    "zmed3": ("... axial median over 3 planes always", {"method": "median_peaks", "median_z": 3}),
    "brightest_choice": ("... brightest candidate instead of the track choice", {"method": "median_peaks", "select": "brightest"}),
    "brightest_zmed3": ("... brightest candidate with the axial median over 3 planes", {"method": "median_peaks", "select": "brightest", "median_z": 3}),
}
SCORE_COLS = ["detection_rate", "n_false_on_dark", "median_error_um", "p95_error_um", "rms_error_um", "frac_far", "rms_error_near_um"]


def run_config(run: Path) -> dict:
    return yaml.safe_load((Path(run) / "stage2_cells" / "params.yaml").read_text())


def truth_with_photons(run: Path) -> pd.DataFrame:
    """The image-space locus truth, with each spot's photon count (rebuilt from the saved ground truth for runs made before it was saved)."""
    run = Path(run)
    s3 = run / "stage3_microscopy"
    ti = s3 / "loci_truth_image.csv"
    truth = pd.read_csv(ti if ti.exists() else run / "stage1_chromatin" / "loci_truth.csv")
    if "photons" not in truth.columns:
        from simlive.stage3_microscopy.render import imager_for_run

        ph = imager_for_run(run, run_config(run)).locus_photons()
        truth = truth.merge(ph, on=["t", "cell_id", "locus_id"], how="left")
    return truth


def compare(run: Path, keys: list[str] | None = None, n_workers: int | None = None, calibrated: bool = True, three_d: bool = True,
            progress=print) -> pd.DataFrame:
    run = Path(run)
    cfg = run_config(run)
    s3, s6 = run / "stage3_microscopy", run / "stage6_analysis"
    out5 = run / "stage5_cells" / "truth"
    if not (out5 / "cell_summary.csv").exists():
        raise FileNotFoundError(f"{out5} has no isolated cells (run the pipeline with isolation.sources including truth)")
    px_um = cfg["optics"]["pixel_size_nm"] / 1000.0
    n_loci = int(cfg["loci"]["n_loci"])
    base = {"dilate_px": (cfg.get("locus_analysis") or {}).get("dilate_px", 2)}
    truth = truth_with_photons(run)
    cal_dir = s3 / "calibration"
    cal_opts = [False] + ([True] if calibrated and (cal_dir / "calibration.json").exists() else [])
    ej = s3 / "imaging_errors.json"
    zinfo = json.loads(ej.read_text()).get("z_stack") if ej.exists() else None
    (s6 / "comparison").mkdir(parents=True, exist_ok=True)
    rows = []
    jobs = [("2d", k, VARIANTS[k]) for k in (keys or VARIANTS) if k in VARIANTS]
    if three_d and zinfo is not None:
        jobs += [("3d", k, VARIANTS_3D[k]) for k in (keys or VARIANTS_3D) if k in VARIANTS_3D]
    for dim, key, (label, over) in jobs:
        for cal in cal_opts:
            t0 = time.time()
            kw = {**base, **over, "n_loci": n_loci}
            if cal:
                kw["calib"] = str(cal_dir)
            if dim == "2d":
                pos = map_cells(localize_loci, out5, n_workers=n_workers, px_um=px_um, **kw)
                err, sc = score_localization(pos, truth, px_um) if len(pos) else (pd.DataFrame(), {"n_detections": 0})
            else:
                from simlive.pipelines.tracking_demo import nominal_axial_sigma

                stacks = {f"locus{k}": str(s3 / f"locus{k}_zstack.tif") for k in range(n_loci)}
                pos = map_cells(localize_loci_3d, out5, n_workers=n_workers, stacks=stacks, z_um=zinfo["z_positions_um"],
                                axial_sigma_um=nominal_axial_sigma(cfg), **kw)
                err, sc = score_localization_3d(pos, truth, px_um) if len(pos) else (pd.DataFrame(), {"n_detections": 0})
            name = f"{key}{'_cal' if cal else ''}{'_3d' if dim == '3d' else ''}"
            if len(err):
                err.merge(pos, on=["cell_id", "t", "locus_id"], how="left").to_csv(s6 / "comparison" / f"{name}.csv", index=False)
            rows.append({"variant": key, "dim": dim, "calibrated": cal, "label": label, "settings": json.dumps(over),
                         **{c: sc.get(c) for c in SCORE_COLS}, "n_detections": sc.get("n_detections"), "seconds": round(time.time() - t0, 1)})
            progress(f"{name:>22}: detected {sc.get('detection_rate', 0):.3f}, median {1000 * (sc.get('median_error_um') or 0):.0f} nm, "
                     f"p95 {1000 * (sc.get('p95_error_um') or 0):.0f} nm, far {sc.get('frac_far') or 0:.3f}  ({time.time() - t0:.0f}s)")
    tab = pd.DataFrame(rows)
    tab.to_csv(s6 / "method_comparison.csv", index=False)
    return tab
