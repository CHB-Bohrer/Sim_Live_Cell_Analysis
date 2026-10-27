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
from simlive.stage6_analysis.localize import localize_loci, photons_in_frame, score_localization
from simlive.stage6_analysis.localize3d import localize_loci_3d, score_localization_3d

# key -> (label, settings that differ from spots.DEFAULTS). Every variant starts from the defaults (not from the run's own locus_analysis block,
# except dilate_px), so the table means the same thing for every run. "default" = dog_viterbi (Chen et al. 2025) with a Gaussian fit and a
# median shell background; the "..." variants change one thing of it.
VARIANTS = {
    "original": ("Original placeholder: brightest pixel after smoothing + 5x5 centroid", {"method": "brightest"}),
    "default": ("Default: Chen et al. 2025 (median 3x3, difference of Gaussians at 3 x noise, Viterbi path with a diffusion potential) "
                "+ Gaussian fit, median shell background for the intensity", {"method": "dog_viterbi"}),
    "chen_literal": ("Chen et al. 2025 as described: spot centroid, mean shell background",
                     {"method": "dog_viterbi", "refine": "centroid", "intensity_background": "mean"}),
    "no_median": ("... without the median filter", {"method": "dog_viterbi", "median_px": 1}),
    "dog_snr5": ("... candidate threshold 5 x noise instead of 3 (fewer candidates)", {"method": "dog_viterbi", "dog_snr": 5.0}),
    "integer": ("... integer pixel position, no sub-pixel step", {"method": "dog_viterbi", "refine": "none"}),
}
VARIANTS_3D = {
    "original": ("Original placeholder (3D): brightest voxel + centroid", {"method": "brightest"}),
    "default": ("Default in 3D: Chen et al. 2025 (3D median auto, 3D DoG, Viterbi with a 3D diffusion potential) + 3D Gaussian fit",
                {"method": "dog_viterbi"}),
    "chen_literal": ("Chen et al. 2025 in 3D as described: spot centroid, mean shell background",
                     {"method": "dog_viterbi", "refine": "centroid", "intensity_background": "mean"}),
    "zmed3": ("... axial median over 3 planes always", {"method": "dog_viterbi", "median_z": 3}),
}
SCORE_COLS = ["detection_rate", "n_false_on_dark", "median_error_um", "p95_error_um", "rms_error_um", "frac_far", "rms_error_near_um",
              "intensity_r_locus0", "intensity_r_locus1"]


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
        truth = truth.merge(ph.drop(columns=[c for c in ("copy",) if c not in truth.columns]), on=["t", "cell_id", "locus_id"]
                            + (["copy"] if "copy" in truth.columns else []), how="left")
    it = s3 / "imaging_truth.csv"
    return photons_in_frame(truth, pd.read_csv(it) if it.exists() else None)


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
    # every variant starts from the run's own locus_analysis settings (copies, pair prior, candidates ...), so "default" is what the pipeline ran;
    # until 2026-10-26 only dilate_px was taken, and on multi-copy runs every row was the one-copy finder
    base = {"dilate_px": 2, **{k: v for k, v in (cfg.get("locus_analysis") or {}).items() if k != "method"}}
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
                                axial_sigma_um=nominal_axial_sigma(cfg), px_um=px_um, **kw)
                err, sc = score_localization_3d(pos, truth, px_um) if len(pos) else (pd.DataFrame(), {"n_detections": 0})
            name = f"{key}{'_cal' if cal else ''}{'_3d' if dim == '3d' else ''}"
            if len(err):
                from simlive.stage6_analysis.localize import match_copies

                err.merge(match_copies(pos, truth, px_um), on=["cell_id", "t", "locus_id"] + (["copy"] if "copy" in err.columns else []),
                          how="left").to_csv(s6 / "comparison" / f"{name}.csv", index=False)
            rows.append({"variant": key, "dim": dim, "calibrated": cal, "label": label, "settings": json.dumps(over),
                         **{c: sc.get(c) for c in SCORE_COLS}, "n_detections": sc.get("n_detections"), "seconds": round(time.time() - t0, 1)})
            progress(f"{name:>22}: detected {sc.get('detection_rate', 0):.3f}, median {1000 * (sc.get('median_error_um') or 0):.0f} nm, "
                     f"p95 {1000 * (sc.get('p95_error_um') or 0):.0f} nm, far {sc.get('frac_far') or 0:.3f}  ({time.time() - t0:.0f}s)")
    tab = pd.DataFrame(rows)
    out = s6 / "method_comparison.csv"
    if out.exists() and len(tab):           # keep the rows of variants not re-run this time
        old = pd.read_csv(out)
        key = set(zip(tab.variant, tab.dim, tab.calibrated))
        old = old[[k not in key for k in zip(old.variant, old.dim, old.calibrated)]]
        order = {k: i for i, k in enumerate(list(VARIANTS) + [k for k in VARIANTS_3D if k not in VARIANTS])}
        tab = pd.concat([old, tab], ignore_index=True)
        tab = tab.assign(_o=tab.variant.map(order).fillna(99)).sort_values(["dim", "_o", "calibrated"]).drop(columns="_o")
    tab.to_csv(out, index=False)
    return tab
