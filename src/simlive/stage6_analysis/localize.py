"""Stage 6 PLACEHOLDER: a per-cell locus localizer, to be replaced by the user's own analysis.

`localize_loci(cell)` takes one isolated cell (CellData). Every locus has its own colour channel (locus0, locus1, ...), so in each frame
the locus is the one spot of that colour inside the (slightly dilated) cell mask. Two methods (`method`, config `locus_analysis.method`):

- `median_peaks` (default): the detection steps of Bennabi et al. 2025 (Sci. Adv. 11, eadv7790), rescaled to diffraction-limited spots:
  median filter -> one threshold per movie (median of the per-frame thresholds) -> area above it, small objects removed -> Gaussian
  smoothing -> local maxima (min distance) = candidates; the brightest is the locus, refined by a 2D Gaussian fit. See `spots.py`.
- `brightest` (the original placeholder): smooth, take the brightest pixel if its SNR passes `min_snr`, refine by a 5 x 5 centroid.

Positions are returned in three frames: crop pixels, lab pixels, and the cell's own frame (origin at the nucleus centroid, u along the
nucleus' long axis, v perpendicular, micrometres). `score_localization` compares against the simulation's true positions (true cell IDs only).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

from simlive.stage6_analysis import spots as SP

METHODS = ("median_peaks", "brightest")
COLUMNS = ["cell_id", "t", "locus_id", "crop_y_px", "crop_x_px", "lab_y_px", "lab_x_px", "u_um", "v_um", "snr", "n_candidates", "fit_ok"]


def _detect_one(img: np.ndarray, region: np.ndarray, sigma: float, min_snr: float, win: int = 2):
    """Original placeholder: brightest sub-pixel spot inside region -> (y, x, snr) or None."""
    vals = img[region]
    if vals.size < 20:
        return None
    bg = np.median(vals)
    sm = gaussian_filter(img - bg, sigma)
    noise = 1.4826 * np.median(np.abs(sm[region] - np.median(sm[region]))) + 1e-9
    work = np.where(region, sm, -np.inf)
    iy, ix = np.unravel_index(np.argmax(work), work.shape)
    if work[iy, ix] / noise < min_snr:
        return None
    y0, y1, x0, x1 = max(iy - win, 0), iy + win + 1, max(ix - win, 0), ix + win + 1
    w = np.clip(img[y0:y1, x0:x1] - bg, 0, None)
    if w.sum() <= 0:
        return None
    yy, xx = np.indices(img.shape)
    return (w * yy[y0:y1, x0:x1]).sum() / w.sum(), (w * xx[y0:y1, x0:x1]).sum() / w.sum(), work[iy, ix] / noise


_CAL = {}


def _load_calibration(path):
    if path not in _CAL:
        from simlive.stage3_microscopy.calibration.model import Calibration

        _CAL[path] = Calibration.load(path)
    return _CAL[path]


def spot_params(kw: dict) -> dict:
    """The median_peaks settings from a locus_analysis config block (defaults for anything missing)."""
    return {k: kw.get(k, v) for k, v in SP.DEFAULTS.items()}


def analyse_channel(stack, region, centre: np.ndarray | None = None, **kw) -> dict:
    """median_peaks on every frame of one colour of one cell: stack (n, S, S), region (n, S, S) bool, centre (n, 2) the nucleus centroid in
    crop pixels (used by select: track). Returns per frame: the frame's own threshold, the threshold used (one per movie by default), the
    candidates (peaks brightest first, smoothed values, SNR), the chosen candidate (None = no locus) and its refined fit."""
    p = spot_params(kw)
    n = len(stack)
    pres = [SP.prefilter(stack[i], region[i], **p) if region[i].sum() >= 20 else None for i in range(n)]
    thr_f = [SP.frame_threshold(q, region[i], p["threshold"], p["min_snr"], p["threshold_on"]) if q is not None else np.nan for i, q in enumerate(pres)]
    thr = SP.movie_threshold(thr_f, p["threshold_scope"])
    cands = [SP.candidates(stack[i], region[i], thr[i], pres[i], **p) if pres[i] is not None else None for i in range(n)]
    del pres
    choice = _choose(cands, centre, p, S=stack.shape[-1])
    fits = [SP.refine_peak(c, j, **p) if c is not None and j is not None else None for c, j in zip(cands, choice)]
    light = [None if c is None else {k: c[k] for k in ("peaks", "values", "snr")} for c in cands]
    return {"params": p, "thr_frame": thr_f, "thr": thr, "candidates": light, "choice": choice, "fits": fits}


def detect_channel(stack, region, method: str = "median_peaks", centre: np.ndarray | None = None, **kw) -> list:
    """Run one method on every frame of one colour of one cell (see `analyse_channel`). Returns one entry per frame: None or
    (y, x, snr, n_candidates, fit_ok) in crop pixels."""
    if method == "brightest":
        out = []
        for i in range(len(stack)):
            d = _detect_one(np.asarray(stack[i], float), region[i], kw.get("smooth_sigma_px", 1.0), kw.get("min_snr", 5.0))
            out.append(None if d is None else (d[0], d[1], d[2], 1, True))
        return out
    if method != "median_peaks":
        raise ValueError(f"unknown locus_analysis.method {method!r}; use one of {METHODS}")
    return [None if f is None else (f["pos"][-2], f["pos"][-1], f["snr"], f["n_candidates"], f["fit_ok"])
            for f in analyse_channel(stack, region, centre, **kw)["fits"]]


def _choose(cands: list, centre, p: dict, S: int) -> list:
    """Which candidate is the locus in each frame: the brightest, or (select: track) the one consistent with the neighbouring frames."""
    if p["select"] == "track":
        c0 = np.full((len(cands), 2), S // 2, float) if centre is None else np.asarray(centre, float)
        return SP.select_track([c["peaks"] if c is not None else np.zeros((0, 2)) for c in cands], c0, int(p["track_window"]),
                               float(p["track_gate_px"]))
    return [0 if c is not None and len(c["peaks"]) else None for c in cands]


def channel_inputs(cell, name: str, dilate_px: int = 2, calib: str | None = None) -> dict:
    """What the localizer works on for one colour of one cell: the crops (calibration-corrected if `calib`), the analysis region, the nucleus
    centroid in crop pixels (the crop is centred on the ROUNDED centroid, so up to 0.5 px off the crop centre) and the chromatic shift to remove."""
    cal = _load_calibration(calib) if calib else None
    fr = cell.frames.reset_index(drop=True)
    stack = cell.channels[name].astype(np.float32)
    if cal is not None:
        stack = np.stack([cal.correct(stack[i], name, int(r.crop_y0), int(r.crop_x0)) for i, r in enumerate(fr.itertuples())])
    return {"stack": stack, "region": cell.analysis_mask(dilate_px), "frames": fr, "centre": centroid_in_crop(fr, cell.crop_size),
            "shift": cal.shift_px(name) if cal is not None else (0.0, 0.0)}


def centroid_in_crop(fr: pd.DataFrame, S: int) -> np.ndarray:
    """(n, 2) nucleus centroid in crop pixels from stage 5's frames.csv (the crop centre if the table has no centroid columns)."""
    if {"cy_px", "cx_px"} <= set(fr.columns):
        return np.c_[fr.cy_px.to_numpy() - fr.crop_y0.to_numpy(), fr.cx_px.to_numpy() - fr.crop_x0.to_numpy()]
    return np.full((len(fr), 2), S // 2, float)


def localize_loci(cell, n_loci: int = 2, method: str = "median_peaks", dilate_px: int = 2, px_um: float = 0.13, calib: str | None = None,
                  **kw) -> pd.DataFrame:
    """Locate every locus (one colour channel each) in every frame of one isolated cell. One row per detection.
    With `calib` (folder of a saved Calibration) every crop is first dark-subtracted, flat-fielded and hot-pixel-corrected, and the measured
    chromatic shift of the colour is subtracted from the detected position. Other keywords: the method's settings (see `spots.DEFAULTS`)."""
    chans = cell.locus_channels[:n_loci]
    if not chans:
        return pd.DataFrame(columns=COLUMNS)
    fr = cell.frames.reset_index(drop=True)
    # nucleus long-axis angle (from mask moments) is only defined modulo 180 deg: unwrap so (u, v) is continuous in time
    theta = 0.5 * np.unwrap(2 * fr.theta_rad.to_numpy())
    rows = []
    for k, name in enumerate(chans):
        ci = channel_inputs(cell, name, dilate_px, calib)
        (sy, sx), cen = ci["shift"], ci["centre"]
        for i, d in enumerate(detect_channel(ci["stack"], ci["region"], method, cen, **kw)):
            if d is None:
                continue
            y, x = d[0] - sy, d[1] - sx
            c, s = np.cos(theta[i]), np.sin(theta[i])
            dy, dx = y - cen[i, 0], x - cen[i, 1]
            r = fr.iloc[i]
            rows.append((cell.cell_id, int(r.t), k, y, x, r.crop_y0 + y, r.crop_x0 + x, (dx * c + dy * s) * px_um, (dy * c - dx * s) * px_um,
                         d[2], d[3], bool(d[4])))
    return pd.DataFrame(rows, columns=COLUMNS)


def locus_distances(pos: pd.DataFrame, px_um: float) -> pd.DataFrame:
    """Distance between locus 0 and locus 1 per cell and frame (frames where both were found), in um."""
    a = pos[pos.locus_id == 0].set_index(["cell_id", "t"])
    b = pos[pos.locus_id == 1].set_index(["cell_id", "t"])
    j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
    d = px_um * np.hypot(j.lab_y_px_a - j.lab_y_px_b, j.lab_x_px_a - j.lab_x_px_b)
    return d.rename("dist_um").reset_index()


def _visible(truth: pd.DataFrame) -> pd.Series:
    """True where the simulation drew a spot at all (the probe model can leave a locus with no probe attached during an exposure)."""
    if "photons" in truth.columns:
        return truth.photons > 0
    return pd.Series(True, index=truth.index)


def score_localization(pos: pd.DataFrame, truth: pd.DataFrame, px_um: float, far_um: float = 0.5) -> tuple[pd.DataFrame, dict]:
    """Compare detections with the true locus positions (same cell, frame and locus). Valid for true cell IDs only.
    detection_rate = detections / visible loci (a locus with no light cannot be found); `n_false_on_dark` = detections where no spot was
    drawn; `frac_far` = share of detections more than `far_um` from the truth (a wrong object was picked: hot pixel, cosmic ray, ghost)."""
    est = pos.assign(est_y_um=pos.lab_y_px * px_um, est_x_um=pos.lab_x_px * px_um)
    tr = truth.assign(visible=_visible(truth))
    j = est.merge(tr, on=["cell_id", "t", "locus_id"], how="inner")
    j["err_um"] = np.hypot(j.est_y_um - j.y_um, j.est_x_um - j.x_um)
    err = j[["cell_id", "t", "locus_id", "err_um", "visible"]]
    if err.empty:
        return err, {"n_detections": 0}
    ok = err[err.visible]
    n_vis = int(tr.visible.sum())
    return err, {"n_detections": len(err), "detection_rate": len(ok) / max(n_vis, 1), "n_visible": n_vis,
                 "n_false_on_dark": int((~err.visible).sum()),
                 "rms_error_um": float(np.sqrt((ok.err_um ** 2).mean())) if len(ok) else np.nan,
                 "median_error_um": float(ok.err_um.median()) if len(ok) else np.nan,
                 "p95_error_um": float(ok.err_um.quantile(0.95)) if len(ok) else np.nan,
                 "frac_far": float((ok.err_um > far_um).mean()) if len(ok) else np.nan,
                 "rms_error_near_um": float(np.sqrt((ok.err_um[ok.err_um <= far_um] ** 2).mean())) if (ok.err_um <= far_um).any() else np.nan}
