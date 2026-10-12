"""Stage 6 PLACEHOLDER: a per-cell locus localizer, to be replaced by the user's own analysis.

`localize_loci(cell)` takes one isolated cell (CellData). Every locus has its own colour channel (locus0, locus1, ...), so in each frame
the locus is the one spot of that colour inside the (slightly dilated) cell mask. Two methods (`method`, config `locus_analysis.method`):

- `dog_viterbi` (default): the spot pipeline of Chen, Levo, Zoller & Gregor 2025 (Nat. Struct. Mol. Biol.) for MS2 spots: median filter ->
  difference of Gaussians with a low threshold (several candidates) -> a Viterbi path through the candidates over the whole movie with a
  diffusion potential estimated from the path itself (at most one spot per cell; frames may have none). See `dog_viterbi.py`.
- `brightest` (the original placeholder): smooth, take the brightest pixel if its SNR passes `min_snr`, refine by a 5 x 5 centroid.
Every detection also gets the ellipsoid intensity of Chen et al. 2025 (`intensity`: mean in a (7, 7) px ellipse at the spot minus the mean
of the (19, 19) px ellipse around it), so spot brightness can be compared with the truth.

Positions are returned in three frames: crop pixels, lab pixels, and the cell's own frame (origin at the nucleus centroid, u along the
nucleus' long axis, v perpendicular, micrometres). `score_localization` compares against the simulation's true positions (true cell IDs only).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

from simlive.stage6_analysis import dog_viterbi as DV
from simlive.stage6_analysis import spots as SP

METHODS = ("dog_viterbi", "brightest")
COLUMNS = ["cell_id", "t", "locus_id", "crop_y_px", "crop_x_px", "lab_y_px", "lab_x_px", "u_um", "v_um", "snr", "n_candidates", "fit_ok",
           "intensity"]


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
    """The spot-finder settings from a locus_analysis config block (defaults for anything missing)."""
    return {k: kw.get(k, v) for k, v in SP.DEFAULTS.items()}


def analyse_channel(stack, region, centre: np.ndarray | None = None, method: str = "dog_viterbi", **kw) -> dict:
    """dog_viterbi on every frame of one colour of one cell: stack (n, S, S), region (n, S, S) bool, centre (n, 2) the nucleus centroid in crop
    pixels (for the Viterbi potential). Returns per frame: the candidate threshold, the candidates (peaks brightest first, DoG values, SNR),
    the chosen candidate (None = no locus), its refined fit and the ellipsoid intensity at the fitted position, plus the diffusion potential D
    (px^2 per frame)."""
    if method != "dog_viterbi":
        raise ValueError(f"unknown locus_analysis.method {method!r}; use one of {METHODS}")
    p = spot_params(kw)
    n = len(stack)
    centre = np.full((n, 2), stack.shape[-1] // 2, float) if centre is None else np.asarray(centre, float)
    cands = [DV.dog_candidates(stack[i], region[i], steps=True, **p) if region[i].sum() >= 20 else None for i in range(n)]
    thr = [c["thr"] if c is not None else np.nan for c in cands]
    choice, D = DV.viterbi_track([c["peaks"] if c is not None else np.zeros((0, 2)) for c in cands],
                                 [c["snr"] if c is not None else np.zeros(0) for c in cands], centre,
                                 **{**p, "search_volume": float(np.mean(region.sum(axis=(-2, -1))))})
    fits = [SP.refine_peak(c, j, **p) if c is not None and j is not None else None for c, j in zip(cands, choice)]
    for i, f in enumerate(fits):
        if f is not None:
            f["intensity"] = DV.ellipsoid_intensity(stack[i], f["pos"], p["intensity_inner_px"], p["intensity_outer_px"], p["intensity_background"])[0]
    light = [None if c is None else {k: c[k] for k in ("peaks", "values", "snr")} for c in cands]
    return {"params": p, "method": method, "thr_frame": thr, "thr": thr, "candidates": light, "choice": choice, "fits": fits, "D_px2": D}


def detect_channel(stack, region, method: str = "dog_viterbi", centre: np.ndarray | None = None, **kw) -> list:
    """Run one method on every frame of one colour of one cell (see `analyse_channel`). Returns one entry per frame: None or
    (y, x, snr, n_candidates, fit_ok, intensity) in crop pixels."""
    if method == "brightest":
        p = spot_params(kw)
        out = []
        for i in range(len(stack)):
            d = _detect_one(np.asarray(stack[i], float), region[i], kw.get("smooth_sigma_px", 1.0), kw.get("min_snr", 5.0))
            out.append(None if d is None else (d[0], d[1], d[2], 1, True, DV.ellipsoid_intensity(
                stack[i], (d[0], d[1]), p["intensity_inner_px"], p["intensity_outer_px"], p["intensity_background"])[0]))
        return out
    return [None if f is None else (f["pos"][-2], f["pos"][-1], f["snr"], f["n_candidates"], f["fit_ok"], f["intensity"])
            for f in analyse_channel(stack, region, centre, method, **kw)["fits"]]


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


def localize_loci(cell, n_loci: int = 2, method: str = "dog_viterbi", dilate_px: int = 2, px_um: float = 0.13, calib: str | None = None,
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
                         d[2], d[3], bool(d[4]), d[5]))
    return pd.DataFrame(rows, columns=COLUMNS)


def locus_distances(pos: pd.DataFrame, px_um: float) -> pd.DataFrame:
    """Distance between locus 0 and locus 1 per cell and frame (frames where both were found), in um."""
    a = pos[pos.locus_id == 0].set_index(["cell_id", "t"])
    b = pos[pos.locus_id == 1].set_index(["cell_id", "t"])
    j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
    d = px_um * np.hypot(j.lab_y_px_a - j.lab_y_px_b, j.lab_x_px_a - j.lab_x_px_b)
    return d.rename("dist_um").reset_index()


def _visible(truth: pd.DataFrame) -> pd.Series:
    """True where the simulation drew a spot at all (the probe model can leave a locus with no probe attached during an exposure) and
    it lies inside the image (cells may walk in and out of the field: geometry.outside_um)."""
    vis = truth.photons > 0 if "photons" in truth.columns else pd.Series(True, index=truth.index)
    if "in_image" in truth.columns:
        vis = vis & truth.in_image.astype(bool)
    return vis


def photons_in_frame(truth: pd.DataFrame, imaging_truth: pd.DataFrame | None) -> pd.DataFrame:
    """Add `photons_frame` = the spot's photons in that frame (photons x that colour's bleaching factor), the quantity a spot intensity
    should be proportional to (up to flat-field and flicker)."""
    if "photons" not in truth.columns:
        return truth
    out = truth.copy()
    out["photons_frame"] = out.photons.astype(float)
    if "bleach" in out.columns:            # per-cell bleaching (each cell by its own light dose): the factor of this spot's cell and colour
        out["photons_frame"] = out.photons * out.bleach.fillna(1.0)
    elif imaging_truth is not None:
        it = imaging_truth.set_index("t")
        for k in out.locus_id.unique():
            col = f"bleach_locus{int(k)}"
            if col in it.columns:
                m = out.locus_id == k
                out.loc[m, "photons_frame"] = out.loc[m, "photons"].to_numpy() * it[col].reindex(out.loc[m, "t"]).to_numpy()
    return out


def intensity_scores(j: pd.DataFrame, near: pd.Series) -> dict:
    """How well the measured spot intensity follows the true photons (correct detections only): Pearson r, per locus colour."""
    out = {}
    if "intensity" not in j.columns or "photons_frame" not in j.columns:
        return out
    for k, g in j[near & j.visible].groupby("locus_id"):
        g = g[np.isfinite(g.intensity) & np.isfinite(g.photons_frame)]
        if len(g) > 5 and g.photons_frame.std() > 0 and g.intensity.std() > 0:
            out[f"intensity_r_locus{int(k)}"] = float(np.corrcoef(g.intensity, g.photons_frame)[0, 1])
            out[f"counts_per_photon_locus{int(k)}"] = float(np.median(g.intensity / g.photons_frame))
    return out


def score_localization(pos: pd.DataFrame, truth: pd.DataFrame, px_um: float, far_um: float = 0.5) -> tuple[pd.DataFrame, dict]:
    """Compare detections with the true locus positions (same cell, frame and locus). Valid for true cell IDs only.
    detection_rate = detections / visible loci (a locus with no light cannot be found); `n_false_on_dark` = detections where no spot was
    drawn; `frac_far` = share of detections more than `far_um` from the truth (a wrong object was picked: hot pixel, cosmic ray, ghost);
    `intensity_r_locusK` = correlation of the measured ellipsoid intensity with the true photons in that frame (if both are known)."""
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
                 "rms_error_near_um": float(np.sqrt((ok.err_um[ok.err_um <= far_um] ** 2).mean())) if (ok.err_um <= far_um).any() else np.nan,
                 **intensity_scores(j, j.err_um <= far_um)}
