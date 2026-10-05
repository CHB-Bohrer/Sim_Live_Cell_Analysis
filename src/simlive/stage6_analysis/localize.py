"""Stage 6 PLACEHOLDER: a simple per-cell locus localizer, to be replaced by the user's own analysis.

`localize_loci(cell)` takes one isolated cell (CellData). Every locus has its own colour channel (locus0, locus1, ...),
so in each frame the brightest spot inside the (slightly dilated) cell mask is that locus; its sub-pixel position is an
intensity-weighted centroid. Positions are returned in three frames: crop pixels, lab pixels, and the cell's own frame
(origin at the nucleus centroid, u along the nucleus' long axis, v perpendicular, micrometres).
`score_localization` compares against the simulation's true positions (valid for true cell IDs only).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


def _detect_one(img: np.ndarray, region: np.ndarray, sigma: float, min_snr: float, win: int = 2):
    """Brightest sub-pixel spot inside region -> (y, x, snr) or None."""
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


def localize_loci(cell, n_loci: int = 2, smooth_sigma_px: float = 1.0, min_snr: float = 5.0, dilate_px: int = 2,
                  px_um: float = 0.13, **_) -> pd.DataFrame:
    """Locate every locus (one colour channel each) in every frame of one isolated cell. One row per detection."""
    chans = cell.locus_channels[:n_loci]
    if not chans:
        return pd.DataFrame()
    reg = cell.analysis_mask(dilate_px)
    fr = cell.frames.reset_index(drop=True)
    S = cell.crop_size
    # nucleus long-axis angle (from mask moments) is only defined modulo 180 deg: unwrap so (u, v) is continuous in time
    theta = 0.5 * np.unwrap(2 * fr.theta_rad.to_numpy())
    rows = []
    for k, name in enumerate(chans):
        stack = cell.channels[name]
        for i, r in enumerate(fr.itertuples()):
            d = _detect_one(stack[i].astype(float), reg[i], smooth_sigma_px, min_snr)
            if d is None:
                continue
            c, s = np.cos(theta[i]), np.sin(theta[i])
            dy, dx = d[0] - S // 2, d[1] - S // 2                 # from the crop centre = nucleus centroid
            rows.append((cell.cell_id, r.t, k, d[0], d[1], r.crop_y0 + d[0], r.crop_x0 + d[1],
                         (dx * c + dy * s) * px_um, (dy * c - dx * s) * px_um, d[2]))
    return pd.DataFrame(rows, columns=["cell_id", "t", "locus_id", "crop_y_px", "crop_x_px", "lab_y_px", "lab_x_px",
                                       "u_um", "v_um", "snr"])


def locus_distances(pos: pd.DataFrame, px_um: float) -> pd.DataFrame:
    """Distance between locus 0 and locus 1 per cell and frame (frames where both were found), in um."""
    a = pos[pos.locus_id == 0].set_index(["cell_id", "t"])
    b = pos[pos.locus_id == 1].set_index(["cell_id", "t"])
    j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
    d = px_um * np.hypot(j.lab_y_px_a - j.lab_y_px_b, j.lab_x_px_a - j.lab_x_px_b)
    return d.rename("dist_um").reset_index()


def score_localization(pos: pd.DataFrame, truth: pd.DataFrame, px_um: float) -> tuple[pd.DataFrame, dict]:
    """Compare detections with the true locus positions (same cell, frame and locus). Valid for true cell IDs only."""
    est = pos.assign(est_y_um=pos.lab_y_px * px_um, est_x_um=pos.lab_x_px * px_um)
    j = est.merge(truth, on=["cell_id", "t", "locus_id"], how="inner")
    j["err_um"] = np.hypot(j.est_y_um - j.y_um, j.est_x_um - j.x_um)
    err = j[["cell_id", "t", "locus_id", "err_um"]]
    if err.empty:
        return err, {"n_detections": 0}
    return err, {"n_detections": len(err), "detection_rate": len(err) / len(truth),
                 "rms_error_um": float(np.sqrt((err.err_um ** 2).mean())), "median_error_um": float(err.err_um.median()),
                 "p95_error_um": float(err.err_um.quantile(0.95))}
