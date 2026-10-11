"""Stage 6 PLACEHOLDER (3D): locate every locus in the full z-stack of one isolated cell.

The 2D crops written by stage 5 (maximum-intensity projections) give each cell's window in every frame (`frames.csv`: crop_y0, crop_x0)
and its mask; this module re-reads the same window from the full z-stack TIFFs (memory-mapped, one frame at a time) and finds, per locus
colour, the brightest voxel inside the (dilated) cell mask after smoothing in 3D, then refines it with an intensity-weighted centroid.
Positions: lab pixels (y, x) and height z_um (from the slice positions). `score_localization_3d` compares with the true 3D positions.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter


def _crop3d(vol, y0: int, x0: int, S: int, pad_value: float) -> np.ndarray:
    nz, ny, nx = vol.shape
    out = np.full((nz, S, S), pad_value, dtype=float)
    ya, yb, xa, xb = max(y0, 0), min(y0 + S, ny), max(x0, 0), min(x0 + S, nx)
    if yb > ya and xb > xa:
        out[:, ya - y0:yb - y0, xa - x0:xb - x0] = vol[:, ya:yb, xa:xb]
    return out


def _detect_one_3d(vol: np.ndarray, region2d: np.ndarray, sigma_xy: float, sigma_z: float, min_snr: float, win: int = 2):
    """Brightest spot inside the region (a 2D mask applied to every plane) -> (z_index, y, x, snr) or None."""
    region = np.broadcast_to(region2d[None], vol.shape)
    vals = vol[region]
    if vals.size < 20:
        return None
    bg = np.median(vals)
    sm = gaussian_filter(vol - bg, (sigma_z, sigma_xy, sigma_xy))
    noise = 1.4826 * np.median(np.abs(sm[region] - np.median(sm[region]))) + 1e-9
    work = np.where(region, sm, -np.inf)
    iz, iy, ix = np.unravel_index(np.argmax(work), work.shape)
    if work[iz, iy, ix] / noise < min_snr:
        return None
    z0, z1 = max(iz - 1, 0), iz + 2
    y0, y1, x0, x1 = max(iy - win, 0), iy + win + 1, max(ix - win, 0), ix + win + 1
    w = np.clip(vol[z0:z1, y0:y1, x0:x1] - bg, 0, None)
    if w.sum() <= 0:
        return None
    zz, yy, xx = np.indices(w.shape)
    return ((w * (zz + z0)).sum() / w.sum(), (w * (yy + y0)).sum() / w.sum(), (w * (xx + x0)).sum() / w.sum(), work[iz, iy, ix] / noise)


def localize_loci_3d(cell, stacks: dict, z_um: list, n_loci: int = 2, smooth_sigma_px: float = 1.0, smooth_sigma_z: float = 0.8,
                     min_snr: float = 5.0, dilate_px: int = 2, calib: str | None = None, **_) -> pd.DataFrame:
    """Locate every locus in every frame of one isolated cell from the full z-stacks. `stacks`: {channel: path of (T, Z, Y, X) TIFF}.
    With `calib` (folder of a saved Calibration) the crops are dark/flat/hot-pixel corrected and the measured chromatic shift (lateral and
    axial) of the colour is subtracted from the detected position."""
    from simlive.stage6_analysis.localize import _load_calibration

    cal = _load_calibration(calib) if calib else None
    chans = cell.locus_channels[:n_loci]
    if not chans:
        return pd.DataFrame()
    reg = cell.analysis_mask(dilate_px)
    fr = cell.frames.reset_index(drop=True)
    S = cell.crop_size
    z = np.asarray(z_um, float)
    step = float(np.mean(np.diff(z))) if len(z) > 1 else 1.0
    rows = []
    for k, name in enumerate(chans):
        st = tifffile.memmap(stacks[name], mode="r")
        pad = float(np.percentile(np.asarray(st[0, len(z) // 2]), 1))
        for i, r in enumerate(fr.itertuples()):
            vol = _crop3d(np.asarray(st[int(r.t)]), int(r.crop_y0), int(r.crop_x0), S, pad)
            if cal is not None:
                vol = cal.correct(vol, name, int(r.crop_y0), int(r.crop_x0)).astype(float)
            d = _detect_one_3d(vol, reg[i], smooth_sigma_px, smooth_sigma_z, min_snr)
            if d is None:
                continue
            y, x, zz = r.crop_y0 + d[1], r.crop_x0 + d[2], z[0] + d[0] * step
            if cal is not None:
                sy, sx = cal.shift_px(name)
                y, x, zz = y - sy, x - sx, zz - cal.shift_z_um(name)
            rows.append((cell.cell_id, r.t, k, y, x, zz, d[3]))
    return pd.DataFrame(rows, columns=["cell_id", "t", "locus_id", "lab_y_px", "lab_x_px", "z_um", "snr"])


def locus_distances_3d(pos: pd.DataFrame, px_um: float) -> pd.DataFrame:
    """3D distance between locus 0 and locus 1 per cell and frame (frames where both were found), in um."""
    a = pos[pos.locus_id == 0].set_index(["cell_id", "t"])
    b = pos[pos.locus_id == 1].set_index(["cell_id", "t"])
    j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
    d = np.sqrt((px_um * (j.lab_y_px_a - j.lab_y_px_b)) ** 2 + (px_um * (j.lab_x_px_a - j.lab_x_px_b)) ** 2 + (j.z_um_a - j.z_um_b) ** 2)
    return d.rename("dist_um").reset_index()


def score_localization_3d(pos: pd.DataFrame, truth: pd.DataFrame, px_um: float) -> tuple[pd.DataFrame, dict]:
    """Compare 3D detections with the true positions (same cell, frame, locus; truth in image coordinates, i.e. drift-corrected).
    err_um = 3D error, err_xy_um = lateral error, err_z_um = signed axial error."""
    est = pos.assign(est_y_um=pos.lab_y_px * px_um, est_x_um=pos.lab_x_px * px_um)
    j = est.merge(truth, on=["cell_id", "t", "locus_id"], how="inner", suffixes=("", "_true"))
    j["err_xy_um"] = np.hypot(j.est_y_um - j.y_um, j.est_x_um - j.x_um)
    j["err_z_um"] = j.z_um - j.z_um_true
    j["err_um"] = np.sqrt(j.err_xy_um**2 + j.err_z_um**2)
    err = j[["cell_id", "t", "locus_id", "err_um", "err_xy_um", "err_z_um"]]
    if err.empty:
        return err, {"n_detections": 0}
    return err, {"n_detections": len(err), "detection_rate": len(err) / len(truth),
                 "rms_error_um": float(np.sqrt((err.err_um ** 2).mean())), "median_error_um": float(err.err_um.median()),
                 "p95_error_um": float(err.err_um.quantile(0.95)), "rms_xy_um": float(np.sqrt((err.err_xy_um ** 2).mean())),
                 "rms_z_um": float(np.sqrt((err.err_z_um ** 2).mean())), "median_z_um": float(err.err_z_um.median())}
