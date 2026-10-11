"""Stage 6 PLACEHOLDER (3D): locate every locus in the full z-stack of one isolated cell.

The 2D crops written by stage 5 (maximum-intensity projections) give each cell's window in every frame (`frames.csv`: crop_y0, crop_x0)
and its mask; this module re-reads the same window from the full z-stack TIFFs (memory-mapped, one frame at a time). Methods as in 2D
(`localize.py`): `dog_viterbi` (default; Chen et al. 2025: 3D median, 3D difference of Gaussians, Viterbi path with a 3D diffusion
potential, z distances scaled by plane spacing / pixel size, 3D Gaussian fit) or `brightest` (the original: brightest voxel after 3D
smoothing, refined by a centroid). Every detection gets the 3D ellipsoid
intensity (inner (7, 7, 5), shell to (19, 19, 7) px / planes). Positions: lab pixels (y, x) and height z_um (nominal slice positions).
`score_localization_3d` compares with the true 3D positions.

Axial median: a cosmic ray hits ONE plane, a real spot spans several planes when they are closer than the axial PSF width. With
median_z = 3 the median over neighbouring planes removes single-plane events but keeps finely sampled spots; with coarse planes (spacing >
axial sigma) it would remove the spot too. median_z = "auto" chooses 3 when the spacing is at most the axial sigma (`axial_sigma_um`, from the
microscope's optics), else 1.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter

from simlive.stage6_analysis import dog_viterbi as DV
from simlive.stage6_analysis import spots as SP

COLUMNS = ["cell_id", "t", "locus_id", "lab_y_px", "lab_x_px", "z_um", "snr", "n_candidates", "fit_ok", "intensity"]


def _crop3d(vol, y0: int, x0: int, S: int, pad_value: float) -> np.ndarray:
    nz, ny, nx = vol.shape
    out = np.full((nz, S, S), pad_value, dtype=np.float32)
    ya, yb, xa, xb = max(y0, 0), min(y0 + S, ny), max(x0, 0), min(x0 + S, nx)
    if yb > ya and xb > xa:
        out[:, ya - y0:yb - y0, xa - x0:xb - x0] = vol[:, ya:yb, xa:xb]
    return out


def _detect_one_3d(vol: np.ndarray, region2d: np.ndarray, sigma_xy: float, sigma_z: float, min_snr: float, win: int = 2):
    """Original placeholder: brightest spot inside the region (a 2D mask applied to every plane) -> (z_index, y, x, snr) or None."""
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


def median_z_for(median_z, step_um: float, axial_sigma_um: float | None) -> int:
    """Resolve median_z = "auto": 3 (axial median removes single-plane events) when the planes sample the axial PSF finely, else 1."""
    if median_z == "auto":
        return 3 if axial_sigma_um and step_um <= axial_sigma_um else 1
    return int(median_z)


def localize_loci_3d(cell, stacks: dict, z_um: list, n_loci: int = 2, method: str = "dog_viterbi", dilate_px: int = 2,
                     calib: str | None = None, axial_sigma_um: dict | None = None, px_um: float = 0.13, **kw) -> pd.DataFrame:
    """Locate every locus in every frame of one isolated cell from the full z-stacks. `stacks`: {channel: path of (T, Z, Y, X) TIFF}.
    With `calib` (folder of a saved Calibration) the crops are dark/flat/hot-pixel corrected and the measured chromatic shift (lateral and
    axial) of the colour is subtracted from the detected position. `axial_sigma_um`: {channel: nominal axial PSF sigma} for median_z: auto."""
    from simlive.stage6_analysis.localize import _load_calibration, centroid_in_crop, spot_params

    cal = _load_calibration(calib) if calib else None
    chans = cell.locus_channels[:n_loci]
    if not chans:
        return pd.DataFrame(columns=COLUMNS)
    reg = cell.analysis_mask(dilate_px)
    fr = cell.frames.reset_index(drop=True)
    S = cell.crop_size
    z = np.asarray(z_um, float)
    step = float(np.mean(np.diff(z))) if len(z) > 1 else 1.0
    p = spot_params(kw)
    rows = []
    for k, name in enumerate(chans):
        st = tifffile.memmap(stacks[name], mode="r")
        pad = float(np.percentile(np.asarray(st[0, len(z) // 2]), 1))

        def volume(i, r):
            v = _crop3d(np.asarray(st[int(r.t)]), int(r.crop_y0), int(r.crop_x0), S, pad)
            return cal.correct(v, name, int(r.crop_y0), int(r.crop_x0)).astype(np.float32) if cal is not None else v

        frames = list(fr.itertuples())
        inner, outer, ibg = p["intensity_inner_px"], p["intensity_outer_px"], p["intensity_background"]
        if method == "brightest":
            dets = []
            for i, r in enumerate(frames):
                v = volume(i, r)
                d = _detect_one_3d(v, reg[i], kw.get("smooth_sigma_px", 1.0), kw.get("smooth_sigma_z", 0.8), kw.get("min_snr", 5.0))
                dets.append(None if d is None else (d[0], d[1], d[2], d[3], 1, True, DV.ellipsoid_intensity(v, d[:3], inner, outer, ibg)[0]))
        elif method == "dog_viterbi":
            pz = {**p, "median_z": median_z_for(p["median_z"], step, (axial_sigma_um or {}).get(name)), "z_scale": step / px_um}
            refined = []
            for i, r in enumerate(frames):
                v = volume(i, r)
                c = DV.dog_candidates(v, reg[i], **pz) if reg[i].sum() >= 20 else None
                rr = [] if c is None else [SP.refine_peak(c, j, **pz) for j in range(len(c["peaks"]))]
                for q in rr:
                    q["intensity"] = DV.ellipsoid_intensity(v, q["pos"], inner, outer, ibg)[0]
                refined.append(rr)
            vol = float(np.mean(reg.sum(axis=(-2, -1)))) * len(z) * pz["z_scale"]      # where candidates can appear, z in pixel units
            choice, _ = DV.viterbi_track([np.array([q["pos"] for q in rr]) if rr else np.zeros((0, 3)) for rr in refined],
                                         [np.array([q["snr"] for q in rr]) for rr in refined], centroid_in_crop(fr, S),
                                         **{**pz, "search_volume": vol})
            dets = [None if j is None else (*rr[j]["pos"], rr[j]["snr"], rr[j]["n_candidates"], rr[j]["fit_ok"], rr[j]["intensity"])
                    for rr, j in zip(refined, choice)]
        else:
            raise ValueError(f"unknown locus_analysis.method {method!r}")
        for i, d in enumerate(dets):
            if d is None:
                continue
            r = frames[i]
            y, x, zz = r.crop_y0 + d[1], r.crop_x0 + d[2], z[0] + d[0] * step
            if cal is not None:
                sy, sx = cal.shift_px(name)
                y, x, zz = y - sy, x - sx, zz - cal.shift_z_um(name)
            rows.append((cell.cell_id, int(r.t), k, y, x, zz, d[3], d[4], bool(d[5]), d[6]))
    return pd.DataFrame(rows, columns=COLUMNS)


def locus_distances_3d(pos: pd.DataFrame, px_um: float) -> pd.DataFrame:
    """3D distance between locus 0 and locus 1 per cell and frame (frames where both were found), in um."""
    a = pos[pos.locus_id == 0].set_index(["cell_id", "t"])
    b = pos[pos.locus_id == 1].set_index(["cell_id", "t"])
    j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
    d = np.sqrt((px_um * (j.lab_y_px_a - j.lab_y_px_b)) ** 2 + (px_um * (j.lab_x_px_a - j.lab_x_px_b)) ** 2 + (j.z_um_a - j.z_um_b) ** 2)
    return d.rename("dist_um").reset_index()


def score_localization_3d(pos: pd.DataFrame, truth: pd.DataFrame, px_um: float, far_um: float = 0.5) -> tuple[pd.DataFrame, dict]:
    """Compare 3D detections with the true positions (same cell, frame, locus; truth in image coordinates, i.e. drift-corrected).
    err_um = 3D error, err_xy_um = lateral error, err_z_um = signed axial error. Rates and far fraction as in `score_localization`."""
    from simlive.stage6_analysis.localize import _visible, intensity_scores

    est = pos.assign(est_y_um=pos.lab_y_px * px_um, est_x_um=pos.lab_x_px * px_um)
    tr = truth.assign(visible=_visible(truth))
    j = est.merge(tr, on=["cell_id", "t", "locus_id"], how="inner", suffixes=("", "_true"))
    j["err_xy_um"] = np.hypot(j.est_y_um - j.y_um, j.est_x_um - j.x_um)
    j["err_z_um"] = j.z_um - j.z_um_true
    j["err_um"] = np.sqrt(j.err_xy_um**2 + j.err_z_um**2)
    err = j[["cell_id", "t", "locus_id", "err_um", "err_xy_um", "err_z_um", "visible"]]
    if err.empty:
        return err, {"n_detections": 0}
    ok = err[err.visible]
    n_vis = int(tr.visible.sum())
    near = ok.err_um <= far_um
    return err, {"n_detections": len(err), "detection_rate": len(ok) / max(n_vis, 1), "n_visible": n_vis,
                 "n_false_on_dark": int((~err.visible).sum()),
                 "rms_error_um": float(np.sqrt((ok.err_um ** 2).mean())), "median_error_um": float(ok.err_um.median()),
                 "p95_error_um": float(ok.err_um.quantile(0.95)), "rms_xy_um": float(np.sqrt((ok.err_xy_um ** 2).mean())),
                 "rms_z_um": float(np.sqrt((ok.err_z_um ** 2).mean())), "median_z_um": float(ok.err_z_um.median()),
                 "frac_far": float((~near).mean()), "rms_error_near_um": float(np.sqrt((ok.err_um[near] ** 2).mean())) if near.any() else np.nan,
                 **intensity_scores(j, j.err_um <= far_um)}
