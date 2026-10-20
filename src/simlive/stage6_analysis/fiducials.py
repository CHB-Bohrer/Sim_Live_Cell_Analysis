"""Stage 6 PLACEHOLDER: use the fiducial beads (and, optionally, the piezo's own position log) to correct locus positions in z-stacks.

Two different things can go wrong with the planes of a z-stack, and they need different cures:

  * the WHOLE picture is displaced in a frame (stage drift, focus error): every bead and every locus moves by the same amount.
    Beads on the slide measure it: track them frame by frame, subtract the median displacement. `track_fiducials`.
  * the planes are not where they should be relative to each other (piezo step errors, non-linearity, hysteresis, settling): the positions
    of planes differ from their nominal heights by different amounts, so two loci at different heights are displaced by different amounts.
    Beads that all sit on the slide see only the planes near the slide, so they cannot measure this. A piezo that logs its own position
    can: replace the nominal plane heights by the logged ones. `correct_with_plane_positions`.

`score_variants` compares the uncorrected positions with each correction, using the locus-locus distance (which a shift common to all
loci cannot change) and the height after removing the constant offset.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_dilation

from simlive.stage3_microscopy.calibration.fit import _fit_bead, detect_beads


def track_fiducials(stack_path, z_um=None, seg_path=None, dilate_px: int = 8, half: int = 6, min_fraction: float = 0.7, z_per_frame: dict | None = None):
    """Follow the fiducial beads through the movie. stack_path: (T, Z, Y, X) (or (T, Y, X)) TIFF of one colour; seg_path: a (T, Y, X) label
    stack whose nuclei are not beads (beads are searched for outside the nuclei in frame 0). Returns (shifts, fits):
    shifts = DataFrame(t, dy_px, dx_px, dz_um, n_beads) = median displacement of the beads relative to frame 0 (z in um, nominal plane heights),
    fits = (T, n_beads, 5) array of the individual fits (y, x, sigma_xy, z, sigma_z). z_per_frame {t: plane heights}: use these (e.g. the piezo's
    logged heights) instead of the nominal ones to turn a bead's profile along the planes into a height."""
    st = tifffile.memmap(stack_path, mode="r")
    three_d = st.ndim == 4
    T = st.shape[0]
    f0 = np.asarray(st[0])
    P0 = detect_beads(f0.max(axis=0) if three_d else f0)
    if seg_path is not None and len(P0):
        seg = tifffile.memmap(seg_path, mode="r")
        m = binary_dilation(np.asarray(seg[0]) > 0, iterations=dilate_px)
        P0 = P0[[not m[int(round(y)), int(round(x))] for y, x in P0]]
    n = len(P0)
    res = np.full((T, n, 5), np.nan)
    disp = np.zeros(2)
    for t in range(T):
        fr = np.asarray(st[t])
        for i in range(n):
            f = _fit_bead(fr, P0[i, 0] + disp[0], P0[i, 1] + disp[1], z_per_frame[t] if z_per_frame is not None else z_um, half)
            if f is not None:
                res[t, i] = f
        d = np.nanmedian(res[t, :, :2] - res[0, :, :2], axis=0) if t > 0 and np.isfinite(res[t, :, 0]).any() else np.zeros(2)
        disp = np.where(np.isfinite(d), d, disp)
    good = np.isfinite(res[:, :, 0]).mean(axis=0) >= min_fraction
    res = res[:, good]
    rel = res - res[:1]
    keep = consistent_beads(rel, three_d)
    res, rel = res[:, keep], rel[:, keep]
    cols = {"dy_px": 0, "dx_px": 1} | ({"dz_um": 3} if three_d else {})
    shifts = pd.DataFrame({"t": np.arange(T), "n_beads": int(res.shape[1])})
    for name, c in cols.items():
        with warnings.catch_warnings():             # frames in which no kept bead was fitted give NaN here (interpolated below)
            warnings.simplefilter("ignore", RuntimeWarning)
            v = np.nanmedian(rel[:, :, c], axis=1) if res.shape[1] else np.zeros(T)
        ok = np.isfinite(v)
        if ok.any() and not ok.all():               # frames where no good bead was fitted: interpolated from the neighbouring frames
            v = np.interp(np.arange(T), np.nonzero(ok)[0], v[ok])
        shifts[name] = np.where(np.isfinite(v), v, 0.0)
    if not three_d:
        shifts["dz_um"] = 0.0
    shifts["n_beads_fitted"] = np.isfinite(res[:, :, 0]).sum(axis=1) if res.shape[1] else 0
    return shifts[["t", "dy_px", "dx_px", "dz_um", "n_beads", "n_beads_fitted"]], res


def consistent_beads(rel: np.ndarray, three_d: bool, max_xy_px: float = 0.5, max_z_um: float = 0.2) -> np.ndarray:
    """Which beads move like the others. All beads lie on the slide, so in every frame they share one displacement; a bead whose own track
    keeps disagreeing with the median of the beads (median absolute difference over frames above max_xy_px or max_z_um, and 3x the most
    consistent bead) is badly fitted: e.g. a bead whose axial profile is cut off by the lowest plane, whose height fit then wanders by
    micrometres (measured: one of three beads did this and, with so few beads, the frame medians jumped by up to 4 um). With fewer than three
    beads there is no majority to judge by, and all are kept."""
    n = rel.shape[1]
    if n < 3:
        return np.ones(n, bool)
    keep = np.ones(n, bool)
    for c, lim in [(0, max_xy_px), (1, max_xy_px)] + ([(3, max_z_um)] if three_d else []):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med = np.nanmedian(rel[:, :, c], axis=1, keepdims=True)
            dev = np.nanmedian(np.abs(rel[:, :, c] - med), axis=0)
        best = np.nanmin(dev) if np.isfinite(dev).any() else 0.0
        keep &= ~((dev > lim) & (dev > 3 * best))
    return keep if keep.any() else np.ones(n, bool)


def correct_with_beads(pos: pd.DataFrame, shifts: pd.DataFrame) -> pd.DataFrame:
    """Subtract the beads' measured whole-picture displacement of each frame from the locus positions (lab_y_px, lab_x_px, z_um)."""
    s = shifts.set_index("t")
    out = pos.copy()
    out["lab_y_px"] = out.lab_y_px - s.loc[out.t, "dy_px"].to_numpy()
    out["lab_x_px"] = out.lab_x_px - s.loc[out.t, "dx_px"].to_numpy()
    out["z_um"] = out.z_um - s.loc[out.t, "dz_um"].to_numpy()
    return out


def correct_with_plane_positions(pos: pd.DataFrame, plane_log: pd.DataFrame, z_nominal) -> pd.DataFrame:
    """Replace the nominal plane heights by the logged ones: a z measured at fractional plane index i is moved to the logged height
    interpolated between the planes. plane_log: columns t, k, readout_um (the piezo's position log)."""
    z0 = np.asarray(z_nominal, float)
    step = float(np.mean(np.diff(z0)))
    logs = {int(t): g.sort_values("k").readout_um.to_numpy() for t, g in plane_log.groupby("t")}
    out = pos.copy()
    idx = (out.z_um.to_numpy() - z0[0]) / step
    new = out.z_um.to_numpy().copy()
    for j, (t, i) in enumerate(zip(out.t.to_numpy(), idx)):
        r = logs.get(int(t))
        if r is not None:
            new[j] = np.interp(i, np.arange(len(r)), r) if 0 <= i <= len(r) - 1 else new[j] + (r[0] - z0[0] if i < 0 else r[-1] - z0[-1])
    out["z_um"] = new
    return out


def _distances(pos: pd.DataFrame, px_um: float) -> pd.DataFrame:
    from simlive.stage1_chromatin.copies import with_copy

    pos = with_copy(pos)
    a = pos[pos.locus_id == 0].set_index(["cell_id", "t", "copy"])
    b = pos[pos.locus_id == 1].set_index(["cell_id", "t", "copy"])
    j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
    d = np.sqrt((px_um * (j.lab_y_px_a - j.lab_y_px_b)) ** 2 + (px_um * (j.lab_x_px_a - j.lab_x_px_b)) ** 2 + (j.z_um_a - j.z_um_b) ** 2)
    dz = j.z_um_a - j.z_um_b
    return pd.DataFrame({"dist_um": d, "dz_um": dz})


def score_variants(variants: dict, truth: pd.DataFrame, px_um: float) -> pd.DataFrame:
    """variants: {name: positions table with cell_id, t, locus_id, lab_y_px, lab_x_px, z_um}. truth: the ORIGINAL locus truth (sample coordinates,
    columns cell_id, t, locus_id, y_um, x_um, z_um). For each variant: median absolute error of the 3D distance between locus 0 and 1, of their
    height difference, and of the height of each locus after the constant offset is removed (a shift common to the
    whole movie is not an error an analysis can or needs to fix)."""
    from simlive.stage1_chromatin.copies import with_copy
    from simlive.stage6_analysis.localize import match_copies

    truth = with_copy(truth)
    tt = truth.assign(lab_y_px=truth.y_um / px_um, lab_x_px=truth.x_um / px_um)[["cell_id", "t", "locus_id", "copy", "lab_y_px", "lab_x_px", "z_um"]]
    td = _distances(tt, px_um)
    rows = []
    for name, pos in variants.items():
        pos = match_copies(pos, truth, px_um)
        d = _distances(pos, px_um).join(td, rsuffix="_true", how="inner")
        j = pos.merge(truth[["cell_id", "t", "locus_id", "copy", "z_um"]], on=["cell_id", "t", "locus_id", "copy"], suffixes=("", "_true"))
        ez = j.z_um - j.z_um_true
        ez = ez - ez.median()
        med = lambda v: 1000 * float(np.median(np.abs(v)))               # typical error: robust against the placeholder detector's rare wrong detections
        rows.append({"correction": name, "n": len(d), "distance error, typical (nm)": med(d.dist_um - d.dist_um_true),
                     "height-difference error, typical (nm)": med(d.dz_um - d.dz_um_true),
                     "height error after removing the constant offset, typical (nm)": med(ez)})
    return pd.DataFrame(rows)
