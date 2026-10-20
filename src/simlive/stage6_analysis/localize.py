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

Several copies of the promoter-enhancer pair per cell (`locus_analysis.max_copies` > 1, or the pair gate `pair_gate_um`): all colours of a
cell go through multicopy.track_copies (several paths per colour, stray spots rejected by their short life, promoter paths paired with
enhancer paths); every detection then carries its `copy` number (0 for the one-copy finder). `match_copies` relabels the found copies to
the true ones before scoring.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

from simlive.stage6_analysis import dog_viterbi as DV
from simlive.stage6_analysis import multicopy as MC
from simlive.stage6_analysis import spots as SP

METHODS = ("dog_viterbi", "brightest")
EXTRA_COPY = 50          # match_copies: a found copy that matches no true copy is numbered EXTRA_COPY + its own number
COLUMNS = ["cell_id", "t", "locus_id", "copy", "crop_y_px", "crop_x_px", "lab_y_px", "lab_x_px", "u_um", "v_um", "snr", "n_candidates",
           "fit_ok", "intensity"]


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


def copy_params(kw: dict) -> dict:
    """The multi-copy linker's settings (multicopy.SETTINGS) from a locus_analysis config block."""
    return {k: kw.get(k, v) for k, v in MC.SETTINGS.items()}


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
    cands = channel_candidates(stack, region, p)
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


def channel_candidates(stack, region, p: dict) -> list:
    """DoG candidates of every frame of one colour of one cell (None where the cell region is too small)."""
    return [DV.dog_candidates(stack[i], region[i], steps=True, **p) if region[i].sum() >= 20 else None for i in range(len(stack))]


def _fit(stack_i, c: dict, j: int, p: dict) -> dict:
    f = SP.refine_peak(c, j, **p)
    f["intensity"] = DV.ellipsoid_intensity(stack_i, f["pos"], p["intensity_inner_px"], p["intensity_outer_px"], p["intensity_background"])[0]
    return f


def analyse_cell_copies(inputs: list[dict], px_um: float = 0.13, **kw) -> dict:
    """The multi-copy linker on all colours of one cell. `inputs[k]` = channel_inputs of colour k. Returns {"candidates": per colour the
    candidates of every frame (as analyse_channel), "copies": [{colour: fits per frame (None = no spot)}], "track": multicopy.track_copies
    result (paths, pairs, choices), "params", "copy_params"}."""
    p, cp = spot_params(kw), copy_params(kw)
    cands = [channel_candidates(ci["stack"], ci["region"], p) for ci in inputs]
    centre = inputs[0]["centre"]
    valid = np.array([c is not None for c in cands[0]])
    area = float(np.mean(inputs[0]["region"].sum(axis=(-2, -1))))
    tr = MC.track_copies([[c["peaks"] if c is not None else np.zeros((0, 2)) for c in ck] for ck in cands],
                         [[c["snr"] if c is not None else np.zeros(0) for c in ck] for ck in cands], centre, valid,
                         px_um=px_um, **cp, **{**p, "search_volume": area})
    copies = []
    for cpy in tr["copies"]:
        copies.append({k: [None if j is None or cands[k][i] is None else _fit(inputs[k]["stack"][i], cands[k][i], j, p)
                           for i, j in enumerate(ch)] for k, ch in cpy.items()})
    light = [[None if c is None else {q: c[q] for q in ("peaks", "values", "snr")} for c in ck] for ck in cands]
    return {"candidates": light, "copies": copies, "track": tr, "params": p, "copy_params": cp}


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


def row_offsets_removed(stack: np.ndarray, mask: np.ndarray, gap_px: int = 6) -> np.ndarray:
    """Camera row (banding) noise from the movie itself, no dark frames: in every frame, each row's offset = the 25th percentile of its
    pixels away from the cell (other nuclei rarely cover a quarter of a row) minus the same over the whole crop; subtracted row by row."""
    from scipy.ndimage import binary_dilation

    out = stack.astype(np.float32).copy()
    st = np.ones((2 * gap_px + 1, 2 * gap_px + 1), bool)
    for i in range(len(out)):
        far = ~binary_dilation(mask[i] > 0, structure=st)
        if far.sum() < 50:
            continue
        v = np.where(far, out[i], np.nan)
        with np.errstate(all="ignore"):
            row = np.nanpercentile(v, 25, axis=1)
        ref = np.nanpercentile(v, 25)
        ok = np.isfinite(row) & (far.sum(axis=1) >= 8)
        out[i] -= np.where(ok, row - ref, 0.0)[:, None]
    return out


def channel_inputs(cell, name: str, dilate_px: int = 2, calib: str | None = None, row_correct: bool = False) -> dict:
    """What the localizer works on for one colour of one cell: the crops (calibration-corrected if `calib`; camera row offsets removed if
    `row_correct`), the analysis region, the nucleus centroid in crop pixels (the crop is centred on the ROUNDED centroid, so up to 0.5 px off
    the crop centre) and the chromatic shift to remove."""
    cal = _load_calibration(calib) if calib else None
    fr = cell.frames.reset_index(drop=True)
    stack = cell.channels[name].astype(np.float32)
    if cal is not None:
        stack = np.stack([cal.correct(stack[i], name, int(r.crop_y0), int(r.crop_x0)) for i, r in enumerate(fr.itertuples())])
        for src, f in ((cal.values.get("crosstalk") or {}).get(name) or {}).items():     # unmixing: remove the bleed-through of `src`
            if src in cell.channels and f:
                from scipy.ndimage import binary_dilation

                s = np.stack([cal.correct(cell.channels[src][i].astype(np.float32), src, int(r.crop_y0), int(r.crop_x0))
                              for i, r in enumerate(fr.itertuples())])
                far = ~binary_dilation(cell.mask, structure=np.ones((1, 17, 17), bool))
                bg = np.array([np.median(s[i][far[i]]) if far[i].any() else 0.0 for i in range(len(s))], np.float32)
                stack = stack - float(f) * (s - bg[:, None, None])
    if row_correct:
        stack = row_offsets_removed(stack, cell.mask)
    return {"stack": stack, "region": cell.analysis_mask(dilate_px), "frames": fr, "centre": centroid_in_crop(fr, cell.crop_size),
            "shift": cal.shift_px(name) if cal is not None else (0.0, 0.0)}


def centroid_in_crop(fr: pd.DataFrame, S: int) -> np.ndarray:
    """(n, 2) nucleus centroid in crop pixels from stage 5's frames.csv (the crop centre if the table has no centroid columns)."""
    if {"cy_px", "cx_px"} <= set(fr.columns):
        return np.c_[fr.cy_px.to_numpy() - fr.crop_y0.to_numpy(), fr.cx_px.to_numpy() - fr.crop_x0.to_numpy()]
    return np.full((len(fr), 2), S // 2, float)


def localize_loci(cell, n_loci: int = 2, method: str = "dog_viterbi", dilate_px: int = 2, px_um: float = 0.13, calib: str | None = None,
                  row_correct: bool = False, **kw) -> pd.DataFrame:
    """Locate every locus (one colour channel each) in every frame of one isolated cell. One row per detection.
    With `calib` (folder of a saved Calibration, from calibration images or from the movie itself: selfcal.py) every crop is first
    dark-subtracted, flat-fielded and hot-pixel-corrected, and the measured chromatic shift of the colour is subtracted from the detected
    position; `row_correct` removes camera row offsets estimated in the crop. Other keywords: the method's settings (see `spots.DEFAULTS`)."""
    chans = cell.locus_channels[:n_loci]
    if not chans:
        return pd.DataFrame(columns=COLUMNS)
    fr = cell.frames.reset_index(drop=True)
    # nucleus long-axis angle (from mask moments) is only defined modulo 180 deg: unwrap so (u, v) is continuous in time
    theta = 0.5 * np.unwrap(2 * fr.theta_rad.to_numpy())
    rows = []

    def emit(k, copy, i, d, shift, cen):
        y, x = d[0] - shift[0], d[1] - shift[1]
        c, s = np.cos(theta[i]), np.sin(theta[i])
        dy, dx = y - cen[i, 0], x - cen[i, 1]
        r = fr.iloc[i]
        rows.append((cell.cell_id, int(r.t), k, copy, y, x, r.crop_y0 + y, r.crop_x0 + x, (dx * c + dy * s) * px_um, (dy * c - dx * s) * px_um,
                     d[2], d[3], bool(d[4]), d[5]))

    inputs = [channel_inputs(cell, name, dilate_px, calib, row_correct) for name in chans]
    if method == "dog_viterbi" and MC.multicopy_on(kw):          # several copies per cell and / or stray spots: all colours together
        a = analyse_cell_copies(inputs, px_um=px_um, **kw)
        for copy, fits in enumerate(a["copies"]):
            for k, ff in fits.items():
                for i, f in enumerate(ff):
                    if f is not None:
                        emit(k, copy, i, (f["pos"][-2], f["pos"][-1], f["snr"], f["n_candidates"], f["fit_ok"], f["intensity"]),
                             inputs[k]["shift"], inputs[k]["centre"])
        return pd.DataFrame(rows, columns=COLUMNS)
    for k, ci in enumerate(inputs):
        for i, d in enumerate(detect_channel(ci["stack"], ci["region"], method, ci["centre"], **kw)):
            if d is not None:
                emit(k, 0, i, d, ci["shift"], ci["centre"])
    return pd.DataFrame(rows, columns=COLUMNS)


def locus_distances(pos: pd.DataFrame, px_um: float) -> pd.DataFrame:
    """Distance between locus 0 and locus 1 of the same copy per cell and frame (frames where both were found), in um."""
    pos = with_copy(pos)
    a = pos[pos.locus_id == 0].set_index(["cell_id", "t", "copy"])
    b = pos[pos.locus_id == 1].set_index(["cell_id", "t", "copy"])
    j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
    d = px_um * np.hypot(j.lab_y_px_a - j.lab_y_px_b, j.lab_x_px_a - j.lab_x_px_b)
    return d.rename("dist_um").reset_index()[["cell_id", "t", "dist_um", "copy"]]


def match_copies(pos: pd.DataFrame, truth: pd.DataFrame, px_um: float, cols=("lab_y_px", "lab_x_px")) -> pd.DataFrame:
    """Relabel the analysis' copies to the true copies, per cell: the Hungarian assignment on the median distance (over all frames and loci
    both have) between every found copy and every true copy. Adds `copy_found` (the analysis' own number, kept if already there); `copy`
    becomes the matched true copy, or EXTRA_COPY + copy_found for a found copy with no true copy left (an extra copy: a stray spot or noise
    taken for a locus). One copy on both sides (every run before copies existed) maps 0 -> 0. Tables without locus_id (MS2 traces) are
    matched to locus 0 (the promoter)."""
    from scipy.optimize import linear_sum_assignment

    pos, truth = with_copy(pos).copy().reset_index(drop=True), with_copy(truth)
    no_locus = "locus_id" not in pos.columns
    if no_locus:
        pos["locus_id"] = 0
    if "copy_found" not in pos.columns:
        pos["copy_found"] = pos["copy"]
    if pos.empty:
        return pos.drop(columns="locus_id") if no_locus else pos
    est_y, est_x = pos[cols[0]] * px_um, pos[cols[1]] * px_um
    pos = pos.assign(_ey=est_y.to_numpy(), _ex=est_x.to_numpy())
    tr = truth[["cell_id", "t", "locus_id", "copy", "y_um", "x_um"]]
    new = pd.Series(index=pos.index, dtype=int)
    for cid, g in pos.groupby("cell_id"):
        tq = tr[tr.cell_id == cid]
        found, true = sorted(g.copy_found.unique()), sorted(tq["copy"].unique())
        if len(found) == 1 and len(true) == 1:
            new[g.index] = true[0]
            continue
        C = np.full((len(found), len(true)), 1e6)
        for a, f in enumerate(found):
            gf = g[g.copy_found == f]
            for b, q in enumerate(true):
                j = gf[["t", "locus_id", "_ey", "_ex"]].merge(tq[tq["copy"] == q][["t", "locus_id", "y_um", "x_um"]], on=["t", "locus_id"])
                if len(j):
                    C[a, b] = float(np.median(np.hypot(j._ey - j.y_um, j._ex - j.x_um)))
        rr, cc = linear_sum_assignment(C)
        m = {found[a]: true[b] for a, b in zip(rr, cc)}
        new[g.index] = [m.get(f, EXTRA_COPY + int(f)) for f in g.copy_found]
    pos["copy"] = new.astype(int)
    return pos.drop(columns=["_ey", "_ex"] + (["locus_id"] if no_locus else []))


def with_copy(df: pd.DataFrame) -> pd.DataFrame:
    from simlive.stage1_chromatin.copies import with_copy as _wc

    return _wc(df)


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


def _nearest_other(j: pd.DataFrame, objects: pd.DataFrame, keys: list) -> np.ndarray:
    """Distance (um) from every detection in j (est_y_um, est_x_um) to the nearest object of `objects` (y_um, x_um) with the same keys."""
    if objects is None or not len(objects) or j.empty:
        return np.full(len(j), np.inf)
    o = objects[keys + ["y_um", "x_um"]].rename(columns={"y_um": "_oy", "x_um": "_ox"})
    m = j[keys + ["est_y_um", "est_x_um"]].reset_index().merge(o, on=keys, how="inner")
    if m.empty:
        return np.full(len(j), np.inf)
    m["_d"] = np.hypot(m.est_y_um - m._oy, m.est_x_um - m._ox)
    d = m.groupby("index")._d.min()
    return d.reindex(j.index).fillna(np.inf).to_numpy()


def copy_scores(pos: pd.DataFrame, truth: pd.DataFrame) -> dict:
    """How well the number of copies per cell was found: share of cells with the right number, true copies never found, extra copies."""
    t = with_copy(truth)
    t = t[_visible(t)].groupby("cell_id")["copy"].nunique()
    f = with_copy(pos).groupby("cell_id")["copy"]
    n_found = f.nunique().reindex(t.index).fillna(0)
    n_matched = pos[pos["copy"] < EXTRA_COPY].groupby("cell_id")["copy"].nunique().reindex(t.index).fillna(0) if "copy" in pos else n_found
    return {"n_cells": int(len(t)), "true_copies": int(t.sum()), "found_copies": int(n_found.sum()),
            "copy_number_right": float((n_found == t).mean()) if len(t) else np.nan,
            "missed_copies": int((t - n_matched).clip(lower=0).sum()), "extra_copies": int((n_found - n_matched).clip(lower=0).sum())}


def score_localization(pos: pd.DataFrame, truth: pd.DataFrame, px_um: float, far_um: float = 0.5,
                       strays: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict]:
    """Compare detections with the true locus positions (same cell, frame, locus and copy: found copies are first matched to true copies,
    `match_copies`). Valid for true cell IDs only.
    detection_rate = detections / visible loci (a locus with no light cannot be found); `n_false_on_dark` = detections where no spot was
    drawn; `frac_far` = share of detections more than `far_um` from the truth (a wrong object was picked: hot pixel, cosmic ray, ghost);
    `frac_on_stray` = share of detections that are far from their locus but within `far_um` of a stray spot of that colour (`strays` =
    stage3_microscopy/stray_spots.csv); `frac_on_other_copy` = far from their locus but within `far_um` of the same locus of another copy;
    `extra_detections` = detections of found copies that match no true copy; copy-number scores (`copy_scores`);
    `intensity_r_locusK` = correlation of the measured ellipsoid intensity with the true photons in that frame (if both are known)."""
    truth = with_copy(truth)
    pos = match_copies(pos, truth, px_um)
    est = pos.assign(est_y_um=pos.lab_y_px * px_um, est_x_um=pos.lab_x_px * px_um)
    est = est.drop(columns=[c for c in truth.columns if c in est.columns and c not in ("cell_id", "t", "locus_id", "copy")])
    tr = truth.assign(visible=_visible(truth))
    j = est.merge(tr, on=["cell_id", "t", "locus_id", "copy"], how="inner")
    j["err_um"] = np.hypot(j.est_y_um - j.y_um, j.est_x_um - j.x_um)
    err = j[["cell_id", "t", "locus_id", "copy", "err_um", "visible"]]
    if err.empty:
        return err, {"n_detections": int(len(pos)), **copy_scores(pos, truth)}
    vis = j.visible.to_numpy(bool)
    far = j.err_um.to_numpy() > far_um
    other = tr[["cell_id", "t", "locus_id", "copy", "y_um", "x_um"]].rename(columns={"copy": "_c"})
    oc = j[["cell_id", "t", "locus_id", "copy", "est_y_um", "est_x_um"]].reset_index().merge(other, on=["cell_id", "t", "locus_id"])
    oc = oc[oc["copy"] != oc["_c"]]
    d_other = (np.hypot(oc.est_y_um - oc.y_um, oc.est_x_um - oc.x_um).groupby(oc["index"]).min().reindex(j.index).fillna(np.inf).to_numpy()
               if len(oc) else np.full(len(j), np.inf))
    if strays is not None and len(strays):
        st = strays[strays.channel.str.startswith("locus")].assign(locus_id=lambda d: d.channel.str[5:].astype(int))
        sy, sx = ("y_drawn_um", "x_drawn_um") if "y_drawn_um" in st else ("y_um", "x_um")
        st = st.assign(y_um=st[sy], x_um=st[sx])
        d_stray = _nearest_other(j, st, ["cell_id", "t", "locus_id"])
    else:
        d_stray = np.full(len(j), np.inf)
    ok = err[vis]
    n_vis = int(tr.visible.sum())
    nv = max(int(vis.sum()), 1)
    out = {"n_detections": len(err), "detection_rate": len(ok) / max(n_vis, 1), "n_visible": n_vis,
           "n_false_on_dark": int((~vis).sum()),
           "rms_error_um": float(np.sqrt((ok.err_um ** 2).mean())) if len(ok) else np.nan,
           "median_error_um": float(ok.err_um.median()) if len(ok) else np.nan,
           "p95_error_um": float(ok.err_um.quantile(0.95)) if len(ok) else np.nan,
           "frac_far": float((ok.err_um > far_um).mean()) if len(ok) else np.nan,
           "frac_on_stray": float((vis & far & (d_stray <= far_um)).sum() / nv),
           "frac_on_other_copy": float((vis & far & (d_other <= far_um)).sum() / nv),
           "extra_detections": int((pos["copy"] >= EXTRA_COPY).sum()),
           "rms_error_near_um": float(np.sqrt((ok.err_um[ok.err_um <= far_um] ** 2).mean())) if (ok.err_um <= far_um).any() else np.nan,
           **copy_scores(pos, truth), **intensity_scores(j, j.err_um <= far_um)}
    err = err.assign(on_stray=far & (d_stray <= far_um), on_other_copy=far & (d_other <= far_um))
    return err, out
