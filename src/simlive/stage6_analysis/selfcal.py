"""Stage 6: calibration from the movie itself, for when the calibration images (dark frames, a flat-field slide, a bead slide) are missing.

It builds the same `Calibration` object that `stage3_microscopy/calibration/run.py` builds from calibration acquisitions, so the spot finder
applies it in exactly the same way (`localize_loci(..., calib=<folder>)`). Each part uses only the movie:

- offset / background map, per channel: for every pixel, the median over the frames in which no nucleus is near it. Cells move, so most
  pixels are free some of the time. The result is camera offset + offset pattern (DSNU) + dark current + background light +
  autofluorescence + the excess of hot pixels; pixels that are never free are filled from their neighbours (normalised Gaussian smoothing).
  Subtracting it removes all fixed additive patterns at once (what dark frames do for the camera part only).
- hot pixels: pixels of that map far above their 3 x 3 neighbourhood (more than 6 standard errors of that pixel's median and 5 x the map's robust spread); `Calibration.correct`
  replaces them in every frame by the median of their neighbours.
- flat-field: the nuclear channel's illumination fitted on the moving nuclei (`flatfield.py`, with its cross-validation safeguard); the
  locus colours use the same shape (spot brightness fluctuates too much to map it from the spots).
- chromatic shift between the two locus colours: the mean enhancer - promoter vector over all cells and frames. The true vector points in
  every direction (cells and loci are randomly oriented), so its mean is the shift between the colours. Only the shift BETWEEN colours is
  identifiable this way (locus0 is the reference), which is what distances need. It is applied only with >= 30 cells and if larger than 3 x its standard
  error (bootstrap over cells): with few cells the true vectors do not average out and the estimate would add an error instead.
- colour cross-talk (bleed-through) of the nuclear channel into each locus colour: inside the nuclei, a locus channel's pixels rise with the
  nuclear channel's by the bleed-through fraction (robust regression; the locus spot itself is rejected as an outlier). `channel_inputs`
  then subtracts fraction x the (background-subtracted) nuclear crop from the locus crop (unmixing). This matters most for bright MS2 spots,
  which sit exactly on the promoter and leak into its colour. Applied only if larger than 2 standard errors.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_dilation, gaussian_filter, median_filter

COMPONENTS = ("offset", "hot", "flat", "shift", "crosstalk")


def background_map(movie: np.ndarray, labels: np.ndarray, gap_px: int = 8, min_free: int = 5, fill_sigma: float = 8.0, stat: str = "low",
                   q: float = 20.0):
    """Per-pixel background from the frames in which the pixel is more than gap_px from any nucleus. Returns (map, number of free frames,
    standard error of each pixel's median: 1.25 x its robust frame-to-frame spread / sqrt(free frames)).

    stat "median": the median of those frames. Biased upwards wherever a CELL BODY (cytoplasm) sits next to a nucleus in some of them: the
    nuclei are excluded, their cytoplasm is not. stat "low" (default): the q-th percentile of those frames plus z_q x the pixel's own
    frame-to-frame noise (measured on differences of consecutive free frames, which a slowly passing cell body barely changes): the
    background is the floor that the pixel returns to whenever nothing is there, so cytoplasm present in fewer than ~(100 - q)% of the
    free frames does not bias it (the percentile of pure noise is corrected back to its centre)."""
    import warnings

    T = len(movie)
    st = np.ones((2 * gap_px + 1, 2 * gap_px + 1), bool)
    stack = np.empty(movie.shape, np.float32)
    for t in range(T):
        free = ~binary_dilation(np.asarray(labels[t]) > 0, structure=st)
        stack[t] = np.where(free, np.asarray(movie[t], np.float32), np.nan)
    n_free = np.isfinite(stack).sum(axis=0)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)          # pixels never free of nuclei: filled below
        m = np.nanmedian(stack, axis=0)
        spread = 1.4826 * np.nanmedian(np.abs(stack - m), axis=0)
        if stat == "low":
            from scipy.stats import norm

            d = np.diff(stack, axis=0)                               # NaN unless both frames are free
            sd_t = 1.4826 * np.nanmedian(np.abs(d - np.nanmedian(d, axis=0)), axis=0) / np.sqrt(2)
            sd_t = np.where(np.isfinite(sd_t), sd_t, spread)
            low = np.nanpercentile(stack, q, axis=0) - norm.ppf(q / 100.0) * sd_t
            m_med = m
            m = np.where(np.isfinite(low), low, m_med)
    spread = np.maximum(spread, np.nanmedian(spread))            # integer camera values: a pixel's MAD can be 0; use at least the typical spread
    se = np.where(n_free >= 10, 1.25 * spread / np.sqrt(np.maximum(n_free, 1)), np.inf)     # too few free frames: not tested for hot pixels
    ok = n_free >= min_free
    if not ok.all():                                 # never free of nuclei: fill from the neighbours (normalised convolution)
        num = gaussian_filter(np.where(ok, m, 0.0), fill_sigma)
        den = gaussian_filter(ok.astype(float), fill_sigma)
        fill = np.where(den > 1e-3, num / np.maximum(den, 1e-3), np.nanmedian(m[ok]) if ok.any() else 0.0)
        m = np.where(ok, m, fill)
        se = np.where(ok, se, np.inf)
    return m.astype(np.float32), n_free, se.astype(np.float32)


def hot_pixel_map(bg: np.ndarray, se: np.ndarray | None = None, z: float = 6.0, k_spread: float = 5.0) -> np.ndarray:
    """Pixels whose background is far above their 3 x 3 neighbourhood: by more than z standard errors of that pixel's own median (so noise
    in the map is not taken for a hot pixel) AND more than k_spread x the map's robust pixel-to-pixel spread (so the camera's ordinary offset
    pattern is not either). Pixels never free of nuclei are not tested."""
    r = bg - median_filter(bg, size=3, mode="nearest")
    sd = max(1.4826 * float(np.median(np.abs(r - np.median(r)))), 0.5)    # at least the rounding step of integer camera values
    hot = r > k_spread * sd
    if se is not None:
        hot &= r > z * se
    return hot


def chromatic_shift_from_positions(pos: pd.DataFrame, n_boot: int = 400, seed: int = 0, min_cells: int = 30) -> dict:
    """Mean enhancer - promoter vector (lab px) over all cells and frames where both were found, its bootstrap standard error over cells
    (frames of one cell are correlated, cells are not) and whether to apply it: at least min_cells cells and > 3 SE in either axis."""
    k = ["cell_id", "t"] + (["copy"] if "copy" in pos.columns else [])          # promoter and enhancer of the same copy
    a = pos[pos.locus_id == 0].set_index(k)[["lab_y_px", "lab_x_px"]]
    b = pos[pos.locus_id == 1].set_index(k)[["lab_y_px", "lab_x_px"]]
    d = (b - a).dropna().reset_index()
    if len(d) < 20 or d.cell_id.nunique() < 3:
        return {"shift_px": [0.0, 0.0], "se_px": [np.nan, np.nan], "significant": False, "n_frames": int(len(d)), "n_cells": int(d.cell_id.nunique())}
    per = d.groupby("cell_id")[["lab_y_px", "lab_x_px"]].agg(["sum", "count"])
    sums = np.c_[per[("lab_y_px", "sum")], per[("lab_x_px", "sum")]]
    cnt = per[("lab_y_px", "count")].to_numpy()
    mean = sums.sum(0) / cnt.sum()
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(cnt), len(cnt))
        boots.append(sums[i].sum(0) / cnt[i].sum())
    n = len(cnt)
    se = np.std(boots, axis=0) * np.sqrt(n / max(n - 1, 1))       # the bootstrap of a mean underestimates the SE for few cells
    # applied only with enough cells: with 12 cells an estimate 2.4 SE off the truth was called significant at 2 SE (2026-10-13), because
    # the cell-to-cell scatter of the enhancer-promoter direction is large and a dozen cells say little about its spread
    ok = n >= min_cells and bool(np.any(np.abs(mean) > 3 * se))
    return {"shift_px": [float(mean[0]), float(mean[1])], "se_px": [float(se[0]), float(se[1])], "significant": ok,
            "reason": "applied" if ok else (f"not applied: {n} cells < {min_cells}" if n < min_cells else "not applied: < 3 standard errors"),
            "n_frames": int(len(d)), "n_cells": int(n)}


def crosstalk_fraction(target, source, labels, frames=None, gap_px: int = 8, erode_px: int = 2) -> dict:
    """Bleed-through of `source` (e.g. the bright nuclear channel, with its MS2 spots) into `target` (a locus colour), from the movie: inside
    the nuclei, the target's background-subtracted pixels against the source's, robust least squares through zero (3 rounds rejecting
    residuals > 4 robust sd: the locus spot itself is an outlier). Background of each = median of the pixels far from any nucleus.
    Returns {fraction, se, n_pixels}."""
    from scipy.ndimage import binary_erosion

    T = len(labels)
    frames = range(0, T, max(1, T // 30)) if frames is None else frames
    st = np.ones((2 * gap_px + 1, 2 * gap_px + 1), bool)
    se_ = np.ones((2 * erode_px + 1, 2 * erode_px + 1), bool)
    xs, ys = [], []
    for t in frames:
        L = np.asarray(labels[t]) > 0
        far = ~binary_dilation(L, structure=st)
        inside = binary_erosion(L, structure=se_)
        if far.sum() < 100 or inside.sum() < 100:
            continue
        s, g = np.asarray(source[t], np.float32), np.asarray(target[t], np.float32)
        xs.append(s[inside] - np.median(s[far]))
        ys.append(g[inside] - np.median(g[far]))
    if not xs:
        return {"fraction": 0.0, "se": np.nan, "n_pixels": 0}
    x, y = np.concatenate(xs).astype(np.float64), np.concatenate(ys).astype(np.float64)
    keep = np.ones(len(x), bool)
    f = 0.0
    for _ in range(3):
        f = float((x[keep] * y[keep]).sum() / max((x[keep] ** 2).sum(), 1e-9))
        r = y - f * x
        sd = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep]))) + 1e-9
        keep = np.abs(r) < 4 * sd
    r = (y - f * x)[keep]
    se = float(np.sqrt(np.mean(r ** 2) / max((x[keep] ** 2).sum(), 1e-9)))
    return {"fraction": f, "se": se, "n_pixels": int(keep.sum())}


def build_selfcal(images: dict, labels, out_dir, dt_s: float, components=COMPONENTS, positions: pd.DataFrame | None = None,
                  exposure_s: float = 0.1, px_um: float = 0.13, background_stat: str = "low") -> dict:
    """Estimate the parts in `components` from the movie (images: {channel: path}, labels: label-stack path) and save a Calibration in
    out_dir (parts not estimated are neutral: constant offset, flat 1, no hot pixels, no shift). positions: a first-pass loci_positions table
    (no calibration) for the chromatic shift. Returns a report (what was estimated, sizes, standard errors)."""
    from simlive.stage3_microscopy.calibration.model import Calibration
    from simlive.stage6_analysis.flatfield import fit_flat_bleach, nuclear_signals

    out_dir = Path(out_dir)
    lab = tifffile.imread(str(labels))
    chans = list(images)
    report = {"components": list(components), "channels": {}, "background_stat": background_stat}
    maps = {}
    shape = lab.shape[-2:]
    flat = np.ones(shape, np.float32)
    if "flat" in components and "nucleus" in images:
        nuc = nuclear_signals(images["nucleus"], labels)
        fit = fit_flat_bleach(nuc[~nuc.touches_border].rename(columns={"signal": "value"}), nuc[["cell_id", "t", "y_px", "x_px"]], shape, dt_s)
        flat = fit.map() if fit.fitted_flat else flat
        report["flat"] = {"used": bool(fit.fitted_flat), "rejected": bool(fit.flat_rejected), "cross_validation": fit.cv,
                          "range": [float(flat.min()), float(flat.max())]}
    for ch in chans:
        mov = tifffile.imread(str(images[ch]))
        bg, n_free, se = background_map(mov, lab, stat=background_stat)
        del mov
        hot = hot_pixel_map(bg, se) if "hot" in components else np.zeros(shape, bool)
        if "offset" in components:
            off = np.where(hot, median_filter(bg, size=3, mode="nearest"), bg) if "hot" in components else bg
        else:
            off = np.full(shape, float(np.median(bg)), np.float32)
        maps[ch] = {"offset": off.astype(np.float32), "flat": flat.astype(np.float32), "hot": hot}
        report["channels"][ch] = {"n_hot": int(hot.sum()), "background_median": float(np.median(bg)),
                                  "background_spread": float(np.percentile(bg, 95) - np.percentile(bg, 5)),
                                  "pixels_never_free": float(np.mean(n_free < 5))}
    beads = {"shift_px": {c: [0.0, 0.0] for c in chans}}
    if "shift" in components and positions is not None and len(positions):
        cs = chromatic_shift_from_positions(positions)
        report["shift"] = cs
        if cs["significant"]:
            beads["shift_px"]["locus1"] = cs["shift_px"]
    values = {}
    if "crosstalk" in components and "nucleus" in images:
        src = tifffile.imread(str(images["nucleus"]))
        report["crosstalk"] = {}
        for ch in chans:
            if not ch.startswith("locus"):
                continue
            ct = crosstalk_fraction(tifffile.imread(str(images[ch])), src, lab)
            ct["applied"] = bool(ct["fraction"] > 2 * (ct["se"] if np.isfinite(ct["se"]) else np.inf))
            report["crosstalk"][ch] = ct
            if ct["applied"]:
                values.setdefault("crosstalk", {})[ch] = {"nucleus": ct["fraction"]}
        del src
    cal = Calibration(chans, tuple(int(v) for v in shape), float(exposure_s), float(px_um), "locus0", maps, values, beads, [])
    cal.save(out_dir)
    (out_dir / "selfcal_report.json").write_text(json.dumps(report, indent=1, default=float))
    return report
