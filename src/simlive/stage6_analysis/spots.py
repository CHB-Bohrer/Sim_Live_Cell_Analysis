"""Stage 6: helpers shared by the spot finders (localize.py, localize3d.py, dog_viterbi.py): robust noise, sub-pixel refinement of a chosen
candidate (2D / 3D Gaussian fit, centroid, or the integer pixel), and the default settings of the `locus_analysis:` config block.

History: a second method, `median_peaks` (filter / threshold / local-maxima steps adapted from Bennabi et al. 2025, Sci. Adv., who detect
Mesp2 expression poles, not spots), lived here until 2026-10-11. It was removed at the user's request because that paper does not locate
spots; its useful parts were generic. Its measurements are kept in docs/SPOT_FINDING.md (git history has the code).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

MAD_TO_SD = 1.4826
# defaults of the locus_analysis block (dog_viterbi, Chen et al. 2025: see dog_viterbi.py); chosen by scoring against the simulation's truth
DEFAULTS = dict(median_px=3, median_z="auto", refine="gaussian", fit_half_px=3, fit_half_z=2, smooth_sigma_px=1.0, smooth_sigma_z=0.8,
                min_sep_px=4,
                dog_sigma_px=(1.0, 3.0), dog_sigma_z=(0.8, 2.0), dog_snr=3.0, max_candidates=5, viterbi_max_gap=5, viterbi_miss_cost=2.0,
                viterbi_D_px2=None, viterbi_iterations=3, z_scale=1.0,
                # ellipsoid intensity of every detection (Chen et al. 2025): diameters (x, y, z) in px / planes
                intensity_inner_px=(7, 7, 5), intensity_outer_px=(19, 19, 7), intensity_background="median")


def robust_bg_noise(vals: np.ndarray) -> tuple[float, float]:
    """Background (median) and noise (1.4826 x median absolute deviation) of a set of pixel values, insensitive to a few bright spots."""
    bg = float(np.median(vals))
    return bg, float(MAD_TO_SD * np.median(np.abs(vals - bg))) + 1e-9

def _reg(region: np.ndarray, shape) -> np.ndarray:
    return np.broadcast_to(region, shape) if len(shape) == 3 and region.ndim == 2 else region

def refine_peak(c: dict, j: int = 0, refine: str = "gaussian", fit_half_px: int = 3, fit_half_z: int = 2, smooth_sigma_px: float = 1.0,
                smooth_sigma_z: float = 0.8, **_) -> dict:
    """Step 7 for candidate j of a `candidates` result: sub-pixel position (array index coordinates) and the fit."""
    pre, peak = c["pre"], c["peaks"][j]
    pos, fit = _refine(c["img"] - pre["bg"], pre["sub"], peak, refine, fit_half_px, fit_half_z, smooth_sigma_px, smooth_sigma_z)
    v = c["values"]
    return {"found": True, "pos": pos, "peak": tuple(int(q) for q in peak), "n_candidates": int(len(v)), "snr": float(c["snr"][j]),
            "second_ratio": float(v[1] / v[0]) if len(v) > 1 and v[0] > 0 else 0.0, **fit}

def _window(shape, centre, halves):
    sl = tuple(slice(max(c - h, 0), min(c + h + 1, n)) for c, h, n in zip(centre, halves, shape))
    return sl, np.meshgrid(*[np.arange(s.start, s.stop) for s in sl], indexing="ij")

def _refine(raw, sub, peak, how: str, half: int, half_z: int, sigma_px: float, sigma_z: float):
    """Step 7: sub-pixel position from the background-subtracted, unfiltered pixels around the chosen peak."""
    nd = raw.ndim
    halves = (half_z, half, half) if nd == 3 else (half, half)
    amp_med = float(sub[tuple(peak)])
    if how == "none":                                                        # the paper: integer pixel coordinates of the peak
        return tuple(float(v) for v in peak), {"fit_ok": True, "amp": amp_med, "sigma_fit_px": np.nan}
    sl, grids = _window(raw.shape, peak, halves)
    w = raw[sl].astype(float)
    if how == "centroid" or w.size < 2 * (nd + 3):
        return _centroid(raw, peak, nd), {"fit_ok": how == "centroid", "amp": amp_med, "sigma_fit_px": np.nan}
    g = [gr.ravel().astype(float) for gr in grids]
    v = w.ravel()
    amp0 = max(float(raw[tuple(peak)]), 1e-3)
    s0 = max(float(sigma_px), 0.5)
    if nd == 2:
        def model(p):
            return p[4] + p[0] * np.exp(-((g[0] - p[1]) ** 2 + (g[1] - p[2]) ** 2) / (2 * p[3] ** 2))
        p0 = [amp0, peak[0], peak[1], s0, 0.0]
        lo = [0.0, peak[0] - 1.5, peak[1] - 1.5, 0.3, -np.inf]
        hi = [np.inf, peak[0] + 1.5, peak[1] + 1.5, 4 * s0, np.inf]
    else:
        def model(p):
            return p[6] + p[0] * np.exp(-((g[1] - p[2]) ** 2 + (g[2] - p[3]) ** 2) / (2 * p[4] ** 2) - (g[0] - p[1]) ** 2 / (2 * p[5] ** 2))
        sz0 = max(float(sigma_z), 0.4)
        p0 = [amp0, peak[0], peak[1], peak[2], s0, sz0, 0.0]
        lo = [0.0, peak[0] - 1.5, peak[1] - 1.5, peak[2] - 1.5, 0.3, 0.2, -np.inf]
        hi = [np.inf, peak[0] + 1.5, peak[1] + 1.5, peak[2] + 1.5, 4 * s0, 6 * sz0, np.inf]
    p0 = np.clip(p0, np.array(lo) + 1e-6, np.array(hi) - 1e-6)
    try:
        r = least_squares(lambda p: model(p) - v, p0, bounds=(lo, hi), x_scale="jac", max_nfev=200)
        p = r.x
        on_edge = any(abs(p[1 + i] - peak[i]) > 1.45 for i in range(nd))      # centre ran to the bound: not a clean spot
        if not r.success or on_edge:
            raise ValueError
    except (ValueError, np.linalg.LinAlgError):
        return _centroid(raw, peak, nd), {"fit_ok": False, "amp": amp_med, "sigma_fit_px": np.nan}
    return tuple(float(p[1 + i]) for i in range(nd)), {"fit_ok": True, "amp": float(p[0]), "sigma_fit_px": float(p[1 + nd])}

def _centroid(raw, peak, nd, win: int = 2):
    """Intensity-weighted centroid of the positive background-subtracted pixels in a (2 win + 1) window (+-1 plane in z)."""
    halves = (1, win, win) if nd == 3 else (win, win)
    sl, grids = _window(raw.shape, peak, halves)
    w = np.clip(raw[sl], 0, None).astype(float)
    if w.sum() <= 0:
        return tuple(float(v) for v in peak)
    return tuple(float((w * gr).sum() / w.sum()) for gr in grids)
