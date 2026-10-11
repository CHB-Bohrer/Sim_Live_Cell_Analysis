"""Stage 6: find one diffraction-limited spot per colour in a cell, following the detection steps of Bennabi et al. 2025.

Source: Bennabi, ..., Gregor, "Size-dependent temporal decoupling of morphogenesis and transcriptional programs in pseudoembryos",
Science Advances 11, eadv7790 (2025); code github.com/gregorlab/SizeProject2024, `modules/raw_live_analysis.py`
(`get_filtered_fluo_movie_mesp`, `compute_movie_threshold`, `get_fluorescent_area`, `compute_multipol_analysis_per_df`).
They locate Mesp2 expression "poles" (hundreds of pixels across) in live gastruloid movies:

    paper                                                        here (spots ~1 px sigma at 130 nm/px)
    1  median rank filter, k x k square footprint (k = 7)          median filter k x k (median_px, default 3; z-stacks (median_z, k, k))
    2  mask with the gastruloid's brightfield contour              mask with the cell's own (dilated) nuclear mask from stage 5
    -  (no background subtraction)                                 background: "local" (default) = a wider median (background_px) of the
                                                                   filtered image, subtracted; "cell" = one value per cell and frame
    3  movie threshold = median over frames of each frame's        same: per-frame threshold, median over frames (threshold_scope "movie"),
       Otsu threshold (or a fixed value, th = 70)                  per-frame value "noise" (min_snr x robust noise, default), "otsu" (the
                                                                   paper's), or a fixed number
    4  area = pixels above it; holes filled; objects < 5000 px     same, objects < min_area_px removed; erosion optional (erode_px, default 0:
       removed; eroded (disk 10, then disk 2)                      a spot is only a few pixels across)
    5  Gaussian filter sigma = 20 px (pole scale)                  Gaussian filter smooth_sigma_px (~ the PSF sigma: a matched filter)
    6  peak_local_max(min_distance = 100) inside the area:         peak_local_max(min_distance = min_sep_px) inside the area: candidate
       each local maximum is a pole (all are kept and counted)     spots. One is the locus (one locus per colour): select "track" (default,
                                                                   OUR addition) = the candidate nearest the locus' position in the
                                                                   neighbouring frames, so one-frame events (cosmic rays) are skipped;
                                                                   "brightest" = the brightest candidate
    7  (integer pixel coordinates)                                 sub-pixel: 2D (3D) Gaussian fit to the unfiltered pixels around the peak
                                                                   ("gaussian", default), centroid, or none (as the paper)

The paper thresholds the median-filtered image (step 3 before step 5); `threshold_on: smoothed` (default) thresholds the smoothed image
instead (the Gaussian is then also the detection filter, which has the better signal-to-noise ratio for spots). The defaults (local
background, threshold on the smoothed image, track choice) were chosen by scoring variants against the simulation's truth on three runs
(docs/SPOT_FINDING.md); the closest version of the paper's own pipeline is the variant `paper_literal` in compare.py.

Why the scales change: in the paper the median kernel (7 px) is far smaller than the object (hundreds of px), so it removes only pixel-level
noise. A spot here is ~2 px across: a 3 x 3 median already halves its peak (sigma 0.8 px) and is the largest that keeps it, while a 7 x 7
median removes the spot entirely, which makes a wider median a good estimate of the background under a spot (step "background"). The
3 x 3 median removes isolated outliers (hot pixels, single-pixel noise); a cosmic ray is as wide as the PSF and survives it in 2D, but in a
finely sampled z-stack it sits in one plane only and an axial median (median_z = 3) removes it.
Otsu needs a two-peaked histogram (expressing vs non-expressing tissue); a cell region holds thousands of background pixels and a handful
of spot pixels, so Otsu splits the background noise in two and rejects nothing. It is kept as an option to show this.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import binary_erosion, binary_fill_holes, gaussian_filter, label, median_filter
from scipy.optimize import least_squares

MAD_TO_SD = 1.4826
# defaults chosen by scoring variants against the simulation's truth (scripts/compare_localizers.py; see docs/SPOT_FINDING.md)
DEFAULTS = dict(median_px=3, median_z="auto", background="local", background_px=9, smooth_sigma_px=1.0, smooth_sigma_z=0.8,
                threshold="noise", threshold_on="smoothed", threshold_scope="movie", min_snr=5.0, min_area_px=3, erode_px=0, min_sep_px=4,
                select="track", track_window=3, track_gate_px=8.0, refine="gaussian", fit_half_px=3, fit_half_z=2,
                # dog_viterbi (Chen et al. 2025; dog_viterbi.py)
                dog_sigma_px=(1.0, 3.0), dog_sigma_z=(0.8, 2.0), dog_snr=3.0, max_candidates=5, viterbi_max_gap=5, viterbi_miss_cost=2.0,
                viterbi_D_px2=None, viterbi_iterations=3, z_scale=1.0,
                # ellipsoid intensity of every detection (Chen et al. 2025): diameters (x, y, z) in px / planes
                intensity_inner_px=(7, 7, 5), intensity_outer_px=(19, 19, 7), intensity_background="median")


def robust_bg_noise(vals: np.ndarray) -> tuple[float, float]:
    """Background (median) and noise (1.4826 x median absolute deviation) of a set of pixel values, insensitive to a few bright spots."""
    bg = float(np.median(vals))
    return bg, float(MAD_TO_SD * np.median(np.abs(vals - bg))) + 1e-9


def otsu(vals: np.ndarray, nbins: int = 256) -> float:
    """Otsu's threshold of a set of values (the value that best splits them into two classes; same as skimage.filters.threshold_otsu)."""
    v = np.asarray(vals, float).ravel()
    if v.size == 0 or v.min() == v.max():
        return float(v.min()) if v.size else 0.0
    hist, edges = np.histogram(v, nbins)
    centres = 0.5 * (edges[1:] + edges[:-1])
    w0 = np.cumsum(hist).astype(float)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centres) / np.maximum(w0, 1e-12)
    m1 = (np.sum(hist * centres) - np.cumsum(hist * centres)) / np.maximum(w1, 1e-12)
    between = w0[:-1] * w1[:-1] * (m0[:-1] - m1[:-1]) ** 2
    return float(centres[:-1][np.argmax(between)])


def _reg(region: np.ndarray, shape) -> np.ndarray:
    return np.broadcast_to(region, shape) if len(shape) == 3 and region.ndim == 2 else region


def prefilter(img: np.ndarray, region: np.ndarray, median_px: int = 3, median_z: int = 1, background: str = "local",
              background_px: int = 9, smooth_sigma_px: float = 1.0, smooth_sigma_z: float = 0.8, **_) -> dict:
    """Steps 1, 2, background and 5 for one image (y, x) or stack (z, y, x); `region` is the cell's 2D mask (used in every plane).
    Returns the median-filtered image, the background (image or number), the background-subtracted filtered image masked to the region
    (`sub`, zero outside like the paper's frame * mask) and its smoothed version, each with its robust noise inside the region."""
    img = np.asarray(img, np.float32)
    three_d = img.ndim == 3
    reg = _reg(region, img.shape)
    k = max(int(median_px), 1)
    kz = int(median_z) if str(median_z).isdigit() else 1             # "auto" is resolved by the 3D localizer; 1 if it was not
    size = (max(kz, 1), k, k) if three_d else (k, k)
    med = median_filter(img, size=size, mode="nearest") if max(size) > 1 else img
    if background == "local":                        # a median wider than the spot: removes the spot, keeps the slowly varying background
        b = max(int(background_px), 3)
        bgv = median_filter(med, size=(1, b, b) if three_d else (b, b), mode="nearest")
    else:
        bgv = float(np.median(med[reg]))
    sub = np.where(reg, med - bgv, 0.0).astype(np.float32)
    sig = (smooth_sigma_z, smooth_sigma_px, smooth_sigma_px) if three_d else smooth_sigma_px
    sm = gaussian_filter(sub, sig) if smooth_sigma_px > 0 else sub
    return {"med": med, "bg": bgv, "sub": sub, "smooth": sm, "noise": robust_bg_noise(sub[reg])[1],
            "noise_smooth": robust_bg_noise(sm[reg])[1]}


def frame_threshold(pre: dict, region: np.ndarray, threshold="noise", min_snr: float = 5.0, threshold_on: str = "filtered") -> float:
    """Step 3 for one frame: a threshold on `sub` (threshold_on "filtered", as the paper) or on `smooth`, in counts above background."""
    im = pre["smooth"] if threshold_on == "smoothed" else pre["sub"]
    if threshold == "otsu":                                                  # the paper's automatic threshold
        return otsu(im[_reg(region, im.shape)])
    if threshold == "noise":
        return float(min_snr) * pre["noise_smooth" if threshold_on == "smoothed" else "noise"]
    return float(threshold)                                                  # a fixed number of counts above background


def candidates(img: np.ndarray, region: np.ndarray, thr: float | None = None, pre: dict | None = None, *, threshold="noise",
               threshold_on: str = "filtered", min_snr: float = 5.0, min_area_px: int = 3, erode_px: int = 0, min_sep_px: int = 4,
               steps: bool = False, **kw) -> dict | None:
    """Steps 1-6 for one 2D image (y, x) or z-stack (z, y, x) of one cell (`region`: its dilated 2D mask). `thr`: the threshold (the movie
    threshold; None = this frame's own). `pre`: `prefilter` result if already computed (else made from kw). Returns None for a region too
    small to analyse, else {peaks (n, ndim) brightest first, values (smoothed), snr, pre, ...}; with steps=True also the area and threshold."""
    from skimage.feature import peak_local_max

    img = np.asarray(img, np.float32)
    three_d = img.ndim == 3
    reg = _reg(region, img.shape)
    if reg.sum() < 20:
        return None
    if pre is None:
        pre = prefilter(img, region, **kw)
    if thr is None:
        thr = frame_threshold(pre, region, threshold, min_snr, threshold_on)
    sm = pre["smooth"]
    area = ((sm if threshold_on == "smoothed" else pre["sub"]) > thr) & reg        # step 4: the "fluorescent area"
    if not three_d:
        area = binary_fill_holes(area)
    lab, n = label(area)
    if n:
        keep = np.bincount(lab.ravel()) >= int(min_area_px)
        keep[0] = False
        area = keep[lab]
    if erode_px > 0 and area.any():
        st = np.ones((1,) + (2 * erode_px + 1,) * 2, bool) if three_d else np.ones((2 * erode_px + 1,) * 2, bool)
        area = binary_erosion(area, structure=st)
    peaks = (peak_local_max(np.where(area, sm, 0.0), min_distance=max(int(min_sep_px), 1), exclude_border=False,
                            labels=area.astype(np.int32)) if area.any() else np.zeros((0, img.ndim), int))   # step 6: candidates
    vals = sm[tuple(peaks.T)] if len(peaks) else np.zeros(0, np.float32)
    order = np.argsort(vals)[::-1]
    out = {"img": img, "pre": pre, "peaks": peaks[order], "values": vals[order], "snr": vals[order] / pre["noise_smooth"]}
    if steps:
        out.update(area=area, thr=thr)
    return out


def refine_peak(c: dict, j: int = 0, refine: str = "gaussian", fit_half_px: int = 3, fit_half_z: int = 2, smooth_sigma_px: float = 1.0,
                smooth_sigma_z: float = 0.8, **_) -> dict:
    """Step 7 for candidate j of a `candidates` result: sub-pixel position (array index coordinates) and the fit."""
    pre, peak = c["pre"], c["peaks"][j]
    pos, fit = _refine(c["img"] - pre["bg"], pre["sub"], peak, refine, fit_half_px, fit_half_z, smooth_sigma_px, smooth_sigma_z)
    v = c["values"]
    return {"found": True, "pos": pos, "peak": tuple(int(q) for q in peak), "n_candidates": int(len(v)), "snr": float(c["snr"][j]),
            "second_ratio": float(v[1] / v[0]) if len(v) > 1 and v[0] > 0 else 0.0, **fit}


def find_spot(img: np.ndarray, region: np.ndarray, thr: float | None = None, pre: dict | None = None, steps: bool = False,
              **kw) -> dict | None:
    """Steps 1-7 on one image or stack, choosing the brightest candidate (no information from other frames). Returns None for a region too
    small to analyse, {"found": False} if nothing passes, else the refined position etc.; with steps=True every intermediate image too."""
    c = candidates(img, region, thr, pre, steps=steps, **kw)
    if c is None:
        return None
    res = refine_peak(c, 0, **kw) if len(c["peaks"]) else {"found": False}
    if steps:
        res["steps"] = {"raw": c["img"], **c["pre"], "area": c["area"], "thr": c["thr"], "candidates": c["peaks"], "candidate_values": c["values"]}
    return res


def select_track(peaks: list, centre: np.ndarray, window: int = 3, gate_px: float = 8.0, n_iter: int = 2) -> list:
    """Choose one candidate per frame using the neighbouring frames (NOT in the paper; our addition against one-frame events such as cosmic
    rays). `peaks[i]`: candidate positions of frame i, brightest first (last two columns = y, x in crop pixels); `centre[i]`: the nucleus
    centroid in crop pixels. Start from the brightest; then, twice, take in every frame the candidate nearest to the median position (relative
    to the centroid) of the choices in the `window` frames on each side, or none if no candidate lies within `gate_px` of it.
    Returns the chosen index per frame (None = no locus in that frame)."""
    n = len(peaks)
    rel = [np.asarray(p, float)[:, -2:] - centre[i] if len(p) else np.zeros((0, 2)) for i, p in enumerate(peaks)]
    choice = [0 if len(p) else None for p in peaks]
    for _ in range(n_iter):
        pos = np.array([rel[i][c] if c is not None else (np.nan, np.nan) for i, c in enumerate(choice)]).reshape(n, 2)
        new = []
        for i in range(n):
            if not len(rel[i]):
                new.append(None)
                continue
            nb = np.r_[pos[max(0, i - window):i], pos[i + 1:i + window + 1]]
            nb = nb[~np.isnan(nb[:, 0])]
            if not len(nb):
                new.append(0)
                continue
            d = np.hypot(*(rel[i] - np.median(nb, axis=0)).T)
            j = int(np.argmin(d))
            new.append(j if d[j] <= gate_px else None)
        choice = new
    return choice


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


def movie_threshold(per_frame: list[float], scope: str = "movie") -> list[float]:
    """Step 3 across frames: the paper uses ONE threshold per movie, the median of the per-frame thresholds ("movie");
    "frame" keeps each frame's own."""
    v = np.asarray(per_frame, float)
    if scope == "movie" and len(v):
        return [float(np.nanmedian(v))] * len(v)
    return v.tolist()
