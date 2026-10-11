"""Blind estimation of the instrument from calibration images. These functions see ONLY the images (counts), never the simulation's
true parameters; run.py compares the estimates with the truth afterwards.

    fit_dark   (dark frames at several exposures)  -> offset map, dark current, read noise, row noise, fixed offset pattern (DSNU), hot pixels,
                                                       cosmic-ray rate and charge
    fit_ptc    (flat slide at several brightnesses) -> gain (ADU per electron) from the photon transfer curve, read-noise cross-check, flicker
    fit_flat   (brightest flat series)              -> illumination pattern (flat-field), vignetting, pixel gain differences (PRNU)
    fit_beads  (multi-colour bead slide, time lapse) -> PSF width per colour, chromatic shift (lateral and axial), stage / focus drift

Not measurable from these images: the quantum efficiency (needs a source of known photon flux). The gain here is ADU per detected electron,
so every "electron" number below is in detected electrons (after QE).
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, label, maximum_filter
from scipy.optimize import curve_fit
from scipy.special import erf


def _linfit(x, Y):
    """Closed-form straight-line fit y = a + b x along axis 0 of Y (no BLAS/LAPACK: safe before torch is imported). Returns (a, b)."""
    x = np.asarray(x, float)
    xm = x.mean()
    sxx = ((x - xm) ** 2).sum()
    Ym = Y.mean(axis=0)
    shape = (len(x),) + (1,) * (Y.ndim - 1)
    b = (((x - xm).reshape(shape)) * (Y - Ym)).sum(axis=0) / sxx
    return Ym - b * xm, b


def _robust_sd(v) -> float:
    v = np.asarray(v).ravel()
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


# ------------------------------------------------------------------------------------------------------------------ dark frames
def _robust_moments(S: np.ndarray, nsigma: float = 6.0):
    """Per-pixel mean and variance over frames with rare outliers (cosmic rays) rejected: values further than nsigma robust standard
    deviations from the pixel's median are ignored. Returns (mean, variance, kept-mask)."""
    med = np.median(S, axis=0)
    dev = S - med
    s = max(1.4826 * float(np.median(np.abs(dev))), 1e-6)
    ok = np.abs(dev) <= nsigma * s
    cnt = ok.sum(axis=0)
    mean = (S * ok).sum(axis=0) / np.maximum(cnt, 1)
    var = (((S - mean) ** 2) * ok).sum(axis=0) / np.maximum(cnt - 1, 1)
    return mean, var, ok


def fit_dark(stacks: dict, run_exposure_s: float, hot_nsigma: float = 6.0) -> dict:
    """stacks: {exposure_s: (n, ny, nx) dark frames}. Returns the estimates and the maps (see keys)."""
    exps = sorted(stacks)
    e0 = exps[0]
    means = np.stack([_robust_moments(stacks[e].astype(np.float64))[0] for e in exps])
    a, b = _linfit(exps, means)                                   # per pixel: counts = a (offset incl. hot pixels) + b x exposure
    offset_run = a + b * run_exposure_s
    S0 = stacks[e0].astype(np.float64)
    n0, ny, nx = S0.shape
    m0 = means[0]
    sd0 = _robust_sd(m0)
    hot = (m0 - np.median(m0)) > max(hot_nsigma * sd0, 3.0)      # a hot pixel is high at every exposure: use the shortest (little dark current)
    quantized = bool(np.mean(np.abs(S0[:2] - np.rint(S0[:2]))) < 1e-6)
    v = S0.var(axis=0, ddof=1)
    v_tot = float(v[v < 8 * np.median(v)].mean())                 # per-pixel temporal variance (read + row noise + digitisation), cosmic hits rejected
    # row noise: every row of every frame carries one common random offset. The median over a row ignores cosmic-ray blobs.
    Rm = np.median(S0 - S0.mean(axis=0), axis=2)                  # (n0, ny)
    var_rows = Rm.var(ddof=0) * n0 / (n0 - 1)
    sig_row2 = max(var_rows - (np.pi / 2) * v_tot / nx, 0.0)
    sig_read2 = max(v_tot - sig_row2 - (1 / 12 if quantized else 0.0), 0.0)
    dsnu2 = max(np.var(m0[~hot]) - v_tot / n0, 0.0)
    # cosmic rays: bright, rare, tiny blobs that are in one frame only
    ev_peak, ev_sum, n_frames = [], [], 0
    thr = max(7.0 * np.sqrt(v_tot), 15.0)
    for e in exps:
        S = stacks[e].astype(np.float64)
        m = np.median(S, axis=0)
        for f in S:
            n_frames += 1
            lab, n = label((f - m) > thr, structure=np.ones((3, 3)))
            for k in range(1, n + 1):
                ys, xs = np.nonzero(lab == k)
                i = np.argmax((f - m)[ys, xs])
                y, x = ys[i], xs[i]
                win = (f - m)[max(y - 3, 0):y + 4, max(x - 3, 0):x + 4]
                ev_peak.append(float((f - m)[y, x])); ev_sum.append(float(win.sum()))
    cosmic_rate = len(ev_peak) / max(n_frames, 1)
    return {
        "offset_run_map": offset_run.astype(np.float32), "offset_zero_map": a.astype(np.float32), "hot_mask": hot,
        "offset_adu": float(np.median(offset_run[~hot])), "dark_rate_adu_s": float(np.median(b[~hot])),
        "read_noise_adu": float(np.sqrt(sig_read2)), "row_noise_adu": float(np.sqrt(sig_row2)), "dsnu_adu": float(np.sqrt(dsnu2)),
        "hot_fraction": float(hot.mean()), "n_hot": int(hot.sum()),
        "hot_excess_adu": float(np.median((m0 - np.median(m0))[hot])) if hot.any() else float("nan"),
        "cosmic_rate_per_frame": cosmic_rate, "cosmic_events": len(ev_peak), "cosmic_frames": n_frames,
        "cosmic_integral_adu": float(np.median(ev_sum)) if ev_sum else float("nan"),
        "quantized": quantized, "temporal_var_adu2": v_tot, "n_frames": int(n0),
    }


# ------------------------------------------------------------------------------------------------------------------ photon transfer
def fit_ptc(levels: dict, offset_map: np.ndarray, roi: int = 256, full_scale_adu: float = 65535.0) -> dict:
    """levels: {nominal brightness: (n, ny, nx) flat-slide frames}. Photon transfer curve: temporal variance V vs mean signal M (both in counts).
    The common brightness flicker of each frame is divided out and cosmic-ray outliers are rejected; slope = gain (ADU per detected electron),
    intercept = noise floor."""
    pts = []
    for key in sorted(levels):
        S = levels[key].astype(np.float64) - offset_map[None]
        n, ny, nx = S.shape
        y0, x0 = max(ny // 2 - roi // 2, 0), max(nx // 2 - roi // 2, 0)
        R = S[:, y0:y0 + roi, x0:x0 + roi]
        mt = R.mean(axis=(1, 2))
        if abs(mt.mean()) > 20:
            R = R / (mt / mt.mean())[:, None, None]               # divide out the common-mode flicker
        mean, var, _ = _robust_moments(R)
        pts.append((float(mean.mean()), float(var.mean())))
    P = np.array(pts)
    ok = (P[:, 0] > -1) & (P[:, 0] < 0.8 * full_scale_adu)
    c, g = _linfit(P[ok, 0], P[ok, 1][:, None])
    return {"gain_adu_per_e": float(g[0]), "ptc_intercept_adu2": float(c[0]), "ptc_points": P.tolist()}


def fit_flicker(stack: np.ndarray, offset_map: np.ndarray, roi: int = 256) -> dict:
    """Brightness flicker from one long series of a uniform slide: relative sd of the frame-to-frame brightness and its one-step memory."""
    S = stack.astype(np.float64) - offset_map[None]
    n, ny, nx = S.shape
    y0, x0 = max(ny // 2 - roi // 2, 0), max(nx // 2 - roi // 2, 0)
    mt = S[:, y0:y0 + roi, x0:x0 + roi].mean(axis=(1, 2))
    eps = mt / mt.mean() - 1
    d = eps - eps.mean()
    return {"flicker_cv": float(np.sqrt(d.var(ddof=1))), "flicker_memory": float((d[:-1] * d[1:]).sum() / (d ** 2).sum()), "n_frames": int(n)}


# ------------------------------------------------------------------------------------------------------------------ flat field
def fit_flat(stack: np.ndarray, offset_map: np.ndarray, temporal_var_adu2: float) -> dict:
    """stack: (n, ny, nx) of the brightest non-saturated flat series. The normalised flat is split into a smooth illumination pattern
    (constant + tilt + quadratic vignetting, fitted) and what is left: the pixel-to-pixel gain pattern (PRNU)."""
    n, ny, nx = stack.shape
    S = stack.astype(np.float64)
    mt = (S - offset_map[None]).mean(axis=(1, 2))
    S = (S - offset_map[None]) / (mt / mt.mean())[:, None, None]              # divide out the flicker
    mean, _, _ = _robust_moments(S)
    M = float(mean.mean())
    N = mean / M                                                                # illumination x pixel gain, mean 1
    yy, xx = np.mgrid[:ny, :nx]
    u, v = (xx - (nx - 1) / 2) / (nx / 2), (yy - (ny - 1) / 2) / (ny / 2)
    A = np.stack([np.ones(N.size), u.ravel(), v.ravel(), (u ** 2 + v ** 2).ravel()], axis=1)[::7]
    coef = np.linalg.lstsq(A, N.ravel()[::7], rcond=None)[0]
    L = coef[0] + coef[1] * u + coef[2] * v + coef[3] * (u ** 2 + v ** 2)
    prnu = N / L
    noise2 = temporal_var_adu2 / (n * M * M)                                    # noise of N itself (relative variance)
    prnu_cv = float(np.sqrt(max(np.var(prnu) - noise2, 0.0)))
    return {"flat_norm": N.astype(np.float32), "illumination_map": L.astype(np.float32), "prnu_cv": prnu_cv,
            "vignette": float(-2 * coef[3] / coef[0]), "tilt": float(np.hypot(coef[1], coef[2]) / coef[0]),
            "mean_signal_adu": M, "flat_noise_cv": float(np.sqrt(noise2))}


# ------------------------------------------------------------------------------------------------------------------ beads
def _pixint(p, c, s):
    r = np.sqrt(2) * s
    return 0.5 * (erf((p + 0.5 - c) / r) - erf((p - 0.5 - c) / r))


def _model2d(coords, bg, A, y0, x0, s):
    ys, xs = coords
    return bg + A * np.outer(_pixint(ys, y0, s), _pixint(xs, x0, s)).ravel()


def fit_gauss2d(patch: np.ndarray, y0: int, x0: int):
    """Fit background + a pixel-integrated 2D Gaussian to a small patch whose top-left pixel is (y0, x0). Returns (y, x, sigma, amplitude) or None."""
    h, w = patch.shape
    ys, xs = np.arange(y0, y0 + h), np.arange(x0, x0 + w)
    bg0 = float(np.median(np.concatenate([patch[0], patch[-1], patch[:, 0], patch[:, -1]])))
    A0 = float(max((patch - bg0).sum(), 1.0))
    iy, ix = np.unravel_index(np.argmax(patch), patch.shape)
    try:
        p, _ = curve_fit(_model2d, (ys, xs), patch.ravel(), p0=[bg0, A0, ys[iy], xs[ix], 1.0],
                         bounds=([-np.inf, 0, ys[0], xs[0], 0.3], [np.inf, np.inf, ys[-1], xs[-1], 6.0]), maxfev=400)
    except Exception:
        return None
    return float(p[2]), float(p[3]), float(p[4]), float(p[1])


def _gauss1d(k, bg, A, k0, s):
    return bg + A * np.exp(-0.5 * ((k - k0) / s) ** 2)


def fit_gauss1d(profile: np.ndarray):
    """Fit background + Gaussian to a z-profile (index units). Returns (centre index, sigma in indices, amplitude) or None."""
    k = np.arange(len(profile), dtype=float)
    bg0 = float(np.percentile(profile, 10))
    k0 = float(np.argmax(profile))
    try:
        p, _ = curve_fit(_gauss1d, k, profile, p0=[bg0, float(profile.max() - bg0), k0, 1.0],
                         bounds=([-np.inf, 0, -2, 0.2], [np.inf, np.inf, len(k) + 1, 12.0]), maxfev=400)
    except Exception:
        return None
    return float(p[2]), float(p[3]), float(p[1])


def detect_beads(img: np.ndarray, sigma: float = 1.5, nsig: float = 10.0, margin: int = 10, min_sep: int = 14) -> np.ndarray:
    """Isolated bright spots of a 2D image -> (n, 2) array of (y, x) pixel positions."""
    sm = gaussian_filter(img.astype(float), sigma)
    bg, mad = np.median(sm), _robust_sd(sm) + 1e-9
    peaks = (sm == maximum_filter(sm, size=2 * min_sep // 2 + 1)) & (sm > bg + nsig * mad)
    ys, xs = np.nonzero(peaks)
    ny, nx = img.shape
    keep = (ys >= margin) & (ys < ny - margin) & (xs >= margin) & (xs < nx - margin)
    ys, xs = ys[keep], xs[keep]
    order = np.argsort(-sm[ys, xs])
    out = []
    for i in order:
        if all(np.hypot(ys[i] - y, xs[i] - x) >= min_sep for y, x in out):
            out.append((ys[i], xs[i]))
    return np.array(out, float).reshape(-1, 2)


def _fit_bead(frame, cy: float, cx: float, z_um, half: int = 6):
    """One bead in one frame of one channel. frame: (ny, nx) or (nz, ny, nx). Returns (y, x, sigma_xy, z or nan, sigma_z or nan) or None."""
    three_d = frame.ndim == 3
    ny, nx = frame.shape[-2:]
    iy, ix = int(round(cy)), int(round(cx))
    y0, y1, x0, x1 = max(iy - half, 0), min(iy + half + 1, ny), max(ix - half, 0), min(ix + half + 1, nx)
    if y1 - y0 < 2 * half + 1 or x1 - x0 < 2 * half + 1:
        return None
    if not three_d:
        r = fit_gauss2d(frame[y0:y1, x0:x1].astype(float), y0, x0)
        return None if r is None else (r[0], r[1], r[2], np.nan, np.nan)
    vol = frame[:, y0:y1, x0:x1].astype(float)
    prof = vol[:, half - 1:half + 2, half - 1:half + 2].mean(axis=(1, 2))
    g = fit_gauss1d(prof)
    if g is None:
        return None
    kc = int(np.clip(round(g[0]), 0, len(prof) - 1))
    lat = vol[max(kc - 1, 0):kc + 2].sum(axis=0)
    r = fit_gauss2d(lat, y0, x0)
    if r is None:
        return None
    step = float(z_um[1] - z_um[0]) if len(z_um) > 1 else 1.0
    zc = float(np.interp(g[0], np.arange(len(z_um)), z_um)) if len(z_um) > 1 else float(z_um[0])    # plane heights need not be equally spaced
    return (r[0], r[1], r[2], zc, float(g[1] * step))


def fit_beads(stacks: dict, px_um: float, z_um=None, ref: str = "nucleus", half: int = 6) -> dict:
    """stacks: {channel: (T, ny, nx) or (T, nz, ny, nx)} images of the same beads (time lapse). The reference channel locates the beads; every
    channel is then fitted at the same beads in every frame. Returns per-channel shift relative to `ref`, PSF width, and drift / focus series."""
    chans = list(stacks)
    three_d = stacks[ref].ndim == 4
    T = stacks[ref].shape[0]

    def flat2d(a):
        return a.max(axis=0) if three_d else a

    P0 = detect_beads(flat2d(stacks[ref][0]))
    n = len(P0)
    res = {c: np.full((T, n, 5), np.nan) for c in chans}                      # y, x, sigma_xy, z, sigma_z
    disp = np.zeros(2)
    ref_pos = np.zeros((T, n, 2))
    for t in range(T):
        for i in range(n):
            f = _fit_bead(stacks[ref][t], P0[i, 0] + disp[0], P0[i, 1] + disp[1], z_um, half)
            if f is not None:
                res[ref][t, i] = f
        ref_pos[t] = res[ref][t, :, :2]
        d = np.nanmedian(ref_pos[t] - ref_pos[0], axis=0) if t > 0 else np.zeros(2)
        disp = np.where(np.isfinite(d), d, disp)
    for c in chans:
        if c == ref:
            continue
        for t in range(T):
            for i in range(n):
                if not np.isfinite(ref_pos[t, i, 0]):
                    continue
                f = _fit_bead(stacks[c][t], ref_pos[t, i, 0], ref_pos[t, i, 1], z_um, half)
                if f is not None:
                    res[c][t, i] = f
    out = {"n_beads": int(n), "channels": chans, "ref": ref, "shift_px": {}, "shift_px_sd": {}, "shift_z_um": {}, "psf_sigma_px": {},
           "psf_sigma_z_um": {}, "n_fits": {}}
    for c in chans:
        r = res[c]
        d = (r[:, :, :2] - res[ref][:, :, :2]).reshape(-1, 2)
        dz = (r[:, :, 3] - res[ref][:, :, 3]).ravel()
        ok = np.isfinite(d).all(axis=1)
        out["n_fits"][c] = int(ok.sum())
        out["shift_px"][c] = [float(np.nanmedian(d[:, 0])), float(np.nanmedian(d[:, 1]))] if c != ref else [0.0, 0.0]
        out["shift_px_sd"][c] = float(_robust_sd(d[ok]) / np.sqrt(max(ok.sum(), 1))) if c != ref and ok.any() else 0.0
        out["shift_z_um"][c] = float(np.nanmedian(dz)) if (three_d and c != ref) else 0.0
        out["psf_sigma_px"][c] = float(np.nanmedian(r[:, :, 2]))
        out["psf_sigma_z_um"][c] = float(np.nanmedian(r[:, :, 4])) if three_d else float("nan")
    drift = np.nanmedian(res[ref][:, :, :2] - res[ref][:1, :, :2], axis=1)       # (T, 2) displacement relative to frame 0
    out["drift_px_series"] = drift.tolist()
    if three_d:
        out["drift_z_um_series"] = (np.nanmedian(res[ref][:, :, 3] - res[ref][:1, :, 3], axis=1)).tolist()
    sig_t = np.nanmedian(res[ref][:, :, 2], axis=1)
    out["psf_sigma_px_series"] = sig_t.tolist()
    out["bead_positions_px"] = P0.tolist()
    out["residual_px"] = {c: np.nanmedian(res[c][:, :, :2] - res[ref][:, :, :2], axis=0).tolist() for c in chans if c != ref}
    return out
