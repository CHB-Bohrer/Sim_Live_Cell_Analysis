"""Stage 6: the uneven illumination (flat-field) and the per-cell photobleaching, estimated from the movie itself (no flat-field image).

Idea (self-calibration, as astronomers calibrate a camera on the stars they observe): the cells move, so the SAME nucleus is seen at many
places in the field. Its background-subtracted brightness is

    m(c, t) = a_c x F(x_c(t)) x b_c(t),        b_c(t) = exp( -k x D_c(t) ),   D_c(t) = sum over its earlier frames of F(x_c(s))

a_c = the cell's own brightness (unknown, free per cell), F = the illumination (a smooth 2D polynomial, mean 1 over the field), b_c = the
cell's bleaching, driven by the light IT received along its own path (so it depends on where the cell was and for how long; a cell that
walks into the field from outside arrives unbleached). Taking logs, log m = log a_c + log F(x) - k D, and the per-cell offset log a_c is
removed exactly by subtracting each cell's mean, so only F's 5 coefficients and k = dt / tau are fitted (robust least squares).

Why it is identifiable: a change of brightness that follows the cell's POSITION is the illumination; a decline that follows the light it has
accumulated is bleaching. A cell resting in a bright spot bleaches faster than one in the dim corner, which separates the two.

What it cannot know: light a cell received before it was first seen (outside the camera's field but inside the lit area, or before the
movie). That only shifts the cell's log a_c, so its corrected trace is relative to its first frame in view; it is absolute if the light only
falls on the imaged field (the usual case for a field stop / laser illumination of the camera's area).

The same model is fitted on spot intensities (locus colours, each with its own illumination), with spot positions for F and the cell's track
for the dose. The standard alternative without cells that move is BaSiC (Peng et al. 2017, Nat. Commun.), which estimates the flat-field from
the statistics of many images; it needs many images with different content at each pixel.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_dilation, find_objects


def _uv(y, x, shape):
    ny, nx = shape
    return (np.asarray(y, float) - (ny - 1) / 2) / (ny / 2), (np.asarray(x, float) - (nx - 1) / 2) / (nx / 2)


def _basis(y, x, shape) -> np.ndarray:
    v, u = _uv(y, x, shape)
    return np.stack([u, v, u * u, u * v, v * v], -1)


def _field_mean_basis(shape) -> np.ndarray:
    yy, xx = np.mgrid[:shape[0]:4, :shape[1]:4]
    return _basis(yy.ravel(), xx.ravel(), shape).mean(0)


@dataclass
class FlatFit:
    """Estimated illumination F (mean 1 over the field) and bleaching rate. theta = coefficients of u, v, u^2, uv, v^2."""
    shape: tuple
    theta: np.ndarray = field(default_factory=lambda: np.zeros(5))
    k_per_frame: float = 0.0               # bleaching per frame at the field's mean illumination
    dt_s: float = 1.0
    n_obs: int = 0
    rms_residual: float = np.nan
    fitted_flat: bool = True
    flat_rejected: bool = False            # an illumination term was fitted but did not predict held-out cells better than even light
    cv: dict = field(default_factory=dict)

    @property
    def tau_s(self) -> float:
        return float(self.dt_s / self.k_per_frame) if self.k_per_frame > 0 else np.inf

    def F(self, y, x) -> np.ndarray:
        raw = 1.0 + _basis(y, x, self.shape) @ self.theta
        return np.clip(raw / (1.0 + _field_mean_basis(self.shape) @ self.theta), 0.05, None)

    def map(self) -> np.ndarray:
        yy, xx = np.mgrid[:self.shape[0], :self.shape[1]]
        return self.F(yy, xx).astype(np.float32)

    def dose(self, track: pd.DataFrame) -> pd.Series:
        """Light received before each frame of each cell's track (columns cell_id, t, y_px, x_px), in frames at mean illumination."""
        return _dose(track, self.F(track.y_px.to_numpy(), track.x_px.to_numpy()))

    def bleach(self, track: pd.DataFrame) -> pd.Series:
        return np.exp(-self.k_per_frame * self.dose(track))

    def summary(self) -> dict:
        m = self.map()
        return {"tau_s": self.tau_s, "theta": [float(v) for v in self.theta], "n_obs": int(self.n_obs), "rms_residual": float(self.rms_residual),
                "flat_min": float(m.min()), "flat_max": float(m.max()), "fitted_flat": self.fitted_flat,
                "flat_rejected": self.flat_rejected, "cross_validation": self.cv}


def _dose(track: pd.DataFrame, f: np.ndarray) -> pd.Series:
    """Exclusive cumulative sum of the illumination over each cell's earlier frames (frames where the cell was not seen give no light)."""
    d = pd.DataFrame({"cell_id": track.cell_id.to_numpy(), "t": track.t.to_numpy(), "f": f}, index=track.index).sort_values(["cell_id", "t"])
    d["D"] = d.groupby("cell_id").f.cumsum() - d.f
    return d.D.reindex(track.index)


def nuclear_signals(images, labels, ring_px: tuple = (4, 14), border_px: int = 2) -> pd.DataFrame:
    """Per cell and frame: background-subtracted mean nuclear brightness (mean inside the mask minus the median of a ring 4-14 px around it,
    other nuclei excluded), centroid and whether the nucleus touches the image edge. One frame in memory at a time."""
    img, lab = tifffile.memmap(str(images), mode="r"), tifffile.memmap(str(labels), mode="r")
    if img.ndim == 4:                          # (T, Z, Y, X) stacks: use their max projection like the rest of the 2D analysis
        img = img.max(axis=1)
    ny, nx = lab.shape[-2:]
    r0, r1 = ring_px
    st0, st1 = np.ones((2 * r0 + 1,) * 2, bool), np.ones((2 * r1 + 1,) * 2, bool)
    rows = []
    for t in range(len(lab)):
        L = np.asarray(lab[t])
        f = np.asarray(img[t], np.float32)
        for i, sl in enumerate(find_objects(L)):
            if sl is None:
                continue
            cid = i + 1
            y0, y1 = max(sl[0].start - r1, 0), min(sl[0].stop + r1, ny)
            x0, x1 = max(sl[1].start - r1, 0), min(sl[1].stop + r1, nx)
            sub, fs = L[y0:y1, x0:x1], f[y0:y1, x0:x1]
            m = sub == cid
            ring = binary_dilation(m, st1) & ~binary_dilation(m, st0) & (sub == 0)
            if ring.sum() < 20:
                continue
            yy, xx = np.nonzero(m)
            edge = sl[0].start < border_px or sl[1].start < border_px or sl[0].stop > ny - border_px or sl[1].stop > nx - border_px
            rows.append((t, cid, float(fs[m].mean() - np.median(fs[ring])), y0 + yy.mean(), x0 + xx.mean(), int(m.sum()), bool(edge)))
    return pd.DataFrame(rows, columns=["t", "cell_id", "signal", "y_px", "x_px", "area_px", "touches_border"])


def fit_flat_bleach(values: pd.DataFrame, track: pd.DataFrame, shape, dt_s: float, fit_flat: bool = True, dose_flat: FlatFit | None = None,
                    prior_sd: float = 0.1, min_frames: int = 5, min_gain: float = 0.02) -> FlatFit:
    """Fit log value = a_c + log F(position) - k D_c(t).

    values: cell_id, t, y_px, x_px, value (> 0 used; the position where F applies: the nucleus' centroid or the spot)
    track:  cell_id, t, y_px, x_px (every frame the cell was seen: the dose integrates F along it)
    fit_flat=False: F = 1 in the intensity term (e.g. intensities already flat-field corrected by a calibration); the dose then uses
    `dose_flat` (another fit's illumination) or a uniform field.
    Two passes: (1) bleaching only, which gives the noise of the log values (robust sd); (2) the full fit with a robust loss scaled to that
    noise and a Gaussian prior of sd `prior_sd` on each illumination coefficient, so that data which say little about the illumination
    (cells that barely move, spot intensities that fluctuate with probe binding) leave F close to even instead of fitting noise."""
    if fit_flat:
        first = fit_flat_bleach(values, track, shape, dt_s, fit_flat=False, min_frames=min_frames)
        if not np.isfinite(first.rms_residual):
            return first
        full = _fit(values, track, shape, dt_s, True, None, first.rms_residual, prior_sd, min_frames)
        full.cv = _cross_validate(values, track, shape, dt_s, first.rms_residual, prior_sd, min_frames)
        if not full.cv.get("gain", 0.0) >= min_gain:   # the illumination term does not predict held-out cells clearly better than even light
            first.cv, first.fitted_flat = full.cv, False
            first.flat_rejected = True
            return first
        return full
    return _fit(values, track, shape, dt_s, False, dose_flat, None, prior_sd, min_frames)


def _held_out_loss(fit: FlatFit, values: pd.DataFrame, track: pd.DataFrame) -> float:
    """Robust loss of the log values of cells NOT used in the fit; each cell's own brightness is still taken out (its mean)."""
    v = values[np.isfinite(values.value) & (values.value > 0)]
    tr = track.drop_duplicates(["cell_id", "t"]).reset_index(drop=True)
    b = pd.Series(fit.bleach(tr).to_numpy(), index=pd.MultiIndex.from_frame(tr[["cell_id", "t"]]))
    bb = b.reindex(pd.MultiIndex.from_arrays([v.cell_id, v.t])).to_numpy()
    ok = np.isfinite(bb)
    f = fit.F(v.y_px.to_numpy(), v.x_px.to_numpy()) if fit.fitted_flat else 1.0
    r = pd.Series(np.log(v.value.to_numpy() / (f * np.where(ok, bb, 1.0))), index=v.index)[ok]
    r = r - r.groupby(v.cell_id[ok]).transform("mean")
    return float(np.median(np.abs(r))) if len(r) else np.nan


def _cross_validate(values, track, shape, dt_s, noise, prior_sd, min_frames, k: int = 4) -> dict:
    """k-fold over cells: does the illumination term lower the held-out cells' robust loss? gain = relative reduction (> 0 = it helps)."""
    cells = np.unique(values.cell_id)
    if len(cells) < k:
        return {"gain": np.nan, "n_cells": int(len(cells))}
    rng = np.random.default_rng(0)
    folds = np.array_split(rng.permutation(cells), k)
    with_f, without = [], []
    for test in folds:
        tr_v, te_v = values[~values.cell_id.isin(test)], values[values.cell_id.isin(test)]
        a = _fit(tr_v, track, shape, dt_s, True, None, noise, prior_sd, min_frames)
        b = _fit(tr_v, track, shape, dt_s, False, None, None, prior_sd, min_frames)
        with_f.append(_held_out_loss(a, te_v, track))
        without.append(_held_out_loss(b, te_v, track))
    w, wo = float(np.nanmean(with_f)), float(np.nanmean(without))
    return {"gain": (wo - w) / wo if wo > 0 else np.nan, "loss_with_flat": w, "loss_even": wo, "n_cells": int(len(cells))}


def _fit(values, track, shape, dt_s, fit_flat, dose_flat, noise, prior_sd, min_frames) -> FlatFit:
    from scipy.optimize import least_squares

    v = values[np.isfinite(values.value) & (values.value > 0)].copy()
    n_per = v.groupby("cell_id").t.transform("size")
    v = v[n_per >= min_frames]
    out = FlatFit(tuple(shape), dt_s=float(dt_s), fitted_flat=fit_flat)
    if len(v) < 20:
        return out
    tr = track.drop_duplicates(["cell_id", "t"]).set_index(["cell_id", "t"]).sort_index()
    key = pd.MultiIndex.from_arrays([v.cell_id, v.t])
    y = np.log(v.value.to_numpy())
    cell_codes, inv = np.unique(v.cell_id.to_numpy(), return_inverse=True)
    counts = np.bincount(inv)
    Bv = _basis(v.y_px.to_numpy(), v.x_px.to_numpy(), shape)
    trk = tr.reset_index()
    Bt = _basis(trk.y_px.to_numpy(), trk.x_px.to_numpy(), shape)
    mb = _field_mean_basis(shape)
    order = trk.sort_values(["cell_id", "t"]).index.to_numpy()
    cid_sorted = trk.cell_id.to_numpy()[order]
    first = np.r_[True, cid_sorted[1:] != cid_sorted[:-1]]
    pos_in_track = pd.Series(np.arange(len(trk)), index=pd.MultiIndex.from_frame(trk[["cell_id", "t"]]))
    idx_v = pos_in_track.reindex(key).to_numpy()
    ok = np.isfinite(idx_v)
    y, inv, Bv, idx_v = y[ok], inv[ok], Bv[ok], idx_v[ok].astype(int)
    counts = np.bincount(inv, minlength=len(cell_codes))

    def dose_of(theta_t):
        f = np.clip((1 + Bt @ theta_t) / (1 + mb @ theta_t), 0.05, None)
        fs = f[order]
        cs = np.cumsum(fs)
        start = np.maximum.accumulate(np.where(first, np.arange(len(fs)), 0))
        excl = cs - fs - np.where(start > 0, cs[start - 1], 0.0)     # light received earlier in the same cell's track
        D = np.empty(len(fs))
        D[order] = excl
        return D

    fixed_D = None
    if not fit_flat:
        fixed_D = dose_of(dose_flat.theta) if dose_flat is not None else dose_of(np.zeros(5))

    def data_resid(p):
        th, k = (p[:5], p[5]) if fit_flat else (np.zeros(5), p[0])
        logF = np.log(np.clip((1 + Bv @ th) / (1 + mb @ th), 0.05, None)) if fit_flat else 0.0
        D = dose_of(th) if fit_flat else fixed_D
        r = y - logF + k * D[idx_v]
        return r - (np.bincount(inv, r, minlength=len(cell_codes)) / np.maximum(counts, 1))[inv]

    sd = float(noise) if noise else 0.05
    # prior: each coefficient ~ N(0, prior_sd); in the units of the data term (residual / noise), one entry per coefficient
    resid = (lambda p: np.r_[data_resid(p) / sd, p[:5] / prior_sd]) if fit_flat else (lambda p: data_resid(p) / sd)
    p0 = np.r_[np.zeros(5), 1e-4] if fit_flat else np.array([1e-4])
    lo = np.r_[np.full(5, -0.9), 0.0] if fit_flat else np.array([0.0])
    hi = np.r_[np.full(5, 0.9), 0.5] if fit_flat else np.array([0.5])
    res = least_squares(resid, p0, bounds=(lo, hi), loss="soft_l1", f_scale=2.0, x_scale="jac")
    if fit_flat:
        out.theta, out.k_per_frame = res.x[:5].copy(), float(res.x[5])
    else:
        out.theta, out.k_per_frame = (dose_flat.theta.copy() if dose_flat is not None else np.zeros(5)), float(res.x[0])
        out.fitted_flat = False
    r = data_resid(res.x)
    out.n_obs, out.rms_residual = int(len(y)), float(1.4826 * np.median(np.abs(r - np.median(r))))     # robust sd of the log values
    return out


def correct_values(df: pd.DataFrame, fit: FlatFit, track: pd.DataFrame, col: str = "intensity", divide_flat: bool = True) -> pd.Series:
    """Background-subtracted intensity / (F(position) x b_cell(t)): F at the measured spot position (lab_y_px, lab_x_px), the bleaching from
    the cell's own track. divide_flat=False when the intensities were already flat-field corrected (calibration applied)."""
    tr = track.drop_duplicates(["cell_id", "t"]).reset_index(drop=True)
    b = pd.Series(fit.bleach(tr).to_numpy(), index=pd.MultiIndex.from_frame(tr[["cell_id", "t"]]))
    bb = b.reindex(pd.MultiIndex.from_arrays([df.cell_id, df.t])).to_numpy()
    bb = np.where(np.isfinite(bb), bb, 1.0)
    f = fit.F(df.lab_y_px.to_numpy(), df.lab_x_px.to_numpy()) if divide_flat else 1.0
    return df[col] / (f * bb)
