"""Stage 3 (nuclear channel only, for now): render cells.csv into an image movie + ground-truth label masks.

Forward model: NLS-GFP nucleus (per-cell brightness x per-nucleus texture: chromatin-like variation + dark
nucleoli, fixed to the nucleus) -> PSF blur (Gaussian, sigma from NA/wavelength)
-> photons (exposure x flux, optional photobleaching) -> Poisson shot noise -> camera gain/read noise/offset.
Ground truth: integer label mask per frame, value = persistent cell ID (0 = background).
Locus channel(s) are added when stage 1 is wired in.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


def psf_sigma_px(optics: dict) -> float:
    """Gaussian approximation to the widefield Airy PSF: sigma ~ 0.21 * lambda / NA."""
    return 0.21 * optics["wavelength_nm"] / optics["NA"] / optics["pixel_size_nm"]


def _nucleus_mask(dy, dx, r, px_um, harm):
    """Pixels inside a deformed ellipse. dy, dx: pixel offsets from the cell centre; r: a cells.csv row."""
    from simlive.stage2_cells.motion import boundary_rho

    c, s = np.cos(r.angle_rad), np.sin(r.angle_rad)
    u, v = dy * s + dx * c, dy * c - dx * s        # rotate into the ellipse frame (u along the long axis)
    a_ax, b_ax = np.sqrt(r.aspect), 1 / np.sqrt(r.aspect)
    uu, vv = u / a_ax, v / b_ax
    amp = np.array([getattr(r, f"amp{int(k)}") for k in harm])
    phase = np.array([getattr(r, f"phase{int(k)}") for k in harm])
    r_px = r.radius_um / px_um
    inside = np.hypot(uu, vv) <= r_px * boundary_rho(np.arctan2(vv, uu), amp, phase, harm)
    return inside, uu / r_px, vv / r_px  # also the position in the nucleus' own frame (units of its radius)


_TEX_GRID, _TEX_EXTENT = 96, 1.5


def make_texture(tex: dict, radius_um: float, rng: np.random.Generator) -> np.ndarray:
    """Per-nucleus intensity texture on a grid covering +-1.5 nuclear radii, mean ~1 inside the nucleus.

    Fixed to the nucleus (rotates/deforms with it), so it is a feature a tracker can use. Components:
    spatially correlated chromatin-like variation (`contrast`, `correlation_um`) and dark nucleoli where nuclear
    NLS-GFP is excluded (`n_nucleoli` on average, `nucleolus_depth`).
    """
    g = _TEX_GRID
    cell = 2 * _TEX_EXTENT / g                          # grid spacing in nuclear radii
    field = gaussian_filter(rng.normal(size=(g, g)), (tex["correlation_um"] / radius_um) / cell, mode="wrap")
    field /= field.std()
    t = np.clip(1 + tex["contrast"] * field, 0.25, None)
    ax = np.linspace(-_TEX_EXTENT, _TEX_EXTENT, g)
    gx, gy = np.meshgrid(ax, ax)
    for _ in range(min(int(rng.poisson(tex.get("n_nucleoli", 0))), 6)):
        ang, rad, s = rng.uniform(0, 2 * np.pi), 0.65 * np.sqrt(rng.uniform()), rng.uniform(0.10, 0.18)
        d2 = (gx - rad * np.cos(ang)) ** 2 + (gy - rad * np.sin(ang)) ** 2
        t *= 1 - tex.get("nucleolus_depth", 0.5) * np.exp(-d2 / (2 * s**2))
    return t / t[np.hypot(gx, gy) <= 1].mean()


def _sample_texture(tex_map: np.ndarray, xn: np.ndarray, yn: np.ndarray) -> np.ndarray:
    from scipy.ndimage import map_coordinates

    to_grid = lambda c: (c + _TEX_EXTENT) / (2 * _TEX_EXTENT) * (_TEX_GRID - 1)
    return map_coordinates(tex_map, [to_grid(yn), to_grid(xn)], order=1, mode="nearest")


def render_loci(loci: pd.DataFrame, cfg: dict, rng: np.random.Generator, locus_id: int = 0) -> np.ndarray:
    """One locus colour channel: diffraction-limited spots (Gaussian PSF integrated over each pixel).

    Each locus has its own colour/channel (`loci.wavelengths_nm[locus_id]`), so only rows of that locus are drawn.
    Photons per spot = loci.photons_per_locus_s x exposure x photobleaching; then Poisson shot noise and camera
    gain/read noise/offset, like the nuclear channel. Pixel centres sit at integer coordinates (pixel i spans i+-0.5).
    """
    from scipy.special import erf

    opt, acq, lc = cfg["optics"], cfg["acquisition"], cfg["loci"]
    loci = loci[loci.locus_id == locus_id]
    px_um = opt["pixel_size_nm"] / 1000.0
    ny, nx = (int(round(v / px_um)) for v in cfg["geometry"]["fov_um"])
    T = acq["n_frames"]
    sig = psf_sigma_px({**opt, "wavelength_nm": lc["wavelengths_nm"][locus_id]})
    sig = max(sig, 0.3)
    bleach_tau = lc.get("bleach_tau_s", opt.get("bleach_tau_s"))
    bleach = np.exp(-acq["frame_interval_s"] * np.arange(T) / bleach_tau) if bleach_tau else np.ones(T)
    photons_per_spot = lc["photons_per_locus_s"] * acq["exposure_s"]
    half = int(np.ceil(4 * sig)) + 1
    s2 = np.sqrt(2) * sig

    img = np.zeros((T, ny, nx), np.float32)
    for t, g in loci.groupby("t"):
        clean = np.zeros((ny, nx), np.float32)
        for r in g.itertuples():
            cy, cx = r.y_um / px_um, r.x_um / px_um
            y0, y1 = max(int(round(cy)) - half, 0), min(int(round(cy)) + half + 1, ny)
            x0, x1 = max(int(round(cx)) - half, 0), min(int(round(cx)) + half + 1, nx)
            if y1 <= y0 or x1 <= x0:
                continue
            ys, xs = np.arange(y0, y1), np.arange(x0, x1)
            py = 0.5 * (erf((ys + 0.5 - cy) / s2) - erf((ys - 0.5 - cy) / s2))
            px = 0.5 * (erf((xs + 0.5 - cx) / s2) - erf((xs - 0.5 - cx) / s2))
            clean[y0:y1, x0:x1] += photons_per_spot * bleach[t] * np.outer(py, px)
        photons = rng.poisson(clean + opt["background_photons"])
        electrons = photons * opt.get("quantum_efficiency", 0.8)
        img[t] = (electrons + rng.normal(0, opt["read_noise_e"], electrons.shape)) * opt.get("gain_adu_per_e", 1.0) \
            + opt.get("offset_adu", 100.0)
    return img.clip(0, 65535).astype(np.float32)


def render_nuclei(cells: pd.DataFrame, cfg: dict, rng: np.random.Generator):
    opt, acq = cfg["optics"], cfg["acquisition"]
    px_um = opt["pixel_size_nm"] / 1000.0
    ny, nx = (int(round(v / px_um)) for v in cfg["geometry"]["fov_um"])
    T = acq["n_frames"]
    yy, xx = np.mgrid[:ny, :nx]
    sig = psf_sigma_px(opt)

    harm = np.array([int(c[3:]) for c in cells.columns if c.startswith("amp")], float)
    ids = cells["cell_id"].unique()
    brightness = dict(zip(ids, rng.lognormal(0.0, opt.get("brightness_cv", 0.3), len(ids))))
    tex_cfg = cfg.get("nucleus_texture")  # None/absent or contrast 0 and no nucleoli -> uniform nuclei
    if tex_cfg and tex_cfg.get("contrast", 0) == 0 and tex_cfg.get("n_nucleoli", 0) == 0:
        tex_cfg = None
    tex_rng = np.random.default_rng(rng.integers(2**63))  # separate stream so texture draws don't shift the noise
    textures: dict = {}
    flux = opt["photons_per_px_s"] * acq["exposure_s"]
    bleach = np.exp(-acq["frame_interval_s"] * np.arange(T) / opt["bleach_tau_s"]) \
        if opt.get("bleach_tau_s") else np.ones(T)

    img = np.zeros((T, ny, nx), np.float32)
    lab = np.zeros((T, ny, nx), np.uint16)
    for t, grp in cells.groupby("t"):
        clean = np.zeros((ny, nx), np.float32)
        for r in grp.itertuples():
            cy, cx = r.y_um / px_um, r.x_um / px_um
            half = int(np.ceil(r.bound_radius_um / px_um)) + 2
            y0, y1, x0, x1 = max(int(cy) - half, 0), min(int(cy) + half + 1, ny), max(int(cx) - half, 0), min(int(cx) + half + 1, nx)
            inside, xn, yn = _nucleus_mask(yy[y0:y1, x0:x1] - cy, xx[y0:y1, x0:x1] - cx, r, px_um, harm)
            lab[t, y0:y1, x0:x1][inside] = r.cell_id
            val = brightness[r.cell_id] * np.ones(inside.shape, np.float32)
            if tex_cfg:
                if r.cell_id not in textures:
                    textures[r.cell_id] = make_texture(tex_cfg, r.radius_um, tex_rng)
                val = val * _sample_texture(textures[r.cell_id], xn, yn)
            clean[y0:y1, x0:x1][inside] = val[inside]
        blurred = gaussian_filter(clean, sig) if sig > 0 else clean
        photons = rng.poisson((blurred * flux * bleach[t] + opt["background_photons"]).clip(0))
        electrons = photons * opt.get("quantum_efficiency", 0.8)
        img[t] = (electrons + rng.normal(0, opt["read_noise_e"], electrons.shape)) * opt.get("gain_adu_per_e", 1.0) \
            + opt.get("offset_adu", 100.0)
    return img.clip(0, 65535).astype(np.float32), lab
