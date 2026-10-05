"""Stage 3 (nuclear channel only, for now): render cells.csv into an image movie + ground-truth label masks.

Forward model: uniform NLS-GFP nucleus (per-cell brightness) -> PSF blur (Gaussian, sigma from NA/wavelength)
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


def render_nuclei(cells: pd.DataFrame, cfg: dict, rng: np.random.Generator):
    opt, acq = cfg["optics"], cfg["acquisition"]
    px_um = opt["pixel_size_nm"] / 1000.0
    ny, nx = (int(round(v / px_um)) for v in cfg["geometry"]["fov_um"])
    T = acq["n_frames"]
    yy, xx = np.mgrid[:ny, :nx]
    sig = psf_sigma_px(opt)

    ids = cells["cell_id"].unique()
    brightness = dict(zip(ids, rng.lognormal(0.0, opt.get("brightness_cv", 0.3), len(ids))))
    flux = opt["photons_per_px_s"] * acq["exposure_s"]
    bleach = np.exp(-acq["frame_interval_s"] * np.arange(T) / opt["bleach_tau_s"]) \
        if opt.get("bleach_tau_s") else np.ones(T)

    img = np.zeros((T, ny, nx), np.float32)
    lab = np.zeros((T, ny, nx), np.uint16)
    for t, grp in cells.groupby("t"):
        clean = np.zeros((ny, nx), np.float32)
        for r in grp.itertuples():
            cy, cx, rad = r.y_um / px_um, r.x_um / px_um, r.radius_um / px_um
            inside = (yy - cy) ** 2 + (xx - cx) ** 2 <= rad**2
            lab[t][inside] = r.cell_id
            clean[inside] = brightness[r.cell_id]
        blurred = gaussian_filter(clean, sig) if sig > 0 else clean
        photons = rng.poisson((blurred * flux * bleach[t] + opt["background_photons"]).clip(0))
        electrons = photons * opt.get("quantum_efficiency", 0.8)
        img[t] = (electrons + rng.normal(0, opt["read_noise_e"], electrons.shape)) * opt.get("gain_adu_per_e", 1.0) \
            + opt.get("offset_adu", 100.0)
    return img.clip(0, 65535).astype(np.float32), lab
