"""What a microscope reports about the E-P vector: 2D widefield, 3D z-stack and 3D biplane (docs/EP_CURVE.md).

A vector-level acquisition model (NOT a rendered movie): for every frame the true E-P vector comes from the polymer-driven population
(`ep_polymer`), and the measured vector is that vector plus Gaussian errors whose sizes are computed from the instrument:

  localization (photon-limited): lateral, Thompson et al. 2002: s2 = (sigma^2 + a^2/12) / N + 8 pi sigma^4 b^2 / (a^2 N^2) per plane, planes combined by inverse
      variance; axial, from the photon count of each plane falling as a Gaussian of the distance to the locus (the repo's PSF model, stage3 `axial_sigma_um`):
      Fisher information  I_z = sum_k (N_k dz_k / sigma_z^2)^2 / (N_k + n_eff b^2)  with n_eff = 4 pi sigma^2 / a^2 pixels;
  motion: the loci move while they are recorded. The simulation stores nothing between its saved blocks, so the single-locus MSD below one block is the power law
      with the small-lag exponent of the simulation (`Calibration.msd_single_s`; an EXTRAPOLATION):
      - exposure and plane spread: a locus is seen over a window of effective length T (exposure + the spread of the planes that carry its signal); its measured
        position is the window average: variance MSD_axis(T) / 12 per axis (exact for Brownian motion, approximate for the subdiffusive polymer);
      - the two colours are not recorded at the same moment (sequential excitation, 0.15 s): the error of the difference vector is the displacement of each
        locus between its own time and the common mean time: 2 MSD_axis(|t_E - t_P| / 2) per axis;
  piezo (z-stack only): every plane's height has an independent repeatability error; planes shared by both loci cancel in the difference
      (variance = repeatability^2 sum_k (a_k^E - a_k^P)^2 with a_k the photon weights of the planes).
  detection: a locus is lost if its photons in all planes together fall below `detect_min_photons`.
The error sizes are returned (`Measured.parts`) so the budget can be read off. All photon numbers, the background and the focus scatter are PLACEHOLDERS.
Modalities: "2d" (one plane, no z), "zstack" (n_planes planes step_um apart, plane_time_s per plane for both colours, scanned over the whole nucleus),
"biplane" (two simultaneous planes plane_sep_um apart, half the photons each, focal pair set per cell with focus_sd_um scatter).
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.ndimage import gaussian_filter1d

from simlive.stage3_microscopy.render import axial_sigma_um
from simlive.stage7_validation import ep_curve as ep


@dataclass(frozen=True)
class Acquisition:
    kind: str = "zstack"                    # "2d" | "zstack" | "biplane"
    exposure_s: float = 0.1
    colour_delay_s: float = 0.15            # between the promoter and the enhancer colour (Vutara preset: sequential excitation)
    photons_per_s: float = 6000.0           # per locus (PLACEHOLDER: configs/loci_3d.yaml)
    background_px: float = 5.0              # photons per pixel per exposure (PLACEHOLDER: background_photons)
    read_noise_e: float = 0.7
    pixel_nm: float = 100.0
    NA: float = 1.3
    n_immersion: float = 1.406
    wavelengths_nm: tuple = (575.0, 665.0)  # promoter label (locus 0), enhancer label (locus 1)
    n_planes: int = 17
    step_um: float = 0.8
    plane_time_s: float = 0.30              # per plane: two colours + piezo move
    plane_sep_um: float = 0.6               # biplane (Bruker bi-plane page)
    focus_sd_um: float = 0.3                # 2D / biplane: scatter of the focal plane(s) around the pair's height
    z_spread_um: float = 2.0                # scatter of the pair's height between cells (PLACEHOLDER; nuclei are ~12 um thick)
    piezo_nm: float = 30.0                  # z-stack: repeatability of every plane
    shift_nm: float = 0.0                   # colour offset left after correction (in plane, x)
    detect_min_photons: float = 50.0
    wrong_pick_rate: float = 0.0            # frames where the spot finder picked something else: distance uniform in 0-2000 nm
    photon_scale: float = 1.0               # multiplies photons_per_s (for the photon scan)

    @property
    def dims(self) -> int:
        return 2 if self.kind == "2d" else 3


@dataclass
class Measured:
    d: np.ndarray                           # measured distance (nm) per frame, NaN where a locus is lost
    sigma_xy: float                         # rms per-axis error of the difference vector (nm), for the errors-in-variables kernel
    sigma_z: float
    dims: int
    shift_nm: float
    detected: float                         # fraction of frames with both loci
    parts: dict = field(default_factory=dict)   # median per-axis error contributions (nm): localization, blur, colour_delay, piezo; mean frame offset (s)


def _plane_heights(acq: Acquisition, z_mid: np.ndarray, rng: np.random.Generator, n_cells: int) -> np.ndarray:
    """Heights of the recorded planes (um), (cells, K): z-stack = fixed grid centred on the nominal nucleus height 0; 2D / biplane = per cell, around the pair."""
    if acq.kind == "zstack":
        zk = (np.arange(acq.n_planes) - (acq.n_planes - 1) / 2.0) * acq.step_um
        return np.broadcast_to(zk, (n_cells, acq.n_planes))
    zf = z_mid + rng.normal(0, acq.focus_sd_um, n_cells)
    if acq.kind == "biplane":
        return np.stack([zf - acq.plane_sep_um / 2, zf + acq.plane_sep_um / 2], 1)
    return zf[:, None]


def _locus_stats(acq: Acquisition, calib, z_c: np.ndarray, zk: np.ndarray, tk: np.ndarray, lam_nm: float, colour_offset_s: float) -> dict:
    """Photon weights and errors of one locus: z_c (F,), planes zk (F, K) um, plane times tk (K,)."""
    optics = {"NA": acq.NA, "immersion_n": acq.n_immersion, "wavelength_nm": lam_nm, "modality": "widefield"}
    sz_nm = axial_sigma_um(optics, lam_nm) * 1000.0
    s_nm = 0.21 * lam_nm / acq.NA
    share = 0.5 if acq.kind == "biplane" else 1.0
    N0 = acq.photons_per_s * acq.photon_scale * acq.exposure_s * share
    dz = (zk - z_c[:, None]) * 1000.0
    Nk = N0 * np.exp(-0.5 * (dz / sz_nm) ** 2)                                       # photons in plane k
    b2 = acq.background_px + acq.read_noise_e ** 2
    Nk_safe = np.maximum(Nk, 1e-3)
    inv_s2 = 1.0 / ((s_nm ** 2 + acq.pixel_nm ** 2 / 12.0) / Nk_safe + 8 * np.pi * s_nm ** 4 * b2 / (acq.pixel_nm ** 2 * Nk_safe ** 2))
    sxy2 = 1.0 / inv_s2.sum(1)
    n_eff = 4 * np.pi * s_nm ** 2 / acq.pixel_nm ** 2
    Iz = ((Nk * dz / sz_nm ** 2) ** 2 / (Nk + n_eff * b2)).sum(1)
    sz2 = 1.0 / np.maximum(Iz, 1e-12)
    Ntot = Nk.sum(1)
    a = Nk / np.maximum(Ntot[:, None], 1e-12)
    tbar = (a * tk[None, :]).sum(1) + colour_offset_s
    var_t = (a * (tk[None, :] - (a * tk[None, :]).sum(1, keepdims=True)) ** 2).sum(1) + acq.exposure_s ** 2 / 12.0
    return {"sxy2": sxy2, "sz2": sz2, "N": Ntot, "a": a, "t": tbar, "var_t": var_t}


def measure(pop: ep.Population, blocks: np.ndarray, acq: Acquisition, seed: int = 0) -> Measured:
    """Measured E-P distance of every frame of the population. blocks: (cells, n_blocks, 3) float32, the true E-P vector (nm) at every saved block (the population
    keeps only the frame-time vector; the block trace is needed for the time offsets of a z-stack)."""
    calib = pop.calib
    rng = np.random.default_rng(seed)
    cells = np.unique(pop.cell_id)
    n_cells, T = len(cells), int((pop.cell_id == cells[0]).sum())
    F = n_cells * T
    ft = pop.t_s
    # the pair's height above the nominal stack centre (um): fixed per cell
    z_mid = rng.normal(0, acq.z_spread_um, n_cells)
    zk_cell = _plane_heights(acq, z_mid, rng, n_cells)
    zk = np.repeat(zk_cell, T, axis=0)
    K = zk.shape[1]
    tk = np.arange(K) * acq.plane_time_s if acq.kind == "zstack" else np.zeros(K)
    vz = pop.vec_nm[:, 2] / 1000.0
    zmid = np.repeat(z_mid, T)
    zP, zE = zmid - vz / 2, zmid + vz / 2                                              # promoter (locus 0), enhancer (locus 1)
    P = _locus_stats(acq, calib, zP, zk, tk, acq.wavelengths_nm[0], 0.0)
    E = _locus_stats(acq, calib, zE, zk, tk, acq.wavelengths_nm[1], acq.colour_delay_s)
    msd_axis = lambda tau: calib.msd_single_s(tau) / 3.0                              # noqa: E731  (nm^2, per axis)
    blur = {c: msd_axis(np.sqrt(12.0 * S["var_t"])) / 12.0 for c, S in (("P", P), ("E", E))}
    dt_ep = np.abs(E["t"] - P["t"])
    offs = 2.0 * msd_axis(dt_ep / 2.0)
    var_xy = P["sxy2"] + E["sxy2"] + blur["P"] + blur["E"] + offs
    var_z = P["sz2"] + E["sz2"] + blur["P"] + blur["E"] + offs
    piezo = np.zeros(F)
    if acq.kind == "zstack":
        piezo = acq.piezo_nm ** 2 * ((P["a"] - E["a"]) ** 2).sum(1)
        var_z = var_z + piezo
    # true vector at the mean time of the two measurements: block k holds for [k, k + 1) x block_s, as for the distance that drives the promoter
    # (transcription.simulate_cell); until 2026-10-26 round() took the NEXT block in most frames (a 0.23 s late measurement on average)
    tmean = ft + 0.5 * (P["t"] + E["t"])
    bi = np.minimum(np.floor(tmean / calib.block_s).astype(int), blocks.shape[1] - 1)
    cidx = np.repeat(np.arange(n_cells), T)
    v = blocks[cidx, bi, :].astype(np.float64)
    noise = rng.normal(size=(F, 3)) * np.sqrt(np.stack([var_xy, var_xy, var_z], 1))
    v = v + noise
    v[:, 0] += acq.shift_nm
    d = np.sqrt((v[:, :acq.dims] ** 2).sum(1))
    ok = (P["N"] >= acq.detect_min_photons) & (E["N"] >= acq.detect_min_photons)
    d[~ok] = np.nan
    if acq.wrong_pick_rate > 0:
        bad = rng.random(F) < acq.wrong_pick_rate
        d[bad] = rng.uniform(0, 2000, bad.sum())
    loc_xy = np.sqrt(np.median((P["sxy2"] + E["sxy2"])[ok]))
    parts = {"localization_xy_nm": float(loc_xy), "localization_z_nm": float(np.sqrt(np.median((P["sz2"] + E["sz2"])[ok]))) if acq.dims == 3 else float("nan"),
             "exposure_blur_nm": float(np.sqrt(np.median(blur["P"] + blur["E"]))), "colour_delay_nm": float(np.sqrt(np.median(offs))),
             "piezo_z_nm": float(np.sqrt(np.median(piezo))) if acq.kind == "zstack" else 0.0,
             "mean_time_offset_s": float(np.mean(0.5 * (P["t"] + E["t"]))), "mean_EP_time_difference_s": float(np.mean(dt_ep))}
    return Measured(d, float(np.sqrt(np.mean(var_xy[ok]))), float(np.sqrt(np.mean(var_z[ok]))) if acq.dims == 3 else 0.0, acq.dims, acq.shift_nm, float(np.mean(ok)), parts)


# ----------------------------------------------------------------------------- kernel for anisotropic noise
_KCACHE: dict = {}


def kernel_aniso(sigma_xy: float, sigma_z: float, dims: int, shift_nm: float, n_draw: int = 30000, seed: int = 7) -> np.ndarray:
    """K[b, g] = P(measured distance in bin b | true distance D_GRID[g]) for an isotropic orientation, noise sigma_xy in x and y, sigma_z in z, a shift along x, and the
    distance taken in the first `dims` axes. Monte Carlo with common random numbers, smoothed along g."""
    key = (round(sigma_xy, 1), round(sigma_z, 1), dims, round(shift_nm, 1))
    if key in _KCACHE:
        return _KCACHE[key]
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(n_draw, 3)); u /= np.linalg.norm(u, axis=1, keepdims=True)
    n = rng.normal(size=(n_draw, 3)) * np.array([sigma_xy, sigma_xy, sigma_z])
    edges, grid = ep.EDGES, ep.D_GRID
    K = np.zeros((len(edges) - 1, len(grid)))
    for j, g in enumerate(grid):
        v = g * u + n
        v[:, 0] += shift_nm
        d = np.sqrt((v[:, :dims] ** 2).sum(1))
        K[:, j] = np.bincount(np.clip(np.digitize(d, edges) - 1, 0, len(edges) - 2), minlength=len(edges) - 1)
    K = gaussian_filter1d(K, 1.0, axis=1)
    K /= np.maximum(K.sum(0, keepdims=True), 1e-300)
    _KCACHE[key] = K
    return K
