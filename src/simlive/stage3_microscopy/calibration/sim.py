"""Simulated calibration acquisitions of the SAME instrument that recorded a movie: dark frames, flat-field (uniform fluorescent slide) series
and multi-colour bead slides.

The instrument (pixel gains, offsets, hot pixels, read-noise map, flat-field, chromatic shifts) comes from the run's `entropy`; every
acquisition has its own `session`, so its noise, drift, focus and flicker are independent of the movie's and of each other. Sample-related
sources (haze, autofluorescence, cross-talk, background, bleaching, locus motion, slice skew) are off: a slide or a bead coverslip has none.

    dark  : shutter closed, several exposure times  -> offset + fixed pattern, read noise, row noise, dark current, hot pixels, cosmic rays
    flat  : uniform slide at several brightness levels -> photon transfer curve (gain, read noise), illumination (flat-field), pixel gains,
            flicker
    beads : sub-resolution beads that emit in every colour -> PSF width, chromatic shift between colours, stage drift / focus drift in time
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from simlive.stage3_microscopy import errors as ERR
from simlive.stage3_microscopy.render import Imager

SAMPLE_SOURCES = ["haze", "autofluorescence", "crosstalk", "background", "bleaching", "exposure_motion", "slice_skew"]
KIND = {"dark": 1, "flat": 2, "beads": 3}


def instrument_errors(E: ERR.ErrorSet) -> ERR.ErrorSet:
    """The run's error set with the sample-related sources switched off (what a dark frame, a slide or a bead slide would show)."""
    return E.with_enabled([s for s in E.enabled_ids() if s not in SAMPLE_SOURCES])


def session_for(entropy: int, kind: str, index: int = 0) -> int:
    return int(np.random.SeedSequence([int(entropy), 7_000_000 + KIND[kind], int(index)]).generate_state(1, np.uint64)[0])


def acquire_dark(cfg: dict, E: ERR.ErrorSet, entropy: int, channel: str, exposure_s: float, n_frames: int, tag: int = 0) -> np.ndarray:
    """(n_frames, ny, nx) camera counts with the shutter closed, at one exposure time."""
    im = Imager(cfg, entropy=entropy, errors=instrument_errors(E), session=session_for(entropy, "dark", tag), n_frames=n_frames, channels=[channel])
    cam = ERR.Camera(im.shape, im.E, cfg["optics"], exposure_s, entropy, im.channel_index[channel], im.session)
    zero = np.zeros(im.shape)
    return np.stack([cam(zero, i) for i in range(n_frames)])


def acquire_flat(cfg: dict, E: ERR.ErrorSet, entropy: int, channel: str, level_photons: float, n_frames: int, tag: int = 0) -> np.ndarray:
    """(n_frames, ny, nx) camera counts of a uniform slide: `level_photons` photons per pixel (before QE) at the brightest part of the
    illumination pattern, with the lamp flicker of this acquisition."""
    Ec = instrument_errors(E)
    im = Imager(cfg, entropy=entropy, errors=Ec, session=session_for(entropy, "flat", tag), n_frames=n_frames, channels=[channel])
    flat = im._flat(channel, Ec)
    cam = im._camera(channel, Ec)
    out = []
    for t in range(n_frames):
        sig = level_photons * (flat if flat is not None else 1.0) * im.flicker(channel, t, Ec)
        out.append(cam(sig, t))
    return np.stack(out)


def _empty_cells(cfg: dict, three_d: bool) -> pd.DataFrame:
    harm = cfg["geometry"]["shape"]["harmonics"]
    cols = (["t", "cell_id", "y_um", "x_um", "radius_um", "aspect", "angle_rad", "bound_radius_um", "parent_id"]
            + [f"amp{int(k)}" for k in harm] + [f"phase{int(k)}" for k in harm] + (["z_um", "z_semi_um"] if three_d else []))
    return pd.DataFrame({c: [] for c in cols})


def make_bead_imager(cfg: dict, E: ERR.ErrorSet, entropy: int, n_beads: int, photons: float, n_frames: int, bead_z_um: float = 1.0,
                     margin_px: int = 16, tag: int = 0) -> tuple[Imager, pd.DataFrame]:
    """An Imager looking at a coverslip with `n_beads` bright sub-resolution beads that emit in every colour (nucleus channel = the green
    colour, loci channels = the others). Returns (imager, truth table of the bead positions: bead, y_um, x_um, z_um)."""
    import copy

    cfgb = copy.deepcopy(cfg)
    three_d = bool((cfg["acquisition"].get("z_stack") or {}).get("enabled"))
    if three_d:
        cfgb["acquisition"]["z_stack"]["center_um"] = bead_z_um
    n_ch = int((cfg.get("loci") or {}).get("n_loci", 0))
    cfgb.setdefault("loci", {"n_loci": 0, "wavelengths_nm": []})
    cfgb["loci"] = {**cfgb["loci"], "photons_per_locus_s": photons / cfg["acquisition"]["exposure_s"], "bleach_tau_s": None}
    cfgb["acquisition"]["n_frames"] = n_frames
    px_um = cfg["optics"]["pixel_size_nm"] / 1000.0
    ny, nx = (int(round(v / px_um)) for v in cfg["geometry"]["fov_um"])
    r = np.random.default_rng(np.random.SeedSequence([int(entropy), 7_100_000, tag]))
    pos = np.column_stack([r.uniform(margin_px, ny - margin_px, n_beads), r.uniform(margin_px, nx - margin_px, n_beads)]) * px_um
    z = np.full(n_beads, bead_z_um) + r.normal(0, 0.03, n_beads)             # beads sit on the coverslip
    beads = pd.DataFrame({"bead": np.arange(1, n_beads + 1), "y_um": pos[:, 0], "x_um": pos[:, 1], "z_um": z})
    rows = [{"t": t, "cell_id": b.bead, "locus_id": k, "y_um": b.y_um, "x_um": b.x_um, **({"z_um": b.z_um} if three_d else {})}
            for t in range(n_frames) for b in beads.itertuples() for k in range(n_ch)]
    truth = pd.DataFrame(rows) if rows else None
    spots = pd.DataFrame([{"t": t, "cell_id": b.bead, "y_um": b.y_um, "x_um": b.x_um, "photons": photons, **({"z_um": b.z_um} if three_d else {})}
                          for t in range(n_frames) for b in beads.itertuples()])
    im = Imager(cfgb, _empty_cells(cfg, three_d), truth, spots, None, entropy, instrument_errors(E), session=session_for(entropy, "beads", tag),
                n_frames=n_frames)
    return im, beads
