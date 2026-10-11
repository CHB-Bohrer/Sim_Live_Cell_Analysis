"""Run the calibration of a run's instrument end to end: simulate dark frames, flat-field series and a bead slide (sim.py), estimate the
instrument blindly from those images (fit.py), then compare every estimate with the simulation's true value.

    cal, table = run_calibration(cfg, entropy)          # table: quantity, channel, estimate, truth, error
    cal.save(run / "stage3_microscopy" / "calibration")
    scripts\\run.cmd python scripts\\run_calibration.py --run <run folder>      (from the command line, for an existing run)

Config block `calibration:` (all optional; these are the defaults, PLACEHOLDERS for what the real acquisitions will be):
    simulate: true, apply: true
    dark:  {n_frames: 60, exposures_s: [0.01, 1.0, 5.0], extra_frames: 8}
    flat:  {levels_photons: [0, 100, 300, 1000, 3000, 8000, 20000], n_frames: 12}
    flicker: {level_photons: 3000, n_frames: 100}
    beads: {n_beads: 60, photons: 30000, n_frames: 10, n_frames_3d: 6, z_um: 1.0}
"""
from __future__ import annotations

import copy

import numpy as np
import pandas as pd

from simlive.stage3_microscopy import errors as ERR
from simlive.stage3_microscopy.calibration import fit as F
from simlive.stage3_microscopy.calibration import sim as S
from simlive.stage3_microscopy.calibration.model import Calibration

DEFAULTS = {
    "simulate": True, "apply": True,
    "dark": {"n_frames": 60, "exposures_s": [0.01, 1.0, 5.0], "extra_frames": 8},
    "flat": {"levels_photons": [0, 100, 300, 1000, 3000, 8000, 20000], "n_frames": 12},
    "flicker": {"level_photons": 3000, "n_frames": 100},
    "beads": {"n_beads": 60, "photons": 30000, "n_frames": 10, "n_frames_3d": 6, "z_um": 1.0},
}


def settings_for(cfg: dict) -> dict:
    out = copy.deepcopy(DEFAULTS)
    for k, v in (cfg.get("calibration") or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k].update(v)
        else:
            out[k] = v
    return out


def channel_names(cfg: dict) -> list[str]:
    return ["nucleus"] + [f"locus{k}" for k in range(int((cfg.get("loci") or {}).get("n_loci", 0)))]


def _vignette_of(m: np.ndarray) -> float:
    ny, nx = m.shape
    h = max(min(ny, nx) // 20, 2)
    p = max(min(ny, nx) // 12, 2)
    center = m[ny // 2 - h:ny // 2 + h, nx // 2 - h:nx // 2 + h].mean()
    corners = np.mean([m[:p, :p].mean(), m[:p, -p:].mean(), m[-p:, :p].mean(), m[-p:, -p:].mean()])
    return float(1 - corners / center)


def run_calibration(cfg: dict, entropy: int, settings: dict | None = None, progress=print, channels: list[str] | None = None):
    """Returns (Calibration, comparison DataFrame). The comparison uses the simulation's true instrument; the Calibration does not."""
    from simlive.stage3_microscopy.render import Imager

    st = settings or settings_for(cfg)
    E = ERR.ErrorSet(cfg)
    Ec = S.instrument_errors(E)
    opt, acq = cfg["optics"], cfg["acquisition"]
    run_exp = float(acq["exposure_s"])
    px_um = opt["pixel_size_nm"] / 1000.0
    shape = tuple(int(round(v / px_um)) for v in cfg["geometry"]["fov_um"])
    chans = channels or channel_names(cfg)
    rows, maps, values = [], {}, {}
    gain = float(opt.get("gain_adu_per_e", 1.0))

    def add(q, c, unit, est, true, note=""):
        err = (est - true) if (true is not None and est is not None and np.isfinite(est) and np.isfinite(true)) else None
        rows.append({"quantity": q, "channel": c, "unit": unit, "estimate": est, "truth": true, "error": err,
                     "error_pct": (100 * err / true) if (err is not None and true) else None, "note": note})

    for ci, ch in enumerate(chans):
        idx = {"nucleus": 0}.get(ch, int(ch[5:]) + 1 if ch.startswith("locus") else ci)
        # ---- the instrument's true camera (only used for the comparison)
        cam_true = ERR.Camera(shape, Ec, opt, run_exp, entropy, idx, S.session_for(entropy, "dark", 99))
        # ---- dark frames
        d = st["dark"]
        exps = sorted(set([float(e) for e in d["exposures_s"]] + [run_exp]))
        progress(f"  calibration [{ch}]: dark frames at {len(exps)} exposures")
        stacks = {}
        for k, e in enumerate(exps):
            n = int(d["n_frames"]) if (e == exps[0] or e == run_exp) else int(d["extra_frames"])
            stacks[e] = S.acquire_dark(cfg, E, entropy, ch, e, n, tag=k)
        dk = F.fit_dark(stacks, run_exp)
        dark_example = stacks[run_exp][0].copy()
        del stacks
        # ---- flat slide: photon transfer curve, illumination, pixel gain, flicker
        fl = st["flat"]
        progress(f"  calibration [{ch}]: flat-field slide at {len(fl['levels_photons'])} brightness levels")
        levels, top_key = {}, None
        for k, lv in enumerate(fl["levels_photons"]):
            stk = S.acquire_flat(cfg, E, entropy, ch, float(lv), int(fl["n_frames"]), tag=k)
            if np.mean(stk == stk.max()) > 1e-3 and stk.max() > 1000:           # saturated level: not usable for the transfer curve
                continue
            levels[float(lv)] = stk
        ptc = F.fit_ptc(levels, dk["offset_run_map"])
        top_key = max(levels)
        pt = max(ptc["ptc_points"], key=lambda p: p[0])
        flat = F.fit_flat(levels[top_key], dk["offset_run_map"], pt[1])
        flat_example = levels[top_key][0].copy()
        del levels
        fk = st["flicker"]
        progress(f"  calibration [{ch}]: {fk['n_frames']}-frame flicker series")
        flick = F.fit_flicker(S.acquire_flat(cfg, E, entropy, ch, float(fk["level_photons"]), int(fk["n_frames"]), tag=100), dk["offset_run_map"])
        maps[ch] = {"offset": dk["offset_run_map"], "flat": flat["flat_norm"], "hot": dk["hot_mask"]}
        g_est = ptc["gain_adu_per_e"]
        values[ch] = {"gain_adu_per_e": g_est, "offset_adu": dk["offset_adu"], "read_noise_adu": dk["read_noise_adu"],
                      "read_noise_e": dk["read_noise_adu"] / g_est if g_est > 0 else float("nan"),
                      "row_noise_adu": dk["row_noise_adu"], "dsnu_adu": dk["dsnu_adu"], "dark_current_e_per_s": dk["dark_rate_adu_s"] / g_est if g_est > 0 else float("nan"),
                      "n_hot_pixels": dk["n_hot"], "hot_fraction": dk["hot_fraction"], "hot_excess_e": dk["hot_excess_adu"] / g_est if g_est > 0 else float("nan"),
                      "cosmic_rate_per_frame": dk["cosmic_rate_per_frame"], "cosmic_events": dk["cosmic_events"],
                      "cosmic_frames": dk["cosmic_frames"], "cosmic_integral_e": dk["cosmic_integral_adu"] / g_est if g_est > 0 else float("nan"),
                      "prnu_cv": flat["prnu_cv"], "vignette": flat["vignette"], "tilt": flat["tilt"],
                      "flicker_cv": flick["flicker_cv"], "flicker_memory": flick["flicker_memory"], "ptc_points": ptc["ptc_points"], "ptc_intercept_adu2": ptc["ptc_intercept_adu2"],
                      "digitised": dk["quantized"], "qe": None}
        maps[ch]["dark_example"] = dark_example.astype(np.float32)
        maps[ch]["flat_example"] = flat_example.astype(np.float32)
        # ---- comparison with the truth
        v = values[ch]
        add("gain", ch, "ADU/e-", v["gain_adu_per_e"], gain)
        add("offset", ch, "ADU", v["offset_adu"], float(opt.get("offset_adu", 100.0)))
        rn_true = (float(opt["read_noise_e"]) * (float(np.sqrt(np.mean(cam_true.read_map ** 2))) if cam_true.read_map is not None else 1.0)
                   * gain) if Ec.on("read_noise") else 0.0
        add("read noise", ch, "ADU", v["read_noise_adu"], rn_true)
        add("row noise", ch, "ADU", v["row_noise_adu"], float(Ec.par("row_noise", "adu")) if Ec.on("row_noise") else 0.0)
        add("offset pattern (DSNU)", ch, "ADU", v["dsnu_adu"], float(np.std(cam_true.dsnu)) if cam_true.dsnu is not None else 0.0)
        add("dark current", ch, "e-/px/s", v["dark_current_e_per_s"], float(Ec.par("dark_current", "e_per_s")) if Ec.on("dark_current") else 0.0)
        hot_true = (cam_true.hot > 0) if cam_true.hot is not None else np.zeros(shape, bool)
        add("hot pixels found", ch, "count", float(v["n_hot_pixels"]), float(hot_true.sum()),
            f"detected {int((dk['hot_mask'] & hot_true).sum())} of {int(hot_true.sum())} true hot pixels, {int((dk['hot_mask'] & ~hot_true).sum())} false")
        if Ec.on("cosmic_rays"):
            add("cosmic-ray rate", ch, "events/frame", v["cosmic_rate_per_frame"], float(Ec.par("cosmic_rays", "rate_per_frame")),
                f"{v['cosmic_events']} events in {v['cosmic_frames']} dark frames: poorly constrained unless many frames are taken")
        add("pixel gain spread (PRNU)", ch, "relative", v["prnu_cv"], float(np.std(cam_true.prnu)) if cam_true.prnu is not None else 0.0)
        add("vignetting (loss at the corners)", ch, "relative", v["vignette"], float(Ec.par("flat_field", "vignette")) if Ec.on("flat_field") else 0.0)
        add("illumination tilt", ch, "relative", v["tilt"], float(Ec.par("flat_field", "tilt")) if Ec.on("flat_field") else 0.0,
            "the tilt direction is random; only its size is compared")
        add("flicker", ch, "relative (cv)", v["flicker_cv"], float(Ec.par("flicker", "cv")) if Ec.on("flicker") else 0.0)
        add("flicker memory", ch, "", v["flicker_memory"], float(Ec.par("flicker", "memory")) if Ec.on("flicker") else None)

    # ---- bead slide: PSF, chromatic shift, drift
    b = st["beads"]
    three_d = bool((acq.get("z_stack") or {}).get("enabled"))
    nb_frames = int(b["n_frames_3d"] if three_d else b["n_frames"])
    progress(f"  calibration: bead slide, {b['n_beads']} beads, {nb_frames} frames" + (" (z-stacks)" if three_d else ""))
    im, beads = S.make_bead_imager(cfg, E, entropy, int(b["n_beads"]), float(b["photons"]), nb_frames, float(b["z_um"]))
    ims = {c: [] for c in chans}
    for t in range(nb_frames):
        fr = im.frame(t, channels=chans)["images"]
        for c in chans:
            ims[c].append(np.rint(fr[c]).clip(0, 65535).astype(np.uint16) if Ec.on("adc") else fr[c].astype(np.float32))
    arr = {c: np.stack(v) for c, v in ims.items()}
    # correct the beads with the dark / flat just estimated, like a real analysis would
    cal0 = Calibration(chans, shape, run_exp, px_um, "nucleus", maps, values, {}, [])
    for c in chans:
        arr[c] = np.stack([cal0.correct(a, c) for a in arr[c]])
    bf = F.fit_beads(arr, px_um, list(im.zk) if three_d else None, ref="nucleus", half=6)
    example_beads = {c: (arr[c][0].max(axis=0) if three_d else arr[c][0]).astype(np.float32) for c in chans}
    del arr, ims
    for c in chans:
        maps[c]["bead_example"] = example_beads[c]
    cal = Calibration(chans, shape, run_exp, px_um, "nucleus", maps, values, bf, list(map(float, im.zk)) if three_d else [])
    # truth for the bead results
    for c in chans:
        if c == "nucleus":
            continue
        ty, tx = im.channel_shift_px(c, Ec)
        add("chromatic shift y", c, "px", bf["shift_px"][c][0], ty, f"± {bf['shift_px_sd'][c]:.4f} px (1 sd of the median)")
        add("chromatic shift x", c, "px", bf["shift_px"][c][1], tx)
        if three_d:
            add("chromatic shift z", c, "um", bf["shift_z_um"][c], im.channel_shift_z_um(c, Ec))
    for c in chans:
        sig_true = float(np.median([im.sigma_px(c, t, Ec) for t in range(nb_frames)]))
        add("PSF width (lateral sigma)", c, "px", bf["psf_sigma_px"][c], sig_true)
        if three_d:
            add("PSF width (axial sigma)", c, "um", bf["psf_sigma_z_um"][c], float(im.sigma_z_um(c, Ec)))
    est = np.array(bf["drift_px_series"])
    true = np.array([im.drift_px(t, Ec) for t in range(nb_frames)])
    true = true - true[0]
    add("stage drift (rms error of the estimated path)", "all", "px", float(np.sqrt(np.mean((est - true) ** 2))), 0.0,
        f"true path rms {float(np.sqrt(np.mean(true ** 2))):.3f} px")
    if three_d:
        ez = np.array(bf["drift_z_um_series"])
        tz = np.array([im.drift_z_um(t, Ec) + im.defocus_um(t, Ec) for t in range(nb_frames)])
        tz = tz - tz[0]
        add("axial drift + focus (rms error of the estimated path)", "all", "um", float(np.sqrt(np.mean((ez - tz) ** 2))), 0.0,
            f"true path rms {float(np.sqrt(np.mean(tz ** 2))):.3f} um")
    return cal, pd.DataFrame(rows)
