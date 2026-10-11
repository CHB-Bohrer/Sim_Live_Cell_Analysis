"""The dashboard's 🎯 Calibration tab: calibrate the imaging instrument from simulated dark frames, a flat-field slide and a bead slide.

The simulation acquires those calibration images with the SAME instrument that recorded the movie (same camera patterns, same
illumination, same chromatic shifts), estimates the instrument from the images alone, and then shows how close each estimate is to the
true value the simulation used. The estimates can be applied to the movie: dark-subtraction, flat-field, hot-pixel replacement and
removal of the chromatic shift, which improves the locus positions (the scores before and after are shown).
Imports numpy / scipy / pandas only, never torch.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

import style
from simlive.stage3_microscopy.calibration.model import Calibration
from simlive.stage3_microscopy.calibration.run import run_calibration, settings_for
from simlive.stage3_microscopy.render import entropy_for

CH_SHORT = {"nucleus": "nucleus", "locus0": "locus 0 (promoter)", "locus1": "locus 1 (enhancer)"}


@st.cache_resource(show_spinner="Loading the calibration…")
def _load(path: str, stamp: float) -> Calibration:
    return Calibration.load(path)


def _tolerance(q: str) -> float:
    """Percent error counted as 'close enough' for each quantity (rough, stated in the table's note)."""
    if "cosmic" in q:
        return 40.0
    if q.startswith("flicker"):
        return 25.0
    if q in ("row noise", "offset pattern (DSNU)", "dark current"):
        return 15.0
    return 8.0


def _status(row) -> str:
    if row.error_pct is None or not np.isfinite(row.error_pct):
        if row.truth is not None and np.isfinite(row.error if row.error is not None else np.nan):
            tol = {"chromatic shift y": 0.03, "chromatic shift x": 0.03, "chromatic shift z": 0.06, "offset": 0.5}.get(row.quantity, None)
            if tol is not None:
                return "✅" if abs(row.error) <= tol else "⚠"
            if row.quantity.startswith(("stage drift", "axial drift")):
                return "✅" if row.estimate < (0.03 if "px" in row.unit else 0.06) else "⚠"
            if "vignetting" in row.quantity or "tilt" in row.quantity:
                return "✅" if abs(row.error) < 0.01 else "⚠"
        return "·"
    if row.quantity in ("hot pixels found",):
        return "✅" if abs(row.error) <= max(1, 0.05 * row.truth) else "⚠"
    return "✅" if abs(row.error_pct) <= _tolerance(row.quantity) else "⚠"


def _tidy(table: pd.DataFrame) -> pd.DataFrame:
    t = table.copy()
    t["ok"] = [_status(r) for r in t.itertuples()]
    t["estimate"] = t.estimate.map(lambda v: None if v is None else float(f"{v:.4g}"))
    t["truth"] = t.truth.map(lambda v: None if v is None or not np.isfinite(v) else float(f"{v:.4g}"))
    t["error %"] = t.error_pct.map(lambda v: None if v is None or not np.isfinite(v) else round(float(v), 1))
    return t[["ok", "quantity", "channel", "unit", "estimate", "truth", "error %", "note"]]


def _simulate(run: Path, R: dict):
    cfg = R["cfg"]
    cdir = run / "stage3_microscopy" / "calibration"
    with st.status("Simulating calibration acquisitions…", expanded=True) as status:
        cal, table = run_calibration(cfg, entropy_for(int(cfg["seed"])), progress=lambda m: st.write(m))
        cal.save(cdir)
        table.to_csv(cdir / "comparison.csv", index=False)
        status.update(label="Done", state="complete", expanded=False)


def _overview() -> None:
    c = st.columns(4)
    with c[0]:
        style.card(st, "Acquire", "Take pictures of things whose true appearance is known: <b>dark frames</b> (shutter closed), a <b>uniform "
                   "slide</b> at several brightnesses, and a slide of <b>tiny multi-colour beads</b>.", 1)
    with c[1]:
        style.card(st, "Estimate", "From those pictures alone, work out the instrument: gain, noise, hot pixels, illumination pattern, blur, "
                   "colour offsets, drift. The estimator never sees the true values.", 2)
    with c[2]:
        style.card(st, "Check", "Here the simulation <b>does</b> know the truth, so every estimate is compared with it. That shows how well the "
                   "calibration method works before it is used on real data.", 3)
    with c[3]:
        style.card(st, "Apply", "Correct the movie with the estimates (subtract the dark offset, divide by the flat-field, replace hot pixels, "
                   "undo the chromatic shift) and see how much the locus positions improve.", 4)


def _dark_section(cal: Calibration) -> None:
    st.markdown("### Dark frames: the camera with the shutter closed")
    st.caption("No light reaches the sensor, so everything seen is the camera itself: a constant **offset**, a fixed per-pixel pattern (**DSNU**), "
               "a few permanently bright **hot pixels**, random **read noise**, stripes from **row noise**, slowly accumulating **dark current** "
               "(measured by changing the exposure time) and occasional **cosmic rays**.")
    ch = st.selectbox("Channel", cal.channels, format_func=CH_SHORT.get, key="cal_dark_ch")
    m = cal.maps[ch]
    c = st.columns(3)
    for col, (title, img, cmap) in zip(c, [("One dark frame", m["dark_example"], "gray"), ("Offset map (average of many frames)", m["offset"], "viridis")]):
        fig, ax = plt.subplots(figsize=(5, 4.4))
        lo, hi = np.percentile(img, [1, 99.7])
        ax.imshow(img, cmap=cmap, vmin=lo, vmax=hi)
        ax.set_title(title, fontsize=11); ax.axis("off")
        col.pyplot(fig, width="stretch"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(5, 4.4))
    ax.imshow(np.zeros_like(m["offset"]), cmap="gray", vmin=0, vmax=1)
    ys, xs = np.nonzero(m["hot"])
    ax.plot(xs, ys, "o", mfc="none", mec="red", ms=7)
    ax.set_title(f"Hot pixels found: {int(m['hot'].sum())}", fontsize=11); ax.axis("off")
    c[2].pyplot(fig, width="stretch"); plt.close(fig)
    v = cal.values[ch]
    style.stat_grid(st, [("Offset", f"{v['offset_adu']:.1f} ADU", "constant counts added by the camera electronics"),
                         ("Read noise", f"{v['read_noise_adu']:.2f} ADU · {v['read_noise_e']:.2f} e⁻", "random electronic noise per pixel per frame"),
                         ("Row noise", f"{v['row_noise_adu']:.2f} ADU", "random offset of each pixel row in each frame"),
                         ("Offset pattern (DSNU)", f"{v['dsnu_adu']:.2f} ADU", "fixed pixel-to-pixel spread of the offset"),
                         ("Dark current", f"{v['dark_current_e_per_s']:.2f} e⁻/px/s", "electrons per pixel per second of exposure in the dark"),
                         ("Hot pixels", f"{v['n_hot_pixels']} ({100 * v['hot_fraction']:.3f} %)", "pixels that read high at every exposure"),
                         ("Cosmic rays", f"{v['cosmic_rate_per_frame']:.3f} per frame", f"{v['cosmic_events']} events in {v['cosmic_frames']} dark frames (rare: needs many frames)")])


def _flat_section(cal: Calibration) -> None:
    st.markdown("### Flat-field slide: a uniform light source")
    st.caption("A slide that glows evenly at a known-to-be-uniform brightness. Taken at several brightnesses it gives the **photon transfer curve**: "
               "light is made of photons, so the noise variance grows in proportion to the signal; the slope is the **gain** (counts per detected "
               "electron) and the intercept the noise floor. Averaged at high brightness it shows the **illumination pattern** (vignetting, tilt) "
               "and the **pixel-to-pixel gain** differences (PRNU). Its frame-to-frame brightness changes give the **flicker**.")
    ch = st.selectbox("Channel", cal.channels, format_func=CH_SHORT.get, key="cal_flat_ch")
    m, v = cal.maps[ch], cal.values[ch]
    c = st.columns(3)
    P = np.array(v["ptc_points"])
    fig, ax = plt.subplots(figsize=(5, 4.4))
    ax.plot(P[:, 0], P[:, 1], "o", color="#2563eb", label="measured")
    xs = np.linspace(0, P[:, 0].max(), 50)
    ax.plot(xs, v["gain_adu_per_e"] * xs + v["ptc_intercept_adu2"], "-", color="#0f766e", label=f"fit: slope {v['gain_adu_per_e']:.3f} ADU/e⁻")
    ax.set_xlabel("mean signal (ADU)"); ax.set_ylabel("variance (ADU²)"); ax.set_title("Photon transfer curve", fontsize=11); ax.legend(fontsize=9)
    c[0].pyplot(fig, width="stretch"); plt.close(fig)
    for col, (title, img) in zip(c[1:], [("Normalised flat (illumination × pixel gain)", m["flat"]), ("One raw flat frame", m["flat_example"])]):
        fig, ax = plt.subplots(figsize=(5, 4.4))
        lo, hi = np.percentile(img, [0.5, 99.5])
        im_ = ax.imshow(img, cmap="magma", vmin=lo, vmax=hi); ax.set_title(title, fontsize=11); ax.axis("off")
        fig.colorbar(im_, ax=ax, fraction=0.04)
        col.pyplot(fig, width="stretch"); plt.close(fig)
    style.stat_grid(st, [("Gain", f"{v['gain_adu_per_e']:.3f} ADU/e⁻", "counts per detected electron (the slope of the photon transfer curve). The quantum efficiency cannot be measured this way."),
                         ("Vignetting", f"{100 * v['vignette']:.1f} % at the corners", "how much dimmer the corners of the field are"),
                         ("Illumination tilt", f"{100 * v['tilt']:.1f} %", "left-right / top-bottom brightness gradient"),
                         ("Pixel gain spread (PRNU)", f"{100 * v['prnu_cv']:.2f} %", "relative spread of the gain between pixels"),
                         ("Flicker", f"{100 * v['flicker_cv']:.2f} % (memory {v['flicker_memory']:.2f})", "random frame-to-frame brightness change of the excitation")])


def _beads_section(cal: Calibration) -> None:
    b = cal.beads
    st.markdown("### Bead slide: tiny sources that emit in every colour")
    st.caption("Beads much smaller than the blur are points of light, so what the microscope records of one is its **point spread function**. "
               "Every bead appears in every colour: if the colours are imaged at slightly different places, the same bead is found at different "
               "positions in different channels. That difference is the **chromatic shift** and it is the same for every bead. Beads that "
               "do not move also reveal **stage drift** (all beads move together over time) and, in z-stacks, the axial focus error.")
    ref = b["ref"]
    c = st.columns([2, 2, 3])
    ch = c[0].selectbox("Show channel", cal.channels, format_func=CH_SHORT.get, key="cal_bead_ch")
    img = cal.maps[ch]["bead_example"]
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    ax.imshow(img, cmap="gray", vmin=np.percentile(img, 1), vmax=np.percentile(img, 99.9))
    P = np.array(b["bead_positions_px"])
    ax.plot(P[:, 1], P[:, 0], "o", mfc="none", mec="cyan", ms=9)
    ax.set_title(f"Beads in {CH_SHORT[ch]} ({len(P)} found)", fontsize=11); ax.axis("off")
    c[0].pyplot(fig, width="stretch"); plt.close(fig)
    others = [x for x in cal.channels if x != ref]
    if others:
        fig, ax = plt.subplots(figsize=(5.2, 4.8))
        for k, x in enumerate(others):
            r = np.array(b["residual_px"][x])
            ax.plot(r[:, 1], r[:, 0], ".", alpha=0.6, label=f"{CH_SHORT[x]}")
            sy, sx = b["shift_px"][x]
            ax.plot(sx, sy, "k+", ms=14, mew=2)
        ax.axhline(0, color="#999", lw=0.8); ax.axvline(0, color="#999", lw=0.8)
        ax.set_xlabel("x offset from the nuclear channel (px)"); ax.set_ylabel("y offset (px)")
        ax.set_title("Each dot = one bead; + = the estimated shift", fontsize=11); ax.legend(fontsize=8); ax.set_aspect("equal")
        c[1].pyplot(fig, width="stretch"); plt.close(fig)
    d = np.array(b["drift_px_series"])
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    ax.plot(d[:, 1], d[:, 0], "o-", color="#b45309")
    ax.set_xlabel("x drift (px)"); ax.set_ylabel("y drift (px)"); ax.set_title("Stage drift measured from the beads", fontsize=11)
    c[2].pyplot(fig, width="stretch"); plt.close(fig)
    rows = []
    for x in cal.channels:
        rows.append({"channel": CH_SHORT[x], "shift y (px)": round(b["shift_px"][x][0], 3), "shift x (px)": round(b["shift_px"][x][1], 3),
                     "shift z (µm)": round(b["shift_z_um"].get(x, 0.0), 3) if cal.z_um else None,
                     "PSF σ (px)": round(b["psf_sigma_px"][x], 3), "PSF σ (nm)": round(b["psf_sigma_px"][x] * cal.pixel_size_um * 1000, 0),
                     "axial σ (µm)": round(b["psf_sigma_z_um"][x], 3) if cal.z_um else None, "fits used": b["n_fits"][x]})
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    st.caption("Shifts are relative to the nuclear channel. The 'fits used' column is the number of bead measurements (beads × frames).")


def _effect_section(run: Path, cal: Calibration) -> None:
    st.markdown("### Does it help? The movie analysed without and with the calibration")
    sj = run / "stage6_analysis" / "scores.json"
    if not sj.exists():
        st.info("This run has no locus analysis to compare.")
        return
    sc = json.loads(sj.read_text())
    rows = []
    for key, label in (("truth", "2D, uncorrected"), ("truth_calibrated", "2D, calibration applied"),
                       ("truth_3d", "3D, uncorrected"), ("truth_3d_calibrated", "3D, calibration applied")):
        if key in sc and "rms_error_um" in sc[key]:
            s = sc[key]
            rows.append({"analysis": label, "detection rate": round(s["detection_rate"], 3), "median error (nm)": round(1000 * s["median_error_um"], 1),
                         "95th percentile (nm)": round(1000 * s["p95_error_um"], 0), "rms error (nm)": round(1000 * s["rms_error_um"], 0)})
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        st.caption("Errors are the distance between each located locus and its true position (true cell IDs only). The calibrated analysis "
                   "subtracts the measured chromatic shift and works on dark-subtracted, flat-fielded images with hot pixels replaced. "
                   "Cosmic rays cannot be removed this way (they are random per frame). The median error is the typical locus; the rms error "
                   "is dominated by the occasional wrong detection.")
    else:
        st.info("This run was made without 'apply calibration'; make a new run with the calibration switched on to compare.")
    ch = st.selectbox("Show the correction on channel", cal.channels, format_func=CH_SHORT.get, key="cal_eff_ch")
    import tifffile
    s3 = run / "stage3_microscopy"
    try:
        with tifffile.TiffFile(s3 / f"{ch}.tif") as tf:
            n_pages = len(tf.pages)
        t = st.slider("Frame", 0, n_pages - 1, 0, key="cal_eff_t")
        raw = tifffile.imread(s3 / f"{ch}.tif", key=int(t)).astype(float)
    except Exception:
        return
    corr = cal.correct(raw, ch)
    y0, x0, h = raw.shape[0] // 2 - 80, raw.shape[1] // 2 - 80, 160
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, (title, im_) in zip(axes, [(f"Recorded (frame {t})", raw), ("Dark-subtracted, flat-fielded, hot pixels replaced", corr), ("What the correction changed", corr - raw)]):
        sub = im_[y0:y0 + h, x0:x0 + h]
        if title.startswith("What"):
            m = max(np.percentile(np.abs(sub), 99.5), 1e-9)
            ax.imshow(sub, cmap="RdBu_r", vmin=-m, vmax=m)
        else:
            ax.imshow(sub, cmap="gray", vmin=np.percentile(sub, 1), vmax=np.percentile(sub, 99.8))
        ax.set_title(title, fontsize=11); ax.axis("off")
    st.pyplot(fig, width="stretch"); plt.close(fig)


def _method_section() -> None:
    with st.expander("How each estimate is made (and what cannot be measured)"):
        st.markdown(
            "- **Offset, DSNU, hot pixels:** average many dark frames. Hot pixels are the ones far above the typical offset at the shortest exposure.\n"
            "- **Dark current:** how the average dark level grows with the exposure time (a straight line per pixel).\n"
            "- **Read noise, row noise:** the frame-to-frame scatter of dark pixels; row noise is the part shared by a whole row.\n"
            "- **Cosmic rays:** blobs in single dark frames far above the noise. Rare, so the rate is only roughly known unless very many frames are taken.\n"
            "- **Gain:** photon transfer curve of the flat slide, after dividing out the common flicker and rejecting cosmic-ray hits.\n"
            "- **Illumination and PRNU:** a smooth model (constant + tilt + quadratic) fitted to the normalised average flat; what is left is the pixel gain pattern.\n"
            "- **Flicker:** frame-to-frame change of the mean brightness over a long series.\n"
            "- **PSF width, chromatic shift, drift, axial focus:** fit a small Gaussian to every bead in every colour and frame; "
            "colours are compared bead by bead.\n"
            "- **Cannot be measured here:** the quantum efficiency (needs a light source of known photon flux), anything that differs between the "
            "calibration slide and the real sample (autofluorescence, haze, cross-talk between colours when both are bright), and cosmic rays in the movie itself.")


def render(st_, run_path: Path, R: dict):
    style.intro(st, "Calibrating the microscope: dark frames, a flat-field slide and a bead slide are simulated with the same instrument that recorded "
                    "the movie, the instrument is estimated from them blindly, and every estimate is checked against the true value.",
                "Press the button to simulate the calibration (about 15 to 60 seconds), read the table of estimate versus truth, look at each calibration "
                "image below, and see how much the locus positions improve when the calibration is applied.")
    run = Path(run_path)
    _overview()
    cdir = run / "stage3_microscopy" / "calibration"
    if not (cdir / "calibration.json").exists():
        st.info("No calibration for this run yet.")
        if not (R["cfg"].get("loci")):
            st.warning("This run has no loci channels; the calibration needs at least the nuclear and one locus channel.")
            return
        if st.button("▶ Simulate the calibration for this run's instrument", type="primary", key="cal_go"):
            _simulate(run, R)
            st.rerun()
        return
    cal = _load(str(cdir), (cdir / "calibration.json").stat().st_mtime)
    table = pd.read_csv(cdir / "comparison.csv").replace({np.nan: None})
    t = _tidy(table)
    n_ok = int((t.ok == "✅").sum())
    st.markdown(f"### Estimated versus true: {n_ok} of {int((t.ok != '·').sum())} within tolerance")
    st.caption("✅ close enough, ⚠ off by more than expected (tolerances: 8% for most, 15% for the small noise terms, 25% for flicker, 40% for the "
               "cosmic-ray rate; 0.03 px / 0.06 µm for registration). 'truth' is what the simulation actually used; the estimate saw only the images.")
    st.dataframe(t, hide_index=True, width="stretch", column_config={"note": st.column_config.TextColumn(width="large")})
    c = st.columns([1, 3])
    if c[0].button("↻ Re-simulate", key="cal_redo", help="Run the calibration again with this run's settings."):
        _simulate(run, R)
        st.rerun()
    _dark_section(cal)
    _flat_section(cal)
    _beads_section(cal)
    _effect_section(run, cal)
    _method_section()
