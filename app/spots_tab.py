"""The dashboard's 🔎 Spot finding tab: how stage 6 finds each locus, step by step, and how well it does against the truth.

Four parts: (1) the method and where it comes from (Bennabi et al. 2025, rescaled from Mesp2 poles to diffraction-limited spots);
(2) every step on one cell and frame, with settings you can change and the error of every frame of that cell; (3) the run's scores and
its worst detections; (4) a comparison of methods on this run (scripts/compare_localizers.py).
Imports numpy / scipy / scikit-image / pandas only, never torch.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from matplotlib.patches import Circle, Rectangle

import style
from simlive.io.runs import REPO_ROOT
from simlive.stage5_linking.celldata import load_cell
from simlive.stage6_analysis import spots as SP
from simlive.stage6_analysis.compare import VARIANTS, VARIANTS_3D, truth_with_photons
from simlive.stage6_analysis.localize import analyse_channel, channel_inputs

PAPER_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC12372888"
CODE_URL = "https://github.com/gregorlab/SizeProject2024"
CH_SHORT = {"locus0": "locus 0 (promoter)", "locus1": "locus 1 (enhancer)"}
FAR_UM = 0.5

STEPS = [
    ("1 Median filter", "rank median, 7 × 7 px", "median, median_px × median_px (3: a spot is ~2 px wide; 5 or more erases it)", "median_px"),
    ("2 Region", "inside the gastruloid (bright-field mask)", "inside the cell's own nuclear mask, grown by dilate_px", "dilate_px"),
    ("· Background", "none", "a wider median (background_px) of the filtered image, subtracted: removes the nuclear cross-talk ghost, "
                              "autofluorescence and haze under the spot ('cell' = one value per cell)", "background"),
    ("3 Threshold", "one per movie: median over frames of each frame's Otsu threshold (or fixed, 70)",
     "one per movie (threshold_scope): median over frames of each frame's min_snr × robust noise ('noise') or Otsu ('otsu')", "threshold"),
    ("4 Area", "pixels above it; holes filled; objects < 5000 px removed; eroded", "pixels above it (threshold_on: filtered image as the "
                                                                                   "paper, or the smoothed image); objects < min_area_px removed", "threshold_on"),
    ("5 Smoothing", "Gaussian, σ = 20 px (pole scale)", "Gaussian, σ = smooth_sigma_px (≈ the PSF: a matched filter)", "smooth_sigma_px"),
    ("6 Local maxima", "peak_local_max, min_distance 100 px: each maximum is a pole",
     "peak_local_max, min_distance = min_sep_px: candidates; one locus per colour", "min_sep_px"),
    ("· Choice", "(all poles are kept and counted)", "'brightest' candidate, or 'track' (ours, not the paper's): the candidate nearest the "
                                                    "locus' position in the neighbouring frames, so one-frame events such as cosmic rays are skipped", "select"),
    ("7 Position", "integer pixel of the maximum", "Gaussian fit to the unfiltered pixels (refine: gaussian), centroid, or none", "refine"),
]


# ----------------------------------------------------------------------------- cached computations
@st.cache_data(show_spinner=False, max_entries=4)
def _truth(run: str, stamp: float) -> pd.DataFrame:
    return truth_with_photons(Path(run))


@st.cache_resource(show_spinner="Loading the cell…", max_entries=4)
def _inputs(sdir: str, cid: int, name: str, calib: str | None, dilate_px: int) -> dict:
    """The crops of one colour of one cell as the localizer sees them (calibration-corrected if asked)."""
    return channel_inputs(load_cell(sdir, cid), name, dilate_px, calib)


@st.cache_data(show_spinner="Running the spot finder on every frame of this cell…", max_entries=12)
def _analyse(sdir: str, cid: int, name: str, params: str, calib: str | None, dilate_px: int):
    """All frames of one colour of one cell with the given settings: one row per frame + the candidate positions."""
    ci = _inputs(sdir, cid, name, calib, dilate_px)
    a = analyse_channel(ci["stack"], ci["region"], ci["centre"], **json.loads(params))
    fr, (sy, sx) = ci["frames"], ci["shift"]
    rows = []
    for i, r in enumerate(fr.itertuples()):
        f, c = a["fits"][i], a["candidates"][i]
        rows.append({"i": i, "t": int(r.t), "crop_y0": int(r.crop_y0), "crop_x0": int(r.crop_x0), "thr_frame": a["thr_frame"][i],
                     "thr": a["thr"][i], "n_candidates": 0 if c is None else len(c["peaks"]), "chosen": a["choice"][i],
                     "found": f is not None, "y_crop": f["pos"][0] - sy if f else np.nan, "x_crop": f["pos"][1] - sx if f else np.nan,
                     "snr": f["snr"] if f else np.nan, "fit_ok": f["fit_ok"] if f else None, "sigma_fit_px": f["sigma_fit_px"] if f else np.nan})
    cands = [None if c is None else c["peaks"] for c in a["candidates"]]
    return pd.DataFrame(rows), cands


def _frame_steps(sdir: str, cid: int, name: str, params: dict, calib: str | None, dilate_px: int, i: int, thr: float):
    ci = _inputs(sdir, cid, name, calib, dilate_px)
    return SP.find_spot(ci["stack"][i], ci["region"][i], thr, steps=True, **params), ci["region"][i]


# ----------------------------------------------------------------------------- parts
def _part_method(cfg: dict) -> dict:
    st.markdown("### 1 · The method and where it comes from")
    la = cfg.get("locus_analysis") or {}
    method = la.get("method", "brightest" if "method" not in la and "median_px" not in la else "median_peaks")
    p = {**SP.DEFAULTS, "dilate_px": 2, **la}
    st.markdown(f"This run located the loci with **`{method}`**. The `median_peaks` method follows the detection steps of "
                f"[Bennabi et al. 2025, *Science Advances* 11, eadv7790]({PAPER_URL}) ([their code]({CODE_URL}), "
                "`modules/raw_live_analysis.py`), who find Mesp2 expression *poles* in live gastruloid movies. Their objects are hundreds of "
                "pixels across; a locus here is about 2 pixels across, so every size had to be rescaled, and two steps were added.")
    lines = ["| step | Bennabi et al. 2025 (Mesp2 poles) | here (one locus per colour) | this run |", "|---|---|---|---|"]
    lines += [f"| **{s}** | {a} | {b} | " + (f"`{k}: {p.get(k)}`" if method == "median_peaks" else "—") + " |" for s, a, b, k in STEPS]
    st.markdown("\n".join(lines))
    with st.expander("Why the sizes change, and what was measured"):
        st.markdown(
            "- **Median filter.** Theirs (7 px) is tiny next to a pole, so it only removes pixel noise. A 3 × 3 median already halves a spot "
            "of σ ≈ 0.8 px; a 7 × 7 median removes it completely. That is why a *wider* median is used here as the **background** under a spot.\n"
            "- **What the median removes:** hot pixels and single-pixel noise. It cannot remove a **cosmic ray**, which is as wide as the PSF; "
            "in a max projection of a z-stack every plane's cosmic rays end up in the picture. The **track** choice skips them because they "
            "last one frame. In a finely sampled z-stack an axial median over 3 planes also removes them (3D settings).\n"
            "- **Otsu** splits a two-peaked histogram (expressing vs non-expressing tissue). A cell holds thousands of background pixels and a "
            "handful of spot pixels, so Otsu splits the background noise and rejects nothing: it finds every locus but also picks hot pixels "
            "and cosmic rays (see part 4).\n"
            "- **Gaussian fit.** The paper's integer pixel positions are fine at their scale but would cost ±65 nm here.\n"
            "- Defaults were chosen by comparing variants against the simulation's truth (part 4; `scripts/compare_localizers.py`) on a "
            "library run with probes, MS2 and all image errors, a toy run and a 3D run. One seed each: re-check on more seeds and, above "
            "all, on real movies before trusting the numbers.")
    return {**p, "method": method}


def _settings(p: dict, has_cal: bool) -> tuple[dict, bool]:
    keys = list(SP.DEFAULTS) + ["dilate_px"]
    if st.session_state.get("spt_params_run") != st.session_state.get("_run_path"):     # new run: start from its own settings
        st.session_state["spt_params"] = {k: p[k] for k in keys}
        st.session_state["spt_params_run"] = st.session_state.get("_run_path")
    cur = st.session_state["spt_params"]
    with st.expander("⚙ Settings of the spot finder (start from this run's; change and press Apply)", expanded=False):
        with st.form("spt_form"):
            c1, c2, c3, c4 = st.columns(4)
            new = dict(cur)
            new["median_px"] = c1.select_slider("Median filter (px)", [1, 3, 5, 7], int(cur["median_px"]), help="1 = off. The paper uses 7 for objects hundreds of px wide.")
            new["background"] = c2.radio("Background", ["local", "cell"], ["local", "cell"].index(cur["background"]), horizontal=True,
                                         help="local: wider median subtracted; cell: one value per cell and frame")
            new["background_px"] = c3.slider("Background window (px)", 5, 21, int(cur["background_px"]), 2)
            new["smooth_sigma_px"] = c4.slider("Gaussian σ (px)", 0.0, 3.0, float(cur["smooth_sigma_px"]), 0.25)
            new["threshold"] = c1.radio("Threshold", ["noise", "otsu"], ["noise", "otsu"].index(cur["threshold"] if cur["threshold"] in ("noise", "otsu") else "noise"),
                                        horizontal=True, help="noise: min_snr × robust noise; otsu: the paper's")
            new["threshold_on"] = c2.radio("Threshold applied to", ["filtered", "smoothed"], ["filtered", "smoothed"].index(cur["threshold_on"]), horizontal=True)
            new["threshold_scope"] = c3.radio("One threshold per", ["movie", "frame"], ["movie", "frame"].index(cur["threshold_scope"]), horizontal=True)
            new["min_snr"] = c4.slider("min_snr", 2.0, 10.0, float(cur["min_snr"]), 0.5)
            new["min_area_px"] = c1.slider("Smallest area (px)", 1, 15, int(cur["min_area_px"]))
            new["min_sep_px"] = c2.slider("Min distance between maxima (px)", 1, 15, int(cur["min_sep_px"]))
            new["select"] = c3.radio("Choice", ["track", "brightest"], ["track", "brightest"].index(cur["select"]), horizontal=True)
            new["refine"] = c4.radio("Position", ["gaussian", "centroid", "none"], ["gaussian", "centroid", "none"].index(cur["refine"]), horizontal=True)
            new["dilate_px"] = c1.slider("Grow the cell mask by (px)", 0, 8, int(cur["dilate_px"]))
            new["track_gate_px"] = c2.slider("Track gate (px)", 2.0, 20.0, float(cur["track_gate_px"]), 1.0, help="choice = track: max distance from the neighbours' median position")
            if st.form_submit_button("Apply", type="primary"):
                st.session_state["spt_params"] = new
                cur = new
        if st.button("Reset to this run's settings"):
            st.session_state["spt_params"] = {k: p[k] for k in keys}
            cur = st.session_state["spt_params"]
    use_cal = st.checkbox("Apply the run's calibration (dark, flat-field, hot pixels, chromatic shift)", has_cal, key="spt_cal",
                          disabled=not has_cal, help="Only for runs made with 'Simulate and apply a calibration'.")
    return cur, use_cal and has_cal


def _part_steps(run: Path, cfg: dict, p_run: dict, truth: pd.DataFrame | None):
    st.markdown("### 2 · Every step on one cell")
    sdir = run / "stage5_cells" / "truth"
    summ = pd.read_csv(sdir / "cell_summary.csv")
    px_um = cfg["optics"]["pixel_size_nm"] / 1000.0
    cal_dir = run / "stage3_microscopy" / "calibration"
    params, use_cal = _settings(p_run, (cal_dir / "calibration.json").exists())
    calib = str(cal_dir) if use_cal else None
    dil = int(params.get("dilate_px", 2))
    params = {k: v for k, v in params.items() if k != "dilate_px"}
    c1, c2 = st.columns([2, 3])
    cid = c1.selectbox("Cell (true identity)", summ.cell_id.tolist(), key="spt_cell")
    first = sdir / "cells" / f"cell_{int(summ.cell_id.iloc[0]):04d}"
    chans = sorted(q.stem for q in first.glob("locus*.tif"))
    name = c2.radio("Colour", chans, format_func=lambda c: CH_SHORT.get(c, c), horizontal=True, key="spt_ch")
    df, cands = _analyse(str(sdir), int(cid), name, json.dumps(params, sort_keys=True), calib, dil)
    lid = int(name[5:])
    if truth is not None:
        tr = truth[(truth.cell_id == cid) & (truth.locus_id == lid)].set_index("t")
        df["y_true"] = tr.reindex(df.t).y_um.to_numpy() / px_um - df.crop_y0
        df["x_true"] = tr.reindex(df.t).x_um.to_numpy() / px_um - df.crop_x0
        df["visible"] = (tr.reindex(df.t).photons.to_numpy() > 0) if "photons" in tr else True
        df["photons"] = tr.reindex(df.t).photons.to_numpy() if "photons" in tr else np.nan
        df["err_nm"] = 1000 * px_um * np.hypot(df.y_crop - df.y_true, df.x_crop - df.x_true)
    n_vis = int(df.visible.sum()) if "visible" in df else len(df)
    found = df[df.found]
    far = (found.err_nm > 1000 * FAR_UM).sum() if "err_nm" in df else 0
    style.stat_grid(st, [("frames", f"{len(df)}", "frames this cell exists in"),
                         ("found", f"{len(found)} of {n_vis}", "frames with a locus found / frames where a spot was drawn"),
                         ("wrong object", f"{far}", f"found more than {FAR_UM} µm from the truth (hot pixel, cosmic ray, noise peak)"),
                         ("median error", f"{found.err_nm.median():.0f} nm" if "err_nm" in df and len(found) else "–", "typical distance to the truth"),
                         ("threshold used", f"{df.thr.iloc[0]:.1f} counts" if np.isfinite(df.thr.iloc[0]) else "–", "above background"),
                         ("frames with >1 candidate", f"{(df.n_candidates > 1).sum()}", "where the choice mattered")])
    # jump list: problem frames first
    probs = []
    if "err_nm" in df:
        for r in df[df.found].sort_values("err_nm", ascending=False).head(8).itertuples():
            probs.append((r.i, f"frame {r.t}: error {r.err_nm:.0f} nm, {r.n_candidates} candidate(s)"))
        for r in df[(~df.found) & df.visible].head(8).itertuples():
            probs.append((r.i, f"frame {r.t}: missed ({r.photons:.0f} photons)"))
    if "spt_i" not in st.session_state or st.session_state["spt_i"] >= len(df):
        st.session_state["spt_i"] = 0
    if probs:
        def _jump():
            st.session_state["spt_i"] = int(st.session_state["spt_jump"][0])
        st.selectbox("Jump to a problem frame of this cell", probs, format_func=lambda x: x[1], key="spt_jump", on_change=_jump, index=None,
                     placeholder="largest errors and missed spots…")
    i = st.slider("Frame (of this cell)", 0, len(df) - 1, key="spt_i")
    row = df.iloc[i]
    res, region = _frame_steps(str(sdir), int(cid), name, params, calib, dil, i, float(row.thr))
    if res is None:
        st.info("The cell region is too small in this frame.")
        return
    _steps_figure(res, region, row, cands[i], params)
    _time_figure(df, i)


def _steps_figure(res: dict, region: np.ndarray, row, cand: np.ndarray | None, params: dict):
    s = res["steps"]
    raw, med, bg, sub, sm, area = s["raw"], s["med"], s["bg"], s["sub"], s["smooth"], s["area"]
    ys, xs = np.nonzero(region)
    y0, y1, x0, x1 = max(ys.min() - 4, 0), min(ys.max() + 5, raw.shape[0]), max(xs.min() - 4, 0), min(xs.max() + 5, raw.shape[1])
    lo, hi = np.percentile(raw[region], 1), np.percentile(raw[region], 99.9)
    ref = (row.y_true, row.x_true) if "y_true" in row and np.isfinite(row.y_true) else ((row.y_crop, row.x_crop) if row.found else None)
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.6))
    bgimg = bg if np.ndim(bg) else np.full(raw.shape, float(bg))
    panels = [("① recorded", raw, (lo, hi)), ("① median filtered", med, (lo, hi)), ("background (wider median)", bgimg, (lo, hi)),
              ("⑤ smoothed, ④ area, ⑥ candidates", sm, (0, max(np.percentile(sm[region], 99.9), 1e-9)))]
    for ax, (title, im, (a, b)) in zip(axes[0], panels):
        ax.imshow(im[y0:y1, x0:x1], cmap="gray", vmin=a, vmax=b, extent=(x0 - 0.5, x1 - 0.5, y1 - 0.5, y0 - 0.5))
        ax.contour(np.arange(x0, x1), np.arange(y0, y1), region[y0:y1, x0:x1], levels=[0.5], colors="#eab308", linewidths=0.8)
        ax.set_title(title, fontsize=13)
        ax.axis("off")
    axa = axes[0, 3]
    if area.any():
        axa.contour(np.arange(x0, x1), np.arange(y0, y1), area[y0:y1, x0:x1], levels=[0.5], colors="#22d3ee", linewidths=1.0)
    for k, pk in enumerate(cand if cand is not None else []):
        chosen = row.chosen == k
        axa.add_patch(Circle((pk[1], pk[0]), 5, fill=False, ec="#ef4444" if chosen else "#f59e0b", lw=2 if chosen else 1.2))
        axa.annotate(str(k + 1), (pk[1] + 5, pk[0] - 5), color="#ef4444" if chosen else "#f59e0b", fontsize=11, weight="bold")
    if ref is not None:
        h = 10
        cy, cx = int(round(ref[0])), int(round(ref[1]))
        for ax in axes[0]:
            ax.add_patch(Rectangle((cx - h - 0.5, cy - h - 0.5), 2 * h + 1, 2 * h + 1, fill=False, ec="#a78bfa", lw=1.2))
        zy0, zx0 = max(cy - h, 0), max(cx - h, 0)
        zy1, zx1 = min(cy + h + 1, raw.shape[0]), min(cx + h + 1, raw.shape[1])
        ext = (zx0 - 0.5, zx1 - 0.5, zy1 - 0.5, zy0 - 0.5)
        zl, zh = np.percentile(raw[zy0:zy1, zx0:zx1], 1), raw[zy0:zy1, zx0:zx1].max()
        zpanels = [("zoom: recorded", raw, (zl, zh)), ("zoom: median − background", sub, (0, max(sub[zy0:zy1, zx0:zx1].max(), 1e-9))),
                   ("zoom: smoothed and threshold", sm, (0, max(sm[zy0:zy1, zx0:zx1].max(), 1e-9))), ("⑦ position (recorded pixels)", raw, (zl, zh))]
        for ax, (title, im, (a, b)) in zip(axes[1], zpanels):
            ax.imshow(im[zy0:zy1, zx0:zx1], cmap="gray", vmin=a, vmax=b, extent=ext, interpolation="nearest")
            ax.set_title(title, fontsize=13)
            ax.axis("off")
        thr_img = sm if params.get("threshold_on") == "smoothed" else sub
        axes[1, 2].contour(np.arange(zx0, zx1), np.arange(zy0, zy1), thr_img[zy0:zy1, zx0:zx1], levels=[s["thr"]], colors="#22d3ee", linewidths=1.2)
        axp = axes[1, 3]
        if "y_true" in row and np.isfinite(row.y_true):
            axp.plot(row.x_true, row.y_true, "o", mfc="none", mec="cyan", ms=16, mew=2, label="true position")
        if row.found:
            axp.plot(row.x_crop, row.y_crop, "+", color="#ef4444", ms=18, mew=2.5, label="found")
            fh = int(params.get("fit_half_px", 3))
            pk = cand[int(row.chosen)] if cand is not None and row.chosen is not None and not pd.isna(row.chosen) else None
            if pk is not None:
                axp.add_patch(Rectangle((pk[1] - fh - 0.5, pk[0] - fh - 0.5), 2 * fh + 1, 2 * fh + 1, fill=False, ec="#f59e0b", lw=1.2, ls="--"))
    else:
        for ax in axes[1]:
            ax.axis("off")
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    bgv = float(np.median(bgimg[region]))
    style.stat_grid(st, [("threshold", f"{s['thr']:.1f}", "counts above background, the same in every frame when one threshold per movie"),
                         ("noise (filtered / smoothed)", f"{s['noise']:.2f} / {s['noise_smooth']:.2f}", "robust noise inside the cell region"),
                         ("background", f"{bgv:.1f}", "median background inside the region"),
                         ("candidates", f"{len(s['candidates'])}", "local maxima inside the area; red circle = chosen"),
                         ("SNR of the choice", f"{row.snr:.1f}" if row.found else "–", "smoothed peak / smoothed noise"),
                         ("error", f"{row.err_nm:.0f} nm" if "err_nm" in row and row.found else "–", "distance to the true position"),
                         ("fitted σ", f"{row.sigma_fit_px:.2f} px" if row.found and np.isfinite(row.sigma_fit_px) else "–", "width of the fitted Gaussian")])
    st.caption("Top: the whole cell (yellow = its mask grown by the dilation). Purple square = the zoom below. Cyan outline = the area above the "
               "threshold; numbered circles = candidates (red = the one chosen as the locus). Bottom right: cyan ○ = true position (drift-corrected), "
               "red + = found (after removing the measured chromatic shift if the calibration is applied), dashed square = the fit window.")


def _time_figure(df: pd.DataFrame, i: int):
    if "err_nm" not in df:
        return
    fig, axes = plt.subplots(1, 2, figsize=(17, 3.6))
    f = df[df.found]
    axes[0].semilogy(f.t, f.err_nm.clip(lower=1), ".", color="#0f766e", ms=6)
    far = f[f.err_nm > 1000 * FAR_UM]
    axes[0].semilogy(far.t, far.err_nm, "o", mfc="none", mec="#b91c1c", ms=9, label="wrong object (> 0.5 µm)")
    miss = df[(~df.found) & df.visible]
    if len(miss):
        axes[0].plot(miss.t, np.full(len(miss), 2.0), "|", color="#b45309", ms=14, mew=2, label="missed")
    axes[0].set_title("Error per frame (nm)")
    axes[0].legend(fontsize=10)
    axes[1].plot(df.t, df.n_candidates, drawstyle="steps-mid", color="#64748b", label="candidates")
    ax2 = axes[1].twinx()
    ax2.plot(df.t, df.snr, ".", color="#2563eb", ms=5)
    ax2.set_ylabel("SNR", color="#2563eb")
    ax2.grid(False)
    axes[1].set_title("Candidates per frame (grey) and SNR of the choice (blue)")
    for ax in axes:
        ax.axvline(df.t.iloc[i], color="orange", lw=1.2)
        ax.set_xlabel("frame")
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)


def _part_scores(run: Path):
    st.markdown("### 3 · This run's result")
    sj = run / "stage6_analysis" / "scores.json"
    if not sj.exists():
        st.info("No scores for this run.")
        return
    sc = json.loads(sj.read_text())
    rows = []
    for key, label in (("truth", "2D"), ("truth_calibrated", "2D, calibrated"), ("truth_3d", "3D"), ("truth_3d_calibrated", "3D, calibrated")):
        s = sc.get(key) or {}
        if "median_error_um" in s:
            rows.append({"analysis": label, "found (of visible)": s.get("detection_rate"), "wrong object": s.get("frac_far"),
                         "median error (nm)": 1000 * s["median_error_um"], "95th pct (nm)": 1000 * s["p95_error_um"],
                         "rms of good ones (nm)": 1000 * s["rms_error_near_um"] if s.get("rms_error_near_um") is not None else None,
                         "spots with no light": s.get("n_false_on_dark")})
    if rows:
        st.dataframe(pd.DataFrame(rows).round(3), hide_index=True, width="stretch")
        st.caption(f"'Wrong object' = share of detections more than {FAR_UM} µm from the truth (the finder picked a hot pixel, cosmic ray or noise "
                   "peak). 'rms of good ones' leaves those out. Runs made before this tab existed lack some columns.")
    ef = run / "stage6_analysis" / "truth" / "localization_errors.csv"
    if ef.exists():
        e = pd.read_csv(ef)
        fig, ax = plt.subplots(figsize=(17, 3.4))
        ax.hist(np.log10(1000 * e.err_um.clip(lower=1e-3)), bins=60, color="#0f766e")
        ax.axvline(np.log10(1000 * FAR_UM), color="#b91c1c", lw=1.2)
        ax.set_xlabel("error (nm, log scale: 1 = 10 nm, 2 = 100 nm, 3 = 1 µm)")
        ax.set_title("All 2D detections of the run (true cell IDs)")
        st.pyplot(fig, width="stretch")
        plt.close(fig)


def _part_compare(run: Path):
    st.markdown("### 4 · Compare methods on this run")
    st.caption("Re-runs only the spot finding (stage 6) on this run's saved single-cell crops with each variant below and scores it against the truth. "
               "Takes about 1 to 3 minutes for a 2D run (longer with z-stacks). One run is one seed: compare several runs before concluding.")
    mc = run / "stage6_analysis" / "method_comparison.csv"
    c1, c2 = st.columns([1, 3])
    no3d = c2.checkbox("Skip the 3D variants (much faster)", True, key="spt_no3d")
    if c1.button("▶ Compare methods", type="primary"):
        cmd = [sys.executable, str(REPO_ROOT / "scripts" / "compare_localizers.py"), "--run", str(run)] + (["--no-3d"] if no3d else [])
        with st.status("Comparing methods…", expanded=True) as status:
            proc = subprocess.Popen(cmd, cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                                    env={**os.environ, "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"})
            log = []
            for line in proc.stdout:
                log.append(line.rstrip())
                if "detected" in line:
                    st.write(line.strip())
            if proc.wait() == 0:
                status.update(label="Done", state="complete", expanded=False)
            else:
                status.update(label="Failed", state="error")
                st.code("\n".join(log[-30:]))
    if not mc.exists():
        st.info("No comparison yet for this run.")
        return
    tab = pd.read_csv(mc)
    labels = {**{k: v[0] for k, v in VARIANTS.items()}, **{f"{k} (3D)": v[0] for k, v in VARIANTS_3D.items()}}
    tab["name"] = [f"{v} (3D)" if d == "3d" else v for v, d in zip(tab.variant, tab.dim)]
    tab["what"] = tab.name.map(labels).fillna(tab.label)
    show = pd.DataFrame({"variant": tab.name, "calibrated": tab.calibrated, "found": tab.detection_rate, "wrong object": tab.frac_far,
                         "median (nm)": 1000 * tab.median_error_um, "95th pct (nm)": 1000 * tab.p95_error_um,
                         "rms of good ones (nm)": 1000 * tab.rms_error_near_um, "what changes": tab.what})
    st.dataframe(show.round(3), hide_index=True, width="stretch")
    for cal in sorted(tab.calibrated.unique()):
        d = tab[(tab.calibrated == cal) & (tab.dim == "2d")]
        if d.empty:
            continue
        fig, axes = plt.subplots(1, 3, figsize=(17, 0.45 * len(d) + 1.6), sharey=True)
        y = np.arange(len(d))
        for ax, col, title, scale in ((axes[0], "detection_rate", "found (of visible spots)", 1), (axes[1], "frac_far", "wrong object (share)", 1),
                                      (axes[2], "median_error_um", "median error (nm)", 1000)):
            ax.barh(y, scale * d[col].fillna(0), color=["#0f766e" if v == "default" else "#94a3b8" for v in d.variant])
            ax.set_title(title, fontsize=13)
        axes[0].set_yticks(y)
        axes[0].set_yticklabels(d.variant)
        axes[0].invert_yaxis()
        fig.suptitle(f"2D, {'calibration applied' if cal else 'no calibration'} (green = current default)", fontsize=13, x=0.01, ha="left")
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)


def render(st_, run_path: Path, R: dict):
    style.intro(st, "How stage 6 finds each locus in each cell: the steps (adapted from Bennabi et al. 2025), every intermediate image for "
                    "any cell and frame, the error of every frame against the simulation's truth, and a comparison of methods.",
                "Part 2: pick a cell, a colour and a frame (or jump to a problem frame), and change the settings to see what each step "
                "does. Part 4: press Compare methods to score every variant on this whole run.")
    run = Path(run_path)
    st.session_state["_run_path"] = str(run)
    cfg = R["cfg"]
    if not cfg.get("loci") or not (run / "stage5_cells" / "truth" / "cell_summary.csv").exists():
        st.info("This run has no loci analysed by true cell identity. Make one with 🔬 New single-cell loci run (sidebar).")
        return
    stamp = (run / "stage7_validation" / "metrics.json").stat().st_mtime
    p_run = _part_method(cfg)
    try:
        truth = _truth(str(run), stamp)
    except Exception as e:  # noqa: BLE001 - old runs may lack files; the steps still work without the truth
        st.warning(f"Could not load the true locus positions ({e}); errors are not shown.")
        truth = None
    _part_steps(run, cfg, p_run, truth)
    _part_scores(run)
    _part_compare(run)
