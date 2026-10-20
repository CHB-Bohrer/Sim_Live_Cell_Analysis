"""The dashboard's 🔎 Spot finding tab: how stage 6 finds each locus, step by step, and how well it does against the truth.

Four parts: (1) the method and where it comes from (Chen, Levo, Zoller & Gregor 2025);
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
import yaml
from matplotlib.patches import Circle, Rectangle

import style
from simlive.io.runs import DATA_ROOT, REPO_ROOT
from simlive.stage5_linking.celldata import load_cell
from simlive.stage6_analysis import dog_viterbi as DV
from simlive.stage6_analysis import spots as SP
from simlive.stage6_analysis.compare import VARIANTS, VARIANTS_3D, truth_with_photons
from simlive.stage6_analysis import multicopy as MC
from simlive.stage6_analysis.localize import analyse_cell_copies, analyse_channel, channel_inputs, copy_params, match_copies, spot_params
from simlive.stage6_analysis.ms2 import ERODE_PX, ms2_inputs, ms2_truth, relabel_to_true_copies, trace_from_analysis

CHEN_URL = "https://doi.org/10.1038/s41594-025-01615-4"
CHEN_ARXIV = "https://arxiv.org/abs/2304.08770"
CH_SHORT = {"locus0": "locus 0 (promoter)", "locus1": "locus 1 (enhancer)", "nucleus": "MS2 spot (nuclear channel)"}
FAR_UM = 0.5
METHOD_LABEL = {"median_peaks": "median_peaks (removed method)", "dog_viterbi": "dog_viterbi (Chen et al. 2025)",
                "brightest": "brightest (original placeholder)"}

STEPS_CHEN = [
    ("1 One stack per nucleus", "nuclei segmented (DoG + watershed) and tracked (Voronoi overlap); one 3D stack per nucleus over time",
     "stage 5's single-cell crops (true or tracked cell IDs), the cell mask grown by dilate_px", "dilate_px"),
    ("2 Median filter", "3D median filter against salt-and-pepper noise", "median_px × median_px (× median_z planes in z-stacks)", "median_px"),
    ("3 Difference of Gaussians", "3D DoG for round objects, tuned to high sensitivity: several candidates per nucleus",
     "DoG with σ = dog_sigma_px (inner ≈ PSF, outer 3×); candidates = local maxima above dog_snr × noise, at most max_candidates", "dog_snr"),
    ("4 Trajectory", "modified Viterbi: at most one spot per nucleus; links scored by distance with a diffusion potential built "
                     "iteratively from the random walk of all spots",
     "Viterbi over the whole movie: link cost = −log of the diffusion propagator, reward = log SNR + log(area / false candidates per frame), "
     "frames may have no spot "
     "(viterbi_miss_cost each, gaps up to viterbi_max_gap); D re-estimated from the path itself, per cell", "viterbi_miss_cost"),
    ("5 Position", "centroid of the spot", "centroid, or the Gaussian fit (refine)", "refine"),
    ("6 Intensity", "mean in a (7, 7, 5) px ellipsoid minus mean of the (19, 19, 7) px shell around it (removes free MCP-GFP)",
     "the same ellipsoids, for every method (column intensity; compared with the true photons); the ring's median by default "
     "(intensity_background: mean = the paper)", "intensity_background"),
]

# ----------------------------------------------------------------------------- cached computations
@st.cache_data(show_spinner=False, max_entries=4)
def _truth(run: str, stamp: float) -> pd.DataFrame:
    return truth_with_photons(Path(run))


@st.cache_resource(show_spinner="Loading the cell…", max_entries=4)
def _inputs(sdir: str, cid: int, name: str, calib: str | None, dilate_px: int) -> dict:
    """The crops of one colour of one cell as the localizer sees them (calibration-corrected if asked). "nucleus" = the MS2 spot: the nuclear
    channel inside the nucleus shrunk by ms2.ERODE_PX (its bright edge would give candidates)."""
    if name == "nucleus":
        return {**ms2_inputs(load_cell(sdir, cid), ERODE_PX, calib), "shift": (0.0, 0.0)}
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
                     "snr": f["snr"] if f else np.nan, "fit_ok": f["fit_ok"] if f else None, "sigma_fit_px": f["sigma_fit_px"] if f else np.nan,
                     "intensity": f["intensity"] if f else np.nan})
    cands = [None if c is None else c["peaks"] for c in a["candidates"]]
    df = pd.DataFrame(rows)
    _, pos, inten = trace_from_analysis(a, ci["stack"], ci["centre"], spot_params(json.loads(params)))
    df["trace_y"], df["trace_x"], df["trace_intensity"] = pos[:, 0], pos[:, 1], inten      # every frame, interpolated where nothing was found
    return df, cands, a.get("D_px2")


def _rows(fr, fits, choice, cands, shift) -> pd.DataFrame:
    """One row per frame of one spot (as _analyse), from fits / chosen candidate indices."""
    sy, sx = shift
    rows = []
    for i, r in enumerate(fr.itertuples()):
        f, c = fits[i], cands[i]
        rows.append({"i": i, "t": int(r.t), "crop_y0": int(r.crop_y0), "crop_x0": int(r.crop_x0), "thr_frame": np.nan, "thr": np.nan,
                     "n_candidates": 0 if c is None else len(c["peaks"]), "chosen": choice[i],
                     "found": f is not None, "y_crop": f["pos"][0] - sy if f else np.nan, "x_crop": f["pos"][1] - sx if f else np.nan,
                     "snr": f["snr"] if f else np.nan, "fit_ok": f["fit_ok"] if f else None, "sigma_fit_px": f["sigma_fit_px"] if f else np.nan,
                     "intensity": f["intensity"] if f else np.nan})
    return pd.DataFrame(rows)


@st.cache_data(show_spinner="Running the multi-copy linker on this cell (all colours together)…", max_entries=8)
def _analyse_copies(sdir: str, cid: int, params: str, calib: str | None, dilate_px: int, px_um: float, n_loci: int):
    """The multi-copy linker (multicopy.py) on all locus colours of one cell. Returns {(colour, found copy): rows per frame}, the candidates of
    every colour, and the paths it tried with their features."""
    kw = json.loads(params)
    ins = [_inputs(sdir, cid, f"locus{k}", calib, dilate_px) for k in range(n_loci)]
    a = analyse_cell_copies(ins, px_um=px_um, **kw)
    fr = ins[0]["frames"]
    tr = a["track"]
    dfs = {}
    for c, cp in enumerate(tr["copies"]):
        for k, ch in cp.items():
            df = _rows(fr, a["copies"][c][k], ch, a["candidates"][k], ins[k]["shift"])
            thr = [cc["values"].min() if cc is not None and len(cc["values"]) else np.nan for cc in a["candidates"][k]]
            df["thr"] = df["thr_frame"] = np.nanmin(thr) if np.isfinite(thr).any() else np.nan
            dfs[(k, c)] = df
    cands = {k: [None if cc is None else cc["peaks"] for cc in a["candidates"][k]] for k in range(n_loci)}
    paths = {k: [{kk: v for kk, v in p.items()} for p in tr["paths"][k]] for k in tr["paths"]}
    used = {k: {} for k in range(n_loci)}
    for c, cp in enumerate(tr["copies"]):
        for k in cp:
            used[k][c] = cp[k]
    return dfs, cands, paths, tr["pairs"], used, fr[["t", "crop_y0", "crop_x0"]].reset_index(drop=True)


def _frame_steps(sdir: str, cid: int, name: str, params: dict, calib: str | None, dilate_px: int, i: int, thr: float):
    """Every intermediate image of one frame, in one form for both methods: raw, med, bg, sub (filtered minus background / DoG), smooth
    (the image the candidates are taken from), area (above the threshold), thr, noise, noise_smooth, candidates."""
    ci = _inputs(sdir, cid, name, calib, dilate_px)
    img, reg = ci["stack"][i], ci["region"][i]
    c = DV.dog_candidates(img, reg, steps=True, **params)
    if c is None:
        return None, reg
    dog = c["pre"]["sub"]
    steps = {"raw": c["img"], "med": c["pre"]["med"], "bg": c["pre"]["bg"], "sub": dog, "smooth": dog, "area": (dog > c["thr"]) & reg,
             "thr": c["thr"], "noise": c["pre"]["noise_smooth"], "noise_smooth": c["pre"]["noise_smooth"], "candidates": c["peaks"]}
    return {"found": len(c["peaks"]) > 0, "steps": steps}, reg


# ----------------------------------------------------------------------------- parts
def _part_method(cfg: dict) -> dict:
    st.markdown("### 1 · The method and where it comes from")
    la = cfg.get("locus_analysis") or {}
    method = la.get("method", "brightest" if "median_px" not in la else "median_peaks")
    p = {**SP.DEFAULTS, "dilate_px": 2, **la}
    st.markdown(f"This run located the loci with **{METHOD_LABEL.get(method, method)}**."
                + (" That method (steps adapted from Bennabi et al. 2025, who detect Mesp2 poles, not spots) has been removed; part 2 shows "
                   "the current method on this run's cells." if method == "median_peaks" else ""))
    with st.expander("dog_viterbi: the MS2 spot pipeline of Chen, Levo, Zoller & Gregor 2025 (transcription sites in fly embryos)",
                     expanded=True):
        st.markdown(f"[Chen et al. 2025, *Nat. Struct. Mol. Biol.*]({CHEN_URL}) (methods from the [arXiv version]({CHEN_ARXIV}), section III.2; no "
                    "code was published, so this follows the text) find faint MS2 transcription spots, one labelled allele per nucleus, on the "
                    "background of free MCP-GFP in the nucleus: the same situation as a locus (or the MS2 spot) in one of our cells. Their 220 nm "
                    "pixels and two-photon PSF are close to our 130 nm pixels and PSF in units of the PSF, so their sizes are kept.")
        lines = ["| step | Chen et al. 2025 (MS2 spots) | here | this run |", "|---|---|---|---|"]
        lines += [f"| **{s}** | {a} | {b} | " + (f"`{k}: {p.get(k)}`" if method == "dog_viterbi" else "—") + " |" for s, a, b, k in STEPS_CHEN]
        st.markdown("\n".join(lines))
        st.caption("What the text leaves open and was chosen here: the DoG widths, the threshold (3 × noise), the Viterbi score (link = "
                   "−log diffusion propagator; reward = log SNR plus the clutter term log(area / false candidates per frame), without which "
                   "the path skipped ~30% of frames holding the right spot; a cost per frame without a spot) and estimating the diffusion "
                   "potential per cell rather than pooled over all spots of the movie. Details: src/simlive/stage6_analysis/dog_viterbi.py "
                   "and docs/SPOT_FINDING.md.")
    st.caption("Settings were checked by comparing variants against the simulation's truth (part 4; scripts/compare_localizers.py) on library, "
               "toy and 3D runs. One seed each: re-check on more seeds and, above all, on real movies before trusting the numbers.")
    return {**p, "method": method}


def _settings(p: dict, has_cal: bool) -> tuple[dict, bool]:
    keys = list(SP.DEFAULTS) + ["dilate_px"]

    def run_settings():
        return {**{k: p[k] for k in keys}, "method": "dog_viterbi"}

    if st.session_state.get("spt_params_run") != st.session_state.get("_run_path"):     # new run: start from its own settings
        st.session_state["spt_params"] = run_settings()
        st.session_state["spt_params_run"] = st.session_state.get("_run_path")
    cur = st.session_state["spt_params"]
    with st.expander("⚙ Settings of the spot finder (start from this run's; change and press Apply)", expanded=False):
        with st.form("spt_form"):
            new = dict(cur)
            st.markdown("**dog_viterbi (Chen et al. 2025)**")
            d1, d2, d3, d4 = st.columns(4)
            ds = list(cur["dog_sigma_px"])
            new["dog_sigma_px"] = [d1.slider("DoG inner σ (px)", 0.5, 3.0, float(ds[0]), 0.25),
                                   d1.slider("DoG outer σ (px)", 1.0, 8.0, float(ds[1]), 0.5)]
            new["dog_snr"] = d2.slider("Candidate threshold (× noise)", 1.5, 10.0, float(cur["dog_snr"]), 0.5,
                                       help="Chen et al. tune it low (high sensitivity): the Viterbi path then picks the real spot")
            new["max_candidates"] = d2.slider("Candidates kept per frame", 1, 15, int(cur["max_candidates"]))
            new["viterbi_miss_cost"] = d3.slider("Cost of a frame without a spot", 0.0, 10.0, float(cur["viterbi_miss_cost"]), 0.5,
                                                 help="Higher = the path takes a weak or distant candidate rather than skip a frame")
            new["viterbi_max_gap"] = d3.slider("Longest gap (frames)", 1, 20, int(cur["viterbi_max_gap"]))
            new["intensity_inner_px"] = [d4.slider("Intensity ellipse (px)", 3, 15, int(cur["intensity_inner_px"][0]), 2)] * 2 + [int(cur["intensity_inner_px"][2])]
            new["intensity_outer_px"] = [d4.slider("Background ellipse (px)", 9, 31, int(cur["intensity_outer_px"][0]), 2)] * 2 + [int(cur["intensity_outer_px"][2])]
            new["intensity_background"] = d4.radio("Background of the ring", ["median", "mean"], ["median", "mean"].index(cur["intensity_background"]),
                                                   horizontal=True, help="mean = Chen et al.; median = robust to a cosmic ray or another spot in the ring")
            c1, c2, c3, c4 = st.columns(4)
            new["median_px"] = c1.select_slider("Median filter (px)", [1, 3, 5], int(cur["median_px"]), help="1 = off")
            new["min_sep_px"] = c2.slider("Min distance between candidates (px)", 1, 15, int(cur["min_sep_px"]))
            new["refine"] = c3.radio("Position", ["gaussian", "centroid", "none"], ["gaussian", "centroid", "none"].index(cur["refine"]), horizontal=True,
                                     help="centroid = as in the paper")
            new["dilate_px"] = c4.slider("Grow the cell mask by (px)", 0, 8, int(cur["dilate_px"]))
            if st.form_submit_button("Apply", type="primary"):
                st.session_state["spt_params"] = new
                cur = new
        if st.button("Reset to this run's settings"):
            st.session_state["spt_params"] = run_settings()
            cur = st.session_state["spt_params"]
    use_cal = st.checkbox("Apply the run's calibration (dark, flat-field, hot pixels, chromatic shift)", has_cal, key="spt_cal",
                          disabled=not has_cal, help="Only for runs made with 'Simulate and apply a calibration'.")
    return cur, use_cal and has_cal


@st.cache_data(show_spinner=False, max_entries=4)
def _ms2_truth(run: str, stamp: float) -> pd.DataFrame | None:
    t = ms2_truth(Path(run))
    return None if t is None else t.assign(locus_id=-1)


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
    mt = _ms2_truth(str(run), (run / "stage1_chromatin").stat().st_mtime) if (run / "stage1_chromatin" / "transcription.csv").exists() else None
    if mt is not None:
        chans.append("nucleus")
        truth = mt if truth is None else pd.concat([truth, mt], ignore_index=True)
    name = c2.radio("Spot", chans, format_func=lambda c: CH_SHORT.get(c, c), horizontal=True, key="spt_ch",
                    help="The two locus colours, and (runs with the transcription model) the MS2 spot at the promoter, in the nuclear channel")
    lid = -1 if name == "nucleus" else int(name[5:])
    cp = copy_params(p_run)
    multi = MC.multicopy_on(cp) and name != "nucleus"
    n_loci = len([c for c in chans if c.startswith("locus")])
    true_copies = sorted(truth[(truth.cell_id == cid)]["copy"].unique()) if truth is not None and "copy" in truth else [0]
    copy_true = 0
    if multi:                  # the run's own linker: several paths per colour, paired into copies; numbered like the true copies they follow
        dfs, cands_all, paths, pairs, used, frc = _analyse_copies(str(sdir), int(cid), json.dumps({**params, **cp}, sort_keys=True), calib, dil,
                                                                   px_um, n_loci)
        mapping = _copy_mapping(dfs, truth, cid, px_um) if truth is not None else {}
        found = sorted({c for (k, c) in dfs if k == lid})
        if not found:
            st.warning("The linker found no copy of this colour in this cell.")
            return
        lab = lambda c: f"found copy {c}" + (f" = true copy {mapping[c]}" if c in mapping else " (no true copy: a stray or noise)" if truth is not None else "")  # noqa: E731
        cf = st.radio("Copy", found, format_func=lab, horizontal=True, key="spt_copy",
                      help="The cell's promoter-enhancer copies as the linker found them (part 2b shows every path it tried)") if len(found) > 1 else found[0]
        df, cands, D = dfs[(lid, cf)], cands_all[lid], None
        copy_true = mapping.get(cf, -1)
    else:
        df, cands, D = _analyse(str(sdir), int(cid), name, json.dumps(params, sort_keys=True), calib, dil)
        if len(true_copies) > 1:          # the MS2 spot (own track, one path) or a one-copy finder: compare with one chosen true copy
            copy_true = st.radio("Compare with true copy", true_copies, horizontal=True, key="spt_copy_true")
    if name == "nucleus":
        st.info("MS2 spot: one spot per nucleus while the promoter is ON, nothing while it is OFF. 'Found' below counts ON frames; the trace "
                "plot underneath shows the intensity measured at the tracked position in every frame against the true MS2 signal.")
    if truth is not None:
        tq = truth[(truth.cell_id == cid) & (truth.locus_id == lid)]
        tr = (tq[tq["copy"] == copy_true] if "copy" in tq else tq).set_index("t")
        df["y_true"] = tr.reindex(df.t).y_um.to_numpy() / px_um - df.crop_y0
        df["x_true"] = tr.reindex(df.t).x_um.to_numpy() / px_um - df.crop_x0
        df["visible"] = (tr.reindex(df.t).photons.to_numpy() > 0) if "photons" in tr else True
        df["photons"] = tr.reindex(df.t).photons.to_numpy() if "photons" in tr else np.nan
        df["photons_frame"] = tr.reindex(df.t).photons_frame.to_numpy() if "photons_frame" in tr else df["photons"]
        df["err_nm"] = 1000 * px_um * np.hypot(df.y_crop - df.y_true, df.x_crop - df.x_true)
    n_vis = int(df.visible.sum()) if "visible" in df else len(df)
    found = df[df.found & df.visible] if "visible" in df else df[df.found]      # only frames where a spot was really drawn (MS2: ON frames)
    far = (found.err_nm > 1000 * FAR_UM).sum() if "err_nm" in df else 0
    good = found[found.err_nm <= 1000 * FAR_UM] if "err_nm" in df else found
    r_int = (np.corrcoef(good.intensity, good.photons_frame)[0, 1] if "photons_frame" in good and len(good) > 5
             and good.photons_frame.std() > 0 and good.intensity.std() > 0 else np.nan)
    style.stat_grid(st, [("frames", f"{len(df)}", "frames this cell exists in"),
                         ("found", f"{len(found)} of {n_vis}", "frames with a locus found / frames where a spot was drawn"),
                         ("wrong object", f"{far}", f"found more than {FAR_UM} µm from the truth (hot pixel, cosmic ray, noise peak)"),
                         ("median error", f"{found.err_nm.median():.0f} nm" if "err_nm" in df and len(found) else "–", "typical distance to the truth"),
                         ("threshold used", f"{df.thr.iloc[0]:.1f} counts" if np.isfinite(df.thr.iloc[0]) else "–", "above background"),
                         ("frames with >1 candidate", f"{(df.n_candidates > 1).sum()}", "where the choice mattered"),
                         ("intensity vs true photons", f"r = {r_int:.2f}" if np.isfinite(r_int) else "–",
                          "correlation of the ellipsoid intensity with the photons the simulation drew (correct detections)")]
                    + ([("diffusion potential D", f"{D:.2f} px²/frame", "estimated by the Viterbi iterations from the path's own steps "
                         "(relative to the nucleus centroid)")] if D is not None else []))
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
    _time_figure(df, i, cands, px_um)
    if name == "nucleus" and "photons_frame" in df:
        prom = _ms2_promoter(str(sdir), int(cid), json.dumps({**params, **cp}, sort_keys=True), calib, str(run), px_um)
        prom = prom[prom["copy"] == copy_true].drop(columns="copy") if "copy" in prom else prom
        _ms2_trace_figure(df.merge(prom, on="t", how="left"), i, float(((cfg.get("loci") or {}).get("transcription") or {}).get("photons_per_loop_s", 100)))
    if multi:
        _part_copies(run, cid, paths, pairs, used, frc, cands_all, truth, px_um, mapping)


def _copy_mapping(dfs: dict, truth: pd.DataFrame, cid: int, px_um: float) -> dict:
    """found copy -> true copy for one cell (localize.match_copies on the found spots)."""
    rows = [df[df.found].assign(locus_id=k, copy=c, cell_id=cid, lab_y_px=lambda d: d.y_crop + d.crop_y0, lab_x_px=lambda d: d.x_crop + d.crop_x0)
            for (k, c), df in dfs.items()]
    rows = [r for r in rows if len(r)]
    if not rows:
        return {}
    m = match_copies(pd.concat(rows)[["cell_id", "t", "locus_id", "copy", "lab_y_px", "lab_x_px"]], truth[truth.cell_id == cid], px_um)
    return {int(f): int(c) for f, c in m.groupby("copy_found")["copy"].first().items() if c < 50}


def _classify(cands: list, frc: pd.DataFrame, truth: pd.DataFrame | None, strays: pd.DataFrame | None, cid: int, k: int, px_um: float,
              tol_px: float = 3.0) -> list:
    """What every candidate of colour k is, frame by frame: 'copy 0', 'copy 1', ... (within tol_px of that true locus), 'stray', or 'noise'."""
    out = []
    tq = truth[(truth.cell_id == cid) & (truth.locus_id == k)] if truth is not None else None
    sq = strays[(strays.cell_id == cid) & (strays.channel == f"locus{k}")] if strays is not None else None
    for i, pk in enumerate(cands):
        if pk is None or not len(pk):
            out.append([])
            continue
        t = int(frc.t[i])
        y = np.asarray(pk)[:, -2] + frc.crop_y0[i]
        x = np.asarray(pk)[:, -1] + frc.crop_x0[i]
        lab = ["noise"] * len(y)
        best = np.full(len(y), tol_px)
        if sq is not None:
            for q in sq[sq.t == t].itertuples():
                d = np.hypot(y - q.y_drawn_um / px_um, x - q.x_drawn_um / px_um)
                for j in np.nonzero(d < best)[0]:
                    lab[j], best[j] = "stray", d[j]
        if tq is not None:
            for q in tq[tq.t == t].itertuples():
                yy, xx = (q.y_drawn_um, q.x_drawn_um) if hasattr(q, "y_drawn_um") else (q.y_um, q.x_um)
                d = np.hypot(y - yy / px_um, x - xx / px_um)
                for j in np.nonzero(d < best)[0]:
                    lab[j], best[j] = f"copy {int(getattr(q, 'copy', 0))}", d[j]
        out.append(lab)
    return out


def _part_copies(run: Path, cid: int, paths: dict, pairs: list, used: dict, frc: pd.DataFrame, cands: dict, truth, px_um: float, mapping: dict):
    """What the multi-copy linker did in this cell: every path it tried in each colour, coloured by what it followed."""
    st.markdown("### 2b · Copies and stray spots: what the linker did in this cell")
    sp_f = run / "stage3_microscopy" / "stray_spots.csv"
    strays = pd.read_csv(sp_f) if sp_f.exists() else None
    st.markdown("The linker finds several paths per colour (each through the candidates the earlier ones left), keeps those that look like a "
                "locus (most frames covered, about as bright as the best path, as slow), and pairs every promoter path with an enhancer path that "
                "stays within 1.5 µm of it. Each row below is one path; every dot is the candidate it took in that frame, coloured by what that "
                "candidate really was in the simulation. A good result: the kept paths are one solid colour from start to end.")
    cat_col = {"copy 0": "#e5e7eb", "copy 1": "#2dd4bf", "copy 2": "#f59e0b", "copy 3": "#a78bfa", "stray": "#f0abfc", "noise": "#64748b"}
    cls = {k: _classify(cands[k], frc, truth, strays, cid, k, px_um) for k in cands}
    nk = len(paths)
    fig, axes = plt.subplots(1, nk, figsize=(17, 0.55 * (max(len(paths[k]) for k in paths) + len(pairs)) + 1.9), squeeze=False)
    rows = []
    for k in sorted(paths):
        ax = axes[0][k]
        copy_of_path = {pr[0 if k == 0 else 1]: c for c, pr in enumerate(pairs)} if k in (0, 1) else {}
        show = [(f"path {r}", pth["choice"], f"→ copy {copy_of_path[r]}" if r in copy_of_path else "dropped", r in copy_of_path)
                for r, pth in enumerate(paths[k])]
        for c, ch in used[k].items():          # the pair prior re-ran the copy's path: show the result where it differs from the first pass
            r0 = next((r for r, cc in copy_of_path.items() if cc == c), None)
            if r0 is None or paths[k][r0]["choice"] != ch:
                show.append((f"copy {c}, final", ch, "after the pair prior", True))
        for y, (_, choice, tag, kept) in enumerate(show):
            for cat in cat_col:
                ts = [int(frc.t[i]) for i, j in enumerate(choice) if j is not None and cls[k][i] and cls[k][i][j] == cat]
                if ts:
                    ax.plot(ts, np.full(len(ts), y), "|", color=cat_col[cat], ms=11, mew=2.2)
            ax.text(1.01, y, tag, transform=ax.get_yaxis_transform(), va="center", fontsize=9, color="#4ade80" if kept else "#94a3b8")
        for r, pth in enumerate(paths[k]):
            whats = pd.Series([cls[k][i][j] for i, j in enumerate(pth["choice"]) if j is not None and cls[k][i]]).value_counts(normalize=True)
            rows.append({"colour": {0: "promoter (locus 0)", 1: "enhancer (locus 1)"}.get(k, f"locus {k}"), "path": r,
                         "frames covered": pth["coverage"], "median SNR": pth.get("snr_median"), "D (px²/frame)": pth["D_px2"],
                         "used as": f"found copy {copy_of_path[r]}" if r in copy_of_path else "dropped",
                         "followed": ", ".join(f"{w} {100 * v:.0f}%" for w, v in whats.items())})
        ax.set_yticks(range(len(show)))
        ax.set_yticklabels([n for n, *_ in show])
        ax.invert_yaxis()
        ax.set_xlabel("frame")
        ax.set_title({0: "Promoter colour", 1: "Enhancer colour"}.get(k, f"locus {k}"))
    from matplotlib.lines import Line2D

    present = [c for c in cat_col if any(c in row for k in cls for row in cls[k])]
    axes[0][-1].legend([Line2D([], [], color=cat_col[c], lw=3) for c in present], present, fontsize=9, loc="upper left",
                       bbox_to_anchor=(1.12, 1.0), title="candidate was", title_fontsize=9)
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    if rows:
        st.dataframe(pd.DataFrame(rows).round(2), hide_index=True, width="stretch")
    if pairs:
        st.caption("Pairs (promoter path, enhancer path, median distance): " + "; ".join(f"({i}, {j}, {d:.2f} µm)" for i, j, d in pairs)
                   + (". Found copy → true copy: " + ", ".join(f"{a} → {b}" for a, b in mapping.items()) if mapping else "")
                   + ". A dropped path that follows a true copy for most frames would be a missed copy; a kept path that turns pink or grey part "
                     "of the time took a stray or a noise peak in those frames.")


def _steps_figure(res: dict, region: np.ndarray, row, cand: np.ndarray | None, params: dict):
    s = res["steps"]
    raw, med, bg, sub, sm, area = s["raw"], s["med"], s["bg"], s["sub"], s["smooth"], s["area"]
    ys, xs = np.nonzero(region)
    y0, y1, x0, x1 = max(ys.min() - 4, 0), min(ys.max() + 5, raw.shape[0]), max(xs.min() - 4, 0), min(xs.max() + 5, raw.shape[1])
    lo, hi = np.percentile(raw[region], 1), np.percentile(raw[region], 99.9)
    ref = (row.y_true, row.x_true) if "y_true" in row and np.isfinite(row.y_true) else ((row.y_crop, row.x_crop) if row.found else None)
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.6))
    bgimg = bg if np.ndim(bg) else np.full(raw.shape, float(bg))
    dog = params.get("method") == "dog_viterbi"
    t = (["recorded", "② median filtered", "③ outer Gaussian (local background)", "③ DoG, area above threshold, candidates",
          "zoom: DoG", "zoom: DoG and threshold", "⑤ position, ⑥ intensity ellipses"] if dog else
         ["① recorded", "① median filtered", "background (wider median)", "⑤ smoothed, ④ area, ⑥ candidates",
          "zoom: median − background", "zoom: smoothed and threshold", "⑦ position, intensity ellipses"])
    panels = [(t[0], raw, (lo, hi)), (t[1], med, (lo, hi)), (t[2], bgimg, (lo, hi)),
              (t[3], sm, (0, max(np.percentile(sm[region], 99.9), 1e-9)))]
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
        zpanels = [("zoom: recorded", raw, (zl, zh)), (t[4], sub, (0, max(sub[zy0:zy1, zx0:zx1].max(), 1e-9))),
                   (t[5], sm, (0, max(sm[zy0:zy1, zx0:zx1].max(), 1e-9))), (t[6], raw, (zl, zh))]
        for ax, (title, im, (a, b)) in zip(axes[1], zpanels):
            ax.imshow(im[zy0:zy1, zx0:zx1], cmap="gray", vmin=a, vmax=b, extent=ext, interpolation="nearest")
            ax.set_title(title, fontsize=13)
            ax.axis("off")
        thr_img = sm
        axes[1, 2].contour(np.arange(zx0, zx1), np.arange(zy0, zy1), thr_img[zy0:zy1, zx0:zx1], levels=[s["thr"]], colors="#22d3ee", linewidths=1.2)
        axp = axes[1, 3]
        if "y_true" in row and np.isfinite(row.y_true):
            axp.plot(row.x_true, row.y_true, "o", mfc="none", mec="cyan", ms=16, mew=2, label="true position")
        if row.found:
            axp.plot(row.x_crop, row.y_crop, "+", color="#ef4444", ms=18, mew=2.5, label="found")
            fh = int(params.get("fit_half_px", 3))
            pk = cand[int(row.chosen)] if cand is not None and row.chosen is not None and not pd.isna(row.chosen) else None
            if pk is not None and params.get("refine") == "gaussian":
                axp.add_patch(Rectangle((pk[1] - fh - 0.5, pk[0] - fh - 0.5), 2 * fh + 1, 2 * fh + 1, fill=False, ec="#f59e0b", lw=1.2, ls="--"))
            from matplotlib.patches import Ellipse

            for diam, col in ((params.get("intensity_inner_px", (7, 7, 5)), "#22c55e"), (params.get("intensity_outer_px", (19, 19, 7)), "#a78bfa")):
                axp.add_patch(Ellipse((row.x_crop, row.y_crop), float(diam[0]), float(diam[1]), fill=False, ec=col, lw=1.4))
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
                         ("fitted σ", f"{row.sigma_fit_px:.2f} px" if row.found and np.isfinite(row.sigma_fit_px) else "–", "width of the fitted Gaussian"),
                         ("intensity", f"{row.intensity:.1f}" if row.found and np.isfinite(row.intensity) else "–",
                          "mean in the green ellipse minus mean of the ring out to the purple one (Chen et al. 2025)"),
                         ("true photons (this frame)", f"{row.photons_frame:.0f}" if "photons_frame" in row and np.isfinite(row.photons_frame) else "–",
                          "photons the simulation drew for this spot, after bleaching")])
    st.caption("Top: the whole cell (yellow = its mask grown by the dilation). Purple square = the zoom below. Cyan outline = the area above the "
               "threshold; numbered circles = candidates (red = the one chosen as the locus; with dog_viterbi the choice is the Viterbi path over "
               "the whole movie, so it need not be the brightest). Bottom right: cyan ○ = true position (drift-corrected), red + = found (after "
               "removing the measured chromatic shift if the calibration is applied), dashed square = the fit window, green / purple ellipses = "
               "where the intensity and its background are measured.")


def _time_figure(df: pd.DataFrame, i: int, cands: list, px_um: float):
    if "err_nm" not in df:
        return
    fig, axes = plt.subplots(1, 3, figsize=(17, 3.8))
    f = df[df.found]
    axes[0].semilogy(f.t, f.err_nm.clip(lower=1), ".", color="#14b8a6", ms=6)
    far = f[f.err_nm > 1000 * FAR_UM]
    axes[0].semilogy(far.t, far.err_nm, "o", mfc="none", mec="#f87171", ms=9, label="wrong object (> 0.5 µm)")
    miss = df[(~df.found) & df.visible]
    if len(miss):
        axes[0].plot(miss.t, np.full(len(miss), 2.0), "|", color="#f59e0b", ms=14, mew=2, label="missed")
    axes[0].set_title("Error per frame (nm)")
    axes[0].legend(fontsize=10)
    # every candidate's distance from the true spot: shows whether the right one was there and whether the choice took it
    for r, pk in zip(df.itertuples(), cands):
        if pk is None or not len(pk) or not np.isfinite(r.y_true):
            continue
        d = 1000 * px_um * np.hypot(pk[:, -2] - r.y_true, pk[:, -1] - r.x_true)
        axes[1].semilogy(np.full(len(d), r.t), np.clip(d, 1, None), ".", color="#94a3b8", ms=5)
        if r.chosen is not None and not pd.isna(r.chosen):
            axes[1].semilogy([r.t], [max(d[int(r.chosen)], 1)], ".", color="#ef4444", ms=7)
    axes[1].axhline(1000 * FAR_UM, color="#f87171", lw=0.8, ls=":")
    axes[1].set_title("Candidates: distance to truth (nm)")
    g = df[df.found & np.isfinite(df.get("photons_frame", np.nan))] if "photons_frame" in df else df.iloc[:0]
    if len(g):
        axes[2].plot(g.t, g.photons_frame / g.photons_frame.max(), color="#94a3b8", lw=1.5, label="true photons")
        axes[2].plot(g.t, g.intensity / max(np.nanmax(g.intensity), 1e-9), ".", color="#14b8a6", ms=5, label="measured intensity")
        axes[2].legend(fontsize=10)
    axes[2].set_title("Brightness (scaled to max)")
    for ax in axes:
        ax.axvline(df.t.iloc[i], color="orange", lw=1.2)
        ax.set_xlabel("frame")
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    st.caption("Left: distance of the found locus from the truth in every frame. Middle: every candidate's distance from the true spot (grey) and "
               "the one chosen (red); a red dot above the dotted 0.5 µm line is a wrong pick, a frame with a grey dot below the line but no red "
               "dot is a spot that was there but not taken. Right: the measured ellipsoid intensity (dots) and the photons the simulation drew "
               "(line), each scaled to its maximum. Orange line = the frame shown above.")


@st.cache_data(show_spinner="Measuring the MS2 signal at the tracked promoter…", max_entries=12)
def _ms2_promoter(sdir: str, cid: int, params: str, calib: str | None, run: str = "", px_um: float = 0.13) -> pd.DataFrame:
    """The MS2 trace measured at the promoter (locus 0, tracked in its own colour) of every copy: t, copy (numbered like the true copies),
    intensity_at_promoter."""
    from simlive.stage6_analysis.ms2 import ms2_trace

    tr = ms2_trace(load_cell(sdir, cid), calib=calib, ms2_position="promoter", px_um=px_um, **json.loads(params))
    mt = ms2_truth(Path(run)) if run else None
    if mt is not None and "copy" in tr:
        tr = relabel_to_true_copies(tr, mt, px_um)
    return tr[["t"] + (["copy"] if "copy" in tr else []) + ["intensity"]].rename(columns={"intensity": "intensity_at_promoter"})


def _ms2_trace_figure(df: pd.DataFrame, i: int, per_loop: float = 100.0):
    """The cell's transcription trace as Chen et al. measure it (intensity at the tracked position in every frame) against the truth, with the
    spot found in the nuclear channel itself and with the MS2 signal read at the promoter tracked in its own colour."""
    fig, ax = plt.subplots(figsize=(17, 3.8))
    on = (df.photons > 0).to_numpy()
    for t0, t1 in _runs(df.t.to_numpy(), on):
        ax.axvspan(t0 - 0.5, t1 + 0.5, color="#a16207", alpha=0.5, lw=0)
    z = lambda v: (v - np.nanmedian(v)) / max(np.nanpercentile(v, 99) - np.nanmedian(v), 1e-9)  # noqa: E731  (scaled: 0 = typical, 1 = high)
    ax.plot(df.t, z(df.intensity_at_promoter), color="#14b8a6", lw=1.6, label="at the tracked promoter (default)")
    ax.plot(df.t, z(df.trace_intensity), color="#a78bfa", lw=1.0, alpha=0.8, label="at the MS2 spot's own track")
    ax.plot(df.t, df.photons_frame / max(df.photons_frame.max(), 1e-9), color="#e5e7eb", lw=1.2, ls="--", label="true MS2 photons")
    ax.axvline(df.t.iloc[i], color="orange", lw=1.2)
    ax.set_xlabel("frame")
    ax.set_ylabel("scaled")
    ax.set_title("MS2 trace of this cell (yellow = promoter truly ON)")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    r = lambda c: np.corrcoef(df[c].fillna(0), df.photons_frame)[0, 1] if df.photons_frame.std() > 0 else np.nan  # noqa: E731
    st.caption(f"Correlation with the true MS2 signal over all frames: measured at the tracked promoter r = {r('intensity_at_promoter'):.2f}; "
               f"at the MS2 spot's own track r = {r('trace_intensity'):.2f}. The MS2 spot sits on the promoter, so the promoter's own colour gives "
               "its position reliably even when the MS2 spot is faint or OFF; the spot's own track can only follow it when it is brighter than "
               f"the nucleus' texture. This run: {per_loop:g} photons per MS2 loop per second (placeholder 100, at which one Pol II gives ~240 "
               "photons on a ~300 photons/px nucleus and neither trace follows the bursts; measured traces need roughly 3-10x that). Traces "
               "are each scaled to their typical (0) and high (1) values.")


def _runs(t: np.ndarray, flag: np.ndarray) -> list:
    """(first, last) frame of every run of True in flag."""
    out, start, prev = [], None, None
    for tt, f in zip(t, flag):
        if f and start is None:
            start = tt
        if not f and start is not None:
            out.append((start, prev))
            start = None
        prev = tt
    if start is not None:
        out.append((start, prev))
    return out


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
                         "spots with no light": s.get("n_false_on_dark"),
                         "on a stray spot": s.get("frac_on_stray"), "on the other copy": s.get("frac_on_other_copy"),
                         "cells with the right number of copies": s.get("copy_number_right"),
                         "copies missed / extra": f"{s.get('missed_copies')} / {s.get('extra_copies')}" if "missed_copies" in s else None,
                         "intensity r (locus 0 / 1)": " / ".join(f"{s[k]:.2f}" for k in ("intensity_r_locus0", "intensity_r_locus1") if k in s) or None})
    if rows:
        st.dataframe(pd.DataFrame(rows).round(3), hide_index=True, width="stretch")
        st.caption("'On a stray spot' / 'on the other copy' = share of detections more than 0.5 µm from their own locus and within 0.5 µm of a "
                   "stray spot of that colour / of the same locus of another copy (copies matched to the true ones per cell).")
    mrows = []
    for k, lab in (("ms2", "2D"), ("ms2_calibrated", "2D, calibrated")):
        m = sc.get(k)
        if m:
            mrows.append({"MS2 spot": lab, "ON frames located (within 0.5 um)": m.get("found_on"),
                          "median position error (nm)": 1000 * m["median_error_um"] if m.get("median_error_um") is not None else None,
                          "trace r": m.get("intensity_r"), "ON called from trace: sensitivity": m.get("call_sensitivity"),
                          "specificity": m.get("call_specificity")})
    if mrows:
        st.dataframe(pd.DataFrame(mrows).round(3), hide_index=True, width="stretch")
        st.caption("MS2 spot: its intensity (ellipse minus ring, Chen et al. 2025) is measured in the nuclear channel in every frame at the "
                   "promoter's position, tracked in the promoter's own colour (the MS2 spot sits on the promoter). 'Located' = ON frames where that "
                   "position is within 0.5 um of the true spot; 'trace r' = correlation of the trace with the true MS2 photons; ON/OFF called "
                   "from the trace alone with a threshold of 3 x its robust noise. Part 2 (spot: MS2) compares this with finding the MS2 spot "
                   "on its own.")
    elif (run / "stage1_chromatin" / "transcription.csv").exists():
        st.info("This run has MS2 spots but no MS2 analysis yet: scripts/compare_localizers.py --run <run> --ms2-only")
        st.caption(f"'Wrong object' = share of detections more than {FAR_UM} µm from the truth (the finder picked a hot pixel, cosmic ray or noise "
                   "peak). 'rms of good ones' leaves those out. Runs made before this tab existed lack some columns.")
    ef = run / "stage6_analysis" / "truth" / "localization_errors.csv"
    if ef.exists():
        def draw():
            e = pd.read_csv(ef)
            e = e[e.visible.astype(bool)] if "visible" in e else e
            nm = 1000 * e.err_um.clip(lower=1e-3)
            fig, ax = plt.subplots(figsize=(17, 3.4))
            ax.hist(nm, bins=np.geomspace(max(float(nm.min()), 1.0), max(float(nm.max()), 2000.0), 70), color="#14b8a6")
            ax.set_xscale("log")
            ax.axvline(1000 * FAR_UM, color="#f87171", lw=1.2)
            far = int((nm > 1000 * FAR_UM).sum())
            extra = ""
            if "on_stray" in e:
                extra = f": {int(e.on_stray.sum())} on a stray spot, {int(e.on_other_copy.sum())} on the other copy"
            ax.text(1000 * FAR_UM, 0.95, f"  wrong object: {far} of {len(e)}{extra}", color="#f87171", transform=ax.get_xaxis_transform(), va="top", fontsize=10)
            ax.set_xlabel("distance of each detection from its true locus (nm, log scale)")
            ax.set_ylabel("detections")
            ax.set_title(f"All 2D detections of the run (true cell IDs): median {float(nm.median()):.0f} nm")
            return fig

        style.cached_figure(st, ("spt_hist", str(ef), ef.stat().st_mtime), draw)


def _part_bleach(run: Path):
    st.markdown("### 3b · Photobleaching correction")
    bf, sj = run / "stage6_analysis" / "bleach.csv", run / "stage6_analysis" / "scores.json"
    if not bf.exists():
        st.info("This run has no bleaching correction yet (runs made from now on get it; for older runs: "
                "`scripts\\run.cmd python scripts\\compare_localizers.py --run <run>`).")
        return
    st.markdown("A pixel records **offset + background + signal × b(t)**: only the fluorophores bleach. So the decay b(t) is measured on the "
                "signal with the background subtracted (nuclei minus the empty area around them; for the loci, their ring-subtracted spot "
                "intensities, median over all cells), and the correction divides only the background-subtracted spot intensity by b(t) "
                "(column `intensity_corr`). Measured on the raw brightness instead (\"naive\"), the constant background dilutes the decay and "
                "bleaching looks weaker than it is.")
    b = pd.read_csv(bf)

    def draw():
        chans = [c for c in ("nucleus", "locus0", "locus1") if f"{c}_est" in b]
        fig, axes = plt.subplots(1, len(chans), figsize=(17, 3.6), squeeze=False)
        for ax, c in zip(axes[0], chans):
            if f"{c}_true" in b:
                ax.plot(b.t, b[f"{c}_true"], color="#e5e7eb", lw=2, label="true (simulation)")
            ax.plot(b.t, b[f"{c}_est"], color="#14b8a6", lw=1.6, ls="--", label="estimated, background removed")
            if f"{c}_naive" in b:
                ax.plot(b.t, b[f"{c}_naive"], color="#f87171", lw=1.4, ls=":", label="naive, background left in")
            ax.set_title({"nucleus": "nuclear channel (also MS2)", "locus0": "locus 0 (promoter)", "locus1": "locus 1 (enhancer)"}[c], fontsize=13)
            ax.set_xlabel("frame")
            ax.legend(fontsize=9)
        axes[0][0].set_ylabel("b(t), 1 = no bleaching")
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("bleach", str(bf), bf.stat().st_mtime), draw)
    sc = (json.loads(sj.read_text()) if sj.exists() else {}).get("bleaching") or {}
    rows = []
    for key, lab in (("locus0_r", "locus 0 intensity"), ("locus1_r", "locus 1 intensity"), ("ms2_r", "MS2 trace")):
        v = sc.get(key) or {}
        if v:
            rows.append({"what": lab, "r with true photons: raw": v.get("r_raw"), "r: one curve": v.get("r_corrected_global"),
                         "r: illumination + per cell": v.get("r_corrected"),
                         "late / early mean: raw": v.get("late_early_raw"), "late / early: one curve": v.get("late_early_corrected_global"),
                         "late / early: illumination + per cell": v.get("late_early_corrected"),
                         "late / early: true photons": v.get("late_early_true")})
    if rows:
        st.dataframe(pd.DataFrame(rows).round(3), hide_index=True, width="stretch")
        st.caption("'late / early' = mean of the last quarter of frames over the first quarter: raw intensities fall with the bleaching, corrected "
                   "ones should match the true photons (which change only through biology: probes, bursts). True photons are counted before bleaching.")


@st.cache_data(show_spinner=False, max_entries=6)
def _true_flats(run: str, stamp: float) -> dict:
    """The simulation's true illumination map of every channel (mean 1; None = even light)."""
    from simlive.stage3_microscopy.render import imager_for_run

    cfg = yaml.safe_load((Path(run) / "stage2_cells" / "params.yaml").read_text())
    im = imager_for_run(Path(run), cfg)
    return {c: im._flat(c, im.E) for c in im.channels}


def _part_flat(run: Path):
    st.markdown("### 3c · Uneven illumination and per-cell bleaching, from the movie itself (no flat-field image)")
    d6 = run / "stage6_analysis" / "truth"
    ff, ns = d6 / "flatfield_est.npz", d6 / "nuclear_signals.csv"
    sc = ((json.loads((run / "stage6_analysis" / "scores.json").read_text()) if (run / "stage6_analysis" / "scores.json").exists() else {})
          .get("bleaching") or {}).get("flatfield") or {}
    if not ff.exists() or not sc:
        st.info("This run has no illumination estimate yet: `scripts\\run.cmd python scripts\\compare_localizers.py --run <run> --ms2-only`.")
        return
    st.markdown(
        "The cells move, so the **same nucleus is seen at many places** in the field. Its background-subtracted brightness is "
        "*its own brightness × illumination at its position × its bleaching*, and its bleaching follows the **light it has received along its "
        "own path**: a cell sitting in the bright centre bleaches faster than one in a dim corner, and a cell that walks in from outside the "
        "field arrives unbleached. Fitting all cells together (per-cell brightness free, illumination = smooth 2D surface, one bleaching rate) "
        "separates the illumination (changes that follow position) from bleaching (decline that follows accumulated light). The locus colours "
        "use the nuclear channel's illumination shape (same light path) and fit only their own bleaching rate: spot brightness fluctuates with "
        "probe binding far too much to map the illumination from the spots (tried: map error 12% vs 4% for doing nothing). A safeguard: the illumination term is kept only if it predicts **held-out cells** at least "
        "2% better than even light (4-fold cross-validation over cells); otherwise even light is assumed and only bleaching is fitted.")
    rows = []
    for c, lab in (("nucleus", "nuclear channel (+ MS2)"), ("locus0", "locus 0"), ("locus1", "locus 1")):
        s = sc.get(c)
        if not s:
            continue
        cv = s.get("cross_validation") or {}
        rows.append({"channel": lab, "illumination used": s.get("flat_source") if s.get("flat_source") and not s.get("fitted_flat") else "fitted" if s.get("fitted_flat") else "even (fit rejected)" if s.get("flat_rejected") else "even",
                     "held-out gain": cv.get("gain"),
                     "map error: estimate": s.get("flat_rms_error_est"), "map error: assume even": s.get("flat_rms_error_none"),
                     "map error: flat-field slide": s.get("flat_rms_error_calibration"),
                     "tau estimated (s)": s.get("tau_s"), "tau true (s)": s.get("tau_true_s")})
    st.dataframe(pd.DataFrame(rows).round(3), hide_index=True, width="stretch")
    st.caption("Map error = rms difference from the true illumination map over the whole field (both scaled to mean 1). 'Assume even' = the error of "
               "doing nothing; 'flat-field slide' = the calibration acquisition's map (needs the extra images). Tau = bleaching time constant at the "
               "field's mean illumination.")
    n = sc.get("nucleus") or {}
    if "within_cell_sd_raw" in n:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Nuclear brightness spread within a cell: raw", f"{100 * n['within_cell_sd_raw']:.1f} %")
        c2.metric("… one bleaching curve for all", f"{100 * n.get('within_cell_sd_global', np.nan):.1f} %")
        c3.metric("… illumination + per-cell bleaching", f"{100 * n['within_cell_sd_corrected']:.1f} %")
        if "per_cell_b_rms_error_est" in n:
            c4.metric("Per-cell bleaching error (rms)", f"{n['per_cell_b_rms_error_est']:.3f}",
                      f"{n['per_cell_b_rms_error_est'] - n['per_cell_b_rms_error_global']:+.3f} vs one curve", delta_color="inverse")
        st.caption("A nucleus' NLS-GFP brightness is constant in the simulation, so whatever spread is left over time (median over cells, sd of "
                   "log brightness) is noise, texture, shape change, or an uncorrected illumination / bleaching effect.")
    est = np.load(ff)
    true = _true_flats(str(run), ff.stat().st_mtime)

    def draw_maps():
        chans = [c for c in ("nucleus", "locus0", "locus1") if c in est.files]
        fig, axes = plt.subplots(2, len(chans), figsize=(17, 8.2), squeeze=False)
        for j, c in enumerate(chans):
            t = true.get(c)
            t = np.ones_like(est[c]) if t is None else t / t.mean()
            lo, hi = float(min(t.min(), est[c].min())), float(max(t.max(), est[c].max()))
            for i, (m, lab) in enumerate(((t, "true"), (est[c], "estimated from the movie"))):
                ax = axes[i][j]
                h = ax.imshow(m, cmap="magma", vmin=lo, vmax=hi)
                ax.set_title(f"{c}: {lab}", fontsize=12)
                ax.set_xticks([]); ax.set_yticks([])
                if i == 1 and np.ptp(m) < 1e-9:      # the self-calibration kept even light: say so instead of showing an unexplained flat colour
                    ax.text(0.5, 0.5, "fit rejected:\nnot better than even light\n(4-fold cross-validation)", transform=ax.transAxes, ha="center", va="center",
                            fontsize=13, color="white", bbox={"boxstyle": "round,pad=0.4", "fc": "black", "ec": "none", "alpha": 0.5})
                fig.colorbar(h, ax=ax, fraction=0.046)
        if ns.exists():
            nn = pd.read_csv(ns)
            for i in range(2):
                axes[i][0].plot(nn.x_px, nn.y_px, ".", ms=1.5, color="#22d3ee", alpha=0.5)
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("flatmaps", str(ff), ff.stat().st_mtime), draw_maps)
    st.caption("Cyan dots (nuclear channel) = where the nuclei were seen: the illumination can only be learnt where the cells went.")
    bc = run / "stage3_microscopy" / "bleach_cells.csv"
    if ns.exists() and bc.exists():
        nn = pd.read_csv(ns).merge(pd.read_csv(bc), on=["t", "cell_id"], how="left")
        bl = pd.read_csv(run / "stage6_analysis" / "bleach.csv").set_index("t")

        def draw_cells():
            fig, ax = plt.subplots(figsize=(17, 4))
            cells = nn.groupby("cell_id").t.min().sort_values()
            pick = list(cells.index[:3]) + list(cells.index[-3:]) if len(cells) > 6 else list(cells.index)
            colors = plt.cm.tab10(np.linspace(0, 1, len(pick)))
            for col, cid in zip(colors, pick):
                g = nn[nn.cell_id == cid].sort_values("t")
                ax.plot(g.t, g.bleach_nucleus, color=col, lw=2.2, alpha=0.5)
                ax.plot(g.t, g.b_est, color=col, lw=1.2, ls="--", label=f"cell {cid} (first seen at frame {int(g.t.min())})")
            ax.plot(bl.index, bl.nucleus_est, color="#e5e7eb", lw=1.2, ls=":", label="one curve for the whole movie")
            ax.set_xlabel("frame"); ax.set_ylabel("bleaching factor, nucleus")
            ax.set_title("Per-cell bleaching: thick = true, dashed = estimated (cells seen first and last)")
            ax.legend(fontsize=8, ncol=2)
            fig.tight_layout()
            return fig

        style.cached_figure(st, ("cellbleach", str(ns), ns.stat().st_mtime), draw_cells)
        st.caption("An estimated curve starts at 1 when the cell is first seen: light it got before that (outside the camera's field but inside "
                   "the lit area, or before the movie) cannot be known; with light only on the imaged field it is exact.")


def _part_ms2_sweep():
    st.markdown("### 3d · How bright must the MS2 spot be? (sweep over many movies)")
    f = DATA_ROOT / "sweeps" / "ms2_brightness" / "results.csv"
    if not f.exists():
        st.info("No sweep yet: `scripts\\run.cmd python scripts\\sweep_ms2_brightness.py` (about 3 minutes per movie).")
        return
    d = pd.read_csv(f)
    st.markdown("**R** = photons of one Pol II carrying all its MS2 loops ÷ photons of one pixel of nuclear GFP (same exposure). The placeholder is "
                "R = 0.8. Each point = mean over seeds (bars = range); calibration applied, true cell IDs.")
    g = d.groupby("ratio")

    def draw():
        fig, axes = plt.subplots(1, 3, figsize=(17, 4))
        for ax, (col_p, col_o, title) in zip(axes, (("promoter_intensity_r", "own_intensity_r", "trace r with the true MS2 signal"),
                                                    ("promoter_found_on", "own_found_on", "ON frames located (within 0.5 µm)"),
                                                    ("promoter_call_sensitivity", "own_call_sensitivity", "ON frames called from the trace"))):
            for col, lab, colr in ((col_p, "read at the tracked promoter (default)", "#14b8a6"), (col_o, "MS2 spot's own track", "#a78bfa")):
                if col not in d:
                    continue
                m, lo, hi = g[col].mean(), g[col].min(), g[col].max()
                ax.errorbar(m.index, m, yerr=[m - lo, hi - m], marker="o", color=colr, capsize=3, label=lab)
            ax.axvline(0.8, color="#f87171", ls=":", lw=1)
            ax.set_xscale("log"); ax.set_xlabel("R (one Pol II / one nuclear pixel)"); ax.set_title(title, fontsize=12)
            ax.set_ylim(-0.05, 1.05)
        axes[0].legend(fontsize=9)
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("ms2sweep", str(f), f.stat().st_mtime), draw)
    tab = g.agg(**{"movies": ("seed", "size"), "r at promoter": ("promoter_intensity_r", "mean"), "r own track": ("own_intensity_r", "mean"),
                   "located, own track": ("own_found_on", "mean"), "ON called (sens.)": ("promoter_call_sensitivity", "mean"),
                   "OFF called (spec.)": ("promoter_call_specificity", "mean")}).reset_index()
    st.dataframe(tab.round(3), hide_index=True, width="stretch")
    st.caption("Red dotted line = the placeholder. To place your experiment on the x axis, measure R on a real movie: background-subtracted spot "
               "photons per Pol II (e.g. the spot of a single nascent transcript, or the spot sum / estimated Pol II number) over the mean nuclear "
               "pixel value with the camera offset removed.")


def _part_copies_sweep():
    st.markdown("### 3e · Several copies and stray spots: does the linker keep them apart? (sweep over many movies)")
    f = DATA_ROOT / "sweeps" / "copies_strays" / "results.csv"
    if not f.exists():
        st.info("No sweep yet: `scripts\\run.cmd python scripts\\sweep_copies_strays.py` (about 3 minutes per movie, 25 movies).")
        return
    d = pd.read_csv(f)
    n_seeds = d.groupby(["stray_rate_per_min", "analysis"]).seed.nunique().min()
    st.markdown("Every movie: polymer-library loci, half the cells with two promoter-enhancer copies (copy 1 bursting twice as often), stray "
                "spots in all three colours at the rate on the x axis (each lives ~1 min, so the rate is also about how many of each colour sit in "
                f"a nucleus at once), true cell IDs. The same movie is analysed four ways. Points = mean over {n_seeds} seeds, bars = range.")
    styles = {"one_path": ("one path per colour (Chen et al.)", "#f87171"), "multi": ("multi-copy linker, no pair prior", "#94a3b8"),
              "multi_gate": ("multi-copy + flat pair gate (default)", "#2dd4bf"), "multi_gauss": ("multi-copy + smooth pair prior", "#a78bfa")}
    panels = [("copy_number_right", "cells with the right number of copies", 1, (0, 1.05)),
              ("detection_rate", "loci found (of visible)", 1, (0, 1.05)),
              ("frac_on_stray", "detections on a stray spot (%)", 100, None),
              ("frac_on_other_copy", "detections on the OTHER copy (%)", 100, None),
              ("dist_rms_nm", "promoter-enhancer distance error, rms (nm)", 1, None)]

    def draw():
        fig, axes = plt.subplots(1, len(panels), figsize=(18, 3.9))
        for ax, (col, title, sc_, ylim) in zip(axes, panels):
            for a, (lab, colr) in styles.items():
                g = d[d.analysis == a].groupby("stray_rate_per_min")[col]
                if not len(g):
                    continue
                m, lo, hi = sc_ * g.mean(), sc_ * g.min(), sc_ * g.max()
                ax.errorbar(m.index + {"one_path": -0.06, "multi": -0.02, "multi_gate": 0.02, "multi_gauss": 0.06}[a], m,
                            yerr=[m - lo, hi - m], marker="o", ms=4, color=colr, capsize=2, lw=1.4, label=lab)
            ax.set_title(title, fontsize=11)
            ax.set_xlabel("stray spots per nucleus per minute")
            if ylim:
                ax.set_ylim(*ylim)
        axes[0].legend(fontsize=8, loc="lower left")
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("copysweep", str(f), f.stat().st_mtime), draw)
    tab = d.groupby(["stray_rate_per_min", "analysis"]).agg(**{"movies": ("seed", "size"), "copies right": ("copy_number_right", "mean"),
                                                             "found": ("detection_rate", "mean"), "on stray": ("frac_on_stray", "mean"),
                                                             "on other copy": ("frac_on_other_copy", "mean"), "distance rms (nm)": ("dist_rms_nm", "mean"),
                                                             "distance bias (nm)": ("dist_bias_nm", "mean"),
                                                             "MS2 copies traced": ("ms2_copies_traced", "mean"), "MS2 trace r": ("ms2_r", "mean")}).reset_index()
    st.dataframe(tab.round(3), hide_index=True, width="stretch")
    st.caption("'One path per colour' follows ONE spot per colour, so in a two-copy cell it finds one copy at best and jumps between the "
               "copies where they come close. The distance error counts every frame where both loci of a copy were found, wrong picks included.")


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
                         "rms of good ones (nm)": 1000 * tab.rms_error_near_um,
                         "intensity r locus 0": tab.get("intensity_r_locus0"), "intensity r locus 1": tab.get("intensity_r_locus1"),
                         "what changes": tab.what})
    st.dataframe(show.round(3), hide_index=True, width="stretch")
    for cal in sorted(tab.calibrated.unique()):
        d = tab[(tab.calibrated == cal) & (tab.dim == "2d")]
        if d.empty:
            continue
        style.cached_figure(st, ("spt_bars", str(mc), mc.stat().st_mtime, bool(cal)), lambda d=d, cal=cal: _bars(d, cal))


def _bars(d: pd.DataFrame, cal: bool):
    """Found / wrong object / median error of every 2D variant (one calibration setting)."""
    fig, axes = plt.subplots(1, 3, figsize=(17, 0.45 * len(d) + 1.6), sharey=True)
    y = np.arange(len(d))
    for ax, col, title, scale in ((axes[0], "detection_rate", "found (of visible spots)", 1), (axes[1], "frac_far", "wrong object (share)", 1),
                                  (axes[2], "median_error_um", "median error (nm)", 1000)):
        ax.barh(y, scale * d[col].fillna(0), color=["#14b8a6" if v == "default" else "#94a3b8" for v in d.variant])
        ax.set_title(title, fontsize=13)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(d.variant)
    axes[0].invert_yaxis()
    fig.suptitle(f"2D, {'calibration applied' if cal else 'no calibration'} (green = current default)", fontsize=13, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def render(st_, run_path: Path, R: dict):
    style.intro(st, "How stage 6 finds each locus in each cell: the steps (Chen et al. 2025), every intermediate image for "
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
    _part_bleach(run)
    _part_flat(run)
    _part_ms2_sweep()
    _part_copies_sweep()
    _part_compare(run)
