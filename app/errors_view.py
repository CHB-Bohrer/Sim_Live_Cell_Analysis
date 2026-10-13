"""The Images page's part 2: pick ONE source of image error and see what it is, what it does to the picture, how it is corrected (by the
analysis itself, with calibration images, or from the movie alone) and - measured - what it does to the RESULTS with and without each
correction (stage7_validation/error_effects.py, scripts/error_effects.py).
"""
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

import style
from simlive.io.runs import REPO_ROOT
from simlive.stage3_microscopy import errors as ERR
from simlive.stage7_validation import error_effects as EE

FIX_COLOR = {"none": "#f87171", "nomedian": "#f59e0b", "global": "#f59e0b", "analysis": "#60a5fa", "movie": "#2dd4bf", "images": "#a78bfa",
             "both": "#e879f9", "movie_median": "#fb923c", "default": "#94a3b8"}
FIX_SHORT = {"none": "not corrected", "nomedian": "path, no median", "global": "one curve", "analysis": "the analysis' own step",
             "movie": "from the movie alone", "images": "calibration images", "both": "calibration images + movie", "movie_median": "movie, median map (naive)", "default": "error-free reference"}
METRIC_INFO = {   # column -> (title, unit, best value: 'low' | 'high' | 'one' | 'zero')
    "frac_far": ("wrong object picked", "% of detections", "low"), "median_err_nm": ("median position error", "nm", "low"),
    "p95_err_nm": ("95th-percentile position error", "nm", "low"), "dist_rms_nm": ("promoter-enhancer distance error (rms)", "nm", "low"),
    "dist_bias_nm": ("promoter-enhancer distance bias", "nm", "zero"), "intensity_err_pct": ("spot intensity error", "% (after one scale)", "low"),
    "intensity_r": ("spot intensity vs true photons (r)", "", "high"), "late_early": ("intensity drift over the movie (late / early)", "1 = none", "one"),
    "ms2_r": ("MS2 trace vs true signal (r)", "", "high"), "found": ("loci found", "fraction", "high")}
SHOW = ["frac_far", "median_err_nm", "dist_rms_nm", "intensity_err_pct", "late_early", "ms2_r"]


def _ee(run: Path) -> pd.DataFrame | None:
    f = run / "stage7_validation" / "error_effects.csv"
    return pd.read_csv(f) if f.exists() else None


def _run_button(run: Path, key: str, sources: str | None = None):
    log = run / "stage7_validation" / "error_effects.log"
    running = log.exists() and (pd.Timestamp.now().timestamp() - log.stat().st_mtime) < 180 and "finished" not in log.read_text(errors="ignore")[-200:]
    if running:
        st.info("The experiment is running in the background (log: stage7_validation/error_effects.log). Its rows appear here as they finish; "
                "come back to this page in a minute.")
        st.code("\n".join(log.read_text(errors="ignore").splitlines()[-6:]))
        return
    label = "▶ Measure this source (about 1-2 minutes, in the background)" if sources else "▶ Measure every source (about 30 minutes, in the background)"
    if st.button(label, key=key):
        (run / "stage7_validation").mkdir(exist_ok=True)
        cmd = [sys.executable, str(REPO_ROOT / "scripts" / "error_effects.py"), "--run", str(run)] + (["--sources", sources] if sources else [])
        with open(log, "w", encoding="utf-8") as fh:
            fh.write("started\n")
        # detached so leaving the page does not stop it; its output goes to the log file
        flags = (subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.DETACHED_PROCESS) if os.name == "nt" else 0
        subprocess.Popen(cmd, cwd=REPO_ROOT, stdout=open(log, "a", encoding="utf-8"), stderr=subprocess.STDOUT, creationflags=flags,
                         env={**os.environ, "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"})
        st.success("Started. Results appear here as each source finishes.")


def _how(sid: str):
    h = EE.HOW.get(sid, {})
    fixes = {k for k, _, _ in EE.FIXES.get(sid, [])}
    c = st.columns(3)
    for col, (key, title, done) in zip(c, (("analysis", "By the analysis itself", bool(fixes & {"analysis", "nomedian", "global"}) or sid in ("bleaching",)),
                                            ("images", "With calibration images", "images" in fixes),
                                            ("movie", "From the movie alone (no calibration images)", "movie" in fixes))):
        txt = h.get(key, "") or "Nothing to correct."
        tag = "✅ applied and measured below" if done else ("➖ nothing to do" if not h.get(key) or txt.startswith("(") else "ℹ described, not applied")
        col.markdown(f"**{title}**  \n{txt}")
        col.caption(tag)


def _images(run: Path, im, R: dict, ctx, sid: str):
    import images_tab as IT
    from simlive.stage3_microscopy.calibration.model import Calibration
    from simlive.stage6_analysis.localize import row_offsets_removed

    t, _, cell, lk, _ = ctx
    chans = im.channels
    default_ch = "nucleus" if sid in ("flat_field", "bleaching", "flicker", "autofluorescence", "background") else next((c for c in chans if c.startswith("locus")), chans[0])
    c = st.columns([2, 2, 3])
    ch = c[0].selectbox("Channel", chans, index=chans.index(default_ch), format_func=IT.CH_SHORT.get, key=f"ee_ch_{sid}")
    mode = c[1].radio("Region", ["cell", "locus"], horizontal=True, key="ee_mode", format_func={"cell": "One cell", "locus": "Zoom on a locus"}.get)
    cells_now = R["cells"][R["cells"].t == t]
    if cells_now.empty:
        st.info("No cells in this frame.")
        return
    if cell not in cells_now.cell_id.tolist():
        cell = int(cells_now.cell_id.iloc[0])
    cell_row = cells_now[cells_now.cell_id == cell].iloc[0]
    lt = IT._csv(run, "stage1_chromatin", "loci_truth.csv")          # true positions without drift (the reference has no drift)
    q = lt[(lt.cell_id == cell) & (lt.t == t)] if lt is not None else pd.DataFrame()
    kk = int(ch[5:]) if ch.startswith("locus") else lk
    qq = q[q.locus_id == kk]
    spot = (qq.y_um.iloc[0] / im.px_um, qq.x_um.iloc[0] / im.px_um) if len(qq) else None
    if mode == "locus" and spot is None:
        mode = "cell"
    bbox = IT._crop_box(mode, im.shape, im.px_um, cell_row, spot)
    y0, y1, x0, x1 = bbox
    c[2].caption(f"Frame {t}, cell {cell} (change them in part 1). Every panel has its own background level removed for display and the same "
                 "grey scale; the difference panels show what changed relative to the error-free picture.")
    stamp = (run / "stage7_validation" / "metrics.json").stat().st_mtime
    ref = tuple(sorted(EE.REFERENCE))
    withs = tuple(sorted(set(EE.REFERENCE) | {sid}))
    a = IT._render(str(run), stamp, t, ch, ref, bbox, -1).astype(float)
    b = IT._render(str(run), stamp, t, ch, withs, bbox, -1).astype(float)
    panels = [("error-free reference", a), (f"+ {ERR.SOURCES[sid]['title'].lower()}", b)]
    for key, label, spec in EE.FIXES.get(sid, []):
        corr = None
        if spec.get("calib") == "images" and (run / "stage3_microscopy" / "calibration" / "calibration.json").exists():
            corr = Calibration.load(run / "stage3_microscopy" / "calibration").restricted(spec["parts"]).correct(b, ch, y0, x0)
        elif spec.get("calib") == "movie":
            d = run / "stage7_validation" / "error_effects_selfcal" / (f"{sid}__{'_'.join(spec['parts'])}"
                                                                         + (f"_{spec['background_stat']}" if spec.get("background_stat") else ""))
            if (d / "calibration.json").exists():
                cal = Calibration.load(d)
                corr = cal.correct(b, ch, y0, x0)
                f = ((cal.values.get("crosstalk") or {}).get(ch) or {}).get("nucleus")
                if f:                              # unmixing: minus fraction x the nuclear picture of the same frame
                    nb = cal.correct(IT._render(str(run), stamp, t, "nucleus", withs, bbox, -1).astype(float), "nucleus", y0, x0)
                    corr = corr - f * (nb - np.median(nb))
        elif spec.get("row_correct"):
            lab = IT._labels_at(run, im, t, -1)[y0:y1, x0:x1]
            corr = row_offsets_removed(b[None], (lab == cell)[None])[0]
        if corr is not None:
            panels.append((f"corrected: {FIX_SHORT.get(key, key)}", corr.astype(float)))
    shown = [(lab, p - np.median(p)) for lab, p in panels]
    lo, hi = np.percentile(shown[1][1], [1, 99.8])
    n = len(shown)
    fig, axes = plt.subplots(2, n, figsize=(4.2 * n, 8.2), squeeze=False)
    ref0 = shown[0][1]
    noise = 1.4826 * float(np.median(np.abs(ref0 - np.median(ref0))))
    dl = max(max(float(np.percentile(np.abs(p - ref0), 99.9)) for _, p in shown[1:]), 0.5 * noise, 1e-3)   # sparse changes (hot pixels) stay visible
    for j, (lab, p) in enumerate(shown):
        axes[0][j].imshow(p, cmap="gray", vmin=lo, vmax=hi)
        axes[0][j].set_title(lab, fontsize=12)
        if j:
            axes[1][j].imshow(p - ref0, cmap="RdBu_r", vmin=-dl, vmax=dl)
            axes[1][j].set_title(f"minus the reference (rms {np.sqrt(np.mean((p - ref0) ** 2)):.1f} counts)", fontsize=11)
        if len(qq) and ch.startswith("locus"):
            for r in ((0, 1) if j else (0,)):          # the bottom-left panel is the explanation
                axes[r][j].plot(qq.x_um.iloc[0] / im.px_um - x0, qq.y_um.iloc[0] / im.px_um - y0, "o", mfc="none", mec="#22d3ee", ms=14, mew=1.6)
    for ax in axes.ravel():
        ax.set_xticks([]); ax.set_yticks([])
    axes[1][0].axis("off")
    axes[1][0].text(0.02, 0.5, "Bottom row: red = brighter than\nthe error-free picture,\nblue = darker.\nCyan circle = the locus' true\nposition.",
                    transform=axes[1][0].transAxes, fontsize=12, va="center")
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)


def _results(ee: pd.DataFrame, sid: str):
    rows = ee[(ee.variant == f"+{sid}") | (ee.variant == "reference")].copy()
    if not (rows.variant == f"+{sid}").any():
        return False
    rows["label"] = [FIX_SHORT.get(f, f) if v != "reference" else "error-free reference" for v, f in zip(rows.variant, rows.fix)]
    rows = pd.concat([rows[rows.variant == "reference"], rows[rows.variant != "reference"]])
    key = EE.KEY_METRIC.get(sid, "median_err_nm")
    cols = [key] + [m for m in SHOW if m != key and m in rows and rows[m].notna().any()]
    cols = cols[:6]

    def draw():
        fig, axes = plt.subplots(1, len(cols), figsize=(17, 0.55 * len(rows) + 1.9), sharey=True, squeeze=False)
        y = np.arange(len(rows))
        for ax, m in zip(axes[0], cols):
            v = rows[m].to_numpy(float) * (100 if m == "frac_far" else 1)
            ax.barh(y, np.nan_to_num(v), color=[FIX_COLOR.get(f, "#94a3b8") for f in rows.fix], edgecolor=["#fde047" if m == key else "none"] * len(v),
                    linewidth=2 if m == key else 0)
            for yi, vi in zip(y, v):
                if np.isfinite(vi):
                    ax.text(vi, yi, f" {vi:.3g}", va="center", fontsize=10)
            t, u, _ = METRIC_INFO[m]
            import textwrap

            ax.set_title(textwrap.fill(("★ " if m == key else "") + t, 24) + (f"\n({u})" if u else ""), fontsize=11)
            if METRIC_INFO[m][2] == "one":
                ax.axvline(1, color="#94a3b8", lw=1)
        axes[0][0].set_yticks(y)
        axes[0][0].set_yticklabels(rows.label)
        axes[0][0].invert_yaxis()
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("ee_bars", sid, tuple(rows.fix), tuple(np.round(rows[cols].to_numpy(float).ravel(), 4))), draw)
    ref = rows[rows.variant == "reference"]
    nc = rows[(rows.variant != "reference") & (rows.fix.isin(["none", "nomedian"]))].head(1)
    best = rows[rows.variant != "reference"].iloc[1:] if (rows.variant != "reference").sum() > 1 else None
    t, u, goal = METRIC_INFO[key]
    if len(ref) and len(nc):
        fmt = lambda v: "n/a" if not np.isfinite(v) else (f"{100 * v:.1f}%" if key == "frac_far" else f"{v:.3g}")  # noqa: E731
        parts = [f"**{t}** — error-free: {fmt(float(ref[key].iloc[0]))}; with this error and no correction: {fmt(float(nc[key].iloc[0]))}"]
        if best is not None:
            parts += [f"{FIX_SHORT.get(r.fix, r.fix)}: {fmt(float(getattr(r, key)))}" for r in best.itertuples()]
        st.markdown("; ".join(parts) + ".")
    st.caption("★ = what this source mainly damages (yellow edge). Bars: red = not corrected, blue = the analysis' own countermeasure, green = "
               "corrected from the movie alone, violet = with calibration images, grey = the error-free reference (same cells, same noise). "
               "One movie: differences of a few % are within the noise of a single run.")
    if "selfcal" in rows and rows.selfcal.notna().any():
        import json

        reps = [json.loads(x) for x in rows.selfcal.dropna()]
        bits = []
        for rep in reps[1:]:                       # several movie-only analyses (e.g. median vs robust map): list each
            for ch, r in rep.get("channels", {}).items():
                if "background_bias" in r and ch == "nucleus":
                    bits.append(f"{ch} ({rep.get('background_stat')} map): bias {r['background_bias']:+.2f}, worst 1% {r['background_max_bias']:+.2f} counts")
        rep = reps[0]
        for ch, r in rep.get("channels", {}).items():
            if "background_bias" in r:
                bits.append(f"{ch}: background map vs the true object-free background: median bias {r['background_bias']:+.2f} counts, "
                            f"worst 1% {r['background_max_bias']:+.2f} ({rep.get('background_stat', 'median')} map)")
            if "hot_true" in r:
                bits.append(f"{ch}: found {r['hot_found_true']} of {r['hot_true']} hot pixels, {r['hot_false']} false")
            if "flat_rms_error" in r:
                bits.append(f"{ch}: illumination map error {100 * r['flat_rms_error']:.1f}% (even light: {100 * r['flat_rms_error_even']:.1f}%)")
        if "flat" in rep:
            bits.append(f"illumination term {'used' if rep['flat'].get('used') else 'rejected by the held-out check (even light kept)'}")
        for ch, ct in (rep.get("crosstalk") or {}).items():
            bits.append(f"{ch}: bleed-through of the nuclear colour {100 * ct['fraction']:.2f}% ± {100 * ct['se']:.2f}% "
                        f"(true {100 * ct.get('true', np.nan):.2f}%){'' if ct.get('applied') else ', not applied'}")
        if "shift" in rep:
            s = rep["shift"]
            bits.append(f"colour shift from the movie ({s['n_cells']} cells, {s['n_frames']} frames): ({s['shift_px'][0]:+.2f}, {s['shift_px'][1]:+.2f}) px "
                        f"± ({s['se_px'][0]:.2f}, {s['se_px'][1]:.2f}); true ({s.get('true_px', [np.nan, np.nan])[0]:+.2f}, "
                        f"{s.get('true_px', [np.nan, np.nan])[1]:+.2f}) px; {'applied' if s['significant'] else 'NOT applied: not larger than 2 standard errors (needs more cells)'}")
        if bits:
            st.markdown("**What the movie-only estimate found** (checked against the truth): " + " · ".join(bits))
    return True


def render_part(st_, run: Path, im, R: dict, ctx):
    st.markdown("### 2 · One error source: what it is, what it does to the results, and how it is corrected")
    st.caption("Pick a source. The page explains it, shows the same frame without it, with it, and corrected, and - measured on this run - how "
               "much it damages the results and how much each correction gives back. 'Error-free reference' = a perfect microscope that still "
               "has the unavoidable blur and photon / read noise.")
    if im.is3d:
        st.info("This part works on 2D runs (the measurements re-render the whole movie per source).")
        return
    ids = [s for s in ERR.ORDER if s not in EE.ONLY_3D]
    on = set(im.E.enabled_ids())
    sid = st.selectbox("Error source", ids, index=ids.index("hot_pixels"), key="ee_src",
                       format_func=lambda s: f"{ERR.SOURCES[s]['title']}  ·  {'on' if s in on else 'off'} in this run")
    spec = ERR.SOURCES[sid]
    c = st.columns(3)
    c[0].markdown(f"**What it is**  \n{spec['what']}")
    c[1].markdown(f"**What it does to the picture**  \n{spec['effect']}")
    import images_tab as IT

    c[2].markdown(f"**This run**  \n{'ON' if sid in on else 'off'}" + (f" · {IT._source_params_text(im.E, sid)}" if spec["params"] else "")
                  + ("  \n_Fundamental: part of every image, including the error-free reference._" if sid in EE.REFERENCE else ""))
    st.markdown("**How it is corrected**")
    _how(sid)
    if sid in EE.REFERENCE:
        st.info("This source is part of the error-free reference itself (it is the physics of light and of reading a camera): there is nothing to "
                "correct, it sets the best possible precision.")
    else:
        if st.toggle("Show the images (renders this frame 2-5 times, a few seconds)", key="ee_img"):
            _images(run, im, R, ctx, sid)
    st.markdown("**What it does to the results** (whole movie, every cell, scored against the truth)")
    ee = _ee(run)
    if sid in EE.REFERENCE:
        pass
    elif ee is None or not _results(ee, sid):
        st.info("Not measured on this run yet." + ("" if sid in on else " (This source is off in this run; it is measured only if it was on.)"))
        if sid in on:
            _run_button(run, f"ee_go_{sid}", f"reference,{sid}")
    st.markdown("**All sources together** (the realistic case)")
    if ee is not None and (ee.variant == "all").any():
        _all(ee)
    else:
        st.caption("Not measured yet: the button below re-renders and analyses the movie for every source (about 30 minutes).")
    _run_button(run, "ee_go_all")


def _all(ee: pd.DataFrame):
    rows = pd.concat([ee[ee.variant == "reference"].assign(fix="default"), ee[ee.variant == "all"]])
    rows["what"] = ["error-free reference"] + [r for r in ee[ee.variant == "all"].fix_label]
    show = rows[["what"] + [m for m in ["found", "frac_far", "median_err_nm", "p95_err_nm", "dist_rms_nm", "intensity_r", "intensity_err_pct", "late_early", "ms2_r"]
                            if m in rows]].copy()
    show.columns = ["analysis"] + [METRIC_INFO[m][0] for m in show.columns[1:]]
    st.dataframe(show.round(3), hide_index=True, width="stretch", column_config={"analysis": st.column_config.TextColumn(width="large")})
    st.caption("Every source this run has, at once. Rows: nothing corrected; the analysis' own countermeasures (median filter, shell background, "
               "Viterbi path); plus everything estimated from the movie alone; or plus the calibration images instead.")
