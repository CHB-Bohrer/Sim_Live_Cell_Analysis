"""The dashboard's 🧮 MS2 → Pol II page: from each cell's MS2 trace to the moments Pol II started on the gene and to its bursts, with the
Bayesian deconvolution of Chen, Levo, Zoller & Gregor 2025 (src/simlive/stage6_analysis/deconvolve.py), step by step and against the truth.

Works on any run with MS2 traces (including sweep runs made without tracking), chosen at the top of the page.
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

import style
from simlive.io.runs import DATA_ROOT, REPO_ROOT
from simlive.stage1_chromatin.copies import trace_label
from simlive.stage6_analysis import deconvolve as DC

RUNS = DATA_ROOT / "runs"
TRUE_C, FIT_C, MEAS_C, ON_C = "#e5e7eb", "#2dd4bf", "#60a5fa", "#a16207"


def _runs_with_ms2() -> list[Path]:
    rs = [p for p in RUNS.glob("*") if (p / "stage6_analysis" / "truth" / "ms2_traces.csv").exists()]
    return sorted(rs, key=lambda p: (p / "stage6_analysis" / "truth" / "ms2_traces.csv").stat().st_mtime, reverse=True)


@st.cache_data(show_spinner=False, max_entries=8)
def _load(path: str, stamp: float, tag: str = ""):
    """tag "" = the default sampler (Gibbs + shift), "_chen_block" = the paper's sampler."""
    run = Path(path)
    d6 = run / "stage6_analysis" / "truth"
    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    dec = pd.read_csv(d6 / f"ms2_deconvolution{tag}.csv")
    info = json.loads((d6 / f"ms2_deconvolution{tag}_info.json").read_text())
    z = np.load(d6 / f"ms2_deconvolution{tag}_samples.npz")
    samples = {int(k): z[k] for k in z.files if not k.startswith("fine_")}
    fine = {int(k[5:]): z[k] for k in z.files if k.startswith("fine_")}      # P(a start) at every 1 s step (newer runs)
    s1 = run / "stage1_chromatin"
    ev = pd.read_csv(s1 / "transcription_events.csv") if (s1 / "transcription_events.csv").exists() else None
    tx = pd.read_csv(s1 / "transcription.csv") if (s1 / "transcription.csv").exists() else None
    # one trace per cell and promoter-enhancer copy: `key` = trace_id (= cell_id for copy 0, so runs made before copies existed are unchanged)
    dec, ev, tx = (None if d is None else d.assign(key=d["trace_id"] if "trace_id" in d else d["cell_id"]) for d in (dec, ev, tx))
    return cfg, dec, info, samples, fine, ev, tx


@st.cache_data(show_spinner="Trying every burst rule against the truth…", max_entries=8)
def _grid(path: str, stamp: float, tag: str):
    cfg, dec, info, samples, fine, ev, tx = _load(path, stamp, tag)
    return DC.burst_rule_grid(samples, tx, float(cfg["acquisition"]["frame_interval_s"]))


def _deconvolve_button(run: Path):
    if st.button("▶ Run the deconvolution on this run (about 1 minute)", type="primary", key="dc_go"):
        cmd = [sys.executable, str(REPO_ROOT / "scripts" / "deconvolve_ms2.py"), "--run", str(run)]
        with st.status("Deconvolving every cell's MS2 trace…", expanded=True) as status:
            proc = subprocess.Popen(cmd, cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
                                    env={**os.environ, "PYTHONUNBUFFERED": "1", "PYTHONIOENCODING": "utf-8"})
            for line in proc.stdout:
                if line.strip().startswith("cell"):
                    st.write(line.strip())
            status.update(label="Done" if proc.wait() == 0 else "Failed", state="complete" if proc.returncode == 0 else "error")
        st.rerun()


def render(st_, run_path: Path, R: dict):
    style.intro(st, "How the MS2 signal of each cell is turned into the moments RNA polymerase II (Pol II) started on the gene and into "
                    "transcriptional bursts, with the Bayesian deconvolution of Chen, Levo, Zoller & Gregor 2025 (Nat. Struct. Mol. Biol.). "
                    "Each step is shown with the simulation's truth next to it.",
                "Pick a run with MS2 spots (top), read steps 1-3 to see what goes in, pick a cell in step 4 to see what comes out, and use "
                "step 5 to see how well it works for the whole run. The burst-calling settings in step 4 can be changed live.")
    runs = _runs_with_ms2()
    if not runs:
        st.info("No run has MS2 traces yet. Make one with 🔬 New single-cell loci run (sidebar) with 'Model distance-dependent transcription' ticked.")
        return
    names = [p.name for p in runs]
    default = Path(run_path).name if Path(run_path).name in names else next((p.name for p in runs if (p / "stage6_analysis" / "truth" / "ms2_deconvolution.csv").exists()), names[0])
    name = st.selectbox("Run (any run with MS2 spots, including sweep runs made without tracking)", names, index=names.index(default), key="dc_run")
    run = RUNS / name
    f = run / "stage6_analysis" / "truth" / "ms2_deconvolution.csv"
    if not f.exists():
        st.info("This run has MS2 traces but no deconvolution yet.")
        _deconvolve_button(run)
        return
    other = run / "stage6_analysis" / "truth" / "ms2_deconvolution_chen_block.csv"
    tag = ""
    if other.exists():
        pick = st.radio("Sampler", ["gibbs_shift", "chen_block"], horizontal=True, key="dc_sampler",
                        format_func={"gibbs_shift": "Gibbs + shift moves (default)", "chen_block": "Chen et al.'s adaptive block sampler"}.get,
                        help="Two ways to draw the same posterior; step 5 compares them. See 'How the sampling works' at the bottom.")
        tag = "" if pick == "gibbs_shift" else "_chen_block"
    ff = run / "stage6_analysis" / "truth" / f"ms2_deconvolution{tag}.csv"
    cfg, dec, info, samples, fine, ev, tx = _load(str(run), ff.stat().st_mtime, tag)
    dt = float(cfg["acquisition"]["frame_interval_s"])
    tc = cfg["loci"]["transcription"]
    kern = DC.kernel_from_config(cfg, info.get("kernel_dt_s", 1.0))
    n_loops = int((tc.get("params") or {}).get("n_loops", 24))
    R_ratio = n_loops * float(tc.get("photons_per_loop_s", 100.0)) / float(cfg["optics"]["photons_per_px_s"])
    if R_ratio < 5:
        st.warning(f"In this run one Pol II is R = {R_ratio:.1f} nuclear pixels bright (MS2 brightness sweep, Spot finding 3d: below ~8 a single "
                   "Pol II is hidden in the noise). The deconvolution then mostly reflects noise; pick a brighter run to see it work.")
    _step1(cfg, info, R_ratio)
    _step2(kern, cfg)
    _step3(dec, info)
    _step4(dec, samples, fine, info, ev, tx, dt, n_loops, kern)
    _step5(dec, samples, ev, tx, dt, info, run)
    _convergence(run)
    _step6(run, ff, tag, dec, samples, tx, dt, info)
    _notes()


# ----------------------------------------------------------------------------- steps
def _h(n: int, title: str, what: str):
    st.markdown(f"### {n} · {title}")
    st.caption(what)


def _step1(cfg, info, R_ratio):
    _h(1, "What goes in: the MS2 trace, in units of one Pol II",
       "The trace is the spot intensity measured at the tracked promoter in every frame (🔎 Spot finding), with the background around the spot "
       "subtracted and photobleaching / uneven illumination divided out. It is then divided by the intensity of ONE Pol II carrying all its MS2 "
       "loops, so 1.0 means 'one fully loaded Pol II on the gene'. Chen et al. get that unit by matching to smFISH counts of single mRNAs; here "
       "it comes from what is known about the microscope (what a single-transcript calibration would measure).")
    opt, tc = cfg["optics"], cfg["loci"]["transcription"]
    n_loops = int((tc.get("params") or {}).get("n_loops", 24))
    items = [("MS2 photons per loop per s", f"{float(tc.get('photons_per_loop_s', 100)):g}", "PLACEHOLDER brightness of one MS2 loop (config)"),
             ("Loops per Pol II", f"{n_loops}", "MS2 stem loops in the cassette"),
             ("One Pol II vs one nuclear pixel (R)", f"{R_ratio:.1f}", "R = loops x photons per loop / nuclear photons per pixel; see the brightness sweep"),
             ("Unit: intensity of one Pol II", f"{info['unit']:.1f} counts/px", "MS2 photons x QE x gain x PSF fraction inside the measuring ellipse / its pixels"),
             ("Trace column used", info.get("column", "?"), "intensity_corr = background-subtracted, bleaching and illumination corrected")]
    if "unit_check" in info:
        items.append(("Unit check against the truth", f"{info['unit_check']:.2f}",
                      "regression of the trace (in units) on the true loaded loops / loops per Pol II: 1.00 = the unit is exact; below 1 = the "
                      "measured spot is dimmer than one Pol II should be (PSF wider, focus, background subtraction)"))
    style.stat_grid(st, items)


def _step2(kern, cfg):
    _h(2, "What one Pol II looks like in the trace (the kernel)",
       "A Pol II that starts at time 0 first copies the MS2 cassette: each loop it makes lights up, so its signal climbs in "
       "steps to 1. It then travels through the rest of the gene at full brightness, waits at the end, and is released: the signal drops to 0. "
       "The trace is the sum of these shapes for every Pol II on the gene (a convolution), plus noise. The deconvolution runs this backwards.")
    p = cfg["loci"]["transcription"].get("params") or {}

    def draw():
        fig, ax = plt.subplots(figsize=(17, 3.2))
        t = np.arange(len(kern.values)) * kern.dt_fine_s
        ax.step(t, kern.values, where="post", color=FIT_C, lw=2.2)
        v = float(p.get("elongation_kb_min", 2.5)) / 60
        tc = float(p.get("cassette_kb", 1.3)) / v
        tg = tc + float(p.get("gene_kb", 5.0)) / v
        ax.axvspan(0, tc, color="#1e3a5f", alpha=0.6, lw=0)
        ax.axvspan(tg, t[-1], color="#3b2a10", alpha=0.6, lw=0)
        ax.text(tc / 2, 1.08, "copies the\nMS2 loops", ha="center", va="bottom", fontsize=11)
        ax.text((tc + tg) / 2, 1.08, "travels through the gene (full brightness)", ha="center", va="bottom", fontsize=11)
        ax.text((tg + t[-1]) / 2, 1.08, "waits, then\nreleased", ha="center", va="bottom", fontsize=11)
        ax.set_ylim(-0.05, 1.45); ax.set_xlabel("seconds after this Pol II started"); ax.set_ylabel("signal (units of one Pol II)")
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("dc_kernel", tuple(np.round(kern.values[::7], 3)), len(kern.values)), draw)
    st.caption(f"Assumptions (as in the paper): constant elongation speed ({float(p.get('elongation_kb_min', 2.5)):g} kb/min), no pausing, and a fixed "
               f"wait at the end ({float(p.get('dwell_s', 116)):g} s, the mean). The simulation's wait is random (exponential), so a Pol II that "
               "stays longer than the mean looks like an extra Pol II to the model: one source of false events.")


def _step3(dec, info):
    _h(3, "How noisy the trace is (noise model)",
       "Before deconvolving, the method needs to know how much of the trace is noise. Camera noise changes from frame to frame, while the "
       "signal changes slowly (a Pol II stays several minutes), so half the squared change between consecutive frames measures the noise. "
       "Brighter spots are noisier (photon noise grows with the light), so the noise is modelled as sigma² = sigma_b² + beta × signal.")
    traces = [g.sort_values("t").trace_cu.to_numpy() for _, g in dec.groupby("key")]
    x, y, w = DC.noise_points(traces)

    def draw():
        fig, ax = plt.subplots(figsize=(9, 3.4))
        ax.scatter(x, y, s=30 + 200 * w / max(w.max(), 1), color=MEAS_C, label="measured (frames binned by level)")
        xx = np.linspace(0, max(float(x.max()) if len(x) else 1, 1), 50)
        ax.plot(xx, info["sigma_b_cu"] ** 2 + info["beta_cu"] * xx, color=FIT_C, label="fit: sigma_b² + beta × signal")
        ax.set_xlabel("signal level (units of one Pol II)"); ax.set_ylabel("noise variance"); ax.legend(fontsize=10)
        fig.tight_layout()
        return fig

    c1, c2 = st.columns([3, 2])
    with c1:
        style.cached_figure(st, ("dc_noise", round(info["sigma_b_cu"], 5), round(info["beta_cu"], 5), len(x)), draw)
    with c2:
        style.stat_grid(st, [("sigma_b (noise without a spot)", f"{info['sigma_b_cu']:.3f}", "in units of one Pol II: below ~0.3 a single Pol II stands out"),
                             ("beta (extra noise per unit of signal)", f"{info['beta_cu']:.3f}", "photon noise of the spot itself")])
        st.caption("Chen et al. estimate it the same way (variance minus the covariance between neighbouring frames, their section IV.4).")


def _step4(dec, samples, fine, info, ev, tx, dt, n_loops, kern):
    _h(4, "One cell: the trace, the Pol II starts found in it, and its bursts",
       "The deconvolution draws 1,000 possible histories of Pol II starts (posterior samples), each one consistent with the trace and its noise. "
       "Where all of them put a start, the start is certain; where they disagree, the probability is spread out. Panels, top to bottom: (1) the "
       "measured trace, the fit (average of the histories) and the true signal; (2) the PROBABILITY OF A POL II START at every second (the "
       "fraction of histories with a start in that second; its area over an interval = the expected number of starts there), true starts as "
       "white ticks; (3) the same summed per frame, with its spread; (4) the probability that the gene is in a burst, from Chen et al.'s rule, "
       "with the true promoter-ON periods shaded.")
    cells = sorted(dec.key.unique())
    T_all = int(dec.t.max()) + 1
    n_true = ev[(ev.init_s > -dt) & (ev.init_s <= T_all * dt)].groupby("key").size() if ev is not None else pd.Series(dtype=int)   # during the movie
    best = int((n_true - n_true.median()).abs().idxmin()) if len(n_true) else int(cells[0])     # a typical cell (median number of starts)
    c = st.columns([2, 2, 2, 3])
    cid = c[0].selectbox("Cell (and copy)", cells, index=cells.index(best) if best in cells else 0, key="dc_cell",
                         format_func=lambda i: f"{trace_label(i)} ({int(n_true.get(i, 0))} true starts)" if len(n_true) else trace_label(i),
                         help="One trace per cell and promoter-enhancer copy (a cell with two copies has two MS2 traces)")
    w = c[1].slider("Burst window w (frames)", 1.0, 30.0, 5.0, 1.0, key="dc_w", help="Chen et al.: 5 frames (~50 s). The rate is smoothed over this window.")
    gb = c[2].slider("Starts needed in the window (g_b)", 0.5, 6.0, 2.0, 0.5, key="dc_gb",
                     help="Chen et al.: 2 initiations within the window = ON. Lower it for genes that start Pol II rarely.")
    c[3].caption("Chen et al.'s rule needs bursts of frequent starts (their fly genes start Pol II every few seconds while ON). For a gene that "
                 "starts one Pol II every few minutes, widen the window and/or lower g_b; step 5 shows the agreement with the truth for the setting chosen.")
    g = dec[dec.key == cid].sort_values("t")
    s = samples[int(cid)]
    p_on = DC.burst_calls(s, dt, w, gb)
    T = len(g)
    t_min = g.t.to_numpy() * dt / 60
    pf = fine.get(int(cid))
    t_end = float(t_min[-1]) if len(t_min) else 1.0
    zoom = st.slider("Show minutes (zoom in to see single Pol II starts)", 0.0, t_end, (0.0, t_end), 0.5, key="dc_zoom")
    fig, ax = plt.subplots(4, 1, figsize=(17, 10.6), sharex=True, gridspec_kw={"height_ratios": [3, 1.8, 1.8, 1.2]})
    if tx is not None:
        q = tx[tx.key == cid].sort_values("t")
        ax[0].plot(q.t * dt / 60, q.ms2_loops / n_loops, color=TRUE_C, lw=1.4, ls="--", label="true signal (loops / loops per Pol II)")
        on = q.promoter_on.to_numpy().astype(bool)
        for a in ax:
            for t0, t1 in _runs(q.t.to_numpy(), on):
                a.axvspan(t0 * dt / 60, (t1 + 1) * dt / 60, color=ON_C, alpha=0.35, lw=0)
    ax[0].plot(t_min, g.trace_cu, ".", color=MEAS_C, ms=5, label="measured trace")
    ax[0].plot(t_min, g.fit_cu, color=FIT_C, lw=2, label="fit (posterior mean)")
    ax[0].set_ylabel("Pol II units"); ax[0].legend(loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=11, ncol=3)    # legends sit above the panels, not on the data
    e = np.zeros(0)
    if ev is not None:
        e = ev[ev.key == cid].init_s.to_numpy() / 60
        e = e[(e >= -dt / 60) & (e <= t_min[-1])]
    if pf is not None:
        tf = (info.get("fine_t0_s", 0.0) + np.arange(len(pf)) * info.get("kernel_dt_s", 1.0)) / 60
        ax[1].fill_between(tf, pf, color="#f472b6", alpha=0.85, lw=0, label="probability of a Pol II start in each second")
        ax[1].plot(e, np.full(len(e), 1.05), "|", color=TRUE_C, ms=14, mew=2, label="true starts")
        ax[1].set_ylim(0, 1.15); ax[1].set_ylabel("P(start) / s")
        ax[1].legend(loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=11, ncol=2)
    else:
        ax[1].text(0.5, 0.5, "(re-run the deconvolution to get the per-second probability)", ha="center", va="center", transform=ax[1].transAxes)
    m, sd = s.mean(0), s.std(0)
    ax[2].bar(t_min, m, width=dt / 60 * 0.9, color=FIT_C, yerr=sd, ecolor="#94a3b8", label="Pol II starts per frame (posterior mean ± sd)")
    ax[2].plot(e, np.full(len(e), max(1.0, float(m.max())) * 1.1), "|", color=TRUE_C, ms=14, mew=2, label="true starts")
    ax[2].set_ylabel("starts / frame"); ax[2].legend(loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=11, ncol=2)
    ax[3].fill_between(t_min, p_on, step="mid", color=FIT_C, alpha=0.8)
    ax[3].set_ylim(0, 1.05); ax[3].set_ylabel("P(burst)"); ax[3].set_xlabel("time (min)")
    ax[3].set_xlim(zoom[0] - 0.05, zoom[1] + 0.05)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.42)                      # room for the legends above the panels
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    true_n = int(n_true.get(cid, 0)) if len(n_true) else None
    st.caption(f"{trace_label(cid).capitalize()}: {s.sum(1).mean():.1f} ± {s.sum(1).std():.1f} Pol II starts found" + (f" (true: {true_n})" if true_n is not None else "")
               + f"; in a burst {100 * (p_on >= 0.5).mean():.0f}% of frames" + (f" (true promoter ON: {100 * tx[tx.key == cid].promoter_on.mean():.0f}%)" if tx is not None else "")
               + ". Yellow = the promoter is truly ON (it can be ON for a while before the first Pol II starts).")


def _step5(dec, samples, ev, tx, dt, info, run):
    _h(5, "The whole run against the truth", "Every cell together. 'Found' = true Pol II starts with an inferred start within 3 frames; 'false' = inferred "
                                             "starts with no true one nearby, per true start; the burst rows use the window and g_b chosen in step 4.")
    sc = info.get("scores") or {}
    if ev is None or tx is None:
        st.info("This run has no saved truth for transcription.")
        return
    w, gb = float(st.session_state.get("dc_w", 5.0)), float(st.session_state.get("dc_gb", 2.0))
    calls, truth = [], []
    for cid, s in samples.items():
        q = tx[tx.key == cid].sort_values("t")
        calls.append(DC.burst_calls(s, dt, w, gb)[: len(q)] >= 0.5)
        truth.append(q.promoter_on.to_numpy().astype(bool)[: len(calls[-1])])
    calls, truth = np.concatenate(calls), np.concatenate(truth)
    per = dec.groupby("key").n_init_mean.sum().rename("found (posterior mean)").to_frame()
    per["true"] = ev[ev.init_s >= -dt].groupby("key").size().reindex(per.index).fillna(0).astype(int)
    style.stat_grid(st, [("True Pol II starts", f"{sc.get('n_true', int(per.true.sum()))}", "in all cells, during the movie"),
                         ("Found (posterior mean)", f"{sc.get('n_estimated', per.iloc[:, 0].sum()):.1f}", "sum over cells of the expected number of starts"),
                         ("True starts found", f"{100 * sc.get('found', np.nan):.0f} %", "a true start with an inferred one within 3 frames"),
                         ("False starts", f"{sc.get('false_per_true', np.nan):.2f} per true", "inferred starts with no true one within 3 frames"),
                         ("Rate correlation (1 min)", f"{sc.get('rate_r_1min', np.nan):.2f}", "starts per minute, inferred vs true"),
                         ("Burst calls agree with the promoter", f"{100 * (calls == truth).mean():.0f} %",
                          f"window {w:g} frames, g_b {gb:g}; ON {100 * truth.mean():.0f}% true vs {100 * calls.mean():.0f}% called"),
                         ("ON frames called ON", f"{100 * calls[truth].mean():.0f} %" if truth.any() else "n/a", "sensitivity of the burst calls"),
                         ("OFF frames called OFF", f"{100 * (~calls[~truth]).mean():.0f} %" if (~truth).any() else "n/a", "specificity")])

    def draw():
        fig, ax = plt.subplots(figsize=(7, 4))
        mx = max(per.true.max(), per.iloc[:, 0].max(), 1)
        ax.plot([0, mx], [0, mx], color="#94a3b8", lw=1)
        ax.scatter(per.true, per.iloc[:, 0], color=FIT_C, s=50)
        ax.set_xlabel("true Pol II starts per cell"); ax.set_ylabel("found (posterior mean)")
        fig.tight_layout()
        return fig

    c1, c2 = st.columns([2, 3])
    with c1:
        style.cached_figure(st, ("dc_scatter", tuple(per.true), tuple(np.round(per.iloc[:, 0], 2))), draw)
    with c2:
        st.dataframe(per.reset_index().round(2), hide_index=True, width="stretch")
    rows = []
    for tg, lab in (("", "Default: Gibbs + shift moves, checked per cell, tempering where the check fails"), ("_chen_block", "Chen et al.'s adaptive block sampler"),
                    ("_gibbs_pt", "Gibbs + shift moves with parallel tempering")):
        fj = run / "stage6_analysis" / "truth" / f"ms2_deconvolution{tg}_info.json"
        if fj.exists():
            ij = json.loads(fj.read_text())
            sc2 = ij.get("scores") or {}
            rows.append({"sampler": lab, "starts found (posterior mean)": sc2.get("n_estimated"), "true starts": sc2.get("n_true"),
                         "true found within 30 s": sc2.get("found"), "false per true": sc2.get("false_per_true"),
                         "1-min rate r": sc2.get("rate_r_1min"), "acceptance (swap rate for tempering)": ij.get("acceptance") if tg else None, "cells tempered": ij.get("n_tempered"),
                         "block (s)": (np.median(ij["block_s"]) if ij.get("block_s") else None)})
    if len(rows) >= 2:
        st.markdown("**The samplers side by side** (same model, same data: they should agree)")
        st.dataframe(pd.DataFrame(rows).round(3), hide_index=True, width="stretch")


def _convergence(run):
    """Has the sampler converged? Four chains from different starts per cell (scripts/check_mcmc.py), R-hat and effective sample size."""
    _h("5b", "Has the sampler converged?",
       "Four chains per cell start from different configurations (no starts, a draw from the prior, 10x too many starts, starts where the "
       "trace rises). If the sampler works, they end up telling the same story. R-hat compares the spread between chains with the spread inside "
       "them (1.00 = same; the limit here is 1.05, fixed in advance); the effective sample size (ESS) says how many independent samples the "
       "1,000 correlated ones are worth (limit 400).")
    f = run / "stage6_analysis" / "truth" / "mcmc_check.csv"
    if not f.exists():
        st.info("Not run for this movie yet. In a terminal: " + f"`scripts\\run.cmd python scripts\\check_mcmc.py --run {run.name}`" + " (about 10 minutes).")
        return
    df = pd.concat([pd.read_csv(x) for x in (f, f.with_name("mcmc_check_pt.csv")) if x.exists()], ignore_index=True)
    ij = run / "stage6_analysis" / "truth" / "ms2_deconvolution_info.json"
    ck = (json.loads(ij.read_text()).get("checks") if ij.exists() else None)
    if ck:
        k = pd.DataFrame(ck)
        st.success(f"The results above (default sampler) were checked as they were made: {int((k.sampler_used == 'gibbs_shift').sum())} of {len(k)} cells "
                   f"passed with Gibbs + shift moves ({int(k.plain_extended.sum())} of them after drawing 3x more samples), "
                   f"{int((k.sampler_used == 'gibbs_pt').sum())} were sampled again with parallel tempering.")
        with st.expander("Which sampler gave each cell, and the check that decided it"):
            st.dataframe(k.round(3), hide_index=True, width="stretch")
    names = {"gibbs_shift": "Gibbs + shift moves (default)", "chen_block": "Chen et al. block sampler", "gibbs": "Gibbs flips only (no shifts)",
             "gibbs_pt": "Gibbs + shift moves with parallel tempering"}
    g = df.groupby("sampler").agg(cells=("converged", "size"), converged=("converged", "sum"), worst_rhat_logpost=("rhat_logpost", "max"),
                                  median_ess_logpost=("ess_logpost", "median"), worst_chain_difference_per_min=("max_chain_diff", "max"),
                                  seconds_per_cell=("seconds", "mean")).reset_index()
    g["sampler"] = g["sampler"].map(names).fillna(g["sampler"])
    st.dataframe(g.round(2), hide_index=True, width="stretch")
    bad = df[(df.sampler == "gibbs_shift") & ~df.converged]
    if len(bad):
        st.warning(f"{len(bad)} of {int((df.sampler == 'gibbs_shift').sum())} cells did NOT pass for the default sampler: their chains settle in different "
                   "local arrangements of the starts (log-posterior levels that differ by far more than the spread inside a chain). The total number "
                   "of starts agrees better than their exact timing. Re-run the movie's deconvolution with `--sampler gibbs_pt` (parallel tempering) and check "
                   "it with `check_mcmc.py --samplers gibbs_pt`; see docs/MCMC_CONVERGENCE.md.")
    with st.expander("Every cell and sampler"):
        cols = ["cell", "sampler", "converged", "rhat_logpost", "rhat_nstarts", "rhat_frames_max", "ess_logpost", "ess_nstarts", "ess_frames_min",
                "max_chain_diff", "seconds", "reasons"]
        st.dataframe(df[[c for c in cols if c in df]].round(3), hide_index=True, width="stretch")
    fb = run / "stage6_analysis" / "truth" / "mcmc_burn_in.csv"
    if fb.exists():
        b = pd.read_csv(fb)
        b = b.groupby("burn_in_sweeps").agg(median_rhat=("rhat_logpost", "median"), worst_rhat=("rhat_logpost", "max"),
                                            spread_of_chain_means=("spread_of_chain_means", "median")).reset_index()
        st.markdown("**How much burn-in?** R-hat of the log-posterior over the 100 sweeps after discarding the first n sweeps (all cells; the "
                    "paper and our default use 500):")
        st.dataframe(b.round(2), hide_index=True, width="stretch")


def _step6(run, ff, tag, dec, samples, tx, dt, info):
    _h(6, "Choosing the burst rule: why Chen et al.'s default misses bursts of a slow gene",
       "A burst is a period in which Pol II starts often. The rule calls 'ON' where the number of starts within a window of w frames reaches "
       "g_b, i.e. where the start RATE is at least g_b / (w x frame interval).")
    rate_thr = 60 * 2 / (5 * dt)
    st.markdown(
        f"- **Chen et al.'s default**: w = 5 frames (here {5 * dt:.0f} s), g_b = 2 starts → ON needs at least {rate_thr:.1f} starts per minute. "
        "Their fly genes start Pol II every few seconds while ON (~10-20 per minute), so an ON period always passes and an OFF period "
        "(no starts) never does.\n"
        "- **A gene that starts Pol II more slowly does not pass, even while ON.** Starts are random (Poisson): with rate λ while ON, the number in "
        "a window of w frames is Poisson with mean λ·w·Δt. With the denser test gene (3 per minute) the mean in 50 s is 2.5, and the chance of "
        "≥ 2 is ~0.7 at any moment, less after smoothing: about half the ON time is missed. With the placeholder gene (one start per 10 minutes "
        "while ON) the mean in 50 s is 0.08: practically never 2, so no burst is ever called.\n"
        "- **The fix is to match the window to the gene**: enough starts must fall in the window while ON (several λ·w·Δt), and the window must "
        "still be shorter than a typical ON and OFF period, or bursts blur together. If one burst holds only one or two starts (the placeholder "
        "gene: ~1.6), its ON period cannot be seen in the MS2 signal at all: only the starts themselves can.\n"
        "- The calculator below gives the expected agreement for a gene you describe; the map uses this run's truth.")
    c = st.columns(4)
    lam = c[0].number_input("Starts per minute while ON (λ)", 0.01, 60.0, 3.0 if dt >= 1 else 3.0, 0.5, key="dc_lam")
    w = c[1].number_input("Window w (frames)", 1, 60, 5, 1, key="dc_wcalc")
    gb = c[2].number_input("Starts needed (g_b)", 1, 20, 2, 1, key="dc_gbcalc")
    from scipy.stats import poisson

    mu = lam / 60 * w * dt
    c[3].metric("Chance an ON window reaches g_b", f"{100 * poisson.sf(gb - 1, mu):.0f} %", help=f"Poisson with mean λ·w·Δt = {mu:.2f} starts")
    st.caption(f"With λ = {lam:g}/min and w = {w} frames ({w * dt:.0f} s), an ON window holds {mu:.2f} starts on average. An OFF window holds none, "
               "so false calls come only from false starts of the deconvolution. The smoothing kernel (exp(-(2t/w)^4)) makes the real numbers a "
               "little lower than this simple count.")
    if tx is None:
        return
    g = _grid(str(run), ff.stat().st_mtime, tag)

    def draw():
        fig, axes = plt.subplots(1, 3, figsize=(17, 4.2))
        for ax, m, title in zip(axes, ("accuracy", "sensitivity", "specificity"),
                                ("all frames right", "ON frames called ON", "OFF frames called OFF")):
            piv = g.pivot(index="g_b", columns="w_frames", values=m)
            h = ax.imshow(piv.to_numpy(), origin="lower", cmap="viridis", vmin=0, vmax=1, aspect="auto")
            ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels([f"{v:g}" for v in piv.columns])
            ax.set_yticks(range(len(piv.index))); ax.set_yticklabels([f"{v:g}" for v in piv.index])
            ax.set_xlabel("window w (frames)"); ax.set_ylabel("starts needed g_b"); ax.set_title(title, fontsize=13); ax.grid(False)
            for i, gv in enumerate(piv.index):
                for j, wv in enumerate(piv.columns):
                    ax.text(j, i, f"{100 * piv.iloc[i, j]:.0f}", ha="center", va="center", fontsize=11, fontweight="bold",
                            color="black" if piv.iloc[i, j] > 0.62 else "white")
            if 5 in list(piv.columns) and 2 in list(piv.index):
                ax.add_patch(plt.Rectangle((list(piv.columns).index(5) - 0.5, list(piv.index).index(2) - 0.5), 1, 1, fill=False,
                                           ec="#f472b6", lw=2.5))
            fig.colorbar(h, ax=ax, fraction=0.046)
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("dc_grid", str(ff), ff.stat().st_mtime), draw)
    best = g.assign(bal=(g.sensitivity + g.specificity) / 2).sort_values("bal", ascending=False).iloc[0]
    d = g[(g.w_frames == 5) & (g.g_b == 2)]
    st.caption("Pink box = Chen et al.'s default (w = 5 frames, g_b = 2). Numbers are % of frames. For this run the best balance of ON and OFF "
               f"frames is w = {best.w_frames:g} frames ({best.w_frames * dt:.0f} s), g_b = {best.g_b:g}: ON {100 * best.sensitivity:.0f}% / OFF "
               f"{100 * best.specificity:.0f}%" + (f", against the default's ON {100 * d.sensitivity.iloc[0]:.0f}% / OFF {100 * d.specificity.iloc[0]:.0f}%."
                                                   if len(d) else ".") + " With real data there is no truth: pick w from the gene's start rate (calculator).")


def _notes():
    with st.expander("How the sampling works (and why there are two samplers)"):
        st.markdown(
            "**What is being sampled.** The unknown is a list of 0/1 values, one per second of the movie (plus the few minutes before it): 1 = a "
            "Pol II started in that second. Every list predicts a trace (each start adds one kernel shape). The posterior probability of a list "
            "= how well its predicted trace matches the measured one (given the noise model) × the prior (few starts unless the activity says "
            "otherwise). A sampler produces many lists, each drawn with that probability; averages over them give the expected number of starts, "
            "the probability of a start in each second, and the burst probability, all with their uncertainty.\n\n"
            "**Gibbs + shift moves (default).** Go through every second in turn and ask: given all other starts, how likely is a start here? That "
            "is exactly computable (the trace with and without a start at this second, times the prior odds), so the start is added or removed "
            "with exactly that probability (a Gibbs update). Then every existing start is offered a small move in time (up to 15 s earlier or "
            "later), accepted by the Metropolis rule; this lets a start slide to a better place without first having to disappear. Both steps "
            "leave the posterior unchanged, so the samples come from the same posterior as the paper's. 500 sweeps are discarded (burn-in), then "
            "every 3rd of 3,000 sweeps is kept (1,000 samples), as in the paper.\n\n"
            "**Chen et al.'s adaptive block sampler.** Pick a random block of consecutive seconds, redraw the whole block from proposal "
            "probabilities q (one per second, learnt from the chain as it runs, never below 0.05 or above 0.95), and accept or reject the new "
            "block with the Metropolis-Hastings ratio. The paper tuned the block length (~6 min for their genes) so that 10-70% of proposals are "
            "accepted. Here the block length is tuned the same way, automatically (largest block with ≥ 10% acceptance in a short pilot run): "
            "for these traces that is ~20 s, because a 6-minute block filled with q ≥ 0.05 proposes ~18 starts at once and is almost always "
            "rejected (1.7% acceptance, stuck). With the tuned block both samplers give the same answer (step 5).\n\n"
            "**Convergence.** Step 5b runs four chains per cell from different starting configurations and compares them (R-hat, effective sample "
            "size). The default sampler passes on constructed traces but not on the dense, almost noise-free gene, where chains get trapped in "
            "different arrangements of the starts. **Parallel tempering** fixes this: 16 copies of the problem run at likelihood strengths from "
            "1 down to 1/1000 and swap configurations, so the hot copies move freely and hand new arrangements down to the real one. Its chains "
            "agree (see docs/MCMC_CONVERGENCE.md); it costs about 20-30x the time.\n\n"
            "**Rest of the method.** Model (Chen et al. Eqs. 8-11): trace = kernel ∗ starts + noise with variance sigma_b² + beta × signal; prior "
            "(Eq. 15): each second has a start with probability R(t) × 1 s, R from the average activity of all cells (activity at t + tau/2 "
            "divided by tau, the kernel's area); rate and bursts (Eq. 18, section V.4) computed on every sample. Not modelled, as in the paper: "
            "pausing, variable speed, a random release time (the kernel assumes a fixed wait at the gene end), two sister chromatids.\n\n"
            "Code: `src/simlive/stage6_analysis/deconvolve.py`; `scripts\\run.cmd python scripts\\deconvolve_ms2.py --run <run> [--sampler chen_block]`.")


def _runs(t: np.ndarray, on: np.ndarray):
    out, start = [], None
    for ti, o in zip(t, on):
        if o and start is None:
            start = ti
        if not o and start is not None:
            out.append((start, prev))
            start = None
        prev = ti
    if start is not None:
        out.append((start, t[-1]))
    return out
