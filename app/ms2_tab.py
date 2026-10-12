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
from simlive.stage6_analysis import deconvolve as DC

RUNS = DATA_ROOT / "runs"
TRUE_C, FIT_C, MEAS_C, ON_C = "#e5e7eb", "#2dd4bf", "#60a5fa", "#a16207"


def _runs_with_ms2() -> list[Path]:
    rs = [p for p in RUNS.glob("*") if (p / "stage6_analysis" / "truth" / "ms2_traces.csv").exists()]
    return sorted(rs, key=lambda p: (p / "stage6_analysis" / "truth" / "ms2_traces.csv").stat().st_mtime, reverse=True)


@st.cache_data(show_spinner=False, max_entries=6)
def _load(path: str, stamp: float):
    run = Path(path)
    d6 = run / "stage6_analysis" / "truth"
    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    dec = pd.read_csv(d6 / "ms2_deconvolution.csv")
    info = json.loads((d6 / "ms2_deconvolution_info.json").read_text())
    z = np.load(d6 / "ms2_deconvolution_samples.npz")
    samples = {int(k): z[k] for k in z.files}
    s1 = run / "stage1_chromatin"
    ev = pd.read_csv(s1 / "transcription_events.csv") if (s1 / "transcription_events.csv").exists() else None
    tx = pd.read_csv(s1 / "transcription.csv") if (s1 / "transcription.csv").exists() else None
    return cfg, dec, info, samples, ev, tx


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
    cfg, dec, info, samples, ev, tx = _load(str(run), f.stat().st_mtime)
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
    _step4(dec, samples, ev, tx, dt, n_loops, kern)
    _step5(dec, samples, ev, tx, dt, info)
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
    traces = [g.sort_values("t").trace_cu.to_numpy() for _, g in dec.groupby("cell_id")]
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


def _step4(dec, samples, ev, tx, dt, n_loops, kern):
    _h(4, "One cell: the trace, the Pol II starts found in it, and its bursts",
       "The deconvolution draws 1,000 possible histories of Pol II starts (posterior samples), each one consistent with the trace and its noise. "
       "Where all of them put a start, the start is certain; where they disagree, the bars are lower and wider. Top: the measured trace, the fit "
       "(average of the histories) and the true signal. Middle: Pol II starts per frame, true ones as white ticks. Bottom: the probability that "
       "the gene is in a burst, from Chen et al.'s rule, with the true promoter-ON periods shaded.")
    cells = sorted(dec.cell_id.unique())
    T_all = int(dec.t.max()) + 1
    n_true = ev[(ev.init_s > -dt) & (ev.init_s <= T_all * dt)].groupby("cell_id").size() if ev is not None else pd.Series(dtype=int)   # during the movie
    best = int(n_true.idxmax()) if len(n_true) else int(cells[0])
    c = st.columns([2, 2, 2, 3])
    cid = c[0].selectbox("Cell", cells, index=cells.index(best) if best in cells else 0, key="dc_cell",
                         format_func=lambda i: f"cell {i} ({int(n_true.get(i, 0))} true starts)" if len(n_true) else f"cell {i}")
    w = c[1].slider("Burst window w (frames)", 1.0, 30.0, 5.0, 1.0, key="dc_w", help="Chen et al.: 5 frames (~50 s). The rate is smoothed over this window.")
    gb = c[2].slider("Starts needed in the window (g_b)", 0.5, 6.0, 2.0, 0.5, key="dc_gb",
                     help="Chen et al.: 2 initiations within the window = ON. Lower it for genes that start Pol II rarely.")
    c[3].caption("Chen et al.'s rule needs bursts of frequent starts (their fly genes start Pol II every few seconds while ON). For a gene that "
                 "starts one Pol II every few minutes, widen the window and/or lower g_b; step 5 shows the agreement with the truth for the setting chosen.")
    g = dec[dec.cell_id == cid].sort_values("t")
    s = samples[int(cid)]
    p_on = DC.burst_calls(s, dt, w, gb)
    T = len(g)
    t_min = g.t.to_numpy() * dt / 60
    fig, ax = plt.subplots(3, 1, figsize=(17, 8.6), sharex=True, gridspec_kw={"height_ratios": [3, 2, 1.3]})
    if tx is not None:
        q = tx[tx.cell_id == cid].sort_values("t")
        ax[0].plot(q.t * dt / 60, q.ms2_loops / n_loops, color=TRUE_C, lw=1.4, ls="--", label="true signal (loops / loops per Pol II)")
        on = q.promoter_on.to_numpy().astype(bool)
        for a in ax:
            for t0, t1 in _runs(q.t.to_numpy(), on):
                a.axvspan(t0 * dt / 60, (t1 + 1) * dt / 60, color=ON_C, alpha=0.35, lw=0)
    ax[0].plot(t_min, g.trace_cu, ".", color=MEAS_C, ms=5, label="measured trace")
    ax[0].plot(t_min, g.fit_cu, color=FIT_C, lw=2, label="fit (posterior mean)")
    ax[0].set_ylabel("Pol II units"); ax[0].legend(loc="upper left", fontsize=10, ncol=3)
    m, sd = s.mean(0), s.std(0)
    ax[1].bar(t_min, m, width=dt / 60 * 0.9, color=FIT_C, yerr=sd, ecolor="#94a3b8", label="Pol II starts per frame (posterior mean ± sd)")
    if ev is not None:
        e = ev[ev.cell_id == cid].init_s.to_numpy() / 60
        e = e[(e >= -dt / 60) & (e <= t_min[-1])]
        ax[1].plot(e, np.full(len(e), max(1.0, float(m.max())) * 1.1), "|", color=TRUE_C, ms=16, mew=2, label="true starts")
    ax[1].set_ylabel("starts"); ax[1].legend(loc="upper left", fontsize=10)
    ax[2].fill_between(t_min, p_on, step="mid", color=FIT_C, alpha=0.8)
    ax[2].set_ylim(0, 1.05); ax[2].set_ylabel("P(burst)"); ax[2].set_xlabel("time (min)")
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    true_n = int(n_true.get(cid, 0)) if len(n_true) else None
    st.caption(f"Cell {cid}: {s.sum(1).mean():.1f} ± {s.sum(1).std():.1f} Pol II starts found" + (f" (true: {true_n})" if true_n is not None else "")
               + f"; in a burst {100 * (p_on >= 0.5).mean():.0f}% of frames" + (f" (true promoter ON: {100 * tx[tx.cell_id == cid].promoter_on.mean():.0f}%)" if tx is not None else "")
               + ". Yellow = the promoter is truly ON (it can be ON for a while before the first Pol II starts).")


def _step5(dec, samples, ev, tx, dt, info):
    _h(5, "The whole run against the truth", "Every cell together. 'Found' = true Pol II starts with an inferred start within 3 frames; 'false' = inferred "
                                             "starts with no true one nearby, per true start; the burst rows use the window and g_b chosen in step 4.")
    sc = info.get("scores") or {}
    if ev is None or tx is None:
        st.info("This run has no saved truth for transcription.")
        return
    w, gb = float(st.session_state.get("dc_w", 5.0)), float(st.session_state.get("dc_gb", 2.0))
    calls, truth = [], []
    for cid, s in samples.items():
        q = tx[tx.cell_id == cid].sort_values("t")
        calls.append(DC.burst_calls(s, dt, w, gb)[: len(q)] >= 0.5)
        truth.append(q.promoter_on.to_numpy().astype(bool)[: len(calls[-1])])
    calls, truth = np.concatenate(calls), np.concatenate(truth)
    per = dec.groupby("cell_id").n_init_mean.sum().rename("found (posterior mean)").to_frame()
    per["true"] = ev[ev.init_s >= -dt].groupby("cell_id").size().reindex(per.index).fillna(0).astype(int)
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


def _notes():
    with st.expander("Method details and how this differs from the paper"):
        st.markdown(
            "- **Model** (Chen et al. Eqs. 8-11): trace = kernel ∗ starts + noise with variance sigma_b² + beta × signal; starts live on a 1 s grid "
            "(the shortest time between two Pol II, footprint / speed), 0 or 1 per step.\n"
            "- **Prior** (Eq. 15): each 1 s step has a start with probability R(t) × 1 s, R from the average activity of all cells "
            "(activity at t + tau/2 divided by tau, tau = the kernel's area).\n"
            "- **Sampling**: the paper uses an adaptive block-Bernoulli Metropolis-Hastings sampler with a floor of 0.05 on the proposal probabilities, "
            "right for genes that start Pol II every few seconds. It cannot propose the sparse starts of a gene that fires every few minutes "
            "(a 6-minute block would propose ~18 starts). The same posterior is sampled here with single-step Gibbs updates (add / remove a start "
            "at each step) and moves that slide a start in time. 500 burn-in sweeps, then 1,000 samples (every 3rd sweep), as in the paper.\n"
            "- **Rate and bursts** (Eq. 18, section V.4): r(t) = starts in the frame / frame interval; bursts = r smoothed with exp(-(2t/w)^4) "
            "≥ g_b / w, on every sample, so P(burst) carries the uncertainty.\n"
            "- **Not modelled** (as in the paper): pausing, variable speed, random release time (the simulation has a random wait), two sister "
            "chromatids (the paper's factor 2 in the grid step).\n"
            "- Code: `src/simlive/stage6_analysis/deconvolve.py`; command line: `scripts\\run.cmd python scripts\\deconvolve_ms2.py --run <run>`.")


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
