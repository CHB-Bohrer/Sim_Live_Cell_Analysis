"""The dashboard's 🔗 E–P curve page: how well is P(burst ON | promoter-enhancer distance) recovered, given measurement errors?

Part 1: simulate a population now, apply the measurement errors you choose, fit the curve (with and without correcting for the errors) and compare
it with the truth. Part 2: the saved scans (scripts/sweep_ep_curve.py -> $SIMLIVE_DATA/sweeps/ep_curve_*/results.csv). Method and numbers: docs/EP_CURVE.md.
Imports numpy / scipy / numba through simlive only, never torch.
"""
from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

import style
from simlive.io.runs import DATA_ROOT
from simlive.stage7_validation import ep_acquisition as aq
from simlive.stage7_validation import ep_curve as ep
from simlive.stage7_validation import ep_history as eh
from simlive.stage7_validation import ep_polymer as pp

SWEEPS = DATA_ROOT / "sweeps"


@st.cache_resource(show_spinner=False, max_entries=6)
def _reference(wd: float, cf: float, fold: float, loop_s: float):
    s = ep.EPSetup(working_distance_nm=wd, contact_fraction=cf, fold=fold, loop_time_s=loop_s)
    return ep.reference_curve(s, n_cells=3000)


@st.cache_resource(show_spinner=False, max_entries=6)
def _population(wd: float, cf: float, fold: float, loop_s: float, n_cells: int, seed: int, coupled: bool):
    s = ep.EPSetup(working_distance_nm=wd, contact_fraction=cf, fold=fold, loop_time_s=loop_s, n_cells=n_cells)
    return ep.simulate_population(s, seed, coupled=coupled, with_trace=False)


def render(st_=st) -> None:
    style.intro(st, "The probability that the promoter is in a burst, as a function of its distance to the enhancer: the true curve (from a very large "
                    "simulated population) against the curve an analysis can recover from a realistic number of cells with measurement errors.",
                "Choose the biology and the measurement errors, press Simulate, and compare the curves. The scans below show the same thing over many "
                "populations and many error sizes.")
    _live_polymer()
    st.divider()
    _live()
    st.divider()
    _scans()


@st.cache_resource(show_spinner=False)
def _pool():
    return pp.load_pool()


@st.cache_resource(show_spinner=False, max_entries=4)
def _poly_reference(lag_s: float, applies_to: str):
    cal = pp.calibrate(_pool(), lag_s=lag_s, applies_to=applies_to)
    pop = pp.simulate_population(ep.EPSetup(n_cells=3000), _pool(), cal, 12345, with_trace=False, keep_blocks=False)
    a = eh.analyse_history(pop, pop.d_true, aq.kernel_aniso(9.0, 9.0, 3, 0.0))
    hist, _ = np.histogram(pop.d_true, bins=np.r_[ep.D_GRID - 5.0, ep.D_GRID[-1] + 5.0])
    a["d_pdf"] = hist / hist.sum()
    return a, cal


def _live_polymer() -> None:
    st.markdown("### 1 · Distances from the polymer simulations, seen through a 2D, z-stack or biplane microscope")
    st.caption("E and P are two monomers 100 kb apart in the 10 validation simulations. Length scale: 55 nm per polymer unit, chosen so that the pair is within 400 nm 10% of the "
               "time (the paper checked gives no distances). Time scale: the single-locus MSD at the chosen lag is matched to 6 D t with D = 0.002 µm²/s. The promoter responds to "
               "the distance over ~45 s but the distance fluctuates within ~10 s, so the curve shown is P(ON) for a distance HELD at d, fitted from the distance history.")
    with st.form("epp_form"):
        c = st.columns(4)
        kind = c[0].selectbox("Microscope", ["biplane", "zstack", "2d"], format_func={"biplane": "3D biplane (2 planes, simultaneous)", "zstack": "3D z-stack (17 planes)", "2d": "2D widefield"}.get)
        n_cells = c[1].number_input("Cells (1 h each)", 100, 2000, 600, 100)
        delay = c[2].number_input("Delay between the two colours (s)", 0.0, 0.5, 0.15, 0.05, help="Sequential excitation (Vutara preset: 0.15 s). 0 = both colours at once.")
        photons = c[3].number_input("Photons (x 6000 per locus per s)", 0.1, 10.0, 1.0, 0.1)
        c = st.columns(4)
        focus = c[0].number_input("Focus scatter, 2D / biplane (um)", 0.0, 2.0, 0.3, 0.1)
        wrong = c[1].number_input("Wrong picks", 0.0, 0.5, 0.0, 0.01)
        lag = c[2].number_input("D = 0.002 um2/s refers to lag (s)", 5.0, 300.0, 60.0, 5.0)
        applies = c[3].radio("D is for", ["locus", "pair"], horizontal=True, format_func={"locus": "one locus", "pair": "the E-P distance"}.get)
        c = st.columns(3)
        coupled = c[0].checkbox("Enhancer acts (untick = null control)", True)
        seed = c[1].number_input("Seed", 0, 10**6, 1, key="epp_seed")
        go = c[2].form_submit_button("Simulate and fit", type="primary")
    if not go and "epp_last" not in st.session_state:
        st.info("Press Simulate and fit (about 40 s).")
        return
    if go:
        with st.spinner("Polymer pool, reference (3000 cells), your population, fits..."):
            ref, cal = _poly_reference(float(lag), applies)
            pop = pp.simulate_population(ep.EPSetup(n_cells=int(n_cells)), _pool(), cal, int(seed), coupled=bool(coupled), with_trace=False)
            acq = aq.Acquisition(kind=kind, colour_delay_s=float(delay), photon_scale=float(photons), focus_sd_um=float(focus), wrong_pick_rate=float(wrong),
                                 plane_time_s=2 * float(delay) if kind == "zstack" and delay > 0 else (0.15 if kind == "zstack" else 0.0))
            m = aq.measure(pop, pop.blocks, acq, seed=int(seed))
            K = aq.kernel_aniso(m.sigma_xy, m.sigma_z, m.dims, m.shift_nm)
            res = eh.analyse_history(pop, m.d, K, n_boot=60, seed=int(seed))
            inst = ep.analyse_measured(pop, m.d, K)
        st.session_state["epp_last"] = (ref, cal, res, inst, m, pop.d_true, bool(coupled))
    ref, cal, res, inst, m, d_true, coupled = st.session_state["epp_last"]
    r, f = ref["fit"], res["fit"]
    lo, hi = res.get("contrast_ci", (np.nan, np.nan))
    c = st.columns(5)
    c[0].metric("Truth: P(ON) at 150 nm / at 1000 nm (held)", f"{r.contrast():.2f}", help="Same fit to the true distance and the true ON state of 3000 cells.")
    c[1].metric("Recovered", f"{res['contrast']:.2f}", f"95% {lo:.2f} to {hi:.2f}", delta_color="off")
    c[2].metric("Interval contains the truth?" if coupled else "Interval excludes 1 (false positive)?",
                ("yes" if lo <= r.contrast() <= hi else "no") if coupled else ("yes" if (lo > 1 or hi < 1) else "no"))
    c[3].metric("Memory of the promoter (s): truth / recovered", f"{r.tau_s:.0f} / {f.tau_s:.0f}")
    c[4].metric("Distance only NOW (instantaneous fit)", f"{inst['contrast']:.2f}", help="P(ON | distance now): nearly flat, because the distance changes within ~10 s while the promoter integrates over ~45 s.")
    g = ep.D_GRID
    bs = res.get("boot", [])
    band = np.array([b.g(g) for b in bs]) if bs else None

    def draw():
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 4.4), gridspec_kw={"width_ratios": [3, 2]})
        ax.plot(g, r.g(g), color="#e5e7eb", lw=2.4, label="truth: distance held at d (3000 cells, true distance and state)")
        if band is not None:
            ax.fill_between(g, np.quantile(band, 0.025, 0), np.quantile(band, 0.975, 0), color="#2dd4bf", alpha=0.25, lw=0)
        ax.plot(g, f.g(g), color="#2dd4bf", lw=2, label="recovered (95% band)")
        ax.set_xlim(0, 1500); ax.set_ylim(0, None); ax.set_xlabel("promoter-enhancer distance held at (nm)"); ax.set_ylabel("P(promoter in a burst)")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        ax2.hist(d_true, bins=np.arange(0, 2500, 50), density=True, color="#94a3b8", alpha=0.6, label="true")
        ax2.hist(m.d[np.isfinite(m.d)], bins=np.arange(0, 2500, 50), density=True, histtype="step", color="#60a5fa", lw=1.8, label="measured")
        ax2.set_xlabel("distance (nm)"); ax2.set_ylabel("share of frames"); ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("epp_live", round(res["contrast"], 6), round(r.contrast(), 6), len(bs)), draw)
    parts = {k.replace("_nm", " (nm per axis)").replace("_s", " (s)").replace("_", " "): round(v, 1) for k, v in m.parts.items()}
    st.caption(f"Calibration: {cal.nm_per_unit:.1f} nm per unit, {cal.block_s:.2f} s per saved block (D {cal.applies_to}, lag {cal.lag_s:.0f} s; MSD exponent {cal.alpha:.2f} below one block, an extrapolation). "
               f"Per-axis error of the measured difference vector: {m.sigma_xy:.0f} nm in the image plane" + (f", {m.sigma_z:.0f} nm in z" if m.dims == 3 else "") + f"; both loci seen in {100 * m.detected:.0f}% of frames.")
    with st.expander("Error budget (median per-axis contributions)"):
        st.dataframe(pd.DataFrame({"term": list(parts), "value": list(parts.values())}), hide_index=True, width="stretch")


def _live() -> None:
    st.markdown("### 1b · The earlier version: an invented slow distance process, noise and 2D / 3D only")
    with st.form("ep_form"):
        c = st.columns(4)
        wd = c[0].number_input("Working distance of the enhancer (nm)", 100.0, 1500.0, 400.0, 50.0, help="Distance at which the enhancer's effect is half of its maximum (hill curve, n = 4). The old sharp contact was 150 nm.")
        cf = c[1].number_input("Contact fraction", 0.02, 0.4, 0.10, 0.01, help="Fraction of time the pair is closer than the working distance.")
        fold = c[2].number_input("Fold on k_on at full contact", 1.0, 50.0, 12.0, 1.0)
        loop_s = c[3].number_input("Mean time per visit inside (s)", 30.0, 3600.0, 240.0, 30.0)
        c = st.columns(4)
        n_cells = c[0].number_input("Cells (1 h each)", 50, 4000, 600, 50)
        sigma = c[1].number_input("Localization precision of the difference, per axis (nm)", 1.0, 500.0, 57.0, 5.0, help="57 nm = 40 nm per locus per axis.")
        dims = c[2].radio("Distance measured in", [2, 3], horizontal=True, format_func=lambda d: "2D image plane" if d == 2 else "3D (z-stacks)")
        shift = c[3].number_input("Colour shift left over (nm)", 0.0, 400.0, 0.0, 10.0)
        c = st.columns(4)
        corrected = c[0].checkbox("Analysis knows the shift", False)
        miss = c[1].number_input("Missed frames", 0.0, 0.9, 0.05, 0.05)
        out = c[2].number_input("Wrong picks", 0.0, 0.5, 0.02, 0.01)
        coupled = c[3].checkbox("Enhancer acts (untick = null control)", True)
        seed = st.number_input("Random seed", 0, 10**6, 1)
        go = st.form_submit_button("Simulate and fit", type="primary")
    if not go and "ep_last" not in st.session_state:
        st.info("Press **Simulate and fit** (about 20 s).")
        return
    key = (wd, cf, fold, loop_s, int(n_cells), float(sigma), int(dims), float(shift), bool(corrected), float(miss), float(out), bool(coupled), int(seed))
    if go:
        with st.spinner("Reference population (3000 cells), your population, fits…"):
            ref = _reference(wd, cf, fold, loop_s)
            pop = _population(wd, cf, fold, loop_s, int(n_cells), int(seed), bool(coupled))
            meas = ep.Measurement(sigma_nm=sigma, dims=int(dims), shift_nm=shift, miss_rate=miss, outlier_rate=out, assumed_shift_nm=shift if corrected else 0.0)
            eiv = ep.analyse(pop, meas, n_boot=80, seed=int(seed))
            nv = ep.analyse(pop, meas, naive=True, n_boot=0, seed=int(seed))
            d_meas = ep.measure_distance(pop.d_true, meas, np.random.default_rng(int(seed)))
            agg = ep.Aggregates.build(pop.cell_id, d_meas, pop.on_true.astype(float))
        st.session_state["ep_last"] = (key, ref, eiv, nv, agg, pop.d_true, d_meas)
    key_, ref, eiv, nv, agg, d_true, d_meas = st.session_state["ep_last"]
    r, f, fn = ref["fit"], eiv["fit"], nv["fit"]
    lo, hi = eiv.get("contrast_ci", (np.nan, np.nan))
    c = st.columns(5)
    c[0].metric("Truth: P(ON) at 150 nm / at 1000 nm", f"{r.contrast():.2f}", help="From 3000 cells with the true distance and the true ON state.")
    c[1].metric("Recovered (errors accounted for)", f"{eiv['contrast']:.2f}", f"95% {lo:.2f}–{hi:.2f}", delta_color="off")
    c[2].metric("Recovered (errors ignored)", f"{nv['contrast']:.2f}")
    c[3].metric("Interval contains the truth?", "yes" if lo <= r.contrast() <= hi else "no")
    c[4].metric("Working distance: truth / recovered", f"{r.d_half:.0f} / {f.d_half:.0f} nm")

    g = ep.D_GRID
    bs = eiv.get("boot", [])
    band = np.array([b.g(g) for b in bs]) if bs else None

    def draw():
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [3, 2]})
        ax.plot(g, r.g(g), color="#e5e7eb", lw=2.4, label="truth (3000 cells, true distance and state)")
        if band is not None:
            ax.fill_between(g, np.quantile(band, 0.025, 0), np.quantile(band, 0.975, 0), color="#2dd4bf", alpha=0.25, lw=0)
        ax.plot(g, f.g(g), color="#2dd4bf", lw=2, label="recovered, errors accounted for (95% band)")
        ax.plot(g, fn.g(g), color="#f59e0b", lw=1.8, ls="--", label="recovered, errors ignored")
        nb, sb = agg.n.sum(0), agg.s.sum(0)
        centres = 0.5 * (ep.EDGES[:-1] + ep.EDGES[1:])
        # data points in 100 nm groups of measured distance
        grp = np.minimum((centres // 100).astype(int), 15)
        xs, ys, es = [], [], []
        for k in range(16):
            m = grp == k
            n_, s_ = nb[m].sum(), sb[m].sum()
            if n_ > 200:
                p_ = s_ / n_
                xs.append(100 * k + 50); ys.append(p_); es.append(np.sqrt(max(p_ * (1 - p_), 1e-9) / n_ * 8))     # x8: frames of a burst are not independent
        ax.errorbar(xs, ys, yerr=es, fmt="o", color="#60a5fa", ms=4, capsize=2, label="data, by MEASURED distance (error bars widened x2.8)")
        ax.set_xlim(0, 1500); ax.set_ylim(0, None)
        ax.set_xlabel("promoter-enhancer distance (nm)"); ax.set_ylabel("P(promoter in a burst)")
        ax.axvline(key_[0], color="#94a3b8", lw=0.8, ls=":")
        ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3)
        ax2.hist(d_true, bins=np.arange(0, 2000, 50), density=True, color="#94a3b8", alpha=0.6, label="true")
        ax2.hist(d_meas[np.isfinite(d_meas)], bins=np.arange(0, 2000, 50), density=True, histtype="step", color="#60a5fa", lw=1.8, label="measured")
        ax2.set_xlabel("distance (nm)"); ax2.set_ylabel("share of frames"); ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("ep_live", key_, len(bs)), draw)
    st.caption(f"{eiv['n_frames']:,} frames with a distance, {eiv['n_on']:.0f} of them ON. The ON state here is the TRUE one (level L2); the MS2-derived "
               "estimate is in the scans below. The 'ignored' curve falls short when the distance is measured in 2D or with large noise: both shorten "
               "or blur the distance, so near and far frames mix.")


def _scans() -> None:
    st.markdown("### 2 · Saved scans (many populations per setting)")
    files = sorted(SWEEPS.glob("ep_polymer_*/results.csv")) + sorted(SWEEPS.glob("ep_curve_*/results.csv"))
    if not files:
        st.info("No scans yet: scripts\\run.cmd python scripts\\sweep_ep_polymer.py --exps modality null wrong cells")
        return
    names = {("polymer: " if f.parent.name.startswith("ep_polymer_") else "invented: ") + f.parent.name.split("_", 2)[2]: f for f in files}
    LABELS = {"polymer: modality": "polymer distances: 2D / z-stack / biplane", "polymer: null": "polymer: no coupling (false positives)",
              "polymer: wrong": "polymer: wrong picks", "polymer: cells": "polymer: number of cells", "polymer: contact": "polymer: contact fraction",
              "polymer: frames": "polymer: frame interval", "polymer: colour": "polymer: delay between the two colours",
              "polymer: photons": "polymer: photons per locus", "polymer: focus": "polymer: focus error (2D, biplane)",
              "polymer: timescale": "polymer: which lag the D = 0.002 refers to", "polymer: compare_sources": "invented vs polymer distances, same measurement",
              "polymer: ms2": "polymer: MS2 noise x frame interval (ON from the MS2 trace)"}
    exp = st.selectbox("Scan", list(names), key="ep_exp", format_func=lambda n: LABELS.get(n, n))
    df = pd.read_csv(names[exp])
    cols = list(df.columns)
    at = cols.index("run")
    settings = [c for c in cols[:at] if c != "seed" and df[c].nunique(dropna=False) > 1]
    st.caption(f"{len(df)} analyses, {df.seed.nunique()} populations per setting. contrast = P(ON at 150 nm) / P(ON at 1000 nm); 'covers' = the 95% interval "
               "contains the reference value; 'excludes 1' = a detected effect.")
    if not settings:
        agg_cols = []
    else:
        agg_cols = settings
    tab = df.groupby(agg_cols).agg(populations=("seed", "size"), contrast=("contrast", "mean"), reference=("contrast_ref", "mean"), sd_over_populations=("contrast", "std"),
                                   interval_width=("contrast_ci_width", "mean"), covers=("contrast_covers", "mean"), excludes_1=("excludes_one", "mean"),
                                   curve_error=("curve_err", "mean")).reset_index() if agg_cols else df.agg({"contrast": "mean"}).to_frame().T
    st.dataframe(tab.round(3), hide_index=True, width="stretch")
    if not settings:
        return
    xcol = st.selectbox("x axis", settings, key="ep_x")
    others = [c for c in settings if c != xcol]

    def draw():
        fig, ax = plt.subplots(figsize=(9, 4.6))
        groups = list(df.groupby(others)) if others else [((), df)]
        xs_all = sorted(df[xcol].unique(), key=lambda v: (isinstance(v, str), v))
        num = all(not isinstance(v, (str, bool)) for v in xs_all)
        for gi, (k, g) in enumerate(groups):
            k = k if isinstance(k, tuple) else (k,)
            m = g.groupby(xcol)["contrast"].agg(["mean", "std"]).reindex(xs_all)
            xv = xs_all if num else list(range(len(xs_all)))
            ax.errorbar(xv, m["mean"], yerr=m["std"].fillna(0), marker="o", capsize=3, lw=1.8, color=style.CYCLE[gi % len(style.CYCLE)],
                        label=", ".join(f"{c}={v}" for c, v in zip(others, k)) or "contrast")
        ref = df.contrast_ref.mean()
        ax.axhline(ref, color="#e5e7eb", lw=1.2, ls="--", label=f"reference {ref:.2f}")
        ax.axhline(1.0, color="#94a3b8", lw=0.8, ls=":")
        if not num:
            ax.set_xticks(range(len(xs_all))); ax.set_xticklabels([str(v) for v in xs_all])
        if num and xs_all[-1] / max(xs_all[0], 1e-9) > 20:
            ax.set_xscale("log")
        ax.set_xlabel(xcol); ax.set_ylabel("recovered contrast (mean ± sd over populations)"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        fig.tight_layout()
        return fig

    style.cached_figure(st, ("ep_scan", exp, xcol, len(df), float(df.contrast.sum())), draw)
