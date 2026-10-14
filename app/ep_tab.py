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
from simlive.stage7_validation import ep_curve as ep

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
    _live()
    st.divider()
    _scans()


def _live() -> None:
    st.markdown("### 1 · One population, your errors")
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
    files = sorted(SWEEPS.glob("ep_curve_*/results.csv"))
    if not files:
        st.info("No scans yet: `scripts\\run.cmd python scripts\\sweep_ep_curve.py --exps sigma shift outliers cells null`")
        return
    names = {f.parent.name.replace("ep_curve_", ""): f for f in files}
    exp = st.selectbox("Scan", list(names), key="ep_exp",
                       format_func=lambda n: {"ladder": "ladder: truth → estimated ON → measured distance → both", "sigma": "localization precision × 2D/3D",
                                              "shift": "colour shift", "outliers": "wrong picks / missed frames", "cells": "number of cells",
                                              "null": "no coupling (false positives)", "dynamics": "time per visit", "contact": "contact fraction",
                                              "frames": "frame interval", "ms2": "MS2 noise × frame interval"}.get(n, n))
    df = pd.read_csv(names[exp])
    cols = list(df.columns)
    at = cols.index("run")
    settings = [c for c in cols[:at] if c != "seed" and df[c].nunique() > 1]
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
