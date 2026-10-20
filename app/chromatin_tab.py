"""The dashboard's 🧬 Chromatin tab: validation report, and a browser for saved chromatin simulations.

Reads only files written by scripts/run_chromatin.py / validate_chromatin.py (DATA_ROOT/chromatin/<library>/<sim>/).
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as _st

import chromatin_view as V
from player import show_player
from simlive.io.runs import DATA_ROOT, REPO_ROOT
from simlive.stage1_chromatin import analysis as A

import style
CHROM = DATA_ROOT / "chromatin"


@_st.cache_resource(show_spinner=False, max_entries=3)
def _watch_frames(sim_dir: str, mode: str, a_kb: int, b_kb: int, show_loops: bool, rotate_deg: float):
    """Pre-render the movie frames for one simulation (cached per setting)."""
    fn = V.snapshot_frames if mode == "snapshots" else V.block_frames
    return fn(sim_dir, a_kb, b_kb, show_loops, rotate_deg)


@_st.cache_resource(show_spinner=False, max_entries=4)
def _maps(sim_dir: str, n_snap: int):
    """Contact map + P(s) from the last n_snap saved snapshots of one simulation."""
    snaps = A.load_snapshots(sim_dir)[-n_snap:]
    m = A.contact_map(snaps, bin_size=10 if len(snaps[0]) < 20000 else 20)
    s, p = A.contact_probability(A.separation_histogram(snaps))
    return m, (s, p), snaps[-1], len(snaps)


def render(st):
    style.intro(st, "Look inside one saved polymer simulation: the chain moving through time, its contact map, and the distance between any two chosen loci.", "Pick a library and simulation, tick the player to watch the chromosome, then choose two genomic positions (kb) to plot their distance over time. Those distances drive the transcription model.")
    libs = sorted(p.name for p in CHROM.glob("*") if p.is_dir()) if CHROM.exists() else []
    if not libs:
        st.info("No chromatin simulations yet. Quick test (1 minute):\n\n"
                "`scripts\\run.cmd python scripts\\run_chromatin.py --config configs\\chromatin\\smoke_test.yaml --seed 1 --library smoke`")
        return

    rep = REPO_ROOT / "docs" / "chromatin_validation.md"
    with st.expander("✅ Validation report (reproducing known Mirny-lab behaviour)", expanded=False):
        if rep.exists():
            txt = rep.read_text(encoding="utf-8")
            st.markdown(txt.split("## Figures")[0])
            for name in ("contact_maps", "contact_probability", "insulation", "msd"):
                f = REPO_ROOT / "docs" / "img" / "chromatin" / f"{name}.png"
                if f.exists():
                    st.image(str(f))
        else:
            st.info("The validation report has not been generated yet (`scripts\\validate_chromatin.py all`).")

    pick = st.session_state.get("lib_sel")              # open on the library chosen in the 🗂 Library tab, if any
    lib = st.selectbox("Library", libs, index=libs.index(pick) if pick in libs else 0)
    sims = sorted(p for p in (CHROM / lib).glob("*_seed*") if (p / "meta.json").exists())
    if not sims:
        st.warning("No finished simulations in this library yet.")
        return
    rows = []
    for d in sims:
        m = json.loads((d / "meta.json").read_text())
        rows.append({"simulation": d.name, "monomers (kb)": m["n_monomers"], "MD steps/s": round(m["md_steps_per_second"]),
                     "minutes": round(m["wall_seconds_total"] / 60, 1), "mean loop (kb)": m.get("mean_loop_size_kb"),
                     "boundary sites": m.get("n_ctcf_sites")})
    st.caption(f"{len(sims)} finished simulations in `{lib}`.")
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")

    name = st.selectbox("Simulation", [d.name for d in sims])
    d = CHROM / lib / name
    meta = json.loads((d / "meta.json").read_text())
    N = meta["n_monomers"]
    cfg_p = d / "params.yaml"
    c1, c2, c3 = st.columns(3)
    c1.metric("Monomers", f"{N:,}", help="1 monomer = 1 kb")
    c2.metric("Speed", f"{meta['md_steps_per_second']:.0f} steps/s")
    c3.metric("Mean loop size", f"{meta['mean_loop_size_kb']:.0f} kb" if "mean_loop_size_kb" in meta else "no loop extrusion")

    has_lef = (d / "lef_positions.npy").exists()
    with st.expander("▶ Watch the chromosome move through time (smooth player)",
                     expanded=st.session_state.get("chrom_play", False)):
        st.caption("**The grey circle is NOT the nucleus.** It is the sphere that confines this simulated region of chromosome (a computational box: the "
                   "chain fills it by construction because the density is fixed). In a movie this territory is scaled by 'nm per polymer unit' and "
                   "placed inside a nucleus several times larger, so it covers only a small part of it. "
                   "Every bead is coloured by genomic position (purple = start of the region, yellow = end) and the chain "
                   "slowly rotates so the 3D shape is visible. Two loci you choose are highlighted (A red, B blue) and "
                   "joined by a line with their distance. White lines = loop-extruding factors holding two monomers together. "
                   "Time is in MD blocks (750 steps each) until the calibration step converts it to seconds.")
        c1, c2, c3 = st.columns([1, 3, 2])
        on = c1.toggle("Prepare player", key="chrom_play", help="Renders every frame once (a few seconds)")
        mode = c2.radio("Time resolution", ["snapshots", "blocks"], horizontal=True, key="chrom_mode",
                        format_func={"snapshots": "Whole chain at each saved snapshot (every monomer)",
                                     "blocks": "Every block, finest time resolution (1 in 10 monomers)"}.get)
        rot = c3.slider("Rotation per frame (degrees)", 0.0, 5.0, 1.5, 0.5, key="chrom_rot")
        wa = st.slider("Locus A position (kb)", 0, int(N) - 1, int(N * 0.40), 10, key="watch_a")
        ws = st.slider("Locus B is this far from A (kb)", 10, int(min(3000, N - 1)), 500, 10, key="watch_sep")
        show_loops = st.checkbox("Draw loop-extruding factors", True, key="watch_loops", disabled=not has_lef)
        if on:
            with st.spinner("Rendering frames…"):
                frames, labels = _watch_frames(str(d), mode, int(wa), int(min(N - 1, wa + ws)), bool(show_loops and has_lef),
                                               float(rot))
            show_player(frames, labels, fps=10, column_width_px=640, key="chrom")

    if st.toggle("Show contact map, P(s) and 3D shape (reads saved snapshots; a few seconds)", key="chrom_maps"):
        n_snap = 20
        with st.spinner("Computing from snapshots…"):
            m, (s, p), last, used = _maps(str(d), n_snap)
        ctcf = pd.read_csv(d / "ctcf_sites.csv") if (d / "ctcf_sites.csv").exists() else None
        fig = plt.figure(figsize=(16, 5))
        ax = fig.add_subplot(1, 3, 1)
        ax.imshow(np.log10(m + 1e-4), cmap="Reds", vmin=-2.5, vmax=0.3)
        ax.set_title(f"Contact map (log10), last {used} snapshots")
        ax.set_xlabel("genomic bin"); ax.set_ylabel("genomic bin")
        ax = fig.add_subplot(1, 3, 2)
        ax.loglog(s, p, "o-", ms=3)
        ax.set_xlabel("genomic separation s (kb)"); ax.set_ylabel("P(s)"); ax.set_title("Contact probability"); ax.grid(alpha=0.3)
        ax = fig.add_subplot(1, 3, 3, projection="3d")
        step = max(1, N // 3000)
        ax.scatter(last[::step, 0], last[::step, 1], last[::step, 2], c=np.arange(0, N, step), cmap="viridis", s=3)
        ax.set_title("Last conformation (colour = genomic position)")
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)
        if ctcf is not None:
            st.caption(f"{len(ctcf)} boundary sites (synthetic positions) at kb: " + ", ".join(map(str, ctcf.pos.tolist()[:20])) +
                       (" …" if len(ctcf) > 20 else ""))

    st.subheader("Two loci: distance over time")
    st.caption("What the microscopy step will use: pick two genomic positions (kb), and see how far apart those two "
               "monomers are in every saved block. Time is in MD blocks until the calibration step converts it to seconds.")
    beads = np.load(d / "tracked_beads.npy")
    traj = np.load(d / "tracked_positions.npy", mmap_mode="r")
    a_kb = st.slider("Locus A position (kb)", 0, int(N) - 1, int(N * 0.40), 10)
    sep = st.slider("Separation to locus B (kb)", 10, 3000, 500, 10)
    b_kb = min(N - 1, a_kb + sep)
    ia, ib = int(np.argmin(np.abs(beads - a_kb))), int(np.argmin(np.abs(beads - b_kb)))
    dist = np.linalg.norm(np.asarray(traj[:, ia]) - np.asarray(traj[:, ib]), axis=1)
    try:                                   # the assumed calibration (configs/calibration/default.json): polymer units -> nm, blocks -> s
        from simlive.calibration.polymer import default_scale

        nm_u, s_blk, _ = default_scale()
    except (OSError, KeyError, ValueError):  # unreadable calibration file: polymer units only
        nm_u, s_blk = None, None
    fig, ax = plt.subplots(1, 2, figsize=(15, 4))
    ax[0].plot(dist, lw=0.8, color="#2dd4bf")
    ax[0].set_xlabel("block"); ax[0].set_ylabel("distance (monomer diameters)")
    ax[0].set_title(f"Loci at {beads[ia]} kb and {beads[ib]} kb ({abs(beads[ib] - beads[ia])} kb apart)")
    ax[1].hist(dist, bins=40, color="#60a5fa", alpha=0.8)
    ax[1].set_xlabel("distance (monomer diameters)")
    ax[1].set_title(f"mean {dist.mean():.1f}, sd {dist.std():.1f} monomer diameters"
                    + (f" = {dist.mean() * nm_u / 1000:.2f} ± {dist.std() * nm_u / 1000:.2f} µm" if nm_u else ""))
    if nm_u:                               # the same axes in physical units (assumed calibration)
        sec = ax[0].secondary_yaxis("right", functions=(lambda v: v * nm_u / 1000, lambda v: v * 1000 / nm_u))
        sec.set_ylabel("µm (assumed calibration)")
        top = ax[0].secondary_xaxis("top", functions=(lambda b: b * s_blk / 60, lambda m: m * 60 / s_blk))
        top.set_xlabel("minutes")
        ax[1].secondary_xaxis("top", functions=(lambda v: v * nm_u / 1000, lambda v: v * 1000 / nm_u)).set_xlabel("µm")
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    lef_f = d / "lef_positions.npy"
    if lef_f.exists():
        _loop_extruders(st, d, N)


@_st.cache_resource(show_spinner=False, max_entries=3)
def _coverage(sim_dir: str, n_monomers: int, bin_kb: int):
    lp = np.load(Path(sim_dir) / "lef_positions.npy")
    return A.loop_coverage(lp, n_monomers, bin_kb), lp.shape


def _loop_extruders(st, d: Path, N: int) -> None:
    """Loop-extruding factors over time, in two readable views: (1) the whole simulation as a picture of how many loops cover each genomic
    position over time (hundreds of overlapping lines were unreadable); (2) a zoom into a genomic window and a stretch of time that draws each
    factor's loop as a shaded band between its two legs."""
    from simlive.calibration.polymer import default_scale

    st.subheader("Loop-extruding factors over time")
    lp = np.load(d / "lef_positions.npy", mmap_mode="r")
    T, K, _ = lp.shape
    _, block_s, _ = default_scale()
    ctcf = pd.read_csv(d / "ctcf_sites.csv").pos.to_numpy() if (d / "ctcf_sites.csv").exists() else np.array([])
    st.caption(f"{K} factors over {T:,} blocks (about {T * block_s / 60:.0f} min with the default calibration of {block_s:.2f} s per block, an assumption). "
               "**Top:** how many loops cover each genomic position (colour) at each time, the whole simulation at a glance: bright horizontal "
               "bands are positions that are often inside a loop, dark gaps are places loops rarely reach (boundary sites are marked on the right). "
               "**Middle:** the average loop size. **Below:** zoom into one region and stretch of time to see individual factors.")
    bin_kb = 20
    cov, _ = _coverage(str(d), int(N), bin_kb)
    sizes = np.asarray(lp[:, :, 1] - lp[:, :, 0]).mean(axis=1)
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(14, 6.2), sharex=True, gridspec_kw={"height_ratios": [3.2, 1], "hspace": 0.12})
    im = ax.imshow(cov.T, aspect="auto", origin="lower", extent=[0, T, 0, N], cmap="magma", vmin=0, vmax=max(1.0, float(np.percentile(cov, 99.5))),
                   interpolation="nearest")
    if len(ctcf):
        ax.scatter(np.full(len(ctcf), T * 1.006), ctcf, marker="<", s=26, color=style.PALETTE["accent"], clip_on=False, zorder=5, label="boundary site")
        ax.legend(loc="upper left", fontsize=9, frameon=True, facecolor=style.PALETTE["fig"], edgecolor="none")
    ax.set_ylabel("genomic position (kb)"); ax.grid(False)
    ax.set_title("Loops covering each genomic position over time", pad=34)       # the minutes axis sits between the plot and the title
    fig.colorbar(im, ax=[ax, ax2], pad=0.02, fraction=0.025, label="loops covering a 20 kb bin")
    ax2.plot(np.arange(T), sizes, lw=0.8, color=style.CYCLE[0])
    ax2.set_ylabel("mean loop (kb)"); ax2.set_xlabel("block"); ax2.set_ylim(0, max(1.0, float(sizes.max()) * 1.15)); ax2.grid(alpha=0.3)
    top = ax.secondary_xaxis("top", functions=(lambda b: b * block_s / 60, lambda m: m * 60 / block_s))
    top.set_xlabel(f"minutes (assumed {block_s:.2f} s per block)")
    st.pyplot(fig, width="stretch")
    plt.close(fig)

    st.markdown("**Zoom: individual factors**")
    c = st.columns(4)
    width = c[0].select_slider("Genomic window (kb)", [200, 300, 400, 600, 1000, 2000], 400, key="lef_w")
    g0 = c[1].slider("Window starts at (kb)", 0, int(max(0, N - width)), int(min(max(0, N - width), N * 0.4)), 50, key="lef_g0")
    span = c[2].select_slider("Time shown (blocks)", [100, 200, 300, 500, 800, 1500], 300, key="lef_span")
    t0 = c[3].slider("Time starts at (block)", 0, int(max(0, T - min(span, T))), 0, 50, key="lef_t0")
    t1 = int(min(T, t0 + span))
    n_all = len(A.loop_traces(lp, t0, t1, g0, g0 + width, max_factors=10 ** 6))
    traces = A.loop_traces(lp, t0, t1, g0, g0 + width, max_factors=16)
    fig, ax = plt.subplots(figsize=(14, 4.4))
    for i, tr in enumerate(traces):
        col = style.CYCLE[i % len(style.CYCLE)]
        ax.fill_between(tr["t"], tr["left"], tr["right"], color=col, alpha=0.16, linewidth=0)
        ax.plot(tr["t"], tr["left"], lw=0.9, color=col)
        ax.plot(tr["t"], tr["right"], lw=0.9, color=col)
    for p in ctcf:
        if g0 <= p <= g0 + width:
            ax.axhline(p, color="#e5e7eb", lw=0.8, ls="--", alpha=0.6)
            ax.text(t1, p, " boundary", va="center", ha="left", fontsize=8, color="#e5e7eb", clip_on=False)
    ax.set_xlim(t0, t1); ax.set_ylim(g0, g0 + width)
    ax.set_xlabel("block"); ax.set_ylabel("genomic position (kb)")
    shown = f"{len(traces)} of {n_all} loop-extruding factors (the busiest)" if n_all > len(traces) else f"{n_all} loop-extruding factors"
    ax.set_title(f"{shown}: each band is one loop between its two legs")
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    st.caption("A band widens as the factor extrudes (both legs move apart), stops growing at a boundary (dashed line) and ends where the factor unbinds "
               "and reloads elsewhere (the line breaks). Narrow the window or time range if bands overlap.")
