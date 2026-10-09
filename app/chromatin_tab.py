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

    lib = st.selectbox("Library", libs)
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
        st.caption("Every bead is coloured by genomic position (purple = start of the region, yellow = end) and the chain "
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
    fig, ax = plt.subplots(1, 2, figsize=(15, 4))
    ax[0].plot(dist, lw=1)
    ax[0].set_xlabel("block"); ax[0].set_ylabel("distance (monomer diameters)")
    ax[0].set_title(f"Loci at {beads[ia]} kb and {beads[ib]} kb (separation {abs(beads[ib] - beads[ia])} kb)")
    ax[1].hist(dist, bins=40, color="tab:blue")
    ax[1].set_xlabel("distance (monomer diameters)"); ax[1].set_title(f"distribution: mean {dist.mean():.1f}, sd {dist.std():.1f}")
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    lef_f = d / "lef_positions.npy"
    if lef_f.exists():
        lp = np.load(lef_f, mmap_mode="r")
        st.subheader("Loop-extruding factors over time")
        fig, ax = plt.subplots(figsize=(14, 3.6))
        for k in range(min(lp.shape[1], 60)):
            ax.plot(np.arange(lp.shape[0]), lp[:, k, 0], lw=0.4, color="tab:blue")
            ax.plot(np.arange(lp.shape[0]), lp[:, k, 1], lw=0.4, color="tab:red")
        ax.set_xlabel("block"); ax.set_ylabel("monomer index (kb)")
        ax.set_title("Left (blue) and right (red) leg of 60 LEFs: loops grow, stall at boundaries, unbind and reload")
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)
