"""The dashboard's 🔬 Cells tab: look at every isolated cell, one cell through time, and its loci.

Reads only files written by stage 5 / stage 6 (data/runs/<run>/stage5_cells/<source>/..., stage6_analysis/...).
Only single TIFF pages (frames) are read, so memory use stays small however many cells or frames a run has.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as _st
import style
import tifffile
from PIL import Image

from fastview import all_cells_frames, cell_composite, single_cell_frames
from player import show_player
from simlive.stage5_linking.celldata import load_cell

SRC_LABEL = {"truth": "True cell IDs (ground truth)", "tracked_gt_masks": "Tracked IDs (tracker on ground-truth masks)",
             "tracked_auto_masks": "Tracked IDs (tracker on automatic segmentation)"}
COLORS = ["tab:red", "tab:green", "tab:blue"]
COLOR_NAMES = ["red", "green", "blue"]


def read_frame(path: Path, i: int) -> np.ndarray:
    """One page (frame) of a TIFF stack, without loading the whole stack."""
    return tifffile.imread(path, key=int(i)).astype(float)


def cell_rgb(cdir: Path, i: int, channels: list) -> np.ndarray:
    """Composite: nucleus in grey, locus channels in red / green / blue."""
    nuc = read_frame(cdir / "nucleus.tif", i)
    nuc = np.clip((nuc - np.percentile(nuc, 1)) / max(np.percentile(nuc, 99.5) - np.percentile(nuc, 1), 1e-9), 0, 1)
    rgb = np.repeat(nuc[..., None] * 0.35, 3, axis=2)
    for k, ch in enumerate(channels[:3]):
        im = read_frame(cdir / f"{ch}.tif", i)
        # a diffraction-limited spot covers only a few pixels, so scale to the true maximum (a percentile lands on noise)
        im = np.clip((im - np.median(im)) / max(im.max() - np.median(im), 1e-9), 0, 1)
        rgb[..., k] = np.clip(rgb[..., k] + im, 0, 1)
    return rgb


@_st.cache_resource(show_spinner=False, max_entries=3)
def _all_cells_movie(sdir: str, stamp: float, loc_names: tuple, tile: int = 150):
    """Every cell as a small RGB tile per frame (computed once, one cell in memory at a time)."""
    summ = pd.read_csv(Path(sdir) / "cell_summary.csv")
    data = []
    for cid in summ.cell_id:
        c = load_cell(sdir, cid)
        tiles = {}
        for i, t in enumerate(c.t):
            rgb = cell_composite(c.channels["nucleus"][i], [c.channels[n][i] for n in loc_names])
            tiles[int(t)] = np.asarray(Image.fromarray(rgb.astype(np.uint8)).resize((tile, tile), Image.BILINEAR))
        data.append((int(cid), tiles))
    return all_cells_frames(data, int(summ.last_t.max()) + 1, tile=tile)


def render(st, run_path: Path, cfg: dict):
    style.intro(st, "Every cell isolated into its own small movie (nucleus plus one colour channel per locus), so each cell can be analysed on its own.", "Choose which cell identities to use (true or tracker), look at all cells at one frame, then pick one cell to watch through time with its located loci.")
    s5 = run_path / "stage5_cells"
    sources = [p.name for p in sorted(s5.glob("*")) if (p / "cell_summary.csv").exists()] if s5.exists() else []
    if not sources:
        st.info("This run has no single-cell data. Make one with the locus config:\n\n"
                "`scripts\\run.cmd python scripts\\run_tracking_demo.py --config configs\\loci_demo.yaml`")
        return
    src = st.radio("Cells identified by", sources, index=sources.index("truth") if "truth" in sources else 0,
                   format_func=lambda s: SRC_LABEL.get(s, s), horizontal=True)
    sdir = s5 / src
    summ = pd.read_csv(sdir / "cell_summary.csv")
    px_um = cfg["optics"]["pixel_size_nm"] / 1000.0
    S = int(summ.crop_size_px.iloc[0])
    first = sdir / "cells" / f"cell_{int(summ.cell_id.iloc[0]):04d}"
    loc_ch = sorted(p.stem for p in first.glob("locus*.tif"))
    pos_f, dist_f = (run_path / "stage6_analysis" / src / n for n in ("loci_positions.csv", "locus_distances.csv"))
    pos = pd.read_csv(pos_f) if pos_f.exists() else pd.DataFrame()
    dist = pd.read_csv(dist_f) if dist_f.exists() else pd.DataFrame()
    tr_f = run_path / "stage1_chromatin" / "loci_truth.csv"
    truth = pd.read_csv(tr_f) if (src == "truth" and tr_f.exists()) else None
    T_all = int(summ.last_t.max()) + 1

    st.caption(f"{len(summ)} cells isolated. Each cell is a fixed {S}×{S} px window ({S * px_um:.1f} µm) that follows "
               "its centroid. Channels: nucleus (grey)" + "".join(
                   f", {c} ({COLOR_NAMES[k]})" for k, c in enumerate(loc_ch[:3])) + ".")
    with st.expander("Cell summary table (lifetime, gaps, border contact, closest neighbour)"):
        st.dataframe(summ.drop(columns=["crop_size_px"]).round(1), hide_index=True, width="stretch")

    # ---- smooth movie of all cells side by side (in-browser player)
    with st.expander("▶ Play all cells together (smooth, runs in your browser)", expanded=st.session_state.get("cells_all_play", False)):
        if st.toggle("Prepare player", key="cells_all_play", help="Takes a few seconds the first time for each run"):
            with st.spinner("Rendering frames…"):
                frames_all, labels_all = _all_cells_movie(str(sdir), (sdir / "cell_summary.csv").stat().st_mtime,
                                                          tuple(loc_ch))
            show_player(frames_all, labels_all, fps=10, column_width_px=700, key="allcells")

    # ---- contact sheet: every cell at one frame
    st.subheader("All cells at one frame")
    ft = st.slider("Frame", 0, T_all - 1, min(10, T_all - 1), key="cells_frame")
    alive = []
    for r in summ.itertuples():
        fr = pd.read_csv(sdir / "cells" / f"cell_{r.cell_id:04d}" / "frames.csv", usecols=["t"])
        idx = np.nonzero(fr.t.to_numpy() == ft)[0]
        if len(idx):
            alive.append((r.cell_id, int(idx[0])))
    if not alive:
        st.info("No cells exist at this frame.")
    else:
        ncol = min(6, len(alive))
        nrow = int(np.ceil(len(alive) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.7 * nrow), squeeze=False)
        for ax in axes.ravel():
            ax.axis("off")
        for ax, (cid, i) in zip(axes.ravel(), alive):
            ax.imshow(cell_rgb(sdir / "cells" / f"cell_{cid:04d}", i, loc_ch), interpolation="nearest")
            ax.set_title(f"cell {cid}", fontsize=10)
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)

    # ---- one cell in detail
    st.subheader("One cell through time")
    by_id = summ.set_index("cell_id")
    cid = st.selectbox("Cell", summ.cell_id.tolist(), format_func=lambda c: f"cell {c}  ({int(by_id.n_frames[c])} frames)")
    row = by_id.loc[cid]
    cdir = sdir / "cells" / f"cell_{cid:04d}"
    frames = pd.read_csv(cdir / "frames.csv")
    flags = []
    if row.n_missing_frames > 0:
        flags.append(f"{int(row.n_missing_frames)} missing frame(s) (gap in the track)")
    if row.border_frames > 0:
        flags.append(f"touches the image border in {int(row.border_frames)} frame(s)")
    if row.n_frames < 0.5 * T_all:
        flags.append(f"short track ({int(row.n_frames)} of {T_all} frames)")
    if flags:
        st.warning("; ".join(flags))
    else:
        st.success(f"Clean track: present in {int(row.n_frames)} frames, no gaps, never touching the border.")
    with st.expander("▶ Play this cell's movie (smooth, runs in your browser)", expanded=st.session_state.get("cell_play", False)):
        if st.toggle("Prepare player", key="cell_play", help="Renders this cell's frames once"):
            with st.spinner("Rendering frames…"):
                cd = load_cell(sdir, cid)
                pc_ = pos[pos.cell_id == cid] if len(pos) else None
                fr_list, lb_list = single_cell_frames(cd.channels, cd.mask, cd.frames, loc_ch, pc_, truth, px_um)
            show_player(fr_list, lb_list, fps=10, column_width_px=800, key="onecell")
    ci = st.slider("Frame (this cell)", 0, len(frames) - 1, 0, key="cell_idx")
    t_now, fr_now = int(frames.t.iloc[ci]), frames.iloc[ci]
    mask = read_frame(cdir / "mask.tif", ci) > 0
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    axes[0].imshow(read_frame(cdir / "nucleus.tif", ci), cmap="gray")
    axes[0].set_title(f"Nucleus (frame {t_now})")
    axes[1].imshow(cell_rgb(cdir, ci, loc_ch), interpolation="nearest")
    axes[1].set_title("Loci (colours) on nucleus (grey)")
    axes[2].imshow(mask, cmap="gray")
    axes[2].set_title("This cell's mask (neighbours removed)")
    for ax in axes:
        ax.contour(mask, levels=[0.5], colors="yellow", linewidths=0.8)
        ax.axis("off")
    if len(pos):
        for q in pos[(pos.cell_id == cid) & (pos.t == t_now)].itertuples():
            axes[1].plot(q.crop_x_px, q.crop_y_px, "+", color="white", ms=16, mew=2)
            axes[1].annotate(f"locus{q.locus_id}", (q.crop_x_px, q.crop_y_px), color="white", fontsize=9,
                             xytext=(6, -8), textcoords="offset points")
    if truth is not None:
        for q in truth[(truth.cell_id == cid) & (truth.t == t_now)].itertuples():
            axes[1].plot(q.x_um / px_um - fr_now.crop_x0, q.y_um / px_um - fr_now.crop_y0, "o", mfc="none",
                         mec="cyan", ms=14, mew=1.5)
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    st.caption("White + = locus found by the analysis" + ("; cyan ○ = true locus position (simulation)." if truth is not None
               else ". (True positions are only available for true cell IDs.)"))

    if len(pos):
        pc = pos[pos.cell_id == cid]
        cen = frames.set_index("t")
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))
        for k in sorted(pc.locus_id.unique()):
            q = pc[pc.locus_id == k]
            ex = (q.lab_x_px.to_numpy() - cen.loc[q.t, "cx_px"].to_numpy()) * px_um
            ey = (q.lab_y_px.to_numpy() - cen.loc[q.t, "cy_px"].to_numpy()) * px_um
            axes[0].plot(q.t, ex, color=COLORS[k % 3], label=f"locus{k} (estimated)")
            axes[1].plot(q.t, ey, color=COLORS[k % 3])
            if truth is not None:
                tq = truth[(truth.cell_id == cid) & (truth.locus_id == k)].set_index("t")
                tt = tq.index.intersection(cen.index)
                axes[0].plot(tt, tq.loc[tt, "x_um"] - cen.loc[tt, "cx_px"] * px_um, "k:", lw=1)
                axes[1].plot(tt, tq.loc[tt, "y_um"] - cen.loc[tt, "cy_px"] * px_um, "k:", lw=1)
        axes[0].set_title("Locus x relative to cell centre (µm)")
        axes[1].set_title("Locus y relative to cell centre (µm)")
        axes[0].legend(fontsize=8)
        dd = dist[dist.cell_id == cid] if len(dist) else dist
        if len(dd):
            axes[2].plot(dd.t, dd.dist_um, color="k", label="estimated")
        if truth is not None:
            a = truth[(truth.cell_id == cid) & (truth.locus_id == 0)].set_index("t")
            b = truth[(truth.cell_id == cid) & (truth.locus_id == 1)].set_index("t")
            j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
            axes[2].plot(j.index, np.hypot(j.y_um_a - j.y_um_b, j.x_um_a - j.x_um_b), "r:", label="true")
        axes[2].set_title("Distance between locus 0 and locus 1 (µm)")
        axes[2].legend(fontsize=8)
        for ax in axes:
            ax.axvline(t_now, color="orange", lw=1)
            ax.set_xlabel("frame")
            ax.grid(alpha=0.3)
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)
        st.caption("Dotted lines = true values from the simulation (true cell IDs only). Orange line = frame shown above.")
    sc_f = run_path / "stage6_analysis" / "scores.json"
    if sc_f.exists():
        with st.expander("Analysis quality and memory use for this run"):
            st.json(json.loads(sc_f.read_text()))
