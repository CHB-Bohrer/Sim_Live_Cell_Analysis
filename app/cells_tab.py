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
from style import PALETTE

PALETTE_WARN = PALETTE["warn"]
import tifffile
from PIL import Image

from fastview import all_cells_frames, cell_composite, movie_limits, single_cell_frames
from player import show_player
from simlive.stage5_linking.celldata import load_cell

SRC_LABEL = {"truth": "True cell IDs (ground truth)", "tracked_gt_masks": "Tracked IDs (tracker on ground-truth masks)",
             "tracked_auto_masks": "Tracked IDs (tracker on automatic segmentation)"}
COLORS = ["#f87171", "#4ade80", "#60a5fa"]          # locus colours as in the images (red, green, blue), light enough for the dark theme
COLOR_NAMES = ["red", "green", "blue"]
COPY_COLORS = ["#e5e7eb", "#2dd4bf", "#f59e0b", "#a78bfa"]       # one colour per promoter-enhancer copy
LOCUS_LETTER = {0: "P", 1: "E"}                                    # promoter, enhancer


def tag(locus_id, copy) -> str:
    """Short name of one spot: P0 = promoter of copy 0, E1 = enhancer of copy 1 (copies >= 50: found copies that match no true copy)."""
    c = int(copy)
    return f"{LOCUS_LETTER.get(int(locus_id), f'L{int(locus_id)}')}{c if c < 50 else '?'}"


def with_copy(df):
    return df if df is None or not len(df) or "copy" in df.columns else df.assign(copy=0)


def read_frame(path: Path, i: int) -> np.ndarray:
    """One page (frame) of a TIFF stack, without loading the whole stack."""
    return tifffile.imread(path, key=int(i)).astype(float)


def cell_rgb(cdir: Path, i: int, channels: list) -> np.ndarray:
    """Composite: nucleus in grey, locus channels in red / green / blue."""
    nuc = read_frame(cdir / "nucleus.tif", i)
    nuc = np.clip((nuc - np.percentile(nuc, 1)) / max(np.percentile(nuc, 99.5) - np.percentile(nuc, 1), 1e-9), 0, 1)
    rgb = np.repeat(nuc[..., None] * 0.35, 3, axis=2)
    from scipy.ndimage import binary_dilation

    region = binary_dilation(read_frame(cdir / "mask.tif", i) > 0, iterations=4)        # scale from the cell's own region: beads, cosmic rays elsewhere are ignored
    for k, ch in enumerate(channels[:3]):
        im = read_frame(cdir / f"{ch}.tif", i)
        v = im[region] if region.any() else im
        lo, hi = np.median(v) + 1.5 * 1.4826 * np.median(np.abs(v - np.median(v))), max(v.max(), np.median(v) + 1)
        im = np.clip((im - lo) / max(hi - lo, 1e-9), 0, 1)
        rgb[..., k] = np.clip(rgb[..., k] + im, 0, 1)
    return rgb


@_st.cache_resource(show_spinner=False, max_entries=3)
def _all_cells_movie(sdir: str, stamp: float, loc_names: tuple, tile: int = 150):
    """Every cell as a small RGB tile per frame (computed once, one cell in memory at a time)."""
    summ = pd.read_csv(Path(sdir) / "cell_summary.csv")
    data = []
    for cid in summ.cell_id:
        c = load_cell(sdir, cid)
        lims = [movie_limits(c.channels[n], c.mask) for n in loc_names]      # one scale per colour for this cell's whole movie, from its own region
        tiles = {}
        for i, t in enumerate(c.t):
            rgb = cell_composite(c.channels["nucleus"][i], [c.channels[n][i] for n in loc_names], lims=lims)
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
    tr_f = run_path / "stage3_microscopy" / "loci_truth_image.csv"       # truth in image coordinates (moved by the stage drift)
    if not tr_f.exists():
        tr_f = run_path / "stage1_chromatin" / "loci_truth.csv"          # runs made before the image-error model
    truth = with_copy(pd.read_csv(tr_f)) if (src == "truth" and tr_f.exists()) else None
    pos, dist = with_copy(pos), with_copy(dist)
    sp_f = run_path / "stage3_microscopy" / "stray_spots.csv"
    strays = pd.read_csv(sp_f) if sp_f.exists() else None             # decoy spots (image error source stray_spots), image coordinates
    err_f = run_path / "stage6_analysis" / src / "localization_errors.csv"
    errs = with_copy(pd.read_csv(err_f)) if (src == "truth" and err_f.exists()) else None
    T_all = int(summ.last_t.max()) + 1
    n_cp = pd.read_csv(run_path / "stage1_chromatin" / "copies.csv").set_index("cell_id").n_copies if (run_path / "stage1_chromatin" / "copies.csv").exists() else None

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
    cid = st.selectbox("Cell", summ.cell_id.tolist(), format_func=lambda c: f"cell {c}  ({int(by_id.n_frames[c])} frames"
                       + (f", {int(n_cp[c])} promoter-enhancer cop{'y' if int(n_cp[c]) == 1 else 'ies'}" if n_cp is not None and src == "truth" and c in n_cp.index else "") + ")")
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
                tc_ = truth[truth.cell_id == cid] if truth is not None else None      # this cell's loci only (neighbours' loci can fall in the crop)
                sc_ = strays[strays.cell_id == cid] if strays is not None and src == "truth" else None
                fr_list, lb_list = single_cell_frames(cd.channels, cd.mask, cd.frames, loc_ch, pc_, tc_, px_um, strays=sc_)
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
            left = int(q.locus_id) % 2 == 0           # alternate sides (and heights) so the labels of two close loci do not overlap
            axes[1].annotate(tag(q.locus_id, q.copy), (q.crop_x_px, q.crop_y_px), color="white", fontsize=10, fontweight="bold",
                             xytext=(-10, 14) if left else (10, -18), textcoords="offset points", ha="right" if left else "left",
                             arrowprops={"arrowstyle": "-", "color": "white", "lw": 0.8, "alpha": 0.7},
                             bbox={"boxstyle": "round,pad=0.15", "fc": "black", "ec": "none", "alpha": 0.55})
    if truth is not None:
        for q in truth[(truth.cell_id == cid) & (truth.t == t_now)].itertuples():
            axes[1].plot(q.x_um / px_um - fr_now.crop_x0, q.y_um / px_um - fr_now.crop_y0, "o", mfc="none",
                         mec="cyan", ms=14, mew=1.5)
            if hasattr(q, "x_drawn_um") and (abs(q.x_drawn_um - q.x_um) + abs(q.y_drawn_um - q.y_um)) > 1e-9:
                axes[1].plot(q.x_drawn_um / px_um - fr_now.crop_x0, q.y_drawn_um / px_um - fr_now.crop_y0, "x", color="red", ms=9, mew=1.5)
        if strays is not None:
            for q in strays[(strays.cell_id == cid) & (strays.t == t_now) & strays.channel.str.startswith("locus")].itertuples():
                axes[1].plot(q.x_drawn_um / px_um - fr_now.crop_x0, q.y_drawn_um / px_um - fr_now.crop_y0, "s", mfc="none",
                             mec="#f0abfc", ms=11, mew=1.3)
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    st.caption("White + = spot found by the analysis, labelled P (promoter) / E (enhancer) + copy number"
               + ("; cyan ○ = true locus position (simulation, drift-corrected); red × = where that colour's spot is really drawn "
                  "(true position + chromatic shift: the gap between ○ and × is an error no analysis can remove without a calibration); "
                  "pink □ = a stray spot (a decoy the finder must ignore)." if truth is not None
                  else ". (True positions are only available for true cell IDs.)"))

    if len(pos):
        pc = pos[pos.cell_id == cid]
        cen = frames.set_index("t")
        tc = truth[truth.cell_id == cid] if truth is not None else None
        dd = dist[dist.cell_id == cid] if len(dist) else dist
        cps = sorted(set(dd["copy"]) if len(dd) else []) or [0]
        # top row: position errors (wide) + offsets from the truth; bottom row: the promoter-enhancer distance, one panel per copy
        fig = plt.figure(figsize=(17, 8.6))
        gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.85])
        ax = {"path": fig.add_subplot(gs[0, 2])}
        if errs is not None:
            ax["err"] = fig.add_subplot(gs[0, :2])
        else:
            ax["path"].remove()
            ax["path"] = fig.add_subplot(gs[0, :])
        sub = gs[1, :].subgridspec(1, len(cps))
        dist_axes = {c: fig.add_subplot(sub[0, i]) for i, c in enumerate(cps)}
        if errs is not None:                       # how far each found spot is from its true locus, frame by frame
            e = errs[(errs.cell_id == cid) & errs.visible.astype(bool)]
            a = ax["err"]
            for (k, c), g in e.groupby(["locus_id", "copy"]):
                g = g.sort_values("t")
                a.plot(g.t, 1000 * g.err_um, color=COLORS[int(k) % 3], ls="-" if int(c) == 0 else "--", lw=1.2, label=tag(k, c))
                if "on_stray" in g:
                    s_ = g[g.on_stray.astype(bool)]
                    a.plot(s_.t, 1000 * s_.err_um, "s", mfc="none", mec="#f0abfc", ms=7)
                if "on_other_copy" in g:
                    o_ = g[g.on_other_copy.astype(bool)]
                    a.plot(o_.t, 1000 * o_.err_um, "D", mfc="none", mec=PALETTE_WARN, ms=7)
            a.axhline(500, color="#94a3b8", lw=0.8, ls=":")
            a.text(0.01, 500, " 0.5 µm: above = a wrong object was picked", color="#94a3b8", fontsize=8, va="bottom", transform=a.get_yaxis_transform())
            a.set_yscale("log")
            a.set_ylim(max(1.0, float(1000 * e.err_um.min()) * 0.7) if len(e) else 1.0, max(3000.0, float(1000 * e.err_um.max()) * 1.5) if len(e) else 3000.0)
            a.set_ylabel("distance from the true locus (nm)")
            a.set_title("Position error of every found spot")
            a.legend(fontsize=8, ncol=4, loc="lower right")
        for c, a in dist_axes.items():               # the experiment's quantity: the promoter-enhancer distance of every copy
            col = COPY_COLORS[int(c) % 4] if c < 50 else "#f87171"
            if tc is not None and c in set(tc["copy"]):  # truth first, thick and faint, so the found curve on top shows where they differ
                pv = tc[tc["copy"] == c].pivot_table(index="t", columns="locus_id", values=["y_um", "x_um"])
                if ("y_um", 1) in pv.columns:
                    a.plot(pv.index, np.hypot(pv[("y_um", 0)] - pv[("y_um", 1)], pv[("x_um", 0)] - pv[("x_um", 1)]), color=col, lw=4, alpha=0.3,
                           label="true")
            g = dd[dd["copy"] == c].sort_values("t") if len(dd) else dd
            if len(g):
                a.plot(g.t, g.dist_um, color=col, lw=1.1, label="found")
            a.set_ylabel("µm")
            a.set_xlabel("frame")
            a.set_title(f"Promoter-enhancer distance, copy {c if c < 50 else '? (matches no true copy)'}")
            a.legend(fontsize=8, loc="upper right")
            a.axvline(t_now, color="#f59e0b", lw=1)
        a = ax["path"]
        if tc is not None:                          # every found spot minus its true position: the cloud's width = precision, its centre = bias
            j = pc.merge(tc[["t", "locus_id", "copy", "y_um", "x_um"]], on=["t", "locus_id", "copy"])
            for (k, c), g in j.groupby(["locus_id", "copy"]):
                dx, dy = 1000 * (g.lab_x_px * px_um - g.x_um), 1000 * (g.lab_y_px * px_um - g.y_um)
                ok = np.hypot(dx, dy) <= 500
                a.plot(dx[ok], dy[ok], "o" if int(c) == 0 else "^", ms=3, alpha=0.45, color=COLORS[int(k) % 3], label=f"{tag(k, c)} ({(~ok).sum()} beyond 500 nm)")
                if ok.any():
                    a.plot(dx[ok].mean(), dy[ok].mean(), "+", color=COLORS[int(k) % 3], ms=16, mew=2.5)
            a.axhline(0, color="#94a3b8", lw=0.6); a.axvline(0, color="#94a3b8", lw=0.6)
            a.set_aspect("equal", adjustable="datalim")
            a.set_xlabel("found − true, x (nm)")
            a.set_ylabel("found − true, y (nm)")
            a.set_title("Offsets from the truth (+ = mean: a colour shift shows here)")
            a.legend(fontsize=8, loc="upper right")
        else:                                       # no truth (tracked IDs): the paths in a short window around the frame shown
            for (k, c), g in pc.groupby(["locus_id", "copy"]):
                g = g[(g.t - t_now).abs() <= 15].sort_values("t")
                g = g[g.t.isin(cen.index)]
                x = (g.lab_x_px.to_numpy() - cen.loc[g.t, "cx_px"].to_numpy()) * px_um
                y = (g.lab_y_px.to_numpy() - cen.loc[g.t, "cy_px"].to_numpy()) * px_um
                a.plot(x, y, color=COLORS[int(k) % 3], ls="-" if int(c) == 0 else "--", lw=1.0, marker=".", ms=4, label=tag(k, c))
            a.set_aspect("equal", adjustable="datalim")
            a.invert_yaxis()
            a.set_xlabel("x from the nucleus centre (µm)"); a.set_ylabel("y (µm)")
            a.set_title("Paths inside the nucleus (±15 frames around the frame shown)")
            a.legend(fontsize=8)
        if "err" in ax:
            ax["err"].axvline(t_now, color="#f59e0b", lw=1)
            ax["err"].set_xlabel("frame")
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)
        st.caption("Solid = copy 0, dashed = copy 1; colours as in the image (red promoter, green enhancer). "
                   + ("Pink □ = the spot found was a stray spot; orange ◇ = it was the other copy's locus. " if errs is not None else "")
                   + "Orange line = the frame shown above. Copies are numbered like the true ones for true cell IDs (matched per cell).")
    sc_f = run_path / "stage6_analysis" / "scores.json"
    if sc_f.exists():
        with st.expander("Analysis quality and memory use for this run"):
            st.json(json.loads(sc_f.read_text()))
