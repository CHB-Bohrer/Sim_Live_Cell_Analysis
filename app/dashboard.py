"""Dashboard: set parameters, run the cell-tracking simulation, and see what happened.

    scripts\\dashboard.cmd          (opens http://localhost:8501)
"""
import colorsys
import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import tifffile
import yaml
from scipy import ndimage as ndi
from skimage.segmentation import find_boundaries

from simlive.io.runs import REPO_ROOT, load_config

# NOTE: the simulation (PyTorch/OpenMM/GPU) runs in a separate process (scripts/run_tracking_demo.py). Importing
# torch into this web server clashes with the plotting libraries' OpenMP runtime on Windows and crashes it.

RUNS = REPO_ROOT / "data" / "runs"
DEFAULT_CFG = REPO_ROOT / "configs" / "tracking_demo.yaml"
VARIANT_LABEL = {"auto_masks": "Automatic segmentation (realistic)", "gt_masks": "Ground-truth masks (linking error only)"}

st.set_page_config(page_title="Cell tracking validation", layout="wide")


# ----------------------------------------------------------------------------- helpers
def color_for(i: int) -> np.ndarray:
    """Stable color per ID (golden-ratio hue spacing), so one ID is always the same color."""
    return np.array(colorsys.hsv_to_rgb((i * 0.61803398875) % 1.0, 0.65, 1.0))


def lut(max_id: int, mapping: dict[int, np.ndarray]) -> np.ndarray:
    t = np.zeros((max_id + 1, 3))
    for k, v in mapping.items():
        t[k] = v
    return t


def list_runs() -> list[Path]:
    rs = [p for p in RUNS.glob("*") if (p / "stage7_validation" / "id_matches_auto_masks.csv").exists()]
    return sorted(rs, key=lambda p: (p / "stage7_validation" / "metrics.json").stat().st_mtime, reverse=True)


@st.cache_data(show_spinner=False)
def load_run(path: str, stamp: float):
    r = Path(path)
    s3, s4, s7 = r / "stage3_microscopy", r / "stage4_segtrack", r / "stage7_validation"
    out = {
        "img": tifffile.imread(s3 / "nucleus.tif"), "gt": tifffile.imread(s3 / "labels.tif"),
        "cells": pd.read_csv(r / "stage2_cells" / "cells.csv"),
        "cfg": yaml.safe_load((r / "stage2_cells" / "params.yaml").read_text()),
        "metrics": json.loads((s7 / "metrics.json").read_text()), "v": {},
    }
    for v in ("gt_masks", "auto_masks"):
        matches = pd.read_csv(s7 / f"id_matches_{v}.csv")
        dom = (matches[matches.pred_id > 0].groupby("pred_id").gt_id.agg(lambda s: s.mode().iloc[0])).to_dict()
        out["v"][v] = {
            "trk": tifffile.imread(s4 / f"tracked_{v}.tif"), "matches": matches, "dominant": dom,
            "events": pd.read_csv(s7 / f"id_switch_events_{v}.csv"), "per_cell": pd.read_csv(s7 / f"id_per_cell_{v}.csv"),
        }
    out["v"]["auto_masks"]["seg"] = tifffile.imread(s4 / "auto_segmentation.tif")
    return out


def draw_labels(ax, lab, rgb_lut, texts=True, alpha=0.5):
    h, w = lab.shape
    rgba = np.zeros((h, w, 4))
    rgba[..., :3] = rgb_lut[lab]
    rgba[..., 3] = np.where(lab > 0, alpha, 0)
    ax.imshow(rgba, interpolation="nearest")
    if texts:
        ids = [i for i in np.unique(lab) if i > 0]
        if ids:
            for i, (cy, cx) in zip(ids, ndi.center_of_mass(lab > 0, lab, ids)):
                ax.text(cx, cy, str(i), color="white", fontsize=8, ha="center", va="center", weight="bold")


def outline(ax, mask, color):
    if mask.any():
        b = find_boundaries(mask, mode="inner")
        o = np.zeros(mask.shape + (4,)); o[b] = (*color, 1.0)
        ax.imshow(o, interpolation="nearest")


def movie_figure(R, variant, t):
    V = R["v"][variant]
    img, gt, trk, m = R["img"][t], R["gt"][t], V["trk"][t], V["matches"]
    mt = m[m.t == t]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6))
    for ax, title in zip(axes, ("Image (what the tracker sees)", "Ground truth (true cell IDs)",
                                 "Tracked (predicted track IDs)")):
        ax.imshow(img, cmap="gray", vmin=np.percentile(R["img"], 1), vmax=np.percentile(R["img"], 99.8))
        ax.set_title(title); ax.axis("off")
    # ground truth: color = true ID; yellow outline = this cell was not found/matched by the tracker's masks
    gmax = int(R["gt"].max())
    draw_labels(axes[1], gt, lut(gmax, {i: color_for(i) for i in range(1, gmax + 1)}))
    lost = mt[mt.pred_id == 0].gt_id.to_numpy()
    outline(axes[1], np.isin(gt, lost) & (gt > 0), (1.0, 0.9, 0.0))
    # tracked: color = the TRUE cell this track mostly follows; red outline = track is on a different cell right now
    pmax = int(max(trk.max(), 1))
    dom = V["dominant"]
    cols = {p: (color_for(dom[p]) if p in dom else np.array([0.55, 0.55, 0.55])) for p in range(1, pmax + 1)}
    draw_labels(axes[2], trk, lut(pmax, cols))
    cur = dict(zip(mt.pred_id, mt.gt_id))  # track -> GT cell it currently overlaps
    wrong = [p for p in np.unique(trk) if p > 0 and p in dom and p in cur and cur[p] != dom[p]]
    nogt = [p for p in np.unique(trk) if p > 0 and p not in cur]
    outline(axes[2], np.isin(trk, wrong) & (trk > 0), (1.0, 0.1, 0.1))
    outline(axes[2], np.isin(trk, nogt) & (trk > 0), (0.8, 0.2, 1.0))
    fig.tight_layout()
    return fig, len(lost), len(wrong), len(nogt)


def trajectory_figure(R, variant):
    cfg = R["cfg"]
    px = cfg["optics"]["pixel_size_nm"] / 1000.0
    V = R["v"][variant]
    fig, ax = plt.subplots(figsize=(8, 8))
    c = R["cells"]
    for cid, g in c.groupby("cell_id"):
        ax.plot(g.x_um, g.y_um, "-", color=color_for(cid), lw=2.2, alpha=0.9)
    trk = V["trk"]
    paths = {}
    for t in range(len(trk)):
        ids = [i for i in np.unique(trk[t]) if i > 0]
        for i, (cy, cx) in zip(ids, ndi.center_of_mass(trk[t] > 0, trk[t], ids)):
            paths.setdefault(i, []).append((cx * px, cy * px))
    for p, pts in paths.items():
        pts = np.array(pts)
        col = color_for(V["dominant"][p]) if p in V["dominant"] else (0.5, 0.5, 0.5)
        ax.plot(pts[:, 0], pts[:, 1], "--", color="k", lw=0.9, alpha=0.8)
    ev = V["events"]
    for r in ev.itertuples():
        row = c[(c.cell_id == r.gt_id) & (c.t == r.t)]
        if len(row):
            ax.plot(row.x_um.iloc[0], row.y_um.iloc[0], "rx", ms=11, mew=2.5)
    ax.set_xlim(0, cfg["geometry"]["fov_um"][1]); ax.set_ylim(cfg["geometry"]["fov_um"][0], 0)
    ax.set_aspect("equal"); ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    ax.set_title("True paths (colored) vs tracked paths (black dashed); red × = identity switch")
    return fig


# ----------------------------------------------------------------------------- sidebar
st.sidebar.title("Cell tracking validation")
base = load_config(DEFAULT_CFG)

with st.sidebar.expander("➕ New simulation", expanded=not list_runs()):
    with st.form("sim"):
        st.caption("Defaults come from configs/tracking_demo.yaml")
        seed = st.number_input("Random seed", 0, 10_000, int(base["seed"]))
        st.markdown("**Acquisition**")
        n_frames = st.slider("Frames", 10, 120, int(base["acquisition"]["n_frames"]))
        dt = st.number_input("Frame interval (s)", 5, 3600, int(base["acquisition"]["frame_interval_s"]))
        st.markdown("**Cells and motion**")
        n_cells = st.slider("Initial cells", 3, 80, int(base["cells"]["n_cells"]))
        fov = st.slider("Field of view (µm, square)", 60, 400, int(base["geometry"]["fov_um"][0]))
        radius = st.slider("Nuclear radius (µm)", 3.0, 12.0, float(base["geometry"]["cell_radius_um"]), 0.5)
        D = st.number_input("Diffusion D (µm²/s)", 0.0, 1.0, float(base["motion"]["D_um2_s"]), 0.005, format="%.4f")
        drift = st.number_input("Directed speed (µm/s)", 0.0, 0.1, float(base["motion"]["drift_um_s"]), 0.002, format="%.4f")
        p_div = st.slider("Division probability per cell per frame", 0.0, 0.1, float(base["cells"]["p_divide_per_frame"]), 0.005)
        st.markdown("**Nuclear shape**")
        aspect = st.slider("Mean aspect ratio (1 = round)", 1.0, 2.5, float(base["geometry"]["shape"]["aspect_mean"]), 0.05)
        deform = st.slider("Deformation (bumps/dents)", 0.0, 0.25, float(base["geometry"]["shape"]["deform_amp"]), 0.01)
        persist = st.slider("Shape persistence", 0.0, 0.99, float(base["geometry"]["shape"]["persistence"]), 0.01)
        st.markdown("**Imaging**")
        photons = st.number_input("Photons/pixel/s (brightness)", 10, 100000, int(base["optics"]["photons_per_px_s"]), 100)
        read_noise = st.number_input("Camera read noise (e-)", 0.0, 20.0, float(base["optics"]["read_noise_e"]), 0.5)
        pix = st.number_input("Pixel size (nm)", 50, 2000, int(base["optics"]["pixel_size_nm"]), 10)
        NA = st.number_input("NA", 0.3, 1.7, float(base["optics"]["NA"]), 0.05)
        st.markdown("**Tracker**")
        mode = st.selectbox("Trackastra mode", ["greedy", "greedy_nodiv", "ilp"], 0)
        go = st.form_submit_button("▶ Run simulation", type="primary", width="stretch")

if go:
    overrides = [
        f"seed={seed}", f"acquisition.n_frames={n_frames}", f"acquisition.frame_interval_s={dt}",
        f"cells.n_cells={n_cells}", f"cells.p_divide_per_frame={p_div}", f"geometry.fov_um=[{fov}, {fov}]",
        f"geometry.cell_radius_um={radius}", f"motion.D_um2_s={D}", f"motion.drift_um_s={drift}",
        f"geometry.shape.aspect_mean={aspect}", f"geometry.shape.deform_amp={deform}",
        f"geometry.shape.persistence={persist}", f"optics.photons_per_px_s={photons}",
        f"optics.read_noise_e={read_noise}", f"optics.pixel_size_nm={pix}", f"optics.NA={NA}",
        f"tracking.mode={mode}"]
    cmd = [sys.executable, str(REPO_ROOT / "scripts" / "run_tracking_demo.py")]
    for o in overrides:
        cmd += ["--set", o]
    with st.status("Running simulation…", expanded=True) as status:
        proc = subprocess.Popen(cmd, cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                env={**os.environ, "PYTHONUNBUFFERED": "1"})
        log, run_name = [], None
        for line in proc.stdout:
            line = line.rstrip()
            log.append(line)
            if line.startswith("Run folder:"):
                run_name = Path(line.split(":", 1)[1].strip()).name
            if line.startswith(("Stage", "Run folder", "Done")):
                st.write(line)
        if proc.wait() == 0 and run_name:
            status.update(label="Done", state="complete", expanded=False)
            st.session_state["run"] = run_name
            st.session_state["frame"] = 0
        else:
            status.update(label="Failed", state="error")
            st.error("The simulation process failed. Last output:")
            st.code("\n".join(log[-40:]))

runs = list_runs()
if not runs:
    st.info("No runs yet. Open **New simulation** in the sidebar and press **Run simulation**.")
    st.stop()
names = [r.name for r in runs]
if st.session_state.get("run") not in names:
    st.session_state["run"] = names[0]
st.sidebar.selectbox("Run to view", names, key="run")
run_path = RUNS / st.session_state["run"]
R = load_run(str(run_path), (run_path / "stage7_validation" / "metrics.json").stat().st_mtime)
T = len(R["gt"])

# ----------------------------------------------------------------------------- headline numbers
st.title("Cell tracking: how well did we do?")
st.caption(f"Run `{run_path.name}` — {R['cells'].cell_id.nunique()} cell IDs, {T} frames, "
           f"{R['cfg']['optics']['pixel_size_nm']} nm/px")
for v in ("auto_masks", "gt_masks"):
    m = R["metrics"][v]
    st.subheader(VARIANT_LABEL[v])
    a, b, c, d = st.columns(4)
    a.metric("CHOTA", f"{m['CHOTAMetric']['CHOTA']:.3f}", help="Detection + trajectory association, 1 = perfect")
    b.metric("LNK", f"{m['CTCMetrics']['LNK']:.3f}", help="Cell Tracking Challenge linking score, 1 = perfect")
    c.metric("ID switches", m["Identity"]["id_switches"], help="Times a true cell changed track ID")
    d.metric("Identity kept", f"{100 * m['Identity']['identity_preserved_fraction']:.1f}%",
             help="Fraction of cell-frames on the cell's main track ID")

tab_movie, tab_metrics, tab_traj, tab_cfg = st.tabs(["🎞 Movie", "📊 Metrics", "🧭 Trajectories", "⚙ Config"])

# ----------------------------------------------------------------------------- movie
with tab_movie:
    variant = st.radio("Tracked on:", list(VARIANT_LABEL), format_func=VARIANT_LABEL.get, horizontal=True)
    if "frame" not in st.session_state:
        st.session_state["frame"] = 0
    st.session_state["frame"] = min(st.session_state["frame"], T - 1)
    c1, c2, c3 = st.columns([1, 1, 10])
    c1.button("◀", on_click=lambda: st.session_state.update(frame=max(0, st.session_state["frame"] - 1)))
    c2.button("▶", on_click=lambda: st.session_state.update(frame=min(T - 1, st.session_state["frame"] + 1)))
    t = c3.slider("Frame", 0, T - 1, key="frame")
    fig, n_lost, n_wrong, n_nogt = movie_figure(R, variant, t)
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    st.caption(
        "**Middle:** each true cell has a fixed color; **yellow outline** = the tracker's masks missed it. "
        "**Right:** each track is colored by the true cell it mostly follows, so *a track whose color differs from "
        "its cell in the middle panel has swapped*. **Red outline** = this track is on a different cell than its "
        "usual one right now (identity swap). A track that merely *breaks* (the cell gets a new track number) keeps its "
        "color, so watch the numbers too. **Purple outline / grey** = a mask with no true cell under it (false positive).")
    st.write(f"This frame: **{n_lost}** true cells missed, **{n_wrong}** tracks on the wrong cell, **{n_nogt}** false-positive masks.")
    ev = R["v"][variant]["events"]
    st.markdown("**Identity switches over time** (each bar = switches at that frame)")
    st.bar_chart(ev.groupby("t").size().reindex(range(T), fill_value=0))
    if (ev.t == t).any():
        st.warning("Switches at this frame (true cell → old track ID → new track ID):")
        st.dataframe(ev[ev.t == t], hide_index=True)

# ----------------------------------------------------------------------------- metrics
with tab_metrics:
    rows = {}
    for v in ("auto_masks", "gt_masks"):
        flat = {}
        for grp, vals in R["metrics"][v].items():
            for k, x in vals.items():
                flat[f"{grp}.{k}"] = x
        rows[VARIANT_LABEL[v]] = flat
    st.dataframe(pd.DataFrame(rows), width="stretch")
    with st.expander("What do these numbers mean?"):
        st.markdown("""
**Matching.** A true cell counts as found if one predicted mask covers >50% of it. *FN nodes* = missed true cells,
*FP nodes* = predicted cells with no true cell, *NS nodes* = one mask covering several true cells.
**Links** connect a cell to itself in the next frame (or a parent to a daughter). *FN edges* = true links the tracker
missed, *FP edges* = links it invented (a swap or wrong merge), *WS edges* = wrong link type (division vs continuation).

**AOGM** = 5·NS + 10·FN nodes + 1·FP nodes + 1.5·FN edges + 1·FP edges + 1·WS edges (Cell Tracking Challenge weights).
**TRA / DET / LNK** = 1 − AOGM-cost ÷ (cost of building the truth from scratch), using all errors / node errors / link
errors only. They saturate near 1, so look at the others too. **CHOTA** blends detection accuracy and
trajectory-association accuracy. **Track purity** = how much of each predicted track follows one true cell;
**target effectiveness** = how much of each true cell's life is covered by one track.

**Identity** (computed here): *id_switches* = times a true cell's matched track ID changed between consecutive frames;
*identity_preserved_fraction* = share of cell-frames carrying that cell's most common track ID.""")
    v = st.selectbox("Per-cell identity table for:", list(VARIANT_LABEL), format_func=VARIANT_LABEL.get)
    st.dataframe(R["v"][v]["per_cell"].sort_values("id_switches", ascending=False), hide_index=True,
                 width="stretch")

# ----------------------------------------------------------------------------- trajectories
with tab_traj:
    v = st.radio("Tracked on:", list(VARIANT_LABEL), format_func=VARIANT_LABEL.get, horizontal=True, key="trajv")
    fig = trajectory_figure(R, v)
    st.pyplot(fig)
    plt.close(fig)

with tab_cfg:
    st.code(yaml.safe_dump(R["cfg"], sort_keys=False), language="yaml")
    st.caption(f"Files for this run: {run_path}")
