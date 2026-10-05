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


@st.cache_resource(show_spinner=False)
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
            "parent_of": pd.read_csv(s4 / f"tracks_{v}.csv").set_index("label").parent.to_dict(),
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


def error_frames(R, variant) -> list[int]:
    """Frames where something went wrong: an identity switch or a true cell the tracker's masks missed."""
    V = R["v"][variant]
    m = V["matches"]
    return sorted(set(V["events"].t) | set(m[m.pred_id == 0].t))


def movie_figure(R, variant, t, show_image=False):
    V = R["v"][variant]
    img, gt, trk, m = R["img"][t], R["gt"][t], V["trk"][t], V["matches"]
    mt = m[m.t == t]
    px = R["cfg"]["optics"]["pixel_size_nm"] / 1000.0
    n = 3 if show_image else 2
    fig, axes = plt.subplots(1, n, figsize=(7.5 * n, 7.6))
    axes = np.atleast_1d(axes)
    titles = (["Image (what the tracker sees)"] if show_image else []) + [
        "TRUE cells (ground truth IDs)", "TRACKED cells (tracker's IDs)"]
    for ax, title in zip(axes, titles):
        ax.imshow(img, cmap="gray", vmin=np.percentile(R["img"], 1), vmax=np.percentile(R["img"], 99.8))
        ax.set_title(title, fontsize=15); ax.axis("off")
    ax_gt, ax_tr = axes[-2], axes[-1]
    gmax = int(R["gt"].max())
    draw_labels(ax_gt, gt, lut(gmax, {i: color_for(i) for i in range(1, gmax + 1)}))
    lost = mt[mt.pred_id == 0].gt_id.to_numpy()
    outline(ax_gt, np.isin(gt, lost) & (gt > 0), (1.0, 0.9, 0.0))  # yellow: tracker's masks missed this cell
    pmax = int(max(trk.max(), 1))
    dom = V["dominant"]
    cols = {p: (color_for(dom[p]) if p in dom else np.array([0.55, 0.55, 0.55])) for p in range(1, pmax + 1)}
    draw_labels(ax_tr, trk, lut(pmax, cols))
    cur = dict(zip(mt.pred_id, mt.gt_id))
    nogt = [p for p in np.unique(trk) if p > 0 and p not in cur]
    outline(ax_tr, np.isin(trk, nogt) & (trk > 0), (0.8, 0.2, 1.0))  # purple: mask with no true cell
    # identity switches happening at THIS frame: big red circle + label on both panels
    ev = V["events"]
    ev_t = ev[ev.t == t]
    for r in ev_t.itertuples():
        row = R["cells"][(R["cells"].cell_id == r.gt_id) & (R["cells"].t == t)]
        if row.empty:
            continue
        cx, cy, rad = row.x_um.iloc[0] / px, row.y_um.iloc[0] / px, row.bound_radius_um.iloc[0] / px
        for ax in (ax_gt, ax_tr):
            ax.add_patch(plt.Circle((cx, cy), rad * 1.7, fill=False, ec="red", lw=3.5))
        ax_tr.annotate(f"track {r.old_pred_id} → {r.new_pred_id}" + (" (called a division)" if V["parent_of"].get(r.new_pred_id, 0) > 0
                                                                   else ""), (cx, cy - rad * 1.7), color="white",
                       fontsize=11, weight="bold", ha="center", va="bottom",
                       bbox=dict(boxstyle="round,pad=0.25", fc="red", ec="none"))
    fig.tight_layout()
    return fig, len(lost), ev_t, len(nogt)


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
        speed = st.slider("Speed scale (1 = default, 0.5 = half as fast)", 0.0, 2.0, float(base["motion"]["speed_scale"]), 0.05)
        D = st.number_input("Diffusion D (µm²/s)", 0.0, 1.0, float(base["motion"]["D_um2_s"]), 0.005, format="%.4f")
        drift = st.number_input("Directed speed (µm/s)", 0.0, 0.1, float(base["motion"]["drift_um_s"]), 0.002, format="%.4f")
        p_div = st.slider("Division probability per cell per frame", 0.0, 0.1, float(base["cells"]["p_divide_per_frame"]), 0.005)
        st.markdown("**Nuclear shape**")
        aspect = st.slider("Mean aspect ratio (1 = round)", 1.0, 2.5, float(base["geometry"]["shape"]["aspect_mean"]), 0.05)
        deform = st.slider("Deformation (bumps/dents)", 0.0, 0.25, float(base["geometry"]["shape"]["deform_amp"]), 0.01)
        persist = st.slider("Shape persistence", 0.0, 0.99, float(base["geometry"]["shape"]["persistence"]), 0.01)
        st.markdown("**Nuclear texture** (chromatin-like structure, fixed to each nucleus)")
        tex_c = st.slider("Texture contrast (0 = uniform blob)", 0.0, 0.8, float(base["nucleus_texture"]["contrast"]), 0.05)
        n_nuc = st.slider("Dark nucleoli per nucleus (average)", 0.0, 5.0, float(base["nucleus_texture"]["n_nucleoli"]), 0.5)
        st.markdown("**Imaging**")
        photons = st.number_input("Photons/pixel/s (brightness)", 10, 100000, int(base["optics"]["photons_per_px_s"]), 100)
        read_noise = st.number_input("Camera read noise (e-)", 0.0, 20.0, float(base["optics"]["read_noise_e"]), 0.5)
        pix = st.number_input("Pixel size (nm)", 50, 2000, int(base["optics"]["pixel_size_nm"]), 10)
        NA = st.number_input("NA", 0.3, 1.7, float(base["optics"]["NA"]), 0.05)
        st.markdown("**Segmentation and tracker**")
        seg_method = st.selectbox("Segmentation", ["cellpose", "threshold_watershed"],
                                  0 if base["segmentation"]["method"] == "cellpose" else 1)
        mode = st.selectbox("Trackastra mode", ["greedy", "greedy_nodiv", "ilp"], 0)
        go = st.form_submit_button("▶ Run simulation", type="primary", width="stretch")

with st.sidebar.expander("🔬 New single-cell loci run"):
    st.caption("Uses configs/loci_demo.yaml: fine pixels, a few cells, two coloured loci per nucleus (stand-in "
               "dynamics), then isolates every cell and locates the loci in each one in parallel (about 2 minutes).")
    with st.form("loci_sim"):
        l_seed = st.number_input("Random seed", 0, 10_000, 1, key="l_seed")
        l_cells = st.slider("Cells in the field", 2, 20, 8, key="l_cells")
        l_frames = st.slider("Frames", 20, 300, 120, key="l_frames")
        l_dt = st.number_input("Frame interval (s)", 1, 600, 10, key="l_dt")
        l_speed = st.slider("Cell speed scale", 0.0, 2.0, 0.4, 0.05, key="l_speed")
        l_phot = st.number_input("Photons per locus per second", 200, 100000, 6000, 200, key="l_phot")
        l_div = st.slider("Division probability per cell per frame", 0.0, 0.05, 0.0, 0.005, key="l_div")
        go_loci = st.form_submit_button("▶ Run single-cell loci simulation", type="primary", width="stretch")


def launch(overrides, config=None):
    """Run scripts/run_tracking_demo.py in a separate process and stream its progress into the page."""
    cmd = [sys.executable, str(REPO_ROOT / "scripts" / "run_tracking_demo.py")]
    if config:
        cmd += ["--config", str(config)]
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


if go:
    launch([
        f"seed={seed}", f"acquisition.n_frames={n_frames}", f"acquisition.frame_interval_s={dt}",
        f"cells.n_cells={n_cells}", f"cells.p_divide_per_frame={p_div}", f"geometry.fov_um=[{fov}, {fov}]",
        f"geometry.cell_radius_um={radius}", f"motion.D_um2_s={D}", f"motion.drift_um_s={drift}",
        f"geometry.shape.aspect_mean={aspect}", f"geometry.shape.deform_amp={deform}",
        f"geometry.shape.persistence={persist}", f"optics.photons_per_px_s={photons}",
        f"optics.read_noise_e={read_noise}", f"optics.pixel_size_nm={pix}", f"optics.NA={NA}",
        f"tracking.mode={mode}", f"nucleus_texture.contrast={tex_c}", f"nucleus_texture.n_nucleoli={n_nuc}",
        f"segmentation.method={seg_method}", f"motion.speed_scale={speed}"])
if go_loci:
    launch([f"seed={l_seed}", f"cells.n_cells={l_cells}", f"acquisition.n_frames={l_frames}",
            f"acquisition.frame_interval_s={l_dt}", f"motion.speed_scale={l_speed}",
            f"loci.photons_per_locus_s={l_phot}", f"cells.p_divide_per_frame={l_div}"],
           config=REPO_ROOT / "configs" / "loci_demo.yaml")

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

tab_movie, tab_metrics, tab_traj, tab_scan, tab_cells, tab_cfg = st.tabs(
    ["🎞 Movie", "📊 Metrics", "🧭 Trajectories", "📈 Scans", "🔬 Cells", "⚙ Config"])

# ----------------------------------------------------------------------------- movie
@st.fragment
def movie_tab():
    variant = st.radio("Tracked on:", list(VARIANT_LABEL), format_func=VARIANT_LABEL.get, horizontal=True)
    if "frame" not in st.session_state:
        st.session_state["frame"] = 0
    st.session_state["frame"] = min(st.session_state["frame"], T - 1)
    errs = error_frames(R, variant)

    def go_prev_err():
        prev = [f for f in errs if f < st.session_state["frame"]]
        if prev:
            st.session_state["frame"] = prev[-1]

    def go_next_err():
        nxt = [f for f in errs if f > st.session_state["frame"]]
        if nxt:
            st.session_state["frame"] = nxt[0]

    b0, b1, b2, b3, b4 = st.columns([1, 1, 2, 2, 3])
    b0.button("◀", help="Previous frame", on_click=lambda: st.session_state.update(frame=max(0, st.session_state["frame"] - 1)))
    b1.button("▶", help="Next frame", on_click=lambda: st.session_state.update(frame=min(T - 1, st.session_state["frame"] + 1)))
    b2.button("⏮ Previous error", on_click=go_prev_err, disabled=not errs)
    b3.button("Next error ⏭", on_click=go_next_err, disabled=not errs, type="primary")
    show_img = b4.checkbox("Also show raw image panel", value=False)
    t = st.slider("Frame", 0, T - 1, key="frame")
    st.caption(f"{len(errs)} of {T} frames contain an error: "
               + (", ".join(map(str, errs[:40])) + (" …" if len(errs) > 40 else "") if errs else "none"))

    fig, n_lost, ev_t, n_nogt = movie_figure(R, variant, t, show_image=show_img)
    st.pyplot(fig, width="stretch")
    plt.close(fig)

    if len(ev_t) or n_lost or n_nogt:
        st.error(f"**Errors in frame {t}:** {len(ev_t)} identity switch(es), {n_lost} true cell(s) missed by the masks, "
                 f"{n_nogt} false-positive mask(s).")
        for r in ev_t.itertuples():
            par = R["v"][variant]["parent_of"].get(r.new_pred_id, 0)
            why = (f" The tracker treated this as a **division** of track {par} (a false division unless the cell really divided)."
                   if par > 0 else " The tracker lost the link and started a new track (or swapped it with a neighbour).")
            st.write(f"🔴 True cell **{r.gt_id}** was track **{r.old_pred_id}** in the previous frame and is now "
                     f"track **{r.new_pred_id}**.{why}")
    else:
        st.success(f"Frame {t}: no errors. Every true cell kept its track ID.")
    with st.expander("How to read the pictures"):
        st.markdown(
            "- **Left (TRUE cells):** the simulation's ground truth. Each true cell has its own color and ID number.\n"
            "- **Right (TRACKED cells):** what the tracker produced, with the tracker's own track numbers. A track is "
            "colored by the true cell it mostly follows, so *matching colors on both sides = tracked correctly*.\n"
            "- **Red circle + red label = an identity switch happening in this frame:** the cell had one track ID in "
            "the previous frame and a different one now (the tracker broke the track or swapped two cells).\n"
            "- **Yellow outline (left):** the segmentation missed this true cell. **Purple outline / grey (right):** a "
            "mask with no true cell under it (false positive).\n"
            "- Use **Next error ⏭** to jump to frames where something went wrong.")
    ev = R["v"][variant]["events"]
    st.markdown("**Identity switches over time** (each bar = switches at that frame)")
    st.bar_chart(ev.groupby("t").size().reindex(range(T), fill_value=0))

# ----------------------------------------------------------------------------- metrics
@st.fragment
def metrics_tab():
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
@st.fragment
def traj_tab():
    v = st.radio("Tracked on:", list(VARIANT_LABEL), format_func=VARIANT_LABEL.get, horizontal=True, key="trajv")
    fig = trajectory_figure(R, v)
    st.pyplot(fig)
    plt.close(fig)

with tab_cfg:
    st.code(yaml.safe_dump(R["cfg"], sort_keys=False), language="yaml")
    st.caption(f"Files for this run: {run_path}")

# ----------------------------------------------------------------------------- scans
SCAN_METRICS = {"CHOTA": "CHOTA (higher is better)", "LNK": "LNK (higher is better)",
                "target_eff": "Target effectiveness (higher is better)", "id_switches": "ID switches (lower is better)",
                "identity_kept": "Identity kept (higher is better)", "DET": "DET (higher is better)",
                "TRA": "TRA (higher is better)"}
@st.fragment
def scan_tab():
    scans = sorted((RUNS.parent / "sweeps").glob("*/results.csv"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not scans:
        st.info("No scans yet. Run one from a terminal, e.g.\n\n"
                "`scripts\\run.cmd python scripts\\sweep.py --name speed --seeds 1-5 --grid motion.speed_scale=1,0.5,0.25`")
    else:
        sname = st.selectbox("Scan", [p.parent.name for p in scans])
        sdf = pd.read_csv(RUNS.parent / "sweeps" / sname / "results.csv")
        params = [c for c in sdf.columns if c not in {"seed", "masks", "run", *SCAN_METRICS}]
        n_seeds = sdf.seed.nunique()
        st.caption(f"{len(sdf) // 2} runs ({n_seeds} seeds per setting). Lines = mean over seeds, bars = ± 1 standard "
                   "deviation across seeds, dots = individual seeds. Differences smaller than the bars are not reliable.")
        c1, c2, c3 = st.columns(3)
        metric = c1.selectbox("Metric", list(SCAN_METRICS), format_func=SCAN_METRICS.get)
        masks_sel = c3.radio("Masks", ["auto_masks", "gt_masks", "both"], horizontal=True,
                             format_func={"auto_masks": "Automatic segmentation", "gt_masks": "Ground-truth masks",
                                          "both": "Both"}.get)
        if not params:
            st.warning("This scan has no varied settings to plot.")
        else:
            xcol = c2.selectbox("X axis", params)
            others = [p for p in params if p != xcol]
            d = sdf if masks_sel == "both" else sdf[sdf.masks == masks_sel]
            group_cols = others + (["masks"] if masks_sel == "both" else [])
            xs = sorted(d[xcol].unique(), key=lambda v: (isinstance(v, str), v))
            xpos = {v: i for i, v in enumerate(xs)}
            numeric = all(not isinstance(v, str) for v in xs)
            fig, ax = plt.subplots(figsize=(9, 5))
            groups = list(d.groupby(group_cols)) if group_cols else [((), d)]
            for gi, (key, g) in enumerate(groups):
                key = key if isinstance(key, tuple) else (key,)
                label = ", ".join(f"{c}={v}" for c, v in zip(group_cols, key)) or "all runs"
                col = plt.cm.tab10(gi % 10)
                m = g.groupby(xcol)[metric].agg(["mean", "std"]).reindex(xs)
                xv = [v if numeric else xpos[v] for v in xs]
                ax.errorbar(xv, m["mean"], yerr=m["std"].fillna(0), marker="o", capsize=4, lw=2, color=col, label=label)
                ax.scatter([v if numeric else xpos[v] for v in g[xcol]], g[metric], s=14, alpha=0.35, color=col)
            if not numeric:
                ax.set_xticks(range(len(xs))); ax.set_xticklabels([str(v) for v in xs])
            ax.set_xlabel(xcol); ax.set_ylabel(SCAN_METRICS[metric])
            ax.grid(alpha=0.3)
            if len(groups) > 1:
                ax.legend(fontsize=9)
            st.pyplot(fig, width="stretch")
            plt.close(fig)
            summary = d.groupby(group_cols + [xcol])[list(SCAN_METRICS)].agg(["mean", "std"]).round(3)
            with st.expander("Table (mean and standard deviation over seeds)"):
                st.dataframe(summary, width="stretch")

# ----------------------------------------------------------------------------- single cells
import cells_tab  # noqa: E402  (app/cells_tab.py)
from fastview import movie_frames  # noqa: E402
from player import show_player  # noqa: E402

@st.fragment
def cells_frag():
    cells_tab.render(st, run_path, R["cfg"])


@st.cache_resource(show_spinner=False, max_entries=4)
def cached_movie(path: str, stamp: float, variant: str, show_image: bool):
    """Pre-render every frame once (numpy/PIL, fast) for the in-browser player."""
    return movie_frames(load_run(path, stamp), variant, show_image)


@st.fragment
def movie_player():
    with st.expander("▶ Play the movie (smooth, runs in your browser)", expanded=st.session_state.get("mv_play", False)):
        st.caption("Pre-renders every frame once, then plays in your browser with no waiting between frames. "
                   "Space = play/pause, arrow keys = step. Red circles mark identity switches.")
        c1, c2, c3 = st.columns([2, 2, 3])
        on = c1.toggle("Prepare player", key="mv_play", help="Takes a few seconds the first time for each run")
        pv = c2.radio("Tracked on", list(VARIANT_LABEL), format_func=lambda v: "Automatic masks" if v == "auto_masks" else "Ground-truth masks",
                      key="mv_play_var", horizontal=True)
        raw = c3.checkbox("Also show the raw image", False, key="mv_play_raw")
        if on:
            with st.spinner("Rendering frames…"):
                frames, labels = cached_movie(str(run_path), (run_path / "stage7_validation" / "metrics.json").stat().st_mtime,
                                              pv, raw)
            show_player(frames, labels, fps=8, key="movie")


with tab_movie:
    movie_player()
    movie_tab()
with tab_metrics:
    metrics_tab()
with tab_traj:
    traj_tab()
with tab_scan:
    scan_tab()
with tab_cells:
    cells_frag()
