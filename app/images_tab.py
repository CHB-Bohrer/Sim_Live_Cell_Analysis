"""The dashboard's 🖼 Images tab: see the recorded images next to the ground truth that produced them.

Four parts: (1) the recorded frame per colour channel with the true masks / true locus positions drawn on top and the numbers that
were behind that frame; (2) the same image built up step by step, one source of image error at a time; (3) the effect of each error
source on its own; (4) a plain-language reference of every source with this run's settings and some derived numbers.

Frames are re-rendered on demand by `Imager` (stage 3) from the run's saved ground truth; the Imager is deterministic per frame, so the
"all errors" version is identical to the saved TIFF (the tab checks and says so). Imports numpy / scipy / pandas only, never torch.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import tifffile
import yaml
from skimage.segmentation import find_boundaries

import style
from simlive.stage3_microscopy import errors as ERR
from simlive.stage3_microscopy.render import imager_for_run, psf_sigma_px

CH_TITLE = {"nucleus": "Nuclear channel (NLS-GFP; also carries the MS2 spots)", "locus0": "Locus 0 channel (promoter)",
            "locus1": "Locus 1 channel (enhancer)"}
CH_SHORT = {"nucleus": "nucleus", "locus0": "locus 0 (promoter)", "locus1": "locus 1 (enhancer)"}


# ----------------------------------------------------------------------------- cached loading
@st.cache_resource(show_spinner="Rebuilding the imaging model of this run…")
def _imager(path: str, stamp: float):
    run = Path(path)
    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    im = imager_for_run(run, cfg)
    im.variants = {}
    return im


def _variant(im, enabled: tuple):
    """ErrorSet with exactly these sources on; cached on the Imager so fixed patterns are built once."""
    if enabled not in im.variants:
        im.variants[enabled] = im.E.with_enabled(enabled)
    return im.variants[enabled]


@st.cache_data(show_spinner=False, max_entries=96)
def _render(path: str, stamp: float, t: int, channel: str, enabled: tuple, bbox: tuple, zsel: int = -1):
    """One 2D image of a frame with exactly the sources in `enabled` on. In a z-stack: slice `zsel`, or the projection if zsel < 0."""
    im = _imager(path, stamp)
    f = im.frame(t, _variant(im, enabled), [channel])["images"][channel]
    if f.ndim == 3:
        f = f[zsel] if zsel >= 0 else im.project(f, im.projection)
    y0, y1, x0, x1 = bbox
    return f[y0:y1, x0:x1].copy()


@st.cache_data(show_spinner=False, max_entries=24)
def _render_vol(path: str, stamp: float, t: int, channel: str, enabled: tuple, bbox: tuple):
    """The whole z-stack (z, y, x) of one channel in a frame, cropped laterally."""
    im = _imager(path, stamp)
    f = im.frame(t, _variant(im, enabled), [channel])["images"][channel]
    y0, y1, x0, x1 = bbox
    return f[:, y0:y1, x0:x1].copy()


@st.cache_data(show_spinner=False)
def _read_page(path: str, stamp: float, name: str, t: int):
    return tifffile.imread(Path(path) / "stage3_microscopy" / f"{name}.tif", key=int(t))


def _recorded(run: Path, stamp: float, im, ch: str, t: int, zsel: int):
    """The saved camera image of a channel: one slice of the z-stack, or the 2D image (projection) the tracker saw."""
    if im.is3d and zsel >= 0:
        return _read_page(str(run), stamp, f"{ch}_zstack", t * im.nz + zsel)
    return _read_page(str(run), stamp, ch, t)


def _labels_at(run: Path, im, t: int, zsel: int):
    if im.is3d and zsel >= 0:
        return _read_page(str(run), 0.0, "labels_3d", t * im.nz + zsel)
    return _read_page(str(run), 0.0, "labels", t)


@st.cache_data(show_spinner=False, max_entries=12)
def _read_stack(path: str, stamp: float, name: str, t: int, nz: int):
    return tifffile.imread(Path(path) / "stage3_microscopy" / f"{name}.tif", key=range(int(t) * nz, (int(t) + 1) * nz))


def _csv(run: Path, *parts) -> pd.DataFrame | None:
    p = run.joinpath(*parts)
    return pd.read_csv(p) if p.exists() else None


# ----------------------------------------------------------------------------- small helpers
def _limits(img: np.ndarray, hi_pct: float) -> tuple[float, float]:
    return float(np.percentile(img, 1)), float(np.percentile(img, hi_pct))


def _crop_box(mode, shape, px_um, cell_row, spot_xy):
    ny, nx = shape
    if mode == "field" or (mode == "cell" and cell_row is None) or (mode == "locus" and spot_xy is None):
        return 0, ny, 0, nx
    if mode == "cell":
        cy, cx = cell_row.y_um / px_um, cell_row.x_um / px_um
        h = int(1.7 * cell_row.bound_radius_um / px_um)
    else:
        (cy, cx), h = spot_xy, 14
    y0, x0 = max(int(cy) - h, 0), max(int(cx) - h, 0)
    return y0, min(int(cy) + h + 1, ny), x0, min(int(cx) + h + 1, nx)


def _centroid(img: np.ndarray, y: float, x: float, half: int = 4):
    """Intensity-weighted centroid of a spot in a small window around (y, x) (local background removed)."""
    y0, x0 = max(int(round(y)) - half, 0), max(int(round(x)) - half, 0)
    w = img[y0:int(round(y)) + half + 1, x0:int(round(x)) + half + 1].astype(float)
    if w.size == 0:
        return np.nan, np.nan
    w = np.clip(w - np.median(img), 0, None)
    if w.sum() <= 0:
        return np.nan, np.nan
    yy, xx = np.indices(w.shape)
    return y0 + (w * yy).sum() / w.sum(), x0 + (w * xx).sum() / w.sum()


def _zcentroid(vol: np.ndarray, y: int, x: int) -> float:
    """Centre of a spot along z (in slice units) from the z-profile through its brightest pixel (background removed)."""
    p = vol[:, max(y - 1, 0):y + 2, max(x - 1, 0):x + 2].mean(axis=(1, 2))
    p = np.clip(p - np.median(p), 0, None)
    return float((p * np.arange(len(p))).sum() / p.sum()) if p.sum() > 0 else float("nan")


def _thompson_nm(sigma_nm: float, a_nm: float, n_photons: float, bg_sd_e: float) -> float:
    """Best-case localization precision (Thompson et al. 2002, with background): sqrt((s^2 + a^2/12)/N + 8 pi s^4 b^2 / (a^2 N^2))."""
    n = max(n_photons, 1e-9)
    return float(np.sqrt((sigma_nm**2 + a_nm**2 / 12) / n + 8 * np.pi * sigma_nm**4 * bg_sd_e**2 / (a_nm**2 * n**2)))


def _source_params_text(E: ERR.ErrorSet, sid: str) -> str:
    spec = ERR.SOURCES[sid]["params"]
    if not spec:
        return "no parameters here (uses the optics / loci settings)"
    return "; ".join(f"{k} = {E.params[sid][k]} {spec[k][1]}".strip() for k in spec)


# ----------------------------------------------------------------------------- part 1
def _ground_truth_numbers(im, run: Path, t: int, cell: int | None, truth_img, tx, occ, cfg):
    px_nm = im.px_um * 1000
    it = _csv(run, "stage3_microscopy", "imaging_truth.csv")
    r = it[it.t == t].iloc[0]
    items = [("Stage drift this frame", f"{np.hypot(r.drift_y_px, r.drift_x_px) * px_nm:.0f} nm", f"({r.drift_y_px:+.2f}, {r.drift_x_px:+.2f}) px in (y, x); everything in the frame moved by this"),
             ("Focus error", f"{r.defocus_um * 1000:+.0f} nm",
              "z-stack: every object appears displaced along z by this much" if im.is3d else "distance of the focal plane from the sample plane in this frame")]
    if im.is3d and "drift_z_um" in r:
        items.insert(1, ("Axial drift this frame", f"{r.drift_z_um * 1000:+.0f} nm", "the sample slid along z by this much since the first frame"))
    for c in im.channels:
        items.append((f"Blur σ, {CH_SHORT[c]}", f"{r[f'sigma_px_{c}'] * px_nm:.0f} nm",
                      "total Gaussian blur this frame (diffraction, defocus, locus motion in quadrature)"))
    for c in im.channels:
        items.append((f"Light level, {CH_SHORT[c]}", f"{100 * r[f'flicker_{c}'] * r[f'bleach_{c}']:.0f} %",
                      "excitation flicker x photobleaching remaining (100% = start of movie, steady lamp)"))
    style.stat_grid(st, items)
    if cell is None or truth_img is None:
        return
    rows = []
    for k in range(len([c for c in im.channels if c.startswith("locus")])):
        q = truth_img[(truth_img.cell_id == cell) & (truth_img.t == t) & (truth_img.locus_id == k)]
        if q.empty:
            continue
        q = q.iloc[0]
        row = {"locus": f"locus {k}", "true x (µm)": round(q.x_um, 3), "true y (µm)": round(q.y_um, 3)}
        if "z_um" in q:
            row["true z (µm)"] = round(q.z_um, 3)
        if "x_drawn_um" in q:
            row["drawn offset from truth, lateral (nm)"] = round(np.hypot(q.x_drawn_um - q.x_um, q.y_drawn_um - q.y_um) * 1000, 0)
        if "z_drawn_um" in q:
            row["drawn offset from truth, axial (nm)"] = round((q.z_drawn_um - q.z_um) * 1000, 0)
        if occ is not None:
            o = occ[(occ.cell_id == cell) & (occ.t == t) & (occ.locus_id == k)]
            if len(o):
                row["probes attached (mean over exposure)"] = round(float(o.mean_attached.iloc[0]), 1)
        rows.append(row)
    if rows:
        st.markdown(f"**Ground truth for cell {cell} in this frame** (positions in image coordinates, i.e. after stage drift)")
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    if tx is not None:
        q = tx[(tx.cell_id == cell) & (tx.t == t)]
        if len(q):
            q = q.iloc[0]
            st.caption(f"Transcription truth: promoter {'ON' if q.promoter_on else 'OFF'}; promoter-enhancer distance {q.d_nm:.0f} nm "
                       f"({'in contact' if q.in_contact else 'not in contact'}); {q.n_polII:.0f} Pol II on the gene; MS2 signal "
                       f"{q.ms2_loops:.1f} loops → this is the bright spot in the nuclear channel at the promoter.")


def _part_view(R, run: Path, im, cfg):
    st.markdown("### 1. The recorded image and the truth behind it")
    st.caption("Left of each pair: what the camera recorded (exactly the saved TIFF). Optionally next to it: the same frame with every error "
               "source switched off, and the difference between the two, which is everything the errors added.")
    T = im.T
    c = st.columns([3, 2, 2, 2])
    t = c[0].slider("Frame", 0, T - 1, min(T // 2, T - 1), key="img_t", help="Which frame of the movie to show.")
    chans = c[1].multiselect("Channels", im.channels, default=im.channels, format_func=CH_SHORT.get, key="img_ch",
                             help="Each colour is a separate camera image.")
    hi = c[2].slider("Brightness ceiling (percentile)", 95.0, 100.0, 99.8, 0.1, key="img_hi",
                     help="Pixels above this percentile of the frame are shown as white. Lower it to see dim structure.")
    lg = c[3].checkbox("Log scale", False, key="img_log", help="Compress bright values so dim and bright objects are both visible.")
    zsel = -1
    if im.is3d:
        cz = st.columns([3, 2, 3])
        proj = cz[1].checkbox("Show the projection", False, key="img_proj",
                              help=f"Show the {im.projection}-intensity projection along z (the 2D image segmentation and tracking used) instead of one slice.")
        zsel = -1 if proj else int(cz[0].slider("z slice", 0, im.nz - 1, im.nz // 2, key="img_z",
                                                help="Each frame is a stack of planes; this chooses which plane to look at."))
        cz[2].caption(f"Stack of {im.nz} planes, {im.zstep:g} µm apart, z = {im.zk[0]:.1f} to {im.zk[-1]:.1f} µm above the coverslip."
                      + ("" if proj else f" Showing plane {zsel + 1} at z = {im.zk[zsel]:.2f} µm."))
    cells_now = R["cells"][R["cells"].t == t]
    c = st.columns([2, 2, 3])
    mode = c[0].radio("View", ["field", "cell", "locus"], horizontal=True, key="img_mode",
                      format_func={"field": "Whole field", "cell": "One cell", "locus": "Zoom on a locus"}.get)
    cell = None
    if mode != "field" and len(cells_now):
        cell = c[1].selectbox("Cell", cells_now.cell_id.tolist(), key="img_cell")
    elif len(cells_now):
        cell = int(cells_now.cell_id.iloc[0])
    lcs = [x for x in im.channels if x.startswith("locus")]
    lk = 0
    if mode == "locus" and lcs:
        lk = c[2].radio("Which locus", list(range(len(lcs))), horizontal=True, key="img_lk",
                        format_func=lambda k: f"locus {k}" + (" (promoter)" if k == 0 else " (enhancer)" if k == 1 else ""))
    truth_img = _csv(run, "stage3_microscopy", "loci_truth_image.csv")
    tx, occ = _csv(run, "stage1_chromatin", "transcription.csv"), _csv(run, "stage1_chromatin", "probe_occupancy.csv")
    show_over = st.multiselect("Draw the ground truth on the images", ["masks", "ids", "loci", "ms2", "fid"],
                               default=["masks", "ids", "loci", "ms2", "fid"], key="img_over",
                               format_func={"masks": "true nucleus outlines", "ids": "true cell ID numbers",
                                            "loci": "true locus positions (○ truth, × where it is drawn)",
                                            "ms2": "promoter with an MS2 spot (◇)", "fid": "fiducial beads on the slide (□)"}.get)
    c = st.columns(2)
    ideal = c[0].checkbox("Also show the error-free version", True, key="img_ideal",
                          help="The same frame rendered with every source of image error switched off.")
    diff = c[1].checkbox("Also show recorded − error-free (what the errors added)", False, key="img_diff")

    cell_row = cells_now[cells_now.cell_id == cell].iloc[0] if cell is not None and (cells_now.cell_id == cell).any() else None
    spot_xy = None
    if truth_img is not None and mode == "locus" and cell is not None:
        q = truth_img[(truth_img.cell_id == cell) & (truth_img.t == t) & (truth_img.locus_id == lk)]
        if len(q):
            q = q.iloc[0]
            spot_xy = (q.get("y_drawn_um", q.y_um) / im.px_um, q.get("x_drawn_um", q.x_um) / im.px_um)
    y0, y1, x0, x1 = bbox = _crop_box(mode if mode != "locus" else "locus", im.shape, im.px_um, cell_row, spot_xy)
    if not chans:
        st.info("Choose at least one channel.")
        return
    lab = _labels_at(run, im, t, zsel)
    rows = 1 + int(ideal) + int(diff)
    fig, axes = plt.subplots(rows, len(chans), figsize=(5.6 * len(chans), 5.4 * rows), squeeze=False)
    stamp = (run / "stage7_validation" / "metrics.json").stat().st_mtime
    for j, ch in enumerate(chans):
        rec_full = _recorded(run, stamp, im, ch, t, zsel).astype(float)
        rec = rec_full[y0:y1, x0:x1]
        lo_, hi_ = _limits(rec_full, hi)
        shown = [("Recorded: " + CH_SHORT[ch], rec)]
        idl = None
        if ideal or diff:
            idl = _render(str(run), stamp, t, ch, (), bbox, zsel).astype(float)
        if ideal:
            shown.append(("Error-free: " + CH_SHORT[ch], idl))
        if diff:
            shown.append(("Recorded − error-free", rec - idl))
        for i, (title, arr) in enumerate(shown):
            ax = axes[i, j]
            if title.startswith("Recorded −"):
                m = max(np.percentile(np.abs(arr), 99.5), 1e-9)
                ax.imshow(arr, cmap="RdBu_r", vmin=-m, vmax=m, extent=(x0 - .5, x1 - .5, y1 - .5, y0 - .5))
            else:
                a = np.log1p(np.clip(arr - lo_, 0, None)) if lg else arr
                lo2, hi2 = (0, np.log1p(max(hi_ - lo_, 1e-9))) if lg else (lo_, hi_)
                ax.imshow(a, cmap="gray", vmin=lo2, vmax=hi2, extent=(x0 - .5, x1 - .5, y1 - .5, y0 - .5))
            ax.set_title(title, fontsize=13)
            ax.set_xticks([]); ax.set_yticks([])
            if not title.startswith("Recorded −"):
                _overlay(ax, lab, truth_img, tx, im, t, ch, cell, show_over, bbox, cells_now, zsel)
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    st.caption("Grey scale: black = low, white = bright (same scale for the recorded and error-free panel of a channel). "
               "Cyan ○ = where the locus truly is in the image; red × = where that colour's spot is actually drawn (differs by the "
               "chromatic shift); yellow outline = true nucleus mask; ◇ = promoter with an active MS2 signal.")
    if (y1 - y0, x1 - x0) != tuple(im.shape):
        st.caption(f"Zoomed region: x {x0}–{x1} px, y {y0}–{y1} px ({(x1 - x0) * im.px_um:.1f} × {(y1 - y0) * im.px_um:.1f} µm).")
    # does the re-rendered frame equal the saved TIFF?
    ch0 = chans[0]
    full = _render(str(run), stamp, t, ch0, tuple(im.E.enabled_ids()), (0, im.shape[0], 0, im.shape[1]), zsel)
    saved = _recorded(run, stamp, im, ch0, t, zsel)
    ok = np.array_equal(np.rint(full).clip(0, 65535).astype(np.uint16), saved.astype(np.uint16))
    (st.success if ok else st.warning)(
        "The rebuilt image with all this run's error sources on is identical to the saved image (so the build-up below is exactly "
        "how this image was made)." if ok else
        "The rebuilt image differs from the saved one (the run was made with other code or settings), so the build-up below is "
        "indicative only.")
    _ground_truth_numbers(im, run, t, cell, truth_img, tx, occ, cfg)
    if im.is3d:
        _side_views(run, stamp, im, t, chans, cell, lk, cells_now, truth_img, hi, ideal)
    return t, ch0, cell, lk, zsel


def _overlay(ax, lab, truth_img, tx, im, t, ch, cell, show_over, bbox, cells_now, zsel=-1):
    y0, y1, x0, x1 = bbox
    if "masks" in show_over and lab.any():
        sub = lab[y0:y1, x0:x1]
        b = find_boundaries(sub, mode="inner")
        o = np.zeros(sub.shape + (4,)); o[b] = (1.0, 0.9, 0.0, 0.9)
        ax.imshow(o, interpolation="nearest", extent=(x0 - .5, x1 - .5, y1 - .5, y0 - .5))
    if "ids" in show_over:
        for r in cells_now.itertuples():
            cx, cy = r.x_um / im.px_um, r.y_um / im.px_um
            if x0 <= cx < x1 and y0 <= cy < y1:
                ax.text(cx, cy, str(int(r.cell_id)), color="yellow", fontsize=11, ha="center", va="center", weight="bold")
    if "fid" in show_over:
        fp = Path(st.session_state.get("_run_path", "")) / "stage3_microscopy" / "fiducials.csv"
        if fp.exists():
            for q in pd.read_csv(fp).itertuples():
                if x0 <= q.x_um / im.px_um < x1 and y0 <= q.y_um / im.px_um < y1:
                    ax.plot(q.x_um / im.px_um, q.y_um / im.px_um, "s", mfc="none", mec="magenta", ms=13, mew=1.5)
    if truth_img is None:
        return
    tt = truth_img[truth_img.t == t]
    if "loci" in show_over:
        for q in tt.itertuples():
            if ch.startswith("locus") and q.locus_id != int(ch[5:]):
                continue
            dz = (q.z_um - im.zk[zsel]) if (im.is3d and zsel >= 0 and hasattr(q, "z_um")) else 0.0
            near = abs(dz) <= im.zstep / 2 + 1e-9                    # in a z-slice, only loci within half a step are "in this plane"
            alpha = 1.0 if near else 0.4
            if x0 <= q.x_um / im.px_um < x1 and y0 <= q.y_um / im.px_um < y1:
                ax.plot(q.x_um / im.px_um, q.y_um / im.px_um, "o", mfc="none", mec="cyan", ms=15, mew=1.6, alpha=alpha, ls="-" if near else "--")
                if not near:
                    ax.annotate(f"z {dz:+.1f} µm", (q.x_um / im.px_um, q.y_um / im.px_um), color="cyan", fontsize=8, xytext=(8, 8), textcoords="offset points")
            if hasattr(q, "x_drawn_um") and ch.startswith("locus"):
                ax.plot(q.x_drawn_um / im.px_um, q.y_drawn_um / im.px_um, "x", color="red", ms=9, mew=1.6, alpha=alpha)
    if "ms2" in show_over and ch == "nucleus" and tx is not None:
        on = tx[(tx.t == t) & (tx.ms2_loops > 0.5)]
        for q in tt[tt.locus_id == 0].itertuples():
            if ((on.cell_id == q.cell_id)).any() and x0 <= q.x_um / im.px_um < x1 and y0 <= q.y_um / im.px_um < y1:
                ax.plot(q.x_um / im.px_um, q.y_um / im.px_um, "D", mfc="none", mec="magenta", ms=14, mew=1.8)
    ax.set_xlim(x0 - .5, x1 - .5); ax.set_ylim(y1 - .5, y0 - .5)


def _side_views(run, stamp, im, t, chans, cell, lk, cells_now, truth_img, hi, ideal):
    """x-z and y-z cuts through the chosen locus (or the cell centre): the axial truth next to the recorded stack."""
    st.markdown("**Side views: x–z and y–z cuts through the chosen locus (or the cell's centre)**")
    st.caption("Horizontal = position along the cut, vertical = z (the planes of the stack). Cyan ○ = the true position, red × = where the spot is "
               "really drawn (true position + chromatic shift + focus error). Blur along z is the axial point spread function.")
    if cell is None:
        st.info("Choose a cell above to see its side views.")
        return
    row = cells_now[cells_now.cell_id == cell]
    if row.empty:
        return
    yc, xc = row.iloc[0].y_um / im.px_um, row.iloc[0].x_um / im.px_um
    zt, zd = float(row.iloc[0].get("z_um", im.zcenter)), None
    xd = yd = None
    if truth_img is not None:
        q = truth_img[(truth_img.cell_id == cell) & (truth_img.t == t) & (truth_img.locus_id == lk)]
        if len(q) and any(c.startswith("locus") for c in chans):
            q = q.iloc[0]
            yc, xc, zt = q.y_um / im.px_um, q.x_um / im.px_um, float(q.get("z_um", zt))
            yd, xd, zd = q.get("y_drawn_um", q.y_um) / im.px_um, q.get("x_drawn_um", q.x_um) / im.px_um, q.get("z_drawn_um", zt)
    w = 26
    iy, ix = int(round(yc)), int(round(xc))
    y0, y1, x0, x1 = max(iy - w, 0), min(iy + w + 1, im.shape[0]), max(ix - w, 0), min(ix + w + 1, im.shape[1])
    bbox = (y0, y1, x0, x1)
    zext = (im.zk[-1] + im.zstep / 2, im.zk[0] - im.zstep / 2)
    per = 1 + int(ideal)                                               # recorded (+ error-free) for each of the two cuts
    fig, axes = plt.subplots(len(chans), 2 * per, figsize=(4.2 * 2 * per, 3.5 * len(chans)), squeeze=False)
    for j, ch in enumerate(chans):
        rec = _read_stack(str(run), stamp, f"{ch}_zstack", t, im.nz)[:, y0:y1, x0:x1].astype(float)
        vols = [("Recorded", rec)]
        if ideal:
            vols.append(("Error-free", _render_vol(str(run), stamp, t, ch, (), bbox).astype(float)))
        lo, hh = float(np.percentile(rec, 1)), float(np.percentile(rec, hi))
        for k, axname in enumerate(("x–z", "y–z")):
            for i, (name, vol) in enumerate(vols):
                cut = vol[:, iy - y0 - 1:iy - y0 + 2, :].mean(axis=1) if k == 0 else vol[:, :, ix - x0 - 1:ix - x0 + 2].mean(axis=2)
                ax = axes[j, k * per + i]
                u0, u1 = ((x0, x1) if k == 0 else (y0, y1))
                ax.imshow(cut, cmap="gray", vmin=lo, vmax=hh, aspect="auto", extent=(u0 * im.px_um, u1 * im.px_um, zext[0], zext[1]))
                ax.set_title(f"{name}: {CH_SHORT[ch]}, {axname}", fontsize=11)
                if ch == "nucleus" or "loci" in st.session_state.get("img_over", ["loci"]):
                    ax.plot((xc if k == 0 else yc) * im.px_um, zt, "o", mfc="none", mec="cyan", ms=14, mew=1.5)
                    if zd is not None and ch.startswith("locus"):
                        ax.plot(((xd if k == 0 else yd)) * im.px_um, zd, "x", color="red", ms=9, mew=1.5)
                ax.set_xlabel("µm"); ax.set_ylabel("z (µm)")
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)


# ----------------------------------------------------------------------------- part 2
def _part_buildup(R, run: Path, im, ctx):
    st.markdown("### 2. How this image was built, one error at a time")
    st.caption("Start from the ground truth drawn with perfect optics and no noise, then add the run's error sources in the physical order "
               "in which light is degraded. Every panel adds one source to the previous panel; the last panel is the recorded image. "
               "Each panel uses the same grey scale, so brightness changes are real.")
    t, ch0, cell, lk, _ = ctx
    if im.channels[-1].startswith("locus"):                      # the loci show the errors best: start there
        st.session_state.setdefault("bu_ch", next(c for c in im.channels if c.startswith("locus")))
        st.session_state.setdefault("bu_mode", "locus")
    c = st.columns([2, 2, 2])
    ch = c[0].selectbox("Channel", im.channels, format_func=CH_SHORT.get, key="bu_ch")
    mode = c[1].radio("Region", ["cell", "locus"], horizontal=True, key="bu_mode",
                      format_func={"cell": "One cell", "locus": "Zoom on a locus"}.get)
    if mode == "locus" and not ch.startswith("locus"):
        mode = "cell"
        c[1].caption("The nuclear channel has no locus; showing the cell.")
    lkk = int(ch[5:]) if ch.startswith("locus") else lk
    cells_now = R["cells"][R["cells"].t == t]
    if cells_now.empty:
        st.info("No cells in this frame.")
        return
    cell = c[2].selectbox("Cell", cells_now.cell_id.tolist(), index=max(0, cells_now.cell_id.tolist().index(cell)) if cell in cells_now.cell_id.tolist() else 0,
                          key="bu_cell")
    cell_row = cells_now[cells_now.cell_id == cell].iloc[0]
    ti = _csv(run, "stage3_microscopy", "loci_truth_image.csv")
    spot = None
    if ti is not None:
        q = ti[(ti.cell_id == cell) & (ti.t == t) & (ti.locus_id == lkk)]
        if len(q):
            spot = (q.iloc[0].get("y_drawn_um", q.iloc[0].y_um) / im.px_um, q.iloc[0].get("x_drawn_um", q.iloc[0].x_um) / im.px_um)
    bbox = _crop_box(mode, im.shape, im.px_um, cell_row, spot)
    y0, y1, x0, x1 = bbox
    stamp = (run / "stage7_validation" / "metrics.json").stat().st_mtime
    zsel = -1
    if im.is3d:
        zref = float(cell_row.get("z_um", im.zcenter))
        if ti is not None and len(q):
            zref = float(q.iloc[0].get("z_um", zref))
        k0 = int(np.argmin(np.abs(im.zk - zref)))
        zsel = int(st.slider("z slice", 0, im.nz - 1, k0, key="bu_z", help="The plane to build up. It starts at the plane nearest the object's true height."))
    if not st.button("▶ Build the image up", type="primary", key="bu_go") and "bu_done" not in st.session_state:
        st.caption("Press the button: it renders about 20 versions of this frame (a few seconds).")
        return
    st.session_state["bu_done"] = True
    order = [s for s in ERR.ORDER if im.E.on(s)]
    stages = [("Ground truth only: perfect optics, no noise, no errors", ())]
    for i, s in enumerate(order):
        stages.append((ERR.SOURCES[s]["title"], tuple(order[: i + 1])))
    with st.spinner("Rendering the stages…"):
        imgs = [_render(str(run), stamp, t, ch, en, bbox, zsel).astype(float) for _, en in stages]
    lo, hi = _limits(imgs[-1], 99.9)
    import textwrap
    ncol = 5
    nrow = int(np.ceil(len(stages) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.9 * nrow), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for k, ((title, en), img) in enumerate(zip(stages, imgs)):
        ax = axes.ravel()[k]
        ax.imshow(img, cmap="gray", vmin=lo, vmax=hi)
        ax.set_title(textwrap.fill(f"{k}  +{title}" if k else "0  Ground truth", 26), fontsize=11)
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    rows = [{"step": 0, "adds": "Ground truth (perfect optics, no noise)", "what it does": "The objects drawn where they truly are, as light intensity.",
             "change vs previous (RMS, ADU)": None}]
    for k in range(1, len(stages)):
        s = order[k - 1]
        rows.append({"step": k, "adds": ERR.SOURCES[s]["title"], "what it does": ERR.SOURCES[s]["effect"],
                     "change vs previous (RMS, ADU)": round(float(np.sqrt(np.mean((imgs[k] - imgs[k - 1]) ** 2))), 2)})
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch", column_config={
        "what it does": st.column_config.TextColumn(width="large")})
    if ch.startswith("locus") or mode == "locus":
        _profiles(stages, imgs, order, ch)


def _profiles(stages, imgs, order, ch):
    st.markdown("**Brightness profile through the brightest point** (left-right line through the peak of the error-free image)")
    opt_end = max([i + 1 for i, s in enumerate(order) if ERR.SOURCES[s]["group"] == "Optics"] or [0])
    illum_end = max([i + 1 for i, s in enumerate(order) if ERR.SOURCES[s]["group"] in ("Optics", "Illumination")] or [0])
    choices = {0: "ground truth", opt_end: "after optics", illum_end: "after illumination", len(stages) - 1: "recorded"}
    pick = st.multiselect("Stages to compare", list(choices), default=list(choices), format_func=lambda i: f"{i}: {choices[i]}", key="bu_prof")
    y, x = np.unravel_index(np.argmax(imgs[0]), imgs[0].shape)
    fig, ax = plt.subplots(figsize=(9, 3.4))
    for i in sorted(pick):
        ax.plot(imgs[i][y], label=choices[i], lw=2.2 if i == len(stages) - 1 else 1.6)
    ax.set_xlabel("pixel along the line"); ax.set_ylabel("camera counts (ADU)"); ax.legend(); ax.set_title("Profile through the brightest spot")
    st.pyplot(fig, width="stretch")
    plt.close(fig)


# ----------------------------------------------------------------------------- part 3
def _part_budget(R, run: Path, im, ctx):
    st.markdown("### 3. How much each error matters on its own")
    st.caption("For the channel and region chosen in part 2: the error-free image plus ONE source at a time, compared with the error-free "
               "image. Size = RMS difference in camera counts; for a locus also how far the spot's measured centre moves.")
    if "bu_ch" not in st.session_state:
        st.caption("Choose a channel and region in part 2 first.")
        return
    t, _, cell, _, _ = ctx
    zsel = int(st.session_state.get("bu_z", -1)) if im.is3d else -1
    ch, mode = st.session_state["bu_ch"], st.session_state.get("bu_mode", "cell")
    if not st.button("▶ Measure every source", key="eb_go"):
        return
    cells_now = R["cells"][R["cells"].t == t]
    cid = st.session_state.get("bu_cell", cell)
    if not (cells_now.cell_id == cid).any():
        return
    cell_row = cells_now[cells_now.cell_id == cid].iloc[0]
    ti = _csv(run, "stage3_microscopy", "loci_truth_image.csv")
    lkk = int(ch[5:]) if ch.startswith("locus") else 0
    spot = None
    if ti is not None:
        q = ti[(ti.cell_id == cid) & (ti.t == t) & (ti.locus_id == lkk)]
        if len(q):
            spot = (q.iloc[0].y_um / im.px_um, q.iloc[0].x_um / im.px_um)          # drift-corrected truth
    mode = mode if ch.startswith("locus") or mode == "cell" else "cell"
    bbox = _crop_box(mode, im.shape, im.px_um, cell_row, (spot[0], spot[1]) if spot else None)
    stamp = (run / "stage7_validation" / "metrics.json").stat().st_mtime
    base = _render(str(run), stamp, t, ch, (), bbox, zsel).astype(float)
    vol0 = _render_vol(str(run), stamp, t, ch, (), bbox).astype(float) if (im.is3d and ch.startswith("locus")) else None
    peak = max(float(base.max() - np.median(base)), 1e-9)
    py, px = np.unravel_index(np.argmax(base), base.shape)
    cy0, cx0 = _centroid(base, py, px)
    rows = []
    prog = st.progress(0.0, text="Measuring…")
    srcs = [s for s in ERR.ORDER if im.E.on(s)]
    for i, s in enumerate(srcs):
        img = _render(str(run), stamp, t, ch, (s,), bbox, zsel).astype(float)
        rms = float(np.sqrt(np.mean((img - base) ** 2)))
        row = {"source": ERR.SOURCES[s]["title"], "group": ERR.SOURCES[s]["group"], "RMS change (ADU)": round(rms, 2),
               "% of the brightest signal": round(100 * rms / peak, 2)}
        if ch.startswith("locus"):
            cy, cx = _centroid(img, py, px)
            row["spot centre moves by (nm)"] = round(float(np.hypot(cy - cy0, cx - cx0) * im.px_um * 1000), 1) if np.isfinite(cy) else None
            if vol0 is not None:
                v1 = _render_vol(str(run), stamp, t, ch, (s,), bbox).astype(float)
                row["spot centre moves along z by (nm)"] = round(float((_zcentroid(v1, py, px) - _zcentroid(vol0, py, px)) * im.zstep * 1000), 1)
        rows.append(row)
        prog.progress((i + 1) / len(srcs), text=f"Measured {ERR.SOURCES[s]['title']}")
    prog.empty()
    df = pd.DataFrame(rows).sort_values("RMS change (ADU)", ascending=False)
    fig, ax = plt.subplots(figsize=(9, 0.38 * len(df) + 1.2))
    ax.barh(df.source[::-1], df["RMS change (ADU)"][::-1], color="#0f766e")
    ax.set_xlabel("RMS difference from the error-free image (ADU)"); ax.set_title("Which sources change this image most?")
    st.pyplot(fig, width="stretch")
    plt.close(fig)
    st.dataframe(df, hide_index=True, width="stretch")
    st.caption("Sources that only move things (stage drift, chromatic shift) do not change the brightness pattern much but they displace "
               "the object; read their effect from 'spot centre moves by'. Effects add up only roughly: noise terms add in quadrature.")


# ----------------------------------------------------------------------------- part 4
def _part_reference(im, cfg):
    st.markdown("### 4. Every source of image error, explained")
    E = im.E
    st.caption(f"This run: {len(E.enabled_ids())} of {len(ERR.ORDER)} sources on "
               + ("(the image-error block was switched on)." if E.master else "(only the original model; the image-error block was off)."))
    px_nm, opt, lc, acq = im.px_um * 1000, cfg["optics"], cfg.get("loci") or {}, cfg["acquisition"]
    items = [("Pixel size", f"{px_nm:.0f} nm", "size of one pixel in the sample")]
    for c in im.channels:
        wl = opt["wavelength_nm"] if c == "nucleus" else lc["wavelengths_nm"][int(c[5:])]
        s_nm = psf_sigma_px({**opt, "wavelength_nm": wl}) * px_nm
        items.append((f"PSF width, {CH_SHORT[c]}", f"σ {s_nm:.0f} nm · FWHM {2.355 * s_nm:.0f} nm",
                      f"diffraction blur at {wl} nm and NA {opt['NA']}; FWHM / pixel = {2.355 * s_nm / px_nm:.1f} samples across the spot"))
    if lc and "photons_per_locus_s" in lc:
        n = lc["photons_per_locus_s"] * acq["exposure_s"] * opt.get("quantum_efficiency", 0.8)
        bg = np.sqrt(opt["background_photons"] * opt.get("quantum_efficiency", 0.8) + opt["read_noise_e"] ** 2)
        s_nm = psf_sigma_px({**opt, "wavelength_nm": lc["wavelengths_nm"][0]}) * px_nm
        items.append(("Best-case locus precision", f"{_thompson_nm(s_nm, px_nm, n, bg):.0f} nm",
                      f"theory limit for a spot of {n:.0f} detected photons (photon shot noise + background + read noise only); "
                      "defocus, drift, chromatic shift etc. make it worse"))
    if E.on("stage_drift"):
        steps = np.array([im.drift_px(t) for t in range(im.T)])
        items.append(("Stage drift over the movie", f"{np.hypot(*(steps[-1] - steps[0])) * px_nm:.0f} nm net",
                      f"largest excursion {np.abs(steps).max() * px_nm:.0f} nm"))
    if E.on("chromatic_shift") and len([c for c in im.channels if c.startswith('locus')]) >= 2:
        a, b = np.array(im.channel_shift_px("locus0")), np.array(im.channel_shift_px("locus1"))
        items.append(("Locus 0 ↔ locus 1 misregistration", f"{np.hypot(*(a - b)) * px_nm:.0f} nm",
                      "constant offset added to every measured distance between the two colours (uncorrected)"))
    if im.is3d:
        zs = [("z-stack", f"{im.nz} planes × {im.zstep:g} µm", f"covers z = {im.zk[0]:.1f} to {im.zk[-1]:.1f} µm; a stack takes about {im.nz * im.slice_time_s:.1f} s")]
        for c in im.channels:
            wl = opt["wavelength_nm"] if c == "nucleus" else lc["wavelengths_nm"][int(c[5:])]
            from simlive.stage3_microscopy.render import axial_sigma_um
            sz = axial_sigma_um({**opt}, wl)
            zs.append((f"Axial PSF, {CH_SHORT[c]}", f"σ {sz * 1000:.0f} nm · FWHM {2.355 * sz * 1000:.0f} nm",
                       f"{str(opt.get('modality', 'widefield'))} imaging; FWHM / plane spacing = {2.355 * sz / im.zstep:.1f} planes across the spot"))
        if E.on("chromatic_shift") and len([c for c in im.channels if c.startswith('locus')]) >= 2:
            zs.append(("Locus 0 ↔ locus 1 axial misregistration", f"{abs(im.channel_shift_z_um('locus0') - im.channel_shift_z_um('locus1')) * 1000:.0f} nm",
                       "constant axial offset between the two colours (uncorrected)"))
        items += zs
    style.stat_grid(st, items)
    st.caption("Hover a box for the explanation.")
    for g in ERR.GROUPS:
        ids = [s for s in ERR.ORDER if ERR.SOURCES[s]["group"] == g]
        with st.expander(f"{g} ({sum(E.on(s) for s in ids)} of {len(ids)} on)"):
            for s in ids:
                spec = ERR.SOURCES[s]
                style.badge(st, "ON" if E.on(s) else "off", "good" if E.on(s) else "idle")
                st.markdown(f"**{spec['title']}**" + (" · _part of the original model_" if spec["baseline"] else ""))
                st.write(spec["what"])
                st.caption(f"Effect on the images: {spec['effect']}")
                st.caption(f"Settings: {_source_params_text(E, s)}")
                st.markdown("---")


# ----------------------------------------------------------------------------- tab
def render(st_, run_path: Path, R: dict):
    style.intro(st, "The recorded images next to the ground truth that produced them: where every nucleus and locus truly was, how the "
                    "imaging conditions (drift, focus, light level) varied, and how each source of image error changed the picture.",
                "Part 1: pick a frame and see the camera images with the truth drawn on top (and the error-free version). Part 2: watch the "
                "image being built up one error at a time. Part 3: see which error matters most. Part 4: read what each error is.")
    run = Path(run_path)
    st.session_state["_run_path"] = str(run)
    if not (run / "stage3_microscopy" / "imaging_truth.csv").exists():
        st.info("This run was made before the image-error model existed, so the build-up and the imaging ground truth are not available. "
                "Make a new run (sidebar) to see them. The recorded images are still in the 🎞 Movie tab (tick 'Also show raw image panel').")
        return
    stamp = (run / "stage7_validation" / "metrics.json").stat().st_mtime
    im = _imager(str(run), stamp)
    cfg = R["cfg"]
    ctx = _part_view(R, run, im, cfg)
    if ctx is None:
        return
    _part_buildup(R, run, im, ctx)
    _part_budget(R, run, im, ctx)
    _part_reference(im, cfg)
