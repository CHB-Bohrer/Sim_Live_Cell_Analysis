"""The dashboard's 📍 Stage & fiducials tab: what the z-piezo really did, what the fiducial beads on the slide can and cannot tell, and what
correcting with them (or with the piezo's own position log) does to the locus heights.

Reads the run's saved ground truth (piezo_truth.csv, fiducials.csv, imaging_truth.csv) and the analysis output (fiducial_shifts.csv,
fiducial_corrections.csv). Imports numpy / pandas / matplotlib only, never torch.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import tifffile

import style


def _csv(p: Path):
    return pd.read_csv(p) if p.exists() else None


def _overview() -> None:
    c = st.columns(4)
    with c[0]:
        style.card(st, "The piezo is asked", "A z-stack is made by stepping a piezo through planes at <b>nominal</b> heights, for example every 0.8 µm. "
                   "Every locus height is read from those nominal heights.", 1)
    with c[1]:
        style.card(st, "…but is imperfect", "Each step has a small <b>random</b> error, the response <b>bows</b> over the range, it differs "
                   "going <b>up and down</b> (hysteresis), and it <b>lags</b> just after starting. So the true plane heights differ from the nominal ones.", 2)
    with c[2]:
        style.card(st, "Beads on the slide", "A few bright beads lie on the slide in every image. They show how far the <b>whole picture</b> moved "
                   "(drift, focus). Because they all sit at one height, they only see the planes near the slide.", 3)
    with c[3]:
        style.card(st, "The piezo's own log", "A piezo with a position sensor <b>logs where each plane really was</b>. That is the only thing "
                   "here that corrects the plane-to-plane errors, which change the distance between loci at different heights.", 4)


def _components(pz: pd.DataFrame) -> dict:
    """Estimate the four ingredients of the piezo error from the saved truth (what was actually applied)."""
    up, dn = pz[pz.scan_dir > 0], pz[pz.scan_dir < 0]
    nz = pz.k.nunique()
    mean_by_k = pz.groupby("k").error_um.mean()
    hyst = float(up.error_um.mean() - dn.error_um.mean()) if len(dn) and len(up) else 0.0
    base = pz.assign(res=pz.error_um - pz.groupby(["k", "scan_dir"]).error_um.transform("mean"))
    ends = float(mean_by_k.iloc[[0, -1]].mean())
    return {"rms_total_nm": 1000 * float(np.sqrt((pz.error_um ** 2).mean())), "repeat_nm": 1000 * float(base.res.std()), "bow_um": float(mean_by_k.max() - ends),
            "hyst_um": hyst, "n_frames": int(pz.t.nunique()), "nz": nz}


def _piezo_section(run: Path, cfg: dict) -> None:
    st.markdown("### 1. Where the planes really were")
    pz = _csv(run / "stage3_microscopy" / "piezo_truth.csv")
    if pz is None:
        st.info("This run has no piezo errors (the 'Piezo z-stage errors' source is off), so every plane was exactly at its nominal height.")
        return
    st.caption("Each row is one stack (one frame of the movie), each column one plane. Colour = how far that plane was from its nominal height, in nm "
               "(red = higher than asked, blue = lower). Rows alternate in scan direction when the piezo is bidirectional.")
    E = pz.pivot(index="t", columns="k", values="error_um") * 1000
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.2), gridspec_kw={"width_ratios": [3, 2]})
    m = float(np.abs(E.values).max())
    im = ax[0].imshow(E.values, aspect="auto", cmap="RdBu_r", vmin=-m, vmax=m, origin="lower")
    ax[0].set_xlabel("plane number (0 = lowest)"); ax[0].set_ylabel("frame"); ax[0].set_title("Error of every plane in every stack (nm)")
    fig.colorbar(im, ax=ax[0], fraction=0.03)
    up, dn = pz[pz.scan_dir > 0].groupby("k").error_um.mean() * 1000, pz[pz.scan_dir < 0].groupby("k").error_um.mean() * 1000
    ax[1].plot(up.index, up.values, "o-", color="#f59e0b", label="up-scans (average)")
    if len(dn):
        ax[1].plot(dn.index, dn.values, "o-", color="#60a5fa", label="down-scans (average)")
    ax[1].axhline(0, color="#999", lw=0.8); ax[1].set_xlabel("plane number"); ax[1].set_ylabel("average error (nm)")
    ax[1].set_title("The systematic part: bow and hysteresis"); ax[1].legend(fontsize=9)
    st.pyplot(fig, width="stretch"); plt.close(fig)
    comp = _components(pz)
    p = (cfg.get("imaging_errors") or {}).get("piezo") or {}
    style.stat_grid(st, [("Total error", f"{comp['rms_total_nm']:.0f} nm rms", "root-mean-square distance of the planes from their nominal heights"),
                         ("Random step error", f"{comp['repeat_nm']:.0f} nm", f"independent in every plane of every stack (set: {p.get('repeatability_nm', '?')} nm)"),
                         ("Bow (middle vs ends)", f"{1000 * comp['bow_um']:.0f} nm", f"systematic non-linearity (set: {1000 * p.get('bow_um', 0):.0f} nm)"),
                         ("Hysteresis (up − down)", f"{1000 * comp['hyst_um']:.0f} nm", f"up-scan planes sit higher than down-scan planes (set: {1000 * p.get('hysteresis_um', 0):.0f} nm)"),
                         ("Position log noise", f"{1000 * float(np.std(pz.readout_um - pz.z_actual_um)):.1f} nm", "the piezo's sensor reports the true height plus this noise")])
    st.caption("What this does to an image: a locus recorded in plane 7 is assigned the height of nominal plane 7, but the plane was really "
               "a few tens of nm away. The height of a locus is therefore wrong by about the error of the planes around it. The common (average) part moves "
               "both loci together and cancels in their distance; the random part does not.")


def _fiducial_section(run: Path, R: dict) -> None:
    st.markdown("### 2. The fiducial beads on the slide")
    fid = _csv(run / "stage3_microscopy" / "fiducials.csv")
    if fid is None:
        st.info("This run has no fiducial beads (set `fiducials.enabled: true`, or tick 'Fiducial beads on the slide' in the sidebar form).")
        return
    ej = run / "stage3_microscopy" / "imaging_errors.json"
    import json
    zs = (json.loads(ej.read_text()).get("z_stack") or {}) if ej.exists() else {}
    zlist = zs.get("z_positions_um", [])
    px_um = R["cfg"]["optics"]["pixel_size_nm"] / 1000.0
    c = st.columns([3, 2])
    img_path = run / "stage3_microscopy" / "nucleus.tif"
    with c[0]:
        img = tifffile.imread(img_path, key=0).astype(float)
        fig, ax = plt.subplots(figsize=(6, 5.6))
        ax.imshow(img, cmap="gray", vmin=np.percentile(img, 1), vmax=np.percentile(img, 99.7))
        ax.plot(fid.x_um / px_um, fid.y_um / px_um, "s", mfc="none", mec="magenta", ms=14, mew=1.6)
        for r in fid.itertuples():
            ax.annotate(f"bead {r.bead}", (r.x_um / px_um, r.y_um / px_um), color="magenta", xytext=(8, -12), textcoords="offset points", fontsize=9)
        ax.set_title("Frame 0, nuclear colour (maximum projection): beads marked", fontsize=11); ax.axis("off")
        st.pyplot(fig, width="stretch"); plt.close(fig)
    with c[1]:
        reach = bool(zlist) and zlist[0] - 0.5 * zs.get("step_um", 0.8) <= fid.z_um.min() + 0.2
        st.markdown(f"**{len(fid)} beads**, lying on the slide at z ≈ {fid.z_um.mean():.2f} µm above the coverslip.")
        if zlist:
            st.markdown(f"The stack covers z = {zlist[0]:.1f} to {zlist[-1]:.1f} µm: " + ("**it reaches the slide, so the beads are in focus in the lowest planes.**" if reach else
                        "**it does NOT reach the slide, so the beads are almost invisible (only faint out-of-focus haze).** Lower the stack's centre or add planes."))
        st.markdown("**What the beads can tell you:** how far the whole picture moved in x, y and z in each frame (stage drift and focus error): "
                    "all beads and all loci moved together, so the median bead displacement is subtracted from every locus.")
        st.markdown("**What they cannot tell you:** the errors of the planes high up near the nuclei. All beads are at one height, so they sample only "
                    "the lowest planes.")
    sh = _csv(run / "stage6_analysis" / "truth" / "fiducial_shifts.csv")
    it = _csv(run / "stage3_microscopy" / "imaging_truth.csv")
    if sh is None or it is None:
        st.info("The bead tracking has not been run for this run.")
        return
    st.markdown("**Whole-picture displacement of each frame: measured from the beads versus what really happened**")
    t = sh.t.to_numpy()
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.6))
    tz = (it.drift_z_um + it.defocus_um).to_numpy() if "drift_z_um" in it else np.zeros(len(it))
    for a, (name, est, true, unit) in zip(ax, [("z (µm)", sh.dz_um, tz - tz[0], "µm"), ("y (px)", sh.dy_px, it.drift_y_px - it.drift_y_px[0], "px"),
                                                 ("x (px)", sh.dx_px, it.drift_x_px - it.drift_x_px[0], "px")]):
        a.plot(t, true, "-", color="#999", lw=3, label="truth (drift + focus)")
        a.plot(t, est, "o-", color="#14b8a6", ms=4, label="measured from the beads")
        a.set_title(f"Displacement in {name}", fontsize=11); a.set_xlabel("frame"); a.legend(fontsize=8)
    st.pyplot(fig, width="stretch"); plt.close(fig)
    st.caption(f"Beads followed: {int(sh.n_beads.iloc[0])}. The z estimate includes the piezo error of the planes near the slide, which is why it can differ a little "
               "from the pure drift + focus truth.")


def _effect_section(run: Path) -> None:
    st.markdown("### 3. Does correcting help?")
    tab = _csv(run / "stage6_analysis" / "truth" / "fiducial_corrections.csv")
    if tab is None:
        st.info("No correction results for this run (needs a z-stack run with fiducial beads).")
        return
    show = tab.copy()
    for c in show.columns[2:]:
        show[c] = show[c].round(1)
    st.dataframe(show, hide_index=True, width="stretch")
    st.caption("Medians of the absolute errors (so a few wrong detections by the placeholder locus finder do not hide the effect). "
               "**Distance error**: error of the 3D distance between locus 0 and locus 1, a shift common to both cannot change it. "
               "**Height-difference error**: error of their height difference. **Height error after removing the constant offset**: how much each "
               "locus' height wanders, with the movie-wide constant taken out. Rows: uncorrected; subtracting the bead-measured picture shift; "
               "replacing the nominal plane heights by the piezo's logged heights; both.")
    base = tab[tab.correction == "uncorrected"].iloc[0]
    notes = []
    for r in tab.itertuples():
        if r.correction == "uncorrected":
            continue
        d0, d1 = base["height-difference error, typical (nm)"], getattr(r, "_4")
        notes.append(f"**{r.correction}**: height-difference error {d0:.0f} → {d1:.0f} nm.")
    st.info("  \n".join(notes) if notes else "Only the uncorrected result exists.")


def _reference() -> None:
    with st.expander("How each piece works (and its assumptions)"):
        st.markdown(
            "- **Piezo error model** (`stage3_microscopy/render.py`, `piezo_error_um`): error of plane k = random N(0, repeatability) + bow × (1 − u²) with u from −1 to +1 across the stack + "
            "direction × hysteresis / 2 − direction × settling × exp(−order of recording / 2). The scan direction alternates every frame if the piezo is bidirectional. "
            "The images are drawn with the planes at their true heights; the analysis is told the nominal heights.\n"
            "- **Fiducial beads** (`stage3_microscopy/fiducials.py`): a few bright sub-resolution beads on the slide, placed away from every nucleus, in every colour, "
            "shifted by the same drift as everything else.\n"
            "- **Bead tracking** (`stage6_analysis/fiducials.py`): beads are found in frame 0 outside the nuclei, fitted in every frame (2D Gaussian laterally, a Gaussian along z "
            "through the planes) and the median displacement from frame 0 is the whole-picture shift.\n"
            "- **Piezo log**: the true plane heights plus the sensor's noise (a few nm). Using it replaces the nominal heights in two places: the bead heights are measured with "
            "the logged heights, and the locus heights are interpolated between the logged plane heights.\n"
            "- **Limits**: beads on the slide only constrain the lowest planes; random plane errors at nuclear height can only be removed with the log (or beads embedded at "
            "nuclear height); the placeholder locus finder's own errors (cosmic rays, hot pixels) are larger than these corrections unless the typical error is looked at.")


def render(st_, run_path: Path, R: dict):
    style.intro(st, "The z-stage that steps through the planes of a stack is imperfect, and the beads on the slide are the tool for measuring part of that.",
                "Read the four cards, then see in order: how far each plane really was from where it was asked to be, the beads and what they can and cannot measure, "
                "and a table of what each correction does to the locus distances and heights.")
    run = Path(run_path)
    _overview()
    if not (R["cfg"].get("acquisition", {}).get("z_stack") or {}).get("enabled"):
        st.info("This run is 2D, so there is no z-piezo. Make a 3D run (sidebar: 'Record 3D z-stacks') to see this tab.")
        return
    _piezo_section(run, R["cfg"])
    _fiducial_section(run, R)
    _effect_section(run)
    _reference()
