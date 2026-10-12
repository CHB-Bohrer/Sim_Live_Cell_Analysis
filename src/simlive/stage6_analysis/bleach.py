"""Stage 6: photobleaching correction, with the background taken out first.

A recorded pixel is  offset + background + signal x b(t):  only the fluorophores bleach (b(t) falls from 1); the camera offset, stray light
and autofluorescence do not (or not with the same time course). Two consequences:

1. Estimating b(t) from raw brightness is biased: the constant background dilutes the decay, so bleaching looks weaker than it is (the
   `naive` estimate below). The background must be subtracted before the decay is measured.
2. Correcting means dividing the background-SUBTRACTED signal by b(t). Spot intensities (ellipse minus the ring around it) are already
   background-subtracted, so intensity / b(t) is the corrected value; dividing raw pixel values would also inflate the background.

Estimates (both from the data alone; the simulation's true b(t) is only used for scoring):
- nuclear channel (NLS-GFP; also the MCP-GFP of the MS2 spots): per frame, mean inside the nuclei minus the mean of the empty area around
  them (pixels more than `gap_px` from any nucleus), then an exponential fit through the frames.
- locus colours: the spots' own background-subtracted intensities, median over all cells per frame (probe binding makes single spots
  fluctuate, the median over cells does not decay for any other reason), then an exponential fit.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_dilation


def fit_exponential(t_s: np.ndarray, s: np.ndarray) -> tuple[float, float]:
    """Least-squares fit of log(s) = log(s0) - t / tau over the frames with s > 0. Returns (s0, tau_s); tau = inf if nothing decays."""
    t_s, s = np.asarray(t_s, float), np.asarray(s, float)
    ok = np.isfinite(s) & (s > 0)
    if ok.sum() < 3:
        return float(np.nanmean(s)) if np.isfinite(s).any() else np.nan, np.inf
    k, c = np.polyfit(t_s[ok], np.log(s[ok]), 1)
    return float(np.exp(c)), (float(-1.0 / k) if k < 0 else np.inf)


def nuclear_bleach(images: str | Path, labels: str | Path, dt_s: float, gap_px: int = 8) -> pd.DataFrame:
    """Per frame: mean in the nuclei, mean of the empty background, their difference, and the bleaching factor from an exponential fit of the
    difference (`b_est`) and, for comparison, of the raw nuclear mean (`b_naive`, background left in). One frame in memory at a time."""
    img, lab = tifffile.memmap(str(images), mode="r"), tifffile.memmap(str(labels), mode="r")
    rows = []
    st = np.ones((2 * gap_px + 1, 2 * gap_px + 1), bool)
    for t in range(len(img)):
        f, m = np.asarray(img[t], np.float32), np.asarray(lab[t]) > 0
        far = ~binary_dilation(m, structure=st)
        rows.append((t, float(f[m].mean()) if m.any() else np.nan, float(f[far].mean()) if far.any() else np.nan))
    d = pd.DataFrame(rows, columns=["t", "inside", "background"])
    d["signal"] = d.inside - d.background
    t_s = d.t.to_numpy() * dt_s
    _, tau = fit_exponential(t_s, d.signal.to_numpy())
    _, tau_naive = fit_exponential(t_s, d.inside.to_numpy())
    d["b_est"], d["b_naive"] = np.exp(-t_s / tau), np.exp(-t_s / tau_naive)
    d.attrs.update(tau_s=tau, tau_naive_s=tau_naive)
    return d


def spot_bleach(intensity: pd.DataFrame, dt_s: float, n_frames: int) -> tuple[np.ndarray, float]:
    """Bleaching factor of one spot colour from its detections (columns t, intensity): exponential fit of the per-frame median over cells.
    Returns (b(t) for t = 0..n_frames-1, tau_s)."""
    med = intensity.groupby("t").intensity.median()
    _, tau = fit_exponential(med.index.to_numpy() * dt_s, med.to_numpy())
    return np.exp(-np.arange(n_frames) * dt_s / tau), tau


def correct(df: pd.DataFrame, b: np.ndarray, col: str = "intensity") -> pd.Series:
    """Background-subtracted intensity / b(t) (the background is not divided: it was already removed by the ring around the spot)."""
    return df[col] / b[np.clip(df.t.to_numpy().astype(int), 0, len(b) - 1)]


def bleach_correct_run(run: Path, sources: tuple = ("truth",), locus_flat: str = "nucleus") -> dict:
    """Estimate b(t) for every channel of a finished run, add `intensity_corr` to its locus and MS2 tables, and return a summary with the
    estimates next to the simulation's true bleaching (imaging_truth.csv). Writes stage6_analysis/bleach.csv."""
    import yaml

    run = Path(run)
    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    dt = float(cfg["acquisition"]["frame_interval_s"])
    s3, s6 = run / "stage3_microscopy", run / "stage6_analysis"
    T = int(cfg["acquisition"]["n_frames"])
    nuc = nuclear_bleach(s3 / "nucleus.tif", s3 / "labels.tif", dt)
    out = pd.DataFrame({"t": np.arange(T)})
    out["nucleus_est"], out["nucleus_naive"] = nuc.b_est.to_numpy()[:T], nuc.b_naive.to_numpy()[:T]
    summary = {"nucleus": {"tau_est_s": nuc.attrs["tau_s"], "tau_naive_s": nuc.attrs["tau_naive_s"]}}
    for src in sources:
        d6 = s6 / src
        for name in ("loci_positions.csv", "loci_positions_cal.csv"):
            f = d6 / name
            if not f.exists():
                continue
            pos = pd.read_csv(f)
            if "intensity" not in pos:
                continue
            pos["intensity_corr"] = np.nan
            for k, g in pos.groupby("locus_id"):
                b, tau = spot_bleach(g, dt, T)
                pos.loc[g.index, "intensity_corr"] = correct(g, b)
                if src == "truth" and name == "loci_positions.csv":
                    out[f"locus{int(k)}_est"] = b
                    summary[f"locus{int(k)}"] = {"tau_est_s": tau}
            pos.to_csv(f, index=False)
        for name in ("ms2_traces.csv", "ms2_traces_cal.csv"):
            f = d6 / name
            if f.exists():
                tr = pd.read_csv(f)
                tr["intensity_corr"] = correct(tr, out.nucleus_est.to_numpy())    # MCP-GFP bleaches like the NLS-GFP of the nucleus
                tr.to_csv(f, index=False)
    it = s3 / "imaging_truth.csv"
    if it.exists():
        tru = pd.read_csv(it)
        for c in ("nucleus", "locus0", "locus1"):
            if f"bleach_{c}" in tru:
                out[f"{c}_true"] = tru[f"bleach_{c}"].to_numpy()[:T]
                tau_true = cfg["optics"].get("bleach_tau_s") if c == "nucleus" else (cfg.get("loci") or {}).get("bleach_tau_s")
                if c in summary:
                    summary[c]["tau_true_s"] = tau_true
                    summary[c]["max_error"] = float(np.nanmax(np.abs(out[f"{c}_est"] - out[f"{c}_true"])))
        if "nucleus_true" in out:
            summary["nucleus"]["max_error_naive"] = float(np.nanmax(np.abs(out.nucleus_naive - out.nucleus_true)))
    s6.mkdir(parents=True, exist_ok=True)
    out.to_csv(s6 / "bleach.csv", index=False)
    for src in sources:        # the global correction above is kept as intensity_corr_global; intensity_corr becomes flat-field + per cell
        try:
            ff = flat_bleach_correct(run, src, cfg, locus_flat)
        except FileNotFoundError:
            continue
        if src == "truth":
            summary["flatfield"] = ff
    summary.update(_scores(run, s6))
    return summary


LABELS = {"truth": ("stage3_microscopy", "labels.tif"), "tracked_gt_masks": ("stage4_segtrack", "tracked_gt_masks.tif"),
          "tracked_auto_masks": ("stage4_segtrack", "tracked_auto_masks.tif")}


def flat_bleach_correct(run: Path, src: str, cfg: dict, locus_flat: str = "nucleus") -> dict:
    """Uneven illumination + per-cell bleaching from the movie itself (flatfield.py), for one identity source: fit on the nuclei (the nuclear
    channel, also the MS2 spots) and on each locus colour's spot intensities, then intensity_corr = intensity / (F(spot) x b_cell(t)) in the
    locus and MS2 tables (tables made with the calibration applied are already flat-field corrected: only b_cell(t) is divided out). The
    previous global correction moves to intensity_corr_global. Returns the fits and, for the truth source, how close they are to the truth."""
    from simlive.stage6_analysis.flatfield import correct_values, fit_flat_bleach, nuclear_signals

    run = Path(run)
    s3, d6 = run / "stage3_microscopy", run / "stage6_analysis" / src
    lab = run / LABELS[src][0] / LABELS[src][1]
    if not lab.exists() or not d6.exists():
        raise FileNotFoundError(lab)
    dt = float(cfg["acquisition"]["frame_interval_s"])
    nuc = nuclear_signals(s3 / "nucleus.tif", lab)
    shape = tifffile.TiffFile(str(lab)).pages[0].shape[-2:]
    track = nuc[["cell_id", "t", "y_px", "x_px"]]
    fits = {"nucleus": fit_flat_bleach(nuc[~nuc.touches_border].rename(columns={"signal": "value"}), track, shape, dt)}
    nuc["b_est"] = fits["nucleus"].bleach(track).to_numpy()
    nuc["flat_est"] = fits["nucleus"].F(nuc.y_px.to_numpy(), nuc.x_px.to_numpy())
    raw_tables = {}
    for name in ("loci_positions.csv", "loci_positions_cal.csv"):
        f = d6 / name
        if not f.exists():
            continue
        pos = pd.read_csv(f)
        if "intensity" not in pos:
            continue
        if "intensity_corr" in pos:
            pos["intensity_corr_global"] = pos.intensity_corr
        cal = name.endswith("_cal.csv")
        for k, g in pos.groupby("locus_id"):
            ch = f"locus{int(k)}"
            if not cal or ch not in fits:
                vals = g.rename(columns={"lab_y_px": "y_px", "lab_x_px": "x_px", "intensity": "value"})
                if locus_flat == "own" and not cal:      # the colour's own illumination, fitted on its spot intensities
                    fits[ch] = fit_flat_bleach(vals, track, shape, dt, fit_flat=True)
                else:      # "nucleus": the nuclear channel's illumination shape (same light path), only the bleaching rate fitted; "none": even light
                    fits[ch] = fit_flat_bleach(vals, track, shape, dt, fit_flat=False, dose_flat=fits["nucleus"] if locus_flat == "nucleus" else None)
            pos.loc[g.index, "intensity_corr"] = correct_values(g, fits[ch], track, divide_flat=not cal).to_numpy()
        pos.to_csv(f, index=False)
        raw_tables[name] = pos
    for name in ("ms2_traces.csv", "ms2_traces_cal.csv"):
        f = d6 / name
        if f.exists():
            tr = pd.read_csv(f)
            if "intensity_corr" in tr:
                tr["intensity_corr_global"] = tr.intensity_corr
            tr["intensity_corr"] = correct_values(tr, fits["nucleus"], track, divide_flat=not name.endswith("_cal.csv")).to_numpy()
            tr.to_csv(f, index=False)
    nuc.to_csv(d6 / "nuclear_signals.csv", index=False)
    np.savez_compressed(d6 / "flatfield_est.npz", **{c: v.map() for c, v in fits.items()})
    summary = {c: v.summary() for c, v in fits.items()}
    for c in summary:
        if c != "nucleus" and locus_flat != "own":
            summary[c]["flat_source"] = {"nucleus": "the nuclear channel's estimate", "none": "even light"}[locus_flat]
    # a nucleus' NLS-GFP brightness is constant in the simulation: what spread over time is left in each cell after each correction?
    ok = (~nuc.touches_border) & (nuc.signal > 0)
    n = nuc[ok]
    bg = pd.read_csv(run / "stage6_analysis" / "bleach.csv").set_index("t").nucleus_est if (run / "stage6_analysis" / "bleach.csv").exists() else None
    spread = lambda v: float(np.nanmedian(pd.Series(np.log(v.to_numpy()), index=n.index).groupby(n.cell_id).std()))  # noqa: E731
    summary["nucleus"]["within_cell_sd_raw"] = spread(n.signal)
    if bg is not None:
        summary["nucleus"]["within_cell_sd_global"] = spread(n.signal / bg.reindex(n.t).to_numpy())
    summary["nucleus"]["within_cell_sd_corrected"] = spread(n.signal / (n.flat_est * n.b_est))
    if src == "truth":
        _flat_truth(run, cfg, fits, nuc, summary)
    return summary


def _flat_truth(run: Path, cfg: dict, fits: dict, nuc: pd.DataFrame, summary: dict) -> None:
    """Compare the estimated illumination maps and per-cell bleaching with the simulation's truth (and with the calibration's flat-field)."""
    from simlive.stage3_microscopy.render import imager_for_run

    s3 = run / "stage3_microscopy"
    try:
        im = imager_for_run(run, cfg)
    except Exception:
        return
    cal = None
    if (s3 / "calibration" / "calibration.json").exists():
        from simlive.stage3_microscopy.calibration.model import Calibration

        cal = Calibration.load(s3 / "calibration")
    rms = lambda a, b: float(np.sqrt(np.mean((a / a.mean() - b / b.mean()) ** 2)))  # noqa: E731
    for c, fit in fits.items():
        true = im._flat(c, im.E)
        true = np.ones(im.shape, np.float32) if true is None else true
        s = summary[c]
        s["flat_rms_error_est"] = rms(fit.map(), true)              # the map used (even light if the fit was rejected)
        s["flat_rms_error_none"] = rms(np.ones_like(true), true)              # error of assuming even illumination
        if cal is not None and c in cal.maps:
            s["flat_rms_error_calibration"] = rms(cal.maps[c]["flat"].astype(np.float64), true)
        s["tau_true_s"] = im._tau(c)
    bc = s3 / "bleach_cells.csv"
    if bc.exists():                         # per-cell truth: how far is each cell's estimated bleaching from its true one?
        b = pd.read_csv(bc).merge(nuc[["t", "cell_id", "b_est"]], on=["t", "cell_id"])
        glob = pd.read_csv(run / "stage6_analysis" / "bleach.csv").set_index("t")
        b["b_global"] = glob.nucleus_est.reindex(b.t).to_numpy()
        summary["nucleus"]["per_cell_b_rms_error_est"] = float(np.sqrt(np.mean((b.b_est - b.bleach_nucleus) ** 2)))
        summary["nucleus"]["per_cell_b_rms_error_global"] = float(np.sqrt(np.mean((b.b_global - b.bleach_nucleus) ** 2)))
        summary["nucleus"]["per_cell_b_true_range"] = [float(b.bleach_nucleus.min()), float(b.bleach_nucleus.max())]


def _scores(run: Path, s6: Path) -> dict:
    """Does the correction help? Correlation of raw and corrected intensities with the photons the simulation drew BEFORE bleaching, and the
    late / early ratio of their means (1 = no trend left; raw intensities fall with the bleaching). Correct detections only (<= 0.5 um)."""
    out = {}
    ti, f = run / "stage3_microscopy" / "loci_truth_image.csv", s6 / "truth" / "loci_positions_cal.csv"
    f = f if f.exists() else s6 / "truth" / "loci_positions.csv"
    if ti.exists() and f.exists():
        tr = pd.read_csv(ti)
        pos = pd.read_csv(f)
        if "photons" in tr and "intensity_corr" in pos:
            import yaml

            px = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())["optics"]["pixel_size_nm"] / 1000.0
            j = pos.merge(tr, on=["cell_id", "t", "locus_id"])
            j = j[np.hypot(j.lab_y_px * px - j.y_um, j.lab_x_px * px - j.x_um) <= 0.5]
            for k, g in j.groupby("locus_id"):
                out[f"locus{int(k)}_r"] = _r(g, "photons")
    for name in ("ms2_traces_cal.csv", "ms2_traces.csv"):
        f = s6 / "truth" / name
        if f.exists():
            from simlive.stage6_analysis.ms2 import ms2_truth

            mt = ms2_truth(run)
            if mt is not None:
                tr = pd.read_csv(f)
                if "intensity_corr" in tr:
                    out["ms2_r"] = _r(tr.merge(mt, on=["cell_id", "t"]), "photons")
            break
    return out


def _r(g: pd.DataFrame, truth_col: str) -> dict:
    g = g[np.isfinite(g.intensity) & np.isfinite(g.intensity_corr)]
    if len(g) < 10 or g[truth_col].std() == 0:
        return {}
    T = g.t.max()
    early, late = g[g.t < T / 4], g[g.t >= 3 * T / 4]
    ratio = lambda c: float(late[c].mean() / early[c].mean()) if len(early) and len(late) and early[c].mean() != 0 else np.nan  # noqa: E731
    out = {"r_raw": float(np.corrcoef(g.intensity, g[truth_col])[0, 1]), "r_corrected": float(np.corrcoef(g.intensity_corr, g[truth_col])[0, 1]),
           "late_early_raw": ratio("intensity"), "late_early_corrected": ratio("intensity_corr"), "late_early_true": ratio(truth_col)}
    if "intensity_corr_global" in g:     # the earlier correction: one b(t) for the whole field, illumination ignored
        gg = g[np.isfinite(g.intensity_corr_global)]
        out["r_corrected_global"] = float(np.corrcoef(gg.intensity_corr_global, gg[truth_col])[0, 1]) if len(gg) > 10 else np.nan
        out["late_early_corrected_global"] = ratio("intensity_corr_global")
    return out
