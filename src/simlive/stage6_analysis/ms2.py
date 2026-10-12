"""Stage 6: the MS2 transcription spot of each cell (in the nuclear channel, the same colour as the NLS-GFP), found and measured as in
Chen, Levo, Zoller & Gregor 2025 (dog_viterbi.py): median filter -> difference of Gaussians (several candidates) -> one Viterbi path per
cell, which may pick the spot up again after an OFF period of any length -> ellipsoid intensity (spot minus the ring around it, which removes
the nuclear MCP/NLS-GFP background) at the path's position in EVERY frame. Frames without a detection (OFF, or too faint) get the position
interpolated in time between the neighbouring detections (relative to the nucleus centroid), so the result is an intensity time series per
cell: the transcription trace the paper analyses.

`ms2_truth(run)` gives the simulation's true MS2 spot per cell and frame (position = the promoter, drift-corrected; photons = MS2 loops x
photons per loop x exposure; `photons_frame` after the nuclear channel's bleaching). `score_ms2` compares detections and the trace with it.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from simlive.stage6_analysis import dog_viterbi as DV
from simlive.stage6_analysis.localize import analyse_channel, centroid_in_crop, spot_params

COLUMNS = ["cell_id", "t", "detected", "crop_y_px", "crop_x_px", "lab_y_px", "lab_x_px", "snr", "intensity"]
ERODE_PX = 3       # the nucleus mask is shrunk by this much: the bright nuclear edge would otherwise give DoG candidates along the rim


def ms2_inputs(cell, erode_px: int = ERODE_PX, calib: str | None = None) -> dict:
    """The nuclear-channel crops of one cell (calibration-corrected if `calib`), its eroded mask and its centroid in crop pixels."""
    from simlive.stage6_analysis.localize import _load_calibration

    fr = cell.frames.reset_index(drop=True)
    stack = cell.channels["nucleus"].astype(np.float32)
    if calib:
        cal = _load_calibration(calib)
        stack = np.stack([cal.correct(stack[i], "nucleus", int(r.crop_y0), int(r.crop_x0)) for i, r in enumerate(fr.itertuples())])
    return {"stack": stack, "region": cell.analysis_mask(-int(erode_px)), "frames": fr, "centre": centroid_in_crop(fr, cell.crop_size)}


def trace_from_analysis(a: dict, stack: np.ndarray, centre: np.ndarray, p: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Position (crop px) and intensity in every frame: the detection where there is one, else interpolated in time (relative to the
    nucleus centroid) between detections. Returns (detected bool (n,), positions (n, 2), intensity (n,))."""
    n = len(stack)
    det = np.array([f is not None for f in a["fits"]])
    pos = np.full((n, 2), np.nan)
    for i, f in enumerate(a["fits"]):
        if f is not None:
            pos[i] = np.asarray(f["pos"][-2:]) - centre[i]
    if det.any():
        t = np.arange(n)
        for k in range(2):
            pos[:, k] = np.interp(t, t[det], pos[det, k])                  # holds the first / last value beyond the ends
    else:
        pos[:] = 0.0                                                         # no detection at all: the nucleus centre
    pos = pos + centre
    inten = np.array([DV.ellipsoid_intensity(stack[i], pos[i], p["intensity_inner_px"], p["intensity_outer_px"], p["intensity_background"])[0]
                      for i in range(n)])
    return det, pos, inten


MS2_POSITIONS = ("promoter", "own")


def ms2_trace(cell, erode_px: int = ERODE_PX, calib: str | None = None, ms2_position: str = "promoter", **kw) -> pd.DataFrame:
    """The MS2 spot of one isolated cell in every frame: whether a spot was detected, where (crop and lab pixels), its SNR and the ellipsoid
    intensity in the nuclear channel at that position (interpolated where nothing was detected). Keywords: the locus_analysis settings.
    ms2_position: "promoter" (default) = measure at the promoter (locus 0), tracked in its own colour with dog_viterbi (the MS2 spot sits on
    the promoter; the measured chromatic shift is removed when calibrated); "own" = find the MS2 spot itself in the nuclear channel with
    dog_viterbi, as Chen et al. do when no other label marks the gene."""
    if "nucleus" not in cell.channels:
        return pd.DataFrame(columns=COLUMNS)
    ci = ms2_inputs(cell, erode_px, calib)
    p = spot_params(kw)
    kw.pop("method", None)                      # the MS2 spot is always found with dog_viterbi
    if ms2_position == "promoter" and "locus0" in cell.channels:
        from simlive.stage6_analysis.localize import channel_inputs

        lc = channel_inputs(cell, "locus0", 2, calib)
        a = analyse_channel(lc["stack"], lc["region"], lc["centre"], "dog_viterbi", **kw)
        det, pos, _ = trace_from_analysis(a, lc["stack"], lc["centre"], p)
        pos = pos - np.asarray(lc["shift"], float)                            # locus-0 colour -> nuclear channel coordinates
        inten = np.array([DV.ellipsoid_intensity(ci["stack"][i], pos[i], p["intensity_inner_px"], p["intensity_outer_px"],
                                                 p["intensity_background"])[0] for i in range(len(pos))])
    else:
        a = analyse_channel(ci["stack"], ci["region"], ci["centre"], "dog_viterbi", **kw)
        det, pos, inten = trace_from_analysis(a, ci["stack"], ci["centre"], p)
    fr = ci["frames"]
    snr = [f["snr"] if f is not None else np.nan for f in a["fits"]]
    return pd.DataFrame({"cell_id": cell.cell_id, "t": fr.t.astype(int).to_numpy(), "detected": det, "crop_y_px": pos[:, 0],
                         "crop_x_px": pos[:, 1], "lab_y_px": fr.crop_y0.to_numpy() + pos[:, 0], "lab_x_px": fr.crop_x0.to_numpy() + pos[:, 1],
                         "snr": snr, "intensity": inten})[COLUMNS]


def ms2_truth(run: Path) -> pd.DataFrame | None:
    """True MS2 spot per cell and frame in image coordinates, or None if the run has no transcription model."""
    run = Path(run)
    tx_f = run / "stage1_chromatin" / "transcription.csv"
    ti = run / "stage3_microscopy" / "loci_truth_image.csv"
    if not tx_f.exists() or not ti.exists():
        return None
    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    tc = cfg["loci"]["transcription"]
    prom = pd.read_csv(ti)
    keep = ["t", "cell_id", "y_um", "x_um"] + [c for c in ("bleach_nucleus", "in_image") if c in prom.columns]
    prom = prom[prom.locus_id == 0][keep]       # the MS2 spot is drawn at the promoter, in the nuclear channel
    tx = pd.read_csv(tx_f)[["t", "cell_id", "ms2_loops", "promoter_on"]]
    m = prom.merge(tx, on=["t", "cell_id"])
    m["photons"] = float(tc.get("photons_per_loop_s", 100.0)) * float(cfg["acquisition"]["exposure_s"]) * m.ms2_loops
    it = run / "stage3_microscopy" / "imaging_truth.csv"
    if "bleach_nucleus" in m.columns:          # per-cell bleaching: this cell's own factor
        m["photons_frame"] = m.photons * m.bleach_nucleus.fillna(1.0)
    elif it.exists():
        b = pd.read_csv(it).set_index("t")
        m["photons_frame"] = m.photons * (b["bleach_nucleus"].reindex(m.t).to_numpy() if "bleach_nucleus" in b else 1.0)
    else:
        m["photons_frame"] = m.photons
    return m


def score_ms2(trace: pd.DataFrame, truth: pd.DataFrame, px_um: float, far_um: float = 0.5, k_call: float = 3.0) -> tuple[pd.DataFrame, dict]:
    """Detections and trace against the truth. A frame is ON when the true spot has any light (an MS2 loop on a Pol II).
    found_on = share of ON frames with a detection within far_um of the true spot; detected_on = with any detection; false_on_off = share of OFF frames with one; err = distance of detections to the true
    spot in ON frames (frac_far > far_um); intensity_r = Pearson r of the trace with the true photons over ALL frames; ON/OFF call from
    the trace alone: intensity > k_call x robust noise of the trace (noise from its median absolute deviation, dominated by OFF frames)."""
    if "in_image" in truth.columns:            # a promoter outside the image (its cell half out of view) cannot be measured
        truth = truth[truth.in_image.astype(bool)]
    j = trace.merge(truth, on=["cell_id", "t"], how="inner")
    if j.empty:
        return j, {}
    j["on"] = j.photons > 0
    j["err_um"] = np.where(j.detected, np.hypot(j.lab_y_px * px_um - j.y_um, j.lab_x_px * px_um - j.x_um), np.nan)
    on, off = j[j.on], j[~j.on]
    dj = on[on.detected]
    out = {"n_frames": len(j), "n_on": len(on), "found_on": float((on.err_um <= far_um).mean()) if len(on) else np.nan,
           "detected_on": float(on.detected.mean()) if len(on) else np.nan,
           "false_on_off": float(off.detected.mean()) if len(off) else np.nan,
           "median_error_um": float(dj.err_um.median()) if len(dj) else np.nan,
           "frac_far": float((dj.err_um > far_um).mean()) if len(dj) else np.nan}
    g = j[np.isfinite(j.intensity)]
    if len(g) > 5 and g.photons_frame.std() > 0:
        out["intensity_r"] = float(np.corrcoef(g.intensity, g.photons_frame)[0, 1])
    med = float(np.median(g.intensity))
    noise = 1.4826 * float(np.median(np.abs(g.intensity - med))) + 1e-9
    called = g.intensity > med + k_call * noise
    out.update(call_threshold=med + k_call * noise, call_sensitivity=float(called[g.on].mean()) if g.on.any() else np.nan,
               call_specificity=float((~called[~g.on]).mean()) if (~g.on).any() else np.nan)
    j["called_on"] = j.intensity > out["call_threshold"]
    return j, out


def analyse_run_ms2(run: Path, n_workers: int | None = None, progress=print) -> dict:
    """MS2 traces of every cell of a finished run (true cell IDs; raw and, if the run has one, calibrated), saved next to the locus results
    (stage6_analysis/truth/ms2_traces[_cal].csv) and scored into stage6_analysis/scores.json (keys ms2, ms2_calibrated)."""
    import json

    from simlive.stage5_linking.celldata import map_cells

    run = Path(run)
    truth = ms2_truth(run)
    if truth is None:
        return {}
    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    ana = {k: v for k, v in (cfg.get("locus_analysis") or {}).items() if k not in ("dilate_px", "method")}
    px_um = cfg["optics"]["pixel_size_nm"] / 1000.0
    out6 = run / "stage6_analysis" / "truth"
    out6.mkdir(parents=True, exist_ok=True)
    cal = run / "stage3_microscopy" / "calibration"
    res = {}
    for tag, cdir in (("", None),) + ((("_cal", str(cal)),) if (cal / "calibration.json").exists() else ()):
        tr = map_cells(ms2_trace, run / "stage5_cells" / "truth", n_workers=n_workers, calib=cdir, **ana)
        tr.to_csv(out6 / f"ms2_traces{tag}.csv", index=False)
        res["ms2" + ("_calibrated" if tag else "")] = sc = score_ms2(tr, truth, px_um)[1]
        progress(f"MS2{tag}: ON frames found {sc.get('found_on', 0):.3f} (any detection {sc.get('detected_on', 0):.3f}), OFF frames with a detection {sc.get('false_on_off', 0):.3f}, "
                 f"median error {1000 * (sc.get('median_error_um') or 0):.0f} nm, trace r {sc.get('intensity_r', float('nan')):.2f}")
    sj = run / "stage6_analysis" / "scores.json"
    scores = json.loads(sj.read_text()) if sj.exists() else {}
    scores.update(res)
    sj.write_text(json.dumps(scores, indent=2, default=str))
    return res
