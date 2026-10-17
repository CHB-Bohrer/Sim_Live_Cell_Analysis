"""Stage 7: what each source of image error does to the RESULTS, and how much each correction gives back.

For one finished 2D run with loci, the movie is rendered again from the same ground truth (same cells, loci, MS2 spots, same noise streams):

- `reference`: a perfect microscope that still obeys physics: diffraction blur, quantum efficiency, photon shot noise, read noise. Nothing
  here can be corrected; it sets the best possible result.
- `+<source>`: the reference plus ONE error source (e.g. `+hot_pixels`).
- `all`: every source the run has (the realistic case).

Every movie is cut into single-cell crops (true cell IDs) and analysed several times: with the correction for that source OFF, and ON in
each available way. The ways are:
- the analysis' own countermeasures: median filter, shell background, Viterbi path;
- calibration images (dark frames, flat-field slide, bead slide; `stage3_microscopy/calibration`);
- the movie alone (`stage6_analysis/selfcal.py`, `flatfield.py`, the row and bleaching corrections).

Each analysis is scored against the truth:
- found, wrong object, median and 95th-percentile position error;
- bias of the promoter-enhancer distance;
- how well spot intensities follow the true photons (correlation, % error after one global scale, late/early drift);
- the MS2 trace correlation, if the run has transcription.

Output: <run>/stage7_validation/error_effects.csv (one row per movie x analysis), resumable. Work files are deleted after each movie.
`scripts/error_effects.py` runs it; the dashboard's Images page shows it per source.
"""
from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import yaml

from simlive.stage3_microscopy import errors as ERR

REFERENCE = ("psf_blur", "qe", "shot_noise", "read_noise")
ONLY_3D = ("slice_skew", "piezo", "biplane_registration")

# per source: the analyses to compare, first = not corrected. spec keys: localize_loci keywords (median_px, method, intensity_background,
# row_correct), "calib": None | "images" | "movie" with "parts" (which calibration components), "post": raw | global | percell (bleaching).
BG_NONE = {"intensity_background": "none"}
FIXES: dict[str, list[tuple[str, str, dict]]] = {
    "stage_drift": [("none", "nothing to correct here: each cell's crop follows the cell, so drift moves the cell and its loci together", {})],
    "chromatic_shift": [("none", "not corrected", {}),
                        ("movie", "shift from the movie (mean enhancer - promoter vector; used only with >= 30 cells and > 3 standard errors)",
                         {"calib": "movie", "parts": ("shift",)}),
                        ("images", "shift from the bead slide", {"calib": "images", "parts": ("shift",)})],
    "sequential_colour": [("none", "no correction possible: the loci really moved between the colour exposures (random motion)", {})],
    "defocus": [("none", "no correction exists in 2D (the Gaussian fit's width follows the blur)", {})],
    "exposure_motion": [("none", "no correction exists (the spot is wider; the fit follows it)", {})],
    "haze": [("none", "no local background (camera offset only)", BG_NONE), ("analysis", "shell background around each spot (Chen et al.)", {})],
    "bleaching": [("none", "not corrected", {"post": "raw"}), ("global", "one decay curve for the whole movie", {"post": "global"}),
                  ("movie", "illumination + each cell's own light dose, from the movie", {"post": "percell"})],
    "flat_field": [("none", "not corrected", {}), ("movie", "illumination map from the moving nuclei", {"calib": "movie", "parts": ("flat",)}),
                   ("images", "flat-field slide", {"calib": "images", "parts": ("flat",)})],
    "flicker": [("none", "not corrected (1% placeholder; no movie-only method for the locus colours yet)", {})],
    "background": [("none", "no local background (camera offset only)", BG_NONE), ("analysis", "shell background", {})],
    "autofluorescence": [("none", "no local background (camera offset only)", BG_NONE), ("analysis", "shell background", {}),
                         ("movie", "background map from the movie + shell", {"calib": "movie", "parts": ("offset",)})],
    "cytoplasm": [("none", "no local background (camera offset only)", BG_NONE), ("analysis", "shell background", {}),
                  ("movie_median", "background map from the movie: MEDIAN of nucleus-free frames (biased by cell bodies) + shell",
                   {"calib": "movie", "parts": ("offset", "flat"), "background_stat": "median"}),
                  ("movie", "background map from the movie: low percentile, noise-corrected (robust to cell bodies) + shell",
                   {"calib": "movie", "parts": ("offset", "flat")})],
    "crosstalk": [("none", "no local background (camera offset only)", BG_NONE), ("analysis", "shell background", {}),
                  ("movie", "bleed-through fraction measured in the movie and subtracted (unmixing) + shell", {"calib": "movie", "parts": ("crosstalk",)})],
    "dark_current": [("none", "no local background (camera offset only)", BG_NONE), ("analysis", "shell background", {})],
    "hot_pixels": [("none", "no median filter, no hot-pixel map", {"median_px": 1}), ("analysis", "median filter 3 x 3 (default)", {}),
                   ("movie", "hot-pixel map from the movie (no median)", {"median_px": 1, "calib": "movie", "parts": ("hot",)}),
                   ("images", "hot-pixel map from dark frames (no median)", {"median_px": 1, "calib": "images", "parts": ("hot",)})],
    "cosmic_rays": [("none", "brightest spot in each frame, frame by frame", {"method": "brightest"}),
                    ("nomedian", "DoG candidates + Viterbi path, no median filter", {"median_px": 1}),
                    ("analysis", "median filter + DoG + Viterbi path (default)", {})],
    "prnu": [("none", "not corrected", {}), ("images", "pixel gains from the flat-field slide", {"calib": "images", "parts": ("prnu",)})],
    "row_noise": [("none", "no row correction", {}), ("movie", "row offsets measured in each crop", {"row_correct": True})],
    "dsnu": [("none", "no local background (camera offset only)", BG_NONE), ("analysis", "shell background", {}),
             ("movie", "offset map from the movie + shell", {"calib": "movie", "parts": ("offset",)}),
             ("images", "dark frames + shell", {"calib": "images", "parts": ("offset",)})],
    "adc": [("none", "nothing to correct (rounding is far below the noise; nothing saturates here)", {})],
}
ALL_FIXES = [("none", "nothing corrected (no median, no shell background, no calibration, no bleaching correction)",
              {"median_px": 1, **BG_NONE, "post": "raw"}),
             ("analysis", "the analysis' own countermeasures only (median, shell background, Viterbi path)", {"post": "raw"}),
             ("movie", "+ everything estimated from the movie (background / hot-pixel / illumination maps, colour shift, cross-talk, rows, "
                       "per-cell bleaching)",
              {"calib": "movie", "parts": ("offset", "hot", "flat", "shift", "crosstalk"), "row_correct": True, "post": "percell"}),
             ("images", "+ calibration images instead (dark, flat-field and bead slides) + rows + per-cell bleaching",
              {"calib": "images", "parts": ("offset", "hot", "flat", "prnu", "shift"), "row_correct": True, "post": "percell"}),
             ("both", "calibration images + what only the movie gives (cross-talk unmixing, rows, per-cell bleaching)",
              {"calib": "images", "parts": ("offset", "hot", "flat", "prnu", "shift"), "movie_parts": ("crosstalk",), "row_correct": True,
               "post": "percell"})]
# how each source is (or is not) corrected, in plain words: by the analysis itself, with calibration images, from the movie alone
HOW = {
    "stage_drift": dict(analysis="Nothing needed for per-cell results: each cell's crop is centred on that cell in every frame, so drift moves "
                                 "the cell and its loci together and cancels in positions inside the nucleus and in distances.",
                        images="Fiducial beads on the slide measure the drift of every frame (📍 Stage & fiducials page; z-stacks).",
                        movie="The cells themselves are the reference here; for absolute stage positions one would register frames on the "
                              "whole nuclear image (not built: no result here needs it)."),
    "chromatic_shift": dict(analysis="Not correctable by the spot finder: each colour is shifted by a constant amount.",
                            images="A bead slide (beads glow in every colour) gives each colour's offset; it is subtracted from every position.",
                            movie="The average enhancer - promoter vector over all cells and frames: true orientations are random, so the "
                                  "average is the shift between the two colours. Used only with >= 30 cells and when larger than 3 standard errors: it needs many "
                                  "cells, because with few cells the true vectors do not average out."),
    "sequential_colour": dict(analysis="Cannot be corrected afterwards: between the colour exposures the chromatin moves randomly.",
                              images="(Nothing to calibrate.)",
                              movie="Only the acquisition can reduce it: shorter exposures / faster switching, the two locus colours next to each "
                                    "other in the order, or simultaneous two-colour detection (a second camera), which brings back some bleed-through."),
    "psf_blur": dict(analysis="Cannot be removed, only modelled: the difference of Gaussians is matched to the blur and the Gaussian fit finds the "
                              "spot's centre to a small fraction of its width. With the photon noise it sets the best possible precision.",
                     images="A bead slide measures the blur (used to check the PSF width).", movie=""),
    "defocus": dict(analysis="Not corrected: the Gaussian fit's width follows the blur, which costs precision and some brightness.",
                    images="(Autofocus hardware or a z-stack avoid it.)", movie="(none)"),
    "exposure_motion": dict(analysis="Not corrected: the spot is slightly wider than the PSF; the fit follows it.", images="(Shorter exposures.)",
                            movie="(none)"),
    "haze": dict(analysis="The shell around each spot (Chen et al.) measures the local glow and subtracts it from the spot.", images="(none needed)",
                 movie="(none needed)"),
    "bleaching": dict(analysis="Spot intensities are divided by the bleaching factor, measured on background-subtracted signal.",
                      images="(none needed)",
                      movie="Each cell's own bleaching from the light it received along its path (illumination map x dwell), fitted on the "
                            "moving nuclei; cells entering the field start unbleached (Spot finding page, 3b-3c)."),
    "flat_field": dict(analysis="Not by itself.", images="A flat-field slide (uniformly fluorescent) gives the illumination map; images are divided by it.",
                       movie="The same nucleus seen at different places: brightness changes that follow position are the illumination. "
                             "Kept only if it predicts held-out cells better than even light (needs cells that move across the field)."),
    "flicker": dict(analysis="Not corrected (1% placeholder).", images="(A power meter in the light path would log it.)",
                    movie="For the nuclear channel the common brightness change of all nuclei in a frame would give it (not built); "
                          "spot brightness fluctuates too much to see 1% flicker."),
    "background": dict(analysis="The shell around each spot measures the local background and subtracts it.", images="Dark frames (camera part only).",
                       movie="Background map from the frames in which no nucleus covers a pixel."),
    "autofluorescence": dict(analysis="The shell around each spot subtracts the local level.", images="(Not in dark frames: it is sample light.)",
                             movie="Background map: per-pixel median over the frames free of nuclei (cells move, so most pixels are free some of the time)."),
    "cytoplasm": dict(analysis="The shell around each spot subtracts the local level, cytoplasm included (spots sit inside the nucleus).",
                      images="(Not in dark frames or a flat slide: it is sample light.)",
                      movie="Careful: 'pixels away from the nuclei' are not empty, a cell body may be there. The movie's background map uses the "
                            "LOW end of each pixel's nucleus-free frames (20th percentile, corrected for the pixel's own noise), i.e. the level the "
                            "pixel returns to when no cell is there, instead of the median, which a passing cell body raises. The illumination "
                            "fit subtracts each nucleus' local surroundings, so the cytoplasm there is removed with the background."),
    "crosstalk": dict(analysis="The shell subtracts the smooth ghost of the nucleus under each spot, but NOT a bright MS2 spot leaking into the "
                                 "promoter's colour: it sits exactly on the promoter.",
                      images="Single-colour samples would give the bleed-through fractions (not simulated here).",
                      movie="Inside the nuclei a locus colour rises with the nuclear colour by the bleed-through fraction; it is measured by "
                            "robust regression and fraction x the nuclear crop is subtracted from the locus crop (unmixing)."),
    "qe": dict(analysis="Nothing to correct: a constant scale; relative intensities are unaffected.", images="", movie=""),
    "shot_noise": dict(analysis="Cannot be corrected: it is the randomness of light itself. It sets the best possible precision.", images="", movie=""),
    "dark_current": dict(analysis="The shell subtracts the extra level; its own shot noise cannot be removed.",
                         images="Dark frames at the movie's exposure give its mean (uniform, so the shell already removes it).",
                         movie="Background map from the movie."),
    "hot_pixels": dict(analysis="A 3 x 3 median filter before spot finding removes single bright pixels; the Viterbi path ignores what does not move "
                                "like a locus.",
                       images="Dark frames: pixels high in every dark frame are flagged and replaced by their neighbours' median.",
                       movie="Pixels high in every frame in which no nucleus covers them (a statistical test per pixel); hot pixels under "
                             "nuclei that never move away cannot be found this way."),
    "cosmic_rays": dict(analysis="The median filter shrinks them and the Viterbi path rejects a spot that is there for one frame and far from "
                                 "the locus' path; a frame-by-frame 'brightest spot' picks them.", images="(Random in time: nothing to calibrate.)",
                        movie="(Their one-frame nature is what the path uses.)"),
    "prnu": dict(analysis="Not by itself.", images="The flat-field slide also holds the pixel-to-pixel gain pattern.",
                 movie="Not from the movie: the illumination fit is smooth and cannot see a 1% speckle."),
    "read_noise": dict(analysis="Cannot be corrected: electronic noise of every read. With the shot noise it sets the precision.", images="", movie=""),
    "row_noise": dict(analysis="The shell subtracts part of it.", images="Dark frames show it but it changes every frame, so they cannot correct it.",
                      movie="Each row's offset measured in every frame from the pixels away from the cell, and subtracted."),
    "dsnu": dict(analysis="The shell subtracts the local level.", images="Dark frames give each pixel's offset.",
                 movie="Background map from the frames free of nuclei holds the offset pattern too."),
    "adc": dict(analysis="Nothing to correct: rounding is far below the noise and nothing saturates at these brightnesses.", images="", movie=""),
}
assert set(HOW) == set(ERR.SOURCES) - set(ONLY_3D)
# what each source mainly damages (the dashboard highlights these columns)
KEY_METRIC = {"sequential_colour": "dist_rms_nm", "cytoplasm": "intensity_err_pct", "chromatic_shift": "dist_rms_nm", "hot_pixels": "frac_far", "cosmic_rays": "frac_far", "bleaching": "late_early",
              "flat_field": "intensity_err_pct", "prnu": "intensity_err_pct", "flicker": "intensity_err_pct", "haze": "intensity_err_pct",
              "background": "intensity_err_pct", "autofluorescence": "intensity_err_pct", "crosstalk": "intensity_err_pct",
              "dark_current": "intensity_err_pct", "dsnu": "intensity_err_pct", "row_noise": "median_err_nm", "defocus": "median_err_nm",
              "exposure_motion": "median_err_nm", "stage_drift": "median_err_nm", "adc": "median_err_nm"}
METRICS = ["found", "frac_far", "median_err_nm", "p95_err_nm", "dist_bias_nm", "dist_rms_nm", "intensity_r", "intensity_err_pct",
           "late_early", "ms2_r"]


def variants(E: ERR.ErrorSet) -> list[tuple[str, str, tuple]]:
    """(name, source, enabled sources) for the reference, every source the run has on (2D sources) and all of them."""
    out = [("reference", "", REFERENCE)]
    for s in ERR.ORDER:
        if s in REFERENCE or s in ONLY_3D or not E.on(s):
            continue
        out.append((f"+{s}", s, tuple(sorted(set(REFERENCE) | {s}))))
    out.append(("all", "all", tuple(E.enabled_ids())))
    return out


def render_variant(im, E, s3: Path, progress=None) -> None:
    """Render every frame of every channel with exactly the sources of E and write <channel>.tif + labels.tif (2D, uint16)."""
    s3.mkdir(parents=True, exist_ok=True)
    out = {c: np.zeros((im.T, *im.shape), np.uint16) for c in im.channels}
    lab = np.zeros((im.T, *im.shape), np.uint16)
    for t in range(im.T):
        f = im.frame(t, E, im.channels)
        for c in im.channels:
            out[c][t] = np.rint(f["images"][c]).clip(0, 65535)
        if f["labels"] is not None:
            lab[t] = f["labels"]
        if progress and t % 60 == 0:
            progress(f"    frame {t + 1}/{im.T}")
    for c, a in out.items():
        tifffile.imwrite(s3 / f"{c}.tif", a, photometric="minisblack")
    tifffile.imwrite(s3 / "labels.tif", lab, photometric="minisblack")


def _distance_errors(j: pd.DataFrame, px_um: float) -> tuple[float, float]:
    """Bias and rms (nm) of the measured promoter-enhancer distance, frames where both were found within 0.5 um of the truth."""
    a = j[j.locus_id == 0].set_index(["cell_id", "t"])
    b = j[j.locus_id == 1].set_index(["cell_id", "t"])
    k = a.index.intersection(b.index)
    if len(k) < 5:
        return np.nan, np.nan
    a, b = a.loc[k], b.loc[k]
    dm = px_um * np.hypot(a.lab_y_px - b.lab_y_px, a.lab_x_px - b.lab_x_px)
    dt = np.hypot(a.y_um - b.y_um, a.x_um - b.x_um)
    e = 1000 * (dm - dt)
    return float(e.mean()), float(np.sqrt((e ** 2).mean()))


def _intensity_metrics(I: np.ndarray, P: np.ndarray, t: np.ndarray) -> dict:
    """How well intensities follow the true photons: Pearson r; rms error (% of the mean) after the best single scale factor through zero;
    late/early = (last quarter / first quarter of the movie) of the intensities divided by the same ratio of the true photons (1 = no drift)."""
    ok = np.isfinite(I) & np.isfinite(P) & (P > 0)
    I, P, t = I[ok], P[ok], t[ok]
    if len(I) < 20 or np.std(P) == 0:
        return {"intensity_r": np.nan, "intensity_err_pct": np.nan, "late_early": np.nan}
    a = float((I * P).sum() / (P * P).sum())
    err = 100 * float(np.sqrt(np.mean((I - a * P) ** 2)) / max(np.mean(a * P), 1e-9))
    T = t.max()
    e, l = t < T / 4, t >= 3 * T / 4
    le = (I[l].mean() / I[e].mean()) / (P[l].mean() / P[e].mean()) if e.any() and l.any() and I[e].mean() != 0 else np.nan
    return {"intensity_r": float(np.corrcoef(I, P)[0, 1]), "intensity_err_pct": err, "late_early": float(le)}


def score(pos: pd.DataFrame, truth: pd.DataFrame, px_um: float, col: str = "intensity", far_um: float = 0.5) -> dict:
    from simlive.stage6_analysis.localize import _visible

    tr = truth.assign(visible=_visible(truth))
    n_vis = int(tr.visible.sum())
    if not len(pos):
        return {"found": 0.0}
    j = pos.merge(tr, on=["cell_id", "t", "locus_id"], how="inner")
    j["err_um"] = np.hypot(j.lab_y_px * px_um - j.y_um, j.lab_x_px * px_um - j.x_um)
    vis = j[j.visible]
    near = vis[vis.err_um <= far_um]
    out = {"found": len(near) / max(n_vis, 1), "frac_far": float((vis.err_um > far_um).mean()) if len(vis) else np.nan,
           "median_err_nm": 1000 * float(near.err_um.median()) if len(near) else np.nan,
           "p95_err_nm": 1000 * float(near.err_um.quantile(0.95)) if len(near) else np.nan}
    out["dist_bias_nm"], out["dist_rms_nm"] = _distance_errors(near, px_um)
    ims = [_intensity_metrics(g[col].to_numpy(float), g.photons.to_numpy(float), g.t.to_numpy()) for _, g in near.groupby("locus_id")] \
        if "photons" in near else []
    for k in ("intensity_r", "intensity_err_pct", "late_early"):
        v = [m[k] for m in ims if np.isfinite(m[k])]
        out[k] = float(np.mean(v)) if v else np.nan
    return out


def _ms2_r(tr: pd.DataFrame, mt: pd.DataFrame, col: str) -> float:
    j = tr.merge(mt, on=["cell_id", "t"])
    j = j[np.isfinite(j[col])]
    return float(np.corrcoef(j[col], j.photons)[0, 1]) if len(j) > 20 and j.photons.std() > 0 else np.nan


def true_background(im, E, ch: str, t: int = 0) -> np.ndarray:
    """The camera image of channel ch in frame t with NO objects (no nuclei, cytoplasm, loci, spots or beads) and no random noise: the true
    background the movie-only map should find (camera offset + offset pattern + hot pixels + background light + autofluorescence)."""
    import copy

    im2 = copy.copy(im)
    if im.nucleus is not None:
        im2.nucleus = copy.copy(im.nucleus)
        im2.nucleus.by_t, im2.nucleus.spots_by_t = {}, {}
    im2.loci = []
    for L in im.loci:
        L2 = copy.copy(L)
        L2.by_t = {}
        im2.loci.append(L2)
    im2.fid = None
    E2 = E.with_enabled([x for x in E.enabled_ids() if x not in ("shot_noise", "read_noise", "row_noise", "cosmic_rays", "dark_current")])
    return im2.frame(t, E2, [ch])["images"][ch].astype(np.float64)


def _compare_selfcal(rep: dict, cal, im, E, cdir: Path) -> None:
    """Add the truth to a movie-only calibration report: hot pixels found / missed / false, the illumination map's error, the true shift,
    and (offset part) the error of the background map against the true object-free background."""
    for ch, r in rep["channels"].items():
        if "offset" in rep["components"]:
            tb = true_background(im, E, ch)
            d = cal.maps[ch]["offset"].astype(np.float64) - tb
            r.update(background_bias=float(np.median(d)), background_p95_abs=float(np.percentile(np.abs(d - np.median(d)), 95)),
                     background_max_bias=float(np.percentile(d, 99)))
        cam = im._camera(ch, E)
        if "hot" in rep["components"] and cam.hot is not None:
            true, found = cam.hot > 0, cal.maps[ch]["hot"]
            r.update(hot_true=int(true.sum()), hot_found_true=int((found & true).sum()), hot_false=int((found & ~true).sum()))
        fl = im._flat(ch, E)
        if "flat" in rep["components"]:
            t = np.ones(im.shape) if fl is None else fl / fl.mean()
            m = cal.maps[ch]["flat"]
            r.update(flat_rms_error=float(np.sqrt(np.mean((m / m.mean() - t) ** 2))), flat_rms_error_even=float(np.sqrt(np.mean((1 - t) ** 2))))
    if "crosstalk" in rep:
        pairs = E.par("crosstalk", "pairs") if E.on("crosstalk") else []
        for ch, r in rep["crosstalk"].items():
            r["true"] = float(sum(f for a, b, f in pairs if a == "nucleus" and b == ch))
    if "shift" in rep and E.on("chromatic_shift"):
        a, b = np.array(im.channel_shift_px("locus0", E)), np.array(im.channel_shift_px("locus1", E))
        rep["shift"]["true_px"] = [float(v) for v in (b - a)]
    (cdir / "selfcal_report.json").write_text(json.dumps(rep, indent=1, default=float))


class _Mini:
    """A throw-away run folder with the same layout as a pipeline run, so the existing analysis functions work on a re-rendered movie."""

    def __init__(self, run: Path, work: Path):
        self.run, self.dir = run, work
        if work.exists():
            shutil.rmtree(work)
        for d in ("stage1_chromatin", "stage2_cells", "stage3_microscopy", "stage6_analysis/truth"):
            (work / d).mkdir(parents=True, exist_ok=True)
        shutil.copy(run / "stage2_cells" / "params.yaml", work / "stage2_cells" / "params.yaml")
        for f in ("transcription.csv", "transcription_events.csv", "loci_truth.csv"):
            if (run / "stage1_chromatin" / f).exists():
                shutil.copy(run / "stage1_chromatin" / f, work / "stage1_chromatin" / f)

    @property
    def s3(self) -> Path:
        return self.dir / "stage3_microscopy"

    @property
    def s5(self) -> Path:
        return self.dir / "stage5_cells" / "truth"

    @property
    def s6(self) -> Path:
        return self.dir / "stage6_analysis" / "truth"

    def images(self, chans) -> dict:
        return {c: self.s3 / f"{c}.tif" for c in chans}


def run_error_effects(run: Path, sources: list[str] | None = None, n_workers: int | None = None, progress=print, keep_work: bool = False) -> pd.DataFrame:
    from simlive.stage3_microscopy.calibration.model import Calibration
    from simlive.stage3_microscopy.render import imager_for_run, save_imaging_truth
    from simlive.stage5_linking.celldata import map_cells
    from simlive.stage5_linking.isolate import isolate_cells
    from simlive.stage6_analysis.bleach import bleach_correct_run
    from simlive.stage6_analysis.localize import localize_loci
    from simlive.stage6_analysis.ms2 import ms2_trace, ms2_truth
    from simlive.stage6_analysis.selfcal import build_selfcal

    run = Path(run)
    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    if not cfg.get("loci"):
        raise ValueError("the run has no loci")
    im = imager_for_run(run, cfg)
    if im.is3d:
        raise ValueError("error_effects works on 2D runs (z-stack runs: not yet)")
    px_um, dt = im.px_um, float(cfg["acquisition"]["frame_interval_s"])
    offset_adu = float(cfg["optics"].get("offset_adu", 0.0))
    cells = pd.read_csv(run / "stage2_cells" / "cells.csv")
    truth0 = pd.read_csv(run / "stage1_chromatin" / "loci_truth.csv")
    ana = {k: v for k, v in (cfg.get("locus_analysis") or {}).items() if k not in ("method",)}
    n_loci = int(cfg["loci"]["n_loci"])
    has_tx = (run / "stage1_chromatin" / "transcription.csv").exists()
    cal_run = Calibration.load(run / "stage3_microscopy" / "calibration") if (run / "stage3_microscopy" / "calibration" / "calibration.json").exists() else None
    s7 = run / "stage7_validation"
    s7.mkdir(parents=True, exist_ok=True)
    out_f = s7 / "error_effects.csv"
    done = pd.read_csv(out_f) if out_f.exists() else pd.DataFrame(columns=["variant", "fix"])
    rows = done.to_dict("records")
    work = s7 / "error_effects_work"
    for name, src, enabled in variants(im.E):
        if sources and src not in sources and name not in sources:
            continue
        fixes = [("default", "default analysis", {})] if name == "reference" else ALL_FIXES if name == "all" else FIXES.get(src, [("none", "not corrected", {})])
        todo = [f for f in fixes if not ((done.variant == name) & (done.fix == f[0])).any()]
        if not todo:
            continue
        t0 = time.time()
        progress(f"{name}: rendering the movie with {len(enabled)} sources on")
        E = im.E.with_enabled(enabled)
        mini = _Mini(run, work / name.replace("+", "plus_"))
        render_variant(im, E, mini.s3, progress)
        save_imaging_truth(im, cells, truth0, mini.s3, E)
        chans = im.channels
        isolate_cells(mini.images(chans), mini.s3 / "labels.tif", mini.s5, margin=(cfg.get("isolation") or {}).get("margin", 0.3), n_workers=n_workers)
        truth = pd.read_csv(mini.s3 / "loci_truth_image.csv")
        mt = ms2_truth(mini.dir) if has_tx else None
        first_pass = None
        for key, label, spec in todo:
            t1 = time.time()
            kw = {**ana, **{k: v for k, v in spec.items() if k not in ("calib", "parts", "post", "movie_parts", "background_stat")}}
            cdir = None
            if spec.get("calib") == "images":
                if cal_run is None:
                    continue
                cdir = mini.dir / f"cal_images_{'_'.join(spec['parts'])}{'_plus_movie' if spec.get('movie_parts') else ''}"
                cal = cal_run.restricted(spec["parts"])
                if spec.get("movie_parts"):             # what calibration images cannot give, taken from the movie
                    mdir = s7 / "error_effects_selfcal" / f"{name.lstrip('+')}__{'_'.join(spec['movie_parts'])}"
                    build_selfcal(mini.images(chans), mini.s3 / "labels.tif", mdir, dt, spec["movie_parts"], None,
                                  float(cfg["acquisition"]["exposure_s"]), px_um)
                    cal.values.update(Calibration.load(mdir).values)
                cal.save(cdir)
            elif spec.get("calib") == "movie":       # kept (small) so the dashboard can show the corrected images
                cdir = s7 / "error_effects_selfcal" / f"{name.lstrip('+')}__{'_'.join(spec['parts'])}{'_' + spec['background_stat'] if spec.get('background_stat') else ''}"
                if first_pass is None and "shift" in spec["parts"]:
                    first_pass = map_cells(localize_loci, mini.s5, n_workers=n_workers, px_um=px_um, n_loci=n_loci, **ana)
                rep = build_selfcal(mini.images(chans), mini.s3 / "labels.tif", cdir, dt, spec["parts"], first_pass,
                                    float(cfg["acquisition"]["exposure_s"]), px_um, spec.get("background_stat", "low"))
                _compare_selfcal(rep, Calibration.load(cdir), im, E, cdir)
            pos = map_cells(localize_loci, mini.s5, n_workers=n_workers, px_um=px_um, n_loci=n_loci,
                            calib=str(cdir) if cdir else None, **kw)
            if kw.get("intensity_background") == "none" and len(pos):       # no shell: take off only the camera's nominal offset
                pos["intensity"] = pos.intensity - offset_adu
            if key == "none" and first_pass is None and not spec.get("calib"):
                first_pass = pos
            tr = None
            if has_tx:
                tr = map_cells(ms2_trace, mini.s5, n_workers=n_workers, calib=str(cdir) if cdir else None,
                               **{k: v for k, v in kw.items() if k not in ("dilate_px", "row_correct")})
                if kw.get("intensity_background") == "none" and len(tr):
                    tr["intensity"] = tr.intensity - offset_adu
            col = "intensity"
            post = spec.get("post", "raw")
            if post in ("global", "percell") and len(pos):
                for f in mini.s6.glob("*.csv"):
                    f.unlink()
                tag = "_cal" if cdir else ""
                pos.to_csv(mini.s6 / f"loci_positions{tag}.csv", index=False)
                if tr is not None:
                    tr.to_csv(mini.s6 / f"ms2_traces{tag}.csv", index=False)
                bleach_correct_run(mini.dir, ("truth",))
                pos = pd.read_csv(mini.s6 / f"loci_positions{tag}.csv")
                if tr is not None:
                    tr = pd.read_csv(mini.s6 / f"ms2_traces{tag}.csv")
                col = "intensity_corr_global" if post == "global" else "intensity_corr"
            sc = score(pos, truth, px_um, col)
            if tr is not None and mt is not None and len(tr):
                sc["ms2_r"] = _ms2_r(tr, mt, col if col in tr else "intensity")
            row = {"variant": name, "source": src, "fix": key, "fix_label": label, **sc, "seconds": round(time.time() - t1, 1)}
            if cdir and spec.get("calib") == "movie" and (cdir / "selfcal_report.json").exists():
                row["selfcal"] = (cdir / "selfcal_report.json").read_text()
            rows.append(row)
            pd.DataFrame(rows).to_csv(out_f, index=False)
            progress(f"  {key:>9}: found {sc.get('found', 0):.3f}, wrong {sc.get('frac_far', np.nan):.3f}, median {sc.get('median_err_nm', np.nan):.0f} nm, "
                     f"distance bias {sc.get('dist_bias_nm', np.nan):+.0f} nm, intensity r {sc.get('intensity_r', np.nan):.3f} "
                     f"err {sc.get('intensity_err_pct', np.nan):.1f}%, late/early {sc.get('late_early', np.nan):.3f}"
                     + (f", MS2 r {sc.get('ms2_r', np.nan):.3f}" if 'ms2_r' in sc else ""))
        if not keep_work:
            shutil.rmtree(mini.dir, ignore_errors=True)
        progress(f"{name}: done in {time.time() - t0:.0f} s")
    return pd.DataFrame(rows)


def summary_json(run: Path) -> dict:
    """The table as {variant: {fix: metrics}} for quick reading."""
    f = Path(run) / "stage7_validation" / "error_effects.csv"
    if not f.exists():
        return {}
    d = pd.read_csv(f)
    return {v: {r.fix: {m: getattr(r, m, None) for m in METRICS} for r in g.itertuples()} for v, g in d.groupby("variant")}


if __name__ == "__main__":       # pragma: no cover
    print(json.dumps(summary_json(Path(".")), indent=1))
