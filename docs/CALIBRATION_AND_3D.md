# Calibration of the imaging instrument, and 3D z-stacks

Both are part of stage 3 (the microscope). All numbers are **placeholders** until measured on the real microscope.

## 1. 3D z-stacks (`acquisition.z_stack.enabled: true`, example `configs/loci_3d.yaml`)

**What changes.** Every frame is a stack of `n_slices` planes `step_um` apart (default 15 planes x 0.8 um). Nuclei are flattened ellipsoids
(`geometry.z.nucleus_semi_axis_ratio`, 1 = sphere; half thickness at a pixel = z_semi x sqrt(1 - q^2), q = relative radius), with a height
(`z_center_um`, scatter `z_spread_um`, random walk `z_walk_um_per_frame`; daughters start at the parent's height). Loci have a height
`z_um` (toy loci: confined diffusion in 3D; library loci: the polymer's z, scaled by the nucleus' z half axis, so keep the ratio at 1 if the
transcription model's polymer distances must equal the imaged ones). The nucleus texture is the same at every height.

**Point spread function.** 3D Gaussian core (lateral sigma 0.21 lambda / NA, x0.8 for confocal; axial sigma from the modality: widefield
FWHM_z = 2 lambda n / NA^2, confocal 0.88 lambda / (n - sqrt(n^2 - NA^2)), or `optics.axial_fwhm_um`) plus the broad out-of-focus haze (source `haze`).
A point source's flux in the in-focus plane equals its photon count, falling as a Gaussian of the axial distance; a nucleus' intensity is its
local average over the axial PSF. This is an approximation of a real widefield PSF (which has long axial tails: use `haze` for them).

**Errors specific to z.** `stage_drift` (axial drift + jitter), `defocus` (in a stack: every object appears displaced along z by the focus
error), `chromatic_shift` (`shift_z_um` per colour), `slice_skew` (planes are recorded one after the other, `slice_time_s` apart; loci
move between planes, positions interpolated between frames). Every plane has independent camera noise; fixed patterns are the same in every plane.

**Files** (`stage3_microscopy/`): `<channel>_zstack.tif` (T, Z, Y, X), `labels_3d.tif`, plus the 2D projections `<channel>.tif` and `labels.tif`
that segmentation, tracking and cell isolation use (maximum projection; nuclei never overlap laterally so the projected labels keep the IDs).
`imaging_errors.json` has the z-stack geometry (`z_positions_um`, axial sigma, channel z shifts); `imaging_truth.csv` has `drift_z_um`;
`loci_truth_image.csv` has `z_um` (drift-corrected truth) and `z_drawn_um` (where the spot is drawn: + chromatic z shift + focus error).

**Analysis.** `stage6_analysis/localize3d.py` (placeholder): re-reads each cell's crop from the full stacks, smooths in 3D, takes the brightest
voxel inside the cell mask and refines it with a centroid. Output `loci_positions_3d.csv`, `localization_errors_3d.csv` (3D, lateral and axial
error), scores `truth_3d` in `stage6_analysis/scores.json`.

**Cost.** frames x planes x pixels^2 x channels x 2 bytes: 40 frames x 15 planes x 385^2 x 3 channels = 0.5 GB. Rendering is about 4x the 2D cost per plane.

## 2. Calibration (`stage3_microscopy/calibration/`, config block `calibration:`, tab 🎯 Calibration, `scripts/run_calibration.py`)

**Idea.** A real calibration takes pictures of things whose true appearance is known and estimates the instrument from them. Here the
simulation does the acquisitions with the SAME instrument as the movie (same camera patterns, illumination, chromatic shifts; independent
noise, drift, focus and flicker), the estimators see only the images, and afterwards every estimate is compared with the true value. That
checks the calibration method before it is used on real data, and gives the corrections to apply to the movie.

| acquisition | what is simulated | what it gives |
|---|---|---|
| dark frames | shutter closed, exposures 0.01 / 1 / 5 s (+ the run's), 60 frames at the shortest and the run's | offset map, read noise, row noise, offset pattern (DSNU), dark current (slope vs exposure), hot pixels, cosmic-ray rate and charge |
| flat-field slide | uniform slide, 7 brightness levels x 12 frames, + a 100-frame series | photon transfer curve -> gain (ADU per detected electron); illumination pattern (constant + tilt + quadratic fit: vignetting, tilt); pixel gain spread (PRNU); flicker |
| bead slide | 60 bright sub-resolution beads in every colour, time lapse (10 frames; 6 stacks in 3D) | PSF width per colour (lateral and axial), chromatic shift (lateral and axial), stage drift and focus drift paths |

Sample-related sources (haze, autofluorescence, cross-talk, background, bleaching, locus motion, slice skew) are off for calibration targets.
Not measurable from these images: the quantum efficiency (needs a source of known photon flux).

**Results on the placeholder instrument** (loci_demo, 40 um field, 3 channels, ~15 s): gain within 0.5%, offset within 0.02 ADU, read noise within 1%, row noise
and DSNU within 1%, dark current within 2%, 100% of the hot pixels found with 0 false alarms, PRNU within 0.2%, vignetting within 0.001 (0.1496 vs 0.15),
flicker within 10%, chromatic shifts within 0.003 px (0.0003 um axially), PSF width within 0.3%, drift path within 0.002 px. The cosmic-ray rate is
only roughly known (a Poisson count of rare events). Check `comparison.csv` for any run.

**Applying it** (`Calibration.correct`, `correct_positions`): counts -> (counts - offset map) / normalised flat, hot pixels replaced by the 3x3 median, detected
positions minus the measured chromatic shift (lateral and axial). The pipeline (`calibration.simulate`/`apply`, run after tracking) analyses the loci a second time this way
(`loci_positions_cal.csv`, scores `truth_calibrated`, `truth_3d_calibrated`). Seed 11, library loci, all errors on: median locus error 54 -> 27 nm, 95th percentile 4.7 um -> 76 nm,
rms 1.57 um -> 0.16 um (hot pixels were fooling the placeholder detector). Cosmic rays are random per frame and are NOT corrected: in z-stacks each plane has its own,
so ~6% of 3D detections are still wrong (3D run: median error 193 -> 90 nm, but the 95th percentile stays ~4 um). A real analysis needs outlier rejection over time.

**Gotcha.** The fits use numpy LAPACK (`lstsq`, `curve_fit`), which loads a second OpenMP runtime; importing torch AFTER that fails on Windows. The pipeline therefore runs the
calibration after the tracker (torch already imported). Do not call it before `import torch` in the same process.

## 3. Piezo plane-position errors and fiducial beads (tab 📍 Stage & fiducials)

**The problem.** A z-stack is made by stepping a piezo through planes at nominal heights; locus heights are read from those nominal heights. The source `piezo`
(`imaging_errors.piezo`) makes the true plane heights differ from the nominal ones: error of plane k = random N(0, `repeatability_nm`) + `bow_um` x (1 - u^2)
(u = -1..+1 across the stack) + direction x `hysteresis_um` / 2 - direction x `settling_um` x exp(-recording order / 2). The scan direction alternates every frame when
`bidirectional`. Images are drawn with the planes at their TRUE heights; analysis is told the nominal heights. Truth: `stage3_microscopy/piezo_truth.csv` (per frame and plane:
nominal height, true height, error, the sensor readout = true + `readout_noise_nm`, scan direction, recording order). The values are placeholders: ask the user for their piezo's
repeatability / hysteresis.

**Fiducial beads** (`fiducials:` block, `stage3_microscopy/fiducials.py`): bright sub-resolution beads ON THE SLIDE (z = `z_um` above the coverslip) in every image and every colour, placed away from
nuclei, drifting with everything else; truth `fiducials.csv`. A stack must reach the slide (configs/loci_3d.yaml: 17 planes x 0.8 um from z = 0.2) or the beads are almost invisible
(only out-of-focus haze). The sidebar lowers the stack automatically when the beads are ticked.

**What each tool can and cannot correct** (`stage6_analysis/fiducials.py`):
- bead tracking (`track_fiducials`: beads found outside the nuclei in frame 0, a Gaussian fit per bead per frame, median displacement relative to frame 0) measures the movement of the
  WHOLE picture (drift, focus error) -> `correct_with_beads` subtracts it. It cannot see plane-to-plane errors at nuclear height, because all beads are at one height.
- the piezo's position log (`readout_um`) gives every plane's real height -> `correct_with_plane_positions` replaces the nominal heights (interpolating the fitted fractional plane index).
  With the log, the beads are fitted with the logged heights, so they measure drift + focus only.
- `score_variants` (uncorrected / beads / piezo log / both) reports MEDIAN absolute errors of: the 3D locus0-locus1 distance, their height difference, and each height after the
  movie-wide constant is removed. Medians, because the placeholder locus finder has rare wrong detections (cosmic rays) that would otherwise hide everything.

**Result with the placeholder values** (loci_3d, seed 22, 4 cells, 8 beads, repeatability 30 nm, bow 0.1 um, hysteresis 0.05 um): height wander 84 -> 65 nm with the beads (the whole-picture shift is removed);
the piezo log changes the height-difference error only 83 -> 80 nm and the distance error 50 -> 46 nm, because at 30 nm the piezo errors are small next to the placeholder detector's own
~80 nm height noise. Increase `repeatability_nm` to see the plane-to-plane errors dominate: only the log helps then. Do not over-read one seed.
