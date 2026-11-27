# Image-error model, z-stacks, calibration, piezo and fiducials, biplane

Working notes moved verbatim from `CLAUDE.md` on 2026-11-27 (section order kept). Add new dated sections at the top; delete statements that are no longer true.

## Image errors (stage 3, built 2026-10-09; the user asked for ALL sources of image error)
- `stage3_microscopy/errors.py` = table `SOURCES` (22 sources: stage drift, chromatic shift, PSF, defocus, haze, exposure motion, bleaching, flat-field,
  flicker, constant background, autofluorescence, cross-talk, QE, shot noise, dark current, hot pixels, cosmic rays, PRNU, read noise, row noise, DSNU,
  ADC/saturation), `ORDER` (physical order), `ErrorSet` (on/off + parameters from the config block `imaging_errors:`), per-source RNG streams
  (`rng_for(entropy, t, source, channel)`), and `Camera`. `render.py` = `Imager`: renders any single frame with any subset of sources; `render_nuclei` /
  `render_loci` are one-channel wrappers. `docs/IMAGE_ERRORS.md` is generated from the table (`scripts/make_image_error_docs.py`; regenerate after editing it).
- `imaging_errors.enabled: false` (or no block) = the ORIGINAL model only (baseline sources: PSF, bleaching, constant background, QE, shot noise, read noise).
  The three configs (`tracking_demo`, `loci_demo`, `loci_library`) now have the block ON with PLACEHOLDER values. Existing tests pass
  `imaging_errors.enabled=false`. QE is now applied before Poisson (electrons ~ Poisson(QE x photons)), slightly more noise than the old Poisson-then-scale.
- Determinism: frame t of (config, truth, entropy) is identical whatever else is rendered, and switching a source never changes another source's noise.
  `entropy_for(seed)` = `SeedSequence(seed).spawn(6)[1]`; the dashboard's Images tab rebuilds the Imager with `imager_for_run` and checks the rebuilt frame equals
  the saved TIFF (it does, bit for bit).
- Ground truth per run in `stage3_microscopy/`: `labels.tif` (moved by the drift), `imaging_truth.csv`, `cells_image.csv`, `loci_truth_image.csv` (drift-corrected
  truth in y_um/x_um + `y_drawn_um/x_drawn_um` = where the spot is drawn, incl. chromatic shift), `imaging_errors.json`. Stage 6 localization is scored against
  `loci_truth_image.csv` (drift cannot be known by an analysis; the chromatic shift IS counted as error). Dashboard Cells tab/trajectories use the image-space files.
- Findings (1 seed each, placeholder error sizes; do not over-read): default hard tracking config, 3 seeds: CHOTA 0.909 with errors vs 0.903 without (id switches
  14/23/23 vs 14/21/25) -> bright nuclear images are robust to these errors. Loci (library smoke, MS2 + probes on, seed 11): the OLD placeholder brightest-spot
  localizer's (replaced 2026-10-11, see "Spot finding") detection rate falls 99.9% -> 85.1% and its RMS error 0.019 -> 2.05 um (median 16 -> 52 nm, p95 5.8 um) with all errors on: hot pixels, cosmic rays and
  the cross-talk ghost of the nucleus fool a "brightest spot in the mask" detector. That is the kind of failure this framework is meant to expose; the user's own
  analysis will replace the placeholder.
- Dashboard: new 🖼 Images tab (`app/images_tab.py`): (1) recorded frame per channel with true outlines / IDs / locus positions (cyan o) / drawn positions (red x) /
  MS2 promoters, optional error-free and difference panels, and the imaging truth numbers of that frame; (2) step-by-step build-up (ground truth -> + each source in
  physical order, same grey scale, with a table of what each adds and a brightness profile); (3) error budget (each source alone, RMS change and spot-centre shift);
  (4) reference of every source with this run's settings and derived numbers (PSF width, Thompson localization limit, drift, chromatic misregistration). Sidebar forms
  have "Image errors" controls (master switch + sources to switch off). The run summary, glossary, app map, intro boxes and tooltips were added the same day.
- Gotcha: sources whose size you change at use time reuse the same unit random series (drift, focus, flicker) so results stay comparable; fixed patterns (PRNU,
  DSNU, hot pixels, flat-field, autofluorescence) are cached on the ErrorSet instance (`E.cache`), so create variants with `ErrorSet.with_enabled`, not by mutating.
- Not done: motion blur of the nuclei (negligible at these exposures), EMCCD gain, rolling shutter, 3D / z-stacks, per-channel focal shift, sample-dependent
  (thick-sample) scattering, measured values for any of the 22 sources (all placeholders; ask the user for bead / dark-frame / flat-field calibration images).

## 3D z-stacks and calibration (built 2026-10-11; full description in docs/CALIBRATION_AND_3D.md)
- z-stacks: `acquisition.z_stack.enabled` (+ `geometry.z`), example `configs/loci_3d.yaml` (toy loci; override the same blocks on the library config). Cells get `z_um`, `z_semi_um`
  (`stage2_cells.motion.add_z`, own random stream `ss[6]`: the lateral simulation is unchanged), loci get `z_um` (toy and library). `Imager` renders (nz, y, x) frames: 3D Gaussian PSF core + haze,
  nuclei = flattened ellipsoids, per-plane noise (`t * 4096 + k` into the camera RNG). New sources: `slice_skew`; z parameters in `stage_drift` (drift_z), `chromatic_shift` (shift_z_um), `defocus` (axial
  shift in a stack). Saved: `<ch>_zstack.tif` (T,Z,Y,X), `labels_3d.tif`, plus 2D projections `<ch>.tif` / `labels.tif` (max) that stages 4-5 use unchanged. `localize3d.py` (placeholder 3D localizer) +
  `score_localization_3d`; scores `truth_3d`. Pipeline stream layout is now `SeedSequence(seed).spawn(10)`: 0 motion, 1 imaging, 2 loci, 4 probes, 5 transcription, 6 heights, 8 copies, 7 and 9 free (entropy_for still uses stream 1).
- calibration (`stage3_microscopy/calibration/`: `sim.py` acquisitions of the SAME instrument with their own session, `fit.py` blind estimators (dark / photon transfer curve / flat-field / flicker / beads),
  `model.py` `Calibration` (maps, values, `correct`, `correct_positions`, save/load), `run.py` `run_calibration` -> (Calibration, estimate-vs-truth table); `scripts/run_calibration.py`; config block `calibration:`).
  Imager now separates `entropy` (instrument: fixed patterns) from `session` (drift, focus, flicker, per-frame noise; equal for a movie). The pipeline runs it AFTER the tracker (numpy LAPACK loads a second
  OpenMP runtime; importing torch after that fails, see the gotcha above) and analyses the loci again with the calibration applied (`loci_positions_cal.csv`, `truth_calibrated`, `truth_3d_calibrated`).
- Results (placeholder instrument; one seed each): 53 of 53 compared quantities within tolerance (gain 0.3%, read noise 0.1%, hot pixels 40/40 and 0 false, vignetting 0.1496 vs 0.15, chromatic shift
  0.001 px / 0.0003 um, PSF width 0.3%, drift 0.001 px); the cosmic-ray rate is only roughly known. Applying it: 2D library run (seed 11) locus median error 54 -> 27 nm, p95 4.7 um -> 76 nm, rms 1.57 -> 0.16 um;
  3D run (seed 5) median 3D error 193 -> 90 nm but p95 stays ~4 um because cosmic rays (independent in every plane, ~0.75 events per stack per channel) fool the brightest-voxel detector (OLD finder;
  the 2026-10-11 spot finder's track choice addresses this: see "Spot finding").
- Dashboard: 🖼 Images has a z-slice slider / projection toggle, truth markers with their z offset, and x-z / y-z side views; new 🎯 Calibration tab (estimate vs truth table, dark / flat / bead panels, photon transfer
  curve, before-after scores); sidebar loci form has field of view, 3D z-stacks (planes, spacing, confocal) and calibration switches (warns about disk size, refuses > 12 GB).
- Piezo + fiducials (2026-10-11, docs/CALIBRATION_AND_3D.md section 3): source `piezo` (true plane heights = nominal + random + bow + hysteresis + settling lag; alternating scan direction; truth
  `piezo_truth.csv` incl. the sensor `readout_um`), `fiducials:` block (beads ON THE SLIDE in every image; the stack must reach z = 0.3 um, loci_3d.yaml is now 17 planes from z = 0.2),
  `stage6_analysis/fiducials.py` (bead tracking -> whole-picture shift; piezo-log correction; `score_variants` with MEDIAN errors), tab 📍 Stage & fiducials. Key point shown in the UI: beads on the slide
  cannot measure plane errors at nuclear height; only the piezo log can. RNG gotcha: `errors.INDEX` is FROZEN (new sources get the next number) so adding a source to ORDER does not change old noise.
  With the placeholder 30 nm repeatability the piezo effect is small next to the placeholder 3D locus finder's ~80 nm noise (beads: height wander 84 -> 65 nm; log: 83 -> 80 nm).
- Not done: real (measured) values for any placeholder; 3D transcription spots in thick tissue; z-dependent chromatic focal shift beyond a constant offset; cell isolation in 3D (stage 5 stays 2D on projections);
  3D tracking; calibration of the quantum efficiency and of the nuclear brightness; lens distortion / field-dependent chromatic shift (shifts are constant over the field).

## Vutara biplane mode (2026-10-17; user: "change it, and why not build in the biplane ... you still need to do a zstack ... calibrate with the beads within the image")
- The sidebar label "Bruker Vutara VXL (2D)" was misleading (the instrument is 3D-capable; the preset modelled only its widefield mode): renamed "(widefield mode)"; new third choice "(biplane 3D mode)" = `configs/instruments/vutara_vxl_biplane.yaml`.
- Facts (Bruker bi-plane page, fetched 2026-10-17): two planes 600 nm apart, 50:50 splitter, two halves of one camera, z from the intensity difference by a 3D-PSF fit, PSF calibrated with a bead imaged as a function of z, z-stacks possible; NO axial range per acquisition and
  no registration step are stated. So a z-stack is still needed for a ~10 um nucleus (the user was right).
- Code: `acquisition.z_stack.biplane: {enabled, separation_um: 0.6, photon_fraction: 0.5}` in `stage3_microscopy/render.py` (`Imager.biplane`; `n_slices` = stage positions, planes = 2x, stage step = 2 x separation so the planes are EVEN 0.6 um and the 3D localizer is unchanged; both planes of a pair share the
  time index and the piezo error; each plane gets `plane_fraction` of the signal, background not halved); new error source `biplane_registration` (second plane shifted `shift_px`, its height `separation_error_um` off; in `_LATER`, off unless a config enables it; listed in `error_effects.ONLY_3D`).
  `stage6_analysis/biplane.py` `estimate_plane_offset` (beads between the lowest pair -> offset) + `localize_loci_3d(plane_shift_px=)`; pipeline `_undo_biplane_registration` writes `truth_3d_registered` and `biplane_registration` in scores.json. `tests/test_biplane.py` (5 tests). 176+ CPU tests.
- Measured once (seed 61, 2 cells, 12 frames): offset (0.321, -0.214) px vs true (0.35, -0.25); locus error unchanged within noise (dominated by the placeholder axial chromatic shift; calibration was off); biplane vs an ordinary 22-plane z-stack: same 3D rms (0.27 um), lateral rms 0.06 vs 0.10 um. ONE seed: not evidence.
- NOT done: true plane-separation and PSF z-calibration from a bead z-scan; registration with rotation / scale; a bead SLIDE calibration of the offset; biplane in the E-P tab's vector model is separate (`ep_acquisition.py`); the Images/Stage tabs have not been looked at with a biplane run.
