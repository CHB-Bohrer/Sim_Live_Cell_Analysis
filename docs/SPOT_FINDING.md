# Spot finding (stage 6): method, source, measured performance

Stage 6 finds one locus per colour in every frame of every isolated cell. It is still a placeholder for the user's own analysis, but it is
now a defensible one: its steps follow a published pipeline, every step can be inspected in the dashboard (🔎 Spot finding tab), and every
choice below was scored against the simulation's ground truth.

Code: `src/simlive/stage6_analysis/spots.py` (the steps), `localize.py` (2D, per cell), `localize3d.py` (z-stacks), `compare.py` +
`scripts/compare_localizers.py` (method comparison). Config: the `locus_analysis:` block (`configs/loci_*.yaml`).

## Source

Bennabi et al., "Size-dependent temporal decoupling of morphogenesis and transcriptional programs in pseudoembryos", *Science Advances* 11,
eadv7790 (2025), [PMC12372888](https://pmc.ncbi.nlm.nih.gov/articles/PMC12372888). Image-analysis code:
[github.com/gregorlab/SizeProject2024](https://github.com/gregorlab/SizeProject2024), `modules/raw_live_analysis.py`
(functions `get_filtered_fluo_movie_mesp`, `compute_movie_threshold`, `get_fluorescent_area`, `compute_multipol_analysis_per_df`), driver
`exe/multipolarity_analysis/live_movie_multipolarity_analysis.py` (`k=7, th=70`). Read from the code on 2026-10-11, not from memory.

Important: the paper does **not** localize diffraction-limited spots. It finds Mesp2 expression *poles* (regions hundreds of pixels
across) in live gastruloid movies and counts them. Its steps are general (filter, threshold, area, smooth, local maxima), but every size
had to be rescaled, and two steps had to be added for single spots.

## The steps

| step | Bennabi et al. 2025 (Mesp2 poles) | here (one locus per colour, ~0.8 px PSF sigma at 130 nm/px) | config key |
|---|---|---|---|
| 1 median filter | `skimage.filters.rank.median`, 7 x 7 square | `scipy.ndimage.median_filter`, 3 x 3 (z-stacks: (median_z, 3, 3)) | `median_px`, `median_z` |
| 2 region | frame x brightfield gastruloid mask | the cell's own nuclear mask from stage 5, dilated 2 px | `dilate_px` |
| + background | none | a wider (9 x 9) median of the filtered image, subtracted | `background`, `background_px` |
| 3 threshold | ONE per movie: median over frames of per-frame Otsu (or fixed `th=70`) | ONE per movie (cell x colour): median over frames of per-frame min_snr x robust noise; Otsu available | `threshold`, `threshold_scope`, `min_snr` |
| 4 area | > threshold, holes filled, objects < 5000 px removed, eroded (disk 10, disk 2) | > threshold, holes filled, objects < 3 px removed, no erosion | `threshold_on`, `min_area_px`, `erode_px` |
| 5 smoothing | Gaussian sigma 20 px | Gaussian sigma 1 px (~PSF: matched filter) | `smooth_sigma_px`, `smooth_sigma_z` |
| 6 local maxima | `peak_local_max(min_distance=100)`; every maximum is a pole | `peak_local_max(min_distance=4)` inside the area = candidates | `min_sep_px` |
| + choice | (all poles kept and counted) | `track`: the candidate nearest the locus' median position in the 3 frames on each side (none if > 8 px away); or `brightest` | `select`, `track_window`, `track_gate_px` |
| 7 position | integer pixel of the maximum | 2D (3D) Gaussian + offset fit to the unfiltered pixels in a 7 x 7 (x 5 planes) window | `refine` |

Why the sizes change: the paper's 7 px median is tiny next to a pole, so it only removes pixel noise. A 3 x 3 median already halves the peak
of a sigma 0.8 px spot and a 5 x 5 median loses 10-25% of the spots (measured below); a 7 x 7 median erases spots, which is why a wider
median is a good estimate of the background *under* a spot. The 3 x 3 median removes hot pixels and single-pixel noise. A cosmic ray is as
wide as the PSF in this camera model and survives any spatial median.

Two additions (not in the paper), both measured below:
- **local background**: the locus channels carry the nucleus' cross-talk ghost (with its texture), autofluorescence and haze. One background
  value per cell inflates the noise estimate and raises the threshold: locus 0 (the channel with the ghost) lost most detections.
- **track choice**: a cosmic ray lasts one frame; the locus stays near where it was. In max projections of z-stacks every plane's cosmic rays
  end up in the 2D image (17 planes: ~17x the rate), and the brightest-candidate rule picked them in ~4.5% of frames.

The default thresholds the *smoothed* image (`threshold_on: smoothed`) rather than the median-filtered one as the paper does: the Gaussian is
then also the detection filter, which has the better signal-to-noise ratio for a spot. `threshold_on: filtered` is the paper's order.

## Measured (stage 6 re-run on saved crops, true cell IDs; `scripts/compare_localizers.py`; 2026-10-11)

"found" = detections / loci that had light; "wrong" = detections > 0.5 um from the truth (a hot pixel, cosmic ray or noise peak was picked);
errors in nm against the drift-corrected truth (the chromatic shift counts as error unless calibrated). One seed per run: differences of a
few percent in "found" are within run-to-run noise.

Library run (seed 11, 8 cells x 120 frames, probes + MS2 on, all image errors; the hardest: spots 60-570 photons):

| finder | calibrated | found | wrong | median | p95 | rms of correct ones |
|---|---|---|---|---|---|---|
| original (brightest pixel + 5x5 centroid) | no | 86.5% | 9.8% | 54 | 4681 | 64 |
| paper literal (cell background, Otsu, brightest, integer) | no | 100% | 15.1% | 79 | 6519 | 86 |
| **default** | no | 89.8% | 0.06% | 48 | 96 | 59 |
| original | yes | 82.2% | 0.13% | 27 | 76 | 42 |
| **default** | yes | 90.7% | 0.06% | 25 | 74 | 41 |

Toy-loci run (seed 12, probes + MS2, all errors): original 96.3% found / 3.8% wrong / p95 104 nm (uncalibrated) and 95.9% / 0.1% / median 25 /
p95 61 nm (calibrated); default 97.0% / 0% / p95 83 nm and 97.9% / 0% / median 21 / p95 54 nm.
Older toy run (seed 1, bright spots): every variant finds 100%; Gaussian fit p95 64 nm vs 72 nm for the centroid; median 44 nm in all (the
uncorrected chromatic shift, ~0.3 px).

3D run (seed 22, 4 cells x 40 frames, 17 planes 0.8 um apart, all errors, fiducial beads):

| finder | analysis | calibrated | found | wrong | median | p95 |
|---|---|---|---|---|---|---|
| original | 2D on max projection | yes | 100% | 4.7% | 40 | 110 |
| **default** | 2D on max projection | yes | 99.4% | 0% | 39 | 80 |
| original | 3D | yes | 99.7% | 4.7% | 101 | 320 |
| **default** (axial median off: auto, spacing > axial sigma) | 3D | yes | 99.1% | 0.3% | 84 | 183 |
| default + axial median always (`median_z: 3`) | 3D | yes | 99.1% | 0% | 85 | 177 |

What each one-at-a-time change does (library run, calibrated, vs default 90.7% found / 0.06% wrong / 25 nm):
no median 89.6% / 0% (but uncalibrated: 4.1% wrong, the hot pixels); 5 x 5 median 86.4% found; **one background per cell 56.0% found**; Otsu
100% found but 1.4% wrong (rejects nothing); threshold on the filtered image 82.9% found; brightest choice 0.7% wrong; per-frame threshold
same as one per movie; integer position median 59 nm (the paper's choice costs ~35 nm at this scale); centroid 29 nm and p95 84 nm.

## Where it still fails / open

- Dim spots (library run, calibrated, photons per exposure after bleaching): 100-200 photons found 71% (original 52%), 200-300 photons 91%
  (original 81%), 300-400 photons 97% (original 95%). By colour: locus 0, the channel with the nucleus' cross-talk ghost, 88% (original 69%);
  locus 1 94% (original 95%). Lower `min_snr` finds more but `track` then has more noise peaks to choose from; re-check with real brightness.
- `median_z: auto` turns the axial median on only when the plane spacing is at most the axial PSF sigma. With bright spots it helped even at
  0.8 um spacing (0% vs 0.3% wrong); with dim spots it would remove the spot. Decide with real brightness and spacing.
- Background structures as small as a spot (< ~0.5 um) cannot be told from a spot by any of these steps (`tests/test_spots.py` uses a 1.3 um
  patch for that reason).
- `track` assumes one locus per colour that moves less than `track_gate_px` (1 um) relative to the nucleus centroid within `track_window`
  frames. Two copies of a locus (replicated sister chromatids), or a fast-moving locus, would need a real tracker.
- 3D z error dominates the 3D error (lateral ~40 nm, 3D ~85 nm): coarse planes, focus error and axial drift (both counted as error; the
  fiducial / piezo corrections in docs/CALIBRATION_AND_3D.md address them).
- All numbers come from 4 simulated runs, one seed each, with PLACEHOLDER brightness and error sizes. A real annotated movie is the test.

## Logic errors found and fixed in the stage 5-6 path (2026-10-11 review)

1. Cell-frame coordinates (`u_um`, `v_um`) were measured from the crop centre, which is the *rounded* centroid: up to 0.5 px (65 nm) off.
   Now from the true centroid (`frames.csv` cy_px/cx_px minus the crop origin). Test: `test_uv_coordinates_are_relative_to_the_true_nucleus_centroid...`.
2. `locus_analysis.min_sep_px` was in every config but silently ignored (swallowed by `**_`). It is now the min distance of step 6.
3. Scoring counted loci with no light at all (probe model: no probe attached for a whole exposure) as misses, and its rms was dominated by a
   few wrong picks without saying so. Now `detection_rate` is over loci with light (`photons` > 0, saved in `loci_truth_image.csv`),
   `n_false_on_dark` counts detections where nothing was drawn, `frac_far` the wrong picks, and `rms_error_near_um` the accuracy of the rest.
4. The single-cell player drew the true positions of *every* cell's loci, so a neighbour's locus inside the crop window appeared as a false
   "true" circle. Now only the cell's own.
5. The player's + and o markers were half a pixel up-left of the spot (PIL coordinates are pixel corners; positions are pixel centres).
6. The old 5 x 5 centroid around the brightest pixel is biased towards the pixel centre ("pixel locking"): replaced by a Gaussian fit
   (centroid vs fit, calibrated library run: median 29 vs 25 nm, p95 84 vs 74 nm).
