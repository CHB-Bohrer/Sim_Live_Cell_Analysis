# Spot finding (stage 6): methods, sources, measured performance

Stage 6 finds one locus per colour in every frame of every isolated cell. It is still a placeholder for the user's own analysis, but it is
now a defensible one: two published pipelines are implemented, every step can be inspected in the dashboard (🔎 Spot finding tab), and every
choice below was scored against the simulation's ground truth.

| `locus_analysis.method` | source | code |
|---|---|---|
| `dog_viterbi` (**default**) | Chen, Levo, Zoller & Gregor 2025, *Nat. Struct. Mol. Biol.* (MS2 transcription spots) | `stage6_analysis/dog_viterbi.py` |
| `median_peaks` | Bennabi et al. 2025, *Sci. Adv.* (Mesp2 poles), rescaled to spots | `stage6_analysis/spots.py` |
| `brightest` | the original placeholder | `localize.py` |

Code around them: `localize.py` (2D, per cell), `localize3d.py` (z-stacks), `compare.py` + `scripts/compare_localizers.py` (method comparison).
Every method also reports Chen et al.'s ellipsoid **intensity** of each detection, scored against the photons the simulation drew.

# Part A: `dog_viterbi`, after Chen et al. 2025 (default)

## Source

Chen, Levo, Zoller & Gregor, "A conserved coupling of transcriptional ON and OFF periods underlies bursting dynamics", *Nat. Struct. Mol.
Biol.* (2025), doi [10.1038/s41594-025-01615-4](https://doi.org/10.1038/s41594-025-01615-4); methods read from the open arXiv version
[2304.08770v3](https://arxiv.org/abs/2304.08770) (section III.2 "Transcription spot tracking and signal integration") on 2026-10-11. No code
was published, so this implementation follows the text. Their setting: MS2 spots of ONE labelled allele per nucleus on the background of free
MCP-GFP in the nucleus, two-photon microscope, 40x/1.3 NA, 220 nm pixels, 12 planes 1 um apart, 10 s per stack. That is the situation of a
locus (or the MS2 spot) in one of our cells; in units of the PSF their pixel sizes match ours, so their numbers are kept.

## The steps

| step | Chen et al. 2025 (text) | here | config key |
|---|---|---|---|
| 1 one stack per nucleus | nuclei segmented (DoG + watershed), tracked by Voronoi overlap; one 3D stack per nucleus over time | stage 5's single-cell crops | `dilate_px` |
| 2 median filter | 3D median filter against salt-and-pepper noise | 3 x 3 (x `median_z` planes in z-stacks; auto) | `median_px`, `median_z` |
| 3 candidates | 3D difference of Gaussians, "tuned to high sensitivity": several candidates per nucleus | DoG sigma 1 and 3 px (z: 0.8 and 2 planes); local maxima above 3 x the DoG noise; at most 5 | `dog_sigma_px`, `dog_snr`, `max_candidates` |
| 4 trajectory | "modified Viterbi": at most one spot per nucleus; links scored by pairwise distance with a diffusion potential built iteratively from the random walk of the spots | Viterbi over the whole movie (below) | `viterbi_*` |
| 5 position | centroid of the spot | Gaussian fit (`refine: gaussian`; `centroid` = as described) | `refine` |
| 6 intensity | mean in a (7, 7, 5) px ellipsoid minus the MEAN of the (19, 19, 7) px shell | same ellipsoids; the shell's MEDIAN by default (`mean` = as described) | `intensity_*` |

The Viterbi score (the text gives no formula; this is the standard track likelihood): nodes are candidates, positions relative to the nucleus
centroid (z scaled to pixels). Link cost over a gap of g frames = -log of the diffusion propagator, |d|^2 / (4 D g) + (ndim/2) log(4 pi D g);
each skipped frame costs `viterbi_miss_cost` (2); gaps up to `viterbi_max_gap` (5). Node reward = log(SNR) + log(region size / false
candidates per frame). The second term is the clutter density of a track likelihood: a candidate left off the path must be explained as a
false one. **Without it the path skipped 547 of 1914 frames that held the right candidate** (the extra link normalisation made skipping
cheaper than keeping a faint spot); with it, 1 of 1914 (`tests/test_spots.py::test_viterbi_keeps_a_faint_spot...`). D starts from the
brightest-candidate path and is re-estimated from the chosen path's own steps (median of |d|^2 / g, robust) up to 3 times, per cell and
colour (the paper pools all spots of a movie; `viterbi_D_px2` fixes it). Estimated D here: 1.4-2.3 px^2 per frame (0.13 um, 10 s frames).

## Measured: dog_viterbi (default) against the alternatives

Same protocol as part B below (stage 6 re-run on the saved crops with `scripts/compare_localizers.py`; "found" over loci with light; "wrong" =
> 0.5 um from the truth; errors in nm; intensity r = Pearson correlation of the measured intensity with the photons drawn in that frame,
correct detections only, locus 0 / locus 1). One seed per run.

| run | finder | cal. | found | wrong | median | p95 | intensity r |
|---|---|---|---|---|---|---|---|
| library (probes + MS2, all errors) | original | yes | 82.2% | 0.13% | 27 | 76 | 0.79 / 0.82 |
| | median_peaks | yes | 90.7% | 0.06% | 25 | 74 | 0.80 / 0.81 |
| | **dog_viterbi** | yes | **98.0%** | **0%** | 26 | 98 | 0.78 / 0.81 |
| | original | no | 86.5% | 9.8% | 54 | 4681 | 0.77 / 0.80 |
| | **dog_viterbi** | no | **97.8%** | 0.1% | 50 | 124 | 0.76 / 0.79 |
| toy (seed 12, probes + MS2, all errors) | original | yes | 95.9% | 0.1% | 25 | 60 | 0.83 / 0.87 |
| | median_peaks | yes | 97.9% | 0% | 21 | 54 | 0.83 / 0.86 |
| | **dog_viterbi** | yes | **99.5%** | 0.05% | 22 | 55 | 0.83 / 0.87 |
| 3D (seed 22), 2D on projections | original | yes | 100% | 4.7% | 40 | 110 | |
| | **dog_viterbi** | yes | 100% | **0%** | 39 | 81 | |
| 3D (seed 22), in 3D | original | yes | 99.7% | 4.7% | 101 | 320 | 0.40 / 0.59 |
| | median_peaks | yes | 99.1% | 0.3% | 84 | 183 | 0.41 / 0.56 |
| | **dog_viterbi** | yes | **100%** | **0%** | 84 | 177 | 0.40 / 0.58 |

A second 3D run (seed 23, made with dog_viterbi as the pipeline default; harder: the original picks the wrong object 7-10% of the time):
calibrated 2D found / wrong / p95 = original 100% / 7.2% / 3258 nm, **dog_viterbi 100% / 0.3% / 91 nm**, median_peaks 98.8% / 0% / 90 nm;
calibrated 3D = original 100% / 6.9% / 5137 nm, **dog_viterbi 100% / 0% / 274 nm** (median 106 nm), median_peaks 99.1% / 0.3% / 275 nm;
uncalibrated 3D wrong picks: original 9.7%, dog_viterbi 1.6%, median_peaks 0.9%. With bright spots the two new methods are about equal.

Intensities use the median of the background ring (the default, for every method). With the paper's MEAN, one frame of the toy run where a
bright object sat in the ring got a negative intensity, and locus 1's r fell from 0.87 to 0.75 (dog_viterbi) or 0.72 (median_peaks).

Reading it: dog_viterbi finds 7-16% more of the spots on the hard (faint, probe-model) run without picking more wrong objects; the median
error is unchanged. With bright spots (the 3D runs) it and median_peaks are about equal (both far better than the original). Its
95th percentile is higher there only because the extra detections are the faint spots (100-300 photons), which no method locates well; filter
on `snr` if a use needs only bright ones. Variants as described in the paper (`chen_literal`: centroid, mean background) and one-change
variants are in each run's `stage6_analysis/method_comparison.csv` (🔎 Spot finding tab, part 4). The intensity follows the true photons only
moderately in 2D (r ~0.8 in the probe runs: shot noise of ~300-photon spots, flat-field). Read r only where the true brightness varies: in the
3D runs the probe model is off, every spot has the same 600 photons apart from slow bleaching, so r is low (0.2-0.6) whatever the method; the
measured intensity there mostly varies with the spot's height between the 0.8 um planes.

Open for dog_viterbi: one path per cell, so an OFF period longer than `viterbi_max_gap` frames ends the track (fine for loci, which are
always there; for MS2 spots, add re-acquisition after long gaps); D per cell rather than pooled; DoG widths and `dog_snr` set by reasoning,
not tuned; the paper's two-photon PSF and our widefield PSF differ in axial shape.

# Part B: `median_peaks`, after Bennabi et al. 2025

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

median_peaks thresholds the *smoothed* image (`threshold_on: smoothed`) rather than the median-filtered one as the paper does: the Gaussian is
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
| **median_peaks** | no | 89.8% | 0.06% | 48 | 96 | 59 |
| original | yes | 82.2% | 0.13% | 27 | 76 | 42 |
| **median_peaks** | yes | 90.7% | 0.06% | 25 | 74 | 41 |

Toy-loci run (seed 12, probes + MS2, all errors): original 96.3% found / 3.8% wrong / p95 104 nm (uncalibrated) and 95.9% / 0.1% / median 25 /
p95 61 nm (calibrated); median_peaks 97.0% / 0% / p95 83 nm and 97.9% / 0% / median 21 / p95 54 nm.
Older toy run (seed 1, bright spots): every variant finds 100%; Gaussian fit p95 64 nm vs 72 nm for the centroid; median 44 nm in all (the
uncorrected chromatic shift, ~0.3 px).

3D run (seed 22, 4 cells x 40 frames, 17 planes 0.8 um apart, all errors, fiducial beads):

| finder | analysis | calibrated | found | wrong | median | p95 |
|---|---|---|---|---|---|---|
| original | 2D on max projection | yes | 100% | 4.7% | 40 | 110 |
| **median_peaks** | 2D on max projection | yes | 99.4% | 0% | 39 | 80 |
| original | 3D | yes | 99.7% | 4.7% | 101 | 320 |
| **median_peaks** (axial median off: auto, spacing > axial sigma) | 3D | yes | 99.1% | 0.3% | 84 | 183 |
| median_peaks + axial median always (`median_z: 3`) | 3D | yes | 99.1% | 0% | 85 | 177 |

What each one-at-a-time change does (library run, calibrated, vs median_peaks 90.7% found / 0.06% wrong / 25 nm):
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
