# Spot finding (stage 6): method, sources, measured performance

Stage 6 finds the promoter (locus 0) and the enhancer (locus 1), one spot per colour, in every frame of every isolated cell, and measures the
MS2 transcription signal at the promoter. It is still a placeholder for the user's own analysis, but a defensible one: it follows a published
pipeline, every step can be inspected in the dashboard (🔎 Spot finding tab), and every choice below was scored against the simulation's truth.

| `locus_analysis.method` | source | code |
|---|---|---|
| `dog_viterbi` (**default**) | Chen, Levo, Zoller & Gregor 2025, *Nat. Struct. Mol. Biol.* (MS2 transcription spots) | `stage6_analysis/dog_viterbi.py` |
| `brightest` | the original placeholder | `localize.py` |

Code around them: `localize.py` (2D, per cell), `localize3d.py` (z-stacks), `ms2.py` (MS2 trace), `bleach.py` (photobleaching correction),
`compare.py` + `scripts/compare_localizers.py` (method comparison; also computes MS2 traces and the bleaching correction for older runs).
Every detection also gets Chen et al.'s ellipse **intensity**, scored against the photons the simulation drew.

# The method: `dog_viterbi`, after Chen et al. 2025

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
After a gap longer than `viterbi_max_gap` frames the path may pick the spot up again anywhere in the cell, at the cost of an unconstrained
position (log of the cell area) plus the per-frame miss cost; so one path can hold several bursts separated by long OFF periods (before this,
a gap longer than 5 frames ended the track: fine for loci, which are always there, fatal for MS2 spots).

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

Open: D per cell rather than pooled; DoG widths and `dog_snr` set by reasoning, not tuned; the paper's two-photon PSF and our widefield
PSF differ in axial shape. (The `median_peaks` rows are the removed second method; see History.)

# The MS2 spot (`ms2.py`)

The MS2 spot is drawn at the promoter, in the nuclear channel (same colour as the NLS-GFP). Its trace is the ellipse intensity (spot minus
ring, removing the nuclear MCP/NLS-GFP background) in EVERY frame at a tracked position, as in Chen et al.; positions in frames without a
detection are interpolated in time. Two ways to get the position (`ms2_trace(..., ms2_position=)`):

- `promoter` (default): the promoter (locus 0), tracked in its own colour with dog_viterbi (the measured chromatic shift removed when
  calibrated). In the user's experiment the promoter is labelled, so this is available and reliable whether the gene is ON or OFF.
- `own`: dog_viterbi on the nuclear channel itself (nucleus eroded by 3 px so its bright edge gives no candidates), as Chen et al. do with no
  other label on the gene.

Measured, seed-12 toy movie (probes + MS2, all errors, calibrated) at three MS2 brightnesses (photons per loop per second; placeholder 100):

| MS2 brightness | own track: ON frames located (<= 0.5 um) / trace r | at the promoter: located / trace r |
|---|---|---|
| 1x (100) | 0% / 0.07 | 99.5% / 0.01 |
| 3x (300) | 0% / 0.07 | 100% / 0.31 |
| 10x (1000) | 74% / 0.37 | 100% / 0.82 |

At the placeholder brightness the MS2 signal is not measurable even at the exact right place (r = 0.01): one Pol II carries ~240 photons per
exposure on a nucleus of ~300 photons per pixel, and in ON frames the true spot's DoG signal-to-noise is ~2 (median rank 29 among the
nucleus' local maxima). The own-track mode then follows a fixed texture feature of the nucleus. The MS2 brightness per loop and the nuclear
MCP/NLS-GFP brightness are placeholders the user should set from data; the result above says transcription becomes measurable at roughly
3-10x the current MS2 brightness (or correspondingly less nuclear background).

**Brightness sweep** (`scripts/sweep_ms2_brightness.py`, 2026-10-11).
- **R** = photons of one Pol II carrying all 24 loops / photons of one nuclear-GFP pixel, same exposure. Placeholder R = 0.8.
- **Setup:** toy loci with probes and transcription on, all image errors, calibration applied, true cell IDs. 3 seeds per value. Mean (min..max).
- **The ON signal:** the gene was ON in 12% of cell-frames (6-19% by seed), with a median of 24 loops (about one Pol II) when ON.
- **What a real movie still needs:** no published ratio of spot to nuclear background was found (2026-10-11 search), so R has to be measured on a
  real movie. Take the background-subtracted spot photons per Pol II over the mean nuclear pixel value, with the camera offset removed.

| R | trace r, at the promoter | ON called from that trace (sensitivity / specificity) | trace r, MS2 spot's own track | ON frames located by its own track |
|---|---|---|---|---|
| 0.25 | 0.07 (-0.17..0.31) | 0.01 / 1.00 | 0.02 | 0% |
| 0.8 (placeholder) | 0.15 (-0.04..0.34) | 0.01 / 1.00 | 0.03 | 0% |
| 2.5 | 0.40 (0.33..0.46) | 0.11 / 1.00 | 0.00 | 5% |
| 8 | 0.76 (0.71..0.86) | 0.64 / 1.00 | 0.30 | 21% |
| 25 | 0.95 (0.93..0.96) | 0.97 / 1.00 | 0.57 | 38% |
| 80 | 0.99 (0.99..0.99) | 0.99 / 1.00 | 0.83 | 69% |

What this means:
- **Reading at the promoter:**
  - single Pol II: needs R of about 8 to be followed (r 0.76) and about 25 to be called reliably ON/OFF (97%);
  - at R <= 0.8 it cannot be seen at all;
  - bursts with several Pol II scale accordingly (n Pol II x R).
- **Finding the MS2 spot by itself** needs roughly 10x more brightness than reading it at the tracked promoter, and is unreliable even at R = 80.
Scores in `scores.json` (`ms2`, `ms2_calibrated`): `found_on`, `median_error_um`, `intensity_r`, ON/OFF called from the trace alone
(`call_sensitivity`, `call_specificity`; threshold 3 x the trace's robust noise).

# Photobleaching correction (`bleach.py`)

A pixel records offset + background + signal x b(t); only the fluorophores bleach. So:
1. b(t) is measured on background-subtracted signal: nuclear channel = mean inside the nuclei minus the mean of the empty area > 8 px away from
   them; locus colours = the spots' own ring-subtracted intensities, median over all cells per frame. Then an exponential fit.
2. The correction divides only the background-subtracted spot intensity by b(t) (`intensity_corr` in loci_positions*.csv and ms2_traces*.csv).

Measured against the true bleaching (10x-MS2 seed-12 run): nucleus tau 18,400 s estimated vs 20,000 true (b at the end 0.937 vs 0.942); the
naive estimate on raw brightness (background left in) gives tau 27,900 s and misses ~30% of the decay. Locus 0: tau 6,001 vs 6,000; locus 1:
7,180 vs 6,000. The spot intensities' late / early ratio goes 0.86-0.89 (bleaching) -> 1.00-1.03 after correction, matching the true photons
(0.97-1.03). In these movies bleaching is mild (6% nucleus, 18% loci over 20 min); the same code handles stronger bleaching.
`tests/test_spots.py::test_bleaching_is_measured_after_removing_the_background...` checks it on a synthetic movie with a known tau.
That single curve is kept as `intensity_corr_global`. `intensity_corr` is now the per-cell correction below.

## Uneven illumination + per-cell bleaching, without a flat-field image (`flatfield.py`, built 2026-10-11)

**The simulation** (`imaging_errors.bleaching.per_cell: true`, on in the loci configs):
- Each cell bleaches by the light it has received: the illumination (flat-field map) at its position, summed over its earlier frames.
- It receives nothing while outside the lit area (the camera's field + `lit_margin_um`).
- Daughters inherit the mother's dose. tau refers to the field's mean illumination.
- `geometry.outside_um` lets cells live outside the imaged field, so cells leave and unbleached cells walk in.
- Truth: `stage3_microscopy/bleach_cells.csv`, plus `bleach` / `bleach_nucleus` / `in_image` columns in loci_truth_image.csv. Spots outside
  the image are not counted as misses.
- The CTC tracking truth writes each stretch a cell spends in view as its own track, linked to the previous stretch.

**The estimate** (self-calibration on the moving cells):
- **Model:** log(background-subtracted nuclear brightness) = log a_cell + log F(position) - k x D_cell(t), where D_cell = the sum of F along
  the cell's own track before frame t. F is a smooth quadratic surface with mean 1.
- **Fit:**
  - Each cell's brightness a_cell is removed exactly (each cell's mean is subtracted).
  - Robust least squares, with a prior of sd 0.1 on each coefficient of F.
  - The illumination term is kept only if it lowers the robust loss of held-out cells by at least 2% (4-fold cross-validation over cells).
- **Locus colours:** use the nuclear channel's illumination shape (same light path) and fit only their own bleaching rate (`locus_flat="nucleus"`;
  `"own"` / `"none"` also exist).
- **What it cannot know:** light received before a cell is first seen. That is absorbed into a_cell, so a cell's correction is relative to
  its first frame in view. It is absolute when light falls only on the imaged field.
- **Correction:** intensity_corr = intensity / (F(spot) x b_cell(t)). Tables made with the calibration applied are already flat-field corrected,
  so only b_cell is divided out there.
- **Outputs:** `stage6_analysis/<src>/nuclear_signals.csv`, `flatfield_est.npz`, `scores.json["bleaching"]["flatfield"]`.
- **Tests:** `tests/test_flatfield.py` (synthetic cells with a known map and tau; renderer dose for a cell entering at frame 10).

**Measured (one seed each; placeholder flat-field: vignette 0.15, tilt 0.03, i.e. 3.8% rms over the field):**

| movie | illumination map error (estimate / assume even / flat-field slide) | tau nucleus (est / true) | nuclear brightness spread within a cell (raw / one curve / new) | per-cell bleaching rms error (new / one curve) |
|---|---|---|---|---|
| seed 12, slow cells (~5 um over the movie), all in view | fit rejected by cross-validation: even light used, 3.8% / 3.8% / 1.0% | 18,300 / 20,000 s | 2.5% / 1.6% / 1.6% | (global bleaching in this run) |
| seed 31, speed 1.0, drift 0.03 um/s, 20 um outside, 16 cells, per-cell bleaching, tau 4000 s | 2.1% / 3.8% / 1.0% | 3,680 / 4,000 s (one curve: 3,170) | 6.6% / 2.1% / 1.5% | 0.044 / 0.091 |

Locus colours (seed 31):
- **Spot intensity vs true photons:** r = 0.784 / 0.801 (raw 0.734 / 0.766, one curve 0.777 / 0.790).
- **Bleaching rate:** their tau stays biased (about 4,200 vs 3,000 s) because spot brightness also changes through probe binding.
- **Locus 0's own illumination fit:** map error 12%, worse than doing nothing. That is why the nuclear shape is used for the loci.

Lesson: the illumination can only be learnt where the cells went, and cells must cross a real part of the field. Slow cells or few
cells give nothing, and the cross-validation then falls back to even light, as it should. Pooling several movies of the same instrument
(one shared F) is the cheap way to make this work. A dark + flat-field slide, which the calibration simulates, is still 2x more accurate.
The standard alternative from still images is BaSiC (Peng et al. 2017, Nat. Commun.).

# From the MS2 trace to Pol II starts: deconvolution (`deconvolve.py`, built 2026-10-12)

Sources: Chen, Levo, Zoller & Gregor 2025 (Nat. Struct. Mol. Biol.). Methods are in arXiv 2304.08770v3, read 2026-10-12:
- IV.3: calibration to C.U.;
- IV.4: noise;
- V.1: deconvolution;
- V.4: burst calling.

1. **Units.** The trace is the bleaching- and illumination-corrected spot intensity at the tracked promoter (`intensity_corr`, calibrated
   traces when available). It is divided by the intensity of one Pol II carrying all its loops.
   - The paper calibrates this unit against smFISH. Here it is computed from the microscope: photons per loop x loops x exposure x QE x
     gain x the PSF fraction inside the 7 px ellipse / the ellipse's pixel count.
   - Check against the truth: the unit is 4-6% too large (`unit_check` 0.94 and 0.96).
2. **Model** (their Eqs. 8-11): trace = kernel * starts + noise.
   - The kernel climbs in 24 steps while Pol II copies the cassette. It is then 1 until release at (cassette + gene) / v + mean dwell.
     The paper adds ~800 bp for retention.
   - Starts are 0/1 per 1 s step.
3. **Noise**: sigma^2 = sigma_b^2 + beta x signal, from half the squared frame-to-frame difference binned by level (their correlation
   estimate).
4. **Prior** (Eq. 15): Bernoulli with p = R(t) x 1 s, where R = population mean activity at t + tau/2, divided by tau (tau = kernel area).
5. **Sampling.** Their adaptive block-Bernoulli Metropolis-Hastings has a proposal floor of 0.05. That suits fly genes firing ~10-20 times
   a minute, but cannot propose the placeholder gene's ~1 start per 10 min.
   - The same posterior is sampled with single-step Gibbs updates plus moves that shift a start in time.
   - 500 burn-in sweeps, then 1000 samples, every 3rd sweep.
6. **Rate and bursts** (Eq. 18, V.4): starts per frame interval for every sample.
   - Bursts: the rate is smoothed with exp(-(2t/w)^4), w = 5 frames, and called ON where it reaches g_b / w (g_b = 2), on every sample.
     This gives P(ON).
   - The dashboard lets w and g_b be changed and shows the agreement with the true promoter state.

Outputs: `stage6_analysis/truth/ms2_deconvolution.csv` (per cell and frame: trace in units, fit, expected starts, sd, rate per min, P(ON)),
plus `_samples.npz` and `_info.json` (unit, noise, scores). The pipeline runs it after the bleaching correction; for an existing run use
`scripts/deconvolve_ms2.py --run <run>`. Tests: `tests/test_deconvolve.py` (kernel, noise model, recovery of known starts). Dashboard
page: 🧮 MS2 → Pol II.

**Measured** (calibration applied, R = 25, i.e. 30x the placeholder MS2 brightness; one seed each):

| Movie | Starts (true / inferred) | True starts with an inferred one within 30 s | False starts per true one | 1-minute rate correlation | Burst calls (Chen rule: agree / ON frames caught) |
|---|---|---|---|---|---|
| Placeholder gene: ~1.6 Pol II per 16-min burst, seed 7, 12 cells, 30 min | 5 / 5.9 | 100% | 0.2 | 0.76 | 75% / 0% (never reaches 2 starts in 50 s) |
| Denser gene: 3 starts/min while ON, 10-min bursts, seed 8 | 389 / 401 | 57% | 0.45 | 0.80 | 81% / 49% |

- **Counting is right; exact timing is weak when many Pol II overlap.** The total is within 3%, but only 57% of starts are placed within
  30 s. The likely cause is the kernel's FIXED release time: the simulation's wait at the end of the gene is random (exponential, mean
  116 s). A Pol II leaving early or late is then explained by moving or adding starts (not yet tested by simulating a fixed wait).
- **The burst rule depends on the gene.** It was designed for fly genes that start Pol II every few seconds, so it misses half the ON
  frames at 3 starts/min and all of them at the placeholder rate. Widening w or lowering g_b helps; the dashboard shows the trade-off.
- **Below R ~ 8 (the brightness sweep) single Pol II are not visible**, and the deconvolution then fits noise. The page warns about it.

# History: `median_peaks` (removed 2026-10-11)

The user first asked for the approach of Bennabi et al. 2025 (*Sci. Adv.* 11, eadv7790; code github.com/gregorlab/SizeProject2024). That paper
finds Mesp2 expression poles hundreds of pixels across, not spots, so its steps (median 7x7, Otsu movie threshold, area, Gaussian sigma 20,
peak_local_max) had to be rescaled and extended; the result (`median_peaks`: median 3x3, local background, threshold on the smoothed image,
choice consistent with neighbouring frames, Gaussian fit) worked (library run, calibrated: 90.7% found vs 82.2% original) but was not
Bennabi's method, and the user had it removed. The paper's pipeline taken literally found every spot but picked the wrong object 15% of the
time (Otsu rejects nothing for a single spot). Code: git history before commit 0142350.

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
