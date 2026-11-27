# Spot finder, MS2 trace, bleaching, deconvolution, samplers, multi-copy linking, isolation

Working notes moved verbatim from `CLAUDE.md` on 2026-11-27 (section order kept). Add new dated sections at the top; delete statements that are no longer true.

## Single-cell isolation + per-cell analysis (stage 5a / 6)
- `isolate_cells(images, labels, out_dir)` writes one fixed-size (S x S, centred on the centroid) movie per cell ID:
  `stage5_cells/<source>/cells/cell_XXXX/{nucleus,locus0,locus1,mask}.tif + frames.csv`, plus `cell_table.csv` and
  `cell_summary.csv` (lifetime, gaps, border contact, nearest-neighbour distance). `mask` holds only that cell's pixels.
- Identity sources: `truth` (ground-truth IDs), `tracked_gt_masks`, `tracked_auto_masks` (tracker IDs). Running the
  same analysis on each is how tracking errors will be propagated into locus results (stage 7b).
- `map_cells(func, source_dir, n_workers)` runs `func(CellData, **kw) -> DataFrame` over all cells in a process
  pool (one cell per worker; `func` must be a top-level function). Serial and parallel results are tested identical.
- Each locus has its OWN colour channel (user's experiment: loci are different colours), so locus identity is unambiguous.
  Same-colour loci that move farther per frame than their separation cannot be told apart (we saw ~140 swaps).
- Memory (measured, loci_demo: 8 cells, 120 frames, 615x615 px, 3 channels): images stored as uint16 (`image_dtype`)
  = 87 MB per channel for the full movie; per-cell crops ~28 MB per cell per identity source on disk; each analysis
  worker peaks ~245 MB (mostly Python/libs); the MAIN process peaks ~2.1 GB during whole-movie segmentation/tracking
  (that is the biggest consumer, not the per-cell work). Dashboard reads single TIFF pages only.
  Scaling caveat: `CellData` loads a whole cell into RAM; for long 3D movies switch to lazy per-frame access.
- Dashboard performance (fixed after the user found movies slow): (1) each tab body is an `@st.fragment`, so a widget
  only re-runs its own tab (Streamlit otherwise re-runs ALL tabs on any click); (2) `load_run` uses
  `st.cache_resource` (cache_data copied the whole movie on every interaction); (3) playback is an in-browser JS
  player (`app/player.py`) fed by pre-rendered JPEG frames from `app/fastview.py` (numpy/PIL, ~44 ms/frame; matplotlib
  was far too slow), gated behind a "Prepare player" toggle and cached; expanders stay open while their toggle is on.
  Do not draw per-frame matplotlib figures for animation.
- Dashboard speed-up (2026-10-11, measured with Streamlit AppTest + a stack-sampling profiler, scripts kept out of the repo): (1) `st.tabs`
  ran ALL 14 tabs' code on every full rerun (42 s cold / 20 s warm); now a `st.segmented_control` ("page", in session state) renders only
  the open page (8.5 s / 3.1 s) and no longer jumps back to the first tab; (2) `style.fast_pyplot` replaces st.pyplot (whose default is
  200 dpi: a 17-inch figure was a 3400 px PNG) with 110 dpi + light PNG compression; (3) the Movie page draws frames with
  `fastview.movie_frame` (PIL) instead of matplotlib, with the movie-wide grey scale cached (`movie_state`; it was recomputed over the
  whole movie per frame): 3.5 -> 0.2 s; (4) tracked-cell paths cached (`tracked_paths`): Trajectories 2.6 -> 0.3 s; (5)
  `style.cached_figure(st, key, draw)` keeps static figures (run histograms, comparison bars, bleaching curves) as PNG bytes per session.
  Per page warm now 0.1-2.4 s (Spot finding and Transcription are the slowest). AppTest note: cProfile only sees the main thread;
  AppTest runs the page in another thread, so sample `sys._current_frames()` instead.
- The stand-in loci (`toy_loci.py`) are INVENTED, not from polychrom or data: confinement 0.1 x nuclear radius, relaxation
  time 120 s, 600 photons/locus/frame are guesses. Never present their statistics as physical results.
- Gotchas: tifffile treats stacks of exactly 3-4 frames as RGB, so always write with `photometric="minisblack"`;
  Cellpose needs the right `diameter_px` at fine pixel sizes (`diameter_px: auto` in loci_demo.yaml; without it DET fell
  to 0.946 at 130 nm/px); Streamlit does not reliably hot-reload helper modules (`app/cells_tab.py`), so restart the
  dashboard after editing them.

## Spot finding (stage 6; 2026-10-11; full description + numbers in docs/SPOT_FINDING.md)
- ONE method now: `dog_viterbi` (`stage6_analysis/dog_viterbi.py`; `brightest` = the old placeholder kept for comparison), after Chen, Levo, Zoller & Gregor
  2025 (Nat. Struct. Mol. Biol., doi 10.1038/s41594-025-01615-4; methods from arXiv 2304.08770v3 section III.2; no code published). It tracks the promoter
  (locus 0) and the enhancer (locus 1), each in its own colour: per-cell crops -> median 3x3 -> difference of Gaussians (1, 3 px) with a LOW threshold (3 x noise,
  up to 5 candidates) -> one Viterbi path per cell and colour (link = -log diffusion propagator; gaps up to 5 frames at 2 nats each; after LONGER gaps the path may
  pick the spot up again anywhere at log(cell area) + the miss cost; reward = log SNR + log(area / false candidates per frame); D re-estimated from the path) ->
  Gaussian fit -> ellipse intensity (mean in (7,7,5) px minus the ring to (19,19,7); ring MEDIAN by default, paper = mean). The clutter term in the reward is
  essential (without it the path skipped 547 of 1914 frames with the right candidate). Calibrated: library run found 98.0% (original 82.2%), 0 wrong picks,
  median 26 nm; toy seed 12 99.5% / 22 nm; 3D run in 3D 100% / 0% / 84 nm (original 4.7% wrong).
- The user asked for Bennabi et al. 2025 (Sci. Adv.) first; that paper finds Mesp2 poles, not spots; the rescaled `median_peaks` method built from it was REMOVED at
  the user's request (it worked, 90.7% found, but was not Bennabi's method). History + numbers in docs/SPOT_FINDING.md; code before commit 0142350.
- MS2 spot (`stage6_analysis/ms2.py`): trace = ellipse intensity in the nuclear channel in EVERY frame (interpolated positions without a detection), as Chen et
  al. Position `ms2_position: promoter` (default) = the promoter tracked in its own colour, or `own` = dog_viterbi on the nuclear channel (nucleus eroded 3 px).
  Measured (seed 12, calibrated; MS2 photons per loop per s 100 / 300 / 1000): own track located 0% / 0% / 74% of ON frames, trace r 0.07 / 0.07 / 0.37;
  at the promoter 99.5-100% located, trace r 0.01 / 0.31 / 0.82. At the PLACEHOLDER MS2 brightness (100) transcription is NOT measurable at all (one Pol II
  ~240 photons on a ~300 photons/px nucleus; DoG SNR ~2): the user must set MS2 / nuclear brightness from data. Runs: loci_demo_seed12_cfcf39 (1x),
  _026a82 (3x), _fc340a (10x). Scores `ms2`, `ms2_calibrated` in scores.json; traces stage6_analysis/<src>/ms2_traces[_cal].csv.
- Photobleaching correction (`stage6_analysis/bleach.py`): b(t) from BACKGROUND-SUBTRACTED signal (nuclei minus the empty area > 8 px away; loci: median of the
  ring-subtracted spot intensities over cells), exponential fit; `intensity_corr` = background-subtracted intensity / b(t) in loci_positions*.csv and
  ms2_traces*.csv; stage6_analysis/bleach.csv (+ true and naive curves), scores["bleaching"]. Measured: nucleus tau 18.4k s vs true 20k (naive, background left in:
  27.9k s, misses ~30% of the decay); locus0 6001 vs 6000 s; spot late/early ratio 0.86-0.89 raw -> 1.00-1.03 corrected (truth 0.97-1.03).
- Uneven illumination + per-cell bleaching (2026-10-11, docs/SPOT_FINDING.md section "Uneven illumination"):
  - **Simulation:** `imaging_errors.bleaching.per_cell` (ON in the loci configs; default off in the table) bleaches each cell by its own light dose: the flat-field at
    its centre (+ stage drift) summed over its earlier frames, 0 outside the lit area (`lit_margin_um`), daughters inherit it (`Imager.cell_bleach`). Truth:
    `bleach_cells.csv`, plus `bleach`, `bleach_nucleus` and `in_image` columns in loci_truth_image.csv. `geometry.outside_um` gives cells room outside the view
    (walls at -m and fov+m), so they leave and enter.
  - **Downstream:** scoring ignores spots outside the image; the CTC truth splits each cell's stretches in view into linked tracks (`ctc_ground_truth`;
    traccuracy rejects tracks with gaps).
  - **Analysis:** `stage6_analysis/flatfield.py` self-calibration fits log nuclear signal = a_cell + log F(x) - k D_cell(t), with a quadratic F, a ridge prior and
    4-fold cross-validation over cells (F kept only if held-out loss improves by >= 2%). Locus colours reuse the nuclear F (`locus_flat="nucleus"`): their own spot
    fit gave a 12% map error. `bleach_correct_run` now writes `intensity_corr` (illumination + per cell) and keeps the old curve as `intensity_corr_global`.
  - **Measured:**
    - Seed 31 (moving cells, entering): map error 2.1% vs 3.8% for even light (flat slide 1.0%); per-cell bleaching error 0.044 vs 0.091 for one curve.
    - Seed 12 (slow cells): the fit was rejected, correctly.
    - Locus tau stays biased (probe brightness changes over time).
- MS2 brightness sweep: `scripts/sweep_ms2_brightness.py` (R = one Pol II's MS2 photons / one nuclear pixel; placeholder R = 0.8; 6 values x 3 seeds, `tracking.skip: true`
  = new config switch that skips stages 4/7 and analyses true IDs only; such runs do not appear in the dashboard's run list). Results in
  `$SIMLIVE_DATA/sweeps/ms2_brightness/results.csv`, Spot finding page part 3d.
  Result (3 seeds each, ~1 Pol II per ON frame): trace r at the promoter 0.07 / 0.15 / 0.40 / 0.76 / 0.95 / 0.99 for R = 0.25 / 0.8 / 2.5 / 8 / 25 / 80; ON called
  from the trace 1% up to R = 0.8, 64% at 8, 97% at 25; own-track MS2 finding needs ~10x more. R must be measured on real data (no published value found).
- NOT built: using promoter + enhancer + MS2 jointly to choose spots, and separating homologous copies: plan in docs/PLAN_three_signals_and_homologs.md.
- MS2 -> Pol II deconvolution (2026-10-12, docs/SPOT_FINDING.md "From the MS2 trace to Pol II starts"; `stage6_analysis/deconvolve.py`, numba; page 🧮 MS2 → Pol II
  `app/ms2_tab.py`; `scripts/deconvolve_ms2.py --run`; run by the pipeline after the bleaching correction):
  - Chen et al. 2025 V.1/V.4: step kernel (24 loops over the cassette, release at (cassette+gene)/v + mean dwell), 1 s grid, Bernoulli prior from the
    population mean, noise sigma_b^2 + beta G from frame-to-frame differences, unit = one fully loaded Pol II computed from the microscope (unit_check vs
    truth 0.94-0.96).
  - DEVIATION: their adaptive block sampler (proposal floor 0.05) cannot propose sparse starts; we Gibbs-sample single steps + shift moves (same
    posterior).
  - Measured (R=25): placeholder gene (seed 7, 47132d) 5 true / 5.9 inferred starts, all found; denser gene (k_init 0.05/s, seed 8, f5c88d)
    389 / 401, 57% within 30 s, 0.45 false per true, 1-min rate r 0.80; Chen's burst rule (2 starts in 50 s) catches 49% of ON frames (0% for
    the placeholder gene).
  - Likely limit: fixed release time vs the simulation's random dwell (untested).
- Image errors -> results, per source (2026-10-12; `stage7_validation/error_effects.py`, `scripts/error_effects.py --run`, Images page part 2
  `app/errors_view.py`):
  - **What it does:** re-renders the run's movie for a reference microscope (PSF + QE + shot + read noise), the reference + ONE source, and all sources.
    It analyses each with the source's correction off / on (FIXES: analysis countermeasure, calibration images via `Calibration.restricted(parts)`, or
    the movie alone) and scores found / wrong / position / distance rms / intensity error / late-early / MS2 r.
  - **Outputs:** `stage7_validation/error_effects.csv` (resumable); movie calibrations kept in `error_effects_selfcal/`. About 2 min per source.
    Only 2D runs.
  - **Movie-only corrections** (`stage6_analysis/selfcal.py` builds a normal `Calibration` from the movie):
    - background/offset map = per-pixel median over nucleus-free frames;
    - hot pixels by a per-pixel test (excess > 6 SE of the pixel's median AND > 5 x map spread, >= 10 free frames). The first version flagged
      11,000 pixels because the MAD of integer camera values is 0; it now finds 100 of 117, 0 false;
    - flat = nuclear self-calibration;
    - colour shift = mean enhancer-promoter vector, applied only if > 2 SE (cell bootstrap; 8-12 cells are NOT enough: SE ~0.4-0.9 px).
  - **Other new options:** `localize_loci(..., row_correct=True)` (row offsets from pixels away from the cell); `intensity_background: none`.
    - cross-talk unmixing: nuclear -> locus bleed-through fraction by robust regression inside nuclei (measured 1.997% vs 2.000%), subtracted in
      `channel_inputs` when the calibration carries values["crosstalk"].
  - `Calibration.restricted(parts)` splits the flat-field slide into "flat" (smooth illumination) and "prnu" (pixel part) and drops the
    hot-pixel excess from "offset" unless "hot" is chosen. Applying the whole slide to a single-source test movie added errors the movie did
    not have.
  - Results (docs/ERROR_EFFECTS.md, showcase run loci_demo_seed8_f5c88d, one seed):
    - **Colour shift dominates positions:** 12 -> 43 nm, distance rms 16 -> 63 nm; beads fix it. The movie estimate is not significant with
      12 cells (needs ~30 to reach 2 SE).
    - **Cross-talk of the bright MS2 spot into the promoter colour dominates intensities:** r 0.95 -> 0.52. Movie unmixing gives r 0.91 but
      positions get noisier (leaked photon noise stays).
    - **Bleaching:** late/early 0.80 -> 0.99.
    - **All sources:** images 28 nm / r 0.47; movie 49 nm / r 0.80; both 26 nm / r 0.80; reference 12 nm / r 0.95.
  - Showcase runs (R = 25, cells moving): loci_demo_seed7_47132d (placeholder gene) and loci_demo_seed8_f5c88d (denser gene, error_effects run on it).
- 2026-10-13 round (user's questions: which optics? cytoplasm? how is the posterior sampled? start probability? burst rule?):
  - **Instrument presets:** `configs/instruments/<name>.yaml`, applied with `instrument: <name>` / `--set instrument=...` (deep-merged by
    `io.runs.load_config`, --set still wins; sidebar "Microscope" radio). `vutara_vxl` = Bruker Vutara VXL in widefield ("bypass") mode; `vutara_vxl_biplane` = its biplane 3D mode (see "Vutara biplane mode"). Every value is
    sourced or ASSUMED (docs/INSTRUMENTS.md):
    - 60x/1.3 NA silicone, ~100 nm px, Orca Fusion BT QE 0.95, top-hat (vignette 3% ASSUMED), active focus (defocus 30 nm), sequential colours;
    - bleed-through FROM SPECTRA (FPbase EGFP: 0.0042 cross-excitation at 555 nm x 0.062/0.717 in-band emission = 0.036%, vs ~8.6% if excited
      simultaneously);
    - biplane (2 planes 600 nm apart, 50:50): built 2026-10-17 (`vutara_vxl_biplane`, see "Vutara biplane mode").
  - **New error sources** (errors.INDEX frozen; `_LATER` list gives new sources the next numbers):
    - `sequential_colour`: later colours drawn at t + delay by interpolation + a Brownian bridge of the locus motion RELATIVE to its nucleus;
      `Imager.locus_positions_at_exposure`, 2D only;
    - `cytoplasm`: cell body = nucleus outline x radius_factor, nuclear channel fraction + locus-colour photons; 2D only; ON in loci_demo/library.
  - **Movie background map:** default `stat="low"` = 20th percentile of nucleus-free frames + 0.84 x the pixel's frame-to-frame noise. The
    median is raised wherever a cell body lingers. Pixels a cell body covers almost always cannot be recovered (the shell background still
    handles spots). `error_effects.true_background()` scores maps against the object-free truth.
  - **Deconvolution:**
    - `sampler="chen_block"` (the paper's adaptive block MH, block auto-tuned to >= 10% acceptance: ~20 s here; the paper's 6 min gives 1.7%
      acceptance and a stuck chain) agrees with gibbs_shift;
    - per-second start probability `p_fine` saved;
    - SCORING FIXED: per posterior sample, one-to-one matching. The old cluster-peak scorer gave 57%; the true value is 82% within 30 s;
    - `dwell_fixed` transcription option: a fixed wait gives 93% / rate r 0.95 vs 82% / 0.80 with the random wait, so the kernel's fixed
      release time is a real limit;
    - burst-rule map (page step 6): on the denser gene w = 30 frames, g_b = 6 gives ON/OFF 94/94% vs Chen default 49/99%.
  - **Sampler convergence (2026-10-13; docs/MCMC_CONVERGENCE.md; `stage6_analysis/mcmc_check.py`, `scripts/check_mcmc.py`, page step 5b):**
    - 4 chains per trace from different starts (`deconvolve.INITS`), rank-normalised split R-hat + bulk ESS on log-posterior, number of starts and
      per-frame counts; limits R-hat < 1.05, ESS >= 400 fixed in advance (do not tune). New sampler `gibbs` = flips only; `deconvolve_trace(init=)`
      returns `log_post`, `n_starts`, `log_post_sweeps`; `prepare_traces()` split out of `deconvolve_run_traces`.
    - All 3 samplers match the exactly enumerated posterior of a 14-step problem (test), so the target is right.
    - **FINDING: the default `gibbs_shift` does NOT converge on the dense, nearly noise-free gene** (seed 8, R = 25): 3 of 12 cells pass; chains trap in
      different local arrangements of starts (log-posterior levels ~100 nats apart), agree on the NUMBER of starts, disagree on timing (1-1.4 starts/min
      in the 1-min rate). More burn-in does not help. Sparse placeholder gene (seed 7): 9 of 12 pass, the rest only miss ESS 400 (320-353). Constructed traces: 5/5.
    - Flips only get stuck (1/5): the shift moves matter. Block sampler 2/5 (needs more samples; 5x cheaper per sample).
    - **FIX BUILT: `sampler="gibbs_pt"` (parallel tempering, 16 replicas, likelihood strength 1 -> 1e-3, adaptive ladder, `n_replicas`, `beta_min`)**:
      constructed 5/5, sparse run 12/12, dense run 11/12 (the miss: one frame's R-hat 1.056 vs 1.05, limit not moved); chain disagreement 0.09-0.26
      starts/min (was 0.12-1.39). ~20-30x slower than gibbs_shift. **DEFAULT NOW `sampler="auto"` (user's choice 2026-10-13)**: per cell, 4 gibbs_shift chains + the check; pass ->
      chain 0; only ESS short -> same sampler with 3x samples; R-hat fails -> gibbs_pt (info["checks"], info["n_tempered"]; dense run 9/12 tempered, sparse 0/12).
      Force one sampler with `scripts/deconvolve_ms2.py --sampler gibbs_pt|gibbs_shift|chen_block` (files `_<name>`; the page's sampler table lists them) and check
      with `scripts/check_mcmc.py --samplers gibbs_pt --tag _pt --workers 6`. Headline scores vs the simulation are unchanged by tempering (found 0.83, false 0.20,
      rate r 0.80; per-frame posterior means differ by about as much as the other samplers differ from each other), so earlier totals stand; intervals and
      single-start timing of dense bursts should come from the tempered chains. Open: fewer replicas / warmer floor than 1e-3 (not explored).
  - Vutara showcase: loci_demo_seed9_1c649d (sequential, tracked) and loci_demo_seed9_17363e (same movie, simultaneous: 8.6% bleed-through,
    no colour delay; tracking skipped). Error effects for both in docs/ERROR_EFFECTS.md.
- **BURSTY DEFAULT GENE** (user 2026-10-13: "the transcription trajectory needs to be way more bursty"). Configs + transcription.DEFAULTS +
  the Transcription page now use:
  - k_on 1/900 s (OFF 15 min), k_off 1/120 s (ON 2 min), k_init 0.2/s → ~24 Pol II per burst, 3.5 bursts/h, ON 12%;
  - the trace spikes and returns to zero. The old TFF1-like values (~1.6 Pol II per burst) are kept in the config comments.
  - Showcase loci_demo_seed10_ebf264 (Vutara, 40 min, PLACEHOLDER MS2 brightness R = 0.8): bursts are measurable even at the placeholder
    brightness. Deconvolution: 86% of starts within 30 s, rate r 0.91; Chen's default burst rule catches 97% of ON frames (92% agreement).
  - The MS2 page now opens on a typical cell (median number of starts).
  - SECOND STEP, same day ("still a ton of accumulated MS2 signal, even more bursty"), now the default:
    - k_on 1/1200 s (OFF 20 min), k_off 1/45 s (ON 45 s), k_init 0.15/s → ~7 Pol II per burst;
    - shorter residence: gene_kb 2.0 (was 5.0), dwell_s 20 (was 116, Coulon 2014's measured value); a Pol II stays ~100 s instead of ~270 s;
    - ~3 bursts/h. The trace is narrow spikes of 4-7 Pol II.
    - check_sampling now WARNS at 10 s frames (only ~4-5 frames per 45 s burst): use faster frames for real data of such a gene.
    - Showcase loci_demo_seed11_24bba2 (Vutara, placeholder MS2): trace r 0.81. Deconvolution: 200 true / 293 inferred starts, 79% within
      30 s, 0.68 false per true (weak spikes in noise ±2 Pol II), rate r 0.88; burst calls 97% agreement, 84% of ON frames.
- scripts/compare_localizers.py --run <run> also computes MS2 traces and the bleaching correction (`--ms2-only` = just those; for older runs). Variants in
  compare.py: original, default, chen_literal, no_median, dog_snr5, integer (+3D: zmed3); tables are merged, so delete method_comparison.csv after renaming variants.
- Scores (`score_localization[_3d]`): detection over loci WITH light (`photons` in loci_truth_image.csv), `n_false_on_dark`, `frac_far` (> 0.5 um), `rms_error_near_um`,
  `intensity_r_locusK` (vs `photons_frame` = photons x bleach).
- Logic errors fixed 2026-10-11: u/v from the rounded crop centre; `min_sep_px` ignored; loci without light counted as misses; single-cell player drew other cells'
  true loci and put markers half a pixel off; centroid pixel-locking; Spot finding per-cell stats counted MS2 OFF frames as wrong picks.
- Dashboard 🔎 Spot finding page: method table, every intermediate image for any cell / spot (locus 0, locus 1, MS2) / frame with editable settings, per-frame
  errors, every candidate's distance to the truth, intensity vs true photons, MS2 trace (promoter vs own track vs truth), run scores, 3b photobleaching (true /
  estimated / naive curves, raw vs corrected), method comparison. Run summary has "Image errors at a glance" (errors.BRIEF: one line per source).
- Gotchas: scripts that call `map_cells` with workers need `if __name__ == "__main__":` (Windows spawn). Keep CLAUDE.md UTF-8 (one stray cp1252 byte was fixed).

## Several promoter-enhancer copies per cell, stray spots, figure audit (2026-10-19/20; user: "multiple sets of E-P pairs and separate bursting dynamics", "stray spots ... the linking algorithms have to distinguish", "make sure the graphs are informative", "fix any problems", "update the tutorial")
- **Simulation:** `loci.copies` (`stage1_chromatin/copies.py`): copies per cell (`n`, `weights`), `min_separation_um`, `ms2: all|first`, per-copy
  transcription overrides; toy + library loci (extra copies = other library simulations, own pair/rotation/offset/territory), probes and transcription
  per copy (`copy` column everywhere, `trace_id` = cell_id + 100000 x copy). Copy 0 uses exactly the old random numbers: with copies and strays off,
  every truth file, image and analysis output was byte-identical to commit 06e4070 (regression run against a git worktree at a SHORT path:
  numba cache files exceed Windows' path limit under the scratchpad). Pipeline seed streams: `SeedSequence(seed).spawn(10)`, [8] = copies.
  New error source `stray_spots` (`Imager.stray_spots`, truth `stage3_microscopy/stray_spots.csv`, in `_LATER`). Configs loci_demo/loci_library: half
  the cells with two copies, copy 1 bursts twice as often (k_on x2), strays 1/min, `max_copies: 2`, `max_candidates: 8`, `pair_prior: gauss`.
- **Analysis:** `stage6_analysis/multicopy.py` (paths per colour, locus test coverage/median SNR/D relative to the best path, `untangle`, `assign_pairs`
  that maximises valid pairs first, pair prior gauss|gate|none, joint untangle of copies); 2D (`localize.analyse_cell_copies`), 3D (`localize3d`), MS2
  per copy (`ms2.relabel_to_true_copies` for the truth source). Scores: `localize.match_copies` (EXTRA_COPY = 50 + found number for unmatched copies),
  `frac_on_stray`, `frac_on_other_copy`, `copy_number_right`, `missed_copies`, `extra_copies`. Truth-source positions are renumbered like the true
  copies before saving (`copy_found` keeps the analysis' number). Deconvolution keyed by `trace_id` (`deconvolve.trace_key`). `locus_analysis.deconvolve:
  false` skips it (sweeps).
- **Measured** (25 movies, docs/SPOT_FINDING.md "Several copies per cell"): one path per colour gets the copy number right in 58% of cells; the linker
  in 100% (0-0.5 strays/nucleus/min), 98% (1-2), 90% (4). Smooth prior: fewest wrong picks (0.2% -> 8% from 0 to 4 strays/min), 0% other-copy picks
  without strays, distance bias 3-30 nm; flat gate cannot separate copies closer than ~1 um. Bias check: the smooth prior drops ~2% of frames, mostly far
  ones (censoring of the tail). Showcase `loci_library_seed12_87a16d` (tracked, calibrated, R = 8): 10/10 copies, 98.6% found, 0.6% wrong picks.
- **Fixed:** provenance cp1252 crash (`io/runs._git`), skimage `max_size`, bead tracking outliers (`fiducials.consistent_beads`: a bead whose height fit
  wandered by um made frame medians jump up to 4.6 um; now within 0.2 um), MS2 deconvolution noise model (sigma_b had collapsed to its floor -> phantom
  Pol II; now measured on the low-signal bins: false starts per true 1.57 -> 1.17 on the showcase, unchanged on seed8/seed11). Per-trace baseline
  (`baseline="auto"`) exists but is OFF (helped one run, hurt two).
- **Dashboard:** Cells (P/E + copy labels, stray squares, figure = position errors with stray/other-copy marks + found-minus-true offsets + one distance
  panel per copy), Spot finding (copy selector, part 2b "what the linker did", 3e sweep, error histogram in nm, MS2 spot starts on the bursting copy),
  Images (stray squares, copy labels, per-copy truth table), MS2 (one trace per cell and copy), Transcription (contact zoom with %, P(ON) in contact vs far
  in the step-4 title, tiles start at the first burst, run view per copy), Trajectories (crowding vs speed at switches), Scans (numeric x axis, no
  lines between categories), Calibration, Probes, Chromatin (um / minutes axes), Stage; home summary lists copies and strays; sidebar: copies share,
  copy 1 rate, MS2 on all copies, stray rate/brightness, copies looked for, pair prior. Every page's figures audited with `scripts/audit_figures.py`
  (AppTest, saves each figure); changed pages also checked in a browser (port 8502).

## Locus brightness in multi-copy linking (2026-10-26; user: "utilize the intensity of the individual loci ... to make the linking better")
- Built: `locus_analysis.viterbi_intensity` (spots.DEFAULTS, default **false**): `dog_viterbi.IntensityModel` = AR(1) model of the path's DoG SNR
  (level + bleaching trend, spread, lag-1 correlation, fitted to the path like D); each spot pays (z^2 - 1) / 2 nats on top of the usual costs.
  Tests `tests/test_viterbi_intensity.py`. `scripts/sweep_copies_strays.py` gained `--variants` (extra: multi_int, multi_gauss_int) and
  `--reanalyse TAG` (saved movies again -> `<TAG>_<variant>.csv`).
- **Result: no net gain, so it stays off** (docs/SPOT_FINDING.md "Locus brightness in the linking"): 25 movies, placeholder probes: found 98 -> 95-98%,
  stray picks slightly lower at 2-4 / min, other-copy picks up at 0.5 and 4 / min, 2 copies lost; slow probes (rates / 10, true lag-1 r 0.9): no gain
  either. Cause: the per-frame MEASURED brightness is noisy (lag-1 r 0.43-0.88 vs true 0.83-0.92, CV ~0.3), copies are equally bright (same probes),
  strays are 0.68 x a locus. A first likelihood-ratio version (vs the leftover candidates) was worse (the other copy counts as clutter) and was dropped.
- Worth re-testing with the user's real numbers: locus photons per frame, probe residence time, stray brightness, and whether copies are labelled
  differently. Not tried: a better per-candidate brightness (fitted amplitude instead of the DoG peak); brightness in `untangle` / `untangle_by_place`.
