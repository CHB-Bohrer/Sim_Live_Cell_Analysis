# Logic audit (2026-10-26) and improvement pass (2026-11-09/10)

Working notes moved verbatim from `CLAUDE.md` on 2026-11-27 (section order kept). Add new dated sections at the top; delete statements that are no longer true.

## Logic audit of the whole analysis (2026-10-26; user: "go through all logic of the live cell analysis and double check everything")
- **Read and checked** (math and bookkeeping, not only "runs"): pipeline wiring and seed streams; cell motion, placement, z; toy / library loci, copies; probe
  Gillespie + stationary theory; transcription telegraph + Pol II integrals + `expected_on`/`expected_ms2`; renderer coordinate maps (nucleus frame <-> lab,
  drift, chromatic shift, truth_in_image), per-cell bleaching dose, stray spots, camera chain; isolation crops; DoG candidates, Viterbi costs (propagator,
  miss / re-acquisition / clutter terms), D estimator (median chi-square), Gaussian fit; multi-copy linker (untangle, untangle_by_place, repair, pair priors:
  Rayleigh 1.18 / Maxwell 1.54 medians correct); 3D localizer; MS2 trace; bleaching + flat-field self-calibration (dose cumsum); calibration estimators (row
  noise pi/2 factor, DSNU, PTC, beads); fiducials; identity matching, stage 7b propagation, CTC gap splitting; E-P kernels (noncentral chi 3D, Rice 2D, shift
  average), EIV fit, history model; polymer calibration; LEF 1D engine. All consistent except the six fixes below.
- **FIXED (each with a regression test that fails on the old code where testable):**
  1. `deconvolve._sweep` shift move was NOT posterior-invariant: only an event AT the scanned step was proposed a move (the reverse move of a pair never
     proposed from the same step; events moved forward are revisited, moved back are not). Alone its stationary distribution was up to 71% off (exact
     transition-matrix check); in the real sampler (flips hide most of it) start probabilities were off by ~0.01 on the enumerable test problem (6 SE),
     and by 0.17 (a 2:1 split instead of 50:50) where only shifts can mix (`tests/test_mcmc_check.py::test_shift_moves_keep_the_posterior_where_only_they_can_mix`).
     Now a symmetric pair move (if exactly one of steps j, j + d holds an event, propose to move it). Affects gibbs_shift, gibbs_pt and "auto". Re-scored
     in memory on the showcase loci_demo_seed11_24bba2 (default spiky gene, 12 cells, auto, none tempered): found 0.786 -> 0.784, false per true 0.677 ->
     0.676, 1-min rate r 0.880 -> 0.881, P(ON) accuracy 0.973 -> 0.975, i.e. the saved totals stand for sparse traces. Re-checked 2026-11-09 (see "Improvement pass"); at the time NOT re-checked: the dense-gene runs
     (seed 8, where shifts carry the mixing), posterior intervals, and the docs/MCMC_CONVERGENCE.md tables (re-run scripts/check_mcmc.py).
  2. `library_loci`: the whole polymer was mapped through the nucleus frame, i.e. stretched by sqrt(aspect) / 1/sqrt(aspect) (and z by z_semi / radius),
     so the drawn promoter-enhancer distance was 0.89-1.14 x (5-95%; 0.82-1.24 extremes) the polymer distance that drives transcription (measured on
     loci_library_seed12_87a16d), contradicting the module's own "TRUE distances" spec. Now only the territory centre deforms with the nucleus; the polymer
     keeps its nm size and turns with the nucleus (`rigid_factors`). Test with elongated, rotated, flattened nuclei. New library movies differ from old ones
     (old runs keep the stretched geometry; their scores against their own truth stay valid).
  3. E-P vector model (`ep_acquisition.measure`, `ep_polymer.simulate_population`): the polymer block was picked by round(t / block) while the distance that
     drives the promoter uses floor (block k holds for [k, k+1) x block_s): a 0.23 s late measurement on average, while the E-P vector moves ~200 nm per axis
     per block. Excess of the measured-vs-true error over the model's own sigma: 2D +32 -> +6 nm, biplane +72 -> +21 nm (150 cells, 20 min). The
     2026-10-14 (b) E-P modality numbers (contrast 2D 4.77 / z-stack 4.25 / biplane 4.45) were made with it; RE-RUN 2026-11-09: 4.78 / 4.34 / 4.39, see "Improvement pass".
  4. `biplane.estimate_plane_offset`: the "bead bright in both planes" filter compared the WHOLE plane's maximum (the same for every bead), so it never
     filtered; now each bead's own peak.
  5. `compare.compare` (method comparison, Spot finding page part 4) ran every variant with only `dilate_px` from the run's settings: on multi-copy runs
     (max_copies 2, pair prior gauss, max_candidates 8 since 2026-10-19) its "default" row was the one-copy finder. Now variants start from the run's own
     `locus_analysis`. Re-run on the two-copy showcase 2026-11-09.
  6. `ms2_trace` tracked the promoter with a hard-coded 2 px mask dilation; now `dilate_px` from locus_analysis (pipeline, error_effects, sweep script pass it).
     No change for the current configs (all use 2).
- **Checked, left as is:** `match_copies` has no distance cap (a junk path could be credited as a found copy): over all 507 matched copies of the saved
  multi-copy runs none is > 1 um from its true copy (worst 0.59 um), so no score is affected. The E-P z-stack measurement is ~2.5 s after the frame time
  (plane scan), so its error against the frame-time truth (314 nm rms) is ~2x the sigma the EIV kernel assumes (141 nm): physical, but the kernel does not
  know it; SETTLED 2026-11-09, no change needed (see "Improvement pass"). `_apply_overrides` failed on a config node that
  is null (`probes:` empty): fixed 2026-11-09.

## Improvement pass: stale results re-run, tutorial updater (2026-11-09/10; user: "make any improvements you see fit ... update the tutorial ... make it easy to update the tutorial in future passes ... make sure the user interface is up")
- **Everything that predated the 2026-10-26 logic fixes was re-run** (CPU only, ~3 h; the old files are kept next to the new ones in folders
  `before_2026-10-26_fixes`). Nothing changed beyond seed scatter, so every earlier conclusion stands:
  - E-P polymer scans (all 12, `scripts/sweep_ep_polymer.py`): contrast 2D 4.77 -> 4.78, z-stack 4.25 -> 4.34, biplane 4.45 -> 4.39 (sd 0.8-1.2; reference 4.19); table old -> new in
    docs/EP_CURVE.md "re-run 2026-11-09".
  - Sampler convergence (`scripts/check_mcmc.py`, with the fixed shift move): constructed traces gibbs_shift 5/5, gibbs_pt 5/5; sparse gene (seed 7) gibbs_shift 12/12 (was 9/12),
    gibbs_pt 12/12; dense gene (seed 8) gibbs_shift 3/12 (unchanged: the trapping is real, not the old bug), gibbs_pt 12/12 (was 11/12).
  - Deconvolution scores (auto): seed 8 found 0.830 -> 0.832, false per true 0.198 -> 0.194, rate r 0.798 both; seed 11 0.786 -> 0.784; seed 7 false per true 0.17 -> 0.10 (5 true starts);
    library seed 12 found 0.889 -> 0.891.
  - Method comparison on the two-copy showcase loci_library_seed12_87a16d now uses the run's own (multi-copy) settings: default finds 98.6% (calibrated 98.8%), 0.6% wrong,
    median 46 nm (calibrated 21 nm); the old brightest-pixel finder 58%, 20% wrong.
- **The z-stack timing question of the audit is settled, no change:** the EIV kernel is the error relative to the distance AT THE MEASUREMENT TIME (correct); the 2.5 s plane-scan
  delay only shifts the distance history against the promoter state, negligible next to the ~45 s the promoter integrates over.
- **Small fixes:** `localize.match_copies` no longer credits a found path as a true copy when its median distance to it exceeds `MATCH_MAX_UM` = 1.5 um (it becomes an extra
  copy and the true one counts as missed; no saved score changes, worst real match 0.59 um); `io.runs._apply_overrides` works on a block left empty in the YAML (`probes:` = None).
- **Tutorial updater (the way to update the tutorial from now on):** `scripts\tutorial\update.cmd` = status (stale screenshots + a check of all 250 sentences in the player);
  `update.cmd all` = re-take stale screenshots, record edited sentences (Kokoro), rewrite NUMBERS.md, check. Guide: `docs/tutorial/UPDATING.md`.
  - `scripts/tutorial/update_tutorial.py` (runs with the shot venv; starts its own dashboard on 8502 and stops only that PID tree), `shot_table.py` (one row per screenshot: capture
    keys, run, source files, sweeps), `docs/tutorial/shots_stamp.json` (what each screenshot was made from: source hashes + newest result file of its run / sweeps; a shot is STALE
    when one changed), `docs/tutorial/NUMBERS.md` (generated: every sentence that states a number, to re-read when results change), `tests/test_tutorial_tables.py`.
  - The check catches script errors and `no rect` / `missing shot` warnings (verified by breaking a rectangle name on purpose). It cannot tell whether a sentence is still TRUE.
  - All 28 screenshots were re-taken on 2026-11-10 from the re-run results; one sentence corrected (stray spots are "nearly as bright as a locus": 0.7 x). The narration quotes
    almost no result numbers, so nothing else needed changing. The MP4 release (2026-10-25) was NOT re-made (pictures differ only in a few numbers; publishing needs the user's OK).
  - Gotcha: writing regex source through a bash heredoc + Python string turned `\b` into backspace bytes; write such files with the Write/Edit tools.
- The tutorial showcase loci_demo_seed11_24bba2 now has its sampler check (gibbs_shift 12/12, gibbs_pt 12/12) and a method comparison (default 99.8% found, 0.07% wrong, 41 nm; calibrated 26 nm; brightest-pixel 9% wrong), so MS2 page 5b and Spot finding part 4 are no longer empty in the tutorial pictures.
- 225 CPU tests.
