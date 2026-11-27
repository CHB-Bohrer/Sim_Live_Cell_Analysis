# Chromatin library, polychrom engine, polymer calibration, stage 7b

Working notes moved verbatim from `CLAUDE.md` on 2026-11-27 (section order kept). Add new dated sections at the top; delete statements that are no longer true.

## Chromatin library (stage 1; decisions agreed with the user 2026-10-05)
- Chromatin is simulated SEPARATELY from the microscopy pipeline as a library of saved polymer simulations
  (polychrom/OpenMM on the GPU, Mirny-lab style). Stage 1 of a movie assigns one saved simulation to each cell
  (`library.assign_simulations`); microscopy labels/fluorophores and mapping into the nucleus come later.
- User's choices: first models = plain confined polymer + loop extrusion with boundary elements (compartments later);
  1 kb per monomer; 100 simulations per library; MSD calibration (nm, s) later (user will supply measured MSD).
  The user first said "chromosome 21 scale" (N = 46,710) and then dropped it ("lets not focus explicitly on chr21
  anymore"): the region size is now a plain setting `polymer.n_monomers`, default 10,000 (the validated 10 Mb).
  Boundary-element positions are SYNTHETIC (seeded random TAD sizes) until real CTCF positions are supplied.
- Code: `stage1_chromatin/lef1d.py` (numba 1D loop extrusion with boundaries; 1 monomer/step/leg, processivity = 2 x
  lifetime), `polymer3d.py` (3D polychrom runner, parameters from the polychrom `loopExtrusion` example), `analysis.py`
  (contact maps, P(s), insulation, corner dots, MSD), `library.py`. Scripts: `run_chromatin.py` (one sim),
  `run_chromatin_library.py` (resumable library runner), `validate_chromatin.py` (reproduce known behaviour ->
  `docs/chromatin_validation.md`). Configs in `configs/chromatin/`. Output: `$SIMLIVE_DATA/chromatin/<library>/<sim>/`.
- Reference numbers (Fudenberg et al. 2016 Cell Reports, fetched from the paper): LEF processivity ~120-240 kb, LEF
  separation ~120 kb, 10 Mb region, impermeable boundaries in the minimal model, contacts ~2-fold lower across TAD
  borders than within. polychrom example: 750 MD steps/block, bond wiggle 0.1, angle k 1.5, repulsion trunc 1.5 x 1.05,
  LEF bonds length 0.5 / wiggle 0.2. Do not quote paper numbers from memory; re-check them.
- Status (2026-10-08): engine built, unit-tested and RE-VALIDATED (`docs/chromatin_validation.md`; 5 conditions x 2
  seeds at 10 Mb; 11 of 12 checks pass, 1 honest CHECK). History: the FIRST validation (b813cf4) was INVALID because the
  3D loop-extruder bonds were frozen at their starting positions (OpenMM `updateParametersInContext` cannot change which
  particles a bond connects; loop legs were ~12 apart instead of ~0.5). Fixed in cc9dad0 (`LEFBondUpdater` pre-registers
  every bond and switches them by stiffness, as the polychrom example does), guarded by `tests/test_chromatin_gpu.py`
  and a leg-distance check in the report; invalid sims are in `$SIMLIVE_DATA/chromatin/validation/_superseded_static_bonds/`.
  Corrected results: leg-to-leg distance 0.61; plain-polymer P(s) slope -1.43; loop extrusion raises P(s) 2.8x
  (peak ~150 kb) and depletes it beyond ~1 Mb; mean loop 101 kb; designed 120-180 kb loops held ~10% of the time with
  anchors ~35x more often in contact than a typical pair (control 0.9); locus MSD exponent 0.49 (COM-subtracted).
  THE ONE CHECK: impermeable boundaries insulate 1.36-fold (seeds 1.42, 1.31; target >= 1.5 set in advance; controls
  1.05 and 0.82-1.11), i.e. a real but weaker effect than the ~2-fold in Fudenberg 2016, whose quantity (between vs
  within TADs) is not the same metric. Leaky (stall 0.99) 1.40 is indistinguishable from impermeable with 2 seeds. More
  seeds (e.g. 5 per condition) would tighten this; do not tune the threshold to pass. The old ~1.9 was a frozen-loop artifact.
- Lesson (bug): never re-point existing OpenMM bonds with updateParametersInContext; and sanity-check couplings
  physically (here: distance between the two legs of every loop must be ~bond length), not only by downstream summary
  statistics, which can look plausible for the wrong reason.
- Lessons for future analysis (do not repeat): a within/across insulation ratio is biased by distance decay (use the
  pooled across-OE fold); dots must be tested on the anchor monomers' own contact probability pooled over loops (10 kb
  bins dilute a single pair 100-fold; per-loop medians are dominated by unformed loops); boundary stall 0.9 per attempt
  gives NO insulation (use >= 0.99); validate metrics on null data (random positions, no-boundary control).
- Speed: N=10,000: ~3,650 MD steps/s (~7.5 min per 2,000-block validation run). N=46,710 (389 LEFs, 52 boundaries):
  ~1,400-1,750 steps/s, i.e. cost scales roughly with region size. Library configs
  (`configs/chromatin/library_*.yaml`: 1000 equilibration + 8000 production blocks x 750 steps since 2026-10-15, was 2000) at 10 Mb. MEASURED with 8000 blocks (loop extrusion): ~2,600 MD steps/s
  (3.5 blocks/s, GPU only ~60% busy) because every loop-extruder bond is pre-registered (about 3.5x more bonds than at 2000 blocks, and `updateParametersInContext` walks all of them each block):
  ~43 min per simulation, ~73 h per 100-simulation library, ~30 GB. The old 2000-block figure was ~11 min / ~18 h. Possible speed-up (not done): register only the bonds of the next few hundred
  blocks. THE LOOP-EXTRUSION LIBRARY IS FINISHED (100/100, 2026-10-18 19:10, ~43.7 min per simulation; see "Library run" below). Whether 1000 equilibration blocks
  suffice is only validated at 10 Mb (P(s) stationarity 0.077 over 2000 blocks); re-check for larger regions.
- Dashboard: the 🧬 Chromatin tab shows the validation report and any saved simulation, and has a smooth in-browser
  player of the polymer through time (`app/chromatin_view.py` + `app/chromatin_tab.py`: whole chain at each snapshot, or
  every block at 1-in-10 monomers; two chosen loci A/B with their distance; loop-extruder lines).
- If a GPU run fails with `CUDA_ERROR_INVALID_VALUE` the GPU/driver was reset (e.g. after a reboot or sleep): check
  `nvidia-smi` and `pytest -m gpu`, then simply re-run; `validate_chromatin.py run` and `run_chromatin_library.py` resume.
- polychrom works on this machine: CUDA platform, ~7,000 steps/s for N=2,000 monomers and ~3,650 steps/s for N=10,000
  (variable Langevin, spherical confinement). In THIS polychrom version the Simulation args are `error_tol`,
  `reporters=[...]` (not `error_tolerance`/`reporter`).

## Polymer calibration, stage 7b and code clean-up (2026-10-15; user said "do all of them" after a code review)
- **Polymer calibration is now ONE shared module** (`src/simlive/calibration/polymer.py`; it was inside `stage7_validation/ep_polymer.py`, which re-exports it): `PolymerPool`, `load_pool`, `calibrate`,
  `Calibration` (JSON save/load, now also stores working distance, contact fraction, separation and the simulations used), `scale_from_config`, `default_scale`. `scripts/calibrate_chromatin.py` writes
  `configs/calibration/default.json` (committed: 55.42 nm per unit, 0.461 s per block, MSD exponent 0.57 below one block). `configs/loci_library.yaml` now has `loci.calibration: default`;
  explicit `loci.nm_per_unit` / `loci.block_duration_s` (config or `--set`) still override that one number; no `calibration` key = the config's own numbers (old placeholders). The dashboard (Library and
  Transcription tabs) starts from the calibration's numbers; RESTART the dashboard (helper modules do not hot-reload).
- **Its inputs are ASSUMPTIONS, not measurements** (the user's 2026-10-14 rules: 100 kb pair within 400 nm 10% of the time; single-locus D = 0.002 um^2/s at 1 min). **CONSEQUENCE found while wiring it:** one saved block is 0.46 s, so a
  2000-production-block simulation (`configs/chromatin/library_*.yaml`) covers only ~15 min of movie, while the default movies are 20-40 min and a transcription study needs hours (bursts ~3/h). Longer movies repeat the
  trajectory back and forth (warning, correlations beyond 15 min are artificial) and loop-extrusion timing is fast (loop lifetime ~46 s). **DECIDED with the user (2026-10-15):** 8,000 production blocks (= 61 min of polymer
  time, a 1 h movie without repeating) at the validated 10 Mb region; every cell draws its own random 100 kb pair (`loci.positions_kb: random`, `pair_separation_kb`, `pair_margin_kb` 1000; `library_loci.draw_pair`; pairs recorded in
  `stage1_chromatin/library_pairs.csv`), so 100 simulations serve hundreds of cells. Movies longer than ~61 min still repeat the trajectory. Still open: whether a surrogate distance process (PLAN_transcription step 5) should take
  over for hours-long movies; loops live only ~46 s at this timescale (real cohesin: minutes), which a slower extrusion per block would fix at the price of even more blocks.
- **Library run (2026-10-15 18:24 to 2026-10-18 19:10, FINISHED 100/100, no failures; detached via PowerShell `Start-Process`, output `C:\Users\cbohr\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data\logs\library_loop_extrusion.out.log`):**
  `scripts\run.cmd python scripts\run_chromatin_library.py --config configs\chromatin\library_loop_extrusion.yaml --library loop_extrusion --n 100`, seeds 1000-1099, folder `$SIMLIVE_DATA/chromatin/loop_extrusion/library_loop_extrusion_seed<seed>`
  (the simulation NAME is `library_loop_extrusion`, so set `loci.library: loop_extrusion`). Watch it in the 🗂 Library tab or `status.json`; resumable: after a reboot / sleep run the same command again. WHILE IT RUNS do not edit
  `configs/chromatin/library_loop_extrusion.yaml`, `polymer3d.py`, `lef1d.py` or `io/runs.py` (each simulation is a new process that re-reads them) and do not run GPU tests or other GPU work. The plain-polymer library
  (`library_polymer.yaml`, also 8000 blocks) has NOT been started.
- **Stage 7b built** (`stage7_validation/propagation.py`, called from `pipelines/tracking_demo._propagate`; tests `tests/test_identity_propagation.py`): for tracked identity sources the 2D locus detections are mapped back to the
  TRUE cell through `identity.match_ids` (largest overlap per frame) and scored with the usual `score_localization` plus locus-distance error, share attributed to a true cell, and share on "mixed" tracks (a tracker label that
  stood for two true cells). In `scores.json` as `propagation.<source>.{uncalibrated,calibrated}`; per-detection errors in `stage6_analysis/<src>/localization_errors_via_tracks[_cal].csv`. Check: with perfect tracking it equals the
  truth-source scores exactly (smoke library, seed 5: rms 0.1187 um both). NOT done: a dashboard view of it, 3D, stitching broken/switched tracks into one cell (the rest of stage 5), a movie with real ID switches (unit tests only).
- **Refactor, verified identical:** `_isolate_and_analyze` (110 lines) is now small functions (`_locate_2d`, `_propagate`, `_trace_ms2`, `_locate_3d`, `_correct_with_fiducials`, one `_write_locus_results` for every positions / distances /
  errors CSV). Same config before and after (2D with probes + MS2, and 3D without beads) gave byte-identical files; only memory readings differ. **Bug found by it and fixed:** the fiducial-bead step needed
  `stage4_segtrack/auto_segmentation.tif`, so every 3D run with beads and `tracking.skip: true` crashed (no scores.json); it now falls back to the true labels.
- Small fixes: `simlive.stage1_chromatin.rotation.random_rotation` replaces two copies; Trackastra loads its model once per tracker (was twice per run); `provenance.json` now has `git_dirty` + `git_diff_sha256` (state at the first
  call in the process); `scripts/sweep.py` writes `failures.csv` and exits non-zero if any run failed (it used to print and carry on); the dead `_stored` is gone; `pyproject.toml` explains why `dependencies` is empty (environment.yml
  pins; do not let pip re-resolve); `docs/ci/syntax.yml` is a ready GitHub Actions workflow that only byte-compiles every file (the tests need the conda env + GPU, so they would NOT be in CI). It is NOT active: the push was refused because the GitHub token lacks the `workflow` scope; to enable it run `gh auth refresh -s workflow` and copy it to `.github/workflows/`. New tests: identity / propagation, calibration, segmentation, tracker registry, provenance.
- Gotchas: (1) `python` in the Bash tool is the Microsoft Store shim: run Python only through `scripts\run.cmd`. (2) Do not edit source while a run is going: the process pool re-imports modules from disk. (3) Python `write_text`
  on this Windows machine writes CRLF (the repo's index is LF, git normalises; fine) and a `"\r"` inside a normal string literal in an edit script silently becomes a carriage return (it garbled a README edit once; reverted before commit): use the
  Edit tool or raw strings for docs. (4) skimage warns that `remove_small_objects(min_size)` is deprecated (removal in skimage 2.0): `segment.threshold_watershed` will need `max_size`.
- From the review, fixed 2026-10-20: Cellpose + Trackastra really run in `tests/test_segtrack_gpu.py` (GPU); every dashboard page is opened by `tests/test_dashboard_pages.py`
  (AppTest, on the saved showcase runs; skipped where they do not exist); the broad `except Exception` blocks now name what they catch (or report why, e.g. the Calibration
  page); the Metrics page shows `propagation` (stage 7b). CLAUDE.md was condensed on 2026-11-27: it now holds only the every-session rules, and its sections were moved verbatim into `docs/notes/`.
