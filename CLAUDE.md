# Sim_Live_Cell_Analysis

Repo for validating live-cell microscopy analysis of chromatin loci against simulations with known ground truth.
Package name: `simlive` (src layout). Owner is a biophysicist (GitHub: CHB-Bohrer; repo
https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis, public). See README.md for how to run everything.

## Session handoff rule (always follow)
**Before ending a work session, or after any significant change, update this file and README.md without being asked:**
what changed, new commands/files, gotchas hit, findings, open questions, next steps. Then commit. A new session starts
from this file and the code alone, so anything not written here is lost. Keep the "Status" and "Findings" sections
below true; delete stale statements rather than appending contradictions. `docs/ARCHITECTURE.md` holds Mermaid diagrams
(pipeline, validation logic, roadmap) whose boxes are coloured built / stand-in / planned: recolour or edit them when a
stage's status changes (GitHub renders them; to check syntax, render the blocks with mermaid.js in a browser).

**Pushing to GitHub (public repo):** the user said (2026-10-09) to push STRAIGHT TO `main` every time: after committing, run
`git push origin HEAD:main` (fast-forward only; never force-push; never open a PR unless asked; the working branch may be pushed too).
Before every push: run `git status -sb`, check the staged diff for secrets (password/token/key strings) and large files, and confirm
`data/runs/` and `data/sweeps/` are still ignored. If a check fails, do not push; tell the user. (Earlier text said "push the working
branch, not main": that was a misreading of the user's intent.)

## The experiment being simulated (from the user; some details still open)
- Loci: promoter (locus 0) and enhancer (locus 1), 100 kb apart (user, 2026-10-09), plus MS2 bursting at the promoter. MS2 is a THIRD
  colour, but the SAME colour as the nuclear NLS-GFP (MS2-MCP has an NLS and accumulates in the nucleus), so MS2 spots are drawn into the
  nuclear channel; the user expects only a small rise in background (we test, not assume, the effect on segmentation).
- Imaging modality (widefield / spinning-disk / lattice light sheet), 2D vs 3D z-stacks, pixel size, NA, wavelengths,
  frame rate, exposure, movie length and cell count must all be easy-to-vary config parameters.
- The user's analysis will quantify each locus's position through time; not written yet (stage 6).
- MSD calibration targets: the user has given only an ASSUMPTION (single locus D = 0.002 um^2/s at a 1 min lag; 100 kb pair within 400 nm 10% of the time); these set
  `configs/calibration/default.json` (see "Polymer calibration" below). Still wanted: a measured MSD (D, alpha, length/time scales).
- Still needed from the user: real frame interval, typical cell speed, cell density, whether cells divide, and ideally
  a real movie with some tracked cells to check how well the simulation matches reality.

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

## Known problems fixed, tutorial link (2026-10-20; user: "Fix all known problems and make a link to the tutorial on the main page please")
- **Tutorial link on the dashboard home page** (both home versions): `style.tutorial_banner` (styled HTML, wraps on narrow windows) links to
  `app/static/tutorial/index.html` (opens in a new tab) and to the video release. `style.publish_tutorial()` copies `docs/tutorial` into
  `app/static/tutorial` (changed files only; gitignored): a page served over http may not link to `file://`, Streamlit refuses symlinks that
  leave `app/static`, and Streamlit 1.65 serves .html/.js with the right content types (checked). docs/tutorial stays the source.
- **Multi-copy linker:** `multicopy.untangle_by_place` + `repair` (docs/SPOT_FINDING.md "Copy swaps through wrong picks and gaps"): copy swaps made
  through wrong picks or gaps were invisible to the step costs (showcase, calibrated, tracked IDs: 7.9% on the other copy -> 0, rms 0.80 -> 0.15 um).
  25 sweep movies re-analysed, all three variants better (default: other copy at most 0.8%, was up to 2.5%; wrong picks 2.3 -> 0.6% at 1 stray/min),
  none worse by > 1 point. Dropped on the way (they made things worse): absolute territories in the nucleus frame (orientation flips by 180 deg on
  round nuclei), the pass on per-colour paths before pairing (stray paths at 4/min). 3D with two copies works (biplane run seed 72: 8/8 copies);
  3D scores now report stray / other-copy picks too.
- **E-P curve:** wrong-pick mixture (`with_outliers`, `assumed_outlier_rate`; corrected 4.11-4.13 vs reference 4.07 at 5-40% wrong picks when the
  share is supplied); bootstrap re-fits the prior per replicate; the "intervals 25% too narrow" claim WITHDRAWN (scan settings reused the same 20
  populations; 60 fresh ones give coverage 0.95). E-P page has the checkbox. docs/EP_CURVE.md "Update 2026-10-20".
- **Crash fixed:** cell placement drew nuclei one at a time and gave up when early ones blocked the rest (4 nuclei in the 40 um biplane field
  crashed; 8 never fit). `motion.place_cells` redraws the whole layout (100 tries; the first try uses the old random numbers, so every run that
  worked is unchanged, tested) and `placement_problem(cfg)` lets the sidebar refuse an impossible layout BEFORE starting, with the reason.
- **Transcription:** default coupling is now hill, half effect at 400 nm (transcription.DEFAULTS, loci_demo/loci_library/loci_3d configs, page):
  the user's E-P choice; the sharp 150 nm contact was reached in 0.75% of frames of the library distances (400 nm: 11.8%), so it barely coupled.
  BUG fixed: the sidebar hand-off fell back to the OLD TFF1-like gene when the Transcription tab had not been opened. Step 5 opens on the first burst.
- **Dashboard:** run summary + tracking scores only above pages that use the run (`NOT_ABOUT_RUN`: Scans, E-P curve, Library, Chromatin); one
  library selector on the Library page (and Chromatin opens on it); the live monitor checked during a real 2-simulation GPU library run, and a
  stalled generator no longer shows "now running" / time left; stage 7b table on the Metrics page (`propagation_view`); guide / glossary updated.
- **Code:** broad `except Exception` blocks narrowed to what they mean; `tests/test_dashboard_pages.py` (every page on the saved showcase runs + the
  empty home page; skipped where the runs do not exist; it caught a syntax slip in library_tab.py that this round introduced); `tests/test_segtrack_gpu.py`
  (Cellpose + Trackastra on a synthetic movie); new tests for placement, untangle_by_place, repair, wrong-pick correction.
- **Lesson (again):** do not edit `src/` while a pipeline run or a re-analysis is going: worker processes import modules from disk (a 3D run crashed
  with "name 'repair' is not defined" between two edits of one file). Do multi-step source edits in one write, or wait.
- **Needs the user (cannot be fixed here):** measured values for every placeholder (MS2 brightness R, instrument / camera, probe kinetics, stray-spot
  rate and brightness, wrong-pick share for the E-P fit), a measured MSD for the polymer calibration, a real movie; which meaning of nucleus
  "steadiness"; whether loops should live longer (more blocks per simulation) and whether a surrogate distance process should cover > 61 min movies;
  the plain-polymer library (73 h of GPU) needs a go-ahead; `gh auth refresh -s workflow` to enable docs/ci/syntax.yml; two stale entries
  `.git/worktrees/head_wt` and `.git/worktrees/wt06e` carry an "Everyone: deny delete" permission (made when worktrees were created inside the
  sandbox) and make every git command print "Permission denied": remove them as the user (e.g. `icacls <dir> /remove:d Everyone` then delete).
- **Planned work, not defects (unchanged):** stage 5 stitching of broken / switched tracks, 3D tracking and 3D cell isolation, Ultrack / TrackMate,
  density and frame-interval scans, choosing spots with promoter + enhancer + MS2 jointly, sister chromatids, E-P pool from the 100-simulation library.

## Side views along z, natural tutorial voice, bug sweep (2026-10-25; user: "visualization of the spot along the X-Z axis in addition to X-Y", "update the tutorial in full", "default speed 1.5x", "still sounds like a robot")
- **Side views (Cells page, `app/cells_tab.side_views`)**: for z-stack runs, below the x-y figure of "One cell through time": x-z and y-z cuts of the chosen spot's own colour stack
  (`<ch>_zstack.tif`, frame of the slider; 3 pixel rows/columns averaged, zoomed 3 um sideways and 4 um up/down, each plane minus its own median background so the spot is not lost in
  the fixed-pattern noise), with found (white +, from `loci_positions_3d[_cal].csv`; radio when both exist), true (cyan o) and drawn (red x) positions, plus the height of every spot
  through time with the true heights as thick faint lines. Works for true and tracked IDs, with copies; plane heights from `imager_for_run(...).zk` (biplane too). 2D runs show a one-line notice.
  Tested on the 3D, biplane and two-copy 3D runs (`tests/test_dashboard_pages.py::test_cells_side_views`). NOT done: x-z view in the Spot finding page's per-step images and in the player.
- **Bugs found and fixed:** Cells page crashed on runs made before the `visible` column of `localization_errors.csv` (loci_demo_seed1_10849d); the tutorial said 26 image-error sources (there are 28; README said 22);
  `build_data.py` dropped every shot not in the folder it was given (it now merges into `shots_data.js`, so single shots can be re-taken).
- **Tutorial voice:** the browser voices sound robotic, so the player plays recorded clips (`docs/tutorial/audio/*.mp3` + `manifest.js`). FIRST version: Microsoft edge-tts Jenny (`make_clips.py`), judged still
  too synthetic by the user; CURRENT: Kokoro-82M, voice `af_heart` (user's pick from 4 samples), run locally on the CPU with NO text leaving the PC: `make_clips.py --export beats.json` (shot venv) then
  `make_clips_kokoro.py beats.json` (kokoro venv `C:/Users/cbohr/kokoro_tmp/v`, model files in `kokoro_tmp/model`, setup + SHA-256s in the script header; downloaded 2026-10-25 with the user's OK, 350 MB). 250 sentences, 22 MB.
  `PRON` in the script respells names the model misreads (Viterbi, Cellpose ...; spoken only, captions unchanged): extend it if a word sounds wrong. A clip is used only if its hash matches the sentence in `engine.js`
  (`hashText`), otherwise the browser voice speaks it, so after editing narration re-run `make_clips.py` (only changed sentences are synthesized). `docs/tutorial/audio/` is gitignored (like the MP4, not in
  the history): a fresh clone has the browser voice until the clips are made or published. Voice list has "Natural voice (recorded)" first. Default speed 1.5x then (2x since 2026-10-26, see "Tutorial playback" below).
  The MP4 was re-made the same day with this voice, the side-view sentences and the 28-source fix: release `tutorial-video-2026-10-25` (49 min, Kokoro af_heart voice; 126 MB + 88 MB 720p; first uploaded with the Jenny voice, assets replaced the same day; published at the user's request; the README, the dashboard home link and docs/tutorial/README point to it; the 2026-10-20 release stays). Local files: `C:\Users\cbohr\SimLiveData\tutorial\video_2026-10-25\`. The video narrates at normal speed (1x); only the web player defaults to 1.5x.
- New shot `pg_cells_3d` (capture_shots.py `cells_3d`, run `loci_3d_seed23_2163f4`); scene "cells" has 4 new sentences (250 sentences in 38 scenes now).

## Tutorial playback (2026-10-26; user: "the movie consistently lags ... the animations tend to go and then stop", "default speed ... 2x")
- **Cause found:** the PC itself. Windows' compositor (`dwm.exe`) kept the GPU's 3D engine at 81-97% with nothing else running (one 3440x1440 60 Hz screen;
  uptime 18 days). A Chrome trace showed the page's own work at ~1% per frame, but `DXGISwapChainImageBacking::Present` holding single frames 0.3-0.8 s
  (~18 frames/s). A blank page with one moving circle stalled the same way, so the tutorial code was not the cause. Check with
  `Get-Counter "\GPU Engine(*)\Utilization Percentage"` (dwm should be a few %); a sign-out or reboot usually clears it. We do NOT change system settings.
- **What the page did wrong:** its animation clock cut every frame's time step to 50 ms, so with late frames the canvas animations ran at ~40% of real
  time (measured: 1.6 s of animation per 4 s), fell behind the voice and looked frozen. Now `engine.js` passes the real elapsed time (capped at 1 s) and
  `Canvas` sub-steps its physics (`step()`) in pieces of <= 50 ms (`STEP`) and paints once: 3.8 s per 4 s at the same frame rate. Draw functions only use dt for
  tweens / accumulators, which are fine with a large dt.
- The next sentence's clip is now preloaded while the current one plays (`preloadAfter`, `takeClip`; the dashboard's static server takes ~0.25 s per file);
  measured readyState 4 (fully loaded) at every play.
- **Default speed 2x** (options now up to 2.5x; storage key `rate20`, so a 1.5x saved under `rate15` does not override the new default).
- **Scene "Distance, promoter, Pol II, MS2" (`A.transcription`) fixed** (user: "the graphs don't show anything"): its seed search found no example
  and silently fell back to a flat distance with no bursts. Now 40 min, distance = smooth random walk (90 s, 750 +- 280 nm, ~10% in contact), the first
  seed with 3-5 bursts (>= 2 starting in contact, >= 1 far, every burst >= 2 Pol II); faint full traces always drawn, contact shaded in every panel,
  burst / Pol II counters, axes in nm / loops / minutes; the sweep restarts when a panel appears and lasts 26 s / speed. Narration unchanged (clips still match).
- Test the tutorial with Playwright from the shot venv (`C:\Users\cbohr\shotenv_tmp\v`) over http (`python -m http.server` in docs/tutorial) with
  `--autoplay-policy=no-user-gesture-required`; `headless=False` uses the real GPU. To measure the animation clock: `T.onFrame(d => sum += d)` vs wall time.

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

## Pinned versions (do not break)
- Exact working versions are saved in `environments/` (see `environments/PINNED.md`): env export, explicit conda list,
  pip freeze, polychrom commit `11a870c` (installer is pinned to it), model-weight SHA-256s, a frozen conda env copy
  `simlive-frozen-2026-10-05`, and git tag `env-2026-10-05`. NEVER run `pip install -U` / `conda update` in `simlive`;
  test upgrades in a separate env and save new lock files first.

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

## Environment
- Windows 11, NVIDIA RTX 3090 (driver supports CUDA 13.2). Conda env `simlive` from `environment.yml`, Python 3.12
  (Trackastra supports 3.10–3.13, Ultrack 3.11–3.13).
- OpenMM CUDA build (conda-forge), PyTorch CUDA wheels (pip), CuPy, polychrom (installed from GitHub, not PyPI), Trackastra.
- The one GPU is shared: OpenMM (stage 1) and PyTorch (stage 4). Run them sequentially, never concurrently.
- Verify setup with `pytest -m gpu`.
- Always run Python through `scripts\run.cmd` (activates the env). Calling the env's python.exe directly fails on
  DLL loading.
- Never import torch into the Streamlit dashboard (`app/dashboard.py`): on Windows it clashes with the plotting libs'
  OpenMP runtime and kills the server. The dashboard shells out to `scripts/run_tracking_demo.py` instead.
- Cellpose runs in a separate worker process (`stage4_segtrack/cellpose_worker.py`) with `skimage` import blocked:
  importing Cellpose after torch otherwise aborts with "OMP: Error #15" (two OpenMP runtimes). Do NOT use
  KMP_DUPLICATE_LIB_OK to work around it (can silently give wrong results). Same rule for any new torch-based tool
  that also pulls in scikit-image: isolate it in its own process.
- `scripts/sweep.py` runs grids x seeds and writes `data/sweeps/<name>/results.csv`. Single runs are not evidence;
  compare conditions over >=5 seeds (seed-to-seed sd of ID switches is large, ~5-10).
- polychrom is installed WITHOUT its Cython extension (`scripts/install_polychrom.ps1`); do not use knot-simplification
  functions.
- Dashboard: `scripts\dashboard.cmd` (port 8501, localhost only). The USER runs it themselves in their terminal; test
  on port 8502 (`.claude/launch.json` config `dashboard-test`) and never kill processes by name pattern (it killed the
  user's server once). Pages (navigation bar, only the open one runs): Movie, Metrics, Trajectories, Scans, Cells, Spot finding, Images, Calibration, Stage & fiducials,
  Transcription, MS2 → Pol II, E–P curve, Probes, Library, Chromatin, Config. DARK theme since 2026-10-12 (user asked for a darker UI): colours in
  `.streamlit/config.toml` + `app/style.PALETTE` (matplotlib figures dark too); hard-coded plot colours were remapped to light ones, so use
  light line colours (#e5e7eb, #2dd4bf, #60a5fa...) in new plots, never near-black.

## Status (update this section every session)
Built: stage 2 (cell motion, collisions, division, irregular time-varying nuclear shapes, `speed_scale`), stage 3
(nuclear channel + one colour channel per locus: per-nucleus texture + nucleoli and 22 switchable image-error sources, see "Image errors"), stage 4
(Cellpose default, threshold-watershed baseline, Trackastra `general_2d` behind the `Tracker` adapter), stage 5a
(per-cell isolation: `stage5_linking/isolate.py`, `celldata.py`), stage 7a (CTC metrics via traccuracy + own
identity-switch analysis), `scripts/run_tracking_demo.py`, `scripts/sweep.py`, dashboard (tabs incl. Scans, Cells).
STAND-INS (to be replaced): stage 1 = `stage1_chromatin/toy_loci.py` (2 loci per nucleus doing confined diffusion in
the nucleus frame, NOT polychrom); stage 6 = `stage6_analysis/spots.py` + `localize.py` / `localize3d.py` (placeholder for the user's
analysis, but now a measured one: see "Spot finding" below). The loci/isolation pipeline runs via `configs/loci_demo.yaml`.
Library wiring (2026-10-08): `stage1_chromatin/library_loci.py` (`loci.source: library`, config `configs/loci_library.yaml`; one simulation per cell, random rotation + time offset, polymer coordinates converted to true nm with `nm_per_unit` (from `loci.calibration`, see "Polymer calibration"; NOT rescaled to the nucleus) and placed at a random territory centre per cell (re-drawn until all loci stay inside the nucleus in every frame), `block_duration_s` likewise) is unit-tested with a fake library on CPU and was run end-to-end (tracked, 2026-10-15) on the 4-simulation `smoke` library; The 100-simulation loop-extrusion library exists since 2026-10-18 (showcase loci_library_seed12_87a16d uses it).
Dashboard redesign + library monitor (2026-10-08; the monitor was checked LIVE on 2026-10-20 with a 2-simulation test library on the GPU: running badge,
1 of 2 finished, current simulation's block and time left, GPU load, then "Stalled" after the generator was stopped): `.streamlit/config.toml` + `app/style.py` (one theme, CSS, matplotlib defaults), new 🗂 Library tab
(`app/library_tab.py`): live monitor (refreshes every 5 s) of `<library>/status.json` (written by `run_chromatin_library.py`) and the
running simulation's `progress.json` (written by `polymer3d.py` every 20 blocks): state badge (running / stalled / finished /
stopped), n done, ETA, GPU load via nvidia-smi, recent log; plus a picker to choose which simulations a movie uses. ONE library selector for the whole
page (key `lib_sel`, opens on a running library, else the one with most simulations; the Chromatin page opens on the same library) since 2026-10-20.
Choosing simulations: `library.assign_simulations(..., pool=, pinned=)`; config `loci.sim_ids` (restrict the pool) and
`loci.pin: {cell_id: sim_id}`; the sidebar's loci form has a source radio (stand-in / library) that uses the Library tab's choice.
The picker's preview uses the movie's own random stream (`SeedSequence(seed).spawn(4)[2]` = `spawn(10)[2]`: a child's numbers do not depend on how many
siblings are spawned), so it matches the real assignment when no cell divides.
Probe binding model (2026-10-08): `stage1_chromatin/probes.py`, stage 1b. Each locus (locus 0 = promoter, 1 = enhancer, separate parameters)
has N probes cycling FREE -(k_bind)-> BOUND -(k_scan)-> SCANNING -(k_off)-> FREE, plus BOUND -(k_unbind)-> FREE; optional `capacity` (binding
sites). Exact Gillespie (numba); exact theory in `summary()` (mean, CV, dwell time, correlation time); tested against each other
(`tests/test_probes.py`, 8 tests). `loci.probes.enabled: true` (configs loci_demo/loci_library, default false) makes locus brightness =
`photons_per_probe_s` x exposure x probes attached (exposure-averaged) in `render_loci(..., occupancy=)`; counts saved to
`stage1_chromatin/probe_occupancy.csv`; own seed stream ss[4]. All rates/brightness are PLACEHOLDERS (user will supply kinetics).
Not modelled: individual-probe bleaching, probe positions, signal spread along the locus while scanning, probe pool shared between loci.
Sweeps: `scripts/sweep_probes.py` (model only, seconds, -> `$SIMLIVE_DATA/probe_sweeps/<name>.csv`) and the 🧪 Probes tab; through the whole movie pipeline
with `scripts/sweep.py --config configs/loci_demo.yaml --set loci.probes.enabled=true --grid loci.probes.promoter.k_off=...` (its metrics are still
the tracking metrics only; locus-analysis metrics wait for stage 6/7b). Dashboard 🧪 Probes tab (`app/probe_tab.py`): explainer, per-locus
parameter panel with presets, live drawing of the model with probe counts, example trace + count histogram, parameter sweep, command/sidebar
hand-off to a movie (sidebar checkbox in the loci form).
Transcription / MS2 model (2026-10-09; stage 1c): `stage1_chromatin/transcription.py` (details and the user's decisions: `docs/PLAN_transcription.md`).
Promoter OFF -(k_on)-> ON -(k_off)-> OFF, Pol II initiation k_init while ON, Pol II ramp over the MS2 cassette then plateau until released
(random dwell); MS2 signal = loops carried by all Pol II. The promoter-enhancer distance d(t) (3D, from the library trajectories at full block
resolution, via `simulate_library_loci(..., return_traces=True)`; toy loci use their 2D lab distance) multiplies ONE chosen rate
(`coupled_rate`: k_on default | k_init | k_off) by 1 + (fold-1) f(d), f = hill (default since 2026-10-20: half effect at 400 nm, n = 4, the user's E-P choice;
the old default, a sharp 150 nm contact, was reached in only 0.75% of frames by the library distances, so it barely coupled) | contact | exp | none (control). ALL rates are per second
of real time, so the frame interval can change freely (library loci are interpolated between blocks; `check_sampling` warns when too coarse).
Exact Gillespie + exact expectations (`expected_on`, `expected_ms2`); 12 tests in `tests/test_transcription.py` (constant-distance theory, two-valued
distance vs exact expectation, frame-interval independence, null/positive correlation with contact, interpolation, MS2 spot in the nuclear
channel). `loci.transcription.enabled: true` (configs loci_demo/loci_library) writes `stage1_chromatin/transcription.csv` + `transcription_events.csv`
and renders the spot at the promoter into the nuclear channel (`render_nuclei(..., spots=)`; pipeline now runs stage 1 BEFORE rendering nuclei; seed
stream ss[5]). If the movie is longer than the saved polymer trajectory, the trajectory is repeated back and forth with a warning (correlations beyond
its length are artificial); the planned surrogate distance process (PLAN step 5) is NOT built. Dashboard ⚡ Transcription tab (`app/transcription_tab.py`): nine
numbered steps (sampling, distance, coupling, promoter, Pol II kymograph, MS2 signal with exact expectation, simulated microscope tiles, parameter sweep,
hand-off to a movie + sidebar checkbox). Its widgets start from `transcription.DEFAULTS` (until 2026-10-20 the hand-off fell back to the OLD TFF1-like gene
when the tab had not been opened: a sidebar movie then used a different gene than the tab showed); step 5 opens on the cell's first Pol II. LITERATURE NUMBERS ARE UNVERIFIED search-summary values (primary papers blocked from the build machine); all
MS2 brightness, coupling strength, initiation rate and nm_per_unit / block_duration_s are placeholders.
NOT built: MSD calibration from MEASURED data (an assumed-D calibration exists, see below), surrogate distance process for long movies, 3D tracking / 3D cell isolation, stage 5 beyond
isolation (stitching the locus results of broken / switched tracks into one cell), the user's stage 6 analysis, 7b for 3D (7b exists for the 2D locus results), Ultrack/TrackMate adapters,
density and frame-interval scans.

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

## Findings so far (cell tracking only; simulated 2D nuclei, 5 seeds per setting; measured BEFORE the image-error model, which changed little for tracking, see "Image errors")
- Errors come from crowding x speed: touching, featureless-ish nuclei; Trackastra (greedy, divisions allowed) often
  calls a touching pair a false division. Texture inside nuclei did NOT improve linking (CHOTA 0.913 vs 0.914 with
  perfect masks) but it breaks the watershed segmenter (0.89 -> 0.78); Cellpose is unaffected.
- Speed scan (speed_scale, default config): CHOTA >= 0.99 up to 0.7x, 0.975 at 0.8x, 0.91 at 1x, 0.85 at 1.25x. At
  <=0.5x linking with perfect masks was perfect (0 ID switches in 15 runs); ~1 switch/run remains from segmentation.
- `ilp` vs `greedy` mode: only helps at >=1x speed (fewer ID switches in 4/5 seeds); no consistent difference slower.
- Trackastra's pretrained model is meant to generalize, but it has NOT been checked on real data; a real annotated
  movie is the proper test.
- Default config (speed_scale 1.0) is deliberately hard; consider setting a slower default once the user's real
  cell speed / frame interval is known.

## Pipeline stages (each is its own subpackage under `src/simlive/`)
1. `stage1_chromatin` — polychrom/OpenMM polymer dynamics -> locus-locus distances vs time (polymer units).
2. `stage2_cells` — multiple cells with persistent IDs, diffusive/directed motion, collisions, optional division.
   Motion is applied *afterward* to saved chromatin coordinates.
3. `stage3_microscopy` — forward model: PSF, exposure blur, camera noise, photon counts, photobleaching, frame rate,
   localization error. Renders locus channel + nuclear/cell channel + per-cell ground-truth label masks.
4. `stage4_segtrack` — segmentation (Cellpose or StarDist) and tracking. Tracker sits behind a swappable adapter
   (`Tracker` protocol); Trackastra first, Ultrack and TrackMate later, all on identical synthetic data.
5. `stage5_linking` — link locus tracks to tracked cells.
6. `stage6_analysis` — the user's analysis pipeline (TBD).
7. `stage7_validation` — ground-truth comparison: tracking metrics (Cell Tracking Challenge: DET, TRA, LNK, etc.)
   and propagation of tracking errors into locus-distance results.
- `calibration` — converts polymer units to nm and seconds by matching locus MSD to measured values.
- `io` — shared readers/writers and schema definitions for stage interchange files.

## Tracking evaluation protocol
Run the tracker twice: (a) on ground-truth masks (isolates linking error), (b) on automatic-segmentation masks
(realistic). Score both against true cell IDs.

## Conventions
- **File contracts:** each stage reads and writes defined files only (no in-memory hand-offs between stages).
  Every stage output directory contains: data files, `params.yaml` (resolved config), `seed.txt`, `provenance.json`
  (git commit, package versions, timestamp), and ground truth where applicable.
- **Run layout:** `$SIMLIVE_DATA/runs/<run_id>/stageN_<name>/...` where `SIMLIVE_DATA` defaults to
  `C:\Users\cbohr\SimLiveData` (set in `scripts\run.cmd`; code uses `simlive.io.runs.DATA_ROOT`). Data was moved out of
  the repo/OneDrive on 2026-10-05. `run_id` encodes the config hash + seed. Outputs are gitignored and must be
  reproducible from config + seed. Chromatin simulations go in `$SIMLIVE_DATA/chromatin/`.
- **Configs:** YAML in `configs/`, one per stage plus sweep files that expand to a list of resolved configs.
  Never hardcode parameters in code.
- **Seeds:** all randomness flows from an explicit seed via `numpy.random.default_rng` (and `torch` / OpenMM seeds
  set from it). Derive per-stage/per-cell seeds with `SeedSequence.spawn`, never reuse global state.
- **Units:** stage 1 in polymer units; calibrated physical units (nm, s) from stage 2 onward. Name columns/vars with
  units (`x_nm`, `t_s`). Calibration parameters are recorded in the run.
- **Formats:** trajectories and tables as HDF5 or Parquet/CSV with a documented schema; images as TIFF/zarr; label
  masks as integer arrays whose value is the persistent cell ID.
- **Tests:** pytest, in `tests/`. GPU tests marked `@pytest.mark.gpu`. They must run in a SEPARATE process from the
  CPU tests (torch vs numba/skimage DLL/OpenMP clash on Windows; neither import order works): plain `pytest` = CPU tests
  only (`addopts = -m 'not gpu'`), `pytest -m gpu` = GPU tests, `scripts\run_tests.cmd` runs both. 218 CPU + 7 GPU tests at 2026-10-26. Prefer small deterministic tests with a
  known analytic answer (e.g. free diffusion MSD).
- Python 3.12, type hints on public functions, no wildcard imports.

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

## E-P distance -> P(burst ON) curve under measurement errors (stage 7b, built 2026-10-14; docs/EP_CURVE.md)
- User's decisions (2026-10-14): readout = P(burst ON); fit the WHOLE curve (not one threshold) with much more data; enhancer with a longer working distance (hill, d_half 400 nm, n = 4,
  fold 12 on k_on; the old sharp contact was 150 nm); contact fraction 0.1 (exact: the stand-in distance process is solved so that P(d < 400 nm) = 0.100).
- Code: `stage7_validation/ep_curve.py` (`EPSetup`, `simulate_population`, `Measurement` / `measure_distance`, `kernel_matrix`, `estimate_on` (MS2 trace -> deconvolution -> burst rule),
  `fit_curve` / `analyse` / `reference_curve`), `scripts/sweep_ep_curve.py` (experiments ladder, sigma, shift, outliers, cells, null, dynamics, contact, frames, ms2, rules ->
  `$SIMLIVE_DATA/sweeps/ep_curve_<exp>/results.csv`, ~60 min for all), page 🔗 E–P curve (`app/ep_tab.py`), `tests/test_ep_curve.py` (6 tests).
- Method: errors-in-variables fit of p_far (1 + (F-1) / (1 + (d/d_half)^4)) integrating out the true distance (measurement kernel + Richardson-Lucy prior), cell bootstrap intervals,
  effect size = contrast P(ON at 150 nm)/P(ON at 1000 nm) (F alone is not identified at the null); reference = same fit to true distance + true ON of 4000 cells (contrast 4.07; 6.7 if the
  distance were fixed: the promoter's response delay flattens it).
- Findings (placeholders: invented distance dynamics, 57 nm per axis, 2D, 5% missed, 2% wrong picks): null false positive 3.3% (60 populations); 2D projection of the distance costs
  ~35% of the effect if ignored (2.7 vs 4.1) and is fully corrected by the EIV fit; localization precision matters only above ~200 nm; colour shift up to 100 nm harmless; missed frames
  harmless to 60%; WRONG PICKS are the bad one (5% -> -6%, 10% -> -13%, 20% -> -24% uncorrected; corrected since 2026-10-20 when their share is supplied, see below); cells: contrast sd 0.77 (50) -> 0.055 (3200), d_half scatter 11% -> 2%
  with a +3.6% offset; contact fraction 0.03 is biased low (2.7 vs 3.5). MS2 -> ON: the balanced-accuracy burst rule has 8% false positives next to real bursts, costs ~10% of the
  contrast, and the simple (sensitivity, false-positive) attenuation correction FAILS for it (9.2 vs 4.07); conservative rules (2,2), (3,3), (8,6) are nearly unbiased uncorrected
  and the correction is then acceptable (4.3-4.9). MS2 noise of 2 x one Pol II with 20 s frames halves the contrast. (An earlier note "intervals ~25% too narrow" was WITHDRAWN 2026-10-20: see below.)
- 2026-10-20 (docs/EP_CURVE.md "Update 2026-10-20"): (1) wrong-pick mixture `ep_curve.with_outliers` / `Measurement.assumed_outlier_rate` (kernel (1 - w) K + w q, q = a wrong
  pick's distance distribution; w must be SUPPLIED, it is not identifiable from distances): contrast 4.11-4.13 at 5-40% wrong picks with the true w (reference 4.07; uncorrected
  3.83 / 3.55 / 3.09 / 2.35); half the share corrects half, twice the share over-corrects and blows up at 40%. E-P page: checkbox "Analysis knows the wrong-pick share".
  (2) bootstrap replicates re-estimate the true-distance prior and refit on the full grid; the "25% too narrow" came from every scan setting reusing the SAME 20 populations:
  with 60 fresh populations coverage of the contrast is 0.95 (null: 6.7% exclude 1, nominal 5%); d_half intervals cover 88% (its +2-4% offset).
- NOT done: other coupling shapes; the real polymer-library distances (needs MSD calibration); linking sigma to photons / the movie's spot finder; a convergence
  check inside `estimate_on` (single chain); the real MS2 noise relative to one Pol II (user must supply).

### E-P curve with POLYMER distances and three microscopes (2026-10-14, second version; docs/EP_CURVE.md "Update 2026-10-14 (b)")
- User's decisions: 3D z-stacks AND biplane, compared (2D kept as a reference); D = 0.002 um^2/s is for a SINGLE LOCUS at a 1 MIN lag (assumed corrected for nucleus motion and drift; if it is for the E-P distance the block is 1.94 s instead of 0.46 s);
  cell type difference (paper = mouse ES cells) ignored for now. The earlier 2D baseline was MY choice (the repo's default configs and the Vutara preset are 2D "bypass"; z-stack config exists; biplane was "not built yet" and is now built at vector level).
- Paper checked (Jusuf et al. bioRxiv 2025, v1 main text only; Methods/supplement and the journal version not readable): mESC, Micro-C calibrated on BILD live imaging of 2 loops, "looping probability" = fraction of time looped ABOVE the random-polymer background; NO distances,
  no nm threshold, no diffusion numbers in what I could read. So nm per unit = 55.4 comes from the user's contact-fraction rule (P(d < 400 nm) = 0.100 for 100 kb pairs, pooled over the 10 validation sims), and block = 0.461 s from the MSD match (exponent 0.50 above one block, 0.57 below).
- Code: `stage7_validation/ep_polymer.py` (pool of 10 sims x 990 pairs, `calibrate`, `simulate_population`), `ep_acquisition.py` (`Acquisition` kinds 2d / zstack / biplane, `measure`, `kernel_aniso`: errors from Thompson / axial Fisher information, exposure blur, the 0.15 s colour delay, piezo,
  all with the extrapolated sub-block MSD), `ep_history.py` (P(ON) from the distance HISTORY), `scripts/sweep_ep_polymer.py` (experiments modality, null, wrong, cells, contact, frames, colour, photons, focus, timescale, compare_sources, ms2 -> `$SIMLIVE_DATA/sweeps/ep_polymer_<exp>`),
  page 🔗 E–P curve part 1 (live polymer + microscope), `tests/test_ep_polymer.py` (5 tests; skipped without the validation sims).
- FINDINGS: polymer E-P distance decorrelates in 9.2 s but the promoter integrates over ~45 s => P(ON | distance NOW) is nearly flat (contrast 1.36 vs 4.1 for the invented slow process); the curve to fit is P(ON) for a distance HELD, from the distance history (reference 4.19, tau 45 s found). Error budget per axis (nm): localization 6-9 in
  plane / 27 (z-stack) or 51 (biplane) in z, exposure motion 38-62, sequential colours 0.15 s apart = 120 (the biggest), piezo 12; total ~126-140. At 600 cells x 1 h: contrast 4.77 +- 1.00 (2D), 4.25 +- 0.72 (z-stack), 4.45 +- 0.80 (biplane), truth 3.80 +- 0.54; 95% interval width 5.1 / 3.9 / 3.1 / 2.3;
  null false positives 3-7%. Need ~10x more cells than with the invented process; frames every <= 10 s (reference falls 4.19 -> 1.78 at 80 s); biplane loses loci when focus is off by > 0.6 um or photons < 0.3x; simultaneous colours cut ~20% of the scatter; the ON state from the MS2 trace is biased high (calls cluster near the enhancer) and unstable.
- NOT done: wrong-pick filter from the locus's past positions (the spot finder already links with a diffusion model; the scan's wrong picks are random draws); nucleus "steadiness" (simulated nuclei: outline changes 0.29 um rms per 10 s frame, long/short radius 1.6, rotation 4.5 deg and centroid step 0.9 um per frame; awaiting which meaning the user intends);
  rendering real z-stack / biplane movies and localizing in them (the acquisition model is vector-level); loop-extrusion timing is too fast at 0.46 s per block (loop lifetime 46 s); the E-P polymer pool still uses the 10 validation simulations (the 100-simulation library exists since 2026-10-18 and could replace it).

## Dashboard walk-through findings (2026-10-09)
- Tested: all 24 coupling shape x rate combinations, frame intervals 5/20/120/600 s (warnings sensible), both sweeps, probe presets,
  Library/Chromatin tabs (Streamlit `AppTest` scripts + the browser), and 4 sidebar movies (toy, library+probes, library+MS2, toy+MS2) all
  finish with DET/CHOTA 1.0 (8 cells, 1 seed; the default MS2 spots, up to ~250 loops x 100 photons, did NOT disturb segmentation; more seeds
  and brighter spots not tested).
- BUG FIXED (important gotcha): every library-source movie died at the tracker with `OSError WinError 127 ... torch\lib\shm.dll`. Cause: the
  FIRST numpy BLAS call (any `@`/dot/`scipy Rotation.apply`) in the main process loads MKL's OpenMP runtimes (libiomp5md + libomp), after which
  `import torch` (Trackastra) fails. Code that runs before the tracker must avoid matmul; `library_loci.py` now uses a broadcast sum and its
  own `_random_rotation`. Diagnose with `psutil.Process().memory_maps()` filtered on "omp". Same class of problem as the Cellpose worker note.
- Fixed: sidebar probe/transcription checkboxes did nothing (never passed on); movie locus positions came from the config's fixed 5050/4950 kb
  (both loci on the last monomer in a 3 Mb library) -> now from the Transcription tab or the middle of the library region; Transcription tab default
  positions follow region size; run log decode crashed on non-ASCII (cp1252) -> UTF-8; figure fonts were tiny.
- `st.rerun(scope="fragment")` buttons work in the browser but raise under AppTest (test artefact only).
- All four "not done" items of this walk-through were fixed on 2026-10-20: the live monitor was checked during a real (test) library run; one library selector per
  page; the run's summary and tracking scores are no longer shown above pages that do not use the run (Scans, E-P curve, Library, Chromatin; `NOT_ABOUT_RUN`);
  Transcription step 5 opens on the cell's first Pol II; the coupling default is the Hill curve at 400 nm (the 150 nm contact was reached in 0.75% of frames).
- Literature checked against PMC full texts, see docs/PLAN_transcription.md: TFF1 ON 16 min (not 6), burst period 66-86 min, 1.5 transcripts/burst;
  Pol II 2.6 kb/min, 3' dwell 116 s (our default 30 s). Defaults not changed: needs the user's choice of gene.

- Scans page (fixed 2026-10-14): it assumed the tracking-sweep layout (`masks` column + CTC metrics) and crashed on opening because the newest sweep
  (`ms2_brightness`, from `scripts/sweep_ms2_brightness.py`) has another layout. Now generic: sweep CSVs are read as "varied settings, seed, [masks], run, results"
  (columns before `run` = settings, after = metrics, `seconds` ignored); tracking sweeps keep their named metrics and the masks switch. Not caused by the UI speed-up.
- Players (app/player.py): every in-browser movie now shows 'Loading frames… n / N' with a progress bar and a disabled button, then 'Ready: press Play'. Chromatin player: the grey circle is the simulation's confining sphere (the territory, NOT the nucleus); caption and in-image label say so. The territory is ~1-1.5 um across at 50 nm/unit vs a 12 um nucleus.

- Display scaling (2026-10-11): single-cell and all-cells players scale each locus colour once for the whole movie from the cell's own mask region (fastview.movie_limits); before, each frame was scaled to the crop maximum, so cosmic rays and fiducial beads in the crop corners made the real spots look missing (spot < 30% of crop max in 70-92% of frames in 3D runs). Measured locus D = MSD/(4t) at 10 s: library loci 0.0017-0.0019 um2/s (falls to 0.0011 at 100 s, nm_per_unit 50 and block_duration_s 10 are PLACEHOLDERS), toy loci 0.0029; the nucleus itself moves 0.43 um per frame.

## Vutara biplane mode (2026-10-17; user: "change it, and why not build in the biplane ... you still need to do a zstack ... calibrate with the beads within the image")
- The sidebar label "Bruker Vutara VXL (2D)" was misleading (the instrument is 3D-capable; the preset modelled only its widefield mode): renamed "(widefield mode)"; new third choice "(biplane 3D mode)" = `configs/instruments/vutara_vxl_biplane.yaml`.
- Facts (Bruker bi-plane page, fetched 2026-10-17): two planes 600 nm apart, 50:50 splitter, two halves of one camera, z from the intensity difference by a 3D-PSF fit, PSF calibrated with a bead imaged as a function of z, z-stacks possible; NO axial range per acquisition and
  no registration step are stated. So a z-stack is still needed for a ~10 um nucleus (the user was right).
- Code: `acquisition.z_stack.biplane: {enabled, separation_um: 0.6, photon_fraction: 0.5}` in `stage3_microscopy/render.py` (`Imager.biplane`; `n_slices` = stage positions, planes = 2x, stage step = 2 x separation so the planes are EVEN 0.6 um and the 3D localizer is unchanged; both planes of a pair share the
  time index and the piezo error; each plane gets `plane_fraction` of the signal, background not halved); new error source `biplane_registration` (second plane shifted `shift_px`, its height `separation_error_um` off; in `_LATER`, off unless a config enables it; listed in `error_effects.ONLY_3D`).
  `stage6_analysis/biplane.py` `estimate_plane_offset` (beads between the lowest pair -> offset) + `localize_loci_3d(plane_shift_px=)`; pipeline `_undo_biplane_registration` writes `truth_3d_registered` and `biplane_registration` in scores.json. `tests/test_biplane.py` (5 tests). 176+ CPU tests.
- Measured once (seed 61, 2 cells, 12 frames): offset (0.321, -0.214) px vs true (0.35, -0.25); locus error unchanged within noise (dominated by the placeholder axial chromatic shift; calibration was off); biplane vs an ordinary 22-plane z-stack: same 3D rms (0.27 um), lateral rms 0.06 vs 0.10 um. ONE seed: not evidence.
- NOT done: true plane-separation and PSF z-calibration from a bead z-scan; registration with rotation / scale; a bead SLIDE calibration of the offset; biplane in the E-P tab's vector model is separate (`ep_acquisition.py`); the Images/Stage tabs have not been looked at with a biplane run.

## Narrated tutorial (2026-10-17; user: "in depth tutorial with animation and voice over ... the person should know how to use the UI")
- `docs/tutorial/index.html` (+ `engine.js`, `anims.js`, `scenes_a.js`, `scenes_b.js`, `tutorial.css`, `shots_data.js`, `img/`): self-contained, no libraries, works from `file://`.
  38 scenes / 246 sentences (~50 min; 2026-10-20). The dashboard home page links to it (served from a copy in app/static/tutorial, refreshed from docs/tutorial on every home-page load: edit docs/tutorial only). Narration = the browser's speech synthesis (Edge "Online (Natural)" voices are best; the picker is on the start screen). It is NOT
  Claude's own voice (not callable from a local page); the MP4 (below) uses a Microsoft neural voice instead.
- Screenshots are real dashboard pages: `scripts/tutorial/capture_shots.py` (Playwright driving the installed Edge, in a separate venv `C:\Users\cbohr\shotenv_tmp\v`, NOT the simlive env) captures TALL page images + the
  rectangle of each control; `scripts/tutorial/build_data.py <shots dir>` copies them to `docs/tutorial/img` and writes `shots_data.js`. Scenes refer to controls by name (`f("next", "Next error")`), so a
  layout change needs only a re-capture. Runs shown: `loci_demo_seed11_24bba2` (most pages), `tracking_demo_seed3_e15523` (tracking errors), `loci_demo_seed8_f5c88d` (Images part 2 measured effects), `loci_3d_seed23_2163f4` (Stage page),
  `loci_library_seed12_87a16d` (two copies + stray spots: shots `pg_cells_copies`, `pg_spots_copies`).
  The capture needs the test dashboard on port 8502 and takes ~10 min; it does not use the GPU. `python capture_shots.py <out> <shot> ...` re-takes single shots.
  Shots of 2026-10-20 in `C:\Users\cbohr\SimLiveData\tutorial\shots_2026-10-20b`. Hand-set rectangles live in `build_data.py` `OVERRIDES` (used where the label text
  also appears higher up, e.g. a dataframe column called "simulation"); check them after a re-capture by drawing them onto the shot (they went stale once when a table was added).
  Check every sentence for script errors / missing rectangles with a Playwright loop over `T.scenes` (each beat: open `index.html#s=<id>&b=<n>`, collect console errors and
  the engine's "no rect" / "missing shot" warnings).
- Gotchas hit: (1) the in-app browser pane gives only 800x500 screenshots and shows `file://` pages as static snapshots (no scripts), so test the tutorial with Playwright (`tut.py`-style: open `index.html#s=<id>&b=<n>`
  and screenshot); (2) Windows long paths: the venv must be at a short path; (3) in Streamlit, `get_by_text` finds the intro boxes before the real widget (e.g. "Next error" appears in the intro text), so rectangles are
  found as the smallest visible container whose text contains the label; (4) a tall (9000 px) viewport stretches the fixed page background, so the capture injects `background-image:none`; (5) text in page JS strings
  written by a Python patch script loses `\n` escapes: use the Edit tool; (6) 2026-10-20: switching runs raced with the page capture (the old run's page was saved; now `select_run` waits for the new run),
  pages taller than the 9000 px viewport were cut (copies shots use 15000 px; `save_page` reads the real viewport height), and clicking an already-open page in the segmented control DESELECTS it
  (`goto_page` now skips when `aria-checked == "true"`); (7) Streamlit keeps the PREVIOUS page's elements, dimmed, until the new page replaces them, so a long page
  saved too early ends in the old page's leftovers (the first Spot finding shot ended in Trajectories plots, parts 3d / 4 were missing): `save_page` now waits until no
  `[data-stale="true"]` element is left. A scene with several pictures (`shot({ sets })`) stays on the last one chosen: switch back explicitly with `fs("main", ...)`.
- MP4s are GitHub releases, NEVER in the git tree (user 2026-10-17: "don't add another 59 MB to the history"; the first 59 MB version is still in history). CURRENT: release `tutorial-video-2026-10-20`
  (https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-10-20 (SUPERSEDED by tutorial-video-2026-10-25, see the 2026-10-25 section); 52 min, with the copies scene; assets `sim_live_cell_analysis_tutorial.mp4` 110 MB and `_720p.mp4` 82 MB;
  published 2026-10-20 at the user's request; README links point to it). Previous: `tutorial-video-2026-10-17` (48 min, kept). Local copies in `C:\Users\cbohr\SimLiveData\tutorial\video_2026-10-20\`
  (+ `release/` with the asset names) and older ones in `C:\Users\cbohr\SimLiveData\tutorial\`. Publishing a new release is outward-facing: ask the user first. `gh release create --target` needs a branch name or full SHA.
  Delete the audio cache before re-running `tts.py` (clips are keyed by scene:sentence, not by text). Found while checking the 2026-10-20 video: the recorder's frames lag the page's sentence log by 0.3-2.2 s, varying, so the voice LED the picture by ~0.7 s (also in the 2026-10-17
  video). `mix.py` now finds each caption change in the recording and places the audio there (voice now 0.17-0.63 s after its caption, median 0.22 s, 49 sentences checked).
  The player also pre-decodes the current and next scene's screenshots (`engine.js` `decodedImg` / `keepDecoded`) so switching to a 9000 px picture does not stall it.
  Voice = Microsoft `en-US-JennyNeural` via `edge-tts` (the user chose it; Claude's own voice cannot be produced). Scripts `scripts/tutorial/video/{tts,record,mix}.py`; the player accepts `#s=<id>&b=<n>&play=1`
  and `window.BEAT_DUR` (fixed sentence lengths) for recording. Regenerating takes ~50 min of real time. Playwright's recorder needed an ffmpeg copied to `%LOCALAPPDATA%\ms-playwright\ffmpeg-1013\ffmpeg-win64.exe`.
- When the dashboard changes: re-run the capture + build, then fix any scene whose control was renamed. Sentences must stay under ~250 characters (some online voices cut off long utterances). NOT checked: how the
  voices actually sound (no audio output was available to the author), and the tutorial on Chrome/Safari.

## Dashboard background (2026-10-16)
- The page background is the Bohrer Lab website's faint chromatin field (user's choice among the site's backgrounds; source: the `.fieldbg` canvas of https://chb-bohrer.github.io, repo CHB-Bohrer.github.io).
  `app/bgfield.py` redraws the same recipe with PIL (seed 67, 1920x1080, transparent) into `app/static/bgfield.png` (committed, 450 KB); `style.CSS` puts it behind `.stApp` with `background-attachment: fixed`;
  `.streamlit/config.toml` has `[server] enableStaticServing = true` (serves app/static at /app/static/). Redraw: `scripts\run.cmd python scripts\make_bgfield.py [--seed --width --height]`. The dashboard colours (navy, teal)
  were NOT changed to the site's near-black/red palette. Restart the dashboard to see it (config change). Checked on port 8502: image served (200), no console errors, body text stays legible.

- Dashboard plot clean-up (2026-10-16; user: "anywhere else this could use improvement? ... fix all"). Method: render every page headlessly with Streamlit AppTest on a real run and save each figure by hooking
  `plt.close` (script kept out of the repo; save with `bbox_inches="tight"` like `style.fast_pyplot` does, otherwise titles/labels look cut off when they are not). Fixed: (1) `app/sweep_plot.py` replaces the one-line-per-combination
  sweep plot on 📈 Scans and the 🔗 E–P scans: each setting has ONE role (x axis / colour of lines / separate panels / held at a value that EXISTS in the sweep, sweeps are often not full grids), defaults by `default_roles`,
  widgets by `controls` (keys include the scan name), legend outside; all 60 (sweep x x-axis) combinations in $SIMLIVE_DATA/sweeps draw; (2) `style.fit_titles` wraps axes titles wider than their axes (called by `fast_pyplot` and
  `cached_figure`); (3) MS2 → Pol II legends above the panels, bigger heatmap numbers; (4) 🔬 Cells locus labels no longer overlap; (5) flat-field maps say "fit rejected" instead of an unexplained flat colour; (6) legends at readable
  size on Transcription/Probes. Tests: `tests/test_sweep_plot.py`. NOT looked at: sections behind toggles (Images build-up, players) and pages that draw no matplotlib figure (Metrics, Stage & fiducials, Library, Movie).
- Chromatin tab "Loop-extruding factors over time" (2026-10-16): the old plot (up to 60 factors x 2 legs x 8000 blocks as thin lines) was unreadable. Now `analysis.loop_coverage` (loops covering each 20 kb bin per
  time bin) as a heatmap + mean loop size (with a minutes axis from the assumed calibration and the boundary sites marked), and a zoom (`analysis.loop_traces`: genomic window x time window, at most 16 factors as shaded bands between the
  two legs, lines broken where a factor reloads; default 400 kb x 300 blocks, ~12-19 factors there). Tested in `tests/test_loop_view.py`; viewed on a real 8000-block library simulation. Streamlit does not
  hot-reload `app/chromatin_tab.py`: restart the dashboard to see edits.

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
- **Library run (2026-10-15 18:24 to 2026-10-18 19:10, FINISHED 100/100, no failures; detached via PowerShell `Start-Process`, output `C:\Users\cbohr\SimLiveData\library_loop_extrusion.out.log`):**
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
  page); the Metrics page shows `propagation` (stage 7b). Still open: CLAUDE.md is long and partly a log (condense when convenient).
