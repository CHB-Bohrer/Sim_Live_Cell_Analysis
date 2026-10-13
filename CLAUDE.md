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
- MSD calibration targets (D, alpha, length/time scales) not yet given.
- Still needed from the user: real frame interval, typical cell speed, cell density, whether cells divide, and ideally
  a real movie with some tracked cells to check how well the simulation matches reality.

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
  (`configs/chromatin/library_*.yaml`: 1000 equilibration + 2000 production blocks x 750 steps) = ~11 min per 10 Mb
  simulation, ~18 h per 100-simulation library, ~7 GB per library at 10 Mb. The libraries have NOT been started; each
  occupies the GPU for many hours, so get the user's go-ahead and region size first. Whether 1000 equilibration blocks
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
  Transcription, MS2 → Pol II, Probes, Library, Chromatin, Config. DARK theme since 2026-10-12 (user asked for a darker UI): colours in
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
Library wiring (2026-10-08): `stage1_chromatin/library_loci.py` (`loci.source: library`, config `configs/loci_library.yaml`; one simulation per cell, random rotation + time offset, polymer coordinates converted to true nm with `nm_per_unit` (PLACEHOLDER 50; NOT rescaled to the nucleus) and placed at a random territory centre per cell (re-drawn until all loci stay inside the nucleus in every frame), `block_duration_s` is a PLACEHOLDER) is unit-tested with a fake library on CPU but NOT yet run end-to-end with a real library or on the GPU machine.
Dashboard redesign + library monitor (2026-10-08, built on the cloud container with fake data, NOT yet seen on the Windows machine
with a real GPU run): `.streamlit/config.toml` + `app/style.py` (one theme, CSS, matplotlib defaults), new 🗂 Library tab
(`app/library_tab.py`): live monitor (refreshes every 5 s) of `<library>/status.json` (written by `run_chromatin_library.py`) and the
running simulation's `progress.json` (written by `polymer3d.py` every 20 blocks): state badge (running / stalled / finished /
stopped), n done, ETA, GPU load via nvidia-smi, recent log; plus a picker to choose which simulations a movie uses.
Choosing simulations: `library.assign_simulations(..., pool=, pinned=)`; config `loci.sim_ids` (restrict the pool) and
`loci.pin: {cell_id: sim_id}`; the sidebar's loci form has a source radio (stand-in / library) that uses the Library tab's choice.
The picker's preview uses the movie's own random stream (`SeedSequence(seed).spawn(4)[2]`), so it matches the real assignment when no cell divides.
Untested: `status.json` writing needs a real library run (GPU).
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
hand-off to a movie (sidebar checkbox in the loci form). Seen only on the cloud container with no real runs.
Transcription / MS2 model (2026-10-09; stage 1c): `stage1_chromatin/transcription.py` (details and the user's decisions: `docs/PLAN_transcription.md`).
Promoter OFF -(k_on)-> ON -(k_off)-> OFF, Pol II initiation k_init while ON, Pol II ramp over the MS2 cassette then plateau until released
(random dwell); MS2 signal = loops carried by all Pol II. The promoter-enhancer distance d(t) (3D, from the library trajectories at full block
resolution, via `simulate_library_loci(..., return_traces=True)`; toy loci use their 2D lab distance) multiplies ONE chosen rate
(`coupled_rate`: k_on default | k_init | k_off) by 1 + (fold-1) f(d), f = sharp contact (default) | hill | exp | none (control). ALL rates are per second
of real time, so the frame interval can change freely (library loci are interpolated between blocks; `check_sampling` warns when too coarse).
Exact Gillespie + exact expectations (`expected_on`, `expected_ms2`); 12 tests in `tests/test_transcription.py` (constant-distance theory, two-valued
distance vs exact expectation, frame-interval independence, null/positive correlation with contact, interpolation, MS2 spot in the nuclear
channel). `loci.transcription.enabled: true` (configs loci_demo/loci_library) writes `stage1_chromatin/transcription.csv` + `transcription_events.csv`
and renders the spot at the promoter into the nuclear channel (`render_nuclei(..., spots=)`; pipeline now runs stage 1 BEFORE rendering nuclei; seed
stream ss[5]). If the movie is longer than the saved polymer trajectory, the trajectory is repeated back and forth with a warning (correlations beyond
its length are artificial); the planned surrogate distance process (PLAN step 5) is NOT built. Dashboard ⚡ Transcription tab (`app/transcription_tab.py`): nine
numbered steps (sampling, distance, coupling, promoter, Pol II kymograph, MS2 signal with exact expectation, simulated microscope tiles, parameter sweep,
hand-off to a movie + sidebar checkbox). Checked in a browser on the cloud container with a FAKE library only; the pipeline was run through stage 4
segmentation (trackastra not installed there). LITERATURE NUMBERS ARE UNVERIFIED search-summary values (primary papers blocked from the build machine); all
MS2 brightness, coupling strength, initiation rate and nm_per_unit / block_duration_s are placeholders.
NOT built: MSD calibration (real nm/s), surrogate distance process for long movies, 3D tracking / 3D cell isolation, the real stage 5 (link
locus tracks to cells beyond isolation), the user's stage 6 analysis, stage 7b (propagation of tracking errors into
locus results; e.g. run stage 6 on 'truth' vs 'tracked_*' identity sources and compare), Ultrack/TrackMate adapters,
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
  only (`addopts = -m 'not gpu'`), `pytest -m gpu` = GPU tests, `scripts\run_tests.cmd` runs both. 130 CPU + 6 GPU tests at 2026-10-13. Prefer small deterministic tests with a
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
  `score_localization_3d`; scores `truth_3d`. Pipeline stream layout is now `SeedSequence(seed).spawn(8)`: 0 motion, 1 imaging, 2 loci, 4 probes, 5 transcription, 6 heights, 7 free (entropy_for still uses stream 1).
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
    `io.runs.load_config`, --set still wins; sidebar "Microscope" radio). `vutara_vxl` = Bruker Vutara VXL in 2D bypass mode. Every value is
    sourced or ASSUMED (docs/INSTRUMENTS.md):
    - 60x/1.3 NA silicone, ~100 nm px, Orca Fusion BT QE 0.95, top-hat (vignette 3% ASSUMED), active focus (defocus 30 nm), sequential colours;
    - bleed-through FROM SPECTRA (FPbase EGFP: 0.0042 cross-excitation at 555 nm x 0.062/0.717 in-band emission = 0.036%, vs ~8.6% if excited
      simultaneously);
    - biplane (2 planes 600 nm apart, 50:50) NOT built.
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
      starts/min (was 0.12-1.39). ~20-30x slower than gibbs_shift, so NOT the pipeline default: `deconvolve_run` still uses gibbs_shift; for dense genes run
      `scripts/deconvolve_ms2.py --sampler gibbs_pt` (files `_gibbs_pt`; the page's sampler table lists all three) and check with
      `scripts/check_mcmc.py --samplers gibbs_pt --tag _pt --workers 6`. Headline scores vs the simulation are unchanged by tempering (found 0.83, false 0.20,
      rate r 0.80; per-frame posterior means differ by about as much as the other samplers differ from each other), so earlier totals stand; intervals and
      single-start timing of dense bursts should come from the tempered chains. Open: make tempering the default for dense genes? fewer replicas?
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
- Not done: the 100-simulation library (waiting for go-ahead; the live monitor is still only tested with fake files + finished libraries);
  one library selector exists in the monitor and a separate one in the picker (confusing; not merged); the header metrics of the selected run
  sit above the Transcription/Probes/Library/Chromatin tabs although unrelated to them; default seed-1 cell in Transcription step 5 shows no Pol II
  (empty plots); with 50 nm/unit and 150 nm contact the smoke distances rarely touch, so coupling looks inert (placeholder scale).
- Literature checked against PMC full texts, see docs/PLAN_transcription.md: TFF1 ON 16 min (not 6), burst period 66-86 min, 1.5 transcripts/burst;
  Pol II 2.6 kb/min, 3' dwell 116 s (our default 30 s). Defaults not changed: needs the user's choice of gene.

- Players (app/player.py): every in-browser movie now shows 'Loading frames… n / N' with a progress bar and a disabled button, then 'Ready: press Play'. Chromatin player: the grey circle is the simulation's confining sphere (the territory, NOT the nucleus); caption and in-image label say so. The territory is ~1-1.5 um across at 50 nm/unit vs a 12 um nucleus.

- Display scaling (2026-10-11): single-cell and all-cells players scale each locus colour once for the whole movie from the cell's own mask region (fastview.movie_limits); before, each frame was scaled to the crop maximum, so cosmic rays and fiducial beads in the crop corners made the real spots look missing (spot < 30% of crop max in 70-92% of frames in 3D runs). Measured locus D = MSD/(4t) at 10 s: library loci 0.0017-0.0019 um2/s (falls to 0.0011 at 100 s, nm_per_unit 50 and block_duration_s 10 are PLACEHOLDERS), toy loci 0.0029; the nucleus itself moves 0.43 um per frame.
