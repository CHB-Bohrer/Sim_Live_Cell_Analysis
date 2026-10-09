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
  user's server once). Tabs: Movie, Metrics, Trajectories, Scans, Cells, Chromatin, Config.

## Status (update this section every session)
Built: stage 2 (cell motion, collisions, division, irregular time-varying nuclear shapes, `speed_scale`), stage 3
(nuclear channel + one colour channel per locus: PSF blur, noise, bleaching, per-nucleus texture + nucleoli), stage 4
(Cellpose default, threshold-watershed baseline, Trackastra `general_2d` behind the `Tracker` adapter), stage 5a
(per-cell isolation: `stage5_linking/isolate.py`, `celldata.py`), stage 7a (CTC metrics via traccuracy + own
identity-switch analysis), `scripts/run_tracking_demo.py`, `scripts/sweep.py`, dashboard (tabs incl. Scans, Cells).
STAND-INS (to be replaced): stage 1 = `stage1_chromatin/toy_loci.py` (2 loci per nucleus doing confined diffusion in
the nucleus frame, NOT polychrom); stage 6 = `stage6_analysis/localize.py` (brightest-spot localizer, placeholder for
the user's analysis). The loci/isolation pipeline runs via `configs/loci_demo.yaml`.
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
NOT built: MSD calibration (real nm/s), surrogate distance process for long movies, 3D z-stacks, the real stage 5 (link
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
  Do not draw per-frame matplotlib figures for animation. Full-app reruns (e.g. sidebar selectbox) reset the active tab.
- The stand-in loci (`toy_loci.py`) are INVENTED, not from polychrom or data: confinement 0.1 x nuclear radius, relaxation
  time 120 s, 600 photons/locus/frame are guesses. Never present their statistics as physical results.
- Gotchas: tifffile treats stacks of exactly 3-4 frames as RGB, so always write with `photometric="minisblack"`;
  Cellpose needs the right `diameter_px` at fine pixel sizes (`diameter_px: auto` in loci_demo.yaml; without it DET fell
  to 0.946 at 130 nm/px); Streamlit does not reliably hot-reload helper modules (`app/cells_tab.py`), so restart the
  dashboard after editing them.

## Findings so far (cell tracking only; simulated 2D nuclei, 5 seeds per setting)
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
  only (`addopts = -m 'not gpu'`), `pytest -m gpu` = GPU tests, `scripts\run_tests.cmd` runs both. 38 tests at 2026-10-08. Prefer small deterministic tests with a
  known analytic answer (e.g. free diffusion MSD).
- Python 3.12, type hints on public functions, no wildcard imports.
