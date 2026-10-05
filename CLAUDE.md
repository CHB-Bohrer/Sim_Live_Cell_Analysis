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

**Pushing to GitHub (public repo):** commit locally as you go, but never push on your own. When something major is
done and verified working (new stage, new feature, important result) and the user has not mentioned pushing, ASK
them whether to push ("This works and is committed locally; push it to GitHub?"). Push immediately only if they
say so. Before every push: run `git status -sb`, check for secrets (password/token/key strings) and large files, and
confirm `data/runs/` and `data/sweeps/` are still ignored.

## The experiment being simulated (from the user; some details still open)
- Loci: LacO/TetO-array-like labels plus MS2 bursting; reported as three different colors. Unclear whether MS2 is a
  third color or one of the two loci — confirm. The nuclear NLS-GFP signal is the same color as the MS2 spots.
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

## Chromatin simulation plan (stage 1, agreed with the user 2026-10-05; not built yet)
- Simulate chromatin SEPARATELY from the microscopy pipeline, as a reusable library of specific polymer simulations
  (Mirny-lab style, with polychrom/OpenMM on the GPU), saved with their parameters, seed and versions. Stage 1 of a
  movie then assigns one specific saved simulation to each cell; the microscopy labels/fluorophores come later.
- Before using them, reproduce known Mirny-lab results to verify the setup (e.g. contact probability P(s) scaling,
  loop-extrusion contact maps, compartments, sub-diffusive locus MSD). Check target numbers against the papers; do not
  rely on memory.
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
  user's server once). Tabs: Movie, Metrics, Trajectories, Scans, Config.

## Status (update this section every session)
Built: stage 2 (cell motion, collisions, division, irregular time-varying nuclear shapes, `speed_scale`), stage 3
(nuclear channel + one colour channel per locus: PSF blur, noise, bleaching, per-nucleus texture + nucleoli), stage 4
(Cellpose default, threshold-watershed baseline, Trackastra `general_2d` behind the `Tracker` adapter), stage 5a
(per-cell isolation: `stage5_linking/isolate.py`, `celldata.py`), stage 7a (CTC metrics via traccuracy + own
identity-switch analysis), `scripts/run_tracking_demo.py`, `scripts/sweep.py`, dashboard (tabs incl. Scans, Cells).
STAND-INS (to be replaced): stage 1 = `stage1_chromatin/toy_loci.py` (2 loci per nucleus doing confined diffusion in
the nucleus frame, NOT polychrom); stage 6 = `stage6_analysis/localize.py` (brightest-spot localizer, placeholder for
the user's analysis). The loci/isolation pipeline runs via `configs/loci_demo.yaml`.
NOT built: real stage 1 (polychrom/OpenMM) and MSD calibration, MS2 bursting, 3D z-stacks, the real stage 5 (link
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
- **Tests:** pytest, in `tests/`. GPU tests marked `@pytest.mark.gpu`. Prefer small deterministic tests with a
  known analytic answer (e.g. free diffusion MSD).
- Python 3.12, type hints on public functions, no wildcard imports.
