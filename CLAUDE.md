# Sim_Live_Cell_Analysis

Repo for validating live-cell microscopy analysis of chromatin loci against simulations with known ground truth.
Package name: `simlive` (src layout). Owner is a biophysicist (GitHub: CHB-Bohrer; repo
https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis, public). See README.md for how to run everything.

## Session handoff rule (always follow)
**Before ending a work session, or after any significant change, update this file and README.md without being asked:**
what changed, new commands/files, gotchas hit, findings, open questions, next steps. Then commit. A new session starts
from this file and the code alone, so anything not written here is lost. Keep the "Status" and "Findings" sections
below true; delete stale statements rather than appending contradictions.

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
nuclear channel only (PSF blur, noise, bleaching, per-nucleus texture + nucleoli), stage 4 (Cellpose default,
threshold-watershed baseline, Trackastra `general_2d` behind the `Tracker` adapter), stage 7a (CTC metrics via
traccuracy + own identity-switch analysis), `scripts/run_tracking_demo.py`, `scripts/sweep.py`, dashboard.
NOT built: stage 1 (chromatin / polychrom), locus + MS2 channels, calibration, stage 5 (locus-to-cell linking),
stage 6 (analysis), stage 7b (propagation of tracking errors into locus results), Ultrack/TrackMate adapters,
density and frame-interval scans.

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
- **Run layout:** `data/runs/<run_id>/stageN_<name>/...`. `run_id` encodes the config hash + seed. Outputs are
  gitignored; they must be reproducible from config + seed.
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
