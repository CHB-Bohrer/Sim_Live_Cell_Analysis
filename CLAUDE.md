# Sim_Live_Cell_Analysis

Repo for validating live-cell microscopy analysis of chromatin loci against simulations with known ground truth.
Package name: `simlive` (src layout). Owner is a biophysicist (GitHub: CHB-Bohrer; repo
https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis, public). See README.md for how to run everything.

## Owner's standing preferences (added 2026-11-27)
- **Short status page.** `STATUS.md` is the one-page view of where the project stands and what is next. This file holds the rules for every session; the detailed record is in `docs/notes/` (index below).
  Rewrite it at the end of a session (under 40 lines; replace, do not append).
- **Parked thoughts.** When the owner mentions an idea that is not for now ("park this", "note for later", an aside outside the
  current task), add it to `THOUGHTS.md` in their words with a concrete **Revisit when** condition, and say in one line which
  condition you recorded. Do not act on it. List the entries with `scripts\run.cmd python scripts\thoughts.py` at the start of a session and at wrap-up;
  when the current work meets an entry's condition, raise it in a sentence or two and let the owner decide, then record it as
  adopted or dropped. Never delete an entry.
- **Do the obvious next thing, at most two.** When the next step is clear, do it instead of asking: it follows directly from the
  task or is tagged `[auto]` in `STATUS.md` (`scripts\run.cmd python scripts\next_steps.py --auto`), git can undo it, it needs no long GPU run, money or
  new dependency, and it does not change the scientific direction or what a result means. After the task the owner asked for, do
  at most TWO such follow-ups per session and write every other one into `STATUS.md` under Next as `- [auto] ...`; they are
  picked up next session, or all at once when the owner says "go ahead with the next steps". Ask first for anything tagged
  `[ask]`, a choice between scientific alternatives, deleting data, publishing, or changing a tolerance or expected value so a
  check passes. Note self-started work in `STATUS.md` under Last session with `(auto)`.
- **Git.** Commit and push without asking for confirmation, following "Pushing to GitHub" below (straight to `main`, the checks before every push still apply). In force until the owner says otherwise.

## Where the project lives (details: `docs/notes/location_and_disk.md`)
- Everything is in `C:\Users\cbohr\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\`: `code\` = this git repo,
  `simulation_data\` = runs, sweeps, chromatin library, tutorial files (`SIMLIVE_DATA`; found relative to the repo by `scripts\run.cmd`).
- `C:\Users\cbohr\SimLiveData` is a JUNCTION to `simulation_data`. Do not delete the junction's target through it.
- The conda envs import the code through an editable install; if the repo moves, fix the `.pth` line (see the notes) or the env imports nothing.
- Open Claude Code sessions in the `code` folder. Dropbox locks files while uploading (`os.replace` can fail with WinError 5: retry).

## Session handoff rule (always follow)
**Before ending a work session, or after any significant change, without being asked:**
1. rewrite `STATUS.md` (where things stand, what was done, what is next);
2. add what changed (new commands / files, gotchas hit, findings, open questions) as a dated section at the top of the matching
   file in `docs/notes/`; a new topic gets a new file and a row in the index below;
3. update README.md if how to run something changed; then commit and push.
A new session starts from this file, `STATUS.md` and the code alone, so anything not written down is lost. Keep THIS file to
rules that apply in every session (aim: under 200 lines); results, history and one-topic gotchas belong in the notes. Delete
stale statements rather than appending contradictions. `docs/ARCHITECTURE.md` holds Mermaid diagrams
(pipeline, validation logic, roadmap) whose boxes are coloured built / stand-in / planned: recolour or edit them when a
stage's status changes (GitHub renders them; to check syntax, render the blocks with mermaid.js in a browser).

**Public status page (added 2026-11-08):** `status.json` at the repo root feeds the page at https://chb-bohrer.github.io/Sim_Live_Cell_Analysis/
(and the Fun Code page on the lab site). When a stage changes state (built / stand-in / planned), or you finish a milestone, update `status.json` to match the
Status section and `docs/ARCHITECTURE.md` section 6: node states, notes, `as_of` date. Keep it short, accurate and public-safe (no unpublished data, no credentials).

**Pushing to GitHub (public repo):** the user said (2026-10-09) to push STRAIGHT TO `main` every time: after committing, run
`git push origin HEAD:main` (fast-forward only; never force-push; never open a PR unless asked; the working branch may be pushed too).
Before every push: run `git status -sb`, check the staged diff for secrets (password/token/key strings) and large files, and confirm
`data/runs/` and `data/sweeps/` are still ignored. If a check fails, do not push; tell the user. (Earlier text said "push the working
branch, not main": that was a misreading of the user's intent.)

## Topic notes (read only the one the task touches)
| File in `docs/notes/` | Covers | Read when |
| --- | --- | --- |
| `location_and_disk.md` | Where the project lives on disk; disk-space rules | moving folders, paths that fail, disk usage, TIFF reading/writing |
| `status_detail.md` | Detailed build status per stage; tracking findings | you need exactly what a stage does or does not do yet |
| `chromatin.md` | Chromatin library, polychrom engine, polymer calibration, stage 7b | touching stage 1, polychrom, the library, nm/second calibration, stage 7b |
| `imaging.md` | Image-error model, z-stacks, calibration, piezo and fiducials, biplane | touching stage 3, the camera/error sources, 3D, calibration images, biplane |
| `spot_finding.md` | Spot finder, MS2 trace, bleaching, deconvolution, samplers, multi-copy linking, isolation | touching stage 5a/6: localization, MS2, bleaching, Pol II deconvolution, copies and strays |
| `ep_curve.md` | E-P distance to P(burst ON) curve under measurement errors | touching the E-P curve fit, its scans or acquisition model |
| `audits.md` | Logic audit (2026-10-26) and improvement pass (2026-11-09/10) | asking whether a result predates a fix, or what was re-run |
| `dashboard.md` | Dashboard pages, performance, plots, walk-through and the 2026-10-20 fix round | editing `app/`, plots, pages, or chasing a dashboard crash |
| `tutorial.md` | Narrated tutorial, voice, playback, video releases | editing the tutorial, its screenshots, voice clips or video |

Longer write-ups with the numbers: `docs/ARCHITECTURE.md`, `SPOT_FINDING.md`, `EP_CURVE.md`, `ERROR_EFFECTS.md`, `IMAGE_ERRORS.md`, `CALIBRATION_AND_3D.md`,
`MCMC_CONVERGENCE.md`, `INSTRUMENTS.md`, `chromatin_validation.md`, `PLAN_transcription.md`, `PLAN_three_signals_and_homologs.md`, `docs/tutorial/UPDATING.md`.

## Rules learned the hard way (the story behind each is in the note named)
- Do not edit `src/`, the library configs or `io/runs.py` while a pipeline run, re-analysis or library run is going: worker
  processes re-import modules from disk (`dashboard.md`, `chromatin.md`).
- Code that runs before the tracker must avoid numpy matmul / BLAS calls: they load a second OpenMP runtime and `import torch`
  then fails (`dashboard.md`). Isolate any torch tool that pulls in scikit-image in its own process (Environment below).
- Never re-point existing OpenMM bonds with `updateParametersInContext`; check couplings physically, not only by summary
  statistics (`chromatin.md`).
- Limits set in advance stay where they are: insulation target 1.5, R-hat < 1.05, ESS >= 400. A miss is reported, not tuned away.
- Single runs are not evidence: compare over >= 5 seeds. Validate a metric on null data before trusting it (`chromatin.md`).
- Every instrument, probe, MS2, stray-spot and polymer-calibration number is a placeholder or an assumption until the owner
  supplies measurements. Never present their statistics as physical results (`status_detail.md`).
- Images: write with `io.images.write_tiff`, read with `open_stack`; always `photometric="minisblack"` (`location_and_disk.md`,
  `spot_finding.md`). `errors.INDEX` is frozen: new error sources get the next number (`imaging.md`).
- Write regex source, JS strings and docs with the Edit / Write tools, not through heredocs or Python string patches (`\b`, `\r`
  and `\n` were silently mangled three times; `audits.md`, `chromatin.md`, `tutorial.md`). Keep this file UTF-8.
- Videos go to GitHub Releases, never into the git tree; publishing a release needs the owner's OK first (`tutorial.md`).
- Do not quote paper numbers from memory; re-check them (`chromatin.md`).

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
  `C:\Users\cbohr\Dropbox\Live_Cell_Microscopy_Investigations\Sim_Live_Cell_Analysis\simulation_data` = the folder `simulation_data` NEXT TO the repo folder (set in `scripts\run.cmd`;
  `simlive.io.runs.DATA_ROOT` does the same when the env var is unset). See "Project location" at the top. `run_id` encodes the config hash + seed. Outputs are gitignored and must be
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
  only (`addopts = -m 'not gpu'`), `pytest -m gpu` = GPU tests, `scripts\run_tests.cmd` runs both. 225 CPU + 7 GPU tests at 2026-11-10. Prefer small deterministic tests with a
  known analytic answer (e.g. free diffusion MSD).
- Python 3.12, type hints on public functions, no wildcard imports.

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

## Pinned versions (do not break)
- Exact working versions are saved in `environments/` (see `environments/PINNED.md`): env export, explicit conda list,
  pip freeze, polychrom commit `11a870c` (installer is pinned to it), model-weight SHA-256s, a frozen conda env copy
  `simlive-frozen-2026-10-05`, and git tag `env-2026-10-05`. NEVER run `pip install -U` / `conda update` in `simlive`;
  test upgrades in a separate env and save new lock files first.
