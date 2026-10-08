# Sim_Live_Cell_Analysis

Validate live-cell chromatin-locus analysis against simulations with known ground truth.
Pipeline: polymer dynamics (polychrom/OpenMM) -> cells + motion -> microscopy forward model ->
segmentation + tracking (Trackastra) -> locus-to-cell linking -> analysis -> validation.
See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for diagrams of how it fits together and the build roadmap, and
[CLAUDE.md](CLAUDE.md) for stages and conventions. **Status:** environment, GPU tests and image viewer are
done; stages 1-7 are not built yet.

## Quick start (Windows)

Everything runs through `scripts\run.cmd`, which uses the `simlive` conda env from any prompt
(no activation needed). Run from the repo root.

| I want to... | Command |
|---|---|
| Check the GPU setup (OpenMM CUDA, PyTorch, CuPy, Trackastra) | `scripts\run.cmd pytest -m gpu -v` |
| Run all tests | `scripts\run.cmd pytest` |
| **Open the dashboard (set parameters, run, see everything)** | `scripts\dashboard.cmd` then open http://localhost:8501 |
| Simulate moving cells and score the tracking (no GUI) | `scripts\run.cmd python scripts\run_tracking_demo.py` |
| ...with different settings | `scripts\run.cmd python scripts\run_tracking_demo.py --set motion.D_um2_s=0.05 --set seed=2` |
| **Single-cell loci run** (cells + two coloured loci, isolate each cell, locate loci per cell in parallel) | `scripts\run.cmd python scripts\run_tracking_demo.py --config configs\loci_demo.yaml` (~2 min), then open the dashboard's **🔬 Cells** tab |
| **Chromatin: run one polymer simulation** (GPU) | `scripts\run.cmd python scripts\run_chromatin.py --config configs\chromatin\smoke_test.yaml --seed 1 --library smoke` (~1 min; browse it in the dashboard's **🧬 Chromatin** tab) |
| Chromatin: re-check against known Mirny-lab results | `scripts\run.cmd python scripts\validate_chromatin.py all` (about 80 min on the GPU) -> [docs/chromatin_validation.md](docs/chromatin_validation.md) |
| Chromatin: build a library of simulations (resumable) | `scripts\run.cmd python scripts\run_chromatin_library.py --config configs\chromatin\library_loop_extrusion.yaml --library loop_extrusion --n 100` (default 10 Mb region, ~11 min each, ~18 h of GPU for 100; change `polymer.n_monomers` for another region size) |
| Compare settings over several seeds | `scripts\run.cmd python scripts\sweep.py --name mytest --seeds 1-5 --grid tracking.mode=greedy,ilp` (results in `%SIMLIVE_DATA%\sweeps\mytest\`) |
| View the tracked cells from that run | `scripts\run.cmd python scripts\view_run.py %SIMLIVE_DATA%\runs\<run_id> --labels stage4_segtrack\tracked_auto_masks.tif` |
| Make a fake demo movie (old, locus-only toy) | `scripts\run.cmd python scripts\make_demo_run.py` |
| **Look at a simulation (interactive window)** | `scripts\run.cmd python scripts\view_run.py %SIMLIVE_DATA%\runs\demo` |
| Save it as a GIF or a PNG instead | `scripts\run.cmd python scripts\view_run.py %SIMLIVE_DATA%\runs\demo --save out.gif` |

Viewer controls: slider or left/right arrows = frame, space = play/pause, `l` = toggle cell-mask overlay
(each persistent cell ID keeps its color), `i` = toggle cell-ID numbers. Works on any run directory containing
`stage3_microscopy/locus.tif`, `nucleus.tif`, `labels.tif` (T,Y,X or T,Z,Y,X).

Or activate the env yourself: `C:\Users\cbohr\miniforge3\condabin\conda.bat activate simlive`, then run the
commands without the `scripts\run.cmd` prefix.

## The dashboard

`scripts\dashboard.cmd` starts a local web app (this PC only, http://localhost:8501). The sidebar's **New
simulation** form sets motion, shape, imaging and tracker parameters and runs the pipeline (~20 s); **Run to view**
reopens any earlier run. Tabs: **Movie** (true cells | tracked cells side by side, large; identity swaps circled in red and
labeled, missed cells outlined in yellow; **Next error ⏭** jumps to the next frame with a mistake; optional raw-image
panel; per-frame error timeline), **Metrics** (all scores, with explanations, and a per-cell
identity table), **Trajectories** (true vs tracked paths, switches marked), **Config**. Stop it with Ctrl-C in its
window. The simulation runs in a separate process, so the dashboard itself never loads PyTorch.

## The tracking test (`run_tracking_demo.py`)

Cells move (diffusion + directed drift, collisions, optional division) -> nuclei are rendered with a
widefield forward model -> segmented -> tracked with Trackastra **twice** (on ground-truth masks = linking error
only; on automatic-segmentation masks = realistic) -> scored against the true cell IDs with Cell Tracking
Challenge metrics (TRA, DET, LNK, CHOTA, track purity / target effectiveness). Prints a table and saves
`stage7_validation/metrics.json`. All parameters (motion, density, optics, noise, frame interval, tracker) are in
[configs/tracking_demo.yaml](configs/tracking_demo.yaml); override any with `--set key=value`. The run folder name
encodes the config hash and seed, so the same config + seed always reproduces the same run.
Segmentation here is a classical threshold + watershed baseline; Cellpose/StarDist plug into
`stage4_segtrack/segment.py`, other trackers (Ultrack, TrackMate) into `stage4_segtrack/trackers.py`.

## Setup from scratch (already done on this PC)

1. Install Miniforge (`C:\Users\<you>\miniforge3`).
2. `conda env create -f environment.yml` (conda part; its pip step builds torch/Trackastra).
3. `powershell -File scripts\install_polychrom.ps1` — polychrom from GitHub without its optional Cython
   extension (needs MSVC; only used for knot simplification).
4. Copy `scripts\cuda_path.bat` to `<env>\etc\conda\activate.d\` so CuPy finds the CUDA headers.
5. `scripts\run.cmd pytest -m gpu -v` should show 5 passed.

## Layout

```
configs/     YAML parameter files and sweeps
data/        only .gitkeep placeholders. Real outputs live OUTSIDE the repo (so OneDrive does not sync gigabytes) in
             %SIMLIVE_DATA%, default C:\Users\<you>\SimLiveData (set by scripts\run.cmd):
             runs/<run_id>/stageN_<name>/, sweeps/<name>/, chromatin/ (polymer simulation library)
scripts/     run.cmd, view_run.py, make_demo_run.py, install scripts
src/simlive/ stage1_chromatin ... stage7_validation, calibration, io
tests/       pytest (GPU tests marked `gpu`)
```


### Chromatin library in the movie
`scripts\run.cmd python scripts\run_tracking_demo.py --config configs\loci_library.yaml` uses saved polymer simulations for the loci (needs a finished library named in `loci.library`; `block_duration_s` is a placeholder until MSD calibration).


### Watching a library run, and choosing simulations (🗂 Library tab)
Start the generator in one window (`scripts\run.cmd python scripts\run_chromatin_library.py --config configs\chromatin\library_loop_extrusion.yaml --library loop_extrusion --n 100`),
the dashboard in another (`scripts\dashboard.cmd`). The **🗂 Library** tab shows live progress (refreshes every 5 s; shows "Stalled" if the PC slept).
Below it you pick how cells get their simulation: random from the whole library (default), random from the ones you tick, or pin cell N to a given simulation.
The same choices work in config files: `loci.sim_ids: [...]`, `loci.pin: {1: loop_extrusion_seed1003}` (see `configs/loci_library.yaml`).


### Probe binding at the promoter and enhancer (🧪 Probes tab)
Each locus gets a stochastic model: free probes land (k_bind), scan (k_scan) and fall off (k_off); the number attached sets the locus
brightness, so binding noise shows up in the movie. Turn it on with `loci.probes.enabled=true` (or the sidebar checkbox in the single-cell loci form).
Explore and sweep parameters in the tab, or from a terminal: `scripts\run.cmd python scripts\sweep_probes.py --name koff --locus promoter --grid k_off=0.02,0.05,0.1,0.2`.
