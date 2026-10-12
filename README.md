# Sim_Live_Cell_Analysis

Validate live-cell chromatin-locus analysis against simulations with known ground truth.
Pipeline: polymer dynamics (polychrom/OpenMM) -> cells + motion -> microscopy forward model ->
segmentation + tracking (Trackastra) -> locus-to-cell linking -> analysis -> validation.
See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for diagrams of how it fits together and the build roadmap, and
[CLAUDE.md](CLAUDE.md) for stages and conventions. **Status:** see the Status section of CLAUDE.md (stages 1-5, 7a built; the chromatin library, image-error model and dashboard are in use; the real stage 6 analysis is still to come).

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


### Transcription and the MS2 spot (⚡ Transcription tab)
The promoter-enhancer distance (from a saved polymer simulation, 100 kb apart by default) sets how readily the promoter switches ON; ON promoters
fire Pol II that carry MS2 loops, which appear as a bright spot at the promoter in the nuclear channel. Explore it step by step in the tab, or run it in
a movie with `loci.transcription.enabled=true` (sidebar checkbox: "Model distance-dependent transcription"), e.g.
`scripts\run.cmd python scripts\run_tracking_demo.py --config configs\loci_library.yaml --set loci.transcription.enabled=true --set acquisition.frame_interval_s=10`.
All kinetic rates are per second, so any frame interval works. Sweep a setting in the tab, or across whole movies with `scripts\sweep.py --grid loci.transcription.params.fold=3,12,48`.

### Sidebar checkboxes and library movies
The single-cell loci form now really passes the **probe binding** and **transcription** settings from their tabs to the run, and library movies place the
two loci at the positions set in the Transcription tab (default: 100 kb apart in the middle of the chosen library's region).


### Image errors and the 🖼 Images tab
The microscope model (stage 3) has 22 switchable sources of image error: stage drift, chromatic shift between colours, diffraction blur, focus drift, out-of-focus haze,
exposure motion, bleaching, uneven illumination, flicker, background, autofluorescence, colour cross-talk, quantum efficiency, shot noise, dark current, hot pixels,
cosmic rays, pixel-gain pattern, read noise, row noise, offset pattern and digitisation/saturation. They are set in the `imaging_errors:` block of the configs (all
values are placeholders), or with the sidebar's **Image errors** controls; `imaging_errors.enabled: false` gives the simple original model. Every source is
explained in [docs/IMAGE_ERRORS.md](docs/IMAGE_ERRORS.md) and in the dashboard's **🖼 Images** tab, which shows each recorded frame next to the ground truth that produced it
(true masks, true locus positions, drift/focus/light level of that frame), builds an image up one error at a time, and measures how much each error matters.


### z-piezo errors and fiducial beads
The z-piezo that steps through the planes has random, bowed, hysteretic and lagging position errors, and beads on the slide track drift and focus in every image. The 📍 **Stage & fiducials** tab explains
each piece, shows where the planes really were, what the beads can and cannot measure, and what correcting with them or with the piezo's own position log does.

### 3D z-stacks and calibration
Tick **Record 3D z-stacks** in the sidebar's loci form (or use `configs\loci_3d.yaml`): every frame becomes a stack of planes, nuclei are 3D ellipsoids, the loci have a height, and
the microscope blurs along z too. The loci are located in 3D and the 🖼 Images tab shows slices, projections and x-z / y-z side views next to the ground truth.
**Calibration**: the simulation also takes dark frames, a flat-field slide and a bead slide with the same instrument, estimates the instrument from them (gain, noise, hot pixels, illumination,
PSF, colour offsets, drift) and compares each estimate with the true value; the movie is then analysed again with the calibration applied. See the 🎯 Calibration tab,
`scripts\run.cmd python scripts\run_calibration.py --run <run folder>` and [docs/CALIBRATION_AND_3D.md](docs/CALIBRATION_AND_3D.md).

### Finding the loci and the MS2 signal: the 🔎 Spot finding page
Stage 6 tracks the promoter and the enhancer (one colour each) in every frame of every isolated cell with the spot pipeline of Chen, Levo, Zoller & Gregor 2025
(*Nat. Struct. Mol. Biol.*): median filter, difference of Gaussians with a low threshold (several candidates), the most likely path through the candidates over the whole
movie (a Viterbi algorithm with a diffusion potential, which can pick the spot up again after long gaps), a Gaussian fit, and their ellipse spot intensity. The MS2
transcription signal is read in the nuclear channel at the tracked promoter in every frame (a trace per cell). Spot intensities are corrected for photobleaching, with
the background removed before the decay is measured. Settings: the `locus_analysis:` block of the configs, or the sidebar's **Spot finding** controls. The page shows
every step on any cell, spot and frame, the error of every frame against the truth, the MS2 trace, the bleaching correction and a comparison of methods. From a
terminal (re-runs only stage 6 on a finished run, ~1-3 min; also adds MS2 traces and the bleaching correction to older runs):
`scripts\run.cmd python scripts\compare_localizers.py --run <run folder or name>`. Method, sources and measured numbers: [docs/SPOT_FINDING.md](docs/SPOT_FINDING.md).

**Uneven illumination, per-cell bleaching, cells entering the field.**
- **Simulation:** with `imaging_errors.bleaching.per_cell: true` (on in the loci configs) every cell bleaches by the light it has received. That dose follows
  the uneven illumination along its path, a cell gets nothing while outside the lit area (`lit_margin_um`), and daughters inherit the dose.
  `geometry.outside_um` lets cells live outside the imaged field and walk in and out (sidebar: *Room outside the field*). Truth per cell is in
  `stage3_microscopy/bleach_cells.csv`.
- **Analysis (`stage6_analysis/flatfield.py`):** estimates the illumination map and the bleaching from the movie itself, with no flat-field image. The same
  nucleus is seen at many places, so brightness changes that follow position reveal the illumination, and decline that follows accumulated light reveals
  bleaching. The illumination term is kept only if it predicts held-out cells better than even light.
- **Results:** shown on the Spot finding page, part 3c.

**How bright must the MS2 spot be?** `scripts\run.cmd python scripts\sweep_ms2_brightness.py` sweeps R = (one Pol II's MS2 photons) / (one nuclear GFP pixel),
3 seeds per value, tracking skipped (`tracking.skip: true`). Results: Spot finding page, part 3d.

**What each image error does, and how it is corrected (🖼 Images page, part 2).** Pick one of the 24 error sources from the list:
- **Explanation:** what it is, what it does to the picture, and how it is corrected. There are three possible routes: by the analysis
  itself, with calibration images (dark frames, flat-field slide, bead slide), or from the movie alone.
- **Images:** the same frame without the error, with it, and corrected.
- **Measured effect on the results**, with and without each correction: wrong picks, position error, distance error, intensity error
  and drift, MS2 trace.

The measurement re-renders the run's movie once per source and analyses it each way (about 30 min for all sources, in the background
from the page, or `scripts\run.cmd python scripts\error_effects.py --run <run>`). The movie-only corrections (`stage6_analysis/selfcal.py`)
build the same calibration object the calibration images give:
- background, offset and hot-pixel maps from the frames in which no nucleus covers a pixel;
- the illumination from the moving nuclei;
- the colour shift from the mean enhancer-promoter vector, applied only when it is significant (needs many cells);
- colour cross-talk: the bleed-through of the bright nuclear / MS2 light into a locus colour, measured and subtracted;
- camera row offsets from the pixels next to each cell;
- per-cell bleaching.

Measured results and what they mean: [docs/ERROR_EFFECTS.md](docs/ERROR_EFFECTS.md). In short:
- **Positions:** the colour shift dominates, and a bead slide fixes it.
- **Spot intensities:** the bleed-through of a bright MS2 spot into the promoter's colour dominates; unmixing from the movie fixes it.
- **Best:** use both, calibration images for what they give and the movie for what only it gives.

**From the MS2 trace to Pol II starts (🧮 MS2 → Pol II page).** The Bayesian deconvolution of Chen et al. 2025 turns each cell's MS2 trace
into the moments Pol II started on the gene (with uncertainty) and into bursts. The page explains it step by step for new users:
units, what one Pol II looks like, noise, one cell against the truth, and the whole run. The burst rule's settings can be changed live.
`scripts\run.cmd python scripts\deconvolve_ms2.py --run <run>`; details and measured accuracy in [docs/SPOT_FINDING.md](docs/SPOT_FINDING.md).

The dashboard is now dark. [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) has updated diagrams of how it all works and of what is built.

**Not built yet:** using the promoter, the enhancer and the MS2 spot together to choose spots, and separating the copies of a locus on different
chromosomes. The plan is in [docs/PLAN_three_signals_and_homologs.md](docs/PLAN_three_signals_and_homologs.md).

The dashboard now runs only the page you have open (a navigation bar replaces the tabs), which made reruns ~5-10x faster; the run summary has a one-line description
of every image-error source.
