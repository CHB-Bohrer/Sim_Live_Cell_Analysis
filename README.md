# Sim_Live_Cell_Analysis

Validate live-cell chromatin-locus analysis against simulations with known ground truth.
Pipeline: polymer dynamics (polychrom/OpenMM) -> cells + motion -> microscopy forward model ->
segmentation + tracking (Trackastra) -> locus-to-cell linking -> analysis -> validation.
See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for diagrams of how it fits together and the build roadmap, and
[CLAUDE.md](CLAUDE.md) for stages and conventions. **Status:** see the Status section of CLAUDE.md (stages 1-5, 7a built; the chromatin library, image-error model and dashboard are in use; the real stage 6 analysis is still to come).

## New here? Take the narrated tutorial
Watch the video **[narrated MP4 on the Releases page](https://github.com/CHB-Bohrer/Sim_Live_Cell_Analysis/releases/tag/tutorial-video-2026-10-25)** (49 min, 720p or full size), or double-click **[docs/tutorial/index.html](docs/tutorial/index.html)** (best in Microsoft Edge): about 50 minutes of animated, narrated scenes covering the
pipeline, starting the dashboard, both simulation forms, and every dashboard page, each ending with a short "Your turn" list. Voice = recorded clips made locally with Kokoro (af_heart) (default speed 2x; without the clips, your browser's
speech synthesis, where Edge's "Natural" voices sound best); the MP4 predates the newest scenes. Details and how to regenerate it in [docs/tutorial/README.md](docs/tutorial/README.md).

## Quick start (Windows)

Everything runs through `scripts\run.cmd`, which uses the `simlive` conda env from any prompt
(no activation needed). Run from the repo root.

| I want to... | Command |
|---|---|
| Check the GPU setup (OpenMM CUDA, PyTorch, CuPy, Trackastra) | `scripts\run.cmd pytest -m gpu -v` |
| Run all tests | `scripts\run.cmd pytest` (CPU; includes opening every dashboard page on the saved showcase runs) and `scripts\run.cmd pytest -m gpu` (GPU: OpenMM, Cellpose, Trackastra), or `scripts\run_tests.cmd` for both |
| **Open the dashboard (set parameters, run, see everything)** | `scripts\dashboard.cmd` then open http://localhost:8501 |
| Simulate moving cells and score the tracking (no GUI) | `scripts\run.cmd python scripts\run_tracking_demo.py` |
| ...with different settings | `scripts\run.cmd python scripts\run_tracking_demo.py --set motion.D_um2_s=0.05 --set seed=2` |
| **Single-cell loci run** (cells + two coloured loci, isolate each cell, locate loci per cell in parallel) | `scripts\run.cmd python scripts\run_tracking_demo.py --config configs\loci_demo.yaml` (~2 min), then open the dashboard's **🔬 Cells** tab |
| **Chromatin: run one polymer simulation** (GPU) | `scripts\run.cmd python scripts\run_chromatin.py --config configs\chromatin\smoke_test.yaml --seed 1 --library smoke` (~1 min; browse it in the dashboard's **🧬 Chromatin** tab) |
| Chromatin: re-check against known Mirny-lab results | `scripts\run.cmd python scripts\validate_chromatin.py all` (about 80 min on the GPU) -> [docs/chromatin_validation.md](docs/chromatin_validation.md) |
| Chromatin: build a library of simulations (resumable) | `scripts\run.cmd python scripts\run_chromatin_library.py --config configs\chromatin\library_loop_extrusion.yaml --library loop_extrusion --n 100` (default 10 Mb region, 8000 production blocks = 61 min of polymer time; ~43 min each, ~73 h of GPU for 100; change `polymer.n_monomers` for another region size) |
| Polymer calibration (nm and seconds per polymer unit) from the saved validation simulations | `scripts\run.cmd python scripts\calibrate_chromatin.py` (writes `configs\calibration\default.json`) |
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

`scripts\dashboard.cmd` starts a local web app (this PC only, http://localhost:8501). The top of its home page links to the narrated tutorial (served by the
dashboard itself, opens in a new tab) and to the tutorial video. The sidebar's **New
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
`scripts\run.cmd python scripts\run_tracking_demo.py --config configs\loci_library.yaml` uses saved polymer simulations for the loci (needs a finished library named in `loci.library`). Its nm-per-unit and seconds-per-block come from `configs/calibration/default.json` (`loci.calibration: default`; remake it with `scripts\run.cmd python scripts\calibrate_chromatin.py`). That calibration rests on ASSUMED numbers (a 100 kb pair is within 400 nm 10% of the time; a single locus has D = 0.002 um^2/s at a 1 min lag) until you supply a measured MSD. One saved block is then 0.46 s, so the 8000-block library simulations cover ~61 min of movie (a longer movie repeats the trajectory back and forth, with a warning). `loci.positions_kb: random` gives every cell its own 100 kb promoter-enhancer pair. `loci.nm_per_unit` / `loci.block_duration_s` still override it.

Stage 7b: when the movie has tracked identity sources, `stage6_analysis/scores.json` has a `propagation` entry per source: the locus results scored through the tracker's cell IDs (detection rate, position and distance errors, share of detections on tracks that mix two true cells).


### Watching a library run, and choosing simulations (🗂 Library tab)
Start the generator in one window (`scripts\run.cmd python scripts\run_chromatin_library.py --config configs\chromatin\library_loop_extrusion.yaml --library loop_extrusion --n 100`),
the dashboard in another (`scripts\dashboard.cmd`). The **🗂 Library** tab shows live progress (refreshes every 5 s; shows "Stalled" if the PC slept). One library selector at the top serves the whole page.
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
The enhancer's effect is a smooth Hill curve, half its strength at 400 nm (the default since 2026-10-20; the old sharp 150 nm contact was reached in under 1% of frames
by the polymer-library distances, so it barely did anything).

### Sidebar checkboxes and library movies
The single-cell loci form now really passes the **probe binding** and **transcription** settings from their tabs to the run, and library movies place the
two loci at the positions set in the Transcription tab (default: 100 kb apart in the middle of the chosen library's region).


### Image errors and the 🖼 Images tab
The microscope model (stage 3) has 28 switchable sources of image error (22 original ones, plus z-slice skew, piezo, biplane registration, sequential colours, cytoplasm and stray spots, which are described in their own sections): stage drift, chromatic shift between colours, diffraction blur, focus drift, out-of-focus haze,
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

**Several promoter-enhancer copies per cell, and stray spots (2026-10-19).**
- **Simulation:** `loci.copies` gives some cells two (or more) promoter-enhancer pairs, e.g. both alleles labelled. Each copy has its own chromatin
  trajectory, probe binding and promoter, so it bursts on its own; `loci.copies.transcription` can give a copy its own rates (the loci configs: half
  the cells have two copies, copy 1 bursts twice as often). The image-error source `stray_spots` adds short-lived decoy spots (probe or MCP
  aggregates) inside the nuclei in every colour. Sidebar: *Promoter-enhancer copies and stray spots* in the loci form.
- **Analysis:** with `locus_analysis.max_copies` > 1 the multi-copy linker (`stage6_analysis/multicopy.py`) follows several paths per colour, keeps the
  ones that look like a locus, untangles swapped copies, pairs each promoter with its enhancer and uses the partner colour (`pair_prior`, default
  `gauss`) to choose spots. Every detection and MS2 trace carries its `copy`. Since 2026-10-20 it also gives every copy back the spots that lie
  where it has recently been (`untangle_by_place`): swaps made through a wrong pick or a gap are undone (picks on the other copy: at most 0.8%,
  was up to 2.5%).
- **Measured** over 25 movies (`scripts\run.cmd python scripts\sweep_copies_strays.py`, Spot finding page part 3e): copy number right in every cell
  without strays and in 90% of cells at 4 strays per nucleus per minute; one path per colour (the one-copy finder) misses the second copy.
  Details and the bias check: [docs/SPOT_FINDING.md](docs/SPOT_FINDING.md).

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
`scripts\run.cmd python scripts\deconvolve_ms2.py --run <run>`; details and measured accuracy in [docs/SPOT_FINDING.md](docs/SPOT_FINDING.md). Has the sampler converged? Four chains per cell, R-hat and effective sample size: `scripts\run.cmd python scripts\check_mcmc.py --run <run>`, results (the plain sampler does not converge on the dense gene; the default now checks every cell and uses parallel tempering where it fails) in [docs/MCMC_CONVERGENCE.md](docs/MCMC_CONVERGENCE.md).

**Promoter-enhancer distance and P(burst ON) (🔗 E–P curve page).** How well can the whole curve of P(promoter in a burst) against the E-P distance be recovered, with localization error, the 2D projection, colour shift, wrong picks, a noisy MS2 trace and a limited number of cells? `scripts\run.cmd python scripts\sweep_ep_curve.py --exps sigma shift outliers cells null --workers 12` (the page also simulates one population live); method and measured numbers in [docs/EP_CURVE.md](docs/EP_CURVE.md). Wrong picks are corrected in
the fit when their share is supplied (`Measurement.assumed_outlier_rate`; 5-40% wrong picks: contrast 4.11-4.13 against 4.07); the cell-bootstrap
intervals cover the truth 95% of the time (checked on 60 populations). Second version (2026-10-14): E–P distances from the polymer simulations (`scripts\run.cmd python scripts\sweep_ep_polymer.py --exps modality null wrong cells --workers 10`), seen through a 2D, z-stack or biplane microscope; the main finding is that the distance changes faster than the promoter responds, so the curve is fitted from the distance history.

**Instruments.** `--set instrument=vutara_vxl` (or the sidebar's *Microscope* choice) switches from the generic placeholder optics to a
Bruker Vutara VXL in 2D mode:
- **Optics:** 60x / 1.3 NA silicone objective, ~100 nm pixels, Orca Fusion BT (QE 95%), flat top-hat illumination.
- **Sequential colours:** bleed-through computed from the FPbase spectra (~0.04%), and the loci move between the colour exposures.

Sources and assumptions are in [docs/INSTRUMENTS.md](docs/INSTRUMENTS.md); the measured sequential-vs-simultaneous trade-off is in
[docs/ERROR_EFFECTS.md](docs/ERROR_EFFECTS.md).

**MS2 → Pol II page:**
- per-second probability of a Pol II start over the trace (with zoom);
- the paper's own sampler next to the default one (they agree);
- a step that explains the burst rule and maps how well each setting matches the truth.

The dashboard is now dark. [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) has updated diagrams of how it all works and of what is built.

**Built since the plan** [docs/PLAN_three_signals_and_homologs.md](docs/PLAN_three_signals_and_homologs.md): two copies per cell and the
promoter-enhancer pairing (above). **Not built yet:** using the MS2 spot as a third witness to choose spots, sister chromatids, quality flags per frame.

The dashboard now runs only the page you have open (a navigation bar replaces the tabs), which made reruns ~5-10x faster; the run summary has a one-line description
of every image-error source.
