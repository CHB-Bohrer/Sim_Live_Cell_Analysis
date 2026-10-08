# Plan: transcription (MS2) model driven by promoter-enhancer distance

Status: PLAN ONLY (2026-10-08). Nothing here is built. Decisions needing the user are marked **ASK**.

## Goal
A stochastic model of transcription whose rate depends on the 3D distance between the promoter and the enhancer, taken from the
saved polymer simulations. Its output is the MS2 signal (a bright spot at the gene) in the simulated movies, so the analysis can later be
tested against known ground truth: true P-E distance, true promoter state, true Pol II count, true MS2 intensity.

## What the polymer simulations already give us
- `tracked_positions.npy`: every 10th monomer (10 kb) at EVERY block, so the P-E distance d(t) between any two tracked beads
  (promoter and enhancer positions must be multiples of 10 kb, or snap to them; `library.locus_trajectories` already does this).
- `lef_positions.npy`: the loop-extruder legs through time, so we also know when a loop bridges P and E (a "loop state" flag).
- Full snapshots every 10 blocks (contact maps, P(s)): the equilibrium distribution of d and the contact probability at the P-E
  genomic separation, used to choose sensible distance thresholds.
- Units: d is in polymer units now; `nm_per_unit` converts to nm, `block_duration_s` to seconds (both PLACEHOLDERS until the MSD calibration).

## The model (stage 1c, `stage1_chromatin/transcription.py`)
Promoter state (telegraph model), with the distance d(t) entering one rate:

    OFF --k_on(d)--> ON --k_off--> OFF            (ON: promoter permissive; Pol II initiates at rate k_init)
    each initiation -> one Pol II travels the MS2 cassette then the gene body at speed v, then is released.

- Coupling forms (config `coupling.kind`): `none` (control), `contact` (k_on = k_basal + k_max if d < d_c, else k_basal),
  `hill` (k_on = k_basal + k_max / (1 + (d/d_half)^n); smooth), `exp` (k_basal + k_max exp(-d/lambda)).
  `coupling.rate` picks WHICH rate depends on d: `k_on` (default), `k_init`, or `k_off` (sticky bursts).
- Exact simulation: d(t) is piecewise constant per polymer block, so within a block the rates are constant and the telegraph
  chain can be simulated exactly (uniformization or per-block exponential waiting times). Pol II elongation is a deterministic delay
  after each initiation event (no extra state).
- MS2 intensity: a Pol II carries n_loops x fraction of the cassette already transcribed (a ramp while it transcribes the cassette),
  then the full signal until it leaves the gene (plateau), then drops. Intensity(t) = sum over Pol II of that profile, averaged over the
  exposure window. (Later: replace the proportionality by the probe model with capacity = loops, so coat-protein binding adds its own noise.)
- Output per cell and frame: d_nm, promoter_on, n_polII, ms2_intensity (as imaged); plus the event list (initiation times) for ground truth.
- One-way coupling only: the polymer does not feel transcription (no transcription-induced compaction), because the polymer
  simulations are pre-computed. State this limitation in the UI.

## The central technical problem: time scales
Polymer trajectories cover (2000 production blocks) x `block_duration_s`; bursts last minutes and movies are 20 min or more. After the
MSD calibration the library may cover only seconds to a few minutes. Plan:
1. First version: use d(t) directly and, if the movie is longer than the trajectory, LOOP it with a flag and a warning
   (time-reversal-symmetric stitching avoids a jump; document that correlations beyond the loop length are artificial).
2. Second version: a surrogate distance process fitted per library to the polymer data (matches the distribution of d, including the
   looped fraction, and the autocorrelation), able to run for arbitrary durations. Validate against held-out polymer simulations.
3. Longer simulations only if the surrogate is not good enough.

## Rendering (stage 3)
MS2 spot drawn at the PROMOTER position (the gene), intensity = photons_per_loop x loops transcribed, same PSF as the loci.
**ASK:** is MS2 its own (third) colour, or the same colour as the nuclear NLS-GFP (as CLAUDE.md says)? Same colour means the spots sit
on top of the nuclear channel and can disturb segmentation and tracking, which would be a real, testable error source.

## Dashboard (new tab "Transcription", same style as 🧪 Probes)
- Choose a library simulation (or a stand-in distance process if there is no library) and the promoter / enhancer positions (kb).
- Plot d(t) with the loop-state shading; plot the rate k_on(d) against d OVER the histogram of d, so it is obvious which distances matter.
- Show promoter state, Pol II count and MS2 intensity through time, with the model drawn like the probe diagram.
- Parameter panel with presets (no coupling, contact-gated, smooth Hill), tooltips, derived numbers (duty cycle, burst frequency,
  mean MS2), and a sweep (e.g. d_half, k_max, loop fraction) with exact-theory overlay where it exists.
- Hand-off to movies like the probes (sidebar checkbox, `--set` command), and `scripts/sweep_transcription.py`.

## Tests (analytic answers first)
- Constant d: ON fraction = k_on/(k_on+k_off); mean Pol II = k_init x ON fraction x gene residence time (Little's law); mean MS2 from the profile.
- Piecewise d with two values: time-average of the ON fraction matches the weighted theory.
- `coupling.kind: none`: no correlation between MS2 and d (null test); `contact`: correlation positive and lagged by the burst/elongation time.
- Reproducible from the seed; looping flag raised when the movie is longer than the trajectory.

## Build order
1. `distance_trace(sim, p_kb, e_kb)` in `library.py` (+ loop-state flag), unit-tested with the fake library.
2. `transcription.py` core + analytic tests.
3. Wire into the pipeline (`loci.transcription.enabled`, stream ss[5]) and the MS2 spot rendering.
4. Dashboard tab + sweep script.
5. Surrogate distance process (time-scale fix) once calibration exists.

## Questions for the user (**ASK**)
1. MS2: third colour, or same colour as nuclear GFP?
2. Positions: which genomic positions / separation for promoter and enhancer inside the region? (multiples of 10 kb)
3. Which rate depends on distance (k_on is the usual choice), and what shape: sharp contact threshold, or smooth?
4. Biology numbers when available: ON/OFF times, initiation rate, Pol II speed, gene and MS2-cassette length, number of loops.
5. Frame interval and movie length (decides how much the time-scale problem matters).
