# Plan: transcription (MS2) model driven by promoter-enhancer distance

Status (2026-10-09): BUILT (steps 1-4 of the build order; the surrogate distance process, step 5, is not). Code: `stage1_chromatin/transcription.py`,
`library_loci.py` (distance traces), `render.py` (MS2 spots), `app/transcription_tab.py`, tests `tests/test_transcription.py`.

Decisions made by the user on 2026-10-09:
- MS2 is a third colour but the SAME colour as the nuclear NLS-GFP (MS2-MCP carries an NLS and accumulates in the nucleus). The spot is
  rendered INTO the nuclear channel at the promoter; the user expects it to raise the background only slightly, not disturb segmentation
  (to be tested, not assumed).
- Promoter (locus 0) and enhancer (locus 1) are 100 kb apart: `loci.positions_kb: [5050, 4950]` in `configs/loci_library.yaml`.
- Coupling: all four shapes are available (contact, hill, exp, none); the default is the sharp contact. Which rate depends on distance: all
  three selectable (default k_on).
- Numbers: from the literature (Larson lab), with the caveat below. Frame interval: must be changeable at any time and the kinetics must follow
  (implemented as: every rate is per SECOND of real time and the model is continuous; the frame interval only sets WHEN it is sampled, and the
  library loci are interpolated between polymer blocks, so nothing has to be re-entered. The dashboard warns when a frame interval is too
  coarse for the chosen kinetics). If the user meant something different by "automatically change kinetic rates", revisit.

**CHECKED AGAINST THE PRIMARY PAPERS (2026-10-09, full text on PMC):**
- Rodriguez et al. 2019 Cell (TFF1-MS2, MCF7, [PMC6331006](https://pmc.ncbi.nlm.nih.gov/articles/PMC6331006)): active period (ON) 16.0 +/- 0.5 min,
  narrowly distributed; the gene bursts every 66 +/- 7 min (complete media) or 86 +/- 18 min (0.5 nM E2), 185 +/- 34 min near the EC50, so mean OFF is
  roughly 50-70 min, with a broad distribution and some alleles inactive > 12 h; burst size 1.5 +/- 0.5 transcripts; nascent-RNA dwell 13.0 +/- 0.8 min;
  initiation 0.5 +/- 0.02 per min. No elongation rate is reported there. (Read only the first 100,000 characters of the page.)
- Coulon et al. 2014 eLife (beta-globin reporter, [PMC4210818](https://pmc.ncbi.nlm.nih.gov/articles/PMC4210818)): elongation ~2.6 kb/min (Table 1: 2.60 +/- 0.16),
  dwell of transcripts at the transcription site (3' end) 116 +/- 6 s. No initiation interval reported. (First 100,000 characters read.)
- So: the ON time of ~6 min used below is the MYC-paper number and is ~2.7x shorter than TFF1's 16 min; our dwell default of 30 s is 4x shorter than
  Coulon's 116 s; 2.5 kb/min matches 2.6. Decide which gene to mimic and set ON (16 min for TFF1), OFF (~50 min), dwell (~116 s) accordingly.
  The initiation rate while ON is not a literature value here (TFF1: ~1.5 transcripts per 16 min burst means fewer than one initiation per 10 min).
  Note that makes our default "a Pol II every 30 s while ON" (about 12 per 6-min burst) far busier than TFF1.

Earlier (unverified) summary, kept for history (SEARCH-RESULT SUMMARIES ONLY; the primary papers could not be opened because PMC / eLife / Europe PMC are blocked from the
build machine; VERIFY against the papers before relying on them): MCF7 cells, MS2 live imaging: mean ON ~5.6-7.1 min, mean OFF ~44-70 min (Cell
Reports 2021, MYC paper, Larson co-author); TFF1 (Rodriguez et al. 2019 Cell, review by Rodriguez and Larson): ~1-2 transcripts per burst, burst
period ~86 min, estradiol ~2x more bursts; human Pol II elongation spans ~1-6 kb/min across studies (Darzacq 2007: ~1.9-4.3), yeast (Larson 2011
Science) ~1.2-2.8 kb/min. Defaults chosen: ON 6 min, OFF 60 min when far, 2.5 kb/min, 24 loops. The coupling strength (fold 12), initiation
interval (30 s), cassette/gene lengths and dwell time are design choices, NOT literature values. All MS2 brightness numbers are placeholders.


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

## Build order (1-4 done 2026-10-09)
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
