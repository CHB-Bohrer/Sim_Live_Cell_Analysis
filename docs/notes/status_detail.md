# Detailed build status per stage; tracking findings

Working notes moved verbatim from `CLAUDE.md` on 2026-11-27 (section order kept). Add new dated sections at the top; delete statements that are no longer true.

## Status (update this section every session)
Built: stage 2 (cell motion, collisions, division, irregular time-varying nuclear shapes, `speed_scale`), stage 3
(nuclear channel + one colour channel per locus: per-nucleus texture + nucleoli and 22 switchable image-error sources, see "Image errors"), stage 4
(Cellpose default, threshold-watershed baseline, Trackastra `general_2d` behind the `Tracker` adapter), stage 5a
(per-cell isolation: `stage5_linking/isolate.py`, `celldata.py`), stage 7a (CTC metrics via traccuracy + own
identity-switch analysis), `scripts/run_tracking_demo.py`, `scripts/sweep.py`, dashboard (tabs incl. Scans, Cells).
STAND-INS (to be replaced): stage 1 = `stage1_chromatin/toy_loci.py` (2 loci per nucleus doing confined diffusion in
the nucleus frame, NOT polychrom); stage 6 = `stage6_analysis/spots.py` + `localize.py` / `localize3d.py` (placeholder for the user's
analysis, but now a measured one: see "Spot finding" below). The loci/isolation pipeline runs via `configs/loci_demo.yaml`.
Library wiring (2026-10-08): `stage1_chromatin/library_loci.py` (`loci.source: library`, config `configs/loci_library.yaml`; one simulation per cell, random rotation + time offset, polymer coordinates converted to true nm with `nm_per_unit` (from `loci.calibration`, see "Polymer calibration"; NOT rescaled to the nucleus) and placed at a random territory centre per cell (re-drawn until all loci stay inside the nucleus in every frame), `block_duration_s` likewise) is unit-tested with a fake library on CPU and was run end-to-end (tracked, 2026-10-15) on the 4-simulation `smoke` library; The 100-simulation loop-extrusion library exists since 2026-10-18 (showcase loci_library_seed12_87a16d uses it).
Dashboard redesign + library monitor (2026-10-08; the monitor was checked LIVE on 2026-10-20 with a 2-simulation test library on the GPU: running badge,
1 of 2 finished, current simulation's block and time left, GPU load, then "Stalled" after the generator was stopped): `.streamlit/config.toml` + `app/style.py` (one theme, CSS, matplotlib defaults), new 🗂 Library tab
(`app/library_tab.py`): live monitor (refreshes every 5 s) of `<library>/status.json` (written by `run_chromatin_library.py`) and the
running simulation's `progress.json` (written by `polymer3d.py` every 20 blocks): state badge (running / stalled / finished /
stopped), n done, ETA, GPU load via nvidia-smi, recent log; plus a picker to choose which simulations a movie uses. ONE library selector for the whole
page (key `lib_sel`, opens on a running library, else the one with most simulations; the Chromatin page opens on the same library) since 2026-10-20.
Choosing simulations: `library.assign_simulations(..., pool=, pinned=)`; config `loci.sim_ids` (restrict the pool) and
`loci.pin: {cell_id: sim_id}`; the sidebar's loci form has a source radio (stand-in / library) that uses the Library tab's choice.
The picker's preview uses the movie's own random stream (`SeedSequence(seed).spawn(4)[2]` = `spawn(10)[2]`: a child's numbers do not depend on how many
siblings are spawned), so it matches the real assignment when no cell divides.
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
hand-off to a movie (sidebar checkbox in the loci form).
Transcription / MS2 model (2026-10-09; stage 1c): `stage1_chromatin/transcription.py` (details and the user's decisions: `docs/PLAN_transcription.md`).
Promoter OFF -(k_on)-> ON -(k_off)-> OFF, Pol II initiation k_init while ON, Pol II ramp over the MS2 cassette then plateau until released
(random dwell); MS2 signal = loops carried by all Pol II. The promoter-enhancer distance d(t) (3D, from the library trajectories at full block
resolution, via `simulate_library_loci(..., return_traces=True)`; toy loci use their 2D lab distance) multiplies ONE chosen rate
(`coupled_rate`: k_on default | k_init | k_off) by 1 + (fold-1) f(d), f = hill (default since 2026-10-20: half effect at 400 nm, n = 4, the user's E-P choice;
the old default, a sharp 150 nm contact, was reached in only 0.75% of frames by the library distances, so it barely coupled) | contact | exp | none (control). ALL rates are per second
of real time, so the frame interval can change freely (library loci are interpolated between blocks; `check_sampling` warns when too coarse).
Exact Gillespie + exact expectations (`expected_on`, `expected_ms2`); 12 tests in `tests/test_transcription.py` (constant-distance theory, two-valued
distance vs exact expectation, frame-interval independence, null/positive correlation with contact, interpolation, MS2 spot in the nuclear
channel). `loci.transcription.enabled: true` (configs loci_demo/loci_library) writes `stage1_chromatin/transcription.csv` + `transcription_events.csv`
and renders the spot at the promoter into the nuclear channel (`render_nuclei(..., spots=)`; pipeline now runs stage 1 BEFORE rendering nuclei; seed
stream ss[5]). If the movie is longer than the saved polymer trajectory, the trajectory is repeated back and forth with a warning (correlations beyond
its length are artificial); the planned surrogate distance process (PLAN step 5) is NOT built. Dashboard ⚡ Transcription tab (`app/transcription_tab.py`): nine
numbered steps (sampling, distance, coupling, promoter, Pol II kymograph, MS2 signal with exact expectation, simulated microscope tiles, parameter sweep,
hand-off to a movie + sidebar checkbox). Its widgets start from `transcription.DEFAULTS` (until 2026-10-20 the hand-off fell back to the OLD TFF1-like gene
when the tab had not been opened: a sidebar movie then used a different gene than the tab showed); step 5 opens on the cell's first Pol II. LITERATURE NUMBERS ARE UNVERIFIED search-summary values (primary papers blocked from the build machine); all
MS2 brightness, coupling strength, initiation rate and nm_per_unit / block_duration_s are placeholders.
NOT built: MSD calibration from MEASURED data (an assumed-D calibration exists, see below), surrogate distance process for long movies, 3D tracking / 3D cell isolation, stage 5 beyond
isolation (stitching the locus results of broken / switched tracks into one cell), the user's stage 6 analysis, 7b for 3D (7b exists for the 2D locus results), Ultrack/TrackMate adapters,
density and frame-interval scans.

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
