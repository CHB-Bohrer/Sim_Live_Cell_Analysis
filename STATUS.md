# Status

**Project:** Sim Live Cell Analysis. Simulated live-cell movies of chromatin loci with known ground truth, used to measure how much each step of the analysis gets wrong.
**Updated:** 2026-11-27 (stage states as of the 2026-11-10 improvement pass; folder moved 2026-11-25)

This is the short front page. The full record, with every finding and gotcha, is `CLAUDE.md`.

## Where things stand
- Built: cell motion (stage 2), microscope image model with 28 switchable error sources, 3D z-stacks, biplane and calibration (stage 3), Cellpose + Trackastra tracking (stage 4), per-cell isolation (stage 5a), tracking metrics and error propagation for 2D (stage 7).
- Stand-ins: stage 6 is a measured placeholder (DoG + Viterbi spot finder, MS2 trace, Pol II deconvolution) until the owner's own analysis replaces it.
- Chromatin: polychrom engine validated (11 of 12 checks; insulation 1.36-fold against a 1.5 target, left as an honest miss). The 100-simulation loop-extrusion library is finished.
- Nearly every instrument, probe, MS2 and polymer-calibration number is a placeholder or an assumption, not a measurement.
- Tests: 225 CPU + 7 GPU at 2026-11-10. Narrated tutorial: 38 scenes, video release `tutorial-video-2026-10-25`.

## Last session
- 2026-11-27: added this file, `THOUGHTS.md`, `scripts/thoughts.py`, `scripts/next_steps.py` and the "Owner's standing preferences" section in `CLAUDE.md`. No code or results changed.

## Next
- [ask] Condense `CLAUDE.md` (about 840 lines, loaded every session) into short rules plus topic files under `docs/`. Largest usage saving available; needs care so no gotcha is lost.
- [ask] Write the paper-style design document (PDF) for this project; today the explanation is spread over `docs/*.md`.
- [ask] Start the plain-polymer library (about 73 h of GPU).
- [ask] Re-make and publish the tutorial MP4 (pictures changed slightly on 2026-11-10).

## Waiting on the owner
- Measured values for the placeholders: MS2 brightness, instrument and camera, probe kinetics, stray-spot rate, wrong-pick share.
- A measured locus MSD for the polymer calibration, and a real movie with some tracked cells.
- Which meaning of nucleus "steadiness"; whether loops should live longer; whether a surrogate distance process should cover movies over 61 min.
- `gh auth refresh -s workflow` to enable `docs/ci/syntax.yml`.

## Where to look
- `CLAUDE.md`: full handoff record. `README.md`: how to run everything.
- `docs/ARCHITECTURE.md`: pipeline diagrams. `docs/tutorial/index.html`: narrated tutorial.
- `src/simlive/stageN_*`: one subpackage per pipeline stage. `status.json`: public status page.
