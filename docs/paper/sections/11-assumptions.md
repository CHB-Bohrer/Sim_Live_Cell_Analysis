<!--
Section: Assumptions, placeholders and limits
Update when: the owner supplies a measured value, or a planned item is built.
Code: configs
Check: none
Sources: STATUS.md, docs/notes/status_detail.md, docs/notes/dashboard.md "Needs the user".
-->

# Assumptions, placeholders and limits

::: {.box}
**Short version.** The machinery is tested. Most of the numbers fed into it are not measured. Until they are, results describe how an analysis behaves under plausible conditions, not what a specific experiment will show.
:::

## Values waiting for measurements

| Quantity | Current basis | What is needed |
| --- | --- | --- |
| Nanometres per polymer unit, seconds per block | two assumed rules (Section 3) | a measured locus MSD |
| MS2 brightness relative to the nucleus ($R$) | placeholder 0.8 | a measurement on real data |
| Camera and optics values | presets, partly sourced, partly assumed | dark, flat and bead images |
| Probe kinetics and photon yield | placeholders | measured kinetics |
| Stray-spot rate and brightness | placeholders | counts from real movies |
| Share of wrong picks, for the E-P fit | must be supplied | an estimate from real data |
| Cell speed, density, frame interval, division | deliberately hard defaults | a real movie with tracked cells |
| Boundary element positions | synthetic | real CTCF positions |
| Transcription rates and coupling strength | chosen to be bursty | the owner's choice of gene |

## Known limits of the model

- Movies longer than 61 minutes replay the polymer trajectory. A surrogate distance process for long movies is planned, not built.
- Loops live about 46 s at the assumed time scale, shorter than real cohesin.
- Tracking and cell isolation are 2D, on projections, even for z-stack runs.
- The acquisition model behind the E-P imaging-mode comparison works on vectors, not on rendered z-stack or biplane movies.
- Not modelled: single-probe bleaching, EMCCD gain, rolling shutter, lens distortion, field-dependent chromatic shift, sister chromatids.

## Planned, not built

Stitching broken or switched tracks, 3D tracking, other trackers (Ultrack, TrackMate), choosing spots jointly from promoter, enhancer and MS2, the plain-polymer library, and the owner's own stage 6 analysis.
