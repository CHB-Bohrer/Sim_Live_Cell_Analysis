# What each image error does to the results, and what each correction gives back

Two sets of measurements:
- **Bruker Vutara VXL** (2026-10-13, below): sequential colours vs simultaneous two-camera detection.
- **Generic placeholder microscope** (2026-10-12, further down).

## Bruker Vutara VXL: sequential colours vs simultaneous detection (2026-10-13)

The instrument and its sources are in [INSTRUMENTS.md](INSTRUMENTS.md). Two versions of the same movie:
- **Common settings:** seed 9, 12 cells, 30 min, the denser gene, MS2 R = 25, cytoplasm on.
- **Sequential** (`loci_demo_seed9_1c649d`, as the Vutara does it): bleed-through of GFP into the promoter's colour 0.036% (from the
  FPbase spectra); the colours are 0.15 s apart.
- **Simultaneous** (`loci_demo_seed9_17363e`, the optional second camera): bleed-through 8.6%, no delay.

| movie | analysis | wrong pick | median nm | p95 nm | distance rms nm | intensity r | intensity err % | late/early | MS2 r |
|---|---|---|---|---|---|---|---|---|---|
| error-free reference | | 0.0% | 10 | 23 | 13 | 0.96 | 7.2 | 1.00 | 1.00 |
| **sequential: colour delay 0.15 s** | nothing can correct it | 0.0% | 17 | 72 | **34** | 0.96 | 7.2 | 1.00 | 1.00 |
| sequential: bleed-through 0.036% | not corrected | 0.0% | 10 | 23 | 13 | 0.95 | 7.8 | 1.00 | 0.99 |
| | unmixing from the movie | 0.0% | 10 | 23 | 13 | 0.96 | 7.3 | 1.00 | 1.00 |
| **simultaneous: bleed-through 8.6%** | not corrected | **9.9%** | 7 | 23 | 11 | **0.47** | 40.8 | 0.90 | 0.98 |
| | shell background | 9.9% | 7 | 23 | 11 | 0.46 | 47.1 | 0.89 | 1.00 |
| | unmixing from the movie | 0.0% | 15 | 69 | **32** | 0.83 | 15.2 | 0.99 | 0.99 |
| cytoplasm | not corrected / shell | 0.0% | 10 | 23 | 13 | 0.96 | 7.2 | 1.00 | 0.99 / 1.00 |
| | movie map: median (naive) | 0.0% | 10 | 23 | 13 | 0.96 | 7.8 | 1.00 | 1.00 |
| | movie map: low percentile (robust) | 0.0% | 10 | 23 | 13 | 0.96 | 7.6 | 1.00 | 1.00 |
| colour shift | not corrected | 0.0% | 34 | 50 | 48 | 0.96 | 7.2 | 1.00 | 1.00 |
| | bead slide | 0.0% | 10 | 23 | 13 | 0.96 | 7.2 | 1.00 | 1.00 |
| | movie alone (12 cells: refused) | 0.0% | 34 | 50 | 48 | 0.96 | 7.2 | 1.00 | 1.00 |
| **all sources, sequential** | not corrected | 0.2% | 38 | 96 | 60 | 0.61 | 20.3 | 0.90 | 0.99 |
| | calibration images | 0.0% | 25 | 80 | 39 | 0.88 | 13.1 | 1.00 | 1.00 |
| | images + movie | 0.0% | 25 | 80 | 39 | 0.88 | 13.1 | 1.00 | 1.00 |
| **all sources, simultaneous** | not corrected | 5.1% | 28 | 64 | 40 | 0.30 | 48.2 | 0.90 | 0.98 |
| | calibration images | 10.1% | 25 | 41 | 27 | 0.42 | 49.3 | 0.95 | 1.00 |
| | images + movie (unmixing) | 0.0% | 23 | 66 | 34 | 0.63 | 27.6 | 1.00 | 1.00 |

What this says (one seed; the colour delay and the dyes are ASSUMED):
1. **The Vutara's sequential excitation removes the bleed-through problem.** At 0.036% it is invisible: intensity r 0.95 vs 0.96. The
   2% of the old placeholder was 55x too pessimistic for this instrument.
2. **Its price is time between colours.** In 0.15 s the loci move, and the promoter-enhancer distance error goes from 13 to 34 nm rms
   (p95 72 nm). Nothing can correct it afterwards.
   - It scales with the square root of the delay and with chromatin mobility. These stand-in loci move fast: 0.6 µm spread, 120 s
     relaxation. Real loci, with D ~ 0.002-0.01 µm²/s, would give roughly 2D·δ ≈ 25-55 nm.
   - **Measure your real delay and MSD.** Keep the two locus colours adjacent in the order (the preset does: promoter → enhancer → GFP).
3. **Simultaneous two-camera detection is not better for distances.**
   - Without unmixing, the 8.6% GFP ghost of the bright nuclei and MS2 spots makes the promoter finder pick the wrong object in 9.9% of
     detections and loses 12% of loci; intensities are useless (r 0.47).
   - Unmixing from the movie fixes detection and intensities (r 0.83) but leaves the leaked photons' noise: distance error 32 nm, about
     the same as the sequential delay.
4. **Cytoplasm:**
   - **The bias is real:** the naive median background map sits +9.5 counts too high in the nuclear channel (median over the field; worst
     1%: +44). The robust low-percentile map sits +1.8 (worst 1%: +36, where a cell body lingered so long that the floor was never seen).
   - **The spots barely notice,** because the shell background around each spot already removes the local cytoplasm. The map matters for
     the flat-field and hot-pixel steps, not for spot intensities.
5. **Colour shift:** the movie-only estimate is now refused with fewer than 30 cells. With 12 cells it had been applied, 2.4 SE off,
   and tripled the position error. The bead slide remains the fix.
6. **All together on the Vutara (sequential):** calibration images bring positions to 25 nm and distances to 39 nm, intensities to
   r 0.88. What is left is mostly the colour delay and the defocus.

## Generic placeholder microscope (2026-10-12)

Measured 2026-10-12 with `scripts/error_effects.py` (code: `src/simlive/stage7_validation/error_effects.py`; dashboard: 🖼 Images page,
part 2, pick a source).

**How it is measured.** The run's movie is rendered again from the same ground truth and the same noise streams:
- **reference**: a perfect microscope that still has diffraction blur, quantum efficiency, photon shot noise and read noise (nothing to
  correct);
- **+ one source**: the reference with one error source added, for every 2D source the run has;
- **all sources**: every source at once, the realistic case.

Every movie is cut into single-cell crops (true cell IDs) and analysed with the correction for that source off and on. The corrections are:
- the analysis' own countermeasure: median filter, shell background, Viterbi path;
- calibration images: dark frames, flat-field slide, bead slide;
- the movie alone (`stage6_analysis/selfcal.py`, `flatfield.py`, the row and bleaching corrections);
- for the realistic case, both together.

**What is scored** (all against the truth):
- wrong-object picks (> 0.5 µm off);
- median / 95th-percentile position error;
- rms error of the promoter-enhancer distance;
- spot intensity against the true photons: Pearson r, rms error after one global scale, and late/early drift (1 = none);
- MS2 trace against the true signal (r).

**This run:** `loci_demo_seed8_f5c88d`.
- 12 cells, 180 frames of 10 s, toy loci.
- Probes and transcription on; the gene starts 3 Pol II per minute while ON.
- MS2 brightness R = 25 (30x the placeholder).
- Cells move enough to map the illumination.

**Caveats:**
- **One seed.** Differences of a few % are within the noise of a single movie.
- **Placeholder error sizes.** They change the ranking.
- **Bright locus spots** (600 photons per exposure x probes): hot pixels and cosmic rays barely matter here. With the old placeholder
  finder or dim spots they did (see CLAUDE.md "Image errors").

| movie | analysis | wrong pick | median nm | p95 nm | distance rms nm | intensity r | intensity err % | late/early | MS2 r |
|---|---|---|---|---|---|---|---|---|---|
| error-free reference | reference | 0.0% | 12 | 28 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
| Stage drift and vibration | not corrected | 0.0% | 12 | 28 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
| Chromatic shift between colour channels | not corrected | 0.0% | 43 | 64 | 63 | 0.95 | 8.5 | 1.00 | 1.00 |
|  | movie alone | 0.0% | 43 | 64 | 63 | 0.95 | 8.5 | 1.00 | 1.00 |
|  | calibration images | 0.0% | 13 | 28 | 16 | 0.95 | 8.5 | 1.00 | 1.00 |
| Focus drift (defocus) | not corrected | 0.0% | 15 | 38 | 21 | 0.95 | 8.8 | 1.00 | 1.00 |
| Locus motion during the exposure | not corrected | 0.0% | 13 | 29 | 16 | 0.95 | 8.4 | 1.00 | 1.00 |
| Out-of-focus haze (widefield) | not corrected | 0.0% | 13 | 29 | 16 | 0.95 | 8.6 | 1.00 | 0.98 |
|  | analysis' own step | 0.0% | 13 | 29 | 16 | 0.95 | 8.6 | 1.00 | 1.00 |
| Photobleaching | not corrected | 0.0% | 14 | 31 | 17 | 0.90 | 12.8 | 0.80 | 1.00 |
|  | one curve | 0.0% | 14 | 31 | 17 | 0.94 | 9.2 | 0.98 | 1.00 |
|  | movie alone (per cell) | 0.0% | 14 | 31 | 17 | 0.94 | 9.2 | 0.99 | 1.00 |
| Uneven illumination (flat-field) | not corrected | 0.0% | 13 | 28 | 16 | 0.95 | 8.8 | 0.99 | 1.00 |
|  | movie alone | 0.0% | 13 | 28 | 16 | 0.95 | 8.8 | 0.99 | 1.00 |
|  | calibration images | 0.0% | 13 | 28 | 16 | 0.95 | 8.4 | 1.00 | 1.00 |
| Lamp / laser flicker | not corrected | 0.0% | 13 | 28 | 16 | 0.95 | 8.4 | 1.00 | 1.00 |
| Constant background light | not corrected | 0.0% | 13 | 31 | 18 | 0.93 | 11.5 | 0.99 | 0.99 |
|  | analysis' own step | 0.0% | 13 | 31 | 18 | 0.93 | 9.9 | 0.99 | 1.00 |
| Autofluorescence | not corrected | 0.0% | 14 | 32 | 17 | 0.77 | 16.5 | 0.97 | 0.99 |
|  | analysis' own step | 0.0% | 14 | 32 | 17 | 0.91 | 11.5 | 0.99 | 1.00 |
|  | movie alone | 0.0% | 14 | 32 | 17 | 0.93 | 10.3 | 0.99 | 1.00 |
| Colour cross-talk (bleed-through) | not corrected | 0.0% | 10 | 27 | 15 | 0.52 | 41.4 | 1.08 | 0.99 |
|  | analysis' own step | 0.0% | 10 | 27 | 15 | 0.52 | 47.9 | 1.11 | 1.00 |
|  | movie alone (unmixing) | 0.0% | 15 | 47 | 24 | 0.91 | 11.2 | 1.00 | 1.00 |
| Dark current | not corrected | 0.0% | 12 | 28 | 16 | 0.95 | 8.2 | 1.00 | 0.99 |
|  | analysis' own step | 0.0% | 12 | 28 | 16 | 0.95 | 8.2 | 1.00 | 1.00 |
| Hot pixels | not corrected | 0.0% | 13 | 29 | 16 | 0.95 | 8.4 | 1.00 | 1.00 |
|  | analysis' own step | 0.0% | 12 | 29 | 16 | 0.95 | 8.4 | 1.00 | 1.00 |
|  | movie alone | 0.0% | 13 | 29 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
|  | calibration images | 0.0% | 13 | 29 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
| Cosmic rays | not corrected (brightest per frame) | 0.1% | 11 | 26 | 14 | 0.95 | 8.1 | 1.00 | 1.00 |
|  | path, no median | 0.0% | 13 | 29 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
|  | analysis' own step | 0.0% | 12 | 28 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
| Pixel gain differences (PRNU) | not corrected | 0.0% | 12 | 28 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
|  | calibration images | 0.0% | 12 | 28 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
| Row (banding) noise | not corrected | 0.0% | 13 | 29 | 16 | 0.94 | 9.0 | 1.00 | 1.00 |
|  | movie alone | 0.0% | 12 | 29 | 16 | 0.95 | 8.7 | 1.00 | 1.00 |
| Offset pattern (DSNU) | not corrected | 0.0% | 13 | 29 | 16 | 0.95 | 8.6 | 1.00 | 0.99 |
|  | analysis' own step | 0.0% | 13 | 29 | 16 | 0.95 | 8.6 | 1.00 | 1.00 |
|  | movie alone | 0.0% | 12 | 29 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
|  | calibration images | 0.0% | 12 | 28 | 16 | 0.95 | 8.6 | 1.00 | 1.00 |
| Digitisation and saturation | not corrected | 0.0% | 12 | 28 | 16 | 0.95 | 8.3 | 1.00 | 1.00 |
| **all sources** | not corrected | 0.1% | 40 | 88 | 59 | 0.38 | 38.3 | 0.96 | 0.98 |
|  | analysis' own steps only | 0.0% | 40 | 87 | 58 | 0.44 | 53.6 | 0.95 | 1.00 |
|  | + movie alone | 0.1% | 49 | 100 | 73 | 0.80 | 18.0 | 1.00 | 0.99 |
|  | + calibration images | 0.0% | 28 | 63 | 35 | 0.47 | 52.0 | 1.09 | 1.00 |
|  | + images + movie | 0.1% | 26 | 85 | 42 | 0.80 | 18.5 | 1.00 | 0.99 |

## What this says

1. **Positions and distances: the colour shift dominates.**
   - Uncorrected, it moves every promoter-enhancer distance by ~60 nm rms; the median position error goes from 12 to 43 nm.
   - A bead slide removes it completely (13 / 16 nm).
   - The movie-only estimate (mean enhancer-promoter vector) was correctly NOT applied: with 12 cells its standard error was 0.41 px,
     against a true shift of (-0.55, +0.35) px. The SE falls as 1 / sqrt(number of cells): just reaching 2 SE needs ~2.2x more cells
     (~27); a precise estimate (SE 0.1 px) needs ~200. **Record a bead slide.**
2. **Spot intensities: colour cross-talk dominates, once the MS2 spot is bright.**
   - 2% of the MS2 light leaks into the promoter's colour, exactly on the promoter. During a burst it is several times the promoter's own
     light, so intensity r drops to 0.52.
   - The shell background cannot remove it.
   - Unmixing from the movie (leak fraction measured at 1.997% vs 2.000% true, then subtracted) restores r 0.91 and removes the false
     drift (1.08 → 1.00).
   - Its price is position precision (median 10 → 15 nm, distance 15 → 24 nm): the leaked photons' noise stays. The fix for positions
     is optical: filters with less bleed-through, or sequential excitation.
   - (Before unmixing the leak even helps positions here, because the simulation puts the MS2 spot exactly on the promoter label.)
3. **Bleaching** is the third intensity error (late/early 0.80). One curve fixes most of it (0.98); per-cell bleaching from the movie
   0.99.
4. **Background-type errors** (constant background, autofluorescence, offset pattern, haze, dark current) need the shell background or
   the movie's background map. Without a local background, autofluorescence alone drops intensity r to 0.77.
5. **Nothing to do / negligible here:** stage drift (the per-cell crops follow the cell), digitisation, flicker (1%), pixel gains (1%),
   row noise, hot pixels, cosmic rays (bright spots, median + path).
   - Uneven illumination matters little for these intensity scores (15% vignetting): each cell stays in a small part of the field, so
     it is mostly a per-cell scale.
   - Its movie-only map was rejected by the held-out check in the single-source movie (gain 1.9% < 2%) but used in the realistic one
     (gain 3.3%).
6. **All together:**
   - Positions: calibration images help (40 → 28 nm, distance 59 → 35 nm). The movie alone cannot replace the bead slide for the
     colour shift.
   - Intensities: only the movie-only steps (unmixing, per-cell bleaching) fix them (r 0.44 → 0.80, error 54 → 18%).
   - **The best is both** (positions 26 nm, intensity r 0.80): images for what images give, the movie for what only the movie gives.
   - What is left against the error-free reference (12 nm, r 0.95) comes from defocus, the leaked photons' noise and the interplay of
     sources.
