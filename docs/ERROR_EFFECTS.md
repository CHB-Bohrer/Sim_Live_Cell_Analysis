# What each image error does to the results, and what each correction gives back

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
