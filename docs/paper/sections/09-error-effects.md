<!--
Section: What each error costs
Update when: the per-source error study, the movie-only corrections or the identity propagation change.
Code: src/simlive/stage7_validation/error_effects.py, src/simlive/stage7_validation/propagation.py, src/simlive/stage6_analysis/selfcal.py
Check: none
Sources: docs/ERROR_EFFECTS.md, docs/notes/spot_finding.md, docs/notes/chromatin.md "Stage 7b built".
-->

# What each error costs

::: {.box}
**Short version.** The same movie is rendered again with one error source at a time, analysed with and without that source's correction, and scored. For positions the colour shift dominates. For brightness it is the bright MS2 spot leaking into the promoter colour. Tracking mistakes are carried into locus results by mapping each detection back to its true cell.
:::

## One source at a time

For a given run the study renders a reference microscope (blur, quantum efficiency, shot noise and read noise only), the reference plus one source, and all sources. Each version is analysed with the correction for that source off and on. A correction can come from the analysis itself, from calibration images, or from the movie alone.

Recorded results (one showcase run, one seed, placeholder error sizes; not to be over-read):

| Situation | Position error | Brightness agreement with truth ($r$) |
| --- | --- | --- |
| Reference microscope | 12 nm | 0.95 |
| All sources, no correction beyond images | 28 nm | 0.47 |
| All sources, corrections from the movie alone | 49 nm | 0.80 |
| All sources, both | 26 nm | 0.80 |

- **Colour shift** dominates position error (12 to 43 nm; distance error 16 to 63 nm). Beads fix it. An estimate from the movie alone is not significant with 12 cells and would need about 30.
- **Cross-talk** of the MS2 spot into the promoter colour dominates brightness error ($r$ from 0.95 to 0.52). Unmixing from the movie restores $r$ to 0.91 but makes positions noisier, because the leaked photons' noise remains.
- **Bleaching**: the late-to-early brightness ratio goes from 0.80 to 0.99 after correction.

## Corrections from the movie alone

`selfcal.py` builds an ordinary calibration object from the movie: a background map from frames without a nucleus over each pixel, hot pixels by a per-pixel test, the flat field from the nuclei, and the colour shift from the mean enhancer-promoter vector (applied only when it exceeds two standard errors).

One lesson is recorded in the code: the first hot-pixel test flagged 11,000 pixels, because the robust spread of integer camera values is zero. The current test finds 100 of 117 with none false.

## Carrying tracking errors into locus results

For identities from the tracker, each 2D locus detection is mapped back to the true cell it overlaps most, then scored as usual, together with the share of detections on "mixed" tracks (a tracker label that stood for two true cells). With perfect tracking this equals the truth-based score exactly, which is the built-in check.

Not built: the same for 3D, stitching broken or switched tracks back into one cell, and a study on a movie with real identity switches.
