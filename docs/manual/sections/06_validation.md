# From image errors to biological answers {#validation}

::: summary
**In one sentence.** Because the truth is saved, we can follow each error from the image to the final result. That covers tracking mistakes, colour shift, cross-talk and wrong spot picks, all the way to the curve linking E–P distance to transcription.
:::

## Error budget, one source at a time

For a chosen movie, `error_effects` redraws it with a reference microscope and adds **one** error source at a time. It analyses each version with and without the matching correction. On the showcase run:

- **Colour shift** dominates the positions: the error rises from 12 to 43 nm, and the distance error from 16 to 63 nm. Bead calibration removes it.
- **Cross-talk** of the bright MS2 spot into the promoter colour dominates the intensities: the correlation with the truth falls from 0.95 to 0.52. Unmixing restores it to 0.91.
- **Bleaching** biases late versus early brightness to 0.80. The correction restores 0.99.

## Tracking errors into locus results

For tracked IDs, every locus detection is mapped back to the true cell through the masks and scored the same way. With perfect tracking, the result equals the true-ID result exactly, which is our check of the bookkeeping.

## The E–P distance → P(burst ON) curve

The main biological readout is the probability that the promoter is bursting, as a function of E–P distance. We fit the whole curve with an errors-in-variables model that integrates over the unknown true distance, and bootstrap over cells. The effect size is the *contrast* P(ON at 150 nm) / P(ON at 1000 nm).

![**E–P curve with polymer distances and three microscopes.** **A**: the polymer distance changes within about 9 s. **B**: the promoter integrates over about 45 s, so P(ON) against the *current* distance is nearly flat (contrast 1.4); against the distance *held* over the past 45 s it is 4.2. **C–D**: from 600 cells × 1 h, widefield, z-stack and biplane imaging all recover the reference contrast (4.19) within their scatter. Biplane is the most precise.](../img/ep_polymer_summary.png){width=100%}

**What we learned.** With placeholder errors:

- **Wrong spot picks are the most harmful error:** 10% of them cut the contrast by about 13%. The fit corrects this if the wrong-pick share is supplied.
- **Projecting the 3D distance onto 2D** costs about 35% of the effect if ignored, but the fit corrects it fully.
- **Frames must come at least every 10 s.**
- **Cell numbers:** about 10× more cells are needed than with a slowly changing distance.
- **False positives:** with no real coupling, 3–7% of populations show a false effect, against a nominal 5%.

## Limits

- **Placeholder values:** all error sizes are placeholders.
- **3D:** the E–P model for 3D is at the level of position vectors (no full 3D movies), and propagation of tracking errors is 2D only.
