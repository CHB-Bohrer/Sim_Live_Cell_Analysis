<!--
Section: Chromatin (stage 1)
Update when: the polymer engine, loop extrusion, the library layout or the nm / second calibration change.
Code: src/simlive/stage1_chromatin/lef1d.py, src/simlive/stage1_chromatin/polymer3d.py, src/simlive/stage1_chromatin/analysis.py, src/simlive/stage1_chromatin/library.py, src/simlive/stage1_chromatin/library_loci.py, src/simlive/stage1_chromatin/toy_loci.py, src/simlive/calibration, configs/chromatin, configs/calibration
Check: none (GPU validation: scripts/validate_chromatin.py)
Sources: docs/chromatin_validation.md, docs/notes/chromatin.md.
-->

# Chromatin: where the loci are

::: {.box}
**Short version.** Chromatin is simulated once, separately, as a library of 100 polymer simulations with loop extrusion. A movie gives each cell one of them, picks a 100 kb pair of monomers as promoter and enhancer, and converts polymer units to nanometres and seconds. That conversion currently rests on two assumptions, not on a measurement.
:::

## The model

A 10 Mb region is a chain of 10,000 monomers, 1 kb each, confined in a sphere and integrated with polychrom on OpenMM (GPU). Loop-extruding factors (LEFs) are simulated first in one dimension (`lef1d.py`): each leg moves one monomer per step, a factor lives for a set number of steps (so its processivity is twice its lifetime), and legs stall at boundary elements with a set probability. The 3D simulation then adds a bond between the two legs of every loop at every block.

Boundary positions are synthetic: random TAD sizes drawn from a seed. **[Assumption]** Real CTCF positions have not been supplied.

## The library

One library is 100 simulations of 1,000 equilibration blocks plus 8,000 production blocks, 750 integration steps each. The loop-extrusion library is finished (100 of 100, about 44 minutes per simulation on an RTX 3090). A plain-polymer library with the same settings has not been run.

For a movie, `library.assign_simulations` gives each cell one simulation. `library_loci.py` then:

- draws a random 100 kb pair of monomers for that cell, at least 1 Mb from the chain ends;
- applies a random rotation and a random time offset;
- converts coordinates to nanometres and places the pair at a random territory centre inside the nucleus;
- lets the territory centre follow the nucleus as it deforms, while the polymer itself keeps its true size and only turns with the nucleus.

So 100 simulations serve hundreds of cells, and the distances that are drawn are the same distances that drive transcription.

## From polymer units to nanometres and seconds

`calibration/polymer.py` produces one committed file, `configs/calibration/default.json`: **55.42 nm per polymer unit** and **0.461 s per saved block**.

Both inputs are assumptions supplied by the project owner, not measurements **[Assumption]**:

- a 100 kb pair is within 400 nm of each other 10% of the time (sets the length);
- a single locus has $D = 0.002\ \mu\text{m}^2/\text{s}$ at a one-minute lag (sets the time, through the mean-squared displacement).

One consequence follows directly. At 0.461 s per block, 8,000 blocks cover 61 minutes. A longer movie replays the trajectory back and forth, with a warning, and correlations beyond 61 minutes are then artificial. Loops also live only about 46 s at this time scale, shorter than real cohesin.

## Does it work?

The engine was validated by reproducing known polymer behaviour (5 conditions, 2 seeds each, at 10 Mb). Recorded results, from `docs/chromatin_validation.md`:

| Quantity | Result | Reference |
| --- | --- | --- |
| Contact probability slope, plain polymer | $-1.43$ | near $-1.5$ for an ideal chain |
| Locus MSD exponent | 0.49 | 0.5 (Rouse) |
| Distance between the two legs of a loop | 0.61 | about the bond length (0.5) |
| Mean loop size | 101 kb | set by processivity and spacing |
| Designed loop anchors in contact | about 35 times a typical pair | far above the control (0.9) |
| Insulation across impermeable boundaries | 1.36-fold | target of at least 1.5, set in advance: **not met** |

Eleven of twelve checks pass. The insulation check does not, and it is reported as a miss. The effect is real (controls give 0.82 to 1.11) but weaker than the target, and the threshold was not lowered to make it pass. More seeds would tighten the number.

::: {.box}
**A bug worth remembering.** The first validation looked plausible and was wrong. OpenMM cannot change which particles an existing bond connects, so the loop bonds stayed frozen where they started and the two legs sat about 12 units apart instead of 0.5. The summary statistics still looked reasonable. The fix registers every bond in advance and switches bonds on and off by stiffness. The lasting rule: check a coupling physically (here, the leg-to-leg distance), not only through downstream statistics.
:::

## The stand-in

`toy_loci.py` gives each nucleus two loci doing confined diffusion. Its numbers are invented and are never to be presented as physical results. It exists so the rest of the pipeline can run without the GPU library.
