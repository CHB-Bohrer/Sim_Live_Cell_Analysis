# Chromatin: the polymer simulation {#chromatin}

::: summary
**In one sentence.** A 10 Mb stretch of chromatin is simulated as a bead-spring polymer in a sphere. Loop-extruding factors (cohesin-like) pull loops and stop at boundary elements (CTCF-like). We save the 3D path of every bead through time.
:::

## What it does

- **The polymer:** 10,000 beads, one per kb, joined by springs, with a little stiffness and soft repulsion. The chain sits in a sphere at 10% volume density. It is run with polychrom and OpenMM on the GPU, using the parameters of the published polychrom loop-extrusion example.
- **Loop extrusion:** a separate 1D model moves loop-extruding factors (LEFs) along the chain, one bead per step per leg. Each LEF stays bound for a set lifetime, and its legs stop at boundary elements. In 3D, each LEF becomes a short stiff bond between its two legs.
- **The library:** 100 simulations with loop extrusion have been run and saved (8,000 time blocks each, about 43 min of GPU time per simulation). Each cell in a movie gets one saved simulation, rotated randomly and started at a random time. Each cell also gets a random pair of loci 100 kb apart, the promoter and the enhancer.

## How we know it works

We ran five conditions with two seeds each and checked them against known polymer physics and published loop-extrusion results (Fudenberg et al. 2016).

![**Contact probability P(s) versus genomic distance.** *Left:* the plain polymer falls with slope −1.43 between 10 and 100 kb, between the equilibrium chain (−1.5) and the fractal globule (−1), as expected for a confined chain. *Right:* loop extrusion raises contacts about 3-fold, with a peak near 150 kb (the loop size), and lowers them beyond about 1 Mb, as in Fudenberg et al.](../img/chromatin/contact_probability.png){width=100%}

![**Locus mean squared displacement (MSD).** A single bead moves sub-diffusively: MSD ∝ lag^0.49, as expected for a Rouse polymer (0.5). The left panel is in the lab frame; the right has the chain's centre of mass removed.](../img/chromatin/msd.png){width=100%}

![**Insulation at boundary elements.** Contacts across a position dip at the boundaries (vertical lines) only when there is loop extrusion and boundaries.](../img/chromatin/insulation.png){width=80%}

| Check | Expected | Measured | |
|---|---|---|---|
| P(s) slope, plain polymer, 10–100 kb | −1.6 to −0.8 | −1.43 | pass |
| Contacts raised at 50–300 kb by loop extrusion | ratio > 1 | 2.83 | pass |
| The two legs of each loop are really held together in 3D | < 1.5 bead sizes | 0.61 | pass |
| Locus MSD exponent | 0.3–0.7 (Rouse 0.5) | 0.49 | pass |
| Designed loop anchors form a "dot" | enrichment ≥ 3 | 34.7 | pass |
| Mean loop size below LEF processivity | 20–200 kb | 101 kb | pass |
| Contacts have settled (first vs second half) | log ratio < 0.1 | 0.077 | pass |
| Insulation across impermeable boundaries | fold ≥ 1.5 | 1.36 | **check** |

Insulation is real but weaker than the target we set in advance. It is 1.36-fold, against 1.05 without boundaries. The roughly 2-fold figure in the paper measures a different quantity, so the two numbers can't be compared directly. More seeds per condition would settle it. We did not move the target so that the check would pass.

::: detail
**A bug this validation caught.** The first run looked plausible, but its loops were frozen in place. OpenMM cannot change which two beads an existing bond connects. Now every possible bond is registered at the start and switched on or off by its stiffness. A GPU test checks that the legs of every loop stay about one bond length apart.
:::

## From polymer units to nanometres and seconds

The simulation has no physical units. We set them by two rules. Both are **assumptions** that wait for a measured MSD:

1. A pair of loci 100 kb apart is within 400 nm of each other 10% of the time. This gives **55.4 nm per bead size**.
2. A single locus has D = 0.002 µm²/s at a 1 min lag. This gives **0.46 s per simulation block**, so a library simulation covers about 61 min of real time.

## Limits

- **Units:** the length and time units are assumed, not measured.
- **Boundaries:** boundary positions are random (synthetic TADs), not real CTCF sites.
- **Loop lifetime:** at 0.46 s per block, loops live only about 46 s, while real cohesin loops last minutes.
- **Long movies:** movies longer than about 61 min replay the trajectory backwards and forwards.
- **Not run yet:** the plain-polymer library, without loop extrusion.
