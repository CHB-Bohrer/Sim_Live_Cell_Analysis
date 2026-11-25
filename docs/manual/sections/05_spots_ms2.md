# Finding loci and reading transcription {#analysis}

::: summary
**In one sentence.** In each cell and each colour, a difference-of-Gaussians filter proposes candidate spots. A Viterbi path then picks the one that moves like a diffusing locus, and a Gaussian fit gives its position. The MS2 trace is turned back into Pol II start times by Bayesian deconvolution.
:::

This chapter is a measured **stand-in** for the lab's own analysis, which is not written yet. It follows Chen, Levo, Zoller & Gregor (2025).

## Locus tracking

1. **Clean and find:** median filter, then a difference of Gaussians (1 and 3 px) with a low threshold, keeping up to 5–8 candidates per frame.
2. **Choose a path through time (Viterbi):** each step is scored by how likely the move is for a diffusing locus. The path may skip short gaps, and is rewarded for bright candidates over clutter. D is re-estimated from the path.
3. **Fit:** a Gaussian fit gives the sub-pixel position, and the brightness is the mean in an ellipse minus the surrounding ring.
4. **Several copies per cell:** paths for each colour are untangled and paired into promoter–enhancer copies, with a prior on their distance. Copy swaps through wrong picks are repaired.

![**Localization precision.** Spots drawn like the microscope draws them (pixelated Gaussian, Poisson noise, background, read noise). **A**: the error falls as 1/√photons and follows the theoretical best possible precision for this background (Mortensen et al. 2010). **B**: the finder is within 0–15% of that limit across 150–10,000 photons. Made by `figure_scripts/fig_spot_precision.py`.](figures/spot_precision.png){width=100%}

On full simulated movies (showcase run, two copies per cell, stray spots at 1 per minute, calibrated), the finder found 98.6% of locus frames with 0.6% wrong picks and a median error of 21 nm. The old brightest-pixel finder found 58% with 20% wrong picks. The linker gets the number of copies per cell right in 100% of cells at up to 0.5 stray spots per nucleus per minute, 98% at 1–2 and 90% at 4. One path per colour, without the linker, gets it right in only 58%.

## MS2 → Pol II starts

The MS2 trace is the intensity at the promoter in every frame, corrected for bleaching and uneven light. Each Pol II adds a known step-shaped pulse: loops build up over the cassette, hold, then drop at release. So the trace is a sum of pulses plus noise. We sample the posterior over start times by Gibbs sampling with shift moves. Four chains per cell are checked for convergence (R-hat < 1.05, effective sample size ≥ 400), and the sampler falls back to parallel tempering when a cell fails.

![**Deconvolution.** **A**: on a problem small enough to add up all 16,384 start patterns exactly, the samplers reproduce the exact start probabilities (largest error 0.005 for the default, 0.03 for the paper's block sampler). **B**: a noisy trace with three known starts (dotted). The probability of a start (red, per 30 s) sits on the true starts. Made by `figure_scripts/fig_deconvolution.py`.](figures/deconvolution.png){width=100%}

On simulated movies with the bursty default gene, 79–86% of true starts were found within 30 s, and the 1-minute rate correlated with the truth at r = 0.88–0.91.

## Limits

- **MS2 brightness:** at the placeholder MS2 brightness (one Pol II ≈ 0.8 nuclear pixels of light), a slow gene can't be measured. Bursty genes can. The real brightness ratio must be measured.
- **Release time:** the kernel assumes a fixed release time, but the simulated wait is random. This costs about 10 points of starts found.
- **Not built:** using promoter, enhancer and MS2 together to choose spots, and separating sister chromatids.
