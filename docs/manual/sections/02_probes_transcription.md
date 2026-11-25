# What makes loci and genes light up {#probes}

::: summary
**In one sentence.** Two small random models sit on top of the polymer. One sets how many fluorescent probes are on each locus, so how bright it is. The other is a promoter that switches ON and OFF, fires RNA polymerase II (Pol II) while ON, and makes an MS2 signal. The E–P distance can speed up switching ON.
:::

## Probe binding (locus brightness)

Each locus has a pool of probes. Each probe cycles through FREE → BOUND → SCANNING → FREE, and a BOUND probe can also let go straight away. Probes are bright while attached, so the locus brightness flickers as probes come and go. The model is simulated exactly, one event at a time (the Gillespie algorithm). All rates are *placeholders* until binding kinetics are measured.

## Transcription and MS2

- **Promoter:** switches OFF → ON at rate k_on and ON → OFF at rate k_off.
- **Pol II:** while the promoter is ON, Pol II starts at rate k_init. Each Pol II moves along the gene, makes MS2 stem-loops as it passes the MS2 cassette, waits at the end of the gene, and is released.
- **MS2 signal:** the signal is the number of loops carried by all Pol II on the gene.
- **Coupling to the enhancer:** the E–P distance d(t), taken from the polymer, multiplies k_on by 1 + (F − 1)·f(d). The default f is a Hill curve, at half effect at 400 nm.

Every rate is per second of real time, so changing the frame rate of the movie does not change the biology. The default gene is very bursty: about 3 bursts per hour, 45 s ON, and about 7 Pol II per burst.

## How we know it works

Both models have exact answers that can be computed without simulating, so we compare the simulation to them.

![**Simulation versus exact theory.** **A**: the time-averaged number of probes on a locus matches the exact binomial distribution. **B**: its autocorrelation matches the exact value from the rate matrix, so the speed of the flicker is right too. **C–D**: the distance switches between contact (80 nm, grey) and apart (400 nm) every 200 s. The mean over 1000 simulated cells of the promoter state (**C**) and of the MS2 signal (**D**) lies on the exact expectation. Made by `figure_scripts/fig_probes_transcription.py`.](figures/probes_transcription.png){width=100%}

Other unit tests check that the MS2 mean is the same at 5, 20 and 60 s frames, and that the signal is correlated with contact only when the coupling is on.

## Limits

- **Placeholder values:** all rates, the MS2 brightness and the coupling strength are placeholders.
- **Probes:** the model leaves out bleaching of single probes and probes shared between the two loci.
- **Literature numbers:** the values quoted from papers for the default gene come from search summaries and have not been checked against the primary papers.
