<!--
Section: The problem and the approach
Update when: the experiment being simulated changes (loci, colours, readout).
Code: configs/loci_demo.yaml, configs/loci_library.yaml
Check: none
Sources: CLAUDE.md "The experiment being simulated".
-->

# The problem and the approach

::: {.box}
**Short version.** An analysis of live-cell movies is a long chain: find nuclei, follow them, find the spots, measure distances and brightness, infer transcription, fit a relationship. Each link adds error. The only way to measure those errors is to analyse data where the answer is known.
:::

## The experiment being imitated

- Two genomic loci, a **promoter** (locus 0) and an **enhancer** (locus 1), 100 kb apart, each in its own colour.
- **MS2** stem loops at the promoter report transcription. MS2 has the same colour as the nuclear marker (NLS-GFP), because the MS2 coat protein carries a nuclear localization signal. MS2 spots are therefore drawn into the nuclear channel.
- Many cells per field. They move, collide and may divide. Some cells carry two copies of the locus pair, and stray spots appear that belong to no locus.
- Imaging mode, pixel size, numerical aperture, frame interval, exposure, movie length and cell number are all configuration settings, never constants in the code.

## The approach

1. Simulate the biology with every hidden quantity saved: where each locus is, whether the promoter is on, when each polymerase started.
2. Turn it into camera images with a forward model of the microscope that can switch each source of image error on or off.
3. Run the analysis on the images only.
4. Compare each analysis output with the saved truth, and trace how an error at one step changes the final scientific result.

**Why simulate the whole experiment and not just the polymer?** A component can be right while the analysis is still biased. For example, projecting a 3D distance onto a 2D image removes about a third of the enhancer effect if ignored (Section 9). Only data with a known answer, passed through every step that touches it, shows whether a reported number can be believed.

## Rules the code follows

- **File contracts.** Each stage reads and writes defined files. Nothing is handed between stages in memory. Every stage folder holds its data, the resolved settings (`params.yaml`), the seed, and `provenance.json` (git commit, package versions, time).
- **Seeds.** All randomness flows from one explicit seed through `numpy.random.SeedSequence.spawn`. Each part has its own stream, so switching one part on or off never changes another part's random numbers.
- **Units.** Stage 1 works in polymer units. From stage 2 on, everything is in nanometres and seconds, and column names carry the unit (`x_nm`, `t_s`).
- **Evidence.** A single run is not evidence. Conditions are compared over at least five seeds. A limit set before a check is not moved afterwards.
