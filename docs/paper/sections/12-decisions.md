<!--
Section: Decisions and why
Update when: a design decision is made or reversed. Append; never delete. A reversal gets a new entry naming the one it replaces.
Code: (none; this is a log)
Check: none
Sources: docs/notes/*.md (each entry names its note).
-->

# Decisions and why

::: {.box}
**Short version.** The choices that shaped the software, each with its date, its reason, and what was rejected. Reasons are taken from the project notes; where the notes give none, none is invented.
:::

**D01 (2026-10-05). Chromatin is simulated separately, as a library of 100 simulations at 1 kb per monomer.**

- *Why:* Agreed with the owner. Reason not recorded beyond the agreement.
- *Rejected or superseded:* A chromosome-21-sized region (dropped by the owner).

**D02 (2026-10-08). Loop bonds are all registered in advance and switched by stiffness.**

- *Why:* Re-pointing existing OpenMM bonds silently left loops frozen; the first validation was invalid.
- *Rejected or superseded:* Updating bond partners in place.

**D03 (2026-10-08). The insulation target (1.5-fold) stays, and the 1.36 result is reported as a miss.**

- *Why:* A limit set in advance is not tuned to pass.
- *Rejected or superseded:* Lowering the threshold.

**D04 (2026-10-11). One spot-finding method, after Chen et al. 2025.**

- *Why:* It tracks through clutter with a diffusion model. An earlier method built from another paper worked but was not that paper's method, and was removed at the owner's request.
- *Rejected or superseded:* Brightest pixel in the mask (kept only for comparison).

**D05 (2026-10-13). The default gene is strongly bursty (ON 45 s, about 7 polymerases per burst).**

- *Why:* The owner asked for a trace of narrow spikes that returns to zero.
- *Rejected or superseded:* The earlier TFF1-like values (kept in config comments).

**D06 (2026-10-13). The deconvolution checks convergence per cell and uses parallel tempering only where needed.**

- *Why:* The simple sampler fails on dense genes; tempering is 20 to 30 times slower.
- *Rejected or superseded:* Tempering everywhere; the simple sampler everywhere.

**D07 (2026-10-14). The E-P readout is the probability of a burst, and the whole curve is fitted.**

- *Why:* The owner's choice.
- *Rejected or superseded:* A single contact threshold.

**D08 (2026-10-15). Library simulations have 8,000 production blocks, and each cell draws its own 100 kb pair.**

- *Why:* 2,000 blocks covered only 15 minutes of movie at the assumed time scale; random pairs let 100 simulations serve hundreds of cells.
- *Rejected or superseded:* 2,000 blocks; one fixed pair.

**D09 (2026-10-17). Tutorial videos are published as GitHub releases, never committed.**

- *Why:* The owner did not want large files added to the history.
- *Rejected or superseded:* Videos in the repository.

**D10 (2026-10-20). The default coupling is a smooth curve with half effect at 400 nm.**

- *Why:* The sharp 150 nm contact was reached in only 0.75% of frames, so it barely coupled.
- *Rejected or superseded:* Sharp contact at 150 nm.

**D11 (2026-10-25). The tutorial voice is a local neural voice (Kokoro).**

- *Why:* The owner judged browser and online voices too synthetic; no text leaves the computer.
- *Rejected or superseded:* Browser speech; an online neural voice.

**D12 (2026-10-26). Spot brightness is not used in multi-copy linking.**

- *Why:* Tested on 25 movies: no net gain, because per-frame brightness is noisy and copies are equally bright.
- *Rejected or superseded:* Brightness term on by default.

**D13 (2026-11-21). All pipeline images are losslessly compressed; tracking-challenge copies are not kept.**

- *Why:* Simulation data fell from 115 GB to 28 GB with identical pixels.
- *Rejected or superseded:* Uncompressed stacks; keeping every copy.

**D14 (2026-11-27). `CLAUDE.md` holds only every-session rules; detail lives in `docs/notes/`, status in `STATUS.md`, and this document explains the design.**

- *Why:* The full record was loaded in every session at 118 KB.
- *Rejected or superseded:* One long handoff file.

