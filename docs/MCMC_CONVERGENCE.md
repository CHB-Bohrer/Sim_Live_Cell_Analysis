# Has the MS2 -> Pol II sampler converged? (2026-10-13)

Code: `src/simlive/stage6_analysis/mcmc_check.py`, runner `scripts/check_mcmc.py`, tests `tests/test_mcmc_check.py`,
dashboard 🧮 MS2 → Pol II step 5b. Samplers: `stage6_analysis/deconvolve.py` (`gibbs_shift` default, `chen_block`, `gibbs` = flips only).

## Why it was needed
The sampler ran a fixed 500 burn-in sweeps + 1000 samples and nothing checked that the chain had reached the posterior. Agreement with the
paper's block sampler and recovery of simulated starts were indirect evidence only: two samplers can agree and both be stuck.

## Method (limits fixed BEFORE looking at results; do not tune them to pass)
1. **Four chains per trace**, each from a different, deliberately over-dispersed start (`deconvolve.INITS`): no starts; a draw from the
   prior; 10x too many starts; starts where the trace rises. Different random seeds.
2. **Scalar summaries** per kept sample: the log-posterior (likelihood + prior log-odds of the starts), the number of starts, and the number
   of starts in every frame interval (what the analysis uses).
3. **Rank-normalised split R-hat** (Vehtari et al. 2021, Bayesian Analysis 16:667): each chain is cut in half so a chain that drifts disagrees
   with itself; the larger of the bulk and the folded version. **Bulk effective sample size** from the autocorrelation (Geyer's initial
   monotone sequence).
4. **Verdict:** converged if R-hat < 1.05 for the log-posterior, the number of starts and every frame, and ESS >= 400 for the log-posterior
   and the number of starts (4 chains x 100). Per-frame ESS is reported, not judged (a frame that almost never holds a start has nothing to
   estimate).
5. **In the units of the result:** the largest between-chain difference of the start rate (starts per minute, smoothed over 1 min).
6. **Burn-in curve:** R-hat of the log-posterior over the 100 sweeps after discarding the first n sweeps, n = 0 ... 500.

## Does the check itself work? (tests/test_mcmc_check.py)
- R-hat ~1 for independent chains, > 1.2 for one displaced chain or for drifting chains; ESS equals the AR(1) theory (842 for phi = 0.9, 4 x 4000).
- **Exact posterior:** on a problem with 14 fine steps (2^14 configurations summed exactly) all three samplers reproduce every start probability to
  within 0.03 (20,000 samples each). So the samplers target the right posterior; what remains is whether they MIX on big problems.
- A run with no burn-in and 16 samples from different starts is flagged as not converged.
- Bug found while building it: a quantity that is constant in theory (log-posterior of the empty configuration) differs by rounding noise (~1e-13)
  between chains, and rank-normalisation turned that into R-hat = 7. Values are rounded to 1e-8 before ranking.

## Results (`scripts/check_mcmc.py`; 4 chains x (500 burn-in + 1000 kept x thin 3); one seed per cell)
Constructed traces with known truth (5 cases):

| sampler | converged | worst R-hat | min ESS | time per case (4 chains) |
|---|---|---|---|---|
| gibbs_shift (default) | **5 / 5** | 1.01 | 3980 | ~10 s |
| chen_block (block auto-tuned to ~20 s) | 2 / 5 | 1.08 | 75 | ~2 s |
| gibbs (flips only) | 1 / 5 (only the empty trace) | inf | 5 | ~10 s |

- **Flips only get stuck**, as expected: chains sit at different log-posterior levels for the whole run. The shift moves are what make the default work.
- The block sampler is ~5x cheaper per sample but its samples are more correlated (ESS 75-650): it would need more samples, not more burn-in.

Showcase runs, every cell:

| run | cells | gibbs_shift converged | what fails |
|---|---|---|---|
| `loci_demo_seed7_47132d` (placeholder gene, 0-3 starts per cell) | 12 | 9 | the 3 misses have R-hat <= 1.02 and ESS 320-353 (limit 400): not enough independent samples, no disagreement between chains |
| `loci_demo_seed8_f5c88d` (denser gene, 8-94 starts per cell, R = 25) | 12 | **3** | 9 cells: **chains settle in different local arrangements of the starts** |

**The dense-gene failure is real and is the main result.** In cell 1 (37 starts) the four chains have log-posterior levels of about -904, -808,
-860 and -815 with a within-chain sd of ~6, ie they differ by ~100 nats: each is trapped. They agree on HOW MANY starts (37.2, 37.2, 38.1,
37.1) much better than on WHEN (largest difference between chains 1.0-1.4 starts per minute in the 1-min rate in 8 of the 9 failing cells).
More burn-in does not help (median R-hat 1.57 at 0 sweeps, 1.53 at 500; the spread of chain means falls from 815 nats at 0 sweeps to 65 at 5 and 24 at 500, never to ~0); a larger shift range
lowers the trapped level (cell 1, best chain's mean log-posterior for max shift 3 / 15 / 60 / 200 s: -867 / -777 / -784 / -743; R-hat of the log-posterior 1.7 / 1.9 /
2.3 / 2.0) but does not produce agreement. The posterior is rugged because the trace is almost noise-free here (sigma_b = 0.01 of one Pol II, R = 25): moving a single start
makes the predicted trace worse by many standard deviations, so rearranging a burst needs several starts to move together.

## The fix: parallel tempering (`sampler="gibbs_pt"`, built 2026-10-13)
16 replicas of the same sampler (Gibbs flips + shift moves) run at likelihood strengths ("inverse temperatures") from 1 down to 0.001; the prior is
never tempered. After every sweep neighbouring replicas swap their configurations with probability min(1, exp((b_k - b_k+1)(loglik_k+1 - loglik_k))), which
leaves the product of the tempered posteriors invariant. The hot replicas move freely between arrangements of the starts (at 1/1000 strength the
noise is ~30x larger) and hand them down; only the strength-1 replica is sampled. The ladder is adapted during the first half of the burn-in towards
30% swap acceptance (ends fixed) and then frozen. All tempered samples pass the same exact-posterior test (4 replicas on the 14-step problem).

| 4 chains per cell, 500 burn-in + 1000 kept | constructed (5) | sparse run (12 cells) | dense run (12 cells) |
|---|---|---|---|
| gibbs_shift (default) | 5 / 5 | 9 / 12 | 3 / 12 |
| **gibbs_pt** | **5 / 5** | **12 / 12** | **11 / 12** |

- Dense run: log-posterior R-hat 1.00 in all 12 cells (was 1.0-2.9), minimum ESS 945, largest between-chain difference of the 1-min rate 0.09-0.26 starts
  per minute (was 0.12-1.39). In cell 1 all four chains now sit at a log-posterior of -741.0 to -741.9 (they were -904 to -808), which equals
  the best level the plain sampler reached with a 200 s shift range (-743). Swap acceptance of the frozen ladder: cell means 0.18-0.40, the weakest pair of each cell 0.10-0.30.
- The one miss (cell 9): the largest per-frame R-hat is 1.056 against the limit 1.05, one frame; everything else passes. Reported as it is, limit not moved.
- Sparse run 12 / 12 (minimum ESS 3156); the three cells that missed ESS with the plain sampler now pass.
- **Cost:** 210-325 s per cell for the 4 chains (6-way parallel load) against ~11 s for the plain sampler: about 20-30x. One chain with 16 replicas and 3,500 sweeps
  is about 1-1.5 min on its own. Not the default in the pipeline for that reason (see "Use").
- **Do the answers change?** The headline scores against the simulation hardly move (dense run, n found 401.1 / 389 true; 30-s match 0.830 vs
  0.823 plain; false per true 0.201 vs 0.207; 1-min rate r 0.797 vs 0.799; sparse run likewise). Per-frame posterior means do change
  (rms 0.12 starts per frame on a mean of 0.19), about as much as the plain and the block sampler differ from each other (0.15), so for this
  truth-scored statistics the trapped modes cost little, but any per-start timing or interval from the dense run should come from the tempered chains.
  What tempering changes for sure is whether the chains agree, and so whether the reported posterior can be believed.

## Update 2026-10-26: the shift move was not exact (fixed)
The shift move of `gibbs_shift` / `gibbs_pt` proposed a move only for an event AT the step being scanned, so the reverse move of a pair was never proposed
from the same step (and the scan revisits an event moved forward, not one moved back). That kernel does not keep the posterior: on its own, for one event on a
flat target, its stationary distribution was up to 71% off; inside the sampler the Gibbs flips hid most of it (start probabilities ~0.01 off on the 14-step
enumerable problem, 6 standard errors; the old test tolerance of 0.03 could not see it), but where only shifts can move a start (two equally likely steps,
flips never accepted) it split the start 2:1 instead of 1:1. It is now a symmetric pair move: at step j, with d drawn, if exactly one of j and j + d holds an
event, propose moving it to the other (`tests/test_mcmc_check.py::test_shift_moves_keep_the_posterior_where_only_they_can_mix`). The convergence tables above
were measured with the old move: R-hat compares chains with each other, so a sampler can pass while sampling a slightly wrong distribution. Re-run
`scripts/check_mcmc.py` to refresh them. On the showcase loci_demo_seed11_24bba2 (default spiky gene, `auto`) the scores did not change: found 0.786 -> 0.784,
false per true 0.677 -> 0.676, 1-min rate r 0.880 -> 0.881; the dense gene (seed 8) has not been re-checked.

## What this means for earlier results
- Totals, matching scores and 1-min rates in the showcase numbers (CLAUDE.md: 389 true / 401 inferred; 1-min r 0.80) stand: they are the same with the
  tempered sampler.
- Posterior intervals and the start times of dense bursts from the plain sampler were not reliable.
- The placeholder (sparse) gene is fine with either sampler.

## Default behaviour: check every cell, temper only where needed (`sampler="auto"`, user's choice 2026-10-13)
`deconvolve_run` / `scripts/deconvolve_ms2.py` / the pipeline now do this per cell: four `gibbs_shift` chains from different starts and the check above;
if it passes, chain 0 is the result. If only the ESS is short (R-hat fine: the chains agree but are too correlated) the cheap sampler is run again with 3x the
samples. If an R-hat fails, or that is still not enough, the cell is sampled again with `gibbs_pt`. `info["checks"]` and the page's step 5 list which
sampler gave which cell. The cheap sampler is kept for the cells where it is right, so a sparse movie costs ~4x one plain run and a dense one about
the tempered cost for the cells that fail.
- Measured: dense run 9 of 12 cells tempered (the same 9 that fail the plain check), 3 plain; sparse run 12 of 12 plain. Scores vs truth with the checked
  default: dense 400.0 found / 389 true, 30-s match 0.830, false per true 0.198, 1-min r 0.798; sparse 5.88 / 5, 1.00, 0.169, 0.758: the same as before.
- Why not "always temper": the cost (20-30x) buys nothing for the cells where the plain sampler already agrees with itself.

## Use
```
scripts\run.cmd python scripts\check_mcmc.py --run <run_id> [--workers 6]                          # plain, block, flips-only
scripts\run.cmd python scripts\check_mcmc.py --run <run_id> --samplers gibbs_pt --workers 6 --tag _pt  # tempered (mcmc_check_pt.csv)
scripts\run.cmd python scripts\deconvolve_ms2.py --run <run_id> --sampler gibbs_pt                  # deconvolve with tempering (files _gibbs_pt)
scripts\run.cmd python scripts\check_mcmc.py --synthetic                                             # constructed hard cases
```
In code: `deconvolve_trace(..., sampler="gibbs_pt", n_replicas=16, beta_min=1e-3)`; `check_convergence(a_cu, kern, frame_s, prior, sigma_b, beta, sampler=..., **kw)`
-> `ChainCheck` (`.converged`, `.reasons`, R-hats, ESS).

## Open
- Fewer replicas / a warmer floor than 1e-3 might do the same job faster; not explored. R, beta_min and the 30% target are settings, not tuned.
- Cell 9's frame R-hat 1.056 and the ESS limit 400 are borderline-sensitive to the seed; one seed per cell so far.
