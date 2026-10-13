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

## What this means for earlier results
- The totals and the 1-min rates in the showcase numbers (CLAUDE.md: 389 true / 401 inferred starts; 1-min rate r 0.80) probably stand (number of
  starts agrees between chains), but the posterior intervals and the exact start times for dense bursts were NOT sampled properly: the reported
  uncertainty is too small, since it reflects one local mode.
- The placeholder (sparse) gene is fine.

## Next step (not built; needs a decision)
Parallel tempering (replicas at lower likelihood weight, swapping configurations) is the standard remedy. The cost is the number of replicas
(the log-posterior range here is ~10^4 nats, so about 15-20): ~15-20x the time per chain. Cheaper first tries: moves that shift two neighbouring starts
together, and a move that replaces a burst's start times with a random re-draw conditional on the count. Neither is built; the check is the
yardstick for either.

## Use
```
scripts\run.cmd python scripts\check_mcmc.py --run <run_id>           # every cell, 3 samplers, ~10 min
scripts\run.cmd python scripts\check_mcmc.py --synthetic              # constructed hard cases
```
In code: `check_convergence(a_cu, kern, frame_s, prior, sigma_b, beta, sampler=...)` -> `ChainCheck` (`.converged`, `.reasons`, R-hats, ESS).
