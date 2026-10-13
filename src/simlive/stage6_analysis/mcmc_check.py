"""Convergence check for the MS2 -> Pol II samplers (stage6_analysis/deconvolve.py).

Method (Gelman-Rubin / Vehtari et al. 2021, "Rank-normalization, folding, and localization", Bayesian Analysis 16:667):
1. Run several chains on the SAME trace from deliberately different starting configurations (`deconvolve.INITS`: empty, prior draw,
   10x too many starts, data-driven) and different seeds.
2. For scalar summaries of every chain compute the split-R-hat (each chain cut in two, so a chain that drifts disagrees with itself) on
   rank-normalised values, and the bulk effective sample size (Geyer's initial monotone sequence on the autocorrelation).
   Scalars: log-posterior, total number of starts, and the number of starts in every frame interval (the quantity the analysis uses).
3. Verdict, with limits fixed in advance (not tuned): converged if every R-hat < `RHAT_MAX` (1.05) and the ESS of the log-posterior and of
   the number of starts is >= `ESS_MIN` (400, ie 100 per chain for 4 chains). Per-frame ESS is reported (minimum, median) but is not a
   criterion, because a frame with almost no starts has no autocorrelation to estimate.
4. Also reported: the largest difference between chains in the start probability per fine step, smoothed over one minute: a direct
   "do the chains tell the same story" number in the units of the result.

What it cannot do: a diagnostic from chains that all miss the same mode looks fine. Over-dispersed starts make that unlikely, not impossible;
`tests/test_mcmc_check.py` also compares both samplers with the exact posterior of a problem small enough to enumerate.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from scipy.stats import norm, rankdata

from simlive.stage6_analysis.deconvolve import INITS, Kernel, deconvolve_trace

RHAT_MAX = 1.05
ESS_MIN = 400.0


def _split(x: np.ndarray) -> np.ndarray:
    """(chains, draws) -> (2 x chains, draws // 2): each chain split in half. Values are rounded to 1e-8 first: a quantity that is
    constant in theory (e.g. the log-posterior of the empty configuration, 210.8 in every chain) differs by rounding noise ~1e-13
    between chains, which rank-normalisation would blow up into an R-hat of 7."""
    x = np.round(np.asarray(x, float), 8)
    n = x.shape[1] // 2
    return np.concatenate([x[:, :n], x[:, n: 2 * n]], axis=0)


def _rhat_plain(x: np.ndarray) -> float:
    m, n = x.shape
    if n < 4:
        return np.nan
    W = x.var(axis=1, ddof=1).mean()
    B = n * x.mean(axis=1).var(ddof=1)
    if W <= 0:
        return 1.0 if B <= 0 else np.inf        # every chain constant: equal -> fine; different constants -> not mixed
    return float(np.sqrt(((n - 1) / n * W + B / n) / W))


def _rank_normalise(x: np.ndarray) -> np.ndarray:
    r = rankdata(x.ravel(), method="average").reshape(x.shape)
    return norm.ppf((r - 0.375) / (x.size + 0.25))


def split_rhat(x: np.ndarray) -> float:
    """Rank-normalised split R-hat of draws x (chains, draws): the larger of the bulk and the folded (tails / scale) version."""
    x = np.asarray(x, float)
    if x.shape[1] < 8:
        return np.nan
    s = _split(x)
    bulk = _rhat_plain(_rank_normalise(s))
    folded = _rhat_plain(_rank_normalise(np.abs(s - np.median(s))))
    return float(np.nanmax([bulk, folded]))


def _autocov(x: np.ndarray) -> np.ndarray:
    n = x.shape[-1]
    f = np.fft.rfft(x - x.mean(-1, keepdims=True), 2 * n)
    return np.fft.irfft(f * np.conj(f))[..., :n] / n


def ess_bulk(x: np.ndarray) -> float:
    """Bulk effective sample size (chains, draws) of the rank-normalised split chains (Geyer's initial monotone positive sequence)."""
    x = np.asarray(x, float)
    if x.shape[1] < 8:
        return np.nan
    z = _rank_normalise(_split(x))
    m, n = z.shape
    ac = _autocov(z)
    W = z.var(axis=1, ddof=1).mean()
    B = n * z.mean(axis=1).var(ddof=1) if m > 1 else 0.0
    if W <= 0:
        return float(m * n)
    var_plus = (n - 1) / n * W + B / n
    rho = 1.0 - (W - ac.mean(0) * n / (n - 1)) / var_plus
    rho[0] = 1.0
    pair, prev = [], np.inf
    for t in range(0, n - 1, 2):
        g = rho[t] + rho[t + 1]
        if g < 0:
            break
        g = min(g, prev)                    # monotone
        pair.append(g)
        prev = g
    tau = max(-1.0 + 2.0 * sum(pair), 1.0 / np.log10(m * n)) if pair else 1.0
    return float(m * n / tau)


@dataclass
class ChainCheck:
    """Result of `check_convergence` for one trace."""
    sampler: str
    inits: tuple
    rhat_logpost: float
    rhat_nstarts: float
    rhat_frames_max: float
    ess_logpost: float
    ess_nstarts: float
    ess_frames_min: float
    ess_frames_median: float
    max_chain_diff: float                   # largest between-chain difference of the smoothed start probability (starts per minute)
    seconds: float
    converged: bool
    reasons: list = field(default_factory=list)
    chains: list = field(default_factory=list, repr=False)         # the raw results of every chain

    def row(self) -> dict:
        d = {k: v for k, v in self.__dict__.items() if k not in ("chains", "reasons", "inits")}
        d["inits"] = ",".join(self.inits)
        d["reasons"] = "; ".join(self.reasons)
        return d


def check_convergence(a_cu: np.ndarray, kern: Kernel, frame_s: float, rate_prior_per_s: np.ndarray, sigma_b: float, beta: float,
                      sampler: str = "gibbs_shift", inits: tuple = INITS, n_burn: int = 500, n_keep: int = 1000, thin: int = 3,
                      max_shift_s: float = 15.0, seed: int = 0, rhat_max: float = RHAT_MAX, ess_min: float = ESS_MIN,
                      block_s="auto") -> ChainCheck:
    """Run one chain per entry of `inits` (different start, different seed) and judge them together. Time is for all chains."""
    chains, t0 = [], time.perf_counter()
    for k, ini in enumerate(inits):
        chains.append(deconvolve_trace(a_cu, kern, frame_s, rate_prior_per_s, sigma_b, beta, n_burn=n_burn, n_keep=n_keep, thin=thin,
                                       max_shift_s=max_shift_s, seed=seed + 1000 * (k + 1), sampler=sampler, block_s=block_s, init=ini))
    secs = time.perf_counter() - t0
    lp = np.array([c["log_post"] for c in chains])
    ns = np.array([c["n_starts"] for c in chains])
    cnt = np.array([c["counts"] for c in chains], float)               # (chains, draws, frames)
    r_lp, r_ns = split_rhat(lp), split_rhat(ns)
    e_lp, e_ns = ess_bulk(lp), ess_bulk(ns)
    live = [n for n in range(cnt.shape[2]) if cnt[:, :, n].std() > 0]  # frames in which any start is ever drawn
    r_fr = [split_rhat(cnt[:, :, n]) for n in live]
    e_fr = [ess_bulk(cnt[:, :, n]) for n in live]
    pf = np.array([c["p_fine"] for c in chains])
    w = max(1, int(round(60.0 / kern.dt_fine_s)))
    sm = np.array([np.convolve(p, np.ones(w), mode="same") for p in pf])        # starts per minute at every fine step
    diff = float((sm.max(0) - sm.min(0)).max())
    r_fr_max = float(np.nanmax(r_fr)) if r_fr else 1.0
    reasons = []
    for name, val, lim in (("R-hat log-posterior", r_lp, rhat_max), ("R-hat number of starts", r_ns, rhat_max), ("R-hat frames", r_fr_max, rhat_max)):
        if not (val < lim):
            reasons.append(f"{name} {val:.3f} >= {lim}")
    for name, val in (("ESS log-posterior", e_lp), ("ESS number of starts", e_ns)):
        if not (val >= ess_min):
            reasons.append(f"{name} {val:.0f} < {ess_min:.0f}")
    return ChainCheck(sampler, tuple(inits), r_lp, r_ns, r_fr_max, e_lp, e_ns,
                      float(np.nanmin(e_fr)) if e_fr else np.nan, float(np.nanmedian(e_fr)) if e_fr else np.nan,
                      diff, secs, not reasons, reasons, chains)


def burn_in_curve(chains: list, sweeps=(0, 5, 10, 25, 50, 100, 250, 500), window: int = 100) -> list[dict]:
    """How many sweeps until the chains agree: split-R-hat of the log-posterior over the `window` sweeps that follow each burn-in length,
    using the per-sweep log-posterior of every chain (needs chains that were run with n_burn >= max(sweeps) + window; shorter ones are skipped)."""
    lp = np.array([c["log_post_sweeps"] for c in chains])
    out = []
    for b in sweeps:
        if b + window <= lp.shape[1]:
            seg = lp[:, b: b + window]
            out.append({"burn_in_sweeps": b, "rhat_logpost": split_rhat(seg), "mean_logpost": float(seg.mean()),
                        "spread_of_chain_means": float(seg.mean(1).max() - seg.mean(1).min())})
    return out
