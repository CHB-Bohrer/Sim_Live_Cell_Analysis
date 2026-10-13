"""Convergence check of the MS2 samplers (stage6_analysis/mcmc_check.py) and exactness of the samplers on a problem small enough to enumerate."""
import itertools

import numpy as np
import pytest

from simlive.stage6_analysis.deconvolve import Kernel, deconvolve_trace, make_kernel
from simlive.stage6_analysis.mcmc_check import burn_in_curve, check_convergence, ess_bulk, split_rhat


def test_rhat_is_one_for_independent_chains_and_large_for_separated_ones():
    rng = np.random.default_rng(0)
    assert split_rhat(rng.normal(size=(4, 1000))) < 1.02
    x = rng.normal(size=(4, 1000))
    x[0] += 3.0                                                      # one chain sits somewhere else
    assert split_rhat(x) > 1.2
    drift = rng.normal(size=(4, 1000)) + np.linspace(0, 4, 1000)       # every chain still drifting: the split halves disagree
    assert split_rhat(drift) > 1.2
    assert split_rhat(np.zeros((4, 100))) == 1.0                       # constant and equal: nothing to disagree about
    assert split_rhat(np.arange(4.0)[:, None] * np.ones((4, 100))) > 100       # constant but different


def test_ess_matches_the_autocorrelation_of_an_ar1_process():
    rng = np.random.default_rng(1)
    n, m, phi = 4000, 4, 0.9
    x = np.zeros((m, n))
    e = rng.normal(size=(m, n))
    for t in range(1, n):
        x[:, t] = phi * x[:, t - 1] + e[:, t]
    # theory: ESS = N (1 - phi) / (1 + phi) = 16000 x 0.0526 = 842
    assert 0.6 * 842 < ess_bulk(x) < 1.5 * 842
    assert ess_bulk(rng.normal(size=(m, n))) > 0.8 * m * n


def _tiny_problem():
    """5 frames x 2 s, kernel of 5 s, so 14 fine steps and 2^14 configurations: the posterior can be summed exactly."""
    kern = Kernel(1.0, np.array([0.0, 0.5, 1.0, 1.0, 0.0]))
    T, s, K = 5, 2, 5
    J0, N = K, K + (T - 1) * s + 1
    rng = np.random.default_rng(5)
    truth = np.zeros(N, int)
    truth[[6, 11]] = 1
    idx = J0 + s * np.arange(T)

    def model(I):
        G = np.zeros(T)
        for n in range(T):
            u = idx[n] - np.flatnonzero(I)
            G[n] = kern.values[u[(u > 0) & (u < K)]].sum()
        return G

    sb, beta = 0.25, 0.05
    G0 = model(truth)
    a = G0 + rng.normal(size=T) * np.sqrt(sb ** 2 + beta * G0)
    p = 0.08
    return kern, a, sb, beta, p, s, N, model


def _exact_marginals(a, sb, beta, p, N, model):
    logw, states = [], []
    for bits in itertools.product((0, 1), repeat=N):
        I = np.array(bits)
        G = model(I)
        v = sb ** 2 + beta * np.clip(G, 0, None)
        logw.append((-0.5 * (a - G) ** 2 / v - 0.5 * np.log(v)).sum() + I.sum() * np.log(p / (1 - p)))
        states.append(I)
    w = np.exp(np.array(logw) - max(logw))
    w /= w.sum()
    return (w[:, None] * np.array(states)).sum(0)


@pytest.mark.parametrize("sampler,kw", [("gibbs_shift", {}), ("gibbs", {}), ("chen_block", {"block_s": 4.0})])
def test_samplers_reproduce_the_exact_posterior_of_a_tiny_problem(sampler, kw):
    kern, a, sb, beta, p, s, N, model = _tiny_problem()
    exact = _exact_marginals(a, sb, beta, p, N, model)
    res = deconvolve_trace(a, kern, 2.0, np.full(len(a), p), sb, beta, n_burn=300, n_keep=20000, thin=3, max_shift_s=3.0, seed=3,
                           sampler=sampler, **kw)
    assert len(res["p_fine"]) == N
    err = np.abs(res["p_fine"] - exact).max()
    assert err < 0.03, f"{sampler}: largest error of a start probability {err:.3f}"


def test_check_passes_for_a_well_mixed_problem_and_flags_chains_that_have_not_mixed():
    kern, a, sb, beta, p, s, N, model = _tiny_problem()
    prior = np.full(len(a), p)
    ok = check_convergence(a, kern, 2.0, prior, sb, beta, n_burn=300, n_keep=1000, thin=3, max_shift_s=3.0, seed=1)
    assert ok.converged, ok.reasons
    assert ok.rhat_logpost < 1.05 and ok.ess_logpost > 400 and ok.max_chain_diff < 0.5
    # a deliberately unconverged run: no burn-in and few samples from very different starts
    bad = check_convergence(a, kern, 2.0, prior, sb, beta, inits=("empty", "crowded"), n_burn=0, n_keep=16, thin=1, max_shift_s=3.0, seed=1)
    assert not bad.converged


def test_burn_in_curve_shrinks_with_burn_in():
    rng = np.random.default_rng(2)
    k = make_kernel(24, 1.3, 5.0, 2.5, 116.0, 1.0)
    T, dt = 60, 10.0
    t = np.arange(T) * dt
    G = np.zeros(T)
    for s0 in (100.0, 300.0, 420.0):
        u = t - s0
        ok = (u > 0) & (u < len(k.values))
        G[ok] += k.values[u[ok].astype(int)]
    a = G + rng.normal(size=T) * 0.12
    chains = [deconvolve_trace(a, k, dt, np.full(T, 3 / (T * dt)), 0.12, 0.0, n_burn=300, n_keep=50, thin=2, seed=i, init=ini)
              for i, ini in enumerate(("empty", "prior", "crowded", "rise"))]
    bc = burn_in_curve(chains, sweeps=(0, 50, 100), window=100)
    assert len(bc) == 3 and bc[0]["spread_of_chain_means"] >= bc[-1]["spread_of_chain_means"]
