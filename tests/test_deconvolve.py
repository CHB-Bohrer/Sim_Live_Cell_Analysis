"""MS2 deconvolution after Chen et al. 2025 (stage6_analysis/deconvolve.py): kernel, noise model and recovery of known initiation events."""
import numpy as np
import pandas as pd

from simlive.stage6_analysis.deconvolve import burst_calls, deconvolve_trace, make_kernel, noise_model, score_deconvolution


def _trace(kern, inits_s, T, dt, sb, beta, rng):
    t = np.arange(T) * dt
    G = np.zeros(T)
    for s0 in inits_s:
        u = t - s0
        ok = u > 0
        d = np.clip(np.floor(u[ok] / kern.dt_fine_s).astype(int), 0, len(kern.values) - 1)
        G[ok] += np.where(u[ok] / kern.dt_fine_s < len(kern.values), kern.values[d], 0.0)
    return G + rng.normal(0, 1, T) * np.sqrt(sb ** 2 + beta * G), G


def test_kernel_has_the_loop_steps_and_ends_at_release():
    k = make_kernel(24, 1.2, 4.8, 2.4, 60.0, 1.0)          # 0.04 kb/s: cassette 30 s, gene 150 s, release at 210 s
    assert k.values[0] == 0 and abs(k.values[15] - 0.5) < 0.05 and k.values[31] == 1.0
    assert k.values[209] == 1.0 and (len(k.values) - 1) * 1.0 >= 210 and k.values[-1] == 0.0
    assert abs(k.area_s - (210 - 15)) < 3                   # ~ release time minus half the cassette


def test_noise_model_recovers_sigma_and_beta():
    rng = np.random.default_rng(0)
    k = make_kernel(24, 1.3, 5.0, 2.5, 116.0, 1.0)
    tr = [_trace(k, rng.uniform(-200, 1800, 8), 180, 10.0, 0.2, 0.05, rng)[0] for _ in range(30)]
    sb, beta = noise_model(tr)
    assert abs(sb - 0.2) < 0.06 and abs(beta - 0.05) < 0.04


def test_deconvolution_finds_known_initiations():
    rng = np.random.default_rng(3)
    k = make_kernel(24, 1.3, 5.0, 2.5, 116.0, 1.0)
    T, dt = 180, 10.0
    inits = np.array([105.0, 140.0, 700.0, 1290.0, 1320.0, 1500.0])
    a, G = _trace(k, inits, T, dt, 0.15, 0.02, rng)
    prior = np.full(T, 6 / (T * dt))
    res = deconvolve_trace(a, k, dt, prior, 0.15, 0.02, n_burn=200, n_keep=300, thin=2)
    assert abs(res["counts"].sum(1).mean() - len(inits)) < 1.0            # number of initiations
    assert np.corrcoef(res["G"], G)[0, 1] > 0.98                          # the fit reproduces the signal
    events = pd.DataFrame({"cell_id": 1, "init_s": inits, "release_s": inits + 300})
    dec = pd.DataFrame({"cell_id": 1, "t": np.arange(T), "n_init_mean": res["counts"].mean(0), "p_on": burst_calls(res["counts"], dt)})
    tx = pd.DataFrame({"cell_id": 1, "t": np.arange(T), "promoter_on": 0})
    sc = score_deconvolution(dec, {1: res["counts"]}, events, tx, dt)
    assert sc["found"] >= 5 / 6 and sc["false_per_true"] <= 0.34
