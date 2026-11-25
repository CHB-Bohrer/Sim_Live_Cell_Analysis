"""Section 6 figure: the MS2 -> Pol II deconvolution samples the exact posterior and finds known starts (CPU, ~1 min)."""
import itertools

import numpy as np

from _style import plt, save, SIM, THEORY, ACCENT, GREY
from simlive.stage6_analysis.deconvolve import Kernel, deconvolve_trace, make_kernel

# --- A: a problem small enough to sum over all 2^14 start patterns (same set-up as tests/test_mcmc_check.py)
kern = Kernel(1.0, np.array([0.0, 0.5, 1.0, 1.0, 0.0]))
T, s, K = 5, 2, 5
N = K + (T - 1) * s + 1
idx = K + s * np.arange(T)


def model(I):
    G = np.zeros(T)
    for n in range(T):
        u = idx[n] - np.flatnonzero(I)
        G[n] = kern.values[u[(u > 0) & (u < K)]].sum()
    return G


rng = np.random.default_rng(5)
truth = np.zeros(N, int); truth[[6, 11]] = 1
sb, beta, p = 0.25, 0.05, 0.08
G0 = model(truth)
a = G0 + rng.normal(size=T) * np.sqrt(sb ** 2 + beta * G0)
logw, states = [], []
for bits in itertools.product((0, 1), repeat=N):
    I = np.array(bits); G = model(I); v = sb ** 2 + beta * np.clip(G, 0, None)
    logw.append((-0.5 * (a - G) ** 2 / v - 0.5 * np.log(v)).sum() + I.sum() * np.log(p / (1 - p))); states.append(I)
w = np.exp(np.array(logw) - max(logw)); w /= w.sum()
exact = (w[:, None] * np.array(states)).sum(0)
samplers = {"gibbs_shift (default)": ("gibbs_shift", {}), "parallel tempering": ("gibbs_pt", {"n_replicas": 4, "beta_min": 0.1}),
            "Chen et al. block sampler": ("chen_block", {"block_s": 4.0})}
est = {n: deconvolve_trace(a, kern, 2.0, np.full(T, p), sb, beta, n_burn=300, n_keep=20000, thin=3, max_shift_s=3.0, seed=3,
                           sampler=sm, **kw)["p_fine"] for n, (sm, kw) in samplers.items()}

# --- B: a realistic trace with three known starts (same set-up as tests/test_deconvolve.py)
rng = np.random.default_rng(8)
k = make_kernel(24, 1.3, 5.0, 2.5, 116.0, 1.0)
TT, dt = 120, 10.0
inits = np.array([203.0, 640.0, 905.0])
t = np.arange(TT) * dt
Gt = np.zeros(TT)
for s0 in inits:
    u = t - s0; ok = u > 0
    d = np.clip(np.floor(u[ok]).astype(int), 0, len(k.values) - 1)
    Gt[ok] += np.where(u[ok] < len(k.values), k.values[d], 0.0)
trace = Gt + rng.normal(0, 1, TT) * np.sqrt(0.1 ** 2 + 0.02 * Gt)
res = deconvolve_trace(trace, k, dt, np.full(TT, 3 / (TT * dt)), 0.1, 0.02, n_burn=200, n_keep=300, thin=2)

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"width_ratios": [1, 1.4]})
x = np.arange(N)
ax[0].bar(x, exact, color=GREY, label="exact (all 16,384 patterns summed)")
for (n, e), m in zip(est.items(), "o^s"):
    ax[0].plot(x, e, m, ms=5, mfc="none", label=n)
ax[0].set(xlabel="time step (s)", ylabel="P(a Pol II starts here)", title="A  Samplers match the exact answer")
ax[0].legend(fontsize=8)
b = ax[1]
b.plot(t / 60, trace, color=GREY, lw=1, label="MS2 trace (units of one Pol II)")
b.plot(t / 60, res["G"], color=SIM, label="fitted signal")
b2 = b.twinx()
edges = np.arange(0, TT * dt + 30, 30.0)
pb, _ = np.histogram(res["t_fine"], edges, weights=res["p_fine"])
b2.bar(edges[:-1] / 60, pb, width=0.5, align="edge", color=ACCENT, alpha=0.45, label="P(start) per 30 s")
b2.set_ylim(0, 2.4); b2.set_ylabel("P(a start) in each 30 s", color=ACCENT); b2.spines["top"].set_visible(False)
b.set_xlim(0, TT * dt / 60)
for s0 in inits:
    b.axvline(s0 / 60, color=THEORY, ls=":", lw=1)
b.set(xlabel="time (min)", ylabel="signal", title="B  Known starts (dotted) are recovered")
b.legend(loc="upper left", fontsize=8)
fig.tight_layout()
save(fig, "deconvolution.png")
print({n: float(np.abs(e - exact).max()) for n, e in est.items()})
