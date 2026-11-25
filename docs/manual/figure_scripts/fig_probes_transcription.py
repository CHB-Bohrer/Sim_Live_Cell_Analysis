"""Section 2 figure: the probe and transcription models agree with their exact theory (CPU, ~1 min)."""
import numpy as np
from scipy.linalg import expm
from scipy.stats import binom

from _style import plt, save, SIM, THEORY, ACCENT, GREY
from simlive.stage1_chromatin import probes as P
from simlive.stage1_chromatin import transcription as T

fig, ax = plt.subplots(2, 2, figsize=(10, 7))

# --- A, B: probes (placeholder promoter kinetics)
p = P.ProbeParams.from_dict(None, "promoter")
th = P.summary(p)
t, F, B, S = P.simulate_trace(p, 200000.0, seed=1)
grid = np.arange(0, 200000.0, 1.0)
att = (B + S)[np.searchsorted(t, grid, side="right") - 1]
k = np.arange(0, p.n_probes + 1)
a = ax[0, 0]
a.bar(k, np.bincount(att, minlength=len(k)) / len(att), color=SIM, alpha=0.6, label="simulation (Gillespie)")
a.plot(k, binom.pmf(k, p.n_probes, th["p_attached"]), "o-", color=THEORY, ms=3, label="exact theory (binomial)")
a.set(xlim=(-0.5, p.n_probes * 0.8), xlabel="probes attached to the locus", ylabel="fraction of time",
      title="A  Probe count at a locus")
a.legend()

lags = np.arange(0, 80, 1.0)
x = att - att.mean()
ac = np.array([np.mean(x[: len(x) - int(l)] * x[int(l):]) for l in lags]) / x.var()
Q, pa = P.generator(p), th["p_attached"]
pi = np.array([1 - pa, th["mean_bound"] / p.n_probes, th["mean_scanning"] / p.n_probes])
on = np.array([0.0, 1.0, 1.0])
ex = np.array([(pi * on) @ expm(Q * l) @ on - pa ** 2 for l in lags]) / (pa * (1 - pa))
a = ax[0, 1]
a.plot(lags, ac, color=SIM, lw=2.5, alpha=0.6, label="simulation")
a.plot(lags, ex, "--", color=THEORY, label="exact theory")
a.set(xlabel="time lag (s)", ylabel="autocorrelation of the count", title="B  How fast the brightness flickers")
a.legend()

# --- C, D: transcription with a distance that switches between contact and far
FAST = dict(k_on=0.004, k_off=0.05, k_init=0.2, coupling="contact", fold=15.0, f_const=0.0, elongation_kb_min=6.0,
            cassette_kb=0.6, gene_kb=2.4, n_loops=24, dwell_s=10.0)
tp = T.TranscriptionParams(**FAST)
d = np.where((np.arange(400) // 20) % 2 == 0, 80.0, 400.0)
ft = np.arange(100, 4000.0, 10.0)
runs = [T.simulate_cell(tp, d, 10.0, ft, 0.1, seed=s).frames for s in range(1000)]
on_sim = np.mean([r.promoter_on for r in runs], axis=0)
ms2_sim = np.mean([r.ms2_loops for r in runs], axis=0)
p_start, _ = T.expected_on(tp, np.tile(d, 10), 10.0)
on_ex = p_start[np.searchsorted(np.arange(len(p_start)) * 10.0, ft, side="right") - 1]
ms2_ex = T.expected_ms2(tp, np.tile(d, 10), 10.0, ft)
contact = d[(ft // 10).astype(int)] < 100
for a in ax[1]:
    for i in np.flatnonzero(np.diff(np.r_[0, contact.astype(int), 0]) == 1):
        j = i + np.argmax(~contact[i:]) if (~contact[i:]).any() else len(ft)
        a.axvspan(ft[i] / 60, ft[min(j, len(ft) - 1)] / 60, color=GREY, alpha=0.2, lw=0)
a = ax[1, 0]
a.plot(ft / 60, on_sim, color=SIM, alpha=0.7, label="simulation, mean of 1000 cells")
a.plot(ft / 60, on_ex, "--", color=THEORY, label="exact expectation")
a.set(xlabel="time (min)", ylabel="P(promoter ON)", title="C  Promoter follows the distance")
a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=8.5); a.set_title(a.get_title(), pad=26)
a = ax[1, 1]
a.plot(ft / 60, ms2_sim, color=SIM, alpha=0.7, label="simulation")
a.plot(ft / 60, ms2_ex, "--", color=THEORY, label="exact expectation")
a.set(xlabel="time (min)", ylabel="mean MS2 signal (loops)", title="D  MS2 signal lags the promoter")
a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=8.5); a.set_title(a.get_title(), pad=26)
fig.text(0.5, -0.01, "Grey bands: enhancer in contact (80 nm); white: apart (400 nm).", ha="center", fontsize=9)
fig.tight_layout()
save(fig, "probes_transcription.png")
