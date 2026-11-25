"""Section 6 figure: the spot finder's position error falls as 1/sqrt(photons) and stays close to the physical limit (CPU, ~1 min)."""
import numpy as np
from scipy.special import erf

from _style import plt, save, SIM, THEORY, GREY
from simlive.stage6_analysis import dog_viterbi as DV
from simlive.stage6_analysis import spots as SP

S, SIGMA, BG, READ = 41, 0.85, 20.0, 1.5          # crop size, PSF sigma (px), background photons/px, read noise


def spot(cy, cx, n):
    s2 = np.sqrt(2) * SIGMA
    i = np.arange(S)
    gy = 0.5 * (erf((i + 0.5 - cy) / s2) - erf((i - 0.5 - cy) / s2))
    gx = 0.5 * (erf((i + 0.5 - cx) / s2) - erf((i - 0.5 - cx) / s2))
    return n * np.outer(gy, gx)


yy, xx = np.indices((S, S))
region = np.hypot(yy - S // 2, xx - S // 2) <= 16
rng = np.random.default_rng(1)
photons = np.array([150, 250, 400, 700, 1200, 2000, 3500, 6000, 10000])
rms, found = [], []
for n in photons:
    err, ok = [], 0
    for _ in range(200):
        cy, cx = rng.uniform(17, 23, 2)
        img = rng.poisson(spot(cy, cx, n) + BG).astype(float) + rng.normal(0, READ, (S, S)) + 100.0
        c = DV.dog_candidates(img, region)
        if c is None or not len(c["peaks"]):
            continue
        r = SP.refine_peak(c, 0, **SP.DEFAULTS)
        e = np.hypot(r["pos"][0] - cy, r["pos"][1] - cx)
        if e < 3:
            ok += 1
            err += [r["pos"][0] - cy, r["pos"][1] - cx]
    rms.append(np.sqrt(np.mean(np.square(err))))
    found.append(ok / 200)

# Mortensen et al. 2010 limit for a pixelated Gaussian with background (per axis), a = 1 px
nn = np.geomspace(photons[0], photons[-1], 100)
sa2 = SIGMA ** 2 + 1 / 12
b = BG + READ ** 2
limit = np.sqrt(sa2 / nn * (16 / 9 + 8 * np.pi * sa2 * b / nn))

fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
a = ax[0]
a.loglog(photons, np.array(rms) * 100, "o", color=SIM, label="spot finder (DoG + Gaussian fit)")
a.loglog(nn, limit * 100, "-", color=THEORY, label="physical limit (Mortensen 2010)")
a.loglog(nn, 100 * rms[-1] * np.sqrt(photons[-1] / nn), ":", color=GREY, label="slope -1/2")
a.set(xlabel="photons in the spot", ylabel="position error per axis (nm at 100 nm/px)", title="A  Precision scales as 1/sqrt(photons)")
a.legend(fontsize=8.5)
from matplotlib.ticker import ScalarFormatter, NullFormatter
for axis in (a.xaxis, a.yaxis):
    axis.set_major_formatter(ScalarFormatter()); axis.set_minor_formatter(NullFormatter())
a.set_yticks([1, 2, 5, 10, 20])
a = ax[1]
a.semilogx(photons, np.array(rms) / np.interp(photons, nn, limit), "o-", color=SIM)
a.axhline(1, color=THEORY, lw=1)
a.set(ylim=(0.8, 2.0), xlabel="photons in the spot", ylabel="error / physical limit", title="B  Close to the limit")
fig.tight_layout()
save(fig, "spot_precision.png")
print(dict(zip(photons.tolist(), np.round(np.array(rms) * 100, 1))), found)
