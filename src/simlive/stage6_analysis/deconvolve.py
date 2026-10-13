"""Stage 6: from an MS2 trace to the moments Pol II started (initiation events) and the bursts, after Chen, Levo, Zoller & Gregor 2025
(Nat. Struct. Mol. Biol.; methods in arXiv 2304.08770v3, sections IV.3-IV.4, V.1 and V.4).

1. Units. The trace (background-subtracted, bleaching/illumination-corrected spot intensity at the tracked promoter, one value per frame)
   is divided by the intensity of ONE Pol II carrying all its MS2 loops ("C.U.", in the paper: one fully tagged transcript, calibrated
   against smFISH). Here the unit comes from the microscope: MS2 photons per loop x loops x exposure x quantum efficiency x gain x the
   fraction of the spot inside the measuring ellipse / the ellipse's pixel count, i.e. what a single-transcript calibration would give.
2. The model (their Eqs. 8-11):  A(t) = G(t) + sigma(G) eta(t),  G = kernel * I,  I = sum of initiation events, and
       kernel(u) = (1/m) sum_i H(u - l_i / v)   for u < L / v,   else 0
   = the fraction of the m loops a Pol II has made u seconds after it started: m steps while it copies the cassette (loop i at position
   l_i), then 1 while it travels the rest of the gene and waits to be released (L = cassette + gene + v x mean dwell, as the paper adds
   ~800 bp for retention), then 0. Assumptions as in the paper: constant elongation speed v, no pausing, fast release.
3. Noise (their IV.4): sigma^2(G) = sigma_b^2 + beta G, estimated from the traces themselves. Imaging noise is white and the signal is
   smooth on the frame scale, so half the mean squared frame-to-frame difference at a given level is the noise variance at that level.
4. Prior (their Eq. 15): independent Bernoulli on a fine time grid (their dt' = Pol II footprint / speed ~1 s: the shortest possible
   interval between two initiations), p_i = R(t_i) dt' with R(t) = population mean activity(t + tau/2) / tau, tau = the kernel's area.
5. Posterior sampling of the initiation configuration. Their sampler is an adaptive block-Bernoulli Metropolis-Hastings with a floor of
   0.05 on the proposal probabilities: right for fly genes firing ~10-20 times a minute, but it cannot propose the sparse initiation of
   the placeholder gene here (~1 per 10 min ON; a 6-minute block would propose ~18 events). The SAME posterior is sampled instead with
   single-step Gibbs updates (add / remove an event at each fine step) plus moves that shift an event in time. Burn-in 500 sweeps, then
   1000 samples (every 3rd sweep), as in the paper (N_s = 3500, N_b = 500).
6. Rate and bursts (their Eq. 18 and V.4): r(t) = initiations in (t - dt, t] / dt for every sample (mean and sd over samples); bursts:
   r smoothed with exp(-(2 t / w)^4), w = 5 frames, ON where it is at least g_b / w (g_b = 2 initiations within the window), applied to
   every sample -> P(ON). With the placeholder gene's sparse initiation the default rule rarely calls ON: g_b and w can be changed.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from numba import njit


@dataclass
class Kernel:
    dt_fine_s: float
    values: np.ndarray            # kernel at u = d x dt_fine for d = 0..len-1 (C.U.)

    @property
    def area_s(self) -> float:
        return float(self.values.sum() * self.dt_fine_s)


def make_kernel(n_loops: int, cassette_kb: float, gene_kb: float, elongation_kb_min: float, dwell_s: float, dt_fine_s: float = 1.0) -> Kernel:
    """Step kernel of Chen et al. Eq. 11: loop i (at l_i = i x cassette / n_loops) adds 1 / n_loops once Pol II has passed it; the signal
    ends when Pol II is released, at (cassette + gene) / v + dwell."""
    v = elongation_kb_min / 60.0                                   # kb per second
    t_end = (cassette_kb + gene_kb) / v + dwell_s
    u = np.arange(int(np.ceil(t_end / dt_fine_s)) + 1) * dt_fine_s
    loops = np.arange(1, n_loops + 1) * cassette_kb / n_loops / v  # time at which each loop has been made
    k = (u[:, None] >= loops[None, :]).sum(1) / n_loops
    k[u >= t_end] = 0.0
    return Kernel(dt_fine_s, k.astype(np.float64))


def kernel_from_config(cfg: dict, dt_fine_s: float = 1.0) -> Kernel:
    p = {**(cfg["loci"]["transcription"].get("params") or {})}
    from simlive.stage1_chromatin.transcription import DEFAULTS

    g = lambda k: float(p.get(k, DEFAULTS[k]))  # noqa: E731
    return make_kernel(int(g("n_loops")), g("cassette_kb"), g("gene_kb"), g("elongation_kb_min"), g("dwell_s"), dt_fine_s)


def unit_intensity(cfg: dict, inner_px=(7, 7)) -> float:
    """Ellipse intensity (camera counts per pixel, background removed) of ONE Pol II carrying all its loops, from the microscope: MS2 photons
    x QE x gain x (fraction of a Gaussian PSF inside the inner ellipse, after the haze fraction) / number of pixels in the ellipse."""
    from simlive.stage3_microscopy.render import psf_sigma_px

    opt, acq, tc = cfg["optics"], cfg["acquisition"], cfg["loci"]["transcription"]
    n_loops = int((tc.get("params") or {}).get("n_loops", 24))
    photons = float(tc.get("photons_per_loop_s", 100.0)) * n_loops * float(acq["exposure_s"])
    counts = photons * float(opt.get("quantum_efficiency", 0.8)) * float(opt.get("gain_adu_per_e", 1.0))
    s = psf_sigma_px(opt)
    haze = float(((cfg.get("imaging_errors") or {}).get("haze") or {}).get("fraction", 0.0)) \
        if (cfg.get("imaging_errors") or {}).get("enabled") and ((cfg.get("imaging_errors") or {}).get("haze") or {}).get("enabled", False) else 0.0
    # pixel-integrated Gaussian centred on a pixel, summed over the pixels of the ellipse
    r = np.arange(-6, 7)
    from scipy.special import erf

    w = 0.5 * (erf((r + 0.5) / (np.sqrt(2) * s)) - erf((r - 0.5) / (np.sqrt(2) * s)))
    yy, xx = np.meshgrid(r, r, indexing="ij")
    inside = (yy / (inner_px[1] / 2)) ** 2 + (xx / (inner_px[0] / 2)) ** 2 <= 1.0
    frac = float(np.outer(w, w)[inside].sum())
    return counts * (1 - haze) * frac / float(inside.sum())


def noise_points(traces: list[np.ndarray], n_bins: int = 8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Half the squared frame-to-frame difference, binned by signal level (robust: median / 0.455, the median of a chi-square with one degree
    of freedom). Returns (level, variance, frames per bin)."""
    d2, lev = [], []
    for a in traces:
        a = np.asarray(a, float)
        ok = np.isfinite(a[1:]) & np.isfinite(a[:-1])
        d2.append(0.5 * (a[1:] - a[:-1])[ok] ** 2)
        lev.append(0.5 * (a[1:] + a[:-1])[ok])
    d2, lev = (np.concatenate(d2), np.concatenate(lev)) if d2 else (np.zeros(0), np.zeros(0))
    if len(d2) < 20:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    edges = np.unique(np.quantile(lev, np.linspace(0, 1, n_bins + 1)))
    x, y, w = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (lev >= a) & (lev <= b)
        if m.sum() >= 10:
            x.append(np.median(lev[m])); y.append(np.median(d2[m]) / 0.455); w.append(m.sum())
    return np.array(x), np.array(y), np.array(w, float)


def noise_model(traces: list[np.ndarray], n_bins: int = 8) -> tuple[float, float]:
    """sigma^2(G) = sigma_b^2 + beta G from the traces (C.U.), fitted by weighted least squares to `noise_points`. Returns (sigma_b, beta)."""
    x, y, w = noise_points(traces, n_bins)
    if len(x) == 0:
        return 0.3, 0.0
    if len(x) < 2:
        return float(np.sqrt(y[0])), 0.0
    A = np.c_[np.ones_like(x), np.clip(x, 0, None)] * np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)
    s2, beta = float(max(coef[0], 1e-4)), float(max(coef[1], 0.0))
    return float(np.sqrt(s2)), beta


@njit(cache=True)
def _delta_ll(A, G, kern, idx, j, sign, sb2, beta, s, J0):
    """Change of the log-likelihood if the event at fine step j is added (sign=+1) or removed (-1). Frames n see it if 0 < idx[n] - j < len."""
    K = len(kern)
    T = len(A)
    n_lo = max(0, (j + 1 - J0 + s - 1) // s)
    n_hi = min(T - 1, (j + K - 1 - J0) // s)
    d = 0.0
    for n in range(n_lo, n_hi + 1):
        if A[n] != A[n]:
            continue
        u = idx[n] - j
        if u <= 0 or u >= K:
            continue
        g0 = G[n]
        g1 = g0 + sign * kern[u]
        v0 = sb2 + beta * max(g0, 0.0)
        v1 = sb2 + beta * max(g1, 0.0)
        d += (-0.5 * (A[n] - g1) ** 2 / v1 - 0.5 * np.log(v1)) - (-0.5 * (A[n] - g0) ** 2 / v0 - 0.5 * np.log(v0))
    return d


@njit(cache=True)
def _apply(G, kern, idx, j, sign, s, J0):
    K = len(kern)
    T = len(G)
    n_lo = max(0, (j + 1 - J0 + s - 1) // s)
    n_hi = min(T - 1, (j + K - 1 - J0) // s)
    for n in range(n_lo, n_hi + 1):
        u = idx[n] - j
        if 0 < u < K:
            G[n] += sign * kern[u]


@njit(cache=True)
def _log_post(A, G, I, logit_p, sb2, beta):
    """Log-posterior of the current configuration up to a constant: Gaussian log-likelihood of the trace + the prior log-odds of the starts."""
    lp = 0.0
    for n in range(len(A)):
        if A[n] == A[n]:
            v = sb2 + beta * max(G[n], 0.0)
            lp += -0.5 * (A[n] - G[n]) ** 2 / v - 0.5 * np.log(v)
    for j in range(len(I)):
        if I[j]:
            lp += logit_p[j]
    return lp


@njit(cache=True)
def _sweep(A, G, I, kern, idx, logit_p, sb2, beta, s, J0, max_shift, invT):
    """One sweep at inverse temperature invT (the likelihood is raised to the power invT, the prior is not): single-step Gibbs updates
    (add / remove an event at every fine step), then a Metropolis shift of every event by up to max_shift steps (max_shift < 1 = none), so an
    event can slide in time without first disappearing. invT = 1 is the posterior."""
    N = len(logit_p)
    for j in range(N):
        # log-odds of "an event at j" vs "none", given everything else: likelihood difference + prior log-odds
        if I[j] == 0:
            lo = invT * _delta_ll(A, G, kern, idx, j, 1.0, sb2, beta, s, J0) + logit_p[j]
            if np.random.random() < 1.0 / (1.0 + np.exp(-lo)):
                I[j] = 1
                _apply(G, kern, idx, j, 1.0, s, J0)
        else:
            lo = -invT * _delta_ll(A, G, kern, idx, j, -1.0, sb2, beta, s, J0) + logit_p[j]
            if np.random.random() >= 1.0 / (1.0 + np.exp(-lo)):
                I[j] = 0
                _apply(G, kern, idx, j, -1.0, s, J0)
    for j in range(N):
        if max_shift < 1 or I[j] == 0:
            continue
        d = np.random.randint(1, max_shift + 1) * (1 if np.random.random() < 0.5 else -1)
        k = j + d
        if k < 0 or k >= N or I[k] == 1:
            continue
        dl = _delta_ll(A, G, kern, idx, j, -1.0, sb2, beta, s, J0)
        _apply(G, kern, idx, j, -1.0, s, J0)
        dl = invT * (dl + _delta_ll(A, G, kern, idx, k, 1.0, sb2, beta, s, J0)) + logit_p[k] - logit_p[j]
        if np.log(np.random.random()) < dl:
            _apply(G, kern, idx, k, 1.0, s, J0)
            I[j] = 0
            I[k] = 1
        else:
            _apply(G, kern, idx, j, 1.0, s, J0)


@njit(cache=True)
def _loglik(A, G, sb2, beta):
    ll = 0.0
    for n in range(len(A)):
        if A[n] == A[n]:
            v = sb2 + beta * max(G[n], 0.0)
            ll += -0.5 * (A[n] - G[n]) ** 2 / v - 0.5 * np.log(v)
    return ll


@njit(cache=True)
def _seed(seed):
    np.random.seed(seed)


@njit(cache=True)
def _pt_sweep(A, kern, logit_p, sb2, beta, s, J0, I2, G2, tix, invT, max_shift, parity, swap_try, swap_acc):
    """Parallel tempering, one round: every replica (row of I2 / G2) does a sweep at the inverse temperature of its rank tix[r]; then
    neighbouring ranks (parity 0: 0-1, 2-3 ...; parity 1: 1-2, 3-4 ...) swap their states with probability min(1, exp((invT_a - invT_b) (ll_b - ll_a))),
    which keeps the product of the tempered posteriors invariant. Returns the row holding the rank-0 (invT = 1) state."""
    R = I2.shape[0]
    idx = J0 + np.arange(len(A)) * s
    for r in range(R):
        _sweep(A, G2[r], I2[r], kern, idx, logit_p, sb2, beta, s, J0, max_shift, invT[tix[r]])
    row_of_rank = np.zeros(R, np.int64)
    for r in range(R):
        row_of_rank[tix[r]] = r
    ll = np.zeros(R)
    for r in range(R):
        ll[r] = _loglik(A, G2[r], sb2, beta)
    for k in range(parity, R - 1, 2):
        ra, rb = row_of_rank[k], row_of_rank[k + 1]
        swap_try[k] += 1
        if np.log(np.random.random()) < (invT[k] - invT[k + 1]) * (ll[rb] - ll[ra]):
            tix[ra], tix[rb] = k + 1, k
            row_of_rank[k], row_of_rank[k + 1] = rb, ra
            swap_acc[k] += 1
    return row_of_rank[0]


def _sample_pt(A, kern, logit_p, sb2, beta, s, J0, n_burn, n_keep, thin, max_shift, seed, I0, n_replicas=16, beta_min=1e-3, target=0.3):
    """Parallel tempering of the Gibbs + shift sampler (the likelihood at inverse temperatures 1 = invT[0] > ... > beta_min; the prior is
    never tempered). Hot replicas move freely between arrangements of the starts and hand them down the ladder, which is what a rugged
    posterior (a nearly noise-free trace: moving one start makes the fit much worse) needs. The ladder is adapted during the first half of the
    burn-in (Robbins-Monro on the log-gaps towards `target` swap acceptance, ends fixed) and then frozen. Samples come from the
    invT = 1 state only. Returns the same as `_sample` plus the mean swap acceptance (last element of the "acceptance" slot) and the ladder."""
    N, T = len(logit_p), len(A)
    R = int(n_replicas)
    _seed(int(seed))
    I2 = np.tile(I0.astype(np.int8), (R, 1))
    G2 = np.zeros((R, T))
    idx = J0 + np.arange(T) * s
    for r in range(R):
        for j in np.flatnonzero(I2[r]):
            _apply(G2[r], kern, idx, int(j), 1.0, s, J0)
    tix = np.arange(R, dtype=np.int64)
    gaps = np.full(R - 1, np.log(1.0 / beta_min) / (R - 1))
    invT = np.exp(-np.r_[0.0, np.cumsum(gaps)])
    swap_try, swap_acc = np.zeros(R - 1), np.zeros(R - 1)
    total = n_burn + n_keep * thin
    adapt_until = n_burn // 2
    counts = np.zeros((n_keep, T), np.uint8)
    lp_kept, nst_kept, lp_sweep = np.zeros(n_keep), np.zeros(n_keep), np.zeros(total)
    Gsum, Isum, kept, rnd = np.zeros(T), np.zeros(N), 0, 0
    hi, lo_ = np.clip(idx + 1, 0, N), np.clip(idx - s + 1, 0, N)
    for sweep in range(total):
        c = _pt_sweep(A, kern, logit_p, sb2, beta, s, J0, I2, G2, tix, invT, max_shift, sweep % 2, swap_try, swap_acc)
        I, G = I2[c], G2[c]
        lp_sweep[sweep] = _log_post(A, G, I, logit_p, sb2, beta)
        if sweep < adapt_until and (sweep + 1) % 25 == 0:           # ladder adaptation: widen gaps that swap too often, narrow the others
            rnd += 1
            acc = swap_acc / np.maximum(swap_try, 1)
            gaps = gaps * np.exp((acc - target) / np.sqrt(rnd))
            gaps = np.clip(gaps, 1e-3, None) * (np.log(1.0 / beta_min) / np.clip(gaps, 1e-3, None).sum())
            invT = np.exp(-np.r_[0.0, np.cumsum(gaps)])
            swap_try[:], swap_acc[:] = 0, 0
        if sweep == max(adapt_until, 1) - 1:
            swap_try[:], swap_acc[:] = 0, 0                          # report the swap rates of the frozen ladder only
        if sweep >= n_burn and (sweep - n_burn) % thin == 0 and kept < n_keep:
            cs = np.r_[0, np.cumsum(I, dtype=np.int64)]
            counts[kept] = np.minimum(cs[hi] - cs[lo_], 255)
            Gsum += G
            Isum += I
            lp_kept[kept] = lp_sweep[sweep]
            nst_kept[kept] = I.sum()
            kept += 1
    swap = swap_acc / np.maximum(swap_try, 1)
    return counts, Gsum / max(kept, 1), Isum / max(kept, 1), float(swap.mean()), lp_kept, nst_kept, lp_sweep, invT, swap


@njit(cache=True)
def _sample(A, kern, logit_p, sb2, beta, s, J0, n_burn, n_keep, thin, max_shift, seed, I0):
    """Gibbs flips + shift moves from the start configuration I0 (max_shift < 1 = flips only). Besides the samples it returns the
    log-posterior and the number of starts of every kept sample and the log-posterior after every sweep (burn-in included)."""
    np.random.seed(seed)
    N = len(logit_p)
    T = len(A)
    idx = J0 + np.arange(T) * s
    I = I0.copy()
    G = np.zeros(T)
    for j in range(N):
        if I[j]:
            _apply(G, kern, idx, j, 1.0, s, J0)
    lp_kept = np.zeros(n_keep)
    nst_kept = np.zeros(n_keep)
    counts = np.zeros((n_keep, T), np.uint8)
    Gsum = np.zeros(T)
    Isum = np.zeros(N)
    kept = 0
    total = n_burn + n_keep * thin
    lp_sweep = np.zeros(total)
    for sweep in range(total):
        _sweep(A, G, I, kern, idx, logit_p, sb2, beta, s, J0, max_shift, 1.0)
        lp_sweep[sweep] = _log_post(A, G, I, logit_p, sb2, beta)
        if sweep >= n_burn and (sweep - n_burn) % thin == 0 and kept < n_keep:
            for n in range(T):          # events in (t_n - dt, t_n]: fine steps idx[n]-s+1 .. idx[n]
                c = 0
                for q in range(idx[n] - s + 1, idx[n] + 1):
                    if 0 <= q < N:
                        c += I[q]
                counts[kept, n] = min(c, 255)
            Gsum += G
            Isum += I
            lp_kept[kept] = lp_sweep[sweep]
            nst_kept[kept] = I.sum()
            kept += 1
    return counts, Gsum / max(kept, 1), Isum / max(kept, 1), 1.0, lp_kept, nst_kept, lp_sweep


@njit(cache=True)
def _sample_block(A, kern, logit_p, sb2, beta, s, J0, n_burn_sweeps, n_keep, w, q0, seed, I0):
    """Chen et al. 2025 (V.1, Eqs. 16-17): Metropolis-Hastings in which a random block of w consecutive steps is redrawn from independent
    Bernoulli(q) proposals; q learns from the chain, q <- q (1 - l) + l I with l = 1 / sqrt(3 j), kept within [0.05, 0.95]. Burn-in
    n_burn_sweeps x N / w iterations, then n_keep samples ~3 sweeps apart. Starts from I0; also returns the log-posterior and the number of
    starts of the kept samples and the log-posterior after every sweep (N / w iterations)."""
    np.random.seed(seed)
    N = len(logit_p)
    T = len(A)
    K = len(kern)
    idx = J0 + np.arange(T) * s
    I = I0.copy()
    G = np.zeros(T)
    for j in range(N):
        if I[j]:
            _apply(G, kern, idx, j, 1.0, s, J0)
    lp_kept = np.zeros(n_keep)
    nst_kept = np.zeros(n_keep)
    q = q0.copy()
    lp1 = np.log(1.0 / (1.0 + np.exp(-logit_p)))          # log p
    lp0 = np.log(1.0 / (1.0 + np.exp(logit_p)))           # log (1 - p)
    counts = np.zeros((n_keep, T), np.uint8)
    Gsum = np.zeros(T)
    Isum = np.zeros(N)
    iters_burn = max(1, n_burn_sweeps * N // w)
    per_keep = max(1, 3 * N // w)
    total = iters_burn + n_keep * per_keep
    per_sweep = max(1, N // w)
    lp_sweep = np.zeros(total // per_sweep + 1)
    kept = 0
    acc = 0
    Inew = np.zeros(w, np.int8)
    for it in range(total):
        b = np.random.randint(0, max(1, N - w + 1))
        e = min(N, b + w)
        for j in range(b, e):
            Inew[j - b] = 1 if np.random.random() < q[j] else 0
        n_lo = max(0, (b + 1 - J0 + s - 1) // s)
        n_hi = min(T - 1, (e - 1 + K - 1 - J0) // s)
        ll_old = 0.0
        for n in range(n_lo, n_hi + 1):
            if A[n] == A[n]:
                v = sb2 + beta * max(G[n], 0.0)
                ll_old += -0.5 * (A[n] - G[n]) ** 2 / v - 0.5 * np.log(v)
        prior = 0.0
        logq_ratio = 0.0
        for j in range(b, e):
            if Inew[j - b] != I[j]:
                _apply(G, kern, idx, j, float(Inew[j - b] - I[j]), s, J0)
            prior += (lp1[j] if Inew[j - b] else lp0[j]) - (lp1[j] if I[j] else lp0[j])
            qj = q[j]
            logq_ratio += (np.log(qj) if I[j] else np.log(1 - qj)) - (np.log(qj) if Inew[j - b] else np.log(1 - qj))
        ll_new = 0.0
        for n in range(n_lo, n_hi + 1):
            if A[n] == A[n]:
                v = sb2 + beta * max(G[n], 0.0)
                ll_new += -0.5 * (A[n] - G[n]) ** 2 / v - 0.5 * np.log(v)
        if np.log(np.random.random()) < ll_new - ll_old + prior + logq_ratio:
            for j in range(b, e):
                I[j] = Inew[j - b]
            acc += 1
        else:                                            # undo the change to G
            for j in range(b, e):
                if Inew[j - b] != I[j]:
                    _apply(G, kern, idx, j, float(I[j] - Inew[j - b]), s, J0)
        lr = 1.0 / np.sqrt(3.0 * (it + 1))
        for j in range(b, e):
            q[j] = min(0.95, max(0.05, q[j] * (1 - lr) + lr * I[j]))
        if (it + 1) % per_sweep == 0:
            lp_sweep[(it + 1) // per_sweep - 1] = _log_post(A, G, I, logit_p, sb2, beta)
        if it >= iters_burn and (it - iters_burn) % per_keep == 0 and kept < n_keep:
            for n in range(T):
                c = 0
                for qq in range(idx[n] - s + 1, idx[n] + 1):
                    if 0 <= qq < N:
                        c += I[qq]
                counts[kept, n] = min(c, 255)
            Gsum += G
            Isum += I
            lp_kept[kept] = _log_post(A, G, I, logit_p, sb2, beta)
            nst_kept[kept] = I.sum()
            kept += 1
    return counts, Gsum / max(kept, 1), Isum / max(kept, 1), acc / max(total, 1), lp_kept, nst_kept, lp_sweep[: total // per_sweep]


SAMPLERS = ("gibbs_shift", "chen_block", "gibbs", "gibbs_pt")   # "gibbs" = flips only (no shifts; for the mixing comparison); "gibbs_pt" = gibbs_shift with parallel tempering


def deconvolve_trace(a_cu: np.ndarray, kern: Kernel, frame_s: float, rate_prior_per_s: np.ndarray, sigma_b: float, beta: float,
                     n_burn: int = 500, n_keep: int = 1000, thin: int = 3, max_shift_s: float = 15.0, seed: int = 0,
                     sampler: str = "gibbs_shift", block_s="auto", init: str = "empty", n_replicas: int = 16, beta_min: float = 1e-3) -> dict:
    """Posterior over the initiation events of one trace (C.U., one value per frame). rate_prior_per_s: prior initiation rate per frame
    (per second). sampler: "gibbs_shift" (default: single-step Gibbs + shift moves) or "chen_block" (the paper's adaptive block-Bernoulli
    Metropolis-Hastings; block_s "auto" = the largest of 360 ... 5 s with >= 10% acceptance in a pilot run, as the paper tuned it). "gibbs" = flips only; "gibbs_pt" = parallel tempering (n_replicas, beta_min; `acceptance` is then the mean swap rate, ladder and swap_rates are returned). init: starting configuration (INITS). Returns counts (n_keep, T) of initiations in each frame interval per
    sample, the posterior mean fit G (C.U.), p_fine = the probability of a start at every fine step, t_fine = their times (s, 0 = first
    frame), the acceptance rate, and for convergence checks log_post / n_starts per kept sample and log_post_sweeps (every sweep)."""
    s = int(round(frame_s / kern.dt_fine_s))
    if abs(s * kern.dt_fine_s - frame_s) > 1e-6:
        raise ValueError("the frame interval must be a whole number of fine steps")
    K = len(kern.values)
    J0 = K                                      # fine steps before the first frame (Pol II started before the movie still show)
    T = len(a_cu)
    N = J0 + (T - 1) * s + 1
    tf = (np.arange(N) - J0) * kern.dt_fine_s   # time of each fine step (s; 0 = first frame)
    r = np.interp(tf, np.arange(T) * frame_s, np.asarray(rate_prior_per_s, float))
    p = np.clip(r * kern.dt_fine_s, 1e-6, 0.5)
    lo = np.log(p / (1 - p))
    A = np.asarray(a_cu, np.float64)
    rise = np.clip(np.r_[np.diff(np.nan_to_num(A)), 0.0], 0, None)           # where the trace rises (the paper: smoothed activity)
    rise_f = np.interp(tf, np.arange(T) * frame_s, rise)
    I0 = initial_state(init, p, rise_f, kern, np.random.default_rng(int(seed) + 7919))
    if sampler == "chen_block":
        q0 = np.clip(rise_f, 0.05, 0.95)         # q starts high where the trace rises
        args = (A, kern.values, lo, float(sigma_b) ** 2, float(beta), s, J0)
        if block_s == "auto":            # the paper tuned the block for an acceptance rate of 10-70%: largest block with >= 10% in a pilot run
            for cand in (360.0, 180.0, 90.0, 45.0, 20.0, 10.0, 5.0):
                w = max(2, int(round(cand / kern.dt_fine_s)))
                if _sample_block(*args, 100, 20, w, q0, int(seed), I0)[3] >= 0.10:
                    break
            block_s = cand
        w = max(2, int(round(float(block_s) / kern.dt_fine_s)))
        counts, G, pf, acc, lp_k, nst_k, lp_sw = _sample_block(*args, int(n_burn), int(n_keep), w, q0, int(seed), I0)
    elif sampler in ("gibbs_shift", "gibbs"):
        shift = max(1, int(round(max_shift_s / kern.dt_fine_s))) if sampler == "gibbs_shift" else 0
        counts, G, pf, acc, lp_k, nst_k, lp_sw = _sample(A, kern.values, lo, float(sigma_b) ** 2, float(beta), s, J0, int(n_burn), int(n_keep),
                                                         int(thin), shift, int(seed), I0)
    elif sampler == "gibbs_pt":
        counts, G, pf, acc, lp_k, nst_k, lp_sw, ladder, swap = _sample_pt(
            A, kern.values, lo, float(sigma_b) ** 2, float(beta), s, J0, int(n_burn), int(n_keep), int(thin),
            max(1, int(round(max_shift_s / kern.dt_fine_s))), int(seed), I0, int(n_replicas), float(beta_min))
        pt = {"ladder": ladder, "swap_rates": swap}
    else:
        raise ValueError(f"unknown sampler {sampler!r}; use one of {SAMPLERS}")
    return {"counts": counts, "G": G, "p_fine": pf, "t_fine": tf, "acceptance": float(acc), "block_s": block_s if sampler == "chen_block" else None,
            "log_post": lp_k, "n_starts": nst_k, "log_post_sweeps": lp_sw, **(pt if sampler == "gibbs_pt" else {})}


INITS = ("empty", "prior", "crowded", "rise")


def initial_state(kind: str, p: np.ndarray, rise: np.ndarray, kern: Kernel, rng: np.random.Generator) -> np.ndarray:
    """Starting configuration of a chain (0/1 per fine step). "empty": no starts; "prior": each step with its prior probability; "crowded": 10x
    the prior probability (a start every few dozen seconds where the prior is high, far too many); "rise": starts where the trace rises
    fastest (a data-driven guess). The chains of a convergence check start from different ones on purpose (over-dispersed)."""
    N = len(p)
    if kind == "empty":
        return np.zeros(N, np.int8)
    if kind == "prior":
        return (rng.random(N) < p).astype(np.int8)
    if kind == "crowded":
        return (rng.random(N) < np.clip(10 * p, 0, 0.5)).astype(np.int8)
    if kind == "rise":
        r = np.asarray(rise, float)
        return (rng.random(N) < np.clip(r / max(r.max(), 1e-9) * 0.2, 0, 0.5)).astype(np.int8)
    raise ValueError(f"unknown init {kind!r}; use one of {INITS}")


def burst_calls(counts: np.ndarray, frame_s: float, w_frames: float = 5.0, g_b: float = 2.0) -> np.ndarray:
    """Chen et al. V.4 on every sample: smooth the per-frame initiation counts with exp(-(2 t / w)^4) (normalised) and call ON where the
    smoothed count per frame >= g_b / w. Returns P(ON) per frame (fraction of samples)."""
    from scipy.ndimage import convolve1d

    h = int(np.ceil(w_frames))
    t = np.arange(-h, h + 1, dtype=float)
    k = np.exp(-(2 * t / w_frames) ** 4)
    k /= k.sum()
    sm = convolve1d(np.asarray(counts, float), k, axis=-1, mode="constant")
    return (sm >= g_b / w_frames - 1e-9).mean(axis=0)


def burst_rule_grid(samples: dict, tx: pd.DataFrame, frame_s: float, ws=(1, 2, 3, 5, 8, 12, 18, 30),
                    gbs=(0.5, 1, 1.5, 2, 3, 4, 6)) -> pd.DataFrame:
    """Agreement of the burst calls with the true promoter state for a grid of windows w (frames) and thresholds g_b (starts in the window):
    accuracy, sensitivity (ON frames called ON) and specificity (OFF frames called OFF), over all cells. Simulation only (needs the truth)."""
    truth = {int(c): g.sort_values("t").promoter_on.to_numpy().astype(bool) for c, g in tx.groupby("cell_id")}
    rows = []
    for w in ws:
        for gb in gbs:
            calls, tru = [], []
            for cid, s in samples.items():
                if cid not in truth:
                    continue
                c = burst_calls(s, frame_s, w, gb) >= 0.5
                n = min(len(c), len(truth[cid]))
                calls.append(c[:n]); tru.append(truth[cid][:n])
            c, t = np.concatenate(calls), np.concatenate(tru)
            rows.append({"w_frames": w, "g_b": gb, "accuracy": float((c == t).mean()), "sensitivity": float(c[t].mean()) if t.any() else np.nan,
                         "specificity": float((~c[~t]).mean()) if (~t).any() else np.nan})
    return pd.DataFrame(rows)


def prepare_traces(tr: pd.DataFrame, cfg: dict, col: str = "intensity_corr", unit: float | None = None):
    """Everything the samplers need from a run's MS2 table: per-cell traces in C.U. {cell_id: array}, the kernel, the frame interval, the
    prior rate per frame (per s), the noise model (sigma_b, beta) and the unit."""
    col = col if col in tr else "intensity"
    dt = float(cfg["acquisition"]["frame_interval_s"])
    kern = kernel_from_config(cfg, 1.0 if dt >= 1 else dt / 4)
    unit = float(unit) if unit else unit_intensity(cfg)
    tr = tr.sort_values(["cell_id", "t"])
    T = int(cfg["acquisition"]["n_frames"])
    traces = {int(c): g.set_index("t")[col].reindex(range(T)).to_numpy(float) / unit for c, g in tr.groupby("cell_id")}
    sb, beta = noise_model(list(traces.values()))
    pop = np.nanmean(np.vstack(list(traces.values())), axis=0)          # population mean activity (C.U.) per frame
    pop = pd.Series(np.clip(pop, 0, None)).rolling(5, center=True, min_periods=1).mean().to_numpy()
    tau = kern.area_s
    shift = int(round(tau / 2 / dt))
    prior = np.clip(np.r_[pop[shift:], np.full(shift, pop[-1] if len(pop) else 0)] / tau, 1e-5, 0.5)   # R(t) = mu(t + tau/2) / tau, per s
    return traces, kern, dt, prior, sb, beta, unit, col


def deconvolve_run_traces(tr: pd.DataFrame, cfg: dict, col: str = "intensity_corr", unit: float | None = None, n_keep: int = 1000,
                          seed: int = 0, progress=None, sampler: str = "gibbs_shift") -> tuple[pd.DataFrame, dict, dict]:
    """Deconvolve every cell's MS2 trace (table with cell_id, t, and the intensity column). Returns (per-frame table, per-cell samples
    {cell_id: counts}, info: unit, noise model, kernel)."""
    traces, kern, dt, prior, sb, beta, unit, col = prepare_traces(tr, cfg, col, unit)
    tau = kern.area_s
    rows, samples, acc, blocks = [], {}, [], []
    for i, (cid, a) in enumerate(traces.items()):
        res = deconvolve_trace(a, kern, dt, prior, sb, beta, n_keep=n_keep, seed=seed + cid, sampler=sampler)
        c = res["counts"]
        samples[cid] = c
        samples[f"fine_{cid}"] = res["p_fine"].astype(np.float32)          # probability of a start at every fine step
        acc.append(res["acceptance"])
        blocks.append(res["block_s"])
        p_on = burst_calls(c, dt)
        p_any = (c > 0).mean(axis=0)
        for n in range(len(a)):
            rows.append((cid, n, a[n], res["G"][n], c[:, n].mean(), c[:, n].std(), 60.0 * c[:, n].mean() / dt, p_on[n], p_any[n]))
        if progress:
            progress(f"  cell {cid}: {c.sum(1).mean():.1f} initiations (posterior mean) [{i + 1}/{len(traces)}]")
    out = pd.DataFrame(rows, columns=["cell_id", "t", "trace_cu", "fit_cu", "n_init_mean", "n_init_sd", "rate_per_min", "p_on", "p_start"])
    info = {"unit": unit, "sigma_b_cu": sb, "beta_cu": beta, "kernel_dt_s": kern.dt_fine_s, "kernel_s": len(kern.values) * kern.dt_fine_s,
            "kernel_area_s": tau, "column": col, "n_samples": int(n_keep), "sampler": sampler,
            "acceptance": float(np.mean(acc)) if acc else None, "block_s": [b for b in blocks if b is not None] or None,
            "fine_t0_s": -len(kern.values) * kern.dt_fine_s}          # time of the first fine step (starts before the movie are allowed)
    return out, samples, info


def deconvolve_run(run, sources: tuple = ("truth",), progress=print, n_keep: int = 1000, sampler: str = "gibbs_shift", tag: str = "") -> dict:
    """Deconvolve the MS2 traces of a finished run (the calibrated traces if there are any, bleaching-corrected column if present).
    Writes stage6_analysis/<src>/ms2_deconvolution.csv (+ _samples.npz, _info.json) and returns the scores of the truth source."""
    import json
    from pathlib import Path

    import yaml

    run = Path(run)
    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    if not (run / "stage1_chromatin" / "transcription.csv").exists():
        return {}
    out = {}
    for src in sources:
        d6 = run / "stage6_analysis" / src
        f = next((d6 / n for n in ("ms2_traces_cal.csv", "ms2_traces.csv") if (d6 / n).exists()), None)
        if f is None:
            continue
        progress(f"MS2 deconvolution ({src}, {f.name})")
        dec, samples, info = deconvolve_run_traces(pd.read_csv(f), cfg, n_keep=n_keep, progress=progress, sampler=sampler)
        info["trace_file"] = f.name
        dec.to_csv(d6 / f"ms2_deconvolution{tag}.csv", index=False)
        np.savez_compressed(d6 / f"ms2_deconvolution{tag}_samples.npz", **{str(c): v for c, v in samples.items()})
        if src == "truth":
            ev = pd.read_csv(run / "stage1_chromatin" / "transcription_events.csv")
            tx = pd.read_csv(run / "stage1_chromatin" / "transcription.csv")
            info["scores"] = out = score_deconvolution(dec, {k: v for k, v in samples.items() if isinstance(k, int)}, ev, tx,
                                                       float(cfg["acquisition"]["frame_interval_s"]))
            # how far the microscope-based unit is from the truth (regression of the trace on the true loops / n_loops, through zero)
            n_loops = int((cfg["loci"]["transcription"].get("params") or {}).get("n_loops", 24))
            j = dec.merge(tx[["cell_id", "t", "ms2_loops"]], on=["cell_id", "t"])
            g = j.ms2_loops.to_numpy() / n_loops
            ok = np.isfinite(j.trace_cu)
            if (g[ok] ** 2).sum() > 0:
                info["unit_check"] = float((j.trace_cu[ok] * g[ok]).sum() / (g[ok] ** 2).sum())   # 1 = unit right; 0.8 = trace 20% too low
        (d6 / f"ms2_deconvolution{tag}_info.json").write_text(json.dumps(info, indent=1, default=float))
    return out


def true_initiations(events: pd.DataFrame, cells: list[int], T: int, frame_s: float) -> pd.DataFrame:
    """Initiations of the simulation counted per frame interval (t - dt, t], per cell."""
    rows = []
    for c in cells:
        e = events[events.cell_id == c].init_s.to_numpy()
        n = np.ceil(e / frame_s).astype(int)
        cnt = np.bincount(n[(n >= 0) & (n < T)], minlength=T)
        rows += [(c, t, int(cnt[t])) for t in range(T)]
    return pd.DataFrame(rows, columns=["cell_id", "t", "n_init_true"])


def _match(est: np.ndarray, tru: np.ndarray, tol: int) -> int:
    """Number of one-to-one matches between two sorted lists of event frames within +-tol frames (greedy in time order)."""
    i = j = m = 0
    while i < len(est) and j < len(tru):
        d = est[i] - tru[j]
        if abs(d) <= tol:
            m += 1; i += 1; j += 1
        elif d < 0:
            i += 1
        else:
            j += 1
    return m


def score_deconvolution(dec: pd.DataFrame, samples: dict, events: pd.DataFrame, tx: pd.DataFrame, frame_s: float,
                        tol_frames: int = 3, n_eval: int = 200) -> dict:
    """Against the simulation, per posterior sample then averaged (so the uncertainty is respected and nothing is 'called'): every
    history's starts are matched one-to-one to the true starts within +-tol_frames. found = matched / true, false_per_true = unmatched
    inferred / true; also the total (posterior mean vs true), the correlation of the 1-minute rate with the true one, and P(ON) against
    the true promoter state. n_eval samples per cell are used."""
    T = int(dec.t.max()) + 1
    tru = true_initiations(events, sorted(dec.cell_id.unique()), T, frame_s)
    j = dec.merge(tru, on=["cell_id", "t"]).merge(tx[["cell_id", "t", "promoter_on"]], on=["cell_id", "t"], how="left")
    matched = inferred = 0.0
    n_true = 0
    for cid, g in j.groupby("cell_id"):
        g = g.sort_values("t")
        te = np.repeat(np.arange(len(g)), g.n_init_true.to_numpy())
        n_true += len(te)
        smp = samples.get(int(cid))
        if smp is None:
            continue
        pick = smp[np.linspace(0, len(smp) - 1, min(n_eval, len(smp))).astype(int)]
        for c in pick:
            ee = np.repeat(np.arange(len(c)), c.astype(int))
            matched += _match(ee, te, tol_frames) / len(pick)
            inferred += len(ee) / len(pick)
    k = max(1, int(round(60 / frame_s)))
    sm = lambda v: pd.Series(v).rolling(k, center=True, min_periods=1).sum().to_numpy()  # noqa: E731
    r_est = np.concatenate([sm(g.sort_values("t").n_init_mean.to_numpy()) for _, g in j.groupby("cell_id")])
    r_true = np.concatenate([sm(g.sort_values("t").n_init_true.to_numpy()) for _, g in j.groupby("cell_id")])
    out = {"n_true": int(n_true), "n_estimated": float(j.n_init_mean.sum()), "found": matched / max(n_true, 1),
           "false_per_true": (inferred - matched) / max(n_true, 1), "tolerance_s": tol_frames * frame_s, "rate_r_1min": float(np.corrcoef(r_est, r_true)[0, 1]) if r_true.std() > 0 and r_est.std() > 0 else np.nan}
    if "promoter_on" in j and j.promoter_on.notna().any():
        on = j.promoter_on.fillna(0).astype(bool).to_numpy()
        call = j.p_on.to_numpy() >= 0.5
        out.update(p_on_accuracy=float((on == call).mean()), on_frames_true=float(on.mean()), on_frames_called=float(call.mean()),
                   on_sensitivity=float(call[on].mean()) if on.any() else np.nan)
    return out
