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
def _sample(A, kern, logit_p, sb2, beta, s, J0, n_burn, n_keep, thin, max_shift, seed):
    np.random.seed(seed)
    N = len(logit_p)
    T = len(A)
    idx = J0 + np.arange(T) * s
    I = np.zeros(N, np.int8)
    G = np.zeros(T)
    counts = np.zeros((n_keep, T), np.uint8)
    Gsum = np.zeros(T)
    kept = 0
    total = n_burn + n_keep * thin
    for sweep in range(total):
        # single-step Gibbs: add or remove an event at every fine step
        for j in range(N):
            # log-odds of "an event at j" vs "none", given everything else: likelihood difference + prior log-odds
            if I[j] == 0:
                lo = _delta_ll(A, G, kern, idx, j, 1.0, sb2, beta, s, J0) + logit_p[j]
                if np.random.random() < 1.0 / (1.0 + np.exp(-lo)):
                    I[j] = 1
                    _apply(G, kern, idx, j, 1.0, s, J0)
            else:
                lo = -_delta_ll(A, G, kern, idx, j, -1.0, sb2, beta, s, J0) + logit_p[j]
                if np.random.random() >= 1.0 / (1.0 + np.exp(-lo)):
                    I[j] = 0
                    _apply(G, kern, idx, j, -1.0, s, J0)
        # shift moves: move every event by a few steps (Metropolis), so an event can slide in time without first disappearing
        for j in range(N):
            if I[j] == 0:
                continue
            d = np.random.randint(1, max_shift + 1) * (1 if np.random.random() < 0.5 else -1)
            k = j + d
            if k < 0 or k >= N or I[k] == 1:
                continue
            dl = _delta_ll(A, G, kern, idx, j, -1.0, sb2, beta, s, J0)
            _apply(G, kern, idx, j, -1.0, s, J0)
            dl += _delta_ll(A, G, kern, idx, k, 1.0, sb2, beta, s, J0)
            dl += logit_p[k] - logit_p[j]
            if np.log(np.random.random()) < dl:
                _apply(G, kern, idx, k, 1.0, s, J0)
                I[j] = 0
                I[k] = 1
            else:
                _apply(G, kern, idx, j, 1.0, s, J0)
        if sweep >= n_burn and (sweep - n_burn) % thin == 0 and kept < n_keep:
            for n in range(T):          # events in (t_n - dt, t_n]: fine steps idx[n]-s+1 .. idx[n]
                c = 0
                for q in range(idx[n] - s + 1, idx[n] + 1):
                    if 0 <= q < N:
                        c += I[q]
                counts[kept, n] = min(c, 255)
            Gsum += G
            kept += 1
    return counts, Gsum / max(kept, 1)


def deconvolve_trace(a_cu: np.ndarray, kern: Kernel, frame_s: float, rate_prior_per_s: np.ndarray, sigma_b: float, beta: float,
                     n_burn: int = 500, n_keep: int = 1000, thin: int = 3, max_shift_s: float = 15.0, seed: int = 0) -> dict:
    """Posterior over the initiation events of one trace (C.U., one value per frame). rate_prior_per_s: prior initiation rate per frame
    (per second). Returns counts (n_keep, T) of initiations in each frame interval per sample and the posterior mean fit G (C.U.)."""
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
    counts, G = _sample(np.asarray(a_cu, np.float64), kern.values, np.log(p / (1 - p)), float(sigma_b) ** 2, float(beta), s, J0,
                        int(n_burn), int(n_keep), int(thin), max(1, int(round(max_shift_s / kern.dt_fine_s))), int(seed))
    return {"counts": counts, "G": G}


def burst_calls(counts: np.ndarray, frame_s: float, w_frames: float = 5.0, g_b: float = 2.0) -> np.ndarray:
    """Chen et al. V.4 on every sample: smooth the per-frame initiation counts with exp(-(2 t / w)^4) (normalised) and call ON where the
    smoothed count per frame >= g_b / w. Returns P(ON) per frame (fraction of samples)."""
    h = int(np.ceil(w_frames))
    t = np.arange(-h, h + 1, dtype=float)
    k = np.exp(-(2 * t / w_frames) ** 4)
    k /= k.sum()
    sm = np.array([np.convolve(c.astype(float), k, mode="same") for c in counts])
    return (sm >= g_b / w_frames).mean(axis=0)


def deconvolve_run_traces(tr: pd.DataFrame, cfg: dict, col: str = "intensity_corr", unit: float | None = None, n_keep: int = 1000,
                          seed: int = 0, progress=None) -> tuple[pd.DataFrame, dict, dict]:
    """Deconvolve every cell's MS2 trace (table with cell_id, t, and the intensity column). Returns (per-frame table, per-cell samples
    {cell_id: counts}, info: unit, noise model, kernel)."""
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
    rows, samples = [], {}
    for i, (cid, a) in enumerate(traces.items()):
        res = deconvolve_trace(a, kern, dt, prior, sb, beta, n_keep=n_keep, seed=seed + cid)
        c = res["counts"]
        samples[cid] = c
        p_on = burst_calls(c, dt)
        for n in range(len(a)):
            rows.append((cid, n, a[n], res["G"][n], c[:, n].mean(), c[:, n].std(), 60.0 * c[:, n].mean() / dt, p_on[n]))
        if progress:
            progress(f"  cell {cid}: {c.sum(1).mean():.1f} initiations (posterior mean) [{i + 1}/{len(traces)}]")
    out = pd.DataFrame(rows, columns=["cell_id", "t", "trace_cu", "fit_cu", "n_init_mean", "n_init_sd", "rate_per_min", "p_on"])
    info = {"unit": unit, "sigma_b_cu": sb, "beta_cu": beta, "kernel_dt_s": kern.dt_fine_s, "kernel_s": len(kern.values) * kern.dt_fine_s,
            "kernel_area_s": tau, "column": col, "n_samples": int(n_keep)}
    return out, samples, info


def deconvolve_run(run, sources: tuple = ("truth",), progress=print, n_keep: int = 1000) -> dict:
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
        dec, samples, info = deconvolve_run_traces(pd.read_csv(f), cfg, n_keep=n_keep, progress=progress)
        info["trace_file"] = f.name
        dec.to_csv(d6 / "ms2_deconvolution.csv", index=False)
        np.savez_compressed(d6 / "ms2_deconvolution_samples.npz", **{str(c): v for c, v in samples.items()})
        if src == "truth":
            ev = pd.read_csv(run / "stage1_chromatin" / "transcription_events.csv")
            tx = pd.read_csv(run / "stage1_chromatin" / "transcription.csv")
            info["scores"] = out = score_deconvolution(dec, samples, ev, tx, float(cfg["acquisition"]["frame_interval_s"]))
            # how far the microscope-based unit is from the truth (regression of the trace on the true loops / n_loops, through zero)
            n_loops = int((cfg["loci"]["transcription"].get("params") or {}).get("n_loops", 24))
            j = dec.merge(tx[["cell_id", "t", "ms2_loops"]], on=["cell_id", "t"])
            g = j.ms2_loops.to_numpy() / n_loops
            ok = np.isfinite(j.trace_cu)
            if (g[ok] ** 2).sum() > 0:
                info["unit_check"] = float((j.trace_cu[ok] * g[ok]).sum() / (g[ok] ** 2).sum())   # 1 = unit right; 0.8 = trace 20% too low
        (d6 / "ms2_deconvolution_info.json").write_text(json.dumps(info, indent=1, default=float))
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


def score_deconvolution(dec: pd.DataFrame, samples: dict, events: pd.DataFrame, tx: pd.DataFrame, frame_s: float,
                        tol_frames: int = 3) -> dict:
    """Against the simulation: total initiations (posterior mean vs true), events found (a true initiation with a called event within
    +-tol_frames: called = frames where the posterior expects >= 0.5 initiations, merged with their neighbours), false events, the
    correlation of the 1-minute rate with the true one, and P(ON) against the true promoter state."""
    T = int(dec.t.max()) + 1
    tru = true_initiations(events, sorted(dec.cell_id.unique()), T, frame_s)
    j = dec.merge(tru, on=["cell_id", "t"]).merge(tx[["cell_id", "t", "promoter_on"]], on=["cell_id", "t"], how="left")
    found = called = false = n_true = 0
    for cid, g in j.groupby("cell_id"):
        g = g.sort_values("t")
        exp = g.n_init_mean.to_numpy()
        # called events: local groups of expected count; each group of consecutive frames with expected count >= 0.15 holding >= 0.5 events
        ev = []
        m = exp >= 0.15
        i = 0
        while i < len(m):
            if m[i]:
                k = i
                while k + 1 < len(m) and m[k + 1]:
                    k += 1
                mass = exp[i:k + 1].sum()
                for _ in range(int(np.floor(mass + 0.5))):
                    ev.append(i + int(np.argmax(exp[i:k + 1])))
                i = k + 1
            else:
                i += 1
        te = np.repeat(np.arange(len(g)), g.n_init_true.to_numpy())
        used = np.zeros(len(ev), bool)
        for t0 in te:
            d = [abs(t0 - e) if not used[q] else 1e9 for q, e in enumerate(ev)]
            if d and min(d) <= tol_frames:
                used[int(np.argmin(d))] = True
                found += 1
        n_true += len(te)
        called += len(ev)
        false += int((~used).sum())
    k = max(1, int(round(60 / frame_s)))
    sm = lambda v: pd.Series(v).rolling(k, center=True, min_periods=1).sum().to_numpy()  # noqa: E731
    r_est = np.concatenate([sm(g.sort_values("t").n_init_mean.to_numpy()) for _, g in j.groupby("cell_id")])
    r_true = np.concatenate([sm(g.sort_values("t").n_init_true.to_numpy()) for _, g in j.groupby("cell_id")])
    out = {"n_true": int(n_true), "n_estimated": float(j.n_init_mean.sum()), "n_called": int(called), "found": found / max(n_true, 1),
           "false_per_true": false / max(n_true, 1), "rate_r_1min": float(np.corrcoef(r_est, r_true)[0, 1]) if r_true.std() > 0 and r_est.std() > 0 else np.nan}
    if "promoter_on" in j and j.promoter_on.notna().any():
        on = j.promoter_on.fillna(0).astype(bool).to_numpy()
        call = j.p_on.to_numpy() >= 0.5
        out.update(p_on_accuracy=float((on == call).mean()), on_frames_true=float(on.mean()), on_frames_called=float(call.mean()),
                   on_sensitivity=float(call[on].mean()) if on.any() else np.nan)
    return out
