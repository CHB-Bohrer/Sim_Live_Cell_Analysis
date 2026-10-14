"""P(burst ON) as a function of the E-P distance when the distance fluctuates faster than the promoter responds (docs/EP_CURVE.md).

With the polymer distances (calibrated to D = 0.002 um^2/s at 1 min) the E-P distance decorrelates in ~10 s, but a burst starts at rate k_on(d) and lasts ~45 s: the
promoter integrates the distance over minutes, so P(ON | d now) is nearly flat (contrast 1.4 for a distance that, held fixed, would give 6.7). What the data can show is the
dependence on the distance HISTORY. Model (linear response):

    P(ON at t) = a + b H_t,      H_t = sum_j w_j E[ h(d_{t-j}) | measured d_{t-j} ],    h(d) = 1 / (1 + (d / d_half)^4),
    w_j = exp(-j dt / tau) (normalised, j >= 0),

a = P(ON) when the pair has been far, a + b when it has been in full contact (h = 1) for longer than tau; (a + b) / a is the effective fold. The curve reported,
P(ON | distance held at d) = a + b h(d; d_half), has the same form as the instantaneous curve of `ep_curve` and is compared with the same reference fitted to the true
distance and true ON state. H is linear in h, so the errors-in-variables step is exact: E[h(d) | d_hat] = sum_g W[bin(d_hat), g] h(g) with the measurement kernel and the
population prior of `ep_curve`. d_half and tau are chosen from a grid by likelihood, (a, b) by Bernoulli maximum likelihood (concave: log of linear functions) on binned H.
Intervals: cells are resampled (the whole grid is refitted every time).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import lfilter

from simlive.stage7_validation import ep_curve as ep

DH = np.geomspace(100.0, 2000.0, 10)
TAU_S = np.array([10.0, 20.0, 45.0, 90.0, 180.0])
NHB = 64


def hill(d, d_half, n=4.0):
    return 1.0 / (1.0 + (np.asarray(d, float) / d_half) ** n)


@dataclass
class HFit:
    a: float
    b: float
    d_half: float
    tau_s: float
    nll: float
    n: float = 4.0

    @property
    def p_far(self) -> float:
        return self.a

    @property
    def fold(self) -> float:
        return (self.a + self.b) / self.a

    def g(self, d):
        """P(ON) for a distance held at d for longer than tau."""
        return self.a + self.b * hill(d, self.d_half, self.n)

    def contrast(self, near: float = 150.0, far: float = 1000.0) -> float:
        return float(self.g(near) / self.g(far))


def _fit_ab(n: np.ndarray, s: np.ndarray, h: np.ndarray, se: float = 1.0, fp: float = 0.0) -> tuple[float, float, float]:
    """Maximum likelihood (a, b) of P_obs = fp + (se - fp)(a + b h) for binned counts: concave in (a, b), Newton with step halving. Returns (a, b, nll)."""
    ok = n > 0
    n, s, h = n[ok], s[ok], h[ok]
    scale = se - fp
    p0 = max(s.sum() / n.sum(), 1e-4)
    th = np.array([(p0 - fp) / scale, 0.0])
    X = np.stack([np.ones_like(h), h], 1) * scale
    lo, hi = 1e-7, 1 - 1e-7

    def nll(t):
        p = np.clip(fp + X @ t, lo, hi)
        return -(s * np.log(p) + (n - s) * np.log1p(-p)).sum(), p

    f, p = nll(th)
    for _ in range(40):
        g = X.T @ (s / p - (n - s) / (1 - p))
        Hm = (X * (s / p ** 2 + (n - s) / (1 - p) ** 2)[:, None]).T @ X
        try:
            step = np.linalg.solve(Hm + 1e-12 * np.eye(2), g)
        except np.linalg.LinAlgError:
            break
        t = 1.0
        while t > 1e-4:
            cand = th + t * step
            pc = fp + X @ cand
            if (pc > lo).all() and (pc < hi).all():
                fc, _ = nll(cand)
                if fc <= f + 1e-10:
                    break
            t *= 0.5
        else:
            break
        th, f_new = cand, fc
        if f - f_new < 1e-9:
            f = f_new
            break
        f = f_new
        _, p = nll(th)
    return float(th[0]), float(th[1]), float(f)


def eiv_weights(K: np.ndarray, prior: np.ndarray) -> np.ndarray:
    num = K * prior[None, :]
    return num / np.maximum(num.sum(1, keepdims=True), 1e-300)


class HistoryData:
    """Per (d_half, tau): per-cell binned counts of (frames, ON, H) so that resampling cells is a sum of rows."""

    def __init__(self, bins: np.ndarray, y: np.ndarray, W: np.ndarray, frame_s: float, grid: np.ndarray = ep.D_GRID):
        self.n_cells, self.T = bins.shape
        self.frame_s = frame_s
        miss = bins < 0
        bsafe = np.where(miss, 0, bins)
        self.stats = {}
        col = np.arange(self.n_cells)[:, None] * NHB
        for dh in DH:
            hb = W @ hill(grid, dh)                                              # E[h | bin]
            fill = hb[bsafe[~miss]].mean() if (~miss).any() else 0.5
            x = np.where(miss, fill, hb[bsafe])
            for tau in TAU_S:
                al = np.exp(-frame_s / tau)
                H, _ = lfilter([1 - al], [1, -al], x, axis=1, zi=al * x[:, :1])
                hbin = np.minimum((H * NHB).astype(int), NHB - 1)
                idx = (col + hbin).ravel()
                n = np.bincount(idx, minlength=self.n_cells * NHB).reshape(self.n_cells, NHB).astype(float)
                s = np.bincount(idx, weights=y.ravel(), minlength=self.n_cells * NHB).reshape(self.n_cells, NHB)
                sh = np.bincount(idx, weights=H.ravel(), minlength=self.n_cells * NHB).reshape(self.n_cells, NHB)
                self.stats[(dh, tau)] = (n, s, sh)

    def fit(self, cells: np.ndarray | None = None, se: float = 1.0, fp: float = 0.0, keys=None) -> HFit:
        best = None
        for key in (keys or self.stats):
            n, s, sh = self.stats[key]
            if cells is not None:
                n, s, sh = n[cells], s[cells], sh[cells]
            n, s, sh = n.sum(0), s.sum(0), sh.sum(0)
            hc = np.where(n > 0, sh / np.maximum(n, 1e-12), 0.0)
            a, b, f = _fit_ab(n, s, hc, se, fp)
            if best is None or f < best.nll:
                best = HFit(a, b, key[0], key[1], f)
        return best


def analyse_history(pop: ep.Population, d_meas: np.ndarray, K: np.ndarray, y: np.ndarray | None = None, correct_calls: tuple | None = None,
                    n_boot: int = 0, seed: int = 0) -> dict:
    """Fit the history model to (measured distance of every frame, ON estimate y; the true ON state if None). K: measurement kernel of the analysis."""
    cells = np.unique(pop.cell_id)
    T = int((pop.cell_id == cells[0]).sum())
    ok = np.isfinite(d_meas)
    bins = np.where(ok, np.clip(np.digitize(np.where(ok, d_meas, 0.0), ep.EDGES) - 1, 0, len(ep.EDGES) - 2), -1).reshape(len(cells), T)
    yy = (pop.on_true.astype(float) if y is None else np.asarray(y, float)).reshape(len(cells), T)
    counts = np.bincount(bins[bins >= 0].ravel(), minlength=len(ep.EDGES) - 1).astype(float)
    prior = ep.population_prior(counts, K)
    W = eiv_weights(K, prior)
    data = HistoryData(bins, yy, W, pop.setup.frame_s)
    se, fp = correct_calls if correct_calls else (1.0, 0.0)
    fit = data.fit(se=se, fp=fp)
    out = {"fit": fit, "contrast": fit.contrast(), "n_frames": int(ok.sum()), "n_on": float(yy[ok.reshape(len(cells), T)].sum()), "frames_used": float(ok.mean()),
           "tau_s": fit.tau_s}
    if n_boot:
        rng = np.random.default_rng(seed + 1)
        bs = [data.fit(cells=rng.integers(0, len(cells), len(cells)), se=se, fp=fp) for _ in range(n_boot)]
        out["boot"] = bs
        for name in ("d_half", "p_far", "fold", "tau_s"):
            out[f"{name}_ci"] = ep.interval([getattr(b, name) for b in bs])
        out["contrast_ci"] = ep.interval([b.contrast() for b in bs])
    return out
