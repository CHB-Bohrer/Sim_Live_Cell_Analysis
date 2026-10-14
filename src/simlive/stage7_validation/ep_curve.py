"""Stage 7b: how well can the relation between promoter-enhancer distance and P(burst ON) be recovered, given measurement errors?

The question (user, 2026-10-14): readout = P(promoter in a burst) as a FUNCTION of the E-P distance (the whole curve, not one threshold), enhancer
with a long working distance, contact fraction 0.1. This module is the machinery; `scripts/sweep_ep_curve.py` runs it over error parameters;
docs/EP_CURVE.md has the method and the measured numbers.

Chain (each step is what a real analysis does, with the simulation's truth kept alongside):
  1. `simulate_population`: many cells; a stand-in distance process (hidden looped / open switch; the fraction of time with d < working distance
     is set EXACTLY to `contact_fraction`) drives the exact transcription model (stage1_chromatin/transcription.py, coupling "hill": f(d) = 1 / (1 +
     (d / working_distance)^n) multiplies k_on by 1 + (fold - 1) f(d)). Truth kept: d(t), promoter ON(t), Pol II starts, the MS2 signal.
     The distance process is INVENTED (like toy_loci); with the library, d(t) comes from polymer simulations once calibrated.
  2. `measure_distance`: what the analysis sees: localization precision per locus, colour shift, 2D projection (images are 2D) or 3D, missed
     frames, wrong picks.
  3. `estimate_on`: the observed MS2 trace (noise model of the real pipeline) -> Bayesian deconvolution (stage6_analysis/deconvolve.py) -> burst
     rule -> per-frame P(ON) estimate, with its sensitivity / false-positive rate against the truth.
  4. `fit_curve`: P(ON | d) = p_far (1 + (F - 1) / (1 + (d / d_half)^n)) fitted to (measured distance, ON estimate) pairs. Errors-in-variables: the
     latent true distance is integrated out with the known measurement kernel and a population prior of true distances (Richardson-Lucy from the
     measured distances), so the localization noise does not flatten the curve; `naive=True` ignores it for comparison. Cells are resampled
     for the confidence intervals (frames of one cell are not independent).
Reference = the same fit to the TRUE distance and TRUE ON state of a very large population (it contains the dynamics: the promoter responds to the
distance with a delay, so P(ON | d now) is not the stationary k_on / (k_on + k_off)).
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import erf, expit, ive

from simlive.stage1_chromatin import transcription as tx

D_GRID = np.arange(10.0, 3010.0, 10.0)             # latent true distances (nm)
EDGES = np.r_[np.arange(0.0, 3000.0, 25.0), 3000.0]  # measured-distance bins (nm); larger values go into the last bin


# ----------------------------------------------------------------------------- the experiment
@dataclass(frozen=True)
class EPSetup:
    working_distance_nm: float = 400.0      # f(d) = 0.5 here ("how far the enhancer acts"; the old sharp contact was 150 nm)
    hill_n: float = 4.0
    contact_fraction: float = 0.10          # fraction of time with d < working distance
    loop_time_s: float = 240.0              # mean time spent inside the working distance per visit
    relax_s: float = 30.0                   # correlation time of the fluctuations inside a state
    fold: float = 12.0
    n_cells: int = 600
    movie_s: float = 3600.0
    frame_s: float = 10.0
    exposure_s: float = 1.0
    block_s: float = 5.0                    # step of the distance process
    tx: dict = field(default_factory=dict)  # overrides of transcription.DEFAULTS (k_on, k_off, k_init, ...)

    def params(self, coupled: bool = True) -> tx.TranscriptionParams:
        base = dict(self.tx)
        base.update({"coupling": "hill" if coupled else "none", "d_half_nm": self.working_distance_nm, "hill_n": self.hill_n,
                     "fold": self.fold if coupled else 1.0})
        if not coupled:
            base["f_const"] = 0.0
        return tx.TranscriptionParams.from_dict(base)

    @property
    def n_frames(self) -> int:
        return int(round(self.movie_s / self.frame_s))


def ep_distance(n: int, dt_s: float, rng: np.random.Generator, working_nm: float, contact_fraction: float, loop_time_s: float,
                relax_s: float, close_factor: float = 0.4, open_factor: float = 2.0) -> np.ndarray:
    """INVENTED E-P distance (nm), n steps of dt_s: hidden looped / open switch, confined fluctuations around close_factor x working distance
    (looped) or open_factor x working distance (open, log-normal). The looped fraction is solved so that the fraction of time with d <
    working distance is `contact_fraction` exactly (stationary): P(d < w) = pl P_loop + (1 - pl) P_open."""
    nrm = lambda z: 0.5 * (1.0 + erf(z / np.sqrt(2.0)))                                  # noqa: E731
    p_open = float(nrm(np.log(1.0 / open_factor) / 0.35))                                # P(d < w | open): ln-normal, sigma 0.35
    p_loop = float(nrm((1.0 / close_factor - 1.0) / 0.3) - nrm((-1.0 / close_factor - 1.0) / 0.3))   # P(|1 + 0.3 x| < w / close)
    pl = (contact_fraction - p_open) / (p_loop - p_open)
    if not 0.0 < pl < 1.0:
        raise ValueError("contact_fraction not reachable with these shape factors")
    k_out = 1.0 / loop_time_s
    k_in = k_out * pl / (1.0 - pl)
    state, a = rng.random() < pl, np.exp(-dt_s / relax_s)
    x, out = rng.normal(), np.empty(n)
    u = rng.random((n, 2))
    z = rng.normal(size=n)
    pin, pout = 1 - np.exp(-k_in * dt_s), 1 - np.exp(-k_out * dt_s)
    for i in range(n):
        if u[i, 0] < (pout if state else pin):
            state = not state
        x = a * x + np.sqrt(1 - a * a) * z[i]
        out[i] = abs(close_factor * working_nm * (1 + 0.3 * x)) if state else open_factor * working_nm * np.exp(0.35 * x)
    return out


@dataclass
class Population:
    setup: EPSetup
    cell_id: np.ndarray                 # per frame
    t_s: np.ndarray
    d_true: np.ndarray
    on_true: np.ndarray                 # 0/1
    trace_cu: dict                      # {cell: noisy MS2 trace, one Pol II = 1}
    starts: dict                        # {cell: true initiation times (s)}
    coupled: bool = True
    vec_nm: np.ndarray | None = None    # true E-P vector per frame (nm, lab frame), polymer-driven populations only
    calib: object = None                # ep_polymer.Calibration of a polymer-driven population
    blocks: np.ndarray | None = None    # (cells, n_blocks, 3) float32: the true E-P vector (nm) at every saved block (polymer-driven populations; for acquisition timing)

    def cells(self) -> np.ndarray:
        return np.unique(self.cell_id)


def simulate_population(setup: EPSetup, seed: int, coupled: bool = True, noise_sb: float = 0.05, noise_beta: float = 0.03,
                        with_trace: bool = True) -> Population:
    """n_cells cells of `movie_s` seconds. noise_sb / noise_beta: the MS2 trace noise (C.U.: sigma^2 = sb^2 + beta G) as the real pipeline would estimate it."""
    p = setup.params(coupled)
    ss = np.random.SeedSequence(seed).spawn(setup.n_cells)
    T = setup.n_frames
    ft = np.arange(T) * setup.frame_s
    n_blocks = int(np.ceil((setup.movie_s + 1000.0) / setup.block_s)) + 2
    ids, ts, ds, ons, traces, starts = [], [], [], [], {}, {}
    for c in range(setup.n_cells):
        rng = np.random.default_rng(ss[c])
        d = ep_distance(n_blocks, setup.block_s, rng, setup.working_distance_nm, setup.contact_fraction, setup.loop_time_s, setup.relax_s)
        cell = tx.simulate_cell(p, d, setup.block_s, ft, setup.exposure_s, int(rng.integers(2**31)))
        f = cell.frames
        ids.append(np.full(T, c)); ts.append(ft); ds.append(f.d_nm.to_numpy()); ons.append(f.promoter_on.to_numpy())
        starts[c] = cell.events.init_s.to_numpy()
        if with_trace:
            g = f.ms2_loops.to_numpy() / p.n_loops
            traces[c] = g + rng.normal(size=T) * np.sqrt(noise_sb ** 2 + noise_beta * np.clip(g, 0, None))
    return Population(setup, np.concatenate(ids), np.concatenate(ts), np.concatenate(ds), np.concatenate(ons), traces, starts, coupled)


# ----------------------------------------------------------------------------- what the microscope reports
@dataclass(frozen=True)
class Measurement:
    sigma_nm: float = 0.0               # localization precision of the E-P DIFFERENCE vector, per axis (= sqrt(2) x per-locus precision)
    shift_nm: float = 0.0               # colour-channel offset left after correction (in plane, x direction)
    dims: int = 3                       # 2 = distance in the image plane (the 2D movies), 3 = from z-stacks
    miss_rate: float = 0.0              # fraction of frames without a distance
    outlier_rate: float = 0.0           # fraction of frames with a wrong pick: distance uniform in [0, 2000] nm
    assumed_shift_nm: float = 0.0       # what the analysis believes the shift to be (correction applied in the fit's kernel)


def measure_distance(d_true: np.ndarray, m: Measurement, rng: np.random.Generator) -> np.ndarray:
    """Measured distance (nm) per frame, NaN where missed. The pair's orientation is redrawn every frame (isotropic), which is exact for the
    magnitude of the 3D noise and a good approximation for the projection."""
    n = len(d_true)
    cos = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    sin = np.sqrt(1 - cos ** 2)
    v = np.stack([d_true * sin * np.cos(phi) + m.shift_nm, d_true * sin * np.sin(phi), d_true * cos], 1)
    v = v + rng.normal(size=(n, 3)) * m.sigma_nm
    out = np.sqrt((v[:, :m.dims] ** 2).sum(1))
    if m.outlier_rate > 0:
        bad = rng.random(n) < m.outlier_rate
        out[bad] = rng.uniform(0, 2000, bad.sum())
    if m.miss_rate > 0:
        out[rng.random(n) < m.miss_rate] = np.nan
    return out


def _kernel_matrix(sigma_nm: float, dims: int, shift_nm: float, edges: np.ndarray = EDGES, grid: np.ndarray = D_GRID) -> np.ndarray:
    """K[b, g] = P(measured distance in bin b | true distance grid[g]) for isotropic orientation, Gaussian noise sigma per axis, in-plane shift.
    3D: noncentral chi (3 dof) with the offset |d u + s| averaged over cos(theta); 2D: Rice distribution of the in-plane component, averaged over
    the orientation. Columns sum to 1 (mass beyond the last edge is put into the last bin)."""
    sigma = max(float(sigma_nm), 9.0)           # the kernel is evaluated at bin centres (25 nm bins): below ~9 nm it would miss mass
    r = 0.5 * (edges[:-1] + edges[1:])
    w = np.diff(edges)
    if dims == 3:
        x, wt = np.polynomial.legendre.leggauss(25)                       # cos theta
        m = np.sqrt(np.maximum(grid[None, :] ** 2 + shift_nm ** 2 + 2 * grid[None, :] * shift_nm * x[:, None], 0.0))   # (node, g)
        wt = wt / 2.0
        mm = np.maximum(m, 1e-6)
        pdf = (r[None, None, :] / (mm[:, :, None] * sigma * np.sqrt(2 * np.pi))) * (
            np.exp(-0.5 * ((r[None, None, :] - mm[:, :, None]) / sigma) ** 2) - np.exp(-0.5 * ((r[None, None, :] + mm[:, :, None]) / sigma) ** 2))
    else:
        x, wx = np.polynomial.legendre.leggauss(14)                       # cos theta -> in-plane size d sin theta
        ph = (np.arange(12) + 0.5) * 2 * np.pi / 12
        sin = np.sqrt(1 - x[:, None] ** 2)
        m = np.sqrt(np.maximum((grid[None, None, :] * sin[:, :, None]) ** 2 + shift_nm ** 2
                               + 2 * grid[None, None, :] * sin[:, :, None] * shift_nm * np.cos(ph)[None, :, None], 0.0)).reshape(-1, len(grid))
        wt = (wx[:, None] / 2.0 * np.ones((1, 12)) / 12.0).reshape(-1)
        z = r[None, None, :] * m[:, :, None] / sigma ** 2
        pdf = (r[None, None, :] / sigma ** 2) * ive(0, z) * np.exp(z - 0.5 * (r[None, None, :] ** 2 + m[:, :, None] ** 2) / sigma ** 2)
    k = (wt[:, None, None] * pdf).sum(0) * w[None, :]                     # (g, b)
    k = k.T
    k /= np.maximum(k.sum(0, keepdims=True), 1e-300)
    return k


_KCACHE: dict = {}


def kernel_matrix(sigma_nm: float, dims: int, shift_nm: float, edges: np.ndarray = EDGES, grid: np.ndarray = D_GRID) -> np.ndarray:
    """Cached `_kernel_matrix` (the 2D one takes ~1.5 s)."""
    key = (round(float(sigma_nm), 3), int(dims), round(float(shift_nm), 3), edges is EDGES and grid is D_GRID)
    if key[3] and key in _KCACHE:
        return _KCACHE[key]
    k = _kernel_matrix(sigma_nm, dims, shift_nm, edges, grid)
    if key[3]:
        _KCACHE[key] = k
    return k


def population_prior(counts: np.ndarray, K: np.ndarray, n_iter: int = 150) -> np.ndarray:
    """Richardson-Lucy / EM estimate of the distribution of TRUE distances on D_GRID from the histogram `counts` of measured ones."""
    pi = np.full(K.shape[1], 1.0 / K.shape[1])
    c = counts / max(counts.sum(), 1.0)
    for _ in range(n_iter):
        pred = K @ pi
        pi = pi * (K.T @ (c / np.maximum(pred, 1e-12)))
        pi /= pi.sum()
    return pi


# ----------------------------------------------------------------------------- P(ON) from the MS2 trace
@dataclass
class OnEstimate:
    """Posterior samples of the Pol II starts of every cell (counts per frame, uint8) and what the burst rule makes of them."""
    pop: Population
    counts: list           # per cell (n_keep, T), in the order of pop.cells()
    noise_sb: float
    noise_beta: float

    def p_on(self, w_frames: float = 5.0, g_b: float = 2.0) -> np.ndarray:
        """Per-frame P(burst ON) (fraction of samples that the rule calls ON), aligned with the population's frames."""
        from simlive.stage6_analysis.deconvolve import burst_calls

        return np.concatenate([burst_calls(c, self.pop.setup.frame_s, w_frames, g_b) for c in self.counts])

    def scores(self, p_on: np.ndarray) -> dict:
        on = self.pop.on_true.astype(bool)
        call = p_on >= 0.5
        return {"sensitivity": float(p_on[on].mean()), "false_positive": float(p_on[~on].mean()),
                "hard_sensitivity": float(call[on].mean()), "hard_false_positive": float(call[~on].mean())}

    def best_rule(self, ws=(1, 2, 3, 5, 8, 12), gbs=(1, 1.5, 2, 3, 4, 6)) -> tuple[float, float]:
        """The (w_frames, g_b) with the best balanced accuracy of the hard call AGAINST THE TRUTH (simulation only: for real data the rule has to come
        from the start rate, see the page's burst-rule calculator)."""
        best, arg = -1.0, (5.0, 2.0)
        for w in ws:
            for g in gbs:
                sc = self.scores(self.p_on(w, g))
                bal = 0.5 * (sc["hard_sensitivity"] + 1.0 - sc["hard_false_positive"])
                if bal > best:
                    best, arg = bal, (float(w), float(g))
        return arg


def estimate_on(pop: Population, n_keep: int = 200, n_burn: int = 150, workers: int = 1, seed: int = 0) -> OnEstimate:
    """Deconvolve every cell's observed MS2 trace: noise model and population prior estimated from the traces as the real pipeline does
    (`prepare_traces`), one gibbs_shift chain per cell (no convergence check: use the pipeline's `auto` for real analyses)."""
    from simlive.stage6_analysis.deconvolve import prepare_traces

    s = pop.setup
    T = s.n_frames
    cfg = {"acquisition": {"frame_interval_s": s.frame_s, "n_frames": T}, "loci": {"transcription": {"params": dict(s.params(pop.coupled).__dict__)}}}
    rows = [(c, t, v) for c, tr in pop.trace_cu.items() for t, v in enumerate(tr)]
    traces, kern, dt, prior, sb, beta, _, _ = prepare_traces(pd.DataFrame(rows, columns=["cell_id", "t", "intensity_corr"]), cfg, unit=1.0)
    cells = list(traces)
    jobs = [(traces[c], kern, dt, prior, sb, beta, n_keep, n_burn, seed + int(c)) for c in cells]
    if workers > 1:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(workers) as ex:
            counts = list(ex.map(_on_job, jobs, chunksize=max(1, len(jobs) // (4 * workers))))
    else:
        counts = [_on_job(j) for j in jobs]
    return OnEstimate(pop, counts, sb, beta)


def _on_job(a):
    from simlive.stage6_analysis.deconvolve import deconvolve_trace

    trace, kern, dt, prior, sb, beta, n_keep, n_burn, seed = a
    return deconvolve_trace(trace, kern, dt, prior, sb, beta, n_keep=n_keep, n_burn=n_burn, thin=2, seed=seed)["counts"]


# ----------------------------------------------------------------------------- the curve
def curve(d, p_far: float, fold: float, d_half: float, n: float) -> np.ndarray:
    return p_far * (1.0 + (fold - 1.0) / (1.0 + (np.asarray(d, float) / d_half) ** n))


@dataclass
class Aggregates:
    """Per cell and measured-distance bin: number of frames and sum of the ON estimate."""
    n: np.ndarray       # (cells, bins)
    s: np.ndarray

    @classmethod
    def build(cls, cell_id: np.ndarray, d_meas: np.ndarray, y: np.ndarray, edges: np.ndarray = EDGES) -> "Aggregates":
        ok = np.isfinite(d_meas) & np.isfinite(y)
        cells, inv = np.unique(cell_id[ok], return_inverse=True)
        b = np.clip(np.digitize(d_meas[ok], edges) - 1, 0, len(edges) - 2)
        nb = len(edges) - 1
        n = np.zeros((len(cells), nb)); s = np.zeros((len(cells), nb))
        np.add.at(n, (inv, b), 1.0)
        np.add.at(s, (inv, b), y[ok])
        return cls(n, s)


_DH = np.geomspace(60.0, 3000.0, 18)
_FG = np.geomspace(0.15, 80.0, 26)
_PG = np.geomspace(0.003, 0.6, 36)
_COARSE = (np.geomspace(60.0, 3000.0, 10), np.geomspace(0.15, 80.0, 12), np.geomspace(0.003, 0.6, 16))      # for the bootstrap refits


def _nll(theta, nb, sb, W, grid, n_shape, se, fp, free_shape):
    p_far, fold, dh = expit(theta[0]), np.exp(theta[1]), np.clip(np.exp(theta[2]), 40.0, 4000.0)
    nsh = np.exp(theta[3]) if free_shape else n_shape
    g = np.clip(curve(grid, p_far, fold, dh, nsh), 1e-6, 1 - 1e-6)
    m = np.clip(fp + (se - fp) * (W @ g), 1e-6, 1 - 1e-6)
    return -(sb * np.log(m) + (nb - sb) * np.log(1 - m)).sum()


@dataclass
class Fit:
    p_far: float
    fold: float
    d_half: float
    n: float
    nll: float
    ok: bool = True

    def g(self, d):
        return curve(d, self.p_far, self.fold, self.d_half, self.n)

    def contrast(self, near: float = 150.0, far: float = 1000.0) -> float:
        """P(ON) at `near` over P(ON) at `far`: the identifiable size of the effect (the fold F alone is not identified when there is no effect:
        any F goes with a tiny d_half)."""
        return float(self.g(near) / self.g(far))


def fit_curve(agg: Aggregates, K: np.ndarray, prior: np.ndarray, n_shape: float = 4.0, naive: bool = False, se: float = 1.0, fp: float = 0.0,
              free_shape: bool = False, start=None, coarse: bool = False, grid: np.ndarray = D_GRID, edges: np.ndarray = EDGES) -> Fit:
    """Maximum-likelihood fit of the P(ON | d) curve to the binned data. EIV: W[b, g] = P(true d = g | measured bin b) from K and the population
    prior; naive: W = identity on the bin centres (the measured distance is taken as the true one). se / fp: sensitivity and false-positive rate
    of the ON estimate (observed P = fp + (se - fp) P(ON)); defaults = the estimate is the truth."""
    nb, sb = agg.n.sum(0), agg.s.sum(0)
    if naive:
        centres = 0.5 * (edges[:-1] + edges[1:])
        W, g = np.eye(len(centres)), centres
    else:
        num = K * prior[None, :] if K.shape[0] == len(nb) and K.shape[1] == len(grid) else None
        W = num / np.maximum(num.sum(1, keepdims=True), 1e-300)
        g = grid
    # 1. coarse grid over (d_half, F, p_far): P(ON | measured bin) = fp + (se - fp) p_far (1 + (F - 1) a_b(d_half)), a_b = E[h(d) | bin b]
    DH, FG, PG = _COARSE if coarse else (_DH, _FG, _PG)
    a_dh = np.stack([W @ (1.0 / (1.0 + (g / dh) ** n_shape)) for dh in DH])                         # (dh, b)
    u = 1.0 + (FG[None, :, None] - 1.0) * a_dh[:, None, :]                                           # (dh, F, b)
    m = fp + (se - fp) * PG[None, None, :, None] * u[:, :, None, :]                                   # (dh, F, p, b)
    m = np.clip(m, 1e-6, 1 - 1e-6)
    nll = -(sb * np.log(m) + (nb - sb) * np.log1p(-m)).sum(-1)
    i, j, k = np.unravel_index(np.argmin(nll), nll.shape)
    x0 = [np.log(PG[k] / (1 - PG[k])), np.log(FG[j]), np.log(DH[i])] + ([np.log(n_shape)] if free_shape else [])
    # 2. local polish (Nelder-Mead from the grid optimum); keep the better of the two
    best = minimize(_nll, x0, args=(nb, sb, W, g, n_shape, se, fp, free_shape), method="Nelder-Mead",
                    options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 400})
    if not best.fun <= nll[i, j, k]:
        best.x, best.fun = np.asarray(x0), nll[i, j, k]
    th = best.x
    th = np.r_[th[:2], np.clip(th[2], np.log(60.0), np.log(3000.0)), th[3:]]
    return Fit(float(expit(th[0])), float(np.exp(th[1])), float(np.exp(th[2])), float(np.exp(th[3])) if free_shape else n_shape, float(best.fun),
               bool(best.success))


def bootstrap_fits(agg: Aggregates, K, prior, n_boot: int = 200, seed: int = 0, **kw) -> list[Fit]:
    """Resample CELLS with replacement (the frames of a cell are correlated) and refit; the prior of true distances is kept fixed."""
    rng = np.random.default_rng(seed)
    n_c = agg.n.shape[0]
    base = fit_curve(agg, K, prior, **kw)
    out = []
    for _ in range(n_boot):
        pick = rng.integers(0, n_c, n_c)
        out.append(fit_curve(Aggregates(agg.n[pick], agg.s[pick]), K, prior, coarse=True, **kw))
    return out


def interval(vals, level: float = 0.95) -> tuple[float, float]:
    a = (1 - level) / 2
    return float(np.quantile(vals, a)), float(np.quantile(vals, 1 - a))


# ----------------------------------------------------------------------------- one experiment
def analyse(pop: Population, meas: Measurement, y: np.ndarray | None = None, naive: bool = False, correct_calls: tuple | None = None,
            n_boot: int = 0, seed: int = 0, n_shape: float = 4.0, free_shape: bool = False) -> dict:
    """Curve fit to (measured distance of `meas`, ON estimate y; the true ON state if None). correct_calls = (sensitivity, false_positive) of y
    for the attenuation correction. Returns the fit, bootstrap intervals (n_boot > 0) and the curve on D_GRID."""
    rng = np.random.default_rng(seed)
    d_meas = measure_distance(pop.d_true, meas, rng)
    K = kernel_matrix(meas.sigma_nm, meas.dims, meas.assumed_shift_nm)
    return analyse_measured(pop, d_meas, K, y=y, naive=naive, correct_calls=correct_calls, n_boot=n_boot, seed=seed, n_shape=n_shape, free_shape=free_shape)


def analyse_measured(pop: Population, d_meas: np.ndarray, K: np.ndarray, y: np.ndarray | None = None, naive: bool = False, correct_calls: tuple | None = None,
                     n_boot: int = 0, seed: int = 0, n_shape: float = 4.0, free_shape: bool = False) -> dict:
    """The instantaneous curve P(ON | distance now) fitted to already measured distances (NaN = no measurement) with the analysis' measurement kernel K."""
    yy = pop.on_true.astype(float) if y is None else np.asarray(y, float)
    agg = Aggregates.build(pop.cell_id, d_meas, yy)
    prior = population_prior(agg.n.sum(0), K)
    se, fp = correct_calls if correct_calls else (1.0, 0.0)
    kw = dict(n_shape=n_shape, naive=naive, se=se, fp=fp, free_shape=free_shape)
    fit = fit_curve(agg, K, prior, **kw)
    out = {"fit": fit, "contrast": fit.contrast(), "n_frames": int(agg.n.sum()), "n_on": float(agg.s.sum()), "frames_used": float(agg.n.sum() / max(len(pop.d_true), 1))}
    if n_boot:
        bs = bootstrap_fits(agg, K, prior, n_boot=n_boot, seed=seed + 1, **kw)
        out["boot"] = bs
        for name in ("fold", "d_half", "p_far"):
            out[f"{name}_ci"] = interval([getattr(b, name) for b in bs])
        out["contrast_ci"] = interval([b.contrast() for b in bs])
    return out


def reference_curve(setup: EPSetup, n_cells: int = 4000, seed: int = 12345) -> dict:
    """The truth the estimators are judged against: the same curve fitted to the TRUE distance and TRUE ON state of a very large population,
    plus the distribution of true distances (for weighting) and the observed contact fraction."""
    pop = simulate_population(replace(setup, n_cells=n_cells), seed, with_trace=False)
    out = analyse(pop, Measurement(sigma_nm=0.5, dims=3))
    hist, _ = np.histogram(pop.d_true, bins=np.r_[D_GRID - 5.0, D_GRID[-1] + 5.0], density=False)
    out["d_pdf"] = hist / hist.sum()
    out["contact_fraction_observed"] = float((pop.d_true < setup.working_distance_nm).mean())
    out["p_on_overall"] = float(pop.on_true.mean())
    return out


def curve_error(fit: Fit, ref: dict) -> float:
    """RMS difference between the fitted and the reference P(ON | d) curve, weighted by how often each true distance occurs."""
    r = ref["fit"]
    return float(np.sqrt((ref["d_pdf"] * (fit.g(D_GRID) - r.g(D_GRID)) ** 2).sum()))
