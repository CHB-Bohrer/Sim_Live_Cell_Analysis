"""E-P distances from the polymer simulations, for the E-P curve analysis (replaces the invented distance process of ep_curve.simulate_population).

Source: the saved validation simulations ($SIMLIVE_DATA/chromatin/validation/val_*, 10 simulations of a 10 Mb region, 1 monomer = 1 kb, the positions of every
10th monomer stored at every one of 1000 production blocks). The E-P pair is two stored monomers `sep_kb` apart (default 100 kb = 10 stored monomers); every such
pair along every chain is used (990 pairs x 10 simulations). A cell gets its distance trace by concatenating 1000-block segments of randomly chosen (simulation,
pair) combinations, each with its own random rotation (the pair's orientation relative to the optical axis is random); the seams are harmless because the
distance decorrelates within ~20 blocks. The conversion to nm and seconds has two numbers (`calibrate`):
  - nm per polymer unit: if the literature gives no usable scale, chosen so that the pair is within the working distance `contact_fraction` of the time
    (user, 2026-10-14: contact fraction 0.1, working distance 400 nm). The paper checked (Jusuf et al. 2025) gives no nm distances.
  - seconds per saved block: the single-locus MSD (centre-of-mass subtracted, 3D) at lag `lag_s` is matched to 6 D lag_s with D = 0.002 um^2/s (user: single locus,
    lag 1 min; the simulated MSD exponent is ~0.5, so D at another lag would give a very different block time).
Both are approximations; see docs/EP_CURVE.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from simlive.io.runs import DATA_ROOT
from simlive.stage1_chromatin import transcription as tx
from simlive.stage7_validation import ep_curve as ep

VAL_DIR = DATA_ROOT / "chromatin" / "validation"
SIMS = ("val_polymer_seed1", "val_polymer_seed2", "val_lef_free_seed1", "val_lef_free_seed2", "val_lef_leaky_seed1", "val_lef_leaky_seed2",
        "val_lef_ctcf_seed1", "val_lef_ctcf_seed2", "val_lef_loops_seed1", "val_lef_loops_seed2")
LAGS = np.array([1, 2, 5, 10, 20, 50, 100, 200, 400], float)        # blocks


@dataclass
class PolymerPool:
    rel: np.ndarray               # (sims, blocks, pairs, 3) float32: vector between the two monomers, polymer units
    msd_single: np.ndarray        # pooled single-locus MSD at LAGS (units^2, 3D, centre of mass of the tracked monomers removed)
    msd_rel: np.ndarray           # pooled MSD of the relative vector (the E-P separation vector)
    names: tuple
    sep_kb: float

    @property
    def n_blocks(self) -> int:
        return self.rel.shape[1]

    def distances(self) -> np.ndarray:
        return np.linalg.norm(self.rel, axis=-1)

    def autocorr_time_blocks(self) -> float:
        """Blocks until the autocorrelation of the E-P distance falls to 1/e (averaged over pairs and simulations)."""
        d = self.distances()
        z = d - d.mean(1, keepdims=True)
        ac = [float((z[:, :-l] * z[:, l:]).mean() / (z * z).mean()) for l in range(1, 400)]
        ac = np.array(ac)
        i = int(np.argmax(ac < 1 / np.e))
        return float(i + 1) if ac.min() < 1 / np.e else float("nan")


def load_pool(sep_kb: float = 100.0, sims: tuple = SIMS, root: Path = VAL_DIR) -> PolymerPool:
    rels, ms, mr = [], [], []
    for s in sims:
        d = root / s
        tp = np.load(d / "tracked_positions.npy").astype(np.float32)                 # (T, beads, 3)
        beads = np.load(d / "tracked_beads.npy")
        step = int(round(sep_kb / float(beads[1] - beads[0])))
        x = tp - tp.mean(1, keepdims=True)
        rel = (tp[:, step:] - tp[:, :-step]).astype(np.float32)
        rels.append(rel)
        dv = x[:, step:] - x[:, :-step]
        ms.append([float(((x[l:] - x[:-l]) ** 2).sum(2).mean()) for l in LAGS.astype(int)])
        mr.append([float(((dv[l:] - dv[:-l]) ** 2).sum(2).mean()) for l in LAGS.astype(int)])
    return PolymerPool(np.stack(rels), np.mean(ms, 0), np.mean(mr, 0), tuple(sims), sep_kb)


@dataclass(frozen=True)
class Calibration:
    nm_per_unit: float
    block_s: float
    alpha: float                  # exponent of the single-locus MSD below one block (extrapolation)
    msd1_nm2: float               # single-locus MSD at one block (nm^2, 3D)
    msd_nm2: tuple                # single-locus MSD at LAGS (nm^2)
    msd_rel_nm2: tuple
    lag_s: float
    D_um2_s: float
    applies_to: str

    def msd_single_s(self, tau_s):
        """Single-locus 3D MSD (nm^2) after tau_s seconds: the simulated table above one block, the power law with the small-lag exponent below it (the simulation
        stores nothing between blocks)."""
        tau = np.asarray(tau_s, float)
        L = tau / self.block_s
        tab = np.interp(np.log(np.maximum(L, 1e-9)), np.log(LAGS), np.log(self.msd_nm2))
        small = self.msd1_nm2 * np.maximum(L, 0.0) ** self.alpha
        return np.where(L < 1.0, small, np.exp(np.where(L > LAGS[-1], np.log(self.msd_nm2[-1]), tab)))


def calibrate(pool: PolymerPool, working_nm: float = 400.0, contact_fraction: float = 0.10, D_um2_s: float = 0.002, lag_s: float = 60.0,
              applies_to: str = "locus", nm_per_unit: float | None = None, block_s: float | None = None) -> Calibration:
    d = pool.distances()
    nm = float(nm_per_unit) if nm_per_unit else float(working_nm / np.quantile(d, contact_fraction))
    msd = np.array(pool.msd_single if applies_to == "locus" else pool.msd_rel, float)
    alpha = float(np.polyfit(np.log(LAGS[:2]), np.log(msd[:2]), 1)[0])
    target = 6.0 * D_um2_s * lag_s * 1e6 / nm ** 2                                  # units^2
    if block_s is None:
        if target < msd[0]:
            L = (target / msd[0]) ** (1.0 / alpha)
        else:
            L = float(np.exp(np.interp(np.log(target), np.log(msd), np.log(LAGS))))
        block_s = lag_s / L
    return Calibration(nm, float(block_s), alpha, float(np.array(pool.msd_single)[0]) * nm ** 2, tuple(np.array(pool.msd_single) * nm ** 2),
                       tuple(np.array(pool.msd_rel) * nm ** 2), lag_s, D_um2_s, applies_to)


def _rotation(rng: np.random.Generator) -> np.ndarray:
    q = rng.normal(size=4)
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def cell_vectors(pool: PolymerPool, calib: Calibration, n_blocks: int, rng: np.random.Generator) -> np.ndarray:
    """E-P vector (nm, lab frame) at every block: independent segments of the saved trajectories, each with a random orientation."""
    out, n = [], 0
    while n < n_blocks:
        s, p = int(rng.integers(pool.rel.shape[0])), int(rng.integers(pool.rel.shape[2]))
        seg = pool.rel[s, :, p, :].astype(np.float64) @ _rotation(rng).T
        out.append(seg * calib.nm_per_unit)
        n += len(seg)
    return np.concatenate(out)[:n_blocks]


def simulate_population(setup: ep.EPSetup, pool: PolymerPool, calib: Calibration, seed: int, coupled: bool = True, noise_sb: float = 0.05,
                        noise_beta: float = 0.03, with_trace: bool = True, keep_blocks: bool = True) -> ep.Population:
    """Like ep_curve.simulate_population with the polymer distances; the Population also carries the true E-P vector at every frame (`vec_nm`)."""
    p = setup.params(coupled)
    ss = np.random.SeedSequence(seed).spawn(setup.n_cells)
    T = setup.n_frames
    ft = np.arange(T) * setup.frame_s
    block = calib.block_s
    n_blocks = int(np.ceil((setup.movie_s + 1000.0) / block)) + 4
    ids, ts, ds, ons, vecs, blocks = [], [], [], [], [], []
    traces, starts = {}, {}
    for c in range(setup.n_cells):
        rng = np.random.default_rng(ss[c])
        v = cell_vectors(pool, calib, n_blocks, rng)
        d = np.linalg.norm(v, axis=1)
        cell = tx.simulate_cell(p, d, block, ft, setup.exposure_s, int(rng.integers(2**31)))
        f = cell.frames
        ids.append(np.full(T, c)); ts.append(ft); ds.append(f.d_nm.to_numpy()); ons.append(f.promoter_on.to_numpy())
        vecs.append(v[np.minimum(np.round(ft / block).astype(int), len(v) - 1)])
        if keep_blocks:
            blocks.append(v.astype(np.float32))
        starts[c] = cell.events.init_s.to_numpy()
        if with_trace:
            g = f.ms2_loops.to_numpy() / p.n_loops
            traces[c] = g + rng.normal(size=T) * np.sqrt(noise_sb ** 2 + noise_beta * np.clip(g, 0, None))
    pop = ep.Population(setup, np.concatenate(ids), np.concatenate(ts), np.concatenate(ds), np.concatenate(ons), traces, starts, coupled)
    pop.vec_nm = np.concatenate(vecs)
    pop.calib = calib
    pop.blocks = np.stack(blocks) if keep_blocks else None
    return pop


def reference_curve(setup: ep.EPSetup, pool: PolymerPool, calib: Calibration, n_cells: int = 4000, seed: int = 12345) -> dict:
    """The truth the estimators are judged against: the curve fitted to the TRUE distance and TRUE ON state of a large polymer-driven population."""
    pop = simulate_population(replace(setup, n_cells=n_cells), pool, calib, seed, with_trace=False)
    out = ep.analyse(pop, ep.Measurement(sigma_nm=0.5, dims=3))
    hist, _ = np.histogram(pop.d_true, bins=np.r_[ep.D_GRID - 5.0, ep.D_GRID[-1] + 5.0])
    out["d_pdf"] = hist / hist.sum()
    out["contact_fraction_observed"] = float((pop.d_true < setup.working_distance_nm).mean())
    out["p_on_overall"] = float(pop.on_true.mean())
    out["distance_autocorr_s"] = distance_autocorr_s(pop, setup.frame_s)
    return out


def distance_autocorr_s(pop: ep.Population, frame_s: float) -> float:
    """Time (s) for the autocorrelation of the E-P distance to fall to 1/e, from the frames (resolution = frame interval; NaN if it is faster than that)."""
    d = pop.d_true
    cells = np.unique(pop.cell_id)
    T = int((pop.cell_id == cells[0]).sum())
    x = d.reshape(len(cells), T)
    x = x - x.mean(1, keepdims=True)
    ac1 = float((x[:, :-1] * x[:, 1:]).mean() / (x * x).mean())
    return float(-frame_s / np.log(ac1)) if ac1 > 0 else float("nan")
