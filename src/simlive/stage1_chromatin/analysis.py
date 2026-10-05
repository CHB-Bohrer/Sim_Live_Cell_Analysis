"""Analysis of chromatin simulations: contact maps, contact probability P(s), boundary insulation, loop 'dots', locus MSD.

Genomic coordinates are in monomers (1 monomer = 1 kb in this project). Contacts: two monomers closer than `cutoff`
monomer diameters (default 3.0). All maps are averaged over the snapshots passed in.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def load_snapshots(sim_dir: str | Path, skip: int = 0, max_n: int | None = None) -> list[np.ndarray]:
    """Full-chain conformations saved by the simulation (polychrom HDF5 blocks)."""
    from polychrom.hdf5_format import list_URIs, load_URI

    uris = list_URIs(str(Path(sim_dir) / "blocks"))[skip:]
    if max_n:
        uris = uris[:max_n]
    return [np.asarray(load_URI(u)["pos"], np.float32) for u in uris]


def _pairs(conf: np.ndarray, cutoff: float) -> np.ndarray:
    return cKDTree(conf).query_pairs(cutoff, output_type="ndarray")


def contact_map(confs: list[np.ndarray], bin_size: int = 10, cutoff: float = 3.0) -> np.ndarray:
    """Coarse (n_bins x n_bins) contact frequency map, symmetric, averaged over conformations. Bin = bin_size monomers."""
    n = len(confs[0])
    nb = int(np.ceil(n / bin_size))
    acc = np.zeros(nb * nb)
    for c in confs:
        p = _pairs(c, cutoff)
        if len(p):
            acc += np.bincount((p[:, 0] // bin_size) * nb + (p[:, 1] // bin_size), minlength=nb * nb)
    m = acc.reshape(nb, nb) / len(confs)
    return m + m.T - np.diag(np.diag(m))


def separation_histogram(confs: list[np.ndarray], cutoff: float = 3.0) -> np.ndarray:
    """Mean number of contacts at each genomic separation s (index = s, in monomers)."""
    n = len(confs[0])
    h = np.zeros(n)
    for c in confs:
        p = _pairs(c, cutoff)
        if len(p):
            h += np.bincount(np.abs(p[:, 1] - p[:, 0]), minlength=n)[:n]
    return h / len(confs)


def contact_probability(hist: np.ndarray, n_bins: int = 40, s_min: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """P(s): contacts at separation s divided by the number of monomer pairs at that separation (n - s), log-binned."""
    n = len(hist)
    s = np.arange(n)
    edges = np.unique(np.round(np.logspace(np.log10(s_min), np.log10(n / 2), n_bins)).astype(int))
    centers, vals = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (s >= a) & (s < b)
        pairs_available = (n - s[sel]).sum()
        if pairs_available > 0:
            centers.append(np.exp(np.log(s[sel]).mean()))
            vals.append(hist[sel].sum() / pairs_available)
    return np.array(centers), np.array(vals)


def loglog_slope(s: np.ndarray, p: np.ndarray, lo: float, hi: float) -> float:
    sel = (s >= lo) & (s <= hi) & (p > 0)
    return float(np.polyfit(np.log(s[sel]), np.log(p[sel]), 1)[0]) if sel.sum() >= 3 else float("nan")


def observed_over_expected(m: np.ndarray) -> np.ndarray:
    """Divide each diagonal by its mean (distance decay removed) so structure stands out."""
    nb = m.shape[0]
    out = np.zeros_like(m, dtype=float)
    for d in range(nb):
        diag = np.diagonal(m, d)
        mean = diag.mean()
        if mean > 0:
            idx = np.arange(nb - d)
            out[idx, idx + d] = diag / mean
            out[idx + d, idx] = diag / mean
    return out


def boundary_insulation(m: np.ndarray, boundaries: list[int], bin_size: int, window_kb: int = 150,
                        use_oe: bool = False) -> list[dict]:
    """For each boundary: mean contact frequency WITHIN the flanking windows vs ACROSS the boundary.
    ratio = within / across; 1 = no insulation (Fudenberg et al. report ~2-fold lower contact across TAD borders).

    Contacts fall steeply with genomic distance and 'within' pixels are closer to the diagonal than 'across' pixels,
    so on a raw map the ratio is >> 1 even with NO boundary. Use use_oe=True (observed/expected per diagonal), which
    compares pixels at equal genomic separation; then 1 really means no insulation."""
    if use_oe:
        m = observed_over_expected(m)
    out, w = [], max(2, window_kb // bin_size)
    nb = m.shape[0]
    for x in boundaries:
        b = x // bin_size
        if b - w < 0 or b + w > nb:
            continue
        left, right = m[b - w:b, b - w:b], m[b:b + w, b:b + w]
        across = m[b - w:b, b:b + w]
        iu = np.triu_indices(w, k=1)
        within = 0.5 * (left[iu].mean() + right[iu].mean())
        out.append({"pos_kb": int(x), "within": float(within), "across": float(across.mean()),
                    "ratio": float(within / max(across.mean(), 1e-12))})
    return out


def corner_peak_enrichment(oe: np.ndarray, boundaries: list[int], bin_size: int, min_tad_kb: int = 300) -> list[float]:
    """Enrichment of observed/expected contacts at the corner (b_i, b_{i+1}) of each domain relative to the ring of
    bins around it (3x3 centre vs the surrounding 9x9 minus 5x5). >1 means a 'dot' at the loop anchors."""
    out = []
    for x0, x1 in zip(boundaries[:-1], boundaries[1:]):
        if x1 - x0 < min_tad_kb:
            continue
        i, j = x0 // bin_size, x1 // bin_size
        if i < 5 or j + 5 >= oe.shape[0]:
            continue
        centre = oe[i - 1:i + 2, j - 1:j + 2].mean()
        ring = oe[i - 4:i + 5, j - 4:j + 5][_RING]
        out.append(float(centre / max(ring.mean(), 1e-12)))
    return out


def loop_dot_enrichment(oe: np.ndarray, pairs: list[tuple[int, int]], bin_size: int) -> list[float]:
    """Enrichment of observed/expected contacts at the anchors (a, b) of designed loops: the 3x3 bins at (a, b) divided
    by the mean of the surrounding ring (outside the central 5x5 of a 9x9 block). >1 = a 'dot' (corner peak)."""
    out = []
    for a, b in pairs:
        i, j = a // bin_size, b // bin_size
        if i < 5 or j + 5 >= oe.shape[0]:
            continue
        centre = oe[i - 1:i + 2, j - 1:j + 2].mean()
        ring = oe[i - 4:i + 5, j - 4:j + 5][_RING]
        out.append(float(centre / max(ring.mean(), 1e-12)))
    return out


_RING = np.ones((9, 9), bool)
_RING[2:7, 2:7] = False     # True for the cells of the 9x9 block that lie outside its central 5x5


def msd(traj: np.ndarray, lags: np.ndarray, subtract_com: bool = False, trim_frac: float = 0.1) -> np.ndarray:
    """Time- and bead-averaged mean squared displacement of tracked beads. traj: (T, n_beads, 3).
    lags in blocks. Beads within trim_frac of either chain end are ignored. subtract_com removes whole-chain motion."""
    t, nb, _ = traj.shape
    if subtract_com:
        traj = traj - traj.mean(1, keepdims=True)
    k = int(nb * trim_frac)
    x = traj[:, k:nb - k]
    return np.array([((x[l:] - x[:-l]) ** 2).sum(-1).mean() for l in lags])


def msd_exponent(lags: np.ndarray, m: np.ndarray, lo: float, hi: float) -> float:
    sel = (lags >= lo) & (lags <= hi)
    return float(np.polyfit(np.log(lags[sel]), np.log(m[sel]), 1)[0])


def radius_of_gyration(conf: np.ndarray) -> float:
    return float(np.sqrt(((conf - conf.mean(0)) ** 2).sum(1).mean()))


def across_insulation(m: np.ndarray, positions: list[int], bin_size: int, window_kb: int = 150) -> list[float]:
    """Insulation fold at each position = 1 / (mean observed/expected contact frequency ACROSS the position).

    The expectation is the map-wide average at the same genomic separation, so 1.0 means 'no insulation' and 2.0 means
    contacts across the position are half of what is typical at those distances. Unlike a within/across ratio this has
    no distance-decay bias: at random positions it is ~1 (checked on plain-polymer and loop-extrusion data)."""
    oe = observed_over_expected(m)
    w = max(2, window_kb // bin_size)
    out = []
    for x in positions:
        b = x // bin_size
        if b - w >= 0 and b + w <= oe.shape[0]:
            out.append(1.0 / max(float(oe[b - w:b, b:b + w].mean()), 1e-9))
    return out


def anchor_contact_enrichment(traj: np.ndarray, beads: np.ndarray, pairs: list[tuple[int, int]], s_curve: tuple,
                              cutoff: float = 3.0) -> list[float]:
    """Contact probability of designed loop anchors relative to the average monomer pair at the same genomic separation.

    traj: (T, n_tracked, 3) positions of tracked monomers at every block; beads: their monomer indices; pairs: (a, b)
    monomer indices (both must be tracked); s_curve: (s, P(s)) of the same simulation (same cutoff). Returns one
    enrichment per pair: P_contact(a, b) / P(b - a). A held loop makes it >> 1; with no loop it is ~1."""
    s, p = s_curve
    out = []
    for a, b in pairs:
        ia, ib = np.nonzero(beads == a)[0], np.nonzero(beads == b)[0]
        if len(ia) == 0 or len(ib) == 0:
            continue
        d = np.linalg.norm(traj[:, ia[0]] - traj[:, ib[0]], axis=1)
        expected = float(np.interp(np.log(b - a), np.log(s), p))
        out.append(float((d < cutoff).mean() / max(expected, 1e-9)))
    return out


def insulation_fold_pooled(m: np.ndarray, positions: list[int], bin_size: int, window_kb: int = 150) -> float:
    """1 / (mean over positions of the observed/expected contact frequency across the position). Averaging the contact
    frequency first and inverting after is robust; averaging per-position folds is inflated by a few near-zero values."""
    oe = observed_over_expected(m)
    w = max(2, window_kb // bin_size)
    vals = [oe[x // bin_size - w:x // bin_size, x // bin_size:x // bin_size + w].mean() for x in positions
            if x // bin_size - w >= 0 and x // bin_size + w <= oe.shape[0]]
    return 1.0 / max(float(np.mean(vals)), 1e-9)


def anchor_contact_pooled(traj: np.ndarray, beads: np.ndarray, pairs: list[tuple[int, int]], s_curve: tuple,
                          cutoff: float = 3.0) -> tuple[float, float]:
    """(observed, expected) contact probability pooled over all designed pairs: observed = fraction of (block, pair)
    samples in which the two anchor monomers are within `cutoff`; expected = the average P(s) at those separations.
    Pooling avoids a few unformed loops (zero contacts in the run) dominating a median."""
    s, p = s_curve
    obs, exp = [], []
    for a, b in pairs:
        ia, ib = np.nonzero(beads == a)[0], np.nonzero(beads == b)[0]
        if len(ia) == 0 or len(ib) == 0:
            continue
        d = np.linalg.norm(traj[:, ia[0]] - traj[:, ib[0]], axis=1)
        obs.append(float((d < cutoff).mean()))
        exp.append(float(np.interp(np.log(b - a), np.log(s), p)))
    return (float(np.mean(obs)), float(np.mean(exp))) if obs else (float("nan"), float("nan"))
