"""1D loop-extruding-factor (LEF) dynamics with CTCF-like boundary elements, in the style of the Mirny lab models
(Fudenberg et al. 2016 Cell Reports; Goloborodko et al. 2016 eLife; the polychrom loopExtrusion example).

Each LEF is two legs holding two monomers together (left leg at l, right leg at r, l < r). In every 1D step each leg
tries to step one monomer outward (left leg to l-1, right leg to r+1). A step is refused when
  * the target monomer is occupied by another LEF leg (LEFs cannot pass each other), or
  * the leg is standing on a boundary element that blocks that direction: a site at monomer x with p_left[x] stops a
    LEFT-moving leg that is at x with probability p_left[x] each attempt (p_right[x] likewise for right-moving legs).
A leg can therefore walk INTO a site and be held there, so a left-blocking site at x1 and a right-blocking site at
x2 > x1 (a *convergent* pair) hold a persistent loop x1..x2.
Each LEF unbinds with probability 1/lifetime per step and immediately rebinds at a random free adjacent pair of
monomers, so the number of bound LEFs is constant (the usual simplification).

Processivity (mean loop size a free LEF extrudes before unbinding) = 2 * speed * lifetime, with speed = 1 monomer
per step per leg.
"""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _seed(s):
    np.random.seed(s)


@njit(cache=True)
def _place(left, right, occ, m, n):
    """Put LEF m on a random free adjacent pair of monomers."""
    while True:
        i = np.random.randint(0, n - 1)
        if occ[i] == 0 and occ[i + 1] == 0:
            left[m] = i
            right[m] = i + 1
            occ[i] = m + 1
            occ[i + 1] = m + 1
            return


@njit(cache=True)
def _steps(left, right, occ, p_left, p_right, unbind_prob, n_steps):
    n = occ.shape[0]
    m_lefs = left.shape[0]
    for _ in range(n_steps):
        for k in range(m_lefs):
            m = np.random.randint(0, m_lefs)  # random order avoids a systematic bias between LEFs
            if np.random.random() < unbind_prob:
                occ[left[m]] = 0
                occ[right[m]] = 0
                _place(left, right, occ, m, n)
                continue
            lp = left[m]
            if lp > 0 and occ[lp - 1] == 0 and np.random.random() >= p_left[lp]:
                occ[lp] = 0
                left[m] = lp - 1
                occ[lp - 1] = m + 1
            rp = right[m]
            if rp < n - 1 and occ[rp + 1] == 0 and np.random.random() >= p_right[rp]:
                occ[rp] = 0
                right[m] = rp + 1
                occ[rp + 1] = m + 1


class LEF1D:
    """A population of LEFs on a chain of n monomers.

    n_lefs: number of bound LEFs (= n / separation).  lifetime_steps: mean bound time in 1D steps.
    ctcf: list of dicts {"pos": int, "blocks_left": float, "blocks_right": float} (stall probabilities 0..1).
    """

    def __init__(self, n: int, n_lefs: int, lifetime_steps: float, ctcf=(), seed: int = 0):
        self.n, self.n_lefs, self.lifetime = int(n), int(n_lefs), float(lifetime_steps)
        self.p_left = np.zeros(n, np.float64)
        self.p_right = np.zeros(n, np.float64)
        for c in ctcf:
            self.p_left[c["pos"]] = max(self.p_left[c["pos"]], c.get("blocks_left", 0.0))
            self.p_right[c["pos"]] = max(self.p_right[c["pos"]], c.get("blocks_right", 0.0))
        self.left = np.zeros(self.n_lefs, np.int64)
        self.right = np.zeros(self.n_lefs, np.int64)
        self.occ = np.zeros(n, np.int64)
        _seed(int(seed))
        for m in range(self.n_lefs):
            _place(self.left, self.right, self.occ, m, self.n)
        self.unbind_prob = 1.0 / self.lifetime

    def step(self, n_steps: int = 1) -> None:
        _steps(self.left, self.right, self.occ, self.p_left, self.p_right, self.unbind_prob, int(n_steps))

    def bonds(self) -> np.ndarray:
        """(n_lefs, 2) array of (left leg, right leg) monomer indices."""
        return np.stack([self.left, self.right], 1).copy()

    def loop_sizes(self) -> np.ndarray:
        return self.right - self.left


def synthetic_tad_boundaries(n: int, mean_tad_kb: float, rng: np.random.Generator, stall: float,
                             min_tad_kb: float = 150.0, margin: int = 100) -> list[dict]:
    """Domain borders with log-normally distributed spacing (human TADs: median ~ 0.5-1 Mb). Each border blocks LEFs
    arriving from both sides (a left-blocking and a right-blocking site at the same monomer)."""
    sites, x = [], margin
    while True:
        x += max(min_tad_kb, rng.lognormal(np.log(mean_tad_kb), 0.35))
        if x > n - margin:
            break
        sites.append({"pos": int(round(x)), "blocks_left": stall, "blocks_right": stall})
    return sites
