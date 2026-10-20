"""STAND-IN for stage 1 (chromatin dynamics): toy loci confined inside each nucleus.

Real stage 1 will be polychrom/OpenMM trajectories calibrated to nm and s. Until then each cell gets `n_loci`
point-like loci that diffuse in the nucleus' own frame as an Ornstein-Uhlenbeck process (confined diffusion:
stationary spread `confinement_std` in units of the nuclear radius, relaxation time `tau_s`). Because the position is
stored in the nucleus frame, loci rotate / deform / move with the cell, exactly the property the per-cell analysis
must handle. Output `loci_truth.csv`: t, cell_id, locus_id, copy, xn, yn (nucleus frame), y_um, x_um (lab frame).
For z-stacks (cells.csv has `z_um`) the loci also move along z (columns zn, z_um) with the same statistics.

Copies (`loci.copies`, see copies.py): copy 0 sits around the nuclear centre and uses exactly the random numbers of a one-copy run; every
further copy is the same process around its own territory centre, at least `min_separation_um` from the others, drawn from `rng_copies`.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from simlive.stage1_chromatin.copies import copy_settings, draw_copy_numbers, extra_centre
from simlive.stage2_cells.motion import nucleus_to_lab_um

R_MAX = 0.7  # loci are kept inside this fraction of the nuclear radius (nuclear boundary wobbles ~ +-25%)


def _confine(q: np.ndarray, room: float = R_MAX) -> np.ndarray:
    rad = np.linalg.norm(q, axis=1)
    return q * (room * np.tanh(rad / room) / np.maximum(rad, 1e-9))[:, None]      # soft confinement within `room` of the centre


def simulate_toy_loci(cells: pd.DataFrame, cfg: dict, rng: np.random.Generator, rng_copies: np.random.Generator | None = None) -> pd.DataFrame:
    lc = cfg["loci"]
    n, std, tau = int(lc["n_loci"]), float(lc["confinement_std"]), float(lc["tau_s"])
    a = np.exp(-cfg["acquisition"]["frame_interval_s"] / tau)
    three_d = "z_um" in cells.columns
    d = 3 if three_d else 2                                   # the 2D draws are unchanged (same random numbers as before z existed)
    cs = copy_settings(cfg)
    n_copies = draw_copy_numbers(sorted(cells.cell_id.unique()), cfg, rng_copies)
    rows = []

    def emit(r, cid, copy, pos):
        y, x = nucleus_to_lab_um(r, pos[:, 0], pos[:, 1])
        for k in range(n):
            row = (r.t, cid, k, copy, pos[k, 0], pos[k, 1], y[k], x[k])
            rows.append(row + ((pos[k, 2], r.z_um + pos[k, 2] * r.z_semi_um) if three_d else ()))

    for cid, g in cells.sort_values("t").groupby("cell_id"):
        pos = rng.normal(0, std, (n, d))
        for r in g.itertuples():
            pos = _confine(a * pos + np.sqrt(1 - a**2) * rng.normal(0, std, (n, d)))     # the confined position feeds the next step
            emit(r, cid, 0, pos)
        if n_copies.get(int(cid), 1) > 1:              # further copies: own territory centre, own random stream
            rmin = float(g.radius_um.min())
            centres = [np.zeros(3)]
            for c in range(1, n_copies[int(cid)]):
                ctr = extra_centre(rng_copies, centres, cs.min_separation_um / rmin, R_MAX - 1.5 * std)[:d]
                centres.append(np.r_[ctr, np.zeros(3 - d)])
                rel = rng_copies.normal(0, std, (n, d))          # position relative to this copy's territory centre
                for r in g.itertuples():
                    # soft confinement around the copy's own centre, within the room left to the nuclear edge (so the centre stays put)
                    rel = _confine(a * rel + np.sqrt(1 - a**2) * rng_copies.normal(0, std, (n, d)), R_MAX - float(np.linalg.norm(ctr)))
                    emit(r, cid, c, ctr + rel)
    cols = ["t", "cell_id", "locus_id", "copy", "xn", "yn", "y_um", "x_um"] + (["zn", "z_um"] if three_d else [])
    return pd.DataFrame(rows, columns=cols)
