"""STAND-IN for stage 1 (chromatin dynamics): toy loci confined inside each nucleus.

Real stage 1 will be polychrom/OpenMM trajectories calibrated to nm and s. Until then each cell gets `n_loci`
point-like loci that diffuse in the nucleus' own frame as an Ornstein-Uhlenbeck process (confined diffusion:
stationary spread `confinement_std` in units of the nuclear radius, relaxation time `tau_s`). Because the position is
stored in the nucleus frame, loci rotate / deform / move with the cell, exactly the property the per-cell analysis
must handle. Output `loci_truth.csv`: t, cell_id, locus_id, xn, yn (nucleus frame), y_um, x_um (lab frame).
For z-stacks (cells.csv has `z_um`) the loci also move along z (columns zn, z_um) with the same statistics.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from simlive.stage2_cells.motion import nucleus_to_lab_um

R_MAX = 0.7  # loci are kept inside this fraction of the nuclear radius (nuclear boundary wobbles ~ +-25%)


def simulate_toy_loci(cells: pd.DataFrame, cfg: dict, rng: np.random.Generator) -> pd.DataFrame:
    lc = cfg["loci"]
    n, std, tau = int(lc["n_loci"]), float(lc["confinement_std"]), float(lc["tau_s"])
    a = np.exp(-cfg["acquisition"]["frame_interval_s"] / tau)
    three_d = "z_um" in cells.columns
    d = 3 if three_d else 2                                   # the 2D draws are unchanged (same random numbers as before z existed)
    rows = []
    for cid, g in cells.sort_values("t").groupby("cell_id"):
        pos = rng.normal(0, std, (n, d))
        for r in g.itertuples():
            pos = a * pos + np.sqrt(1 - a**2) * rng.normal(0, std, (n, d))
            rad = np.linalg.norm(pos, axis=1)
            pos = pos * (R_MAX * np.tanh(rad / R_MAX) / np.maximum(rad, 1e-9))[:, None]  # soft confinement
            y, x = nucleus_to_lab_um(r, pos[:, 0], pos[:, 1])
            for k in range(n):
                row = (r.t, cid, k, pos[k, 0], pos[k, 1], y[k], x[k])
                rows.append(row + ((pos[k, 2], r.z_um + pos[k, 2] * r.z_semi_um) if three_d else ()))
    cols = ["t", "cell_id", "locus_id", "xn", "yn", "y_um", "x_um"] + (["zn", "z_um"] if three_d else [])
    return pd.DataFrame(rows, columns=cols)
