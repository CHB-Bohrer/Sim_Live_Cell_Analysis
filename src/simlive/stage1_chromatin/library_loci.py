"""Stage 1 -> movie: loci taken from saved polymer simulations (replaces `toy_loci.py` when `loci.source: library`).

Each cell gets its own simulation (`library.assign_simulations`). Its loci are monomers at `positions_kb`; their polymer
coordinates are mapped into the nucleus frame (so they rotate / deform / move with the nucleus, as in `toy_loci`):

  - space: the simulated region is only a small part of the genome (10 Mb of ~6 Gb diploid is ~1/600 of the DNA, i.e.
    ~0.12 of the nuclear radius if compact), so its confinement sphere (radius R, polymer units) is scaled to
    `region_radius` x nuclear radius and placed at a random territory centre inside the nucleus (fixed in the nucleus
    frame, drawn per cell). The polymer is given a random rotation per cell; (x, y) is the 2D image-plane position
    (z is kept as `zn` for later 3D work). PLACEHOLDER until the length calibration (nm per polymer unit) exists;
  - time: one saved block = `block_duration_s` seconds. PLACEHOLDER until the MSD calibration exists (the user will supply
    measured MSD); frame k reads block `offset + k * frame_interval_s / block_duration_s` (nearest block), where `offset`
    is random per cell so cells do not all start at the same polymer state.

Output has the same columns as `toy_loci` plus `sim_id` and `zn`: t, cell_id, locus_id, xn, yn, zn, y_um, x_um, sim_id.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from simlive.stage1_chromatin import library as lib
from simlive.stage1_chromatin.toy_loci import R_MAX
from simlive.stage2_cells.motion import nucleus_to_lab_um


def simulate_library_loci(cells: pd.DataFrame, cfg: dict, rng: np.random.Generator) -> pd.DataFrame:
    lc = cfg["loci"]
    positions = list(lc["positions_kb"])
    reg, block_s = float(lc.get("region_radius", 0.12)), float(lc["block_duration_s"])
    if not 0 < reg < R_MAX:
        raise ValueError(f"region_radius must be in (0, {R_MAX}) nuclear radii")
    dt = float(cfg["acquisition"]["frame_interval_s"])
    cell_ids = sorted(cells.cell_id.unique())
    assigned = lib.assign_simulations(lc["library"], cell_ids, rng, replace=bool(lc.get("replace", False)))
    rows = []
    for cid, g in cells.sort_values("t").groupby("cell_id"):
        sim = lib.library_dir(lc["library"]) / assigned[int(cid)]
        traj, _ = lib.locus_trajectories(sim, positions)                       # (T, n_loci, 3) polymer units
        radius = float(json.loads((sim / "meta.json").read_text())["confinement_radius"])
        steps = np.rint((g.t.to_numpy() - g.t.min()) * dt / block_s).astype(int)
        offset = int(rng.integers(0, max(1, len(traj) - steps.max())))
        if offset + steps.max() >= len(traj):
            raise ValueError(f"movie needs {steps.max() + 1} blocks of {assigned[int(cid)]} but it has {len(traj)}; "
                             "raise block_duration_s or simulate longer")
        centre = rng.normal(size=3); centre *= (R_MAX - reg) * rng.random() ** (1 / 3) / np.linalg.norm(centre)
        p = Rotation.random(random_state=int(rng.integers(2**31))).apply(
            traj[offset + steps].reshape(-1, 3)).reshape(len(steps), -1, 3) * (reg / radius) + centre
        for i, r in enumerate(g.itertuples()):
            y, x = nucleus_to_lab_um(r, p[i, :, 0], p[i, :, 1])
            for k in range(p.shape[1]):
                rows.append((r.t, cid, k, p[i, k, 0], p[i, k, 1], p[i, k, 2], y[k], x[k], assigned[int(cid)]))
    return pd.DataFrame(rows, columns=["t", "cell_id", "locus_id", "xn", "yn", "zn", "y_um", "x_um", "sim_id"])
