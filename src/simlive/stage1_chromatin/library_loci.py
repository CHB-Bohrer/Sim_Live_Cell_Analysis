"""Stage 1 -> movie: loci taken from saved polymer simulations (replaces `toy_loci.py` when `loci.source: library`).

Each cell gets its own simulation (`library.assign_simulations`: random draw from the library, optionally restricted to
`loci.sim_ids`, with `loci.pin: {cell_id: sim_id}` forcing particular cells onto particular simulations). Its loci are monomers at `positions_kb`; their polymer
coordinates are mapped into the nucleus frame (so they rotate / deform / move with the nucleus, as in `toy_loci`):

  - space: polymer coordinates are converted to TRUE distances with one calibration number, `nm_per_unit` (nm per polymer
    length unit); the region keeps whatever physical size that gives (no rescaling to the nucleus). The polymer gets a random
    rotation per cell and is placed at a random territory centre inside the nucleus (fixed in the nucleus frame, drawn per
    cell). A position whose loci would leave the nucleus (0.7 of its radius) in any frame is rejected and a new one is
    drawn; an error is raised only if no position works.
    (x, y) is the 2D image-plane position; z is kept as `zn` for later 3D work. `nm_per_unit` is a PLACEHOLDER until the
    MSD / compaction calibration exists;
  - time: one saved block = `block_duration_s` seconds. PLACEHOLDER until the MSD calibration exists (the user will supply
    measured MSD); frame k reads block `offset + k * frame_interval_s / block_duration_s` (nearest block), where `offset`
    is random per cell so cells do not all start at the same polymer state.

Output has the same columns as `toy_loci` (xn, yn, zn in units of the nuclear radius) plus `sim_id` and `zn`: t, cell_id, locus_id, xn, yn, zn, y_um, x_um, sim_id.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from simlive.stage1_chromatin import library as lib
from simlive.stage1_chromatin.toy_loci import R_MAX
from simlive.stage2_cells.motion import nucleus_to_lab_um


def _territory_centre(p: np.ndarray, rng: np.random.Generator, cell_id, tries: int = 1000) -> np.ndarray:
    """Random centre (nuclear radii) such that every locus stays inside R_MAX of the nuclear centre in every frame;
    candidates that do not fit are rejected and a new position is drawn. Only a region whose loci cannot fit anywhere
    (excursion larger than the nucleus) fails, after `tries` draws."""
    for _ in range(tries):
        c = rng.normal(size=3); c *= R_MAX * rng.random() ** (1 / 3) / np.linalg.norm(c)
        if np.linalg.norm(p + c, axis=-1).max() <= R_MAX:
            return c
    raise ValueError(f"cell {cell_id}: the loci' excursion ({np.linalg.norm(p, axis=-1).max():.2f} nuclear radii from "
                     f"the region centre) cannot fit inside the nucleus; lower nm_per_unit or pick loci nearer the centre")


def simulate_library_loci(cells: pd.DataFrame, cfg: dict, rng: np.random.Generator) -> pd.DataFrame:
    lc = cfg["loci"]
    positions = list(lc["positions_kb"])
    nm_per_unit, block_s = float(lc["nm_per_unit"]), float(lc["block_duration_s"])
    dt = float(cfg["acquisition"]["frame_interval_s"])
    cell_ids = sorted(cells.cell_id.unique())
    assigned = lib.assign_simulations(lc["library"], cell_ids, rng, replace=bool(lc.get("replace", False)),
                                      pool=lc.get("sim_ids"), pinned=lc.get("pin"))
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
        nuc_um = float(g.radius_um.min())                                       # smallest nuclear radius of this cell
        p = Rotation.random(random_state=int(rng.integers(2**31))).apply(
            traj[offset + steps].reshape(-1, 3)).reshape(len(steps), -1, 3) * (nm_per_unit / 1000.0 / nuc_um)  # nuclear radii
        p += _territory_centre(p, rng, cid)
        for i, r in enumerate(g.itertuples()):
            y, x = nucleus_to_lab_um(r, p[i, :, 0], p[i, :, 1])
            for k in range(p.shape[1]):
                rows.append((r.t, cid, k, p[i, k, 0], p[i, k, 1], p[i, k, 2], y[k], x[k], assigned[int(cid)]))
    return pd.DataFrame(rows, columns=["t", "cell_id", "locus_id", "xn", "yn", "zn", "y_um", "x_um", "sim_id"])
