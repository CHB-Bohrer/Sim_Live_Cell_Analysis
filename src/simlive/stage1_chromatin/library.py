"""The chromatin simulation LIBRARY: many independent saved simulations that movies draw from.

Layout:  <DATA_ROOT>/chromatin/<library>/<name>_seed<seed>/   (one folder per simulation, see polymer3d.py)
         <DATA_ROOT>/chromatin/<library>/library.json         what the library is (config, seeds, versions)

A movie picks a library and gives every cell its own simulation (`assign_simulations`), then reads that simulation's
locus trajectories (`locus_trajectories`). Mapping the polymer's coordinates and time into a nucleus, nm and seconds is
a separate step (calibration), done when the microscopy labels are introduced.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from simlive.io.runs import DATA_ROOT


def library_dir(library: str) -> Path:
    return DATA_ROOT / "chromatin" / library


def list_library(library: str) -> pd.DataFrame:
    """One row per finished simulation in the library."""
    rows = []
    for d in sorted(library_dir(library).glob("*_seed*")):
        mp = d / "meta.json"
        if not mp.exists():
            continue
        m = json.loads(mp.read_text())
        rows.append({"sim_id": d.name, "seed": m["seed"], "n_monomers": m["n_monomers"],
                     "extrusion": "mean_loop_size_kb" in m, "mean_loop_kb": m.get("mean_loop_size_kb"),
                     "md_steps_per_s": m["md_steps_per_second"], "minutes": m["wall_seconds_total"] / 60, "path": str(d)})
    return pd.DataFrame(rows)


def assign_simulations(library: str, cell_ids, rng: np.random.Generator, replace: bool = False) -> dict[int, str]:
    """Give every cell its own simulation from the library (without replacement while simulations last).

    Cells are independent in reality, so each should have its own independent trajectory. If there are more cells than
    simulations, `replace=True` is required and some cells will share a trajectory (documented, not silent).
    """
    lib = list_library(library)
    cell_ids = list(cell_ids)
    if len(lib) == 0:
        raise ValueError(f"library '{library}' has no finished simulations")
    if len(cell_ids) > len(lib) and not replace:
        raise ValueError(f"{len(cell_ids)} cells but only {len(lib)} simulations; pass replace=True to reuse some")
    picks = rng.choice(lib.sim_id.to_numpy(), size=len(cell_ids), replace=replace)
    return {int(c): str(p) for c, p in zip(cell_ids, picks)}


def locus_trajectories(sim_path: str | Path, positions_kb) -> tuple[np.ndarray, np.ndarray]:
    """Polymer coordinates of loci at the given genomic positions (kb = monomer index), every block of the production
    run. Returns (traj (T, n_loci, 3), bead indices actually used). Only every `track_stride`-th monomer is stored
    each block, so each position snaps to the nearest stored monomer (within track_stride/2 kb)."""
    d = Path(sim_path)
    traj = np.load(d / "tracked_positions.npy", mmap_mode="r")
    beads = np.load(d / "tracked_beads.npy")
    idx = np.array([int(np.argmin(np.abs(beads - p))) for p in positions_kb])
    return np.asarray(traj[:, idx, :]), beads[idx]
