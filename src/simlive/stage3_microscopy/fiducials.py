"""Fiducial beads: a few bright sub-resolution beads lying on the slide in every image of the movie, emitting in every colour.

They are not part of the sample. They exist so that the analysis can measure, frame by frame, how far the whole picture has moved
(stage drift, focus error) and subtract it. They sit on the slide (`fiducials.z_um` above the coverslip), so a z-stack must reach
down to the slide for them to be seen (centre the stack low enough: see configs/loci_3d.yaml). Config block:

    fiducials: {enabled: true, n_beads: 8, photons: 20000, z_um: 0.3, min_distance_um: 3.0, min_separation_um: 8.0}

`make_fiducials` places them away from every nucleus (over the whole movie) so they do not sit on top of cells; the table is saved as
`stage3_microscopy/fiducials.csv` (bead, y_um, x_um, z_um) and is part of the ground truth.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

DEFAULTS = {"enabled": False, "n_beads": 8, "photons": 20000.0, "z_um": 0.3, "min_distance_um": 3.0, "min_separation_um": 8.0}


def fiducial_settings(cfg: dict) -> dict:
    return {**DEFAULTS, **(cfg.get("fiducials") or {})}


def make_fiducials(cfg: dict, cells: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame | None:
    """Positions of the fiducial beads (None if switched off): random, at least `min_distance_um` from the edge of every nucleus at every
    time, and at least `min_separation_um` from each other (as many as fit, up to n_beads)."""
    fc = fiducial_settings(cfg)
    if not fc["enabled"] or int(fc["n_beads"]) <= 0:
        return None
    fov = np.array(cfg["geometry"]["fov_um"], float)
    cy, cx, br = cells.y_um.to_numpy(), cells.x_um.to_numpy(), cells.bound_radius_um.to_numpy()
    pts: list = []
    for _ in range(20000):
        if len(pts) >= int(fc["n_beads"]):
            break
        p = rng.uniform(0.05 * fov, 0.95 * fov)
        if np.min(np.hypot(cy - p[0], cx - p[1]) - br) < fc["min_distance_um"]:
            continue
        if any(np.hypot(*(p - q)) < fc["min_separation_um"] for q in pts):
            continue
        pts.append(p)
    if not pts:
        return None
    pts = np.array(pts)
    z = float(fc["z_um"]) + rng.normal(0, 0.03, len(pts))
    return pd.DataFrame({"bead": np.arange(1, len(pts) + 1), "y_um": pts[:, 0], "x_um": pts[:, 1], "z_um": z})
