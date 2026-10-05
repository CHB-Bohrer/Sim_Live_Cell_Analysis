"""Stage 2: multiple cells with persistent IDs and a motion model.

Cells are discs. Motion per frame = diffusive step (D) + persistent directed drift (speed, with heading
diffusion), hard-core collisions (overlaps are pushed apart), reflecting walls, optional division.
Output `cells.csv`: one row per cell per frame (physical units: um, s).
Division follows the Cell Tracking Challenge convention: the parent track ends the frame before, and two new
daughter IDs (with `parent_id`) start at the division frame.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def simulate_cells(cfg: dict, rng: np.random.Generator) -> pd.DataFrame:
    m, g = cfg["motion"], cfg["geometry"]
    n_frames, dt = cfg["acquisition"]["n_frames"], cfg["acquisition"]["frame_interval_s"]
    fov = np.array(g["fov_um"], float)  # (y, x)
    sigma = np.sqrt(2 * m["D_um2_s"] * dt)  # per-axis diffusive step (um)
    drift = m.get("drift_um_s", 0.0) * dt  # directed step length per frame (um)

    n0 = cfg["cells"]["n_cells"]
    radius = np.full(n0, g["cell_radius_um"], float)
    pos = np.zeros((n0, 2))
    for i in range(n0):  # random non-overlapping initial placement
        for _ in range(1000):
            p = rng.uniform(radius[i], fov - radius[i])
            if all(np.hypot(*(p - pos[j])) >= radius[i] + radius[j] for j in range(i)):
                pos[i] = p
                break
        else:
            raise RuntimeError("Could not place cells without overlap; reduce n_cells or radius")
    heading = rng.uniform(0, 2 * np.pi, n0)
    ids = list(range(1, n0 + 1))
    parent = [0] * n0
    next_id = n0 + 1
    rows = []

    for t in range(n_frames):
        for k, cid in enumerate(ids):
            rows.append((t, cid, pos[k, 0], pos[k, 1], radius[k], parent[k]))
        if t == n_frames - 1:
            break
        # --- division (daughters appear next frame) ---
        p_div, max_cells = cfg["cells"].get("p_divide_per_frame", 0.0), cfg["cells"].get("max_cells", 10**9)
        dividing, born = set(), []  # born: (pos, radius, heading, parent_id)
        if p_div > 0:
            for k, cid in enumerate(ids):
                if len(ids) - len(dividing) + len(born) < max_cells and rng.random() < p_div:
                    dividing.add(k)
                    r_d = radius[k] / np.sqrt(2)
                    ang = rng.uniform(0, np.pi)
                    off = r_d * np.array([np.sin(ang), np.cos(ang)])
                    for sgn in (+1, -1):
                        born.append((pos[k] + sgn * off, r_d, rng.uniform(0, 2 * np.pi), cid))
        if dividing:
            keep = [k for k in range(len(ids)) if k not in dividing]
            ids = [ids[k] for k in keep] + list(range(next_id, next_id + len(born)))
            next_id += len(born)
            pos = np.vstack([pos[keep]] + [b[0][None] for b in born])
            radius = np.concatenate([radius[keep], [b[1] for b in born]])
            heading = np.concatenate([heading[keep], [b[2] for b in born]])
            parent = [parent[k] for k in keep] + [b[3] for b in born]
        # --- motion ---
        heading = heading + rng.normal(0, m.get("heading_diffusion_rad", 0.3), len(ids))
        step = rng.normal(0, sigma, pos.shape) + drift * np.stack([np.sin(heading), np.cos(heading)], 1)
        pos = pos + step
        _resolve_collisions(pos, radius, fov, rng)

    return pd.DataFrame(rows, columns=["t", "cell_id", "y_um", "x_um", "radius_um", "parent_id"])


def _resolve_collisions(pos: np.ndarray, radius: np.ndarray, fov: np.ndarray, rng, iters: int = 30) -> None:
    n = len(pos)
    for _ in range(iters):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                d = pos[j] - pos[i]
                dist = np.hypot(*d)
                gap = radius[i] + radius[j] - dist
                if gap > 0:
                    u = d / dist if dist > 1e-9 else rng.normal(size=2)
                    u = u / np.hypot(*u)
                    pos[i] -= u * gap / 2; pos[j] += u * gap / 2
                    moved = True
        np.clip(pos, radius[:, None], fov - radius[:, None], out=pos)  # reflecting walls
        if not moved:
            break
