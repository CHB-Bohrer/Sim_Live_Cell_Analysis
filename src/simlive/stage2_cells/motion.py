"""Stage 2: multiple cells with persistent IDs, a motion model and time-varying nuclear shape.

Motion per frame = diffusive step (D) + persistent directed drift (speed, with heading diffusion), hard-core
collisions between bounding circles, reflecting walls, optional division.

Shape: each nucleus is an ellipse (aspect ratio, orientation) whose boundary is perturbed by low-order Fourier
harmonics (bumps and dents):  rho(theta) = 1 + sum_k amp_k cos(k theta + phase_k)  in the ellipse's own frame.
Aspect, orientation, amplitudes and phases evolve slowly in time (AR(1) / random walk), so nuclei rotate and
deform. `radius_um` is the equivalent-circle radius; `bound_radius_um` is the true max extent (used for collisions).

Output `cells.csv`: one row per cell per frame (um, s). Division follows the Cell Tracking Challenge convention:
the parent track ends the frame before; two new daughter IDs (with `parent_id`) start at the division frame.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

_THETA = np.linspace(0, 2 * np.pi, 72, endpoint=False)


def boundary_rho(theta, amp, phase, harmonics):
    """Relative boundary radius at angle(s) theta (cell frame). amp/phase: (K,), harmonics: (K,)."""
    th = np.asarray(theta)[..., None]
    rho = 1.0 + (amp * np.cos(harmonics * th + phase)).sum(-1)
    return np.clip(rho, 0.4, None)


def _bound_radius(radius, aspect, amp, phase, harmonics):
    """Max distance from the centre to the boundary, per cell (n,)."""
    out = np.empty(len(radius))
    a_ax, b_ax = np.sqrt(aspect), 1 / np.sqrt(aspect)
    for i in range(len(radius)):
        rho = boundary_rho(_THETA, amp[i], phase[i], harmonics)
        out[i] = radius[i] * np.max(rho * np.hypot(a_ax[i] * np.cos(_THETA), b_ax[i] * np.sin(_THETA)))
    return out


def _new_shapes(n, sh, rng):
    k = len(sh["harmonics"])
    aspect = np.maximum(1.0, rng.normal(sh["aspect_mean"], sh["aspect_std"], n))
    return (aspect, rng.uniform(0, np.pi, n), rng.normal(0, sh["deform_amp"], (n, k)),
            rng.uniform(0, 2 * np.pi, (n, k)))


def simulate_cells(cfg: dict, rng: np.random.Generator) -> pd.DataFrame:
    m, g = cfg["motion"], cfg["geometry"]
    sh = g["shape"]
    harm = np.array(sh["harmonics"], float)
    n_frames, dt = cfg["acquisition"]["n_frames"], cfg["acquisition"]["frame_interval_s"]
    fov = np.array(g["fov_um"], float)  # (y, x)
    sigma = np.sqrt(2 * m["D_um2_s"] * dt)
    drift = m.get("drift_um_s", 0.0) * dt
    rho_p = sh.get("persistence", 0.9)  # AR(1) memory of deformation per frame

    n0 = cfg["cells"]["n_cells"]
    radius = np.full(n0, g["cell_radius_um"], float)
    aspect, angle, amp, phase = _new_shapes(n0, sh, rng)
    bound = _bound_radius(radius, aspect, amp, phase, harm)
    pos = np.zeros((n0, 2))
    for i in range(n0):  # random non-overlapping initial placement
        for _ in range(2000):
            p = rng.uniform(bound[i], fov - bound[i])
            if all(np.hypot(*(p - pos[j])) >= bound[i] + bound[j] for j in range(i)):
                pos[i] = p
                break
        else:
            raise RuntimeError("Could not place cells without overlap; reduce n_cells or radius")
    heading = rng.uniform(0, 2 * np.pi, n0)
    ids = list(range(1, n0 + 1))
    parent = [0] * n0
    next_id = n0 + 1
    rows = []
    cols = (["t", "cell_id", "y_um", "x_um", "radius_um", "aspect", "angle_rad", "bound_radius_um", "parent_id"]
            + [f"amp{int(k)}" for k in harm] + [f"phase{int(k)}" for k in harm])

    for t in range(n_frames):
        for k, cid in enumerate(ids):
            rows.append((t, cid, pos[k, 0], pos[k, 1], radius[k], aspect[k], angle[k], bound[k], parent[k],
                         *amp[k], *phase[k]))
        if t == n_frames - 1:
            break
        # --- division (daughters appear next frame) ---
        p_div, max_cells = cfg["cells"].get("p_divide_per_frame", 0.0), cfg["cells"].get("max_cells", 10**9)
        dividing, born = set(), []  # born: (pos, radius, heading, parent_id)
        if p_div > 0:
            for k, cid in enumerate(ids):
                if len(ids) + len(dividing) < max_cells and rng.random() < p_div:
                    dividing.add(k)
                    r_d = radius[k] / np.sqrt(2)
                    ang = angle[k]  # daughters separate along the long axis
                    off = 0.5 * bound[k] * np.array([np.sin(ang), np.cos(ang)])
                    for sgn in (+1, -1):
                        born.append((pos[k] + sgn * off, r_d, rng.uniform(0, 2 * np.pi), cid))
        if dividing:
            keep = [k for k in range(len(ids)) if k not in dividing]
            nb = len(born)
            ids = [ids[k] for k in keep] + list(range(next_id, next_id + nb))
            next_id += nb
            pos = np.vstack([pos[keep]] + [b[0][None] for b in born])
            radius = np.concatenate([radius[keep], [b[1] for b in born]])
            heading = np.concatenate([heading[keep], [b[2] for b in born]])
            parent = [parent[k] for k in keep] + [b[3] for b in born]
            a2, an2, am2, ph2 = _new_shapes(nb, sh, rng)
            aspect, angle = np.concatenate([aspect[keep], a2]), np.concatenate([angle[keep], an2])
            amp, phase = np.vstack([amp[keep], am2]), np.vstack([phase[keep], ph2])
        # --- shape dynamics: rotation, aspect and bump amplitudes drift slowly (stationary AR(1)) ---
        n = len(ids)
        angle = angle + rng.normal(0, sh.get("rot_diffusion_rad", 0.1), n)
        aspect = np.maximum(1.0, sh["aspect_mean"] + rho_p * (aspect - sh["aspect_mean"])
                            + np.sqrt(1 - rho_p**2) * rng.normal(0, sh["aspect_std"], n))
        amp = rho_p * amp + np.sqrt(1 - rho_p**2) * rng.normal(0, sh["deform_amp"], amp.shape)
        phase = phase + rng.normal(0, 0.15, phase.shape)
        bound = _bound_radius(radius, aspect, amp, phase, harm)
        # --- motion ---
        heading = heading + rng.normal(0, m.get("heading_diffusion_rad", 0.3), n)
        pos = pos + rng.normal(0, sigma, pos.shape) + drift * np.stack([np.sin(heading), np.cos(heading)], 1)
        _resolve_collisions(pos, bound, fov, rng)

    return pd.DataFrame(rows, columns=cols)


def _resolve_collisions(pos: np.ndarray, radius: np.ndarray, fov: np.ndarray, rng, iters: int = 100) -> None:
    n = len(pos)
    for _ in range(iters):
        np.clip(pos, radius[:, None], fov - radius[:, None], out=pos)  # reflecting walls (before the overlap
        moved = False                                                  # pass, so a clip can't leave an overlap)
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
        if not moved:
            break
