"""Stage 1 -> movie: loci taken from saved polymer simulations (replaces `toy_loci.py` when `loci.source: library`).

Each cell gets its own simulation (`library.assign_simulations`: random draw from the library, optionally restricted to
`loci.sim_ids`, with `loci.pin: {cell_id: sim_id}` forcing particular cells onto particular simulations). Its loci are monomers at `positions_kb`; their polymer
coordinates are mapped into the nucleus frame (so they rotate / deform / move with the nucleus, as in `toy_loci`):

  - space: polymer coordinates are converted to TRUE distances with one calibration number, `nm_per_unit` (nm per polymer
    length unit; both scale numbers come from `loci.calibration` = a file in configs/calibration/, see simlive.calibration.polymer, unless the
    config gives `nm_per_unit` / `block_duration_s` itself); the region keeps whatever physical size that gives (no rescaling to the nucleus). The polymer gets a random
    rotation per cell and is placed at a random territory centre inside the nucleus (fixed in the nucleus frame, drawn per
    cell). Only the territory centre deforms with the nucleus (it is stretched with the nucleus' aspect ratio and half thickness);
    the polymer around it keeps its true size and shape and turns with the nucleus (`rigid_factors`; before 2026-10-26 the whole polymer
    was stretched, so drawn promoter-enhancer distances were 0.89-1.14 x the polymer's, 5-95%). A position whose loci would leave the
    nucleus (0.7 of its radius) in any frame is rejected and a new one is drawn; an error is raised only if no position works.
    (x, y) is the 2D image-plane position; z is kept as `zn` for later 3D work. `nm_per_unit` is only as good as the calibration's
    assumptions (contact fraction and working distance);
  - time: one saved block = `block_duration_s` seconds (from the MSD calibration: an ASSUMED single-locus D until the user supplies a measured MSD); frame k reads polymer time `offset + k * frame_interval_s / block_duration_s` (blocks, interpolated), where
    `offset` is random per cell so cells do not all start at the same polymer state.

Output has the same columns as `toy_loci` (xn, yn, zn in units of the nuclear radius) plus `sim_id` and `zn`: t, cell_id, locus_id, copy, xn, yn, zn,
y_um, x_um, sim_id. `copy` numbers the promoter-enhancer pairs of a cell (`loci.copies`; 0 only, unless copies are switched on).
"""
from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd

from simlive.calibration.polymer import scale_from_config
from simlive.stage1_chromatin import library as lib
from simlive.stage1_chromatin.copies import copy_settings, draw_copy_numbers, extra_centre
from simlive.stage1_chromatin.rotation import random_rotation
from simlive.stage1_chromatin.toy_loci import R_MAX
from simlive.stage1_chromatin.transcription import pingpong
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


def rigid_factors(g: pd.DataFrame, radius_um: float, three_d: bool) -> np.ndarray:
    """(frames, 3) factors that turn polymer offsets (in units of `radius_um`) into nucleus-frame offsets that `nucleus_to_lab_um` draws at their
    TRUE size: that mapping stretches the nucleus frame by the nucleus' shape (x by radius x sqrt(aspect), y by radius / sqrt(aspect), z by the
    half thickness z_semi_um), which is right for the territory centre (it moves with the deforming nucleus) but would stretch a promoter-enhancer
    pair by -18% / +24% (aspect 1.0-1.57; found 2026-10-26) away from the polymer distance that drives transcription. g: the cell's cells.csv rows."""
    r = g.radius_um.to_numpy(float)
    s = np.sqrt(g.aspect.to_numpy(float))
    fz = radius_um / g.z_semi_um.to_numpy(float) if three_d else radius_um / r
    return np.stack([radius_um / (r * s), radius_um * s / r, fz], axis=1)


def draw_pair(rng: np.random.Generator, n_monomers: int, beads: np.ndarray, separation_kb: float = 100.0, margin_kb: float = 1000.0) -> list[int]:
    """Genomic positions (kb = monomer index) of [promoter, enhancer] for one cell: a pair `separation_kb` apart, on the stored monomer grid,
    anywhere at least `margin_kb` from both ends of the region (loop extrusion and boundaries behave differently at a chain end), the promoter
    on either side with equal probability. With every cell drawing its own pair, a library of 100 simulations gives hundreds of different pairs."""
    stride = int(beads[1] - beads[0])
    sep = int(round(separation_kb / stride)) * stride
    lo = -(-int(margin_kb) // stride) * stride                       # ceil to the grid
    hi = (int(n_monomers) - 1 - int(margin_kb)) // stride * stride - sep   # last allowed left position
    if hi < lo:
        raise ValueError(f"a {separation_kb:g} kb pair does not fit {margin_kb:g} kb from both ends of a {n_monomers} kb region")
    left = lo + stride * int(rng.integers(0, (hi - lo) // stride + 1))
    pair = [left, left + sep]
    if rng.random() < 0.5:
        pair.reverse()
    return pair


def simulate_library_loci(cells: pd.DataFrame, cfg: dict, rng: np.random.Generator, return_traces: bool = False,
                          rng_copies: np.random.Generator | None = None):
    """Loci of every cell from its own saved simulation. Frame k is read at polymer time k x frame_interval_s (any frame interval
    works: positions are interpolated between saved blocks). `loci.positions_kb` is a list [promoter, enhancer] (the same pair in every cell) or the
    word "random": every cell then draws its own pair (`draw_pair`: `loci.pair_separation_kb` default 100, `loci.pair_margin_kb` default 1000); the
    pairs are in the result's `df.attrs["pairs_kb"]` = {cell_id: [promoter_kb, enhancer_kb]} (copy 0) and `df.attrs["pairs"]` = rows
    (cell_id, copy, sim_id, promoter_kb, enhancer_kb) of every copy. If the movie is longer than the saved trajectory the trajectory is
    repeated forwards/backwards (no jump at the seams) and a warning is raised.
    Copies (`loci.copies`): copy 0 is exactly the one-copy result; every further copy of a cell gets another simulation of the library (a
    different chromosome, so a different trajectory; drawn from `rng_copies`), its own pair, rotation, time offset and territory, at least
    `min_separation_um` from the cell's other copies in the image plane.
    With return_traces=True also returns {cell_id: {"dt_block_s", "d_nm", "looped"}} (copy 0; key (cell_id, copy) for further copies): the 3D
    distance between locus 0 (promoter) and locus 1 (enhancer) at full polymer-block resolution, starting at the cell's first frame (for the
    transcription model)."""
    lc = cfg["loci"]
    random_pair = isinstance(lc["positions_kb"], str)
    if random_pair and lc["positions_kb"] != "random":
        raise ValueError(f"loci.positions_kb must be a list [promoter_kb, enhancer_kb] or 'random', not {lc['positions_kb']!r}")
    positions = None if random_pair else list(lc["positions_kb"])
    pairs, pair_rows = {}, []
    nm_per_unit, block_s, scale_source = scale_from_config(lc)
    dt = float(cfg["acquisition"]["frame_interval_s"])
    cell_ids = sorted(cells.cell_id.unique())
    assigned = lib.assign_simulations(lc["library"], cell_ids, rng, replace=bool(lc.get("replace", False)),
                                      pool=lc.get("sim_ids"), pinned=lc.get("pin"))
    cs = copy_settings(cfg)
    n_copies = draw_copy_numbers(cell_ids, cfg, rng_copies)
    pool = sorted(str(x) for x in lc.get("sim_ids") or lib.list_library(lc["library"]).sim_id) if cs.active else []
    three_d = "z_um" in cells.columns
    rows, traces = [], {}

    def place(cid, g, copy, sim_name, r, others):
        """One copy of one cell: its polymer trajectory mapped into the nucleus. Returns (positions (frames, n_loci, 3) in nuclear radii,
        territory centre, beads used, base trajectory, looped)."""
        sim = lib.library_dir(lc["library"]) / sim_name
        meta = json.loads((sim / "meta.json").read_text())
        pos_kb = positions
        if random_pair:
            pos_kb = draw_pair(r, int(meta["n_monomers"]), np.load(sim / "tracked_beads.npy"),
                               float(lc.get("pair_separation_kb", 100.0)), float(lc.get("pair_margin_kb", 1000.0)))
        traj, beads_used = lib.locus_trajectories(sim, pos_kb)                 # (T, n_loci, 3) polymer units
        tau = (g.t.to_numpy() - g.t.min()) * dt / block_s                       # polymer time of each frame, in blocks
        need = int(np.ceil(tau.max())) + 2
        offset = int(r.integers(0, max(1, len(traj) - need)))
        base = traj[offset:]
        ext, looped = pingpong(base, need)
        if looped:
            warnings.warn(f"cell {cid}: the movie ({tau.max() * block_s:.0f} s) is longer than the saved trajectory "
                          f"({(len(base) - 1) * block_s:.0f} s from the chosen start); it is repeated back and forth. "
                          "Raise block_duration_s or use longer simulations.")
        i0 = np.floor(tau).astype(int)
        w = (tau - i0)[:, None, None]
        pos = (1 - w) * ext[i0] + w * ext[i0 + 1]                               # linear interpolation between blocks
        nuc_um = float(g.radius_um.min())                                       # smallest nuclear radius of this cell
        p = (pos[..., None, :] * random_rotation(r)).sum(-1) * (nm_per_unit / 1000.0 / nuc_um)   # nuclear radii; no BLAS call here
        p = p * rigid_factors(g, nuc_um, three_d)[:, None, :]                   # undo the nucleus' stretch for the polymer's own offsets
        if copy == 0:
            c = _territory_centre(p, r, cid)
        else:
            c = extra_centre(r, others, cs.min_separation_um / nuc_um, R_MAX, fits=lambda q: np.linalg.norm(p + q, axis=-1).max() <= R_MAX)
        return p + c, c, beads_used, base, looped

    for cid, g in cells.sort_values("t").groupby("cell_id"):
        sims, centres = [assigned[int(cid)]], []
        for copy in range(n_copies.get(int(cid), 1)):
            r = rng if copy == 0 else rng_copies
            if copy > 0:                        # another chromosome: a simulation this cell does not use yet (if the library has one)
                free = [s for s in pool if s not in sims] or pool
                sims.append(str(free[int(r.integers(len(free)))]))
            p, c, beads_used, base, looped = place(cid, g, copy, sims[copy], r, centres)
            centres.append(c)
            if copy == 0:
                pairs[int(cid)] = [int(b) for b in beads_used]
            pair_rows.append((int(cid), copy, sims[copy], *[int(b) for b in beads_used]))
            for i, rr in enumerate(g.itertuples()):
                y, x = nucleus_to_lab_um(rr, p[i, :, 0], p[i, :, 1])
                for k in range(p.shape[1]):
                    zk = rr.z_um + p[i, k, 2] * rr.z_semi_um if three_d else np.nan   # z-stacks: polymer z scaled into the (flattened) nucleus
                    rows.append((rr.t, cid, k, copy, p[i, k, 0], p[i, k, 1], p[i, k, 2], y[k], x[k], sims[copy], zk))
            if return_traces and base.shape[1] >= 2:
                traces[int(cid) if copy == 0 else (int(cid), copy)] = {
                    "dt_block_s": block_s, "looped": bool(looped), "d_nm": np.linalg.norm(base[:, 0] - base[:, 1], axis=1) * nm_per_unit}
    df = pd.DataFrame(rows, columns=["t", "cell_id", "locus_id", "copy", "xn", "yn", "zn", "y_um", "x_um", "sim_id", "z_um"])
    if not three_d:
        df = df.drop(columns="z_um")
    df.attrs["pairs_kb"] = pairs
    df.attrs["pairs"] = pair_rows
    df.attrs["scale"] = {"nm_per_unit": nm_per_unit, "block_duration_s": block_s, "source": scale_source}
    return (df, traces) if return_traces else df
