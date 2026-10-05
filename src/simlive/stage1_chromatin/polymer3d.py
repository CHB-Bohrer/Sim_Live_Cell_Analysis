"""Stage 1: one chromatin simulation = a confined polymer (polychrom / OpenMM on the GPU), optionally driven by
loop-extruding factors with boundary elements (1D dynamics in lef1d.py), following the polychrom loopExtrusion example
from the Mirny lab (variable Langevin, harmonic bonds, stiffness 1.5, soft-core repulsion, LEF bonds of rest length 0.5).

Units are polychrom's: 1 length unit = one monomer diameter, time = MD steps. With 1 monomer = 1 kb here. Converting to
nm and seconds is the calibration step (match locus MSD to measurements) and is NOT done in this file.

Output folder (one per simulation):
  meta.json                config, seed, timings, versions, LEF statistics
  ctcf_sites.csv           boundary elements used
  tracked_positions.npy    (n_production_blocks, n_tracked, 3) float32: every `track_stride`-th monomer, EVERY block
  tracked_beads.npy        which monomers those are
  lef_positions.npy        (n_production_blocks, n_lefs, 2) int32: (left leg, right leg) of every LEF, every block
  blocks/                  polychrom HDF5 snapshots of the whole chain every `full_every_blocks` production blocks
  final_conformation.npy   last full conformation (float32)
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from simlive.io.runs import write_provenance
from simlive.stage1_chromatin.lef1d import LEF1D, synthetic_tad_boundaries


def build_ctcf(cfg: dict, n: int) -> list[dict]:
    """Boundary elements from config: explicit list, or synthetic TAD boundaries (seeded, so every simulation of the
    same configuration shares the same 'genome')."""
    ex = cfg.get("extrusion") or {}
    c = ex.get("ctcf")
    if not c:
        return []
    if isinstance(c, list):
        return [dict(s) for s in c]
    if c.get("mode") == "tads":
        return synthetic_tad_boundaries(n, c["mean_tad_kb"], np.random.default_rng(c.get("seed", 0)), c["stall"])
    raise ValueError(f"unknown ctcf spec: {c}")


class LEFBondUpdater:
    """Loop-extruder bonds that really move from block to block.

    IMPORTANT (a bug in an earlier version of this file): OpenMM's `updateParametersInContext` can change a bond's
    length and stiffness but NOT which two particles it connects. Re-pointing existing bonds at new monomers is silently
    ignored, which froze the loops at their starting positions. The correct method (the one in the polychrom
    loopExtrusion example) is to register EVERY bond that will ever occur up front, with stiffness 0, and then switch
    them on and off by changing only the stiffness.

    lef_traj: (n_blocks, n_lefs, 2) int array of (left leg, right leg) for every MD block of the whole run.
    """

    def __init__(self, sim, lef_traj: np.ndarray, length: float = 0.5, wiggle: float = 0.2):
        from polychrom import forces

        flat = lef_traj.reshape(-1, 2)
        self.pairs, inv = np.unique(flat, axis=0, return_inverse=True)
        self.index = inv.reshape(lef_traj.shape[:2])                  # bond number active for each (block, LEF)
        self.force = forces.harmonic_bonds(sim, bonds=[(int(a), int(b)) for a, b in self.pairs],
                                           bondWiggleDistance=0, bondLength=length, name="lef_bonds",
                                           override_checks=True)       # wiggle 0 -> stiffness 0 (all inactive)
        self.sim = sim
        self.length = length * sim.length_scale
        self.k_active = float(sim.kbondScalingFactor / (wiggle * sim.length_scale) ** 2)
        self._active = np.array([], dtype=int)

    def step(self, block: int) -> None:
        """Make exactly the bonds of this block active (call before sim.do_block)."""
        cur = np.unique(self.index[block])
        for i in np.setdiff1d(self._active, cur):
            a, b = self.pairs[i]
            self.force.setBondParameters(int(i), int(a), int(b), self.length, 0.0)
        for i in np.setdiff1d(cur, self._active):
            a, b = self.pairs[i]
            self.force.setBondParameters(int(i), int(a), int(b), self.length, self.k_active)
        self.force.updateParametersInContext(self.sim.context)
        self._active = cur


def run_chromatin_simulation(cfg: dict, out_dir: str | Path, seed: int, progress=print) -> dict:
    from polychrom import forcekits, forces, simulation, starting_conformations
    from polychrom.hdf5_format import HDF5Reporter

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    P, D = cfg["polymer"], cfg["dynamics"]
    N = int(P["n_monomers"])
    ex = cfg.get("extrusion")
    t_start = time.time()

    reporter = HDF5Reporter(folder=str(out / "blocks"), max_data_length=50, overwrite=True)
    sim = simulation.Simulation(platform="CUDA", integrator="variableLangevin", error_tol=0.01, GPU="0",
                                collision_rate=D.get("collision_rate", 0.03), N=N, PBCbox=False,
                                reporters=[reporter], precision="mixed")
    sim.integrator.setRandomNumberSeed(int(seed) % (2 ** 31))      # thermal noise
    np.random.seed(int(seed) % (2 ** 32))                            # starting conformation + tiny initial offsets
    radius = (3 * N / (4 * np.pi * P["density"])) ** (1 / 3)
    box = int(2 * radius / np.sqrt(3) * 0.98)                        # the whole starting cube fits inside the sphere
    sim.set_data(starting_conformations.grow_cubic(N, box), center=True)
    sim.add_force(forces.spherical_confinement(sim, density=P["density"], k=5.0))
    sim.add_force(forcekits.polymer_chains(
        sim, chains=[(0, N, False)],
        bond_force_func=forces.harmonic_bonds,
        bond_force_kwargs={"bondLength": 1.0, "bondWiggleDistance": P.get("bond_wiggle", 0.1)},
        angle_force_func=forces.angle_force, angle_force_kwargs={"k": P.get("angle_k", 1.5)},
        nonbonded_force_func=forces.polynomial_repulsive,
        nonbonded_force_kwargs={"trunc": P.get("repulsion_trunc", 1.5), "radiusMult": P.get("radius_mult", 1.05)},
        except_bonds=True))

    eq, prod = int(D["equilibration_blocks"]), int(D["production_blocks"])
    lef, updater, lef_traj, ctcf = None, None, None, []
    if ex:
        ctcf = build_ctcf(cfg, N)
        lifetime_steps = ex["processivity_kb"] / 2.0       # processivity = 2 * speed(=1/step) * lifetime
        n_lefs = int(round(N / ex["separation_kb"]))
        lef = LEF1D(N, n_lefs, lifetime_steps, ctcf=ctcf, seed=seed + 1)
        lef.step(int(ex.get("warmup_lifetimes", 10) * lifetime_steps))   # reach the steady state in 1D first
        lef_traj = np.zeros((eq + prod, n_lefs, 2), np.int32)            # the whole 1D history, computed up front
        for b in range(eq + prod):
            lef.step(ex.get("lef_steps_per_block", 1))
            lef_traj[b] = lef.bonds()
        updater = LEFBondUpdater(sim, lef_traj, ex.get("bond_length", 0.5), ex.get("bond_wiggle", 0.2))
        sim.add_force(updater.force)
        import pandas as pd
        pd.DataFrame(ctcf).to_csv(out / "ctcf_sites.csv", index=False)
        progress(f"loop extrusion: {n_lefs} LEFs, {len(updater.pairs):,} distinct bonds over {eq + prod} blocks")

    progress(f"minimizing energy (N={N}, density={P['density']}, confinement radius {radius:.1f})")
    sim.local_energy_minimization()

    stride = int(D.get("track_stride", 10))
    tracked = np.arange(0, N, stride)
    traj = np.zeros((prod, len(tracked), 3), np.float32)
    lef_pos = lef_traj[eq:] if lef else None
    steps = int(D["steps_per_block"])
    full_every = int(D.get("full_every_blocks", 10))
    t0 = time.time()
    for b in range(eq + prod):
        if updater:
            updater.step(b)
        p = b - eq
        sim.do_block(steps, save=(p >= 0 and p % full_every == 0))
        if p >= 0:
            traj[p] = sim.get_data()[tracked]
        if b % 100 == 0:
            progress(f"block {b}/{eq + prod}  ({(b + 1) * steps / (time.time() - t0):.0f} MD steps/s)")
    reporter.dump_data()
    wall = time.time() - t0

    np.save(out / "tracked_positions.npy", traj)
    np.save(out / "tracked_beads.npy", tracked)
    if lef:
        np.save(out / "lef_positions.npy", lef_pos)
    final = np.asarray(sim.get_data(), np.float32)
    np.save(out / "final_conformation.npy", final)

    meta = {"seed": int(seed), "n_monomers": N, "confinement_radius": float(radius), "n_tracked": int(len(tracked)),
            "md_steps_total": (eq + prod) * steps, "md_steps_per_second": (eq + prod) * steps / wall,
            "wall_seconds_md": wall, "wall_seconds_total": time.time() - t_start, "n_ctcf_sites": len(ctcf)}
    if lef:
        sizes = (lef_pos[:, :, 1] - lef_pos[:, :, 0]).ravel()
        meta.update({"lef_count": lef.n_lefs, "mean_loop_size_kb": float(sizes.mean()),
                     "processivity_kb": ex["processivity_kb"], "separation_kb": ex["separation_kb"]})
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    write_provenance(out, cfg, int(seed))
    return meta
