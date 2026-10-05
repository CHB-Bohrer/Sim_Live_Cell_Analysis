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


def _update_lef_bonds(sim, force, bonds, length, wiggle):
    k = sim.kbondScalingFactor / (wiggle * sim.length_scale) ** 2
    for i, (a, b) in enumerate(bonds):
        force.setBondParameters(i, int(a), int(b), length * sim.length_scale, float(k))
    force.updateParametersInContext(sim.context)


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

    lef, lef_force, ctcf = None, None, []
    if ex:
        ctcf = build_ctcf(cfg, N)
        speed = ex.get("lef_steps_per_block", 1)
        lifetime_steps = ex["processivity_kb"] / 2.0       # processivity = 2 * speed(=1/step) * lifetime
        n_lefs = int(round(N / ex["separation_kb"]))
        lef = LEF1D(N, n_lefs, lifetime_steps, ctcf=ctcf, seed=seed + 1)
        lef.step(int(ex.get("warmup_lifetimes", 10) * lifetime_steps))   # reach the steady state in 1D first
        lef_force = forces.harmonic_bonds(sim, bonds=lef.bonds(), bondWiggleDistance=ex.get("bond_wiggle", 0.2),
                                          bondLength=ex.get("bond_length", 0.5), name="lef_bonds",
                                          override_checks=True)
        sim.add_force(lef_force)
        import pandas as pd
        pd.DataFrame(ctcf).to_csv(out / "ctcf_sites.csv", index=False)

    progress(f"minimizing energy (N={N}, density={P['density']}, confinement radius {radius:.1f})")
    sim.local_energy_minimization()

    eq, prod = int(D["equilibration_blocks"]), int(D["production_blocks"])
    stride = int(D.get("track_stride", 10))
    tracked = np.arange(0, N, stride)
    traj = np.zeros((prod, len(tracked), 3), np.float32)
    lef_pos = np.zeros((prod, lef.n_lefs, 2), np.int32) if lef else None
    steps = int(D["steps_per_block"])
    full_every = int(D.get("full_every_blocks", 10))
    t0 = time.time()
    for b in range(eq + prod):
        if lef:
            lef.step(ex.get("lef_steps_per_block", 1))
            _update_lef_bonds(sim, lef_force, lef.bonds(), ex.get("bond_length", 0.5), ex.get("bond_wiggle", 0.2))
        p = b - eq
        sim.do_block(steps, save=(p >= 0 and p % full_every == 0))
        if p >= 0:
            traj[p] = sim.get_data()[tracked]
            if lef:
                lef_pos[p] = lef.bonds()
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
