"""GPU regression test for loop-extruder bonds (run with `pytest -m gpu`).

History: an earlier version re-pointed existing OpenMM bonds at new monomers with updateParametersInContext, which
OpenMM silently ignores, so the loops stayed frozen at their starting positions for the entire simulation and the
first 'validation' measured frozen loops instead of dynamic extrusion. This test fails if bonds ever stop moving.
"""
import os
import shutil
import tempfile

import numpy as np
import pytest

pytestmark = pytest.mark.gpu


def test_moving_loop_bond_pulls_new_pair_together_and_releases_the_old_one():
    from polychrom import forcekits, forces, simulation, starting_conformations
    from polychrom.hdf5_format import HDF5Reporter
    from simlive.stage1_chromatin.polymer3d import LEFBondUpdater

    N = 1000
    out = tempfile.mkdtemp()
    try:
        sim = simulation.Simulation(platform="CUDA", integrator="variableLangevin", error_tol=0.01, GPU="0",
                                    collision_rate=0.03, N=N, PBCbox=False, precision="mixed",
                                    reporters=[HDF5Reporter(folder=os.path.join(out, "b"), max_data_length=10, overwrite=True)])
        np.random.seed(2)
        box = int(2 * (3 * N / (4 * np.pi * 0.1)) ** (1 / 3) / np.sqrt(3) * 0.98)
        sim.set_data(starting_conformations.grow_cubic(N, box), center=True)
        sim.add_force(forces.spherical_confinement(sim, density=0.1, k=5.0))
        sim.add_force(forcekits.polymer_chains(
            sim, chains=[(0, N, False)], bond_force_func=forces.harmonic_bonds,
            bond_force_kwargs={"bondLength": 1.0, "bondWiggleDistance": 0.1}, angle_force_func=forces.angle_force,
            angle_force_kwargs={"k": 1.5}, nonbonded_force_func=forces.polynomial_repulsive,
            nonbonded_force_kwargs={"trunc": 1.5, "radiusMult": 1.05}, except_bonds=True))
        n_blocks = 40
        traj = np.zeros((n_blocks, 1, 2), np.int32)
        traj[:20, 0] = (100, 600)       # one LEF holding (100, 600) ...
        traj[20:, 0] = (300, 400)       # ... that then moves to (300, 400)
        upd = LEFBondUpdater(sim, traj, length=0.5, wiggle=0.2)
        sim.add_force(upd.force)
        sim.local_energy_minimization()

        def dist(a, b):
            x = sim.get_data()
            return float(np.linalg.norm(x[a] - x[b]))

        for blk in range(n_blocks):
            upd.step(blk)
            sim.do_block(750, save=False)
            if blk == 19:
                first_pair = dist(100, 600)
        assert first_pair < 1.5, f"held pair should be ~0.5 apart, got {first_pair:.2f}"
        assert dist(300, 400) < 1.5, "the bond did not move to the new pair (OpenMM ignores re-pointed bonds)"
        assert dist(100, 600) > 4.0, "the old pair was not released"
    finally:
        shutil.rmtree(out, ignore_errors=True)
