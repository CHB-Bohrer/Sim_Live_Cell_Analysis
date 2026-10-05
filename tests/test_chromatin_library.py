"""Library bookkeeping tests (no GPU): listing, one-simulation-per-cell assignment, locus extraction."""
import json

import numpy as np
import pytest

import simlive.stage1_chromatin.library as lib


@pytest.fixture()
def fake_library(tmp_path, monkeypatch):
    monkeypatch.setattr(lib, "DATA_ROOT", tmp_path)
    d = tmp_path / "chromatin" / "demo"
    for seed in (1, 2, 3):
        s = d / f"demo_seed{seed}"
        s.mkdir(parents=True)
        (s / "meta.json").write_text(json.dumps({"seed": seed, "n_monomers": 1000, "md_steps_per_second": 3000.0,
                                                 "wall_seconds_total": 600.0, "mean_loop_size_kb": 90.0}))
        beads = np.arange(0, 1000, 10)
        traj = np.zeros((5, len(beads), 3), np.float32)
        traj[:, :, 0] = beads[None, :]                    # x coordinate = monomer index, so we can check the snapping
        np.save(s / "tracked_positions.npy", traj)
        np.save(s / "tracked_beads.npy", beads)
    (d / "demo_seed9").mkdir()                            # unfinished simulation (no meta.json) must be ignored
    return d


def test_list_library_ignores_unfinished(fake_library):
    df = lib.list_library("demo")
    assert sorted(df.seed) == [1, 2, 3] and df.extrusion.all()


def test_each_cell_gets_its_own_simulation(fake_library):
    a = lib.assign_simulations("demo", [10, 11, 12], np.random.default_rng(0))
    assert len(set(a.values())) == 3
    assert a == lib.assign_simulations("demo", [10, 11, 12], np.random.default_rng(0))   # reproducible from the seed


def test_more_cells_than_simulations_must_be_explicit(fake_library):
    with pytest.raises(ValueError):
        lib.assign_simulations("demo", range(5), np.random.default_rng(0))
    a = lib.assign_simulations("demo", range(5), np.random.default_rng(0), replace=True)
    assert len(a) == 5


def test_locus_positions_snap_to_nearest_stored_bead(fake_library):
    traj, beads = lib.locus_trajectories(fake_library / "demo_seed1", [203, 500])
    assert list(beads) == [200, 500] and traj.shape == (5, 2, 3)
    assert np.all(traj[:, 0, 0] == 200)
