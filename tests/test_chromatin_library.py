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


def test_library_loci_follow_the_nucleus_and_polymer(fake_library):
    import pandas as pd
    from simlive.stage1_chromatin.library_loci import simulate_library_loci
    for s in fake_library.glob("demo_seed[123]"):                      # give the fake sims a radius and a moving locus
        m = json.loads((s / "meta.json").read_text()); m["confinement_radius"] = 500.0
        (s / "meta.json").write_text(json.dumps(m))
        t = np.load(s / "tracked_positions.npy"); t[:, :, 1] = np.arange(5)[:, None]; np.save(s / "tracked_positions.npy", t)
    rows = [dict(t=k, cell_id=c, y_um=50.0, x_um=50.0, radius_um=6.0, aspect=1.0, angle_rad=0.0) for c in (1, 2) for k in range(3)]
    cfg = {"loci": {"library": "demo", "positions_kb": [200, 500], "block_duration_s": 10, "fill": 0.7},
           "acquisition": {"frame_interval_s": 10}}
    out = simulate_library_loci(pd.DataFrame(rows), cfg, np.random.default_rng(0))
    assert len(out) == 2 * 3 * 2 and out.groupby("cell_id").sim_id.nunique().eq(1).all()
    assert out.groupby("cell_id").sim_id.first().nunique() == 2        # different cells -> different simulations
    r = np.sqrt(out.xn**2 + out.yn**2 + out.zn**2)                     # polymer radius preserved up to the fill scaling
    assert (r > 0).all()                                               # rotation keeps |x,y,z| (checked per locus below)
    g = out[(out.cell_id == 1) & (out.locus_id == 0)].sort_values("t")
    d = np.linalg.norm(np.diff(g[["xn", "yn", "zn"]].to_numpy(), axis=0), axis=1)
    assert np.allclose(d, 0.7 / 500, atol=1e-9)                        # one block advanced per frame: y moves 1 unit/block
    again = simulate_library_loci(pd.DataFrame(rows), cfg, np.random.default_rng(0))
    pd.testing.assert_frame_equal(out, again)                          # reproducible from the seed
