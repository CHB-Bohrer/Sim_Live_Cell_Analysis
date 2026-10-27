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
        m = json.loads((s / "meta.json").read_text()); m["confinement_radius"] = 20.0
        (s / "meta.json").write_text(json.dumps(m))
        t = np.load(s / "tracked_positions.npy"); t[:, :, 1] = np.arange(5)[:, None]; np.save(s / "tracked_positions.npy", t)
    rows = [dict(t=k, cell_id=c, y_um=50.0, x_um=50.0, radius_um=6.0, aspect=1.0, angle_rad=0.0) for c in (1, 2) for k in range(3)]
    cfg = {"loci": {"library": "demo", "positions_kb": [20, 50], "block_duration_s": 10, "nm_per_unit": 50.0},
           "acquisition": {"frame_interval_s": 10}}
    out = simulate_library_loci(pd.DataFrame(rows), cfg, np.random.default_rng(0))
    assert len(out) == 2 * 3 * 2 and out.groupby("cell_id").sim_id.nunique().eq(1).all()
    assert out.groupby("cell_id").sim_id.first().nunique() == 2        # different cells -> different simulations
    d = out[(out.cell_id == 1) & (out.locus_id == 0)].sort_values("t")
    steps_nm = np.linalg.norm(np.diff(d[["xn", "yn", "zn"]].to_numpy(), axis=0), axis=1) * 6.0 * 1000   # radius 6 um
    assert np.allclose(steps_nm, 50.0, atol=1e-6)                      # one block = 1 polymer unit = 50 nm: TRUE distances
    two = out[(out.cell_id == 1) & (out.t == 0)].sort_values("locus_id")   # loci at 20 and 50 kb: x differs by 30 units
    sep = np.linalg.norm(two[["xn", "yn", "zn"]].to_numpy()[0] - two[["xn", "yn", "zn"]].to_numpy()[1]) * 6000
    assert np.isclose(sep, 30 * 50.0, rtol=1e-6)                     # separation in nm independent of rotation / placement
    assert (np.linalg.norm(out[["xn", "yn", "zn"]].to_numpy(), axis=1) <= 0.7 + 1e-9).all()   # loci always inside the nucleus
    with pytest.raises(ValueError):                                    # loci excursion larger than the whole nucleus
        simulate_library_loci(pd.DataFrame(rows), {**cfg, "loci": {**cfg["loci"], "nm_per_unit": 5000.0}}, np.random.default_rng(0))
    again = simulate_library_loci(pd.DataFrame(rows), cfg, np.random.default_rng(0))
    pd.testing.assert_frame_equal(out, again)                          # reproducible from the seed


@pytest.mark.parametrize("three_d", [False, True])
def test_library_loci_keep_true_distances_in_elongated_flattened_nuclei(fake_library, three_d):
    """The nucleus' shape moves and deforms the territory, not the polymer: in the drawn (lab) coordinates the promoter-enhancer distance
    equals the polymer's own distance in nm, whatever the aspect ratio, orientation and thickness (before 2026-10-26 it was stretched by
    sqrt(aspect) along the long axis)."""
    import pandas as pd
    from simlive.stage1_chromatin.library_loci import simulate_library_loci
    for s in fake_library.glob("demo_seed[123]"):
        m = json.loads((s / "meta.json").read_text()); m["confinement_radius"] = 20.0
        (s / "meta.json").write_text(json.dumps(m))
    rows = [dict(t=k, cell_id=c, y_um=50.0, x_um=50.0, radius_um=6.0, aspect=1.2 + 0.2 * k, angle_rad=0.3 * c + 0.1 * k,
                 **({"z_um": 6.0, "z_semi_um": 3.0} if three_d else {})) for c in (1, 2) for k in range(3)]
    cfg = {"loci": {"library": "demo", "positions_kb": [20, 50], "block_duration_s": 10, "nm_per_unit": 50.0},
           "acquisition": {"frame_interval_s": 10}}
    out = simulate_library_loci(pd.DataFrame(rows), cfg, np.random.default_rng(0))
    for _, g in out.groupby(["cell_id", "t"]):
        g = g.sort_values("locus_id")
        d = np.hypot(*np.diff(g[["y_um", "x_um"]].to_numpy(), axis=0)[0])
        dz = np.diff(g.z_um.to_numpy())[0] if three_d else 6.0 * np.diff(g.zn.to_numpy())[0]   # 2D movies keep z in nuclear radii (6 um)
        assert np.isclose(np.hypot(d, dz) * 1000, 30 * 50.0, rtol=1e-6)    # 30 polymer units apart = 1500 nm
    assert (np.linalg.norm(out[["xn", "yn", "zn"]].to_numpy(), axis=1) <= 0.7 + 1e-9).all()   # still inside the nucleus


def test_pool_restricts_and_pin_forces(fake_library):
    ids = sorted(lib.list_library("demo").sim_id)
    a = lib.assign_simulations("demo", [1, 2], np.random.default_rng(0), pool=ids[:2])
    assert set(a.values()) == set(ids[:2])
    b = lib.assign_simulations("demo", [1, 2, 3], np.random.default_rng(0), pinned={2: ids[2]})
    assert b[2] == ids[2] and len(set(b.values())) == 3                    # pinned sim is not handed to anyone else
    assert b == lib.assign_simulations("demo", [1, 2, 3], np.random.default_rng(0), pinned={2: ids[2]})
    with pytest.raises(ValueError):
        lib.assign_simulations("demo", [1], np.random.default_rng(0), pinned={1: "nope_seed9"})
    with pytest.raises(ValueError):
        lib.assign_simulations("demo", [1, 2, 3], np.random.default_rng(0), pool=ids[:2])   # 3 cells, 2 sims, no replace


def test_random_rotation_is_a_proper_rotation():
    from simlive.stage1_chromatin.rotation import random_rotation as _random_rotation
    rng = np.random.default_rng(3)
    for _ in range(20):
        R = _random_rotation(rng)
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-12) and np.isclose(np.linalg.det(R), 1.0)
    # mean of many random rotations applied to a fixed vector is ~0 (uniform over the sphere)
    v = np.mean([_random_rotation(rng)[:, 0] for _ in range(4000)], axis=0)
    assert np.abs(v).max() < 0.05


def test_draw_pair_is_on_the_grid_inside_the_margins_and_on_both_orientations():
    from simlive.stage1_chromatin.library_loci import draw_pair
    beads = np.arange(0, 10000, 10)
    rng = np.random.default_rng(3)
    pairs = [draw_pair(rng, 10000, beads, 100.0, 1000.0) for _ in range(400)]
    for p, e in pairs:
        assert abs(p - e) == 100 and p % 10 == 0 and e % 10 == 0
        assert min(p, e) >= 1000 and max(p, e) <= 10000 - 1 - 1000
    assert {p < e for p, e in pairs} == {True, False}                     # promoter on either side
    assert len({min(p, e) for p, e in pairs}) > 100                       # a wide spread of positions
    with pytest.raises(ValueError):
        draw_pair(rng, 1500, beads[:150], 100.0, 1000.0)


def test_random_pairs_differ_between_cells_are_reproducible_and_recorded(fake_library):
    import pandas as pd
    from simlive.stage1_chromatin.library_loci import simulate_library_loci
    for s in fake_library.glob("demo_seed[123]"):
        m = json.loads((s / "meta.json").read_text()); m["confinement_radius"] = 20.0
        (s / "meta.json").write_text(json.dumps(m))
    rows = [dict(t=k, cell_id=c, y_um=50.0, x_um=50.0, radius_um=6.0, aspect=1.0, angle_rad=0.0) for c in (1, 2, 3) for k in range(3)]
    cfg = {"loci": {"library": "demo", "positions_kb": "random", "pair_separation_kb": 100, "pair_margin_kb": 200, "block_duration_s": 10, "nm_per_unit": 1.0},
           "acquisition": {"frame_interval_s": 10}}
    a = simulate_library_loci(pd.DataFrame(rows), cfg, np.random.default_rng(5))
    b = simulate_library_loci(pd.DataFrame(rows), cfg, np.random.default_rng(5))
    pairs = a.attrs["pairs_kb"]
    assert pairs == b.attrs["pairs_kb"] and sorted(pairs) == [1, 2, 3]
    assert all(abs(p - e) == 100 and min(p, e) >= 200 and max(p, e) <= 799 for p, e in pairs.values())
    assert len({tuple(v) for v in pairs.values()}) > 1                    # different cells, different pairs
    pd.testing.assert_frame_equal(a, b)
    with pytest.raises(ValueError, match="positions_kb"):
        simulate_library_loci(pd.DataFrame(rows), {**cfg, "loci": {**cfg["loci"], "positions_kb": "middle"}}, np.random.default_rng(0))
