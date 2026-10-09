"""Transcription model: simulation vs exact theory, frame-interval independence, coupling behaviour (no GPU needed)."""
import numpy as np
import pandas as pd
import pytest

from simlive.stage1_chromatin import transcription as T

FAST = dict(k_on=0.02, k_off=0.05, k_init=0.2, coupling="none", f_const=0.0, elongation_kb_min=6.0, cassette_kb=0.6, gene_kb=2.4,
            n_loops=24, dwell_s=10.0)          # transit 30 s, cassette ramp 6 s


def test_constant_distance_matches_theory():
    p = T.TranscriptionParams(**FAST)
    th = T.summary(p, 100.0)
    c = T.simulate_cell(p, np.full(2000, 100.0), 100.0, np.arange(0, 200000.0, 10.0), 0.1, seed=1)
    f = c.frames
    assert abs(f.promoter_on.mean() - th["on_fraction"]) < 0.02
    assert abs(f.n_polII.mean() - th["mean_polII"]) < 0.06 * th["mean_polII"]
    assert abs(f.ms2_loops.mean() - th["mean_ms2_loops"]) < 0.06 * th["mean_ms2_loops"]


def test_two_valued_distance_matches_exact_expectation():
    """Contact coupling with a distance that alternates: mean ON probability and MS2 over many cells equal the exact expectation."""
    p = T.TranscriptionParams(**{**FAST, "coupling": "contact", "fold": 15.0, "k_on": 0.004})
    d = np.where((np.arange(400) // 20) % 2 == 0, 80.0, 400.0)        # 20 steps of 10 s in contact, then 20 apart
    ft = np.arange(100, 4000.0, 25.0)
    runs = [T.simulate_cell(p, d, 10.0, ft, 0.1, seed=s).frames for s in range(300)]
    on = np.mean([r.promoter_on for r in runs], axis=0)
    ms2 = np.mean([r.ms2_loops for r in runs], axis=0)
    p_start, _ = T.expected_on(p, np.tile(d, 10), 10.0)
    exp_on = p_start[np.searchsorted(np.arange(len(p_start)) * 10.0, ft, side="right") - 1]
    assert np.abs(on - exp_on).max() < 0.12                           # binomial noise of 300 cells ~ 0.03 sd
    assert abs(on.mean() - exp_on.mean()) < 0.02
    # expected MS2 for the start-at-stationary distance series (the simulation includes a warm-up, so compare after the transit time)
    e = T.expected_ms2(p, np.tile(d, 10), 10.0, ft)
    ok = ft > 200
    assert abs(ms2[ok].mean() - e[ok].mean()) < 0.1 * e[ok].mean()


def test_changing_the_frame_interval_does_not_change_the_biology():
    p = T.TranscriptionParams(**FAST)
    d = np.full(3000, 100.0)
    means = []
    for dt in (5.0, 20.0, 60.0):
        ft = np.arange(0, 120000.0, dt)
        means.append(T.simulate_cell(p, d, 100.0, ft, 0.1, seed=7).frames.ms2_loops.mean())
    th = T.summary(p, 100.0)["mean_ms2_loops"]
    assert all(abs(m - th) < 0.08 * th for m in means)                # same physics whatever the sampling


def test_no_coupling_is_uncorrelated_and_contact_coupling_is_correlated():
    d = T.standin_distance(6000, 10.0, np.random.default_rng(3), contact_fraction_=0.2, loop_time_s=600.0)
    ft = np.arange(0, 60000.0, 10.0)
    out = {}
    for kind in ("none", "contact"):
        p = T.TranscriptionParams(**{**FAST, "coupling": kind, "f_const": 0.2, "fold": 20.0, "k_on": 0.002,
                                     "d_contact_nm": 150.0})
        f = T.simulate_cell(p, d, 10.0, ft, 0.1, seed=5).frames
        out[kind] = np.corrcoef(f.ms2_loops, (f.d_nm <= 150).astype(float))[0, 1]
    assert abs(out["none"]) < 0.12 and out["contact"] > 0.3


def test_coupled_rate_options_and_shapes():
    p = T.TranscriptionParams(coupling="hill", coupled_rate="k_off", fold=4.0, d_half_nm=200.0, hill_n=4.0)
    kon, koff, kinit = T.rates_from_distance(p, [1.0, 200.0, 1e5])
    assert koff[0] == pytest.approx(p.k_off / 4.0, rel=1e-3) and koff[2] == pytest.approx(p.k_off, rel=1e-6)   # sticky when close
    assert kon[0] == p.k_on and kinit[0] == p.k_init
    assert T.contact_fraction(p, [200.0])[0] == pytest.approx(0.5)
    assert T.contact_fraction(T.TranscriptionParams(coupling="exp", d_decay_nm=100.0), [100.0])[0] == pytest.approx(np.exp(-1))
    with pytest.raises(ValueError):
        T.TranscriptionParams(coupling="bogus")


def test_short_distance_series_is_repeated_and_flagged():
    a, looped = T.pingpong(np.array([1.0, 2.0, 3.0]), 8)
    assert looped and len(a) == 8 and list(a[:5]) == [1, 2, 3, 2, 1]        # no jump at the seam
    assert T.pingpong(np.arange(5.0), 3)[1] is False
    c = T.simulate_cell(T.TranscriptionParams(**FAST), np.full(10, 100.0), 10.0, np.arange(0, 1000.0, 10.0), 0.1, seed=1)
    assert c.looped


def test_pipeline_table_is_reproducible():
    cells = pd.DataFrame({"t": np.tile(np.arange(30), 2), "cell_id": np.repeat([1, 2], 30)})
    cfg = {"acquisition": {"frame_interval_s": 30, "exposure_s": 0.1}, "loci": {"transcription": {"params": FAST}}}
    traces = {c: {"dt_block_s": 10.0, "d_nm": np.full(100, 100.0)} for c in (1, 2)}
    a, ea = T.simulate_transcription(cells, traces, cfg, np.random.default_rng(1))
    b, eb = T.simulate_transcription(cells, traces, cfg, np.random.default_rng(1))
    pd.testing.assert_frame_equal(a, b); pd.testing.assert_frame_equal(ea, eb)
    assert len(a) == 60 and {"ms2_loops", "n_polII", "promoter_on", "d_nm"} <= set(a.columns)


def test_sampling_warnings():
    p = T.TranscriptionParams()
    assert T.check_sampling(p, 10.0, 0.1)[0][0] == "good"
    assert T.check_sampling(p, 600.0, 0.1)[0][0] in ("bad", "warn")
    assert any(l == "bad" for l, _ in T.check_sampling(p, 10.0, 20.0))


# ----------------------------------------------------------------------------- integration with the polymer library and the renderer
@pytest.fixture()
def fake_lib(tmp_path, monkeypatch):
    import json
    import simlive.stage1_chromatin.library as lib
    monkeypatch.setattr(lib, "DATA_ROOT", tmp_path)
    s = tmp_path / "chromatin" / "demo" / "demo_seed1"
    s.mkdir(parents=True)
    (s / "meta.json").write_text(json.dumps({"seed": 1, "n_monomers": 1000, "md_steps_per_second": 3000.0, "wall_seconds_total": 60.0,
                                             "confinement_radius": 20.0}))
    beads = np.arange(0, 1000, 10)
    traj = np.zeros((50, len(beads), 3), np.float32)
    traj[:, :, 0] = np.arange(50)[:, None]                         # every monomer moves +1 unit in x per block
    traj[:, 50, 1] = 1.0                                           # monomer 500 sits 1 unit away in y: distance(500, 510) = sqrt(1+...) known
    np.save(s / "tracked_positions.npy", traj); np.save(s / "tracked_beads.npy", beads)
    return tmp_path


def _cells(n_frames, radius=6.0):
    return pd.DataFrame({"t": np.arange(n_frames), "cell_id": 1, "y_um": 50.0, "x_um": 50.0, "radius_um": radius, "aspect": 1.0,
                         "angle_rad": 0.0})


def test_library_loci_interpolate_between_blocks_and_return_distance_traces(fake_lib):
    from simlive.stage1_chromatin.library_loci import simulate_library_loci
    cfg = {"loci": {"library": "demo", "positions_kb": [500, 510], "nm_per_unit": 100.0, "block_duration_s": 10.0},
           "acquisition": {"frame_interval_s": 4.0}}                # 0.4 block per frame: needs interpolation
    df, tr = simulate_library_loci(_cells(20), cfg, np.random.default_rng(0), return_traces=True)
    g = df[df.locus_id == 0].sort_values("t")
    steps_nm = np.linalg.norm(np.diff(g[["xn", "yn", "zn"]].to_numpy(), axis=0), axis=1) * 6000.0
    assert np.allclose(steps_nm, 0.4 * 100.0, atol=1e-6)            # smooth motion: 0.4 block x 100 nm per frame
    d = tr[1]["d_nm"]
    assert tr[1]["dt_block_s"] == 10.0 and np.allclose(d, 100.0)    # monomers 500 and 510 are 1 unit apart (fixture) = 100 nm
    # a different frame interval reads the same polymer faster/slower; the distance trace does not depend on it
    cfg2 = {**cfg, "acquisition": {"frame_interval_s": 30.0}}
    _, tr2 = simulate_library_loci(_cells(15), cfg2, np.random.default_rng(0), return_traces=True)
    assert np.allclose(tr2[1]["d_nm"][:5], 100.0)


def test_long_movie_repeats_the_trajectory_with_a_warning(fake_lib):
    from simlive.stage1_chromatin.library_loci import simulate_library_loci
    cfg = {"loci": {"library": "demo", "positions_kb": [500, 510], "nm_per_unit": 10.0, "block_duration_s": 10.0},
           "acquisition": {"frame_interval_s": 60.0}}               # 40 frames x 6 blocks = 240 blocks > 50 saved
    with pytest.warns(UserWarning, match="repeated back and forth"):
        df, tr = simulate_library_loci(_cells(40, radius=60.0), cfg, np.random.default_rng(1), return_traces=True)
    assert len(df) == 80 and tr[1]["looped"] is True


def test_ms2_spot_appears_in_the_nuclear_channel_and_scales_with_loops():
    from simlive.io.runs import REPO_ROOT, load_config
    from simlive.stage3_microscopy.render import render_nuclei
    cfg = load_config(REPO_ROOT / "configs" / "loci_demo.yaml", ["geometry.fov_um=[20, 20]", "acquisition.n_frames=2",
                                                                "cells.n_cells=1", "optics.read_noise_e=0", "optics.bleach_tau_s=null",
                                                                "optics.photons_per_px_s=0", "optics.background_photons=0",
                                                                "imaging_errors.enabled=false"])
    cells = pd.DataFrame({"t": [0, 1], "cell_id": 1, "y_um": 10.0, "x_um": 10.0, "radius_um": 5.0, "aspect": 1.0, "angle_rad": 0.0,
                          "bound_radius_um": 5.5, "parent_id": 0, "amp2": 0.0, "phase2": 0.0, "amp3": 0.0, "phase3": 0.0,
                          "amp4": 0.0, "phase4": 0.0, "amp5": 0.0, "phase5": 0.0})
    spots = pd.DataFrame({"t": [0, 1], "y_um": 10.0, "x_um": 10.0, "photons": [1000.0, 3000.0]})
    img, _ = render_nuclei(cells, cfg, np.random.default_rng(0), spots=spots)
    off = cfg["optics"]["offset_adu"]
    sig = [(img[t] - off).sum() for t in (0, 1)]
    assert abs(sig[1] / sig[0] - 3.0) < 0.2 and abs(sig[0] - 1000 * cfg["optics"]["quantum_efficiency"]) < 0.1 * 1000
    peak = np.unravel_index(np.argmax(img[1]), img[1].shape)
    px = cfg["optics"]["pixel_size_nm"] / 1000
    assert abs(peak[0] * px - 10.0) < 0.2 and abs(peak[1] * px - 10.0) < 0.2        # the spot sits at the requested position


def test_ms2_spots_table_links_promoter_position_and_signal():
    truth = pd.DataFrame({"t": [0, 0], "cell_id": 1, "locus_id": [0, 1], "y_um": [5.0, 6.0], "x_um": [7.0, 8.0]})
    tx = pd.DataFrame({"t": [0], "cell_id": 1, "ms2_loops": [48.0]})
    cfg = {"loci": {"transcription": {"photons_per_loop_s": 100.0}}, "acquisition": {"exposure_s": 0.1}}
    s = T.ms2_spots(tx, truth, cfg)
    assert s.iloc[0].y_um == 5.0 and s.iloc[0].x_um == 7.0 and s.iloc[0].photons == pytest.approx(480.0)
