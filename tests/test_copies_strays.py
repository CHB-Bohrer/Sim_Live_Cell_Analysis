"""Several promoter-enhancer copies per cell (stage1_chromatin/copies.py), stray spots (image error source `stray_spots`) and the multi-copy
linker that must tell them apart (stage6_analysis/multicopy.py), plus the copy-aware scores."""
import numpy as np
import pandas as pd
import pytest

from simlive.io.runs import REPO_ROOT, load_config
from simlive.stage1_chromatin import transcription as T
from simlive.stage1_chromatin.copies import copy_settings, split_trace_id, trace_id
from simlive.stage1_chromatin.toy_loci import simulate_toy_loci
from simlive.stage2_cells.motion import simulate_cells
from simlive.stage3_microscopy.errors import ErrorSet
from simlive.stage3_microscopy.render import Imager
from simlive.stage6_analysis import multicopy as MC
from simlive.stage6_analysis.localize import EXTRA_COPY, locus_distances, match_copies, score_localization

CFG = REPO_ROOT / "configs" / "loci_demo.yaml"
SMALL = ["geometry.fov_um=[40, 40]", "cells.n_cells=3", "acquisition.n_frames=30", "motion.speed_scale=0.3"]


def _cfg(*extra):
    return load_config(CFG, SMALL + list(extra))


def test_new_random_streams_leave_the_old_ones_unchanged():
    old, new = np.random.SeedSequence(7).spawn(8), np.random.SeedSequence(7).spawn(10)
    for a, b in zip(old, new):
        assert np.array_equal(np.random.default_rng(a).integers(0, 2**31, 5), np.random.default_rng(b).integers(0, 2**31, 5))


def test_one_copy_and_copy_zero_reproduce_the_old_toy_loci():
    cfg1 = _cfg("loci.copies=null")
    cells = simulate_cells(cfg1, np.random.default_rng(1))
    one = simulate_toy_loci(cells, cfg1, np.random.default_rng(2), np.random.default_rng(3))
    assert set(one["copy"]) == {0}
    two = simulate_toy_loci(cells, _cfg("loci.copies.n=[2]", "loci.copies.weights=[1]"), np.random.default_rng(2), np.random.default_rng(3))
    a = two[two["copy"] == 0].reset_index(drop=True)
    pd.testing.assert_frame_equal(a, one.reset_index(drop=True))          # copy 0 = exactly the one-copy run
    assert (two.groupby("cell_id")["copy"].nunique() == 2).all()


def test_extra_copies_sit_apart_and_inside_their_nucleus():
    cfg = _cfg("loci.copies.n=[2]", "loci.copies.weights=[1]", "loci.copies.min_separation_um=2.5", "imaging_errors.enabled=false",
               "geometry.fov_um=[60, 60]")
    cells = simulate_cells(cfg, np.random.default_rng(4))
    truth = simulate_toy_loci(cells, cfg, np.random.default_rng(5), np.random.default_rng(6))
    centres = truth.groupby(["cell_id", "copy"])[["xn", "yn"]].mean()                    # nucleus frame (radii)
    r = cells.groupby("cell_id").radius_um.min()
    for cid in truth.cell_id.unique():
        d_um = np.hypot(*(centres.loc[(cid, 0)] - centres.loc[(cid, 1)]).to_numpy()) * r[cid]
        assert d_um > 1.8, f"cell {cid}: copies only {d_um:.2f} um apart"
    im = Imager(cfg, cells, truth)
    lab = im.frame(0)["labels"]
    px = im.px_um
    for q in truth[truth.t == 0].itertuples():
        assert lab[int(round(q.y_um / px)), int(round(q.x_um / px))] == q.cell_id


def test_every_copy_bursts_on_its_own_with_its_own_parameters():
    cfg = _cfg("loci.transcription.enabled=true", "loci.copies.transcription=[{}, {k_on: 0.01}]")
    cells = pd.DataFrame({"t": np.arange(2000), "cell_id": 1})
    rng_d = np.random.default_rng(0)
    d = rng_d.uniform(400, 900, 4000)                                           # far: no enhancer contact, plain k_on
    traces = {1: {"dt_block_s": 10.0, "d_nm": d}, (1, 1): {"dt_block_s": 10.0, "d_nm": d}}
    tx, ev = T.simulate_transcription(cells, traces, cfg, np.random.default_rng(1), np.random.default_rng(2))
    one, _ = T.simulate_transcription(cells, {1: traces[1]}, cfg, np.random.default_rng(1), np.random.default_rng(2))
    c0, c1 = tx[tx["copy"] == 0], tx[tx["copy"] == 1]
    assert np.array_equal(c0.promoter_on.to_numpy(), one.promoter_on.to_numpy())          # copy 0 is the one-copy run
    p = T.params_from_config(cfg, 1)
    assert p.k_on == 0.01 and T.params_from_config(cfg, 0).k_on == T.params_from_config(cfg).k_on
    on0, on1 = c0.promoter_on.mean(), c1.promoter_on.mean()
    exp1 = 0.01 / (0.01 + p.k_off)
    assert on1 > 3 * on0 and abs(on1 - exp1) < 0.5 * exp1                       # copy 1 follows ITS k_on
    assert np.corrcoef(c0.promoter_on, c1.promoter_on)[0, 1] < 0.1              # independent promoters
    assert set(tx.trace_id) == {trace_id(1, 0), trace_id(1, 1)} and split_trace_id(trace_id(1, 1)) == (1, 1)


def test_ms2_cassette_on_the_first_copy_only():
    cfg = _cfg("loci.transcription.enabled=true", "loci.copies.ms2=first")
    tx = pd.DataFrame({"t": [0, 0], "cell_id": [1, 1], "copy": [0, 1], "ms2_loops": [24.0, 24.0],
                       "ms2_labelled": [copy_settings(cfg).labelled(0), copy_settings(cfg).labelled(1)]})
    truth = pd.DataFrame({"t": [0, 0, 0, 0], "cell_id": 1, "locus_id": [0, 1, 0, 1], "copy": [0, 0, 1, 1], "y_um": 1.0, "x_um": 1.0})
    s = T.ms2_spots(tx, truth, cfg)
    assert s.set_index("copy").photons[0] > 0 and s.set_index("copy").photons[1] == 0


def test_stray_spots_follow_their_rate_lifetime_and_stay_in_the_nucleus():
    cfg = _cfg("acquisition.n_frames=120", "imaging_errors.stray_spots.rate_per_min=3.0", "imaging_errors.stray_spots.lifetime_s=40")
    cells = simulate_cells(cfg, np.random.default_rng(1))
    im = Imager(cfg, cells, None)
    st = im.stray_spots()
    assert len(st) and set(st.channel) <= {"nucleus"}          # no locus truth here: only the nuclear colour has a reference
    life = st.groupby(["cell_id", "stray_id"]).t.agg(["min", "max", "size"])
    dt = cfg["acquisition"]["frame_interval_s"]
    inner = life[(life["min"] > 0) & (life["max"] < cfg["acquisition"]["n_frames"] - 1)]
    assert abs(inner["size"].mean() * dt - 40 - dt / 2) < 15                     # ceil(exponential) frames
    births = life[life["min"] > 0].groupby(level=0).size().mean() / ((cfg["acquisition"]["n_frames"] - 1) * dt / 60)
    assert 1.5 < births < 4.5
    for q in st[st.t == 5].itertuples():
        assert im.frame(5)["labels"][int(round(q.y_um / im.px_um)), int(round(q.x_um / im.px_um))] == q.cell_id
    pd.testing.assert_frame_equal(st, Imager(cfg, cells, None).stray_spots())    # deterministic


def test_switching_stray_spots_on_adds_only_the_spots():
    cfg = _cfg("loci.transcription.enabled=false")
    cells = simulate_cells(cfg, np.random.default_rng(1))
    truth = simulate_toy_loci(cells, cfg, np.random.default_rng(2))
    im = Imager(cfg, cells, truth)
    E_on = im.E
    E_off = E_on.with_enabled(set(E_on.enabled_ids()) - {"stray_spots"})
    st = im.stray_spots()
    t = int(st[st.channel == "locus0"].t.iloc[0])
    a, b = im.signal("locus0", t, E_on)[0], im.signal("locus0", t, E_off)[0]
    diff = a - b
    q = st[(st.channel == "locus0") & (st.t == t)].iloc[0]
    dy, dx = im.drift_px(t)
    cy, cx = im.channel_shift_px("locus0")
    y, x = q.y_um / im.px_um + dy + cy, q.x_um / im.px_um + dx + cx
    assert diff.max() > 0 and diff.min() >= -1e-3
    iy, ix = np.unravel_index(np.argmax(diff), diff.shape)
    assert np.hypot(iy - y, ix - x) < 1.5 or len(st[(st.channel == "locus0") & (st.t == t)]) > 1


# ---------------------------------------------------------------- the linker on synthetic candidates
def _scene(rng, n=80, copies=((40, 40, 46, 43), (85, 80, 80, 86)), n_stray=12, noise=2):
    """Candidates of two colours per frame: each copy's promoter and enhancer (slow walks, SNR ~12), stray spots (live 6 frames, SNR ~10)
    in either colour anywhere in a 120 px cell, and weak noise peaks. Returns peaks, snr, truth index per (colour, copy) and frame."""
    walk = lambda: np.cumsum(rng.normal(0, 0.6, (n, 2)), axis=0)       # noqa: E731
    obj = {}
    for c, (py, px, ey, ex) in enumerate(copies):
        obj[(0, c)] = np.array([py, px]) + walk()
        obj[(1, c)] = np.array([ey, ex]) + walk()
    strays = []
    for _ in range(n_stray):
        t0 = int(rng.integers(0, n - 6))
        strays.append((int(rng.integers(0, 2)), t0, rng.uniform(10, 110, 2)))
    peaks = [[[] for _ in range(n)] for _ in range(2)]
    snr = [[[] for _ in range(n)] for _ in range(2)]
    true_idx = {}
    for i in range(n):
        for k in range(2):
            for c in range(len(copies)):
                if rng.random() < 0.95:                     # a few frames without the locus
                    true_idx[(k, c, i)] = len(peaks[k][i])
                    peaks[k][i].append(obj[(k, c)][i]); snr[k][i].append(rng.uniform(8, 16))
            for kk, t0, p in strays:
                if kk == k and t0 <= i < t0 + 6:
                    peaks[k][i].append(p + rng.normal(0, 0.5, 2)); snr[k][i].append(rng.uniform(8, 14))
            for _ in range(noise):
                peaks[k][i].append(rng.uniform(5, 115, 2)); snr[k][i].append(rng.uniform(3, 4.5))
    peaks = [[np.array(f, float).reshape(-1, 2) for f in pk] for pk in peaks]
    snr = [[np.array(f, float) for f in s] for s in snr]
    return peaks, snr, true_idx, n


def _agreement(res_copies, true_idx, n, k, c_found, c_true):
    hits = [res_copies[c_found][k][i] == true_idx[(k, c_true, i)] for i in range(n) if (k, c_true, i) in true_idx]
    return float(np.mean(hits))


def test_linker_finds_both_copies_pairs_them_and_ignores_stray_spots():
    rng = np.random.default_rng(11)
    peaks, snr, true_idx, n = _scene(rng)
    centre = np.full((n, 2), 60.0)
    res = MC.track_copies(peaks, snr, centre, np.ones(n, bool), px_um=0.13, max_copies=2, pair_gate_um=1.5, search_volume=120.0 ** 2)
    assert len(res["copies"]) == 2 and not res["fallback"]
    first = [np.nanmedian([peaks[0][i][j][0] if j is not None else np.nan for i, j in enumerate(cp[0])]) for cp in res["copies"]]
    order = np.argsort(first)                                  # found copy whose promoter sits near y = 40 is true copy 0
    for c_true, c_found in enumerate(order):
        for k in (0, 1):
            assert _agreement(res["copies"], true_idx, n, k, c_found, c_true) > 0.93


def test_one_copy_cell_gives_one_copy_although_two_are_allowed():
    rng = np.random.default_rng(12)
    peaks, snr, true_idx, n = _scene(rng, copies=((50, 50, 55, 53),))
    res = MC.track_copies(peaks, snr, np.full((n, 2), 60.0), np.ones(n, bool), px_um=0.13, max_copies=2, pair_gate_um=1.5,
                          search_volume=120.0 ** 2)
    assert len(res["copies"]) == 1
    for k in (0, 1):
        assert _agreement(res["copies"], true_idx, n, k, 0, 0) > 0.93


def test_pair_gate_keeps_the_path_off_a_stray_after_the_locus_goes_dark():
    """The promoter goes dark after frame 14 (no probe attached, or out of focus) and a stray spot appears later, far away: one path
    alone re-acquires the stray (a spot is worth more than a run of misses); with the pair gate (the enhancer stays visible next to where
    the promoter was) it does not."""
    n = 40
    rng = np.random.default_rng(3)
    prom = np.array([50.0, 50.0]) + np.cumsum(rng.normal(0, 0.3, (n, 2)), axis=0)
    enh = prom + np.array([4.0, 3.0])
    peaks = [[], []]
    snr = [[], []]
    for i in range(n):
        p, s = [], []
        if i < 15:
            p.append(prom[i]); s.append(10.0)
        if 20 <= i < 35:
            p.append(np.array([100.0, 20.0])); s.append(30.0)     # stray, far from the enhancer, brighter
        peaks[0].append(np.array(p, float).reshape(-1, 2)); snr[0].append(np.array(s, float))
        peaks[1].append(enh[i][None]); snr[1].append(np.array([10.0]))
    kw = dict(px_um=0.13, max_copies=1, search_volume=120.0 ** 2)
    plain = MC.track_copies(peaks, snr, np.full((n, 2), 60.0), np.ones(n, bool), pair_gate_um=None, **kw)
    gated = MC.track_copies(peaks, snr, np.full((n, 2), 60.0), np.ones(n, bool), pair_gate_um=1.5, **kw)
    on_stray = lambda ch: sum(1 for i, j in enumerate(ch) if j is not None and peaks[0][i][j][0] > 90)    # noqa: E731
    assert on_stray(plain["copies"][0][0]) >= 1
    assert on_stray(gated["copies"][0][0]) == 0


def test_untangle_undoes_two_paths_that_swapped_copies():
    n = 30
    a = np.c_[np.full(n, 20.0), np.linspace(20, 25, n)]          # copy A drifts slowly at y = 20
    b = np.c_[np.full(n, 34.0), np.linspace(40, 36, n)]          # copy B at y = 34, 14 px away
    peaks = [[np.vstack([a[i], b[i]]) for i in range(n)]]
    # greedy extraction result: path 0 = A until frame 17, then B; path 1 = the complement
    p0 = [0 if i < 18 else 1 for i in range(n)]
    p1 = [1 if i < 18 else 0 for i in range(n)]
    fixed = MC.untangle([{0: p0}, {0: p1}], peaks, np.zeros((n, 2)))
    assert fixed[0][0] == [0] * n and fixed[1][0] == [1] * n
    # a copy = promoter + enhancer swapped together is untangled by their summed steps
    fixed2 = MC.untangle([{0: p0, 1: p0}, {0: p1, 1: p1}], peaks + peaks, np.zeros((n, 2)))
    assert fixed2[0][0] == [0] * n and fixed2[0][1] == [0] * n
    same = MC.untangle([{0: [0] * n}, {0: [1] * n}], peaks, np.zeros((n, 2)))   # nothing to undo: unchanged
    assert same[0][0] == [0] * n and same[1][0] == [1] * n


def test_pairing_keeps_a_good_pair_instead_of_shortening_an_invalid_one():
    # promoter paths (rows) x enhancer paths (columns): (0, 0) is a real pair; path 1 of each colour wanders far from everything
    C = np.array([[0.6, 15.0], [15.0, 40.0]])
    assert MC.assign_pairs(C, 11.5) == [(0, 0)]          # raw costs would choose (0, 1) + (1, 0) = 30 < 40.6 and lose the real pair
    assert MC.assign_pairs(np.array([[3.0, 1.0], [1.0, 3.0]]), 11.5) == [(0, 1), (1, 0)]


# ---------------------------------------------------------------- copy-aware scores
def _truth_two_copies(n=10):
    rows = []
    for t in range(n):
        for c, base in enumerate(((5.0, 5.0), (9.0, 9.0))):
            for k in range(2):
                rows.append((t, 1, k, c, base[0] + 0.3 * k, base[1], 100.0))
    return pd.DataFrame(rows, columns=["t", "cell_id", "locus_id", "copy", "y_um", "x_um", "photons"])


def test_found_copies_are_matched_to_the_true_ones_whatever_their_numbers():
    truth = _truth_two_copies()
    px = 0.1
    pos = truth.assign(copy=1 - truth["copy"], lab_y_px=truth.y_um / px, lab_x_px=truth.x_um / px)      # numbers swapped
    m = match_copies(pos, truth, px)
    assert (m["copy"] == truth["copy"]).all() and (m.copy_found == 1 - truth["copy"]).all()
    err, sc = score_localization(pos, truth, px)
    assert sc["frac_far"] == 0 and sc["copy_number_right"] == 1.0 and sc["missed_copies"] == 0 and sc["extra_copies"] == 0
    d = locus_distances(m, px)
    assert np.allclose(d.dist_um, 0.3)                                     # promoter-enhancer of the SAME copy


def test_an_extra_copy_and_picks_of_a_stray_are_counted():
    truth = _truth_two_copies()
    px = 0.1
    one = truth[truth["copy"] == 0]
    pos = one.assign(lab_y_px=one.y_um / px, lab_x_px=one.x_um / px)
    stray_pick = pos.iloc[[0]].assign(lab_y_px=20.0 / px, lab_x_px=2.0 / px)        # frame 0 promoter taken at a stray spot
    pos = pd.concat([stray_pick, pos.iloc[1:]])
    extra = pos.assign(copy=1, lab_y_px=1.0 / px, lab_x_px=15.0 / px)               # a third "copy" at nothing
    truth1 = one
    strays = pd.DataFrame({"t": [0], "cell_id": [1], "channel": ["locus0"], "y_um": [20.05], "x_um": [2.0]})
    err, sc = score_localization(pd.concat([pos, extra]), truth1, px, strays=strays)
    assert sc["extra_copies"] == 1 and sc["extra_detections"] == len(extra)
    assert sc["frac_on_stray"] == pytest.approx(1 / len(pos))
    assert (match_copies(pd.concat([pos, extra]), truth1, px)["copy"] >= EXTRA_COPY).sum() == len(extra)
