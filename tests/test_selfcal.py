"""Corrections from the movie alone (stage6_analysis/selfcal.py, the row correction) and the error-effects scoring (stage7_validation)."""
import numpy as np
import pandas as pd

from simlive.stage6_analysis.localize import row_offsets_removed
from simlive.stage6_analysis.selfcal import background_map, chromatic_shift_from_positions, hot_pixel_map


def _movie(rng, T=60, n=96, hot=12):
    """Integer camera frames: offset pattern + noise + a few hot pixels + 'nuclei' (bright discs) that move across the field."""
    pattern = 100 + rng.normal(0, 1.0, (n, n))
    hot_mask = np.zeros((n, n), bool)
    hot_mask.flat[rng.choice(n * n, hot, replace=False)] = True
    yy, xx = np.mgrid[:n, :n]
    mov = np.empty((T, n, n), np.float32)
    lab = np.zeros((T, n, n), np.uint16)
    for t in range(T):
        f = pattern + rng.normal(0, 2.0, (n, n)) + 40 * hot_mask
        for k, (cy, cx) in enumerate(((20 + t, 30), (70, 10 + t))):       # two nuclei sweeping across
            disc = (yy - cy) ** 2 + (xx - cx) ** 2 < 10 ** 2
            f = f + 300 * disc
            lab[t][disc] = k + 1
        mov[t] = np.rint(f)
    return mov, lab, pattern, hot_mask


def test_hot_pixels_and_background_from_the_movie():
    rng = np.random.default_rng(0)
    mov, lab, pattern, hot = _movie(rng)
    bg, n_free, se = background_map(mov, lab)
    found = hot_pixel_map(bg, se)
    testable = n_free >= 10
    assert (found & ~hot).sum() == 0                                 # no false hot pixels despite integer camera values
    assert (found & hot).sum() >= 0.8 * (hot & testable).sum()        # finds those that are ever free of nuclei
    free = testable & ~hot
    assert np.abs((bg - pattern)[free]).mean() < 0.5                  # the offset pattern itself


def test_colour_shift_from_random_orientations():
    rng = np.random.default_rng(1)
    rows = []
    shift = np.array([0.6, -0.4])
    for c in range(200):                                              # many cells: the true vectors average out
        for t in range(5):
            p = rng.uniform(50, 100, 2)
            v = rng.normal(0, 3.0, 2)
            rows += [(c, t, 0, *p), (c, t, 1, *(p + v + shift))]
    pos = pd.DataFrame(rows, columns=["cell_id", "t", "locus_id", "lab_y_px", "lab_x_px"])
    r = chromatic_shift_from_positions(pos)
    assert r["significant"] and np.allclose(r["shift_px"], shift, atol=0.25)
    few = chromatic_shift_from_positions(pos[pos.cell_id < 4])        # 4 cells: not enough, must not be applied
    assert not few["significant"]


def test_row_offsets_are_removed_using_pixels_away_from_the_cell():
    rng = np.random.default_rng(2)
    S = 60
    row = rng.normal(0, 3, (1, S, 1))
    img = 100 + rng.normal(0, 0.5, (1, S, S)) + row
    mask = np.zeros((1, S, S), bool)
    mask[0, 20:40, 20:40] = True
    img[0][mask[0]] += 200
    out = row_offsets_removed(img, mask)
    resid = (out[0] - 100)[~mask[0]]
    assert resid.std() < 0.8 < (img[0] - 100)[~mask[0]].std()


def test_restricted_calibration_keeps_only_the_chosen_parts():
    from simlive.stage3_microscopy.calibration.model import Calibration

    rng = np.random.default_rng(3)
    maps = {"locus0": {"offset": 100 + rng.normal(0, 1, (8, 8)).astype(np.float32), "flat": 1 + rng.normal(0, 0.1, (8, 8)).astype(np.float32),
                       "hot": rng.random((8, 8)) < 0.1}}
    cal = Calibration(["locus0"], (8, 8), 0.1, 0.13, "locus0", maps, {}, {"shift_px": {"locus0": [0.5, 0.2]}}, [])
    r = cal.restricted(("hot",))
    assert (r.maps["locus0"]["flat"] == 1).all() and np.ptp(r.maps["locus0"]["offset"]) == 0 and r.shift_px("locus0") == (0.0, 0.0)
    assert (r.maps["locus0"]["hot"] == maps["locus0"]["hot"]).all()


def test_error_effects_score_is_perfect_for_perfect_detections():
    from simlive.stage7_validation.error_effects import score

    rng = np.random.default_rng(4)
    rows, det = [], []
    for c in range(5):
        for t in range(40):
            for k in (0, 1):
                y, x, ph = rng.uniform(5, 50), rng.uniform(5, 50), rng.uniform(100, 600)
                rows.append((t, c, k, y, x, ph))
                det.append((c, t, k, y / 0.13, x / 0.13, 0.01 * ph))
    truth = pd.DataFrame(rows, columns=["t", "cell_id", "locus_id", "y_um", "x_um", "photons"])
    pos = pd.DataFrame(det, columns=["cell_id", "t", "locus_id", "lab_y_px", "lab_x_px", "intensity"])
    s = score(pos, truth, 0.13)
    assert s["found"] == 1 and s["frac_far"] == 0 and s["median_err_nm"] < 1e-6 and abs(s["dist_bias_nm"]) < 1e-6
    assert s["intensity_r"] > 0.999 and s["intensity_err_pct"] < 1e-6 and abs(s["late_early"] - 1) < 1e-9


def test_crosstalk_fraction_is_measured_inside_the_nuclei():
    from simlive.stage6_analysis.selfcal import crosstalk_fraction

    rng = np.random.default_rng(5)
    T, n = 20, 120
    yy, xx = np.mgrid[:n, :n]
    nuc = np.zeros((T, n, n), np.float32)
    loc = np.zeros((T, n, n), np.float32)
    lab = np.zeros((T, n, n), np.uint16)
    for t in range(T):
        disc = (yy - 60) ** 2 + (xx - 40 - t) ** 2 < 25 ** 2
        bright = 100 + 300 * disc * (1 + 0.3 * np.sin(xx / 3.0))            # textured nucleus
        spot = 800 * np.exp(-((yy - 55) ** 2 + (xx - 45 - t) ** 2) / 2.0)    # the locus spot (an outlier for the regression)
        nuc[t] = rng.poisson(bright)
        loc[t] = rng.poisson(50 + 0.02 * (bright - 100) + spot)
        lab[t][disc] = 1
    ct = crosstalk_fraction(loc, nuc, lab)
    assert abs(ct["fraction"] - 0.02) < 0.004 and ct["fraction"] > 2 * ct["se"]
    assert abs(crosstalk_fraction(rng.poisson(50 + 0 * nuc).astype(np.float32), nuc, lab)["fraction"]) < 0.004    # no leak -> ~0
