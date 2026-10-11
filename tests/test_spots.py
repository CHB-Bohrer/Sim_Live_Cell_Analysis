"""Stage 6 spot finder (dog_viterbi, Chen et al. 2025) and the localizer around it: known answers on synthetic images."""
import numpy as np
import pandas as pd
import pytest
from scipy.special import erf

from simlive.stage5_linking.celldata import CellData
from simlive.stage6_analysis import dog_viterbi as DV
from simlive.stage6_analysis import spots as SP
from simlive.stage6_analysis.localize import detect_channel, localize_loci, score_localization

S = 41


def _spot(cy, cx, photons, sigma=0.85, shape=(S, S)):
    """Gaussian PSF integrated over each pixel (as the renderer draws it), centred at array index coordinates (cy, cx)."""
    s2 = np.sqrt(2) * sigma
    ys, xs = np.arange(shape[0]), np.arange(shape[1])
    gy = 0.5 * (erf((ys + 0.5 - cy) / s2) - erf((ys - 0.5 - cy) / s2))
    gx = 0.5 * (erf((xs + 0.5 - cx) / s2) - erf((xs - 0.5 - cx) / s2))
    return photons * np.outer(gy, gx)


def _region(shape=(S, S), r=16):
    yy, xx = np.indices(shape)
    return np.hypot(yy - shape[0] // 2, xx - shape[1] // 2) <= r


def _noisy(img, rng, bg=20.0, read=1.5):
    return rng.poisson(img + bg).astype(float) + rng.normal(0, read, img.shape) + 100.0


def _first(img, reg, **kw):
    """Brightest DoG candidate of one image, refined (spots.refine_peak)."""
    c = DV.dog_candidates(img, reg, **kw)
    return SP.refine_peak(c, 0, **{**SP.DEFAULTS, **kw}) if c is not None and len(c["peaks"]) else None


def test_gaussian_fit_finds_a_noise_free_spot_to_a_hundredth_of_a_pixel():
    img = _spot(20.3, 18.7, 2000) + 100.0
    r = _first(img, _region())
    assert abs(r["pos"][0] - 20.3) < 0.01 and abs(r["pos"][1] - 18.7) < 0.01
    assert _first(img, _region(), refine="none")["pos"] == (20.0, 19.0)      # integer pixel of the maximum


def test_noisy_spots_are_found_with_sub_pixel_accuracy():
    rng = np.random.default_rng(1)
    err = []
    for _ in range(40):
        cy, cx = rng.uniform(16, 24, 2)
        r = _first(_noisy(_spot(cy, cx, 500), rng), _region())
        err.append(np.hypot(r["pos"][0] - cy, r["pos"][1] - cx))
    assert np.sqrt(np.mean(np.square(err))) < 0.15        # 500 photons, sigma 0.85 px: Thompson limit ~0.05 px; generous bound


def test_median_filter_removes_a_hot_pixel_that_fools_the_brightest_pixel_finder():
    rng = np.random.default_rng(2)
    img = _noisy(_spot(14.2, 25.6, 400), rng)
    img[27, 12] += 600.0                                    # a single hot pixel, brighter than the spot's peak pixel
    stack, regs = img[None], _region()[None]
    old = detect_channel(stack, regs, "brightest", smooth_sigma_px=0.0)[0]
    new = detect_channel(stack, regs, "dog_viterbi")[0]
    assert np.hypot(old[0] - 27, old[1] - 12) < 1.0                     # without smoothing the old method picks the hot pixel
    assert np.hypot(new[0] - 14.2, new[1] - 25.6) < 0.3                  # the 3 x 3 median removes it


def test_uv_coordinates_are_relative_to_the_true_nucleus_centroid_not_the_crop_centre():
    img = _spot(20.4 + 0.0, 20.0 + 0.3, 3000) + 100.0           # the spot sits exactly on the centroid
    fr = pd.DataFrame({"t": [0], "cy_px": [120.4], "cx_px": [80.3], "theta_rad": [0.0], "crop_y0": [100], "crop_x0": [60]})
    cell = CellData(1, fr, _region()[None], {"locus0": img[None].astype(np.uint16)}, S)
    pos = localize_loci(cell, n_loci=1, px_um=0.1)
    assert abs(pos.u_um.iloc[0]) < 0.005 and abs(pos.v_um.iloc[0]) < 0.005
    assert abs(pos.lab_y_px.iloc[0] - 120.4) < 0.02 and abs(pos.lab_x_px.iloc[0] - 80.3) < 0.02


def test_scores_leave_out_loci_that_had_no_light_and_count_detections_there():
    truth = pd.DataFrame({"t": [0, 1, 2, 3], "cell_id": 1, "locus_id": 0, "y_um": 1.0, "x_um": 1.0, "photons": [500, 500, 0, 500]})
    pos = pd.DataFrame({"t": [0, 1, 2], "cell_id": 1, "locus_id": 0, "lab_y_px": [10.0, 10.0, 30.0], "lab_x_px": [10.0, 18.0, 10.0]})
    err, sc = score_localization(pos, truth, 0.1)
    assert sc["n_visible"] == 3 and sc["detection_rate"] == pytest.approx(2 / 3) and sc["n_false_on_dark"] == 1
    assert sc["frac_far"] == pytest.approx(0.5) and sc["rms_error_near_um"] == pytest.approx(0.0)


def test_3d_axial_median_removes_a_single_plane_cosmic_ray():
    rng = np.random.default_rng(7)
    nz, zc = 9, 4.3
    vol = np.stack([_spot(20.2, 19.6, 800 * np.exp(-0.5 * ((k - zc) / 1.3) ** 2)) for k in range(nz)])   # planes finer than the axial PSF
    vol[2] += _spot(9.0, 28.0, 5000, sigma=0.8)                      # a cosmic ray in one plane
    vol = rng.poisson(vol + 20.0).astype(float) + 100.0
    lateral = _first(vol, _region(), median_z=1)
    axial = _first(vol, _region(), median_z=3)
    assert np.hypot(lateral["pos"][1] - 9, lateral["pos"][2] - 28) < 1.0          # lateral median only: the ray is the brightest
    assert np.hypot(axial["pos"][1] - 20.2, axial["pos"][2] - 19.6) < 0.4 and abs(axial["pos"][0] - zc) < 0.5


# ----------------------------------------------------------------------------- dog_viterbi (Chen et al. 2025)


def test_ellipsoid_intensity_is_inner_mean_minus_shell_mean():
    img = np.full((41, 41), 100.0)
    yy, xx = np.indices(img.shape)
    inner = ((yy - 20) / 3.5) ** 2 + ((xx - 20) / 3.5) ** 2 <= 1
    img[inner] += 50.0
    val, n = DV.ellipsoid_intensity(img, (20, 20), (7, 7, 5), (19, 19, 7))
    assert n == inner.sum() and val == pytest.approx(50.0)
    assert DV.ellipsoid_intensity(np.full((9, 41, 41), 7.0), (4, 20, 20))[0] == pytest.approx(0.0)


def test_diffusion_potential_is_recovered_from_a_random_walk():
    rng = np.random.default_rng(8)
    D = 1.5                                                     # px^2 per frame: per-axis step variance 2 D
    walk = np.cumsum(rng.normal(0, np.sqrt(2 * D), (2000, 2)), axis=0)
    est = DV.estimate_D([w[None] for w in walk], [0] * len(walk))
    assert abs(est - D) < 0.1 * D


def test_viterbi_follows_the_spot_past_a_cosmic_ray_and_a_dark_frame():
    rng = np.random.default_rng(9)
    n, centre = 12, np.full((12, 2), 20.0)
    truth = np.cumsum(rng.normal(0, 0.7, (n, 2)), axis=0) + 20
    peaks, snr = [], []
    for t in range(n):
        p, s = [truth[t]], [15.0]
        if t == 5:
            p, s = [np.array([5.0, 33.0])], [400.0]                  # spot not detected; only a bright cosmic ray elsewhere
        p += [rng.uniform(5, 35, 2) for _ in range(2)]               # noise candidates just above the low threshold
        s += [3.5, 3.2]
        order = np.argsort(s)[::-1]
        peaks.append(np.array(p)[order]); snr.append(np.array(s)[order])
    choice, D = DV.viterbi_track(peaks, snr, centre)
    for t in range(n):
        if t == 5:
            assert choice[t] is None or np.hypot(*(peaks[t][choice[t]] - [5, 33])) > 1, "the cosmic ray was taken"
        else:
            assert np.hypot(*(peaks[t][choice[t]] - truth[t])) < 1e-9, t
    assert 0.25 <= D < 3


def test_dog_viterbi_localizes_a_moving_spot_through_noise_and_a_cosmic_ray():
    rng = np.random.default_rng(10)
    n, reg = 9, _region()
    truth = [(20 + 0.2 * np.sin(i), 18 + 0.2 * np.cos(i)) for i in range(n)]
    stack = np.stack([_noisy(_spot(y, x, 600), rng) for y, x in truth])
    stack[4] += _spot(10.0, 27.0, 4000, sigma=0.8)
    out = detect_channel(stack, np.repeat(reg[None], n, 0), "dog_viterbi", np.full((n, 2), S // 2, float))
    for i in range(n):
        assert out[i] is not None and np.hypot(out[i][0] - truth[i][0], out[i][1] - truth[i][1]) < 0.3, i
        assert out[i][5] > 0                                        # the ellipsoid intensity of a real spot is positive


def test_viterbi_keeps_a_faint_spot_among_noise_candidates_instead_of_skipping_frames():
    """Regression: without the clutter term of the track likelihood the path skipped ~30% of frames that held the right (faint) candidate."""
    rng = np.random.default_rng(11)
    n = 40
    truth = 20 + np.cumsum(rng.normal(0, np.sqrt(2 * 1.5), (n, 2)), axis=0) * 0.5
    peaks, snr = [], []
    for t in range(n):
        p = np.vstack([truth[t] + rng.normal(0, 0.3, 2), rng.uniform(0, 40, (3, 2))])
        s = np.r_[4.0, rng.uniform(3.0, 4.5, 3)]
        order = np.argsort(s)[::-1]
        peaks.append(p[order]); snr.append(s[order])
    choice, _ = DV.viterbi_track(peaks, snr, np.full((n, 2), 20.0), search_volume=40 * 40)
    right = sum(c is not None and np.hypot(*(peaks[t][c] - truth[t])) < 1.5 for t, c in enumerate(choice))
    assert right >= n - 1, right                 # without the clutter term: 37 of 40


def test_viterbi_picks_the_spot_up_again_after_an_off_period_longer_than_the_gap_limit():
    """MS2-like: the spot is ON, OFF for 20 frames (> viterbi_max_gap = 5), then ON again elsewhere in the cell; both bursts must be kept."""
    rng = np.random.default_rng(12)
    peaks, snr, truth = [], [], {}
    for t in range(40):
        on = t < 10 or t >= 30
        p, s = [], []
        if on:
            truth[t] = np.array([12.0, 12.0]) if t < 10 else np.array([28.0, 26.0])
            p.append(truth[t] + rng.normal(0, 0.3, 2)); s.append(12.0)
        p += list(rng.uniform(0, 40, (2, 2))); s += [3.2, 3.1]
        o = np.argsort(s)[::-1]
        peaks.append(np.array(p)[o]); snr.append(np.array(s)[o])
    choice, _ = DV.viterbi_track(peaks, snr, np.full((40, 2), 20.0), search_volume=40 * 40)
    for t, x in truth.items():
        assert choice[t] is not None and np.hypot(*(peaks[t][choice[t]] - x)) < 1.5, t


def test_bleaching_is_measured_after_removing_the_background_and_the_naive_estimate_is_biased(tmp_path):
    import tifffile

    from simlive.stage6_analysis.bleach import correct, fit_exponential, nuclear_bleach

    T, dt, tau = 30, 10.0, 600.0
    yy, xx = np.indices((60, 60))
    nuc = np.hypot(yy - 30, xx - 30) < 12
    b = np.exp(-np.arange(T) * dt / tau)
    imgs = np.stack([100.0 + 50.0 + 300.0 * b[t] * nuc for t in range(T)]).astype(np.float32)     # offset 100 + background 50 + bleaching signal
    tifffile.imwrite(tmp_path / "n.tif", imgs, photometric="minisblack")
    tifffile.imwrite(tmp_path / "l.tif", np.repeat(nuc[None], T, 0).astype(np.uint16), photometric="minisblack")
    d = nuclear_bleach(tmp_path / "n.tif", tmp_path / "l.tif", dt)
    assert abs(d.attrs["tau_s"] - tau) < 0.01 * tau
    assert d.attrs["tau_naive_s"] > 1.3 * tau                       # background left in: bleaching looks much weaker
    assert fit_exponential(np.arange(T) * dt, 5 * b)[1] == pytest.approx(tau, rel=1e-6)
    spots = pd.DataFrame({"t": np.arange(T), "intensity": 80.0 * b})
    assert np.allclose(correct(spots, b), 80.0)
