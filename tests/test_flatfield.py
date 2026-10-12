"""Flat-field + per-cell bleaching estimated from moving cells alone (stage6_analysis/flatfield.py), and the per-cell bleaching of the renderer."""
import numpy as np
import pandas as pd

from simlive.stage6_analysis.flatfield import FlatFit, correct_values, fit_flat_bleach


def _synthetic(seed=0, shape=(300, 300), T=120, n=12, tau_frames=300.0, noise=0.02):
    rng = np.random.default_rng(seed)
    true = FlatFit(shape, theta=np.array([0.04, -0.03, -0.08, 0.01, -0.07]), k_per_frame=1 / tau_frames, dt_s=10.0)
    rows = []
    for c in range(1, n + 1):
        t0 = 0 if c <= n - 3 else int(rng.integers(30, 60))      # the last three cells walk in later, unbleached
        p = rng.uniform(30, shape[0] - 30, 2)
        a = np.exp(rng.normal(0, 0.3))
        for t in range(t0, T):
            p = np.clip(p + rng.normal(0, 4, 2), 5, shape[0] - 5)
            rows.append((c, t, p[0], p[1], a))
    tr = pd.DataFrame(rows, columns=["cell_id", "t", "y_px", "x_px", "a"])
    f = true.F(tr.y_px, tr.x_px)
    b = true.bleach(tr).to_numpy()
    tr["value"] = tr.a * f * b * np.exp(rng.normal(0, noise, len(tr)))
    return tr, true


def test_flat_and_bleaching_are_recovered_from_moving_cells():
    tr, true = _synthetic()
    fit = fit_flat_bleach(tr, tr, true.shape, dt_s=10.0)
    d = fit.map() - true.map()
    assert np.sqrt(np.mean(d ** 2)) < 0.01 and np.abs(d).max() < 0.04, (np.sqrt(np.mean(d ** 2)), np.abs(d).max())  # true range ~0.85-1.07
    assert abs(fit.tau_s / true.tau_s - 1) < 0.1, (fit.tau_s, true.tau_s)
    # corrected values are flat in time and space: only the cell's own brightness (a) is left
    df = tr.rename(columns={"y_px": "lab_y_px", "x_px": "lab_x_px"}).assign(intensity=tr.value)
    corr = correct_values(df, fit, tr)
    resid = np.log(corr / tr.a)
    assert resid.std() < 0.03 and abs(resid.mean()) < 0.03
    # a uniform-illumination, global-bleaching model leaves the illumination in
    naive = np.log(tr.value / tr.a / np.exp(-tr.t / 300.0))
    assert naive.std() > 2 * resid.std()


def test_without_flat_term_the_fixed_dose_flat_is_used():
    tr, true = _synthetic(seed=1)
    tr["value"] = tr.value / true.F(tr.y_px, tr.x_px)              # as if a calibration had already divided out the flat-field
    fit = fit_flat_bleach(tr, tr, true.shape, 10.0, fit_flat=False, dose_flat=true)
    assert abs(fit.tau_s / true.tau_s - 1) < 0.1


def test_renderer_bleaches_each_cell_by_its_own_dose():
    from simlive.stage3_microscopy.errors import ErrorSet
    from simlive.stage3_microscopy.render import Imager

    cfg = {"optics": {"pixel_size_nm": 1000, "photons_per_px_s": 100, "bleach_tau_s": 100.0, "NA": 1.2, "wavelength_nm": 509,
                      "background_photons": 0, "quantum_efficiency": 1, "read_noise_e": 0, "gain_adu_per_e": 1, "offset_adu": 0},
           "acquisition": {"n_frames": 20, "frame_interval_s": 10, "exposure_s": 0.1}, "geometry": {"fov_um": [40, 40]},
           "imaging_errors": {"enabled": True, "bleaching": {"enabled": True, "per_cell": True},
                              **{s: {"enabled": False} for s in ("flat_field", "stage_drift", "flicker")}}}
    rows = []
    for t in range(20):        # cell 1 sits in view all the time; cell 2 is outside the field until frame 10
        rows.append((t, 1, 20.0, 20.0))
        rows.append((t, 2, -30.0 if t < 10 else 10.0, 10.0))
    cells = pd.DataFrame(rows, columns=["t", "cell_id", "y_um", "x_um"]).assign(
        radius_um=3.0, aspect=1.0, angle_rad=0.0, bound_radius_um=3.0, parent_id=0)
    im = Imager(cfg, cells, errors=ErrorSet(cfg))
    d = im.cell_bleach("nucleus")
    assert np.isclose(d[(15, 1)][0], np.exp(-10 * 15 / 100.0))     # lit in all 15 earlier frames
    assert np.isclose(d[(15, 2)][0], np.exp(-10 * 5 / 100.0))      # lit only since it entered at frame 10
    assert d[(5, 2)] == (1.0, False)
    # with the switch off the global model is unchanged
    cfg2 = {**cfg, "imaging_errors": {**cfg["imaging_errors"], "bleaching": {"enabled": True}}}
    assert np.isclose(Imager(cfg2, cells, errors=ErrorSet(cfg2)).bleach("nucleus", 15), np.exp(-1.5))
