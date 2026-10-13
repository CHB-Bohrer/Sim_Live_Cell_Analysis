"""Instrument presets (configs/instruments/) and the cost of sequential colours (the loci move between colour exposures)."""
import numpy as np
import pandas as pd

from simlive.io.runs import REPO_ROOT, load_config


def test_instrument_preset_is_merged_and_explicit_settings_win():
    c = load_config(REPO_ROOT / "configs" / "loci_demo.yaml", ["instrument=vutara_vxl", "optics.NA=1.25"])
    assert c["optics"]["pixel_size_nm"] == 100 and c["optics"]["quantum_efficiency"] == 0.95
    assert c["optics"]["NA"] == 1.25                                          # --set beats the preset
    assert c["imaging_errors"]["sequential_colour"]["enabled"] and c["imaging_errors"]["crosstalk"]["pairs"][0][2] < 0.001
    assert c["imaging_errors"]["hot_pixels"]["enabled"]                      # the config's other settings are kept
    assert "instrument_info" not in c


def test_sequential_colours_move_the_later_colours_by_the_chromatin_motion():
    from simlive.stage3_microscopy.errors import ErrorSet
    from simlive.stage3_microscopy.render import Imager

    rng = np.random.default_rng(0)
    T, dt = 200, 10.0
    cfg = {"optics": {"pixel_size_nm": 100, "photons_per_px_s": 100, "NA": 1.2, "wavelength_nm": 509, "background_photons": 0,
                      "quantum_efficiency": 1, "read_noise_e": 0, "gain_adu_per_e": 1, "offset_adu": 0},
           "acquisition": {"n_frames": T, "frame_interval_s": dt, "exposure_s": 0.1}, "geometry": {"fov_um": [20, 20]},
           "loci": {"n_loci": 2, "wavelengths_nm": [575, 665], "photons_per_locus_s": 1000},
           "imaging_errors": {"enabled": True, "sequential_colour": {"enabled": True, "order": ["locus0", "nucleus", "locus1"], "delay_s": 1.0}}}
    cells = pd.DataFrame([(t, 1, 10.0, 10.0 + 0.01 * t, 3.0, 1.0, 0.0, 3.0, 0) for t in range(T)],
                         columns=["t", "cell_id", "y_um", "x_um", "radius_um", "aspect", "angle_rad", "bound_radius_um", "parent_id"])
    walk = np.cumsum(rng.normal(0, 0.05, (T, 2)), axis=0)                 # locus random walk relative to the cell, 50 nm / frame
    truth = pd.DataFrame([(t, 1, k, 10.0 + walk[t, 0], 10.0 + 0.01 * t + walk[t, 1]) for t in range(T) for k in (0, 1)],
                         columns=["t", "cell_id", "locus_id", "y_um", "x_um"])
    im = Imager(cfg, cells, truth, errors=ErrorSet(cfg))
    assert im.colour_delay_s("locus0") == 0 and im.colour_delay_s("locus1") == 2.0
    d0 = np.array([np.r_[im.locus_positions_at_exposure("locus0", t)] for t in range(T - 1)])
    d1 = np.array([np.r_[im.locus_positions_at_exposure("locus1", t)] for t in range(T - 1)])
    exact = truth[truth.locus_id == 0].sort_values("t")[["y_um", "x_um"]].to_numpy()[: T - 1] / 0.1
    assert np.allclose(d0, exact)                                           # the first colour: at the frame time
    extra = d1 - exact
    f = 2.0 / dt
    sd_expected = 0.5 * np.sqrt(f * (1 - f))                               # 0.05 um = 0.5 px per frame per axis, Brownian bridge
    assert abs((extra - f * np.diff(np.r_[exact, exact[-1:]], axis=0)).std() / sd_expected - 1) < 0.25
