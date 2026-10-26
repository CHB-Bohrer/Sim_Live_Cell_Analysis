"""Brightness in the Viterbi spot path (dog_viterbi.IntensityModel, option viterbi_intensity)."""
import numpy as np

from simlive.stage6_analysis import dog_viterbi as DV


def _ar1(n, mu, s, rho, rng):
    e = np.zeros(n)
    e[0] = rng.normal(0, s)
    for i in range(1, n):
        e[i] = rho * e[i - 1] + rng.normal(0, s * np.sqrt(1 - rho ** 2))
    return mu + e


def test_model_recovers_level_spread_and_correlation():
    rng = np.random.default_rng(1)
    t = np.arange(3000)
    I = _ar1(3000, 10.0, 2.0, 0.6, rng) - 0.001 * t          # slow bleaching
    m = DV.IntensityModel(t, I)
    assert abs(m.b + 0.001) < 0.0005
    assert abs(m.s - 2.0) < 0.2
    assert abs(m.rho - 0.6) < 0.05


def test_link_cost_prefers_continuing_brightness():
    m = DV.IntensityModel(np.arange(50), np.r_[np.full(25, 8.0), np.full(25, 12.0)])
    m.a, m.b, m.s, m.rho = 10.0, 0.0, 2.0, 0.9
    c = m.link_cost(0, np.array([13.0]), 1, np.array([13.0, 7.0]))[0]
    assert c[0] < c[1]


def _scene(seed, n=80):
    """One locus (SNR ~ N(10, 1.5)) wandering near the centre, a dim noise peak somewhere in every frame, and in every 4th frame a much
    brighter spot (SNR 40, e.g. a cosmic ray or a bright stray) 1.5 px from the locus."""
    rng = np.random.default_rng(seed)
    pos = np.cumsum(rng.normal(0, 0.7, (n, 2)), axis=0) + 30
    peaks, snr, decoy = [], [], []
    for i in range(n):
        pk = [pos[i] + rng.normal(0, 0.2, 2), rng.uniform(5, 55, 2)]
        sn = [rng.normal(10, 1.5), rng.lognormal(np.log(3.5), 0.3)]
        if i % 4 == 2:
            pk.append(pos[i] + np.array([1.5, 0.0]))
            sn.append(40.0)
            decoy.append(i)
        peaks.append(np.array(pk))
        snr.append(np.array(sn))
    return peaks, snr, np.zeros((n, 2)), decoy


def test_intensity_keeps_the_path_off_a_much_brighter_neighbour():
    off_old, off_new = 0, 0
    for seed in range(5):
        peaks, snr, centre, decoy = _scene(seed)
        kw = dict(search_volume=60 * 60)
        old, _ = DV.viterbi_track(peaks, snr, centre, **kw)
        info = {}
        new, _ = DV.viterbi_track(peaks, snr, centre, viterbi_intensity=True, model_out=info, **kw)
        off_old += sum(old[i] == 2 for i in decoy)
        off_new += sum(new[i] == 2 for i in decoy)
        assert sum(c == 0 for c in new) >= 0.95 * len(new)          # the locus is still followed
        assert 8 < info["a"] < 12
    assert off_old >= 10                                             # the log-SNR reward takes the bright neighbour (18 of 100 decoy frames) ...
    assert off_new <= 3                                              # ... the brightness model hardly (2: frames where the locus was 3 sd dim)


def test_off_by_default_is_unchanged():
    peaks, snr, centre, _ = _scene(3)
    a, Da = DV.viterbi_track(peaks, snr, centre, search_volume=3600)
    b, Db = DV.viterbi_track(peaks, snr, centre, search_volume=3600, viterbi_intensity=False)
    assert a == b and Da == Db
