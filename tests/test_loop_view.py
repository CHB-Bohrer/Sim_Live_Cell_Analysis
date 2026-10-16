"""Helpers behind the dashboard's loop-extruder graph: coverage per genomic bin over time, and per-factor traces for a zoom window."""
import numpy as np

from simlive.stage1_chromatin.analysis import loop_coverage, loop_traces


def lef_array():
    # 3 blocks, 2 factors: factor 0 spans 0..40 kb, 0..60, 0..80; factor 1 has not separated (size 0) at block 0, then 100..120
    lef = np.zeros((3, 2, 2), np.int32)
    lef[:, 0, 1] = [40, 60, 80]
    lef[:, 1] = [[500, 500], [100, 120], [100, 120]]
    return lef


def test_coverage_counts_loops_per_bin_and_ignores_zero_size_loops():
    cov = loop_coverage(lef_array(), 200, bin_kb=20, max_time_bins=3)
    assert cov.shape == (3, 10)
    assert list(cov[0]) == [1, 1, 1, 0, 0, 0, 0, 0, 0, 0]          # 0..40 kb = bins 0-2; the zero-size loop at 500 is out of range / not counted
    assert list(cov[1][:7]) == [1, 1, 1, 1, 0, 1, 1]               # 0..60 plus the second factor over bins 5-6
    assert cov.sum(axis=1).tolist() == [3, 6, 7 + 0]               # bins covered per block


def test_coverage_averages_blocks_into_time_bins():
    lef = np.tile(lef_array()[:1], (100, 1, 1))
    cov = loop_coverage(lef, 200, bin_kb=20, max_time_bins=10)
    assert cov.shape == (10, 10) and np.allclose(cov[:, :3], 1.0)


def test_traces_pick_factors_in_the_window_and_break_lines_at_reloads():
    lef = np.zeros((6, 3, 2), np.int32)
    lef[:, 0] = [[10, 12], [9, 13], [8, 14], [300, 300], [299, 301], [298, 302]]      # factor 0 reloads from ~10 kb to 300 kb at block 3
    lef[:, 1] = [[500, 520]] * 6                                                      # never near the window
    lef[:, 2] = [[0, 2000]] * 6                                                       # spans the whole window
    tr = loop_traces(lef, 0, 6, 0, 100, jump=3)
    assert sorted(t["factor"] for t in tr) == [0, 2]
    f0 = next(t for t in tr if t["factor"] == 0)
    assert np.isnan(f0["left"][3]) and np.isnan(f0["right"][3]) and not np.isnan(f0["left"][2])
    assert np.nanmax(f0["left"][:3]) == 10 and np.nanmin(f0["left"][4:]) == 298      # both stretches kept, only the jump itself is blanked
    assert len(loop_traces(lef, 0, 6, 0, 100, max_factors=1)) == 1
