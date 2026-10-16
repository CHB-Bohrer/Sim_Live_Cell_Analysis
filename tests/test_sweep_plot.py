"""Roles of the settings in a sweep plot, and that the figure has one panel / line per value and a legend outside the axes."""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "app"))
import sweep_plot as sp  # noqa: E402


def sweep(seeds=3):
    rows = []
    rng = np.random.default_rng(0)
    for sb in (0.05, 0.4, 1.0, 2.0):
        for mod in ("biplane", "zstack"):
            for fr in (10.0, 20.0):
                for cor in (False, True):
                    for s in range(seeds):
                        rows.append(dict(ms2_sb=sb, modality=mod, frame_s=fr, corrected=cor, seed=s, run="r", sensitivity=0.5 + 0.1 * sb + 0.02 * rng.normal()))
    return pd.DataFrame(rows)


SETTINGS = ["modality", "ms2_sb", "frame_s", "corrected"]


def test_default_roles_put_the_most_varied_setting_on_x_and_keep_the_grid_small():
    r = sp.default_roles(sweep(), SETTINGS)
    assert r["x"] == "ms2_sb" and r["lines"] == ["modality"]
    assert set(r["panels"]) == {"frame_s", "corrected"} and r["held"] == {}      # 2 x 2 = 4 panels (<= 6)


def test_roles_cover_every_setting_exactly_once():
    df = sweep()
    for x in SETTINGS:
        r = sp.default_roles(df, SETTINGS, x)
        used = [r["x"], *r["lines"], *r["panels"], *r["held"]]
        assert sorted(used) == sorted(SETTINGS) and len(set(used)) == len(SETTINGS)


def test_too_many_panels_are_held_fixed_instead():
    rows = [dict(a=a, b=b, c=c, d=d, seed=0, run="r", y=1.0) for a in range(2) for b in range(5) for c in range(5) for d in range(5)]
    r = sp.default_roles(pd.DataFrame(rows), ["a", "b", "c", "d"])
    n = {"a": 2, "b": 5, "c": 5, "d": 5}
    assert np.prod([n[c] for c in r["panels"]]) <= sp.MAX_PANELS
    assert r["held"]


def test_figure_has_one_panel_per_combination_and_the_legend_outside():
    df = sweep()
    r = sp.default_roles(df, SETTINGS)
    fig = sp.draw(df, "ms2_sb", "sensitivity", "sensitivity", lines=r["lines"], panels=r["panels"], held=r["held"], hlines=[(0.6, "reference", "--")])
    axes = [a for a in fig.axes if a.axison]
    assert len(axes) == 4
    assert all(len(a.lines) >= 2 for a in axes)                       # two modality lines (+ the reference line)
    leg = fig.legends[0]
    fig.canvas.draw()
    assert leg.get_window_extent().x0 >= max(a.get_window_extent().x1 for a in axes) - 1       # to the right of every panel
    assert {t.get_text() for t in leg.get_texts()} == {"modality=biplane", "modality=zstack", "reference"}


def test_held_settings_filter_the_data_and_are_named_in_the_title():
    df = sweep()
    fig = sp.draw(df, "ms2_sb", "sensitivity", "sensitivity", lines=["modality"], panels=[], held={"frame_s": 10.0, "corrected": True})
    assert "frame_s=10.0" in fig._suptitle.get_text() and "corrected=True" in fig._suptitle.get_text()
    ax = fig.axes[0]
    n_points = sum(len(c.get_offsets()) for c in ax.collections if type(c).__name__ == "PathCollection")
    assert n_points == 4 * 2 * 3                                      # 4 x values x 2 modalities x 3 seeds


def test_sweeps_that_are_not_full_grids_get_held_values_that_exist_together():
    rows = [dict(a=a, b=b, seed=s, run="r", y=1.0) for a, b in ((1, 10), (1, 20), (2, 10)) for s in range(3)]    # (2, 20) was never run
    rows += [dict(a=2, b=30, seed=s, run="r", y=1.0) for s in range(3)]
    df = pd.DataFrame(rows)
    r = sp.default_roles(df, ["a", "b"], x="b")
    assert r["lines"] == ["a"] or r["panels"] == ["a"] or r["held"]
    fig = sp.draw(df, "b", "y", "y", lines=r["lines"], panels=r["panels"], held=r["held"])      # must not raise
    assert fig.axes
    fig2 = sp.draw(df, "b", "y", "y", lines=[], panels=[], held={"a": 99})                       # nothing matches: a message, not an error
    assert not fig2.axes[0].axison
