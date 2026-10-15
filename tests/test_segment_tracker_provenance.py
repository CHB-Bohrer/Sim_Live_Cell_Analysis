"""Segmentation baseline, the tracker registry/adapter protocol, and run provenance (cheap CPU tests; Trackastra and Cellpose themselves need the GPU stack)."""
import json

import numpy as np
import pytest

from simlive.io import runs
from simlive.stage4_segtrack import trackers
from simlive.stage4_segtrack.segment import SEGMENTERS, segment


def disks(T=3, centres=((20, 20), (20, 60)), r=9, shape=(40, 80)):
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    im = np.zeros((T, *shape), np.float32)
    for cy, cx in centres:
        im += ((yy - cy) ** 2 + (xx - cx) ** 2 < r * r)
    rng = np.random.default_rng(0)
    return im * 100 + rng.normal(size=im.shape).astype(np.float32) * 2


def test_threshold_watershed_finds_each_separate_nucleus_in_every_frame():
    lab = segment(disks(), "threshold_watershed", smooth_px=1.0, min_area_px=30, min_distance_px=5)
    assert lab.shape == (3, 40, 80)
    for frame in lab:
        assert len(np.unique(frame[frame > 0])) == 2
        assert frame[20, 20] > 0 and frame[20, 60] > 0 and frame[20, 20] != frame[20, 60]
        assert frame[0, 0] == 0


def test_unknown_segmenter_is_rejected():
    with pytest.raises(ValueError, match="Unknown segmenter"):
        segment(disks(), "nope")
    assert {"threshold_watershed", "cellpose"} <= set(SEGMENTERS)


def test_unknown_tracker_is_rejected_and_the_registry_has_trackastra():
    with pytest.raises(ValueError, match="Unknown tracker"):
        trackers.make_tracker("nope")
    assert "trackastra" in trackers.TRACKERS


def test_a_new_tracker_registers_and_runs_through_make_tracker(monkeypatch):
    class Identity:
        name = "identity"

        def __init__(self, **_):
            pass

        def track(self, images, masks):
            import pandas as pd
            return trackers.TrackResult(masks.copy(), pd.DataFrame({"label": [1], "t1": [0], "t2": [len(masks) - 1], "parent": [0]}))

    monkeypatch.setitem(trackers.TRACKERS, "identity", Identity)
    res = trackers.make_tracker("identity").track(np.zeros((2, 4, 4)), np.ones((2, 4, 4), np.uint16))
    assert res.masks.shape == (2, 4, 4) and list(res.table.columns) == ["label", "t1", "t2", "parent"]


def test_trackastra_adapter_keeps_its_settings_without_importing_torch():
    t = trackers.make_tracker("trackastra", model="general_2d", mode="ilp")
    assert (t.model_name, t.mode) == ("general_2d", "ilp")


def test_provenance_records_commit_and_whether_the_tree_was_dirty(tmp_path):
    runs.write_provenance(tmp_path, {"a": 1}, 7)
    p = json.loads((tmp_path / "provenance.json").read_text())
    assert (tmp_path / "seed.txt").read_text() == "7" and (tmp_path / "params.yaml").exists()
    assert "git_commit" in p and isinstance(p["git_dirty"], bool)
    if p["git_dirty"]:
        assert len(p["git_diff_sha256"]) == 64
    else:
        assert p["git_diff_sha256"] == ""
