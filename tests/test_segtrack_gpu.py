"""Cellpose and Trackastra really run (until 2026-10-20 only their registry was tested; they ran only inside full pipeline runs).

A small synthetic movie: four bright, well separated discs drifting a little each frame, with camera-like noise. Cellpose (in its worker
process) must find every disc in every frame, and Trackastra must keep one track per disc on the true masks and on Cellpose's masks.
"""
import numpy as np
import pytest

pytestmark = pytest.mark.gpu

T, H, W, R = 8, 192, 192, 16
CENTRES = np.array([[50, 50], [50, 140], [140, 50], [140, 140]], float)
STEP = np.array([[1.0, 0.5], [0.5, -1.0], [-1.0, 0.5], [-0.5, -1.0]])


def _movie(seed: int = 0):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[:H, :W]
    img = np.zeros((T, H, W), np.float32)
    lab = np.zeros((T, H, W), np.uint16)
    for t in range(T):
        for k, (cy, cx) in enumerate(CENTRES + t * STEP):
            inside = (yy - cy) ** 2 + (xx - cx) ** 2 <= R ** 2
            lab[t][inside] = k + 1
            img[t][inside] += 400.0
    img = rng.poisson(img + 20.0).astype(np.float32) + rng.normal(0, 2.0, img.shape).astype(np.float32)
    return img, lab


def _best_iou(pred: np.ndarray, true: np.ndarray, k: int) -> tuple[float, int]:
    m = true == k
    ids, n = np.unique(pred[m], return_counts=True)
    keep = ids > 0
    if not keep.any():
        return 0.0, 0
    p = int(ids[keep][np.argmax(n[keep])])
    inter = np.logical_and(m, pred == p).sum()
    return inter / np.logical_or(m, pred == p).sum(), p


def _one_track_per_disc(tracked: np.ndarray, lab: np.ndarray) -> list[set]:
    return [{_best_iou(tracked[t], lab[t], k)[1] for t in range(T)} for k in range(1, len(CENTRES) + 1)]


def test_cellpose_finds_every_disc_and_trackastra_keeps_identities():
    from simlive.stage4_segtrack.segment import cellpose_seg
    from simlive.stage4_segtrack.trackers import make_tracker

    img, lab = _movie()
    seg = cellpose_seg(img, diameter_px=2 * R)
    assert seg.shape == lab.shape
    ious = [_best_iou(seg[t], lab[t], k)[0] for t in range(T) for k in range(1, len(CENTRES) + 1)]
    assert min(ious) > 0.7, f"Cellpose missed or badly cut a disc (IoU min {min(ious):.2f})"
    assert all(len(np.unique(seg[t])) - 1 == len(CENTRES) for t in range(T)), "Cellpose found extra objects"

    tracker = make_tracker("trackastra", model="general_2d", mode="greedy")
    for masks in (lab, seg):
        res = tracker.track(img, masks)
        assert res.masks.shape == lab.shape
        ids = _one_track_per_disc(res.masks, lab)
        assert all(len(s) == 1 and 0 not in s for s in ids), f"identity switch or lost disc: {ids}"
        assert len(set().union(*ids)) == len(CENTRES), f"two discs share a track: {ids}"
