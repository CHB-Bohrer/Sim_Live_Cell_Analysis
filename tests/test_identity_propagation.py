"""Stage 7 identity matching and stage 7b propagation of tracking errors into locus results, on tiny hand-made label movies with known answers."""
import numpy as np
import pandas as pd

from simlive.stage6_analysis.localize import COLUMNS
from simlive.stage7_validation.identity import identity_summary, match_ids
from simlive.stage7_validation.propagation import propagation_scores, remap_to_truth, true_distances


def movie(T=6, n=2):
    """Two 8x8 squares (true cells 1 and 2) that stay put for T frames."""
    gt = np.zeros((T, 20, 40), np.uint16)
    gt[:, 4:12, 4:12] = 1
    gt[:, 4:12, 24:32] = 2
    return gt


def test_perfect_tracking_has_no_switches_and_full_identity():
    gt = movie()
    m = match_ids(gt, gt)
    totals, per, events = identity_summary(m)
    assert totals["id_switches"] == 0 and totals["frames_unmatched"] == 0
    assert totals["identity_preserved_fraction"] == 1.0 and events.empty


def test_a_label_swap_is_one_switch_per_cell_and_a_lost_cell_is_unmatched():
    gt = movie()
    pred = gt.copy()
    pred[3:] = np.where(gt[3:] == 1, 7, gt[3:])            # cell 1 gets a new label from frame 3
    pred[:, :, 24:32][4:] = 0                              # cell 2 vanishes from frame 4 on
    totals, per, events = identity_summary(match_ids(gt, pred))
    p = per.set_index("gt_id")
    assert p.loc[1, "id_switches"] == 1 and p.loc[2, "id_switches"] == 0
    assert p.loc[2, "frames_unmatched"] == 2
    assert totals["cells_with_any_switch"] == 1
    assert list(events[["gt_id", "old_pred_id", "new_pred_id", "t"]].iloc[0]) == [1, 1, 7, 3]


def test_a_label_covering_under_half_of_the_cell_does_not_count():
    gt = movie(T=1)
    pred = np.zeros_like(gt)
    pred[0, 4:12, 4:7] = 5                                  # 3 of 8 columns = 37.5% of cell 1
    m = match_ids(gt, pred)
    assert (m[m.gt_id == 1].pred_id == 0).all()
    assert np.isclose(m[m.gt_id == 1].overlap_frac.iloc[0], 0.375)


def det(cell_id, t, locus, y, x):
    row = dict.fromkeys(COLUMNS, 0.0)
    row.update(cell_id=cell_id, t=t, locus_id=locus, lab_y_px=y, lab_x_px=x, fit_ok=True)
    return row


def truth_table(T=6):
    rows = []
    for t in range(T):
        for cid, x0 in ((1, 8.0), (2, 28.0)):
            for k, dx in ((0, -1.0), (1, 1.0)):
                rows.append(dict(cell_id=cid, t=t, locus_id=k, y_um=8.0 * 0.1, x_um=(x0 + dx) * 0.1))
    return pd.DataFrame(rows)


def perfect_positions(label_of, T=6):
    rows = []
    for t in range(T):
        for cid, x0 in ((1, 8.0), (2, 28.0)):
            for k, dx in ((0, -1.0), (1, 1.0)):
                rows.append(det(label_of(cid, t), t, k, 8.0, x0 + dx))
    return pd.DataFrame(rows, columns=COLUMNS)


def test_remap_is_the_identity_when_labels_equal_true_ids():
    gt = movie()
    pos = perfect_positions(lambda c, t: c)
    out, book = remap_to_truth(pos, match_ids(gt, gt))
    assert len(out) == len(pos) and book["attributed_fraction"] == 1.0 and book["mixed_track_fraction"] == 0.0
    pd.testing.assert_frame_equal(out.reset_index(drop=True), pos.reset_index(drop=True), check_dtype=False)


def test_identity_switch_is_mixed_but_scores_are_still_correct_per_true_cell():
    gt = movie()
    pred = gt.copy()
    pred[3:] = np.where(gt[3:] == 1, 2, gt[3:])            # label 2 first follows cell 2, then also takes over cell 1... both labelled 2 from frame 3 on
    pred[3:, :, 24:32] = np.where(gt[3:, :, 24:32] == 2, 9, 0)   # ... so give cell 2 label 9 to keep the labels distinct
    matches = match_ids(gt, pred)
    label_of = lambda c, t: {(1, 0): 1, (1, 1): 1, (1, 2): 1, (2, 0): 2, (2, 1): 2, (2, 2): 2}.get((c, t), 2 if c == 1 else 9)
    pos = perfect_positions(label_of)
    err, sc = propagation_scores(pos, matches, truth_table(), 0.1)
    # label 2 stood for cell 2 (frames 0-2) and for cell 1 (frames 3-5): that track is mixed
    assert sc["mixed_track_fraction"] > 0
    assert sc["attributed_fraction"] == 1.0
    # mapped back to the true cells every detection is on its own spot
    assert sc["detection_rate"] == 1.0 and sc["rms_error_um"] < 1e-9
    assert sc["distance_rms_um"] < 1e-9


def test_detections_on_untracked_debris_are_not_attributed():
    gt = movie()
    pos = pd.concat([perfect_positions(lambda c, t: c), pd.DataFrame([det(99, 0, 0, 15.0, 15.0)], columns=COLUMNS)], ignore_index=True)
    out, book = remap_to_truth(pos, match_ids(gt, gt))
    assert len(out) == len(pos) - 1
    assert np.isclose(book["attributed_fraction"], (len(pos) - 1) / len(pos))


def test_cells_lost_by_the_tracker_lower_the_detection_rate():
    gt = movie()
    pred = gt.copy()
    pred[:, :, 24:32] = 0                                   # cell 2 never tracked
    pos = perfect_positions(lambda c, t: c)
    pos = pos[pos.cell_id == 1]                             # (a lost cell has no isolated movie, hence no detections)
    err, sc = propagation_scores(pos, match_ids(gt, pred), truth_table(), 0.1)
    assert np.isclose(sc["detection_rate"], 0.5)


def test_true_distance_from_the_truth_table():
    d = true_distances(truth_table(2))
    assert np.allclose(d.dist_true_um, 0.2)
