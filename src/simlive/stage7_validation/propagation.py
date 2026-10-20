"""Stage 7b: how tracking errors reach the locus results.

The locus analysis runs once per identity source: `truth` (the simulation's own cell IDs) and `tracked_*` (the tracker's IDs). A tracker ID is
not a cell: after an identity switch one tracked "cell" carries the loci of two true cells, and a cell the tracker lost has no locus track at all.
`remap_to_truth` turns the tracked results into per-TRUE-cell results by following the frame-by-frame identity matching of stage 7
(`identity.match_ids`: for every true cell and frame, the tracked label that covers most of it), so the usual scores
(`localize.score_localization`) can be computed through the tracker's IDs and compared with the same scores on the true IDs.

What is reported (`propagation_scores`):
  - `attributed_fraction`: share of the tracked detections that sit on a tracked label matched to a true cell (the rest are debris / false cells);
  - `mixed_track_fraction`: share of the detections whose tracked label covers more than one true cell over the movie (an identity switch or a
    merged pair: such a "cell" has a locus trace made of two cells, which a per-cell analysis will treat as one);
  - the localization scores of the remapped detections (detection rate, errors, intensity correlation), so cells the tracker lost show up as
    a lower detection rate;
  - `distance_rms_um` / `distance_median_um`: error of the locus-locus distance per cell and frame against the true distance, where both loci were found.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from simlive.stage1_chromatin.copies import with_copy
from simlive.stage6_analysis.localize import locus_distances, match_copies, score_localization

VARIANT_OF_SOURCE = {"tracked_gt_masks": "gt_masks", "tracked_auto_masks": "auto_masks"}


def _pred_to_gt(matches: pd.DataFrame) -> pd.DataFrame:
    """One true cell per (frame, tracked label): the one that label overlaps most (a label can be the best match of two merged cells)."""
    m = matches[matches.pred_id > 0].sort_values("overlap_frac", ascending=False)
    return m.drop_duplicates(["t", "pred_id"])[["t", "pred_id", "gt_id"]]


def remap_to_truth(pos: pd.DataFrame, matches: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Locus detections of tracked cells -> detections of true cells. `pos.cell_id` is the tracker's label, `matches` the table of `match_ids`
    (t, gt_id, pred_id, overlap_frac). Returns (detections with cell_id = the true cell, bookkeeping)."""
    if pos.empty:
        return pos.copy(), {"n_detections": 0, "attributed_fraction": float("nan"), "mixed_track_fraction": float("nan")}
    p2g = _pred_to_gt(matches)
    j = pos.merge(p2g, left_on=["t", "cell_id"], right_on=["t", "pred_id"], how="left")
    n_gt_per_label = p2g.groupby("pred_id").gt_id.nunique()                # true cells a tracked label ever stood for
    mixed = j.cell_id.map(n_gt_per_label).fillna(0) > 1
    ok = j.gt_id.notna()
    out = j[ok].drop(columns=["pred_id", "cell_id"]).rename(columns={"gt_id": "cell_id"})
    out["cell_id"] = out.cell_id.astype(int)
    out = out[pos.columns]
    return out, {"n_detections": int(len(pos)), "attributed_fraction": float(ok.mean()),
                 "mixed_track_fraction": float((mixed & ok).sum() / max(ok.sum(), 1)),
                 "n_tracked_labels": int(pos.cell_id.nunique()), "n_mixed_labels": int((n_gt_per_label[n_gt_per_label.index.isin(pos.cell_id.unique())] > 1).sum())}


def true_distances(truth: pd.DataFrame) -> pd.DataFrame:
    """True locus-locus distance (um) per cell, copy and frame from the truth table (y_um, x_um of locus 0 and 1)."""
    truth = with_copy(truth)
    a = truth[truth.locus_id == 0].set_index(["cell_id", "t", "copy"])
    b = truth[truth.locus_id == 1].set_index(["cell_id", "t", "copy"])
    j = a.join(b, lsuffix="_a", rsuffix="_b", how="inner")
    return np.hypot(j.y_um_a - j.y_um_b, j.x_um_a - j.x_um_b).rename("dist_true_um").reset_index()


def propagation_scores(pos: pd.DataFrame, matches: pd.DataFrame, truth: pd.DataFrame, px_um: float) -> tuple[pd.DataFrame, dict]:
    """Scores of tracked-ID locus detections `pos` against the true loci, through the identity matching. Returns (errors, scores)."""
    remapped, book = remap_to_truth(pos, matches)
    if remapped.empty:
        return pd.DataFrame(columns=["cell_id", "t", "locus_id", "err_um", "visible"]), book
    remapped = match_copies(remapped, truth, px_um)          # the tracked cell's copies -> the true cell's copies
    err, sc = score_localization(remapped, truth, px_um)
    sc = {**book, **sc}
    d = locus_distances(remapped, px_um)
    d = d.merge(true_distances(truth), on=["cell_id", "t", "copy"], how="inner")
    if len(d):
        e = d.dist_um - d.dist_true_um
        sc.update({"distance_rms_um": float(np.sqrt((e ** 2).mean())), "distance_median_um": float(e.abs().median()), "n_distances": int(len(d))})
    return err, sc
