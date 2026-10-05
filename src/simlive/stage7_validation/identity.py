"""Stage 7b (part 1): frame-by-frame identity matching between ground-truth and tracked masks.

For every GT cell in every frame, find the tracked label with the largest overlap. A GT cell is "matched" if that
overlap covers >50% of it (same rule as the CTC matcher). From the sequence of matched track IDs per GT cell we
count identity switches and the fraction of frames the cell carries its majority track ID. These are the quantities
that matter when locus tracks are later attached to cells.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def match_ids(gt: np.ndarray, pred: np.ndarray) -> pd.DataFrame:
    """Rows: t, gt_id, pred_id (0 = unmatched), overlap_frac (of the GT cell's area)."""
    rows = []
    for t in range(len(gt)):
        g, p = gt[t].ravel().astype(np.int64), pred[t].ravel().astype(np.int64)
        sel = g > 0
        g, p = g[sel], p[sel]
        key = g * (int(p.max()) + 1 if p.size else 1) + p
        uniq, cnt = np.unique(key, return_counts=True)
        k = int(p.max()) + 1 if p.size else 1
        gid, pid = uniq // k, uniq % k
        area = dict(zip(*np.unique(g, return_counts=True)))
        best = {}
        for a, b, c in zip(gid, pid, cnt):
            if b > 0 and c > best.get(a, (0, 0))[1]:
                best[a] = (b, c)
        for a, ar in area.items():
            b, c = best.get(a, (0, 0))
            frac = c / ar
            rows.append((t, int(a), int(b) if frac > 0.5 else 0, float(frac)))
    return pd.DataFrame(rows, columns=["t", "gt_id", "pred_id", "overlap_frac"])


def identity_summary(matches: pd.DataFrame) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Returns (totals, per-GT-cell table, switch events table)."""
    per, events = [], []
    for gid, g in matches.sort_values("t").groupby("gt_id"):
        seq = g[g.pred_id > 0]
        ids = seq.pred_id.to_numpy()
        ts = seq.t.to_numpy()
        sw = np.nonzero(ids[1:] != ids[:-1])[0] if len(ids) > 1 else []
        for i in sw:
            events.append((int(ts[i + 1]), int(gid), int(ids[i]), int(ids[i + 1])))
        majority = pd.Series(ids).mode().iloc[0] if len(ids) else 0
        per.append((int(gid), len(g), int((g.pred_id == 0).sum()), len(sw),
                    float((g.pred_id == majority).mean()), int(majority)))
    per = pd.DataFrame(per, columns=["gt_id", "n_frames", "frames_unmatched", "id_switches",
                                     "frac_frames_majority_track", "majority_pred_id"])
    events = pd.DataFrame(events, columns=["t", "gt_id", "old_pred_id", "new_pred_id"])
    totals = {
        "id_switches": int(per.id_switches.sum()),
        "frames_unmatched": int(per.frames_unmatched.sum()),
        "n_gt_cell_frames": int(per.n_frames.sum()),
        "identity_preserved_fraction": float((per.frac_frames_majority_track * per.n_frames).sum() / per.n_frames.sum()),
        "cells_with_any_switch": int((per.id_switches > 0).sum()),
        "n_gt_cells": int(len(per)),
    }
    return totals, per, events
