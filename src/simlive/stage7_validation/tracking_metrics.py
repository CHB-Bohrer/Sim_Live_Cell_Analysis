"""Stage 7a: score tracking against ground truth with Cell Tracking Challenge metrics (via traccuracy).

Ground truth and results are written in CTC format (mask%04d.tif + man_track.txt) then matched with the CTC
matcher (a GT cell matches a predicted cell if the prediction covers >50% of it).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import tifffile


def gt_table(cells: pd.DataFrame, masks: np.ndarray | None = None) -> pd.DataFrame:
    """CTC track table of the true cells. With `masks`, only the frames in which a cell is actually in the image count (cells can walk
    out of and into the field: geometry.outside_um); cells never in view are dropped and a parent never in view is set to 0."""
    if masks is not None:
        seen = [(t, int(c)) for t, m in enumerate(masks) for c in np.unique(m) if c]
        vis = pd.DataFrame(seen, columns=["t", "cell_id"])
        cells = cells.merge(vis, on=["t", "cell_id"])
    g = cells.groupby("cell_id").agg(t1=("t", "min"), t2=("t", "max"), parent=("parent_id", "first")).reset_index()
    g.loc[~g.parent.isin(g.cell_id), "parent"] = 0
    return g.rename(columns={"cell_id": "label"})[["label", "t1", "t2", "parent"]]


def ctc_ground_truth(masks: np.ndarray, cells: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    """True masks + track table in CTC form when cells may leave the field and come back: every stretch of frames in which a cell is in
    view becomes its own track; a later stretch gets a new label whose parent is the stretch before (the CTC way to write a gap).
    Without gaps this is gt_table(cells, masks) and the masks unchanged."""
    tab = gt_table(cells, masks)
    present = {}
    for t, m in enumerate(masks):
        for c in np.unique(m):
            if c:
                present.setdefault(int(c), []).append(t)
    if all(fr[-1] - fr[0] + 1 == len(fr) for fr in present.values()):
        return masks, tab
    out = masks.copy()
    rows, next_id = [], int(max(present)) + 1
    parent_of = dict(zip(tab.label, tab.parent))
    for c, fr in sorted(present.items()):
        fr = np.array(fr)
        starts = np.r_[0, np.nonzero(np.diff(fr) > 1)[0] + 1]
        ends = np.r_[starts[1:], len(fr)]
        prev = None
        for s, e in zip(starts, ends):
            lab = c if prev is None else next_id
            if prev is not None:
                next_id += 1
                for t in fr[s:e]:
                    out[t][masks[t] == c] = lab
            rows.append((lab, int(fr[s]), int(fr[e - 1]), int(parent_of.get(c, 0)) if prev is None else prev))
            prev = lab
    return out, pd.DataFrame(rows, columns=["label", "t1", "t2", "parent"])


def write_ctc(masks: np.ndarray, table: pd.DataFrame, out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    for t, m in enumerate(masks):
        tifffile.imwrite(out / f"mask{t:04d}.tif", m.astype(np.uint16))
    table.astype(int).to_csv(out / "man_track.txt", index=False, header=False, sep=" ")
    return out


def evaluate(gt_dir: Path, pred_dir: Path) -> dict:
    import logging

    logging.getLogger("traccuracy").setLevel(logging.CRITICAL)  # its log calls use a format stdlib rejects
    from traccuracy import run_metrics
    from traccuracy.loaders import load_ctc_data
    from traccuracy.matchers import CTCMatcher
    from traccuracy.metrics import CHOTAMetric, CTCMetrics, TrackOverlapMetrics

    gt = load_ctc_data(str(gt_dir), str(gt_dir / "man_track.txt"), name="gt")
    pred = load_ctc_data(str(pred_dir), str(pred_dir / "man_track.txt"), name="pred")
    results, _ = run_metrics(gt, pred, CTCMatcher(), [CTCMetrics(), CHOTAMetric(), TrackOverlapMetrics()])
    out = {}
    for r in results:
        out[r["metric"]["name"]] = r["results"]
    return out
