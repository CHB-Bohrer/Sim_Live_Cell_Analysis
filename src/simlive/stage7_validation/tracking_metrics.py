"""Stage 7a: score tracking against ground truth with Cell Tracking Challenge metrics (via traccuracy).

Ground truth and results are written in CTC format (mask%04d.tif + man_track.txt) then matched with the CTC
matcher (a GT cell matches a predicted cell if the prediction covers >50% of it).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import tifffile


def gt_table(cells: pd.DataFrame) -> pd.DataFrame:
    g = cells.groupby("cell_id").agg(t1=("t", "min"), t2=("t", "max"), parent=("parent_id", "first"))
    return g.reset_index().rename(columns={"cell_id": "label"})[["label", "t1", "t2", "parent"]]


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
