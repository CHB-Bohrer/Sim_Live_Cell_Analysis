"""Stage 5a: isolate every cell through time (one self-contained single-cell movie per cell ID).

Input : image stacks (nucleus.tif, optional locus.tif; (T,Y,X)) + a label stack whose value is the cell/track ID
        (ground-truth IDs, or tracker output) -> the *identity source*.
Output: <out>/cell_table.csv   one row per (cell, frame): centroid, area, orientation, crop origin, border/neighbour flags
        <out>/cell_summary.csv one row per cell: lifetime, gaps, border frames, nearest-neighbour distance
        <out>/cells/cell_XXXX/{nucleus,locus,mask}.tif + frames.csv   crops of one fixed size S x S centred on the
        cell's centroid in every frame it exists. `mask` marks only that cell's pixels (neighbours are zero), so
        analysis can ignore neighbouring cells' loci that fall inside the crop.

Memory: stacks are memory-mapped (never loaded whole), each worker reads only its own cell's pixels, crops are saved
as uint16 on disk, and the pool size is capped (`n_workers`). Peak per-worker memory is a few crops, not a movie.
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile


def _memmap(path):
    try:
        return tifffile.memmap(path, mode="r")
    except Exception:  # compressed / non-contiguous files can't be memory-mapped
        return tifffile.imread(path)


def build_cell_table(labels_path: str | Path) -> pd.DataFrame:
    """Per-frame geometry of every labeled cell, computed one frame at a time (low memory)."""
    labels = _memmap(labels_path)
    T, ny, nx = labels.shape
    rows = []
    for t in range(T):
        lab = np.asarray(labels[t])
        ys, xs = np.nonzero(lab)
        if ys.size == 0:
            continue
        ids = lab[ys, xs].astype(np.int64)
        uid, inv = np.unique(ids, return_inverse=True)
        area = np.bincount(inv).astype(float)
        cy, cx = np.bincount(inv, ys) / area, np.bincount(inv, xs) / area
        dy, dx = ys - cy[inv], xs - cx[inv]
        mu20, mu02, mu11 = (np.bincount(inv, dx * dx) / area, np.bincount(inv, dy * dy) / area,
                            np.bincount(inv, dx * dy) / area)
        theta = 0.5 * np.arctan2(2 * mu11, mu20 - mu02)  # major axis, from +x towards +y (same as renderer)
        rmax = np.zeros(len(uid)); np.maximum.at(rmax, inv, np.hypot(dy, dx))
        ymin = np.full(len(uid), ny); ymax = np.full(len(uid), -1); xmin = np.full(len(uid), nx); xmax = np.full(len(uid), -1)
        np.minimum.at(ymin, inv, ys); np.maximum.at(ymax, inv, ys); np.minimum.at(xmin, inv, xs); np.maximum.at(xmax, inv, xs)
        border = (ymin == 0) | (xmin == 0) | (ymax == ny - 1) | (xmax == nx - 1)
        if len(uid) > 1:
            d = np.hypot(cy[:, None] - cy[None], cx[:, None] - cx[None]); np.fill_diagonal(d, np.inf); nn = d.min(1)
        else:
            nn = np.full(1, np.inf)
        for k, cid in enumerate(uid):
            rows.append((t, int(cid), cy[k], cx[k], area[k], theta[k], rmax[k], bool(border[k]), nn[k]))
    return pd.DataFrame(rows, columns=["t", "cell_id", "cy_px", "cx_px", "area_px", "theta_rad", "rmax_px",
                                       "touches_border", "nn_dist_px"])


def _crop(frame, y0: int, x0: int, S: int, pad_value):
    """S x S window with top-left corner (y0, x0); outside the image is filled with pad_value."""
    ny, nx = frame.shape
    out = np.full((S, S), pad_value, dtype=frame.dtype)
    ya, yb, xa, xb = max(y0, 0), min(y0 + S, ny), max(x0, 0), min(x0 + S, nx)
    if yb > ya and xb > xa:
        out[ya - y0:yb - y0, xa - x0:xb - x0] = frame[ya:yb, xa:xb]
    return out


def _extract_cell(job: dict) -> dict:
    cid, S, out = job["cell_id"], job["S"], Path(job["out"])
    g = job["rows"]
    stacks = {name: _memmap(p) for name, p in job["images"].items()}
    labels = _memmap(job["labels"])
    d = out / "cells" / f"cell_{cid:04d}"
    d.mkdir(parents=True, exist_ok=True)
    crops = {name: np.zeros((len(g), S, S), np.uint16) for name in stacks}
    mask = np.zeros((len(g), S, S), np.uint8)
    y0s, x0s = [], []
    for i, r in enumerate(g.itertuples()):
        y0, x0 = int(round(r.cy_px)) - S // 2, int(round(r.cx_px)) - S // 2
        y0s.append(y0); x0s.append(x0)
        for name, st in stacks.items():
            crops[name][i] = np.clip(np.rint(_crop(np.asarray(st[r.t]), y0, x0, S, job["pad"][name])), 0, 65535)
        mask[i] = _crop(np.asarray(labels[r.t]), y0, x0, S, 0) == cid
    for name, arr in crops.items():
        tifffile.imwrite(d / f"{name}.tif", arr, photometric="minisblack")  # 3-4 frame stacks would otherwise be read as RGB
    tifffile.imwrite(d / "mask.tif", mask, photometric="minisblack")
    fr = g.copy(); fr["crop_y0"], fr["crop_x0"] = y0s, x0s
    fr.to_csv(d / "frames.csv", index=False)
    t = g.t.to_numpy()
    summary = dict(cell_id=cid, first_t=int(t.min()), last_t=int(t.max()), n_frames=len(t),
                   n_missing_frames=int(t.max() - t.min() + 1 - len(t)), border_frames=int(g.touches_border.sum()),
                   min_nn_dist_px=float(g.nn_dist_px.min()), mean_area_px=float(g.area_px.mean()))
    try:
        import psutil
        summary["_peak_mb"] = psutil.Process().memory_info().peak_wset / 2**20
    except Exception:
        summary["_peak_mb"] = float("nan")
    return summary


def isolate_cells(images: dict[str, str | Path], labels: str | Path, out_dir: str | Path, margin: float = 0.3,
                  n_workers: int | None = None, min_frames: int = 1) -> pd.DataFrame:
    """Write one single-cell movie per cell ID. Returns the per-cell summary.

    images: {"nucleus": path, "locus": path, ...}; labels: label stack path; margin: extra crop room as a fraction of
    the largest cell radius; n_workers: processes (default: CPUs-2, max 8); None/1 runs in-process.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    table = build_cell_table(labels)
    counts = table.groupby("cell_id").t.transform("count")
    table = table[counts >= min_frames].reset_index(drop=True)
    S = int(2 * np.ceil(table.rmax_px.max() * (1 + margin)) + 1)
    pad = {}
    for name, p in images.items():
        pad[name] = float(np.percentile(np.asarray(_memmap(p)[0]), 1))
    table["crop_size_px"] = S
    table.to_csv(out / "cell_table.csv", index=False)
    jobs = [dict(cell_id=int(cid), rows=g.sort_values("t"), S=S, out=str(out), images={k: str(v) for k, v in images.items()},
                 labels=str(labels), pad=pad) for cid, g in table.groupby("cell_id")]
    nw = n_workers if n_workers is not None else max(1, min(8, (os.cpu_count() or 2) - 2))
    if nw <= 1:
        res = [_extract_cell(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=nw) as ex:
            res = list(ex.map(_extract_cell, jobs, chunksize=max(1, len(jobs) // (4 * nw))))
    summ = pd.DataFrame(res).sort_values("cell_id").reset_index(drop=True)
    summ["crop_size_px"] = S
    summ.drop(columns="_peak_mb").to_csv(out / "cell_summary.csv", index=False)
    summ.attrs["peak_worker_mb"] = float(np.nanmax(summ["_peak_mb"])) if "_peak_mb" in summ else float("nan")
    return summ
