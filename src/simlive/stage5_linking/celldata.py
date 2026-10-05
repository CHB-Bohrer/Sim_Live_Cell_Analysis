"""Load isolated single-cell movies (written by isolate.py) and run a function on every cell in parallel."""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_dilation


@dataclass
class CellData:
    cell_id: int
    frames: pd.DataFrame           # one row per frame the cell exists in (t, cy_px, cx_px, theta_rad, crop_y0, ...)
    mask: np.ndarray               # (n, S, S) bool: this cell's pixels only
    channels: dict                 # name -> (n, S, S) uint16, e.g. "nucleus", "locus0", "locus1" (one per colour)
    crop_size: int

    @property
    def t(self) -> np.ndarray:
        return self.frames.t.to_numpy()

    @property
    def nucleus(self):
        return self.channels.get("nucleus")

    @property
    def locus_channels(self) -> list[str]:
        return sorted(k for k in self.channels if k.startswith("locus"))

    def analysis_mask(self, dilate_px: int = 2) -> np.ndarray:
        """This cell's pixels grown by dilate_px (to include the PSF halo of loci near the nuclear edge)."""
        if dilate_px <= 0:
            return self.mask.copy()
        st = np.ones((1, 2 * dilate_px + 1, 2 * dilate_px + 1), bool)
        return binary_dilation(self.mask, structure=st)


def list_cells(source_dir: str | Path) -> pd.DataFrame:
    return pd.read_csv(Path(source_dir) / "cell_summary.csv")


def load_cell(source_dir: str | Path, cell_id: int) -> CellData:
    d = Path(source_dir) / "cells" / f"cell_{int(cell_id):04d}"
    fr = pd.read_csv(d / "frames.csv")
    mask = tifffile.imread(d / "mask.tif").astype(bool)
    channels = {p.stem: tifffile.imread(p) for p in sorted(d.glob("*.tif")) if p.stem != "mask"}
    return CellData(int(cell_id), fr, mask, channels, mask.shape[-1])


def _run_one(args):
    func, source_dir, cell_id, kw = args
    return func(load_cell(source_dir, cell_id), **kw)


def map_cells(func: Callable[..., pd.DataFrame], source_dir: str | Path, cell_ids=None, n_workers: int | None = None,
              **kw) -> pd.DataFrame:
    """Run func(cell, **kw) -> DataFrame for every cell, one process per core, and concatenate (cell_id column kept).

    `func` must be a top-level (picklable) function. Each worker loads only its own cell, so memory stays small.
    """
    source_dir = Path(source_dir)
    ids = list(cell_ids) if cell_ids is not None else list_cells(source_dir).cell_id.tolist()
    jobs = [(func, str(source_dir), int(i), kw) for i in ids]
    nw = n_workers if n_workers is not None else max(1, min(8, (os.cpu_count() or 2) - 2))
    if nw <= 1 or len(jobs) <= 1:
        parts = [_run_one(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=nw) as ex:
            parts = list(ex.map(_run_one, jobs, chunksize=max(1, len(jobs) // (4 * nw))))
    parts = [p for p in parts if p is not None and len(p)]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
