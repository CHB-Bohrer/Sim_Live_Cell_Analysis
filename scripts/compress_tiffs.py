"""Losslessly compress every uncompressed TIFF under a folder (default: all runs in $SIMLIVE_DATA).

    scripts\\run.cmd python scripts\\compress_tiffs.py                 # every run
    scripts\\run.cmd python scripts\\compress_tiffs.py <folder> --dry-run

Each file is rewritten with simlive.io.images.write_tiff (zlib), read back and compared with the original array; only
if they are identical is the original replaced. New runs are written compressed already, so this is for old runs.
Do not run it while a pipeline run or the dashboard is reading the same runs.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import tifffile

from simlive.io.images import write_tiff
from simlive.io.runs import DATA_ROOT


def is_compressed(path: Path) -> bool:
    with tifffile.TiffFile(path) as tif:
        return tif.pages[0].compression != 1          # 1 = COMPRESSION.NONE


def compress_one(path: str) -> tuple[str, int, int, str]:
    p = Path(path)
    before = p.stat().st_size
    try:
        if is_compressed(p):
            return path, before, before, "already"
        with tifffile.TiffFile(p) as tif:
            series_shape = tif.series[0].shape
        a = tifffile.imread(p)
        if a.shape != series_shape:
            return path, before, before, f"skipped: shape {a.shape} vs series {series_shape}"
        tmp = p.with_name(p.name + ".ztmp")
        write_tiff(tmp, a)
        if not np.array_equal(tifffile.imread(tmp), a):
            tmp.unlink()
            return path, before, before, "FAILED: read-back differs"
        for attempt in range(5):                       # Dropbox / the indexer may hold the file for a moment
            try:
                os.replace(tmp, p); break
            except PermissionError:
                if attempt == 4:
                    raise
                time.sleep(2.0)
        return path, before, p.stat().st_size, "compressed"
    except (OSError, ValueError) as e:                   # unreadable / locked file: leave the original as it is
        Path(path + ".ztmp").unlink(missing_ok=True)
        return path, before, before, f"FAILED (original kept; run again): {e}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", nargs="?", default=str(DATA_ROOT / "runs"))
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    files = sorted(str(p) for p in Path(a.folder).rglob("*.tif"))
    print(f"{len(files)} TIFF files under {a.folder}", flush=True)
    if a.dry_run:
        todo = [f for f in files if not is_compressed(Path(f))]
        print(f"{len(todo)} uncompressed, {sum(Path(f).stat().st_size for f in todo) / 1e9:.1f} GB"); return
    t0, tot_b, tot_a, n_fail = time.time(), 0, 0, 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for i, (f, b, after, status) in enumerate(ex.map(compress_one, files, chunksize=4), 1):
            tot_b += b; tot_a += after
            if status.startswith(("FAILED", "skipped")):
                n_fail += 1; print(f"{status}: {f}", flush=True)
            if i % 200 == 0 or i == len(files):
                print(f"{i}/{len(files)} files, {tot_b / 1e9:.1f} GB -> {tot_a / 1e9:.1f} GB, {time.time() - t0:.0f} s", flush=True)
    print(f"done: {tot_b / 1e9:.1f} GB -> {tot_a / 1e9:.1f} GB; {n_fail} not compressed")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
