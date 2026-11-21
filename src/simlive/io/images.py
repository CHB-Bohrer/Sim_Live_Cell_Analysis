"""Image files: every TIFF the pipeline writes goes through `write_tiff` (lossless zlib compression), and every stack that
used to be memory-mapped is opened with `open_stack`.

Why compressed: label masks shrink ~150x and noisy camera movies ~2.5x, bit-for-bit identical (measured 2026-11-21 on
loci_demo_seed11_24bba2: labels 307 -> 2 MB, nucleus 307 -> 140 MB, locus0 307 -> 110 MB). Each frame is its own
TIFF page, so reading one frame stays a few milliseconds. Compressed files cannot be memory-mapped, so `open_stack`
falls back to `PageStack`, which decodes only the frames that are indexed (old uncompressed runs still memory-map).
"""
from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import numpy as np
import tifffile

COMPRESSION = "zlib"          # readable everywhere (Fiji / ImageJ, napari, MATLAB, tifffile); zstd is only ~10% smaller
LEVEL = 6


def write_tiff(path: str | Path, arr: np.ndarray, compress: bool = True) -> None:
    """Write a (T, Y, X) / (T, Z, Y, X) stack as single-channel pages (never read back as RGB), losslessly compressed."""
    arr = np.asarray(arr)
    kw = {}
    if compress:
        kw = dict(compression=COMPRESSION, compressionargs={"level": LEVEL},
                  predictor=True if arr.dtype.kind in "ui" else None)
    tifffile.imwrite(path, arr, photometric="minisblack", **kw)


class PageStack:
    """Read-only array-like view of a TIFF stack that decodes pages on demand (for files that cannot be memory-mapped).

    Supports what the analysis uses: `shape`, `dtype`, `ndim`, `len()`, and indexing whose first index is an integer
    or a slice of frames (`st[t]`, `st[t, z]`, `st[t, y0:y1, x0:x1]`, `st[a:b]`). A few recently read frames are cached.
    """

    def __init__(self, path: str | Path, cache_frames: int = 4):
        self.path = str(path)
        with tifffile.TiffFile(self.path) as tif:
            s = tif.series[0]
            self.shape, self.dtype = tuple(s.shape), np.dtype(s.dtype)
            n_pages = len(s.pages)
        self.ndim = len(self.shape)
        self._per_frame = int(np.prod(self.shape[1:-2])) if self.ndim > 2 else 1   # pages per frame (z planes)
        if self.ndim > 2 and n_pages != self.shape[0] * self._per_frame:
            raise ValueError(f"{self.path}: {n_pages} pages do not match shape {self.shape}")
        self._tif = None
        self._cache: OrderedDict[int, np.ndarray] = OrderedDict()
        self._cache_n = cache_frames

    def __len__(self) -> int:
        return self.shape[0]

    def __getstate__(self):                       # picklable (process pools): reopen the file in the worker
        d = dict(self.__dict__); d["_tif"] = None; d["_cache"] = OrderedDict(); return d

    def _file(self) -> tifffile.TiffFile:
        if self._tif is None:
            self._tif = tifffile.TiffFile(self.path)
        return self._tif

    def _page(self, i: int) -> np.ndarray:
        return self._file().pages[i].asarray()

    def _frame(self, t: int) -> np.ndarray:
        if t < 0:
            t += self.shape[0]
        if not 0 <= t < self.shape[0]:
            raise IndexError(f"frame {t} out of range for {self.shape[0]} frames")
        if t in self._cache:
            self._cache.move_to_end(t); return self._cache[t]
        if self.ndim <= 2:
            f = self._file().asarray()[t]
        else:
            n = self._per_frame
            f = np.stack([self._page(t * n + k) for k in range(n)]).reshape(self.shape[1:]) if n > 1 else self._page(t)
        self._cache[t] = f
        if len(self._cache) > self._cache_n:
            self._cache.popitem(last=False)
        return f

    def __getitem__(self, key):
        key = key if isinstance(key, tuple) else (key,)
        first, rest = key[0], key[1:]
        if isinstance(first, (int, np.integer)):
            t = int(first)
            if self.ndim == 4 and rest and isinstance(rest[0], (int, np.integer)) and t in self._cache:
                return self._cache[t][rest]
            if self.ndim == 4 and rest and isinstance(rest[0], (int, np.integer)):   # one plane: decode one page only
                z = int(rest[0]) % self.shape[1]
                return self._page((t % self.shape[0]) * self._per_frame + z)[rest[1:]]
            return self._frame(t)[rest] if rest else self._frame(t)
        if isinstance(first, slice):
            ts = range(*first.indices(self.shape[0]))
            out = np.stack([self._frame(t) for t in ts]) if len(ts) else np.empty((0,) + self.shape[1:], self.dtype)
            return out[(slice(None),) + rest] if rest else out
        return np.asarray(self)[key]                 # anything else (fancy indexing): read everything

    def __array__(self, dtype=None, copy=None):
        a = tifffile.imread(self.path)
        return a.astype(dtype) if dtype is not None else a

    def close(self) -> None:
        if self._tif is not None:
            self._tif.close(); self._tif = None

    def __del__(self):
        try:
            self.close()
        except Exception:                            # interpreter shutdown: nothing useful to do
            pass


def open_stack(path: str | Path):
    """Memory-map an uncompressed TIFF stack, or return a `PageStack` for a compressed one. Read-only either way."""
    try:
        return tifffile.memmap(str(path), mode="r")
    except ValueError:                               # compressed / non-contiguous: tifffile cannot memory-map it
        return PageStack(path)
