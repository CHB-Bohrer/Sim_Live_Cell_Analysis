"""Stage 6: biplane registration from the beads that are in the image.

A biplane microscope puts the two focal planes of a pair on two halves of one camera. The halves are never registered perfectly, so a
point seen in both planes appears at slightly different (y, x) positions. A fiducial bead on the slide whose height lies between the two
planes of the lowest pair is bright in BOTH planes, so the difference of its two lateral positions measures the offset, in every frame,
from the same images that hold the cells (no separate calibration slide needed). The measured offset is then undone when the second plane
of each pair is cut out for the locus finder (`localize3d.localize_loci_3d(..., plane_shift_px=...)`).

What this does NOT measure: the true plane separation (a bead of unknown height gives one intensity ratio, not two unknowns), and the
z calibration of the PSF (a z-scan of a bead, as the manufacturer does). Both stay at their nominal values.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_dilation

from simlive.stage3_microscopy.calibration.fit import detect_beads


def _centroid(img: np.ndarray, y: float, x: float, half: int = 6) -> tuple[float, float, float] | None:
    """Background-subtracted centroid of the bead near (y, x) in one plane, and its peak above that background."""
    y0, x0 = int(round(y)), int(round(x))
    ya, xa = max(y0 - half, 0), max(x0 - half, 0)
    w = img[ya:y0 + half + 1, xa:x0 + half + 1].astype(float)
    if w.size < 9:
        return None
    w = np.clip(w - np.percentile(w, 20), 0, None)
    if w.sum() <= 0:
        return None
    yy, xx = np.indices(w.shape)
    return ya + float((w * yy).sum() / w.sum()), xa + float((w * xx).sum() / w.sum()), float(w.max())


def estimate_plane_offset(stack_path, seg_path=None, dilate_px: int = 8, pair: int = 0, min_beads: int = 2) -> tuple[tuple[float, float], pd.DataFrame]:
    """Offset (dy, dx) in pixels of the second plane of a biplane pair relative to the first, from the fiducial beads (default: the lowest pair).
    Returns ((dy, dx), table of every bead x frame measurement). (0, 0) if fewer than `min_beads` beads are usable."""
    st = tifffile.memmap(stack_path, mode="r")
    f0 = np.asarray(st[0])
    P0 = detect_beads(f0.max(axis=0))
    if seg_path is not None and len(P0):
        seg = tifffile.memmap(seg_path, mode="r")
        m = binary_dilation(np.asarray(seg[0]) > 0, iterations=dilate_px)
        P0 = P0[[not m[int(round(y)), int(round(x))] for y, x in P0]]
    rows = []
    ka, kb = 2 * pair, 2 * pair + 1
    for t in range(st.shape[0]):
        fr = np.asarray(st[t])
        for i, (y, x) in enumerate(P0):
            a, b = _centroid(fr[ka], y, x), _centroid(fr[kb], y, x)
            if a is not None and b is not None:      # peaks of THIS bead (until 2026-10-26 the whole plane's maximum, the same for every bead)
                rows.append((t, i, b[0] - a[0], b[1] - a[1], a[2], b[2]))
    tab = pd.DataFrame(rows, columns=["t", "bead", "dy_px", "dx_px", "peak_a", "peak_b"])
    if tab.empty or tab.bead.nunique() < min_beads:
        return (0.0, 0.0), tab
    good = tab[(tab[["peak_a", "peak_b"]].min(axis=1) > 0.25 * tab[["peak_a", "peak_b"]].max(axis=1).median())]   # the bead must be bright in both planes
    if good.bead.nunique() < min_beads:
        return (0.0, 0.0), tab
    return (float(good.dy_px.median()), float(good.dx_px.median())), tab
