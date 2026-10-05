"""Stage 4a: automatic segmentation. Registry of segmenters; each maps (T,Y,X) image -> (T,Y,X) instance labels.

`threshold_watershed` is a classical baseline needing no model download. `cellpose` downloads its model weights
on first use. StarDist can be added to SEGMENTERS with the same signature.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label
from skimage.morphology import remove_small_objects
from skimage.segmentation import watershed


def threshold_watershed(images: np.ndarray, smooth_px: float = 2.0, min_area_px: int = 30,
                        min_distance_px: int = 8, **_) -> np.ndarray:
    out = np.zeros(images.shape, np.uint16)
    for t, im in enumerate(images):
        sm = gaussian(im, smooth_px)
        mask = remove_small_objects(sm > threshold_otsu(sm), min_area_px)
        dist = ndi.distance_transform_edt(mask)
        peaks = peak_local_max(gaussian(dist, 1.0), min_distance=min_distance_px, labels=label(mask),
                               exclude_border=False)
        markers = np.zeros(mask.shape, int)
        markers[tuple(peaks.T)] = np.arange(1, len(peaks) + 1)
        out[t] = watershed(-dist, markers, mask=mask)
    return out


def cellpose_seg(images: np.ndarray, diameter_px: float | None = None, flow_threshold: float = 0.4,
                 cellprob_threshold: float = 0.0, model: str = "cpsam", gpu: bool = True, batch_size: int = 8,
                 **_) -> np.ndarray:
    """Cellpose (v4 default model 'cpsam') on every frame; labels are unique per frame.

    Runs in a separate process (see cellpose_worker.py for why) and exchanges data through temporary TIFF files.
    """
    import json
    import subprocess
    import sys
    import tempfile
    from pathlib import Path

    import tifffile

    params = dict(diameter_px=diameter_px, flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold,
                  model=model, gpu=gpu, batch_size=batch_size)
    with tempfile.TemporaryDirectory() as d:
        src, dst = Path(d) / "in.tif", Path(d) / "out.tif"
        tifffile.imwrite(src, images.astype(np.float32))
        proc = subprocess.run([sys.executable, "-m", "simlive.stage4_segtrack.cellpose_worker", str(src), str(dst),
                               json.dumps(params)], capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"Cellpose worker failed (exit {proc.returncode}):\n{proc.stdout[-1500:]}\n{proc.stderr[-1500:]}")
        return tifffile.imread(dst)


SEGMENTERS = {"threshold_watershed": threshold_watershed, "cellpose": cellpose_seg}


def segment(images: np.ndarray, method: str, **params) -> np.ndarray:
    if method not in SEGMENTERS:
        raise ValueError(f"Unknown segmenter {method!r}; available: {sorted(SEGMENTERS)}")
    return SEGMENTERS[method](images, **params)
