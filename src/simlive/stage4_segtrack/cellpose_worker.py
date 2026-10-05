"""Run Cellpose on a TIFF stack in a clean process:  python -m simlive.stage4_segtrack.cellpose_worker IN.tif OUT.tif '{json params}'

Why a separate process: on Windows, importing Cellpose after PyTorch makes cellpose.plot pull in scikit-image, whose
compiled OpenMP extension loads a second OpenMP runtime next to PyTorch's (libomp.dll vs libiomp5md.dll) and aborts the
process ("OMP: Error #15"). Cellpose only uses scikit-image for drawing, so we block that import; Cellpose already
handles its absence. We deliberately do NOT use KMP_DUPLICATE_LIB_OK (can silently give wrong results).
"""
import json
import sys

sys.modules["skimage"] = None  # makes `from skimage import ...` raise ImportError, which cellpose.plot tolerates

import numpy as np  # noqa: E402
import tifffile  # noqa: E402


def main(in_path: str, out_path: str, params_json: str) -> None:
    from cellpose import models

    p = json.loads(params_json)
    images = tifffile.imread(in_path)
    model = models.CellposeModel(gpu=p.get("gpu", True), pretrained_model=p.get("model", "cpsam"))
    masks, _flows, _styles = model.eval(
        list(images), diameter=p.get("diameter_px"), flow_threshold=p.get("flow_threshold", 0.4),
        cellprob_threshold=p.get("cellprob_threshold", 0.0), batch_size=p.get("batch_size", 8))
    tifffile.imwrite(out_path, np.stack(masks).astype(np.uint16))


if __name__ == "__main__":
    main(*sys.argv[1:4])
