"""Write a small FAKE run (diffusing cells + blinking spots) so you can try view_run.py before stage 3 exists.

    scripts\\run.cmd python scripts\\make_demo_run.py            # -> data\\runs\\demo
Replaced later by the real stage-3 forward model; the file layout below is the contract it will follow:
    data/runs/<run_id>/stage3_microscopy/{locus.tif, nucleus.tif, labels.tif}   (T,Y,X); labels = persistent cell ID
"""
import argparse
from pathlib import Path

import numpy as np
import tifffile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("data/runs/demo"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--frames", type=int, default=60)
    ap.add_argument("--cells", type=int, default=5)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    N, T, R = 256, a.frames, 22
    pos = rng.uniform(40, N - 40, (a.cells, 2))
    yy, xx = np.mgrid[:N, :N]
    loc = np.zeros((T, N, N), np.float32); nuc = np.zeros_like(loc); lab = np.zeros((T, N, N), np.uint16)
    for t in range(T):
        pos += rng.normal(0, 2.0, pos.shape)
        pos = np.clip(pos, R, N - R)
        for i, (y, x) in enumerate(pos, start=1):
            m = (yy - y) ** 2 + (xx - x) ** 2 < R**2
            lab[t][m] = i
            nuc[t][m] += 200
            for _ in range(2):  # two "loci" jittering inside the nucleus
                oy, ox = y + rng.normal(0, 6), x + rng.normal(0, 6)
                loc[t] += 800 * np.exp(-((yy - oy) ** 2 + (xx - ox) ** 2) / (2 * 1.5**2))
    nuc += rng.poisson(20, nuc.shape); loc += rng.poisson(10, loc.shape)
    d = a.out / "stage3_microscopy"; d.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(d / "locus.tif", loc); tifffile.imwrite(d / "nucleus.tif", nuc.astype(np.float32))
    tifffile.imwrite(d / "labels.tif", lab)
    print(f"Wrote demo run to {d}")


if __name__ == "__main__":
    main()
