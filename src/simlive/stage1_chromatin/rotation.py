"""Uniformly random 3D rotation, with plain arithmetic only (shared by `library_loci` and the E-P polymer study)."""
from __future__ import annotations

import numpy as np


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    """Uniformly random 3x3 rotation matrix (from a normalised Gaussian quaternion). Not scipy's Rotation: its `.apply` loads a
    second OpenMP runtime (libiomp5md + libomp), after which `import torch` in the tracking step fails with WinError 127 on shm.dll
    (every library-source movie crashed at tracking until this was replaced). The same holds for ANY numpy matmul (`@`, dot):
    numpy here is MKL-linked and the first BLAS call loads those runtimes, so the code that runs before the tracker must use
    elementwise arithmetic (the loci are applied with a broadcast sum for that reason)."""
    q = rng.normal(size=4)
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
