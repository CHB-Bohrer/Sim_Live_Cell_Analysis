"""Compressed TIFFs (io.images): lossless, smaller, and readable frame by frame exactly like a memory map."""
import pickle

import numpy as np
import pytest
import tifffile

from simlive.io.images import PageStack, open_stack, write_tiff


@pytest.mark.parametrize("shape", [(6, 20, 30), (5, 3, 20, 30), (3, 20, 30)])   # 3 frames: must not become RGB
def test_compressed_stack_reads_back_bit_for_bit_and_frame_by_frame(tmp_path, shape):
    rng = np.random.default_rng(1)
    a = rng.poisson(100, shape).astype(np.uint16)
    write_tiff(tmp_path / "c.tif", a)
    write_tiff(tmp_path / "u.tif", a, compress=False)
    assert np.array_equal(tifffile.imread(tmp_path / "c.tif"), a)
    st, mm = open_stack(tmp_path / "c.tif"), open_stack(tmp_path / "u.tif")
    assert isinstance(st, PageStack) and isinstance(mm, np.memmap)       # old uncompressed runs keep memory-mapping
    assert st.shape == a.shape and st.dtype == a.dtype and len(st) == shape[0]
    for key in [0, -1, (1,), (2, slice(3, 9), slice(4, 12)), slice(1, 3), (slice(None, None, 2), 1)]:
        assert np.array_equal(np.asarray(st[key]), a[key]), key
    if len(shape) == 4:
        assert np.array_equal(st[3, 2], a[3, 2]) and np.array_equal(st[1, 0, 5:9], a[1, 0, 5:9])
    assert np.array_equal(np.asarray(st), a)
    st2 = pickle.loads(pickle.dumps(st))                                   # process pools pickle it
    assert np.array_equal(st2[2], a[2])
    st.close(); st2.close()


def test_label_masks_shrink_a_lot(tmp_path):
    lab = np.zeros((20, 128, 128), np.uint16)
    yy, xx = np.mgrid[:128, :128]
    for t in range(20):
        lab[t][(yy - 40 - t) ** 2 + (xx - 60) ** 2 < 400] = 3
    write_tiff(tmp_path / "c.tif", lab)
    assert (tmp_path / "c.tif").stat().st_size < lab.nbytes / 30
