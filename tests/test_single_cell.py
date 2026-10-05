"""CPU tests for stand-in loci, locus rendering, per-cell isolation, and the parallel per-cell runner."""
import numpy as np
import pandas as pd
import pytest
import tifffile

from simlive.io.runs import REPO_ROOT, load_config
from simlive.stage1_chromatin.toy_loci import simulate_toy_loci
from simlive.stage2_cells.motion import simulate_cells
from simlive.stage3_microscopy.render import render_loci, render_nuclei
from simlive.stage5_linking.celldata import load_cell, map_cells
from simlive.stage5_linking.isolate import isolate_cells
from simlive.stage6_analysis.localize import localize_loci, locus_distances, score_localization

CFG = REPO_ROOT / "configs" / "loci_demo.yaml"
SMALL = ["geometry.fov_um=[30, 30]", "geometry.cell_radius_um=3.0", "cells.n_cells=3", "acquisition.n_frames=8",
         "motion.speed_scale=0.3", "loci.photons_per_locus_s=8000"]


@pytest.fixture(scope="module")
def sim(tmp_path_factory):
    d = tmp_path_factory.mktemp("sim")
    cfg = load_config(CFG, SMALL)
    cells = simulate_cells(cfg, np.random.default_rng(3))
    img, lab = render_nuclei(cells, cfg, np.random.default_rng(4))
    truth = simulate_toy_loci(cells, cfg, np.random.default_rng(5))
    tifffile.imwrite(d / "nucleus.tif", np.rint(img).astype(np.uint16))
    tifffile.imwrite(d / "labels.tif", lab)
    for k in range(2):
        tifffile.imwrite(d / f"locus{k}.tif", np.rint(render_loci(truth, cfg, np.random.default_rng(10 + k), k)).astype(np.uint16))
    return cfg, cells, truth, lab, d


def test_loci_stay_inside_their_nucleus(sim):
    cfg, cells, truth, lab, d = sim
    px = cfg["optics"]["pixel_size_nm"] / 1000
    for r in truth.itertuples():
        y, x = int(round(r.y_um / px)), int(round(r.x_um / px))
        assert lab[r.t, y, x] == r.cell_id, "a locus fell outside its own nucleus"


def test_rendered_spot_is_at_the_true_position(sim):
    cfg, cells, truth, lab, d = sim
    px = cfg["optics"]["pixel_size_nm"] / 1000
    im = tifffile.imread(d / "locus0.tif", key=0).astype(float)
    t0 = truth[(truth.t == 0) & (truth.locus_id == 0)].iloc[0]
    cy, cx = t0.y_um / px, t0.x_um / px
    y0, x0 = int(round(cy)) - 4, int(round(cx)) - 4
    w = np.clip(im[y0:y0 + 9, x0:x0 + 9] - 104, 0, None)
    yy, xx = np.indices(w.shape)
    ey, ex = y0 + (w * yy).sum() / w.sum(), x0 + (w * xx).sum() / w.sum()
    assert abs(ey - cy) < 0.15 and abs(ex - cx) < 0.15, "spot centroid should match the true position (sub-pixel)"


def test_isolation_crops_are_centred_and_parallel_equals_serial(sim, tmp_path):
    cfg, cells, truth, lab, d = sim
    imgs = {"nucleus": d / "nucleus.tif", "locus0": d / "locus0.tif", "locus1": d / "locus1.tif"}
    s1 = isolate_cells(imgs, d / "labels.tif", tmp_path / "serial", n_workers=1)
    s2 = isolate_cells(imgs, d / "labels.tif", tmp_path / "parallel", n_workers=2)
    pd.testing.assert_frame_equal(s1.drop(columns="_peak_mb"), s2.drop(columns="_peak_mb"))
    S = int(s1.crop_size_px.iloc[0])
    assert S % 2 == 1
    for cid in s1.cell_id:
        a, b = load_cell(tmp_path / "serial", cid), load_cell(tmp_path / "parallel", cid)
        for k in a.channels:
            assert np.array_equal(a.channels[k], b.channels[k])
        # mask area equals the labeled area, and the mask is centred in the crop (centroid within 1 px of the middle)
        for i, t in enumerate(a.t):
            assert a.mask[i].sum() == (lab[t] == cid).sum()
            ys, xs = np.nonzero(a.mask[i])
            assert abs(ys.mean() - S // 2) < 1.0 and abs(xs.mean() - S // 2) < 1.0
        assert not a.mask[:, :2].any() and not a.mask[:, -2:].any(), "crop is big enough to contain the whole cell"


def test_per_cell_localization_matches_truth_and_parallel_equals_serial(sim, tmp_path):
    cfg, cells, truth, lab, d = sim
    imgs = {"nucleus": d / "nucleus.tif", "locus0": d / "locus0.tif", "locus1": d / "locus1.tif"}
    isolate_cells(imgs, d / "labels.tif", tmp_path, n_workers=1)
    px = cfg["optics"]["pixel_size_nm"] / 1000
    p1 = map_cells(localize_loci, tmp_path, n_workers=1, px_um=px, n_loci=2)
    p2 = map_cells(localize_loci, tmp_path, n_workers=2, px_um=px, n_loci=2)
    pd.testing.assert_frame_equal(p1.sort_values(["cell_id", "t", "locus_id"]).reset_index(drop=True),
                                  p2.sort_values(["cell_id", "t", "locus_id"]).reset_index(drop=True))
    err, sc = score_localization(p1, truth, px)
    assert sc["detection_rate"] > 0.95 and sc["rms_error_um"] < 0.05, sc
    dist = locus_distances(p1, px)
    assert len(dist) > 0 and (dist.dist_um > 0).all()


def test_isolation_handles_a_cell_that_leaves_the_image(tmp_path):
    lab = np.zeros((4, 60, 60), np.uint16)
    nuc = np.full((4, 60, 60), 100, np.uint16)
    for t in range(4):
        lab[t, 20:40, 2 - t if t < 3 else 0:20 - t if t < 3 else 18] = 1  # cell slides against the left border
    tifffile.imwrite(tmp_path / "lab.tif", lab, photometric="minisblack")
    tifffile.imwrite(tmp_path / "nuc.tif", nuc, photometric="minisblack")
    s = isolate_cells({"nucleus": tmp_path / "nuc.tif"}, tmp_path / "lab.tif", tmp_path / "out", n_workers=1)
    assert int(s.border_frames.iloc[0]) >= 1
    assert load_cell(tmp_path / "out", 1).channels["nucleus"].shape[1:] == (s.crop_size_px.iloc[0],) * 2


def test_three_frame_cell_is_stored_as_frames_not_colour_channels(tmp_path):
    lab = np.zeros((3, 40, 40), np.uint16)
    lab[:, 10:25, 10:25] = 1  # a cell that exists for exactly 3 frames (tifffile would treat 3 pages as RGB)
    nuc = np.full((3, 40, 40), 500, np.uint16)
    tifffile.imwrite(tmp_path / "lab.tif", lab, photometric="minisblack")
    tifffile.imwrite(tmp_path / "nuc.tif", nuc, photometric="minisblack")
    isolate_cells({"nucleus": tmp_path / "nuc.tif"}, tmp_path / "lab.tif", tmp_path / "out", n_workers=1)
    cell = load_cell(tmp_path / "out", 1)
    assert cell.channels["nucleus"].shape[0] == 3 and cell.mask.shape[0] == 3
    assert tifffile.imread(tmp_path / "out" / "cells" / "cell_0001" / "nucleus.tif", key=2).shape == cell.mask.shape[1:]
