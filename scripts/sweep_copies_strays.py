"""Can the spot finder keep the promoter-enhancer copies apart and ignore stray spots? A sweep over the stray-spot rate.

    scripts\\run.cmd python scripts\\sweep_copies_strays.py
    scripts\\run.cmd python scripts\\sweep_copies_strays.py --rates 0,1,4 --seeds 51,52 --name quick
    scripts\\run.cmd python scripts\\sweep_copies_strays.py --name copies_strays_bright --rates 1 --set imaging_errors.stray_spots.brightness=1.5

Every movie: configs/loci_library.yaml (polymer-library loci; half the cells with two promoter-enhancer copies, copy 1 bursting twice as
often, probes and MS2 on; --config configs/loci_demo.yaml for the stand-in loci),
stray spots in both locus colours and the nuclear (MS2) colour at `rate` new spots per nucleus per minute (mean life 60 s, so `rate` strays
of each colour are present per nucleus on average), true cell IDs (tracking is not what is tested), no calibration, no deconvolution.
The SAME crops are then analysed three ways:
  one_path     the one-copy finder (Chen et al.: one Viterbi path per colour; no coverage test, no pairing, no pair gate),
  multi        the multi-copy linker without the pair gate (paths, coverage test, promoter-enhancer pairing),
  multi_gate   the multi-copy linker with the flat pair gate (1.5 um),
  multi_gauss  the multi-copy linker with the smooth (gauss) pair prior (the configs' default since 2026-10-19).
Scores per analysis (stage6_analysis.localize.score_localization, with the stray-spot truth): found, wrong picks (frac_far), picks of a
stray spot, picks of the other copy, copies found right per cell, missed / extra copies, promoter-enhancer distance error, and the MS2
trace read at each found promoter (copies matched to the true ones): share of labelled copies traced, ON frames found, trace r.
Writes $SIMLIVE_DATA/sweeps/<name>/results.csv (resumable); the dashboard's Spot finding page (part 3e) and Scans page show it.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from simlive.io.runs import DATA_ROOT, REPO_ROOT, load_config
from simlive.pipelines.tracking_demo import run_tracking_demo

VARIANTS = {"one_path": {"max_copies": 1, "pair_prior": "none"}, "multi": {"pair_prior": "none"},
            "multi_gate": {"pair_prior": "gate", "pair_gate_um": 1.5}, "multi_gauss": {"pair_prior": "gauss"}}
KEEP = ("detection_rate", "frac_far", "frac_on_stray", "frac_on_other_copy", "median_error_um", "copy_number_right", "missed_copies",
        "extra_copies", "true_copies")


def analyse(run: Path, over: dict, workers: int | None) -> dict:
    """Localize the loci and read the MS2 traces of every cell of a finished run again with the settings `over`, and score both."""
    from simlive.stage1_chromatin.copies import with_copy
    from simlive.stage5_linking.celldata import map_cells
    from simlive.stage6_analysis.localize import localize_loci, locus_distances, match_copies, photons_in_frame, score_localization
    from simlive.stage6_analysis.ms2 import ms2_trace, ms2_truth, relabel_to_true_copies, score_ms2
    from simlive.stage7_validation.propagation import true_distances

    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    px = cfg["optics"]["pixel_size_nm"] / 1000.0
    ana = {**(cfg.get("locus_analysis") or {}), **over}
    s3, out5 = run / "stage3_microscopy", run / "stage5_cells" / "truth"
    truth = photons_in_frame(pd.read_csv(s3 / "loci_truth_image.csv"), pd.read_csv(s3 / "imaging_truth.csv"))
    strays = pd.read_csv(s3 / "stray_spots.csv") if (s3 / "stray_spots.csv").exists() else None
    pos = map_cells(localize_loci, out5, n_workers=workers, px_um=px, n_loci=int(cfg["loci"]["n_loci"]), **ana)
    _, sc = score_localization(pos, truth, px, strays=strays)
    out = {k: sc.get(k) for k in KEEP}
    m = match_copies(pos, truth, px)
    d = locus_distances(m[m["copy"] < 50], px).merge(true_distances(with_copy(truth)), on=["cell_id", "t", "copy"])
    if len(d):
        e = 1000 * (d.dist_um - d.dist_true_um)
        out.update(dist_bias_nm=float(e.mean()), dist_rms_nm=float(np.sqrt((e ** 2).mean())), dist_p95_abs_nm=float(np.abs(e).quantile(0.95)))
    mt = ms2_truth(run)
    if mt is not None:
        kw = {k: v for k, v in ana.items() if k not in ("dilate_px", "method")}
        tr = relabel_to_true_copies(map_cells(ms2_trace, out5, n_workers=workers, px_um=px, **kw), mt, px)
        _, ms = score_ms2(tr[tr["copy"] < 50], mt, px)
        lab = mt[mt.ms2_labelled.astype(bool)] if "ms2_labelled" in mt else mt
        n_lab = lab.groupby(["cell_id", "copy"]).ngroups
        traced = tr[tr["copy"] < 50].groupby(["cell_id", "copy"]).ngroups
        out.update(ms2_copies_traced=traced / max(n_lab, 1), ms2_found_on=ms.get("found_on"), ms2_r=ms.get("intensity_r"))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rates", default="0,0.5,1,2,4", help="new stray spots per nucleus per minute (each colour)")
    ap.add_argument("--seeds", default="51,52,53,54,55")
    ap.add_argument("--name", default="copies_strays")
    ap.add_argument("--config", default=str(REPO_ROOT / "configs" / "loci_library.yaml"))
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL", help="extra overrides for every movie")
    ap.add_argument("--workers", type=int, default=None)
    a = ap.parse_args()
    out = DATA_ROOT / "sweeps" / a.name
    out.mkdir(parents=True, exist_ok=True)
    res_f = out / "results.csv"
    done = pd.read_csv(res_f) if res_f.exists() else pd.DataFrame()
    rows = done.to_dict("records")
    for rate in [float(v) for v in a.rates.split(",")]:
        for seed in [int(v) for v in a.seeds.split(",")]:
            if len(done) and ((done.stray_rate_per_min == rate) & (done.seed == seed)).any():
                continue
            sets = [f"seed={seed}", "tracking.skip=true", "loci.probes.enabled=true", "loci.transcription.enabled=true",
                    "loci.transcription.photons_per_loop_s=1000", "calibration.simulate=false", "locus_analysis.deconvolve=false",
                    f"imaging_errors.stray_spots.rate_per_min={rate}", f"imaging_errors.stray_spots.enabled={'true' if rate > 0 else 'false'}", *a.set]
            cfg = load_config(a.config, sets)
            t0 = time.time()
            run = run_tracking_demo(cfg, DATA_ROOT / "runs", progress=lambda s: None)
            base = {"stray_rate_per_min": rate, "seed": seed}
            for name, over in VARIANTS.items():
                rows.append({**base, "analysis": name, "run": run.name, **analyse(run, over, a.workers)})
            for r in rows[-len(VARIANTS):]:
                r["seconds"] = round(time.time() - t0)
            pd.DataFrame(rows).to_csv(res_f, index=False)
            g = {r["analysis"]: r for r in rows[-len(VARIANTS):]}
            print(f"rate {rate:g} seed {seed}: " + "; ".join(f"{k}: stray {v['frac_on_stray']:.3f} other copy {v['frac_on_other_copy']:.3f} "
                                                              f"copies right {v['copy_number_right']:.2f}" for k, v in g.items())
                  + f" ({time.time() - t0:.0f} s)", flush=True)
    print(f"Wrote {res_f}")


if __name__ == "__main__":
    main()
