"""How bright must the MS2 spot be, relative to the nuclear MCP/NLS-GFP, for transcription to be measurable?

    scripts\\run.cmd python scripts\\sweep_ms2_brightness.py
    scripts\\run.cmd python scripts\\sweep_ms2_brightness.py --ratios 0.8,8 --seeds 41,42 --name quick

The ratio R = photons of ONE Pol II carrying all its MS2 loops / photons of one pixel of nuclear GFP (both per exposure, before bleaching)
= n_loops x photons_per_loop_s / optics.photons_per_px_s. The placeholder (100 photons per loop per s, 3000 per pixel per s, 24 loops) is
R = 0.8. A spot's brightest pixel holds ~25-30% of its photons (PSF sigma ~0.7 px), so its peak contrast over the nucleus is ~0.25 x R per
Pol II. To measure R on a real movie: (sum of the background-subtracted spot photons) / (mean nuclear pixel value minus camera offset), divided
by the number of Pol II on the gene (or report it per spot and compare with R x Pol II per spot).

Every movie: the toy loci config (configs/loci_demo.yaml) with probes and transcription on, tracking skipped (true cell IDs; tracking is not
what is tested here), calibration applied. Per movie: the MS2 trace read at the tracked promoter (default) and at the MS2 spot's own track.
Writes $SIMLIVE_DATA/sweeps/<name>/results.csv; the dashboard's Spot finding page shows it.
"""
import argparse
import json
import time
from pathlib import Path

import pandas as pd

from simlive.io.runs import DATA_ROOT, REPO_ROOT, load_config
from simlive.pipelines.tracking_demo import run_tracking_demo

KEEP = ("found_on", "detected_on", "false_on_off", "median_error_um", "frac_far", "intensity_r", "call_sensitivity", "call_specificity")


def own_track_scores(run: Path, workers: int | None) -> dict:
    """The same trace with the MS2 spot found in the nuclear channel itself (ms2_position: own), calibrated."""
    import yaml

    from simlive.stage5_linking.celldata import map_cells
    from simlive.stage6_analysis.ms2 import ms2_trace, ms2_truth, score_ms2

    cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    ana = {k: v for k, v in (cfg.get("locus_analysis") or {}).items() if k not in ("dilate_px", "method")}
    cal = run / "stage3_microscopy" / "calibration"
    tr = map_cells(ms2_trace, run / "stage5_cells" / "truth", n_workers=workers, ms2_position="own",
                   calib=str(cal) if (cal / "calibration.json").exists() else None, **ana)
    tr.to_csv(run / "stage6_analysis" / "truth" / "ms2_traces_own_cal.csv", index=False)
    return score_ms2(tr, ms2_truth(run), cfg["optics"]["pixel_size_nm"] / 1000.0)[1]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ratios", default="0.25,0.8,2.5,8,25,80", help="R = one Pol II's MS2 photons / one nuclear pixel")
    ap.add_argument("--seeds", default="41,42,43")
    ap.add_argument("--name", default="ms2_brightness")
    ap.add_argument("--config", default=str(REPO_ROOT / "configs" / "loci_demo.yaml"))
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL", help="extra overrides for every movie")
    ap.add_argument("--workers", type=int, default=None)
    a = ap.parse_args()
    out = DATA_ROOT / "sweeps" / a.name
    out.mkdir(parents=True, exist_ok=True)
    res_f = out / "results.csv"
    done = pd.read_csv(res_f) if res_f.exists() else pd.DataFrame()
    rows = done.to_dict("records")
    for R in [float(v) for v in a.ratios.split(",")]:
        for seed in [int(v) for v in a.seeds.split(",")]:
            if len(done) and ((done.ratio == R) & (done.seed == seed)).any():
                continue
            base = load_config(a.config, [])
            n_loops = int(base["loci"]["transcription"]["params"]["n_loops"])
            ppl = R * float(base["optics"]["photons_per_px_s"]) / n_loops
            cfg = load_config(a.config, [f"seed={seed}", "tracking.skip=true", "loci.probes.enabled=true", "loci.transcription.enabled=true",
                                         f"loci.transcription.photons_per_loop_s={ppl}", *a.set])
            t0 = time.time()
            run = run_tracking_demo(cfg, DATA_ROOT / "runs", progress=lambda s: None)
            sc = json.loads((run / "stage6_analysis" / "scores.json").read_text())
            row = {"ratio": R, "seed": seed, "photons_per_loop_s": ppl, "run": run.name, "seconds": round(time.time() - t0)}
            for key, tag in (("ms2_calibrated", "promoter"), ("ms2", "promoter_raw")):
                for k in KEEP:
                    row[f"{tag}_{k}"] = (sc.get(key) or {}).get(k)
            own = own_track_scores(run, a.workers)
            for k in KEEP:
                row[f"own_{k}"] = own.get(k)
            tx = pd.read_csv(run / "stage1_chromatin" / "transcription.csv")
            row["frac_frames_on"] = float((tx.ms2_loops > 0).mean())
            row["median_loops_when_on"] = float(tx.ms2_loops[tx.ms2_loops > 0].median()) if (tx.ms2_loops > 0).any() else 0.0
            rows.append(row)
            pd.DataFrame(rows).to_csv(res_f, index=False)
            print(f"R={R:g} seed={seed}: at the promoter r={row['promoter_intensity_r']:.2f}, own track found {row['own_found_on'] or 0:.2f} "
                  f"r={row['own_intensity_r'] if row['own_intensity_r'] is not None else float('nan'):.2f} ({row['seconds']} s)", flush=True)
    print(f"Wrote {res_f}")


if __name__ == "__main__":
    main()
