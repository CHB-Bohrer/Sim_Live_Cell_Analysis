"""Run the tracking pipeline over a grid of settings x seeds and summarize accuracy.

    scripts\\run.cmd python scripts\\sweep.py --name texture_vs_seg --seeds 1-5 \\
        --grid nucleus_texture.contrast=0,0.3 --grid segmentation.method=threshold_watershed,cellpose

Each --grid is KEY=v1,v2,... (all combinations are run, each with every seed). Results go to
data/sweeps/<name>/results.csv (one row per run x mask variant) and a mean +- sd summary is printed.
Individual runs are kept in data/runs/ and can be opened in the dashboard.
"""
import argparse
import itertools
import json
import logging
import time
import traceback
from pathlib import Path

import pandas as pd
import yaml

from simlive.io.runs import DATA_ROOT, REPO_ROOT, load_config
from simlive.pipelines.tracking_demo import run_tracking_demo

METRICS = {"CHOTA": ("CHOTAMetric", "CHOTA"), "LNK": ("CTCMetrics", "LNK"), "TRA": ("CTCMetrics", "TRA"),
           "DET": ("CTCMetrics", "DET"), "target_eff": ("TrackOverlapMetrics", "target_effectiveness"),
           "id_switches": ("Identity", "id_switches"), "identity_kept": ("Identity", "identity_preserved_fraction")}


def parse_seeds(s: str) -> list[int]:
    if "-" in s:
        a, b = s.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(",")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--config", default=str(REPO_ROOT / "configs" / "tracking_demo.yaml"))
    ap.add_argument("--grid", action="append", default=[], metavar="KEY=v1,v2")
    ap.add_argument("--seeds", default="1-5")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL", help="fixed override for every run")
    a = ap.parse_args()
    logging.getLogger("traccuracy").setLevel(logging.CRITICAL)

    axes = {}
    for g in a.grid:
        k, _, v = g.partition("=")
        axes[k] = v.split(",")
    out_dir = DATA_ROOT / "sweeps" / a.name
    out_dir.mkdir(parents=True, exist_ok=True)
    combos = list(itertools.product(*axes.values())) if axes else [()]
    seeds = parse_seeds(a.seeds)
    total, rows, n = len(combos) * len(seeds), [], 0
    for combo in combos:
        for seed in seeds:
            n += 1
            ov = a.set + [f"{k}={v}" for k, v in zip(axes, combo)] + [f"seed={seed}"]
            t0 = time.time()
            try:
                run = run_tracking_demo(load_config(a.config, ov), DATA_ROOT / "runs", progress=lambda s: None)
                m = json.loads((run / "stage7_validation" / "metrics.json").read_text())
                for variant in ("gt_masks", "auto_masks"):
                    row = {**dict(zip(axes, combo)), "seed": seed, "masks": variant, "run": run.name}
                    row.update({k: m[variant][g][x] for k, (g, x) in METRICS.items()})
                    rows.append(row)
                print(f"[{n}/{total}] {dict(zip(axes, combo))} seed={seed}  {time.time() - t0:.0f}s", flush=True)
            except Exception:
                print(f"[{n}/{total}] FAILED {dict(zip(axes, combo))} seed={seed}\n{traceback.format_exc()[-800:]}", flush=True)
            pd.DataFrame(rows).to_csv(out_dir / "results.csv", index=False)
    df = pd.DataFrame(rows)
    (out_dir / "sweep.yaml").write_text(yaml.safe_dump({"grid": axes, "seeds": seeds, "fixed": a.set}))
    if df.empty:
        raise SystemExit("No successful runs.")
    keys = list(axes) + ["masks"]
    summ = df.groupby(keys)[list(METRICS)].agg(["mean", "std"]).round(3)
    print("\n=== mean (sd) over", len(seeds), "seeds ===")
    for k, g in df.groupby(keys):
        parts = [f"{m}={g[m].mean():.3f}±{g[m].std():.3f}" for m in ("CHOTA", "LNK", "target_eff", "id_switches")]
        print(f"{dict(zip(keys, k if isinstance(k, tuple) else (k,)))}  " + "  ".join(parts))
    summ.to_csv(out_dir / "summary.csv")
    print(f"\nSaved {out_dir / 'results.csv'} and summary.csv")


if __name__ == "__main__":
    main()
