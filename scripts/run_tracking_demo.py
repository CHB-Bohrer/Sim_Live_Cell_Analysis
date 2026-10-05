"""Cells-only tracking test: simulate moving cells, render nuclei, segment, track, score against ground truth.

    scripts\\run.cmd python scripts\\run_tracking_demo.py
    scripts\\run.cmd python scripts\\run_tracking_demo.py --set motion.D_um2_s=0.05 --set seed=2

Tracking is run twice: on GROUND-TRUTH masks (isolates linking error) and on AUTOMATIC-segmentation masks
(realistic). Both are scored against the true cell IDs. Prefer the dashboard (scripts\\dashboard.cmd) to see it.
Outputs in data/runs/<run_id>/ (one folder per stage; each has params.yaml, seed.txt, provenance.json).
"""
import argparse
import json
from pathlib import Path

from simlive.io.runs import DATA_ROOT, REPO_ROOT, load_config
from simlive.pipelines.tracking_demo import run_tracking_demo


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(REPO_ROOT / "configs" / "tracking_demo.yaml"))
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL", help="override, e.g. motion.D_um2_s=0.05")
    ap.add_argument("--out", type=Path, default=DATA_ROOT / "runs")
    a = ap.parse_args()

    run = run_tracking_demo(load_config(a.config, a.set), a.out)
    metrics = json.loads((run / "stage7_validation" / "metrics.json").read_text())
    print("\n=== Tracking accuracy vs ground truth ===")
    for variant, m in metrics.items():
        print(f"[{variant}]")
        for group, vals in m.items():
            print(f"  {group}: " + ", ".join(f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}"
                                              for k, v in vals.items()))
    print(f"\nSee it:  scripts\\dashboard.cmd")


if __name__ == "__main__":
    main()
