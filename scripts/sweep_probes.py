"""Sweep the probe-binding model's parameters (fast: seconds, no movie, no GPU) and compare simulation with theory.

    scripts\\run.cmd python scripts\\sweep_probes.py --name koff --locus promoter --grid k_off=0.02,0.05,0.1,0.2
    scripts\\run.cmd python scripts\\sweep_probes.py --name two --locus enhancer --grid k_off=0.02,0.1 --grid n_probes=10,30,100 --set k_bind=0.1

Each --grid is PARAM=v1,v2,... (all combinations run); --set fixes a parameter for every run; everything else comes from
--config (default configs/loci_demo.yaml, section loci.probes). Results: <SIMLIVE_DATA>/probe_sweeps/<name>.csv, one row per
combination x replicate, with exact theory (th_*) next to simulated numbers (sim_*). The 🧪 Probes tab does the same interactively.

To sweep the probe parameters through the WHOLE movie pipeline (how binding noise changes locus-analysis accuracy), use
scripts\\sweep.py with keys like loci.probes.promoter.k_off (and --config configs\\loci_demo.yaml --set loci.probes.enabled=true).
"""
import argparse
from pathlib import Path

import yaml

from simlive.io.runs import DATA_ROOT, REPO_ROOT
from simlive.stage1_chromatin import probes as P


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--locus", choices=P.LOCUS_NAMES, default="promoter")
    ap.add_argument("--config", default=str(REPO_ROOT / "configs" / "loci_demo.yaml"))
    ap.add_argument("--grid", action="append", default=[], metavar="PARAM=v1,v2", required=True)
    ap.add_argument("--set", action="append", default=[], metavar="PARAM=VAL")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--duration", type=float, default=3000.0, help="seconds simulated per replicate")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    cfg = yaml.safe_load(Path(a.config).read_text())
    base = {**P.DEFAULTS[a.locus], **((cfg.get("loci", {}).get("probes") or {}).get(a.locus) or {})}
    for kv in a.set:
        k, _, v = kv.partition("=")
        base[k] = yaml.safe_load(v)
    grid = {}
    for g in a.grid:
        k, _, v = g.partition("=")
        if k not in P.PARAM_HELP:
            raise SystemExit(f"unknown parameter '{k}'; choose from {list(P.PARAM_HELP)}")
        grid[k] = [yaml.safe_load(x) for x in v.split(",")]
    df = P.sweep(base, grid, a.locus, a.reps, a.duration, a.seed)
    out = DATA_ROOT / "probe_sweeps" / f"{a.name}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    cols = list(grid) + ["th_mean_attached", "sim_mean_attached", "th_cv_attached", "sim_cv", "th_residence_s"]
    print(df.groupby(list(grid))[cols[len(grid):]].mean().round(3).to_string())
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
