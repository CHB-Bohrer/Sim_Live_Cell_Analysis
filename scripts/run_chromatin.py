"""Run ONE chromatin simulation from a config (stage 1, standalone).

    scripts\\run.cmd python scripts\\run_chromatin.py --config configs\\chromatin\\smoke_test.yaml --seed 1
    scripts\\run.cmd python scripts\\run_chromatin.py --config ... --seed 1 --set dynamics.production_blocks=200

Output: %SIMLIVE_DATA%\\chromatin\\<library>\\<name>_seed<seed>\\ (see simlive/stage1_chromatin/polymer3d.py).
"""
import argparse
from pathlib import Path

from simlive.io.runs import DATA_ROOT, REPO_ROOT, load_config
from simlive.stage1_chromatin.polymer3d import run_chromatin_simulation


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL")
    ap.add_argument("--library", default="scratch", help="subfolder of chromatin/ to write into")
    ap.add_argument("--out", type=Path, help="exact output folder (overrides --library)")
    a = ap.parse_args()
    cfg = load_config(a.config, a.set)
    out = a.out or DATA_ROOT / "chromatin" / a.library / f"{cfg['name']}_seed{a.seed}"
    print(f"Output: {out}")
    meta = run_chromatin_simulation(cfg, out, a.seed)
    print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in meta.items()})


if __name__ == "__main__":
    main()
