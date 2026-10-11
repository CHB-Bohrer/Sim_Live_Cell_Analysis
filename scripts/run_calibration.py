"""Simulate the calibration of an existing run's instrument (dark frames, flat-field slide, bead slide), estimate the instrument from those
images, compare with the truth, and save everything next to the run's images.

    scripts\\run.cmd python scripts\\run_calibration.py --run %SIMLIVE_DATA%\\runs\\<run folder>
    scripts\\run.cmd python scripts\\run_calibration.py --config configs\\loci_demo.yaml --seed 1      (no run needed; just prints the table)

Settings (number of frames, brightness levels, number of beads...) come from the `calibration:` block of the run's config; see
src/simlive/stage3_microscopy/calibration/run.py for the defaults.
"""
import argparse
from pathlib import Path

import pandas as pd
import yaml

from simlive.io.runs import REPO_ROOT, load_config
from simlive.stage3_microscopy.calibration.run import run_calibration
from simlive.stage3_microscopy.render import entropy_for


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", help="a run folder (uses its resolved config and seed)")
    ap.add_argument("--config", help="a config file (instead of --run)")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--set", action="append", default=[], help="override, e.g. --set imaging_errors.flat_field.vignette=0.3")
    a = ap.parse_args()
    if a.run:
        run = Path(a.run)
        cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
    elif a.config:
        run = None
        cfg = load_config(Path(a.config), a.set)
    else:
        ap.error("give --run or --config")
    seed = int(a.seed if a.seed is not None else cfg["seed"])
    cal, table = run_calibration(cfg, entropy_for(seed))
    pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 60)
    print(table[["quantity", "channel", "unit", "estimate", "truth", "error_pct"]].round(4).to_string())
    if run is not None:
        d = run / "stage3_microscopy" / "calibration"
        cal.save(d)
        table.to_csv(d / "comparison.csv", index=False)
        print("saved to", d)


if __name__ == "__main__":
    main()
