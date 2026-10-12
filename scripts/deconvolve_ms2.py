"""MS2 trace -> Pol II initiation events and bursts (stage6_analysis/deconvolve.py, after Chen, Levo, Zoller & Gregor 2025), on a finished run.

    scripts\\run.cmd python scripts\\deconvolve_ms2.py --run <run folder or name>

Needs the run's MS2 traces (made by the pipeline, or by scripts/compare_localizers.py --ms2-only for older runs). Writes
stage6_analysis/truth/ms2_deconvolution.csv (+ _samples.npz, _info.json with the scores against the simulation). The dashboard's
🧮 MS2 → Pol II page shows it.
"""
import argparse
import json
from pathlib import Path

from simlive.io.runs import DATA_ROOT
from simlive.stage6_analysis.deconvolve import deconvolve_run


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--samples", type=int, default=1000, help="posterior samples kept per cell (the paper: 1000)")
    a = ap.parse_args()
    run = Path(a.run) if Path(a.run).exists() else DATA_ROOT / "runs" / a.run
    sc = deconvolve_run(run, ("truth",), progress=lambda s: print(s, flush=True), n_keep=a.samples)
    sj = run / "stage6_analysis" / "scores.json"
    if sj.exists():
        scores = json.loads(sj.read_text())
        scores["ms2_deconvolution"] = sc
        sj.write_text(json.dumps(scores, indent=2, default=str))
    print(json.dumps(sc, indent=1))


if __name__ == "__main__":
    main()
