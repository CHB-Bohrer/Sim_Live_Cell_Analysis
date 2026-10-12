"""Compare locus-finding methods on a finished run (stage 6 only, on the saved single-cell crops with the true cell IDs).

    scripts\\run.cmd python scripts\\compare_localizers.py --run %SIMLIVE_DATA%\\runs\\<run folder>
    scripts\\run.cmd python scripts\\compare_localizers.py --run <folder> --variants original,paper --no-cal --workers 4

Variants (see src/simlive/stage6_analysis/compare.py): original placeholder, the Bennabi et al. 2025 steps (median filter, movie threshold,
area, Gaussian smoothing, local maxima) and one-change-at-a-time versions of it; 3D variants for z-stack runs. Writes
<run>/stage6_analysis/method_comparison.csv, shown in the dashboard's 🔎 Spot finding tab. One run is one seed: compare over several runs.
"""
import argparse
import json
from pathlib import Path

import pandas as pd

from simlive.io.runs import DATA_ROOT
from simlive.stage6_analysis.compare import compare
from simlive.stage6_analysis.bleach import bleach_correct_run
from simlive.stage6_analysis.ms2 import analyse_run_ms2


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="run folder, or a run name under $SIMLIVE_DATA/runs")
    ap.add_argument("--variants", default=None, help="comma-separated keys (default: all)")
    ap.add_argument("--workers", type=int, default=None, help="processes (default: CPU cores - 2, max 8)")
    ap.add_argument("--no-cal", action="store_true", help="skip the calibrated versions")
    ap.add_argument("--no-3d", action="store_true", help="skip the 3D versions of z-stack runs")
    ap.add_argument("--ms2-only", action="store_true", help="only (re)compute the MS2 spot traces (runs with transcription)")
    a = ap.parse_args()
    run = Path(a.run) if Path(a.run).exists() else DATA_ROOT / "runs" / a.run
    analyse_run_ms2(run, a.workers)                     # MS2 spot traces, if the run has the transcription model
    sj = run / "stage6_analysis" / "scores.json"         # photobleaching correction of the saved intensities (background removed first)
    scores = json.loads(sj.read_text()) if sj.exists() else {}
    scores["bleaching"] = bleach_correct_run(run)
    from simlive.stage6_analysis.deconvolve import deconvolve_run

    scores["ms2_deconvolution"] = deconvolve_run(run)          # MS2 trace -> Pol II initiation events (uses the bleaching-corrected trace)
    sj.write_text(json.dumps(scores, indent=2, default=str))
    print("bleaching:", {k: (round(v["tau_est_s"]) if isinstance(v, dict) and "tau_est_s" in v else v) for k, v in scores["bleaching"].items()})
    if a.ms2_only:
        return
    tab = compare(run, a.variants.split(",") if a.variants else None, a.workers, not a.no_cal, not a.no_3d)
    pd.set_option("display.width", 250)
    show = tab.drop(columns=["label", "settings"]).copy()
    for c in ("median_error_um", "p95_error_um", "rms_error_um", "rms_error_near_um"):
        show[c.replace("_um", "_nm")] = (1000 * show.pop(c)).round(1)
    print(show.round(4).to_string(index=False))
    print(f"\nsaved {run / 'stage6_analysis' / 'method_comparison.csv'}")


if __name__ == "__main__":
    main()
