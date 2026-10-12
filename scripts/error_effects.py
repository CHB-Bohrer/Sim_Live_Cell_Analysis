"""What each source of image error does to the results, and how much each correction gives back (stage7_validation/error_effects.py).

    scripts\\run.cmd python scripts\\error_effects.py --run <run folder or name>
    scripts\\run.cmd python scripts\\error_effects.py --run <run> --sources hot_pixels,chromatic_shift,reference

Re-renders the run's movie for the reference microscope, for the reference + each error source the run has, and for all sources together;
analyses each with the correction off and on (calibration images / the movie alone / the analysis' own countermeasures) and scores it
against the truth. About 1-2 minutes per error source. Resumable: rows already in <run>/stage7_validation/error_effects.csv are skipped
(delete the file to start again). The dashboard's 🖼 Images page shows the result per source.
"""
import argparse
from pathlib import Path

from simlive.io.runs import DATA_ROOT
from simlive.stage7_validation.error_effects import run_error_effects


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--sources", default=None, help="comma-separated source ids (and/or reference, all); default: every source the run has")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--keep-work", action="store_true", help="keep the re-rendered movies (several GB)")
    a = ap.parse_args()
    run = Path(a.run) if Path(a.run).exists() else DATA_ROOT / "runs" / a.run
    run_error_effects(run, a.sources.split(",") if a.sources else None, a.workers, keep_work=a.keep_work,
                      progress=lambda s: print(s, flush=True))
    print("finished", flush=True)


if __name__ == "__main__":
    main()
