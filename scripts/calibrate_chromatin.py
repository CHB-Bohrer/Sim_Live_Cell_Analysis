"""Make the polymer -> nm / seconds calibration from the saved validation simulations and save it as configs/calibration/<name>.json.

    scripts\\run.cmd python scripts\\calibrate_chromatin.py
    scripts\\run.cmd python scripts\\calibrate_chromatin.py --D 0.004 --lag 60 --name fast_locus
    scripts\\run.cmd python scripts\\calibrate_chromatin.py --nm-per-unit 60        # a measured length scale instead of the contact-fraction rule

Library movies read it through `loci.calibration: default` (see src/simlive/calibration/polymer.py). The assumptions it rests on (a 100 kb pair is within
400 nm 10% of the time; a single locus has D = 0.002 um^2/s at a 1 min lag) are the user's of 2026-10-14 and are stored in the file; re-run with measured numbers.
"""
import argparse

from simlive.calibration.polymer import CALIBRATION_DIR, SIMS, VAL_DIR, calibrate, calibration_path, load_pool


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="default")
    ap.add_argument("--sep-kb", type=float, default=100.0, help="genomic separation of the promoter-enhancer pair used for the contact fraction")
    ap.add_argument("--working-nm", type=float, default=400.0)
    ap.add_argument("--contact-fraction", type=float, default=0.10)
    ap.add_argument("--D", type=float, default=0.002, help="single-locus D (um^2/s) at --lag")
    ap.add_argument("--lag", type=float, default=60.0, help="lag (s) at which D applies")
    ap.add_argument("--nm-per-unit", type=float, default=None, help="measured length scale (skips the contact-fraction rule)")
    ap.add_argument("--block-s", type=float, default=None, help="measured block duration (skips the MSD match)")
    a = ap.parse_args()

    missing = [s for s in SIMS if not (VAL_DIR / s / "tracked_positions.npy").exists()]
    if missing:
        raise SystemExit(f"validation simulations missing in {VAL_DIR}: {missing}\nRun scripts\\validate_chromatin.py run first.")
    pool = load_pool(a.sep_kb)
    cal = calibrate(pool, a.working_nm, a.contact_fraction, a.D, a.lag, nm_per_unit=a.nm_per_unit, block_s=a.block_s)
    out = cal.save(calibration_path(a.name))
    print(f"nm per polymer unit : {cal.nm_per_unit:.2f}")
    print(f"seconds per block   : {cal.block_s:.3f}   (MSD exponent {cal.alpha:.2f} below one block)")
    print(f"the 2000-block saved trajectory covers {2000 * cal.block_s / 60:.1f} min")
    print(f"saved {out.relative_to(CALIBRATION_DIR.parents[1])}")


if __name__ == "__main__":
    main()
