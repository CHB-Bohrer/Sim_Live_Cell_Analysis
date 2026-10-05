"""Generate a LIBRARY of independent chromatin simulations (resumable; one simulation at a time on the GPU).

    scripts\\run.cmd python scripts\\run_chromatin_library.py --config configs\\chromatin\\chr21_loop_extrusion.yaml ^
        --library chr21_loop_extrusion --n 100

Each simulation runs in its own process (so GPU memory is released), seeds are first_seed .. first_seed+n-1, and
finished simulations are skipped, so you can stop (Ctrl-C / close the window) and run the same command later to
continue. Progress is printed and appended to <library>\\progress.log.
"""
import argparse
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

from simlive.io.runs import DATA_ROOT, REPO_ROOT, load_config


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--library", required=True)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--first-seed", type=int, default=1000)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VAL")
    ap.add_argument("--max-hours", type=float, default=None, help="stop starting new simulations after this long")
    a = ap.parse_args()

    cfg = load_config(a.config, a.set)
    lib = DATA_ROOT / "chromatin" / a.library
    lib.mkdir(parents=True, exist_ok=True)
    (lib / "library.json").write_text(json.dumps({"config": cfg, "n_requested": a.n, "first_seed": a.first_seed,
                                                  "created": datetime.now().isoformat(timespec="seconds")}, indent=2))
    log = (lib / "progress.log").open("a")

    def say(msg):
        line = f"{datetime.now():%Y-%m-%d %H:%M:%S}  {msg}"
        print(line, flush=True); log.write(line + "\n"); log.flush()

    t_start, done_times = time.time(), []
    seeds = list(range(a.first_seed, a.first_seed + a.n))
    todo = [s for s in seeds if not (lib / f"{cfg['name']}_seed{s}" / "meta.json").exists()]
    say(f"library '{a.library}': {len(seeds) - len(todo)} done, {len(todo)} to run")
    for k, s in enumerate(todo):
        if a.max_hours and (time.time() - t_start) / 3600 > a.max_hours:
            say(f"stopping: reached --max-hours {a.max_hours}")
            break
        t0 = time.time()
        cmd = [sys.executable, str(REPO_ROOT / "scripts" / "run_chromatin.py"), "--config", a.config, "--seed", str(s),
               "--library", a.library] + [x for kv in a.set for x in ("--set", kv)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            say(f"FAILED seed {s}: {r.stderr[-600:]}")
            continue
        done_times.append(time.time() - t0)
        left = len(todo) - k - 1
        say(f"seed {s} done in {done_times[-1] / 60:.1f} min  |  {left} left, about "
            f"{left * sum(done_times) / len(done_times) / 3600:.1f} h remaining")
    say("finished")


if __name__ == "__main__":
    main()
