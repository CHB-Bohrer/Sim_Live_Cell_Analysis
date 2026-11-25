"""List the manual chapters whose code changed on this branch (compared with origin/main, or a ref you give).

    python docs/manual/stale.py            # vs origin/main
    python docs/manual/stale.py HEAD~3     # vs any git ref

Only these chapters need a look; everything else stays as it is. A figure script is rerun only if the
validation result it shows has changed (not for every code edit)."""
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def changed_files(ref: str) -> list[str]:
    out = subprocess.run(["git", "diff", "--name-only", f"{ref}...HEAD"], cwd=ROOT, capture_output=True, text=True, check=True)
    work = subprocess.run(["git", "diff", "--name-only", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True)
    return sorted(set(out.stdout.split()) | set(work.stdout.split()))


def main():
    ref = sys.argv[1] if len(sys.argv) > 1 else "origin/main"
    files = changed_files(ref)
    chapters = yaml.safe_load((HERE / "manifest.yaml").read_text(encoding="utf-8"))["chapters"]
    hits = []
    for ch in chapters:
        touched = [f for f in files if any(f == c or f.startswith(c.rstrip("/") + "/") for c in ch["code"])]
        if touched:
            hits.append((ch, touched))
    if not hits:
        print("No chapter of the manual describes the changed code: nothing to update.")
    for ch, touched in hits:
        print(f"{ch['file']}: {len(touched)} changed file(s), e.g. {', '.join(touched[:3])}")
        if ch["figures"]:
            print(f"    figures (rerun only if the result they show changed): {', '.join(ch['figures'])}")


if __name__ == "__main__":
    main()
