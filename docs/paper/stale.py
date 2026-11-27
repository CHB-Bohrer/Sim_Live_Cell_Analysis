r"""Say which paper sections a code change made stale (so only those are opened).

    scripts\run.cmd python docs\paper\stale.py            # changes since the last documented commit
    scripts\run.cmd python docs\paper\stale.py --since REF

Each section file starts with a comment whose `Code:` line lists the paths it describes
(folders or files, comma separated) and whose `Check:` line names the checks to re-run.
The documented commit is in data/paper_state.json (written by `build.py --mark-done`).
"""
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    state = HERE / "data" / "paper_state.json"
    since = sys.argv[sys.argv.index("--since") + 1] if "--since" in sys.argv else (
        json.loads(state.read_text())["documented_commit"] if state.exists() else None)
    if not since:
        sys.exit("No documented commit yet: run build.py --mark-done once.")
    run = lambda *a: subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.split("\n")  # noqa: E731
    changed = sorted({f for f in run("diff", "--name-only", f"{since}..HEAD") + run("diff", "--name-only") if f
                      and not f.startswith("docs/paper/")})
    hits, covered = [], set()
    for sec in sorted((HERE / "sections").glob("*.md")):
        head = sec.read_text(encoding="utf-8").split("-->")[0]
        code = re.search(r"^Code:\s*(.+)$", head, re.M)
        check = re.search(r"^Check:\s*(.+)$", head, re.M)
        paths = [p.strip().rstrip("/") for p in code.group(1).split(",")] if code else []
        mine = [f for f in changed if any(f == p or f.startswith(p + "/") for p in paths if p)]
        if mine:
            covered.update(mine)
            hits.append((sec.name, mine, check.group(1).strip() if check else "none"))
    print(f"Changes since {since}:")
    for name, files, check in hits:
        shown = ", ".join(files[:6]) + (f" (+{len(files) - 6} more)" if len(files) > 6 else "")
        print(f"  edit sections/{name}   re-run checks: {check}\n      because of: {shown}")
    rest = [f for f in changed if f not in covered and f.startswith(("src/", "configs/"))]
    if rest:
        print("  NOT COVERED by any section (add to a Code: line): " + ", ".join(rest[:12]))
    if not hits and not rest:
        print("  nothing in the paper is stale.")


if __name__ == "__main__":
    main()
