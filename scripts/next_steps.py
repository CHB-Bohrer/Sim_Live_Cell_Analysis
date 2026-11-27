#!/usr/bin/env python3
"""List recorded next steps.

Sources: STATUS.md lines like "- [auto] ..." / "- [ask] ...", and, if the project
has a design document, \\nextstep[auto|ask]{...} in docs/sections/*.tex.

  scripts\run.cmd python scripts\next_steps.py          # everything
  scripts\run.cmd python scripts\next_steps.py --auto   # only those the assistant may just do
"""
import pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
TEX = re.compile(r"\\nextstep(?:\[(\w+)\])?\{((?:[^{}]|\{[^{}]*\})*)\}", re.S)
MD = re.compile(r"^\s*-\s*\[(auto|ask)\]\s+(.+)$", re.M)
only_auto = "--auto" in sys.argv
found = []

status = ROOT / "STATUS.md"
if status.exists():
    found += [(k, "STATUS.md", b) for k, b in MD.findall(status.read_text()) if not b.startswith("(")]
for f in sorted((ROOT / "docs" / "sections").glob("*.tex")):
    text = "\n".join(l for l in f.read_text().splitlines() if not l.lstrip().startswith("%"))
    found += [(k or "ask", f.name, b) for k, b in TEX.findall(text)]

found = [x for x in found if not only_auto or x[0] == "auto"]
for kind, src, body in found:
    print(f"[{kind}] {src}: {' '.join(body.split())}")
if not found:
    print("No open next steps.")
