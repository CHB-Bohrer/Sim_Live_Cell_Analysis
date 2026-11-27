#!/usr/bin/env python3
"""List parked thoughts and the condition for revisiting each (one line apiece).

  scripts\run.cmd python scripts\thoughts.py          # still parked
  scripts\run.cmd python scripts\thoughts.py --all
"""
import pathlib, re, sys

text = (pathlib.Path(__file__).resolve().parents[1] / "THOUGHTS.md").read_text()
text = re.sub(r"<!--.*?-->", "", text, flags=re.S)
n = 0
for block in re.split(r"^## ", text, flags=re.M)[1:]:
    head = block.splitlines()[0].strip()
    get = lambda k: (re.search(rf"\*\*{k}:\*\*\s*(.+)", block) or [None, "?"])[1].strip()
    status = get("Status")
    if "--all" in sys.argv or status.startswith("parked"):
        n += 1
        print(f"{head}  |  revisit when: {get('Revisit when')}" + ("" if status.startswith("parked") else f"  [{status}]"))
if not n:
    print("No parked thoughts.")
