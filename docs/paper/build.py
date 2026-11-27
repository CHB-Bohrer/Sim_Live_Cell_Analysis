r"""Build docs/paper/sim-live-paper.pdf from the section files.

    scripts\run.cmd python docs\paper\build.py              # redraw stale figures, then build the PDF
    scripts\run.cmd python docs\paper\build.py --no-figures # build from the figures already on disk
    scripts\run.cmd python docs\paper\build.py --mark-done  # also record this commit as "documented"

Needs pandoc on PATH, matplotlib, and Playwright with Chromium (the tutorial's shot venv has
Playwright; or set CHROMIUM_PATH to an existing Chrome / Edge executable).
The build never runs a simulation or a check: figures come from docs/paper/data/.
Sections are the files in sections/, in file-name order. A line `@include data/<file>.md`
is replaced by that file (generated tables). Same system as chromatin-studio's paper.
"""
import argparse
import datetime
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_PDF = HERE / "sim-live-paper.pdf"
BUILD = HERE / "build"
STATE = HERE / "data" / "paper_state.json"


def commit():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    except OSError:
        return ""


def assemble():
    parts = []
    for path in sorted((HERE / "sections").glob("*.md")):
        lines = []
        for line in path.read_text(encoding="utf-8").splitlines():
            m = re.match(r"^@include\s+(\S+)\s*$", line)
            if m:
                target = HERE / m.group(1)
                lines.append(target.read_text(encoding="utf-8") if target.exists() else f"*({m.group(1)} has not been generated yet.)*")
            else:
                lines.append(line)
        parts.append("\n".join(lines))
    return "\n\n".join(parts)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--html-only", action="store_true")
    parser.add_argument("--mark-done", action="store_true")
    args = parser.parse_args()
    if not args.no_figures:
        subprocess.run([sys.executable, str(HERE / "figures" / "make_figures.py")], check=True)
    BUILD.mkdir(exist_ok=True)
    source = BUILD / "paper.md"
    source.write_text(assemble(), encoding="utf-8")
    html = BUILD / "paper.html"
    today = datetime.date.today().isoformat()
    subprocess.run([
        "pandoc", str(source), "--from", "markdown+fenced_divs+header_attributes+implicit_figures",
        "--to", "html5", "--standalone", "--mathml", "--number-sections", "--toc", "--toc-depth=1",
        "--resource-path", str(HERE), "--css", str(HERE / "style.css"), "--embed-resources",
        "--metadata", "title=Sim Live Cell Analysis: how the pipeline works",
        "--metadata", "subtitle=A plain-language technical description, with the checks that each part works",
        "--metadata", f"date=Built {today} from commit {commit() or 'unknown'}",
        "--metadata", "toc-title=Contents", "--output", str(html)], check=True, cwd=HERE)
    if args.html_only:
        print(html)
        return
    from playwright.sync_api import sync_playwright
    footer = ('<div style="font-size:7.5pt;color:#777;width:100%;text-align:center;font-family:sans-serif">'
              'Sim Live Cell Analysis · how the pipeline works · page <span class="pageNumber"></span> of '
              '<span class="totalPages"></span></div>')
    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=os.environ.get("CHROMIUM_PATH") or None)
        page = browser.new_page()
        page.goto(html.as_uri(), wait_until="load")
        page.pdf(path=str(OUT_PDF), format="A4", print_background=True, display_header_footer=True,
                 header_template="<div></div>", footer_template=footer,
                 margin={"top": "18mm", "bottom": "18mm", "left": "17mm", "right": "17mm"})
        browser.close()
    if args.mark_done:
        STATE.write_text(json.dumps({"documented_commit": commit(), "date": today}, indent=1) + "\n", newline="\n")
    print(OUT_PDF)


if __name__ == "__main__":
    main()
