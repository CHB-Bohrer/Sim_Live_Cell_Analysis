"""Redraw the dashboard's page background (the lab website's faint chromatin field) into app/static/bgfield.png.

    scripts\\run.cmd python scripts\\make_bgfield.py
    scripts\\run.cmd python scripts\\make_bgfield.py --seed 41 --width 2560 --height 1440

Restart the dashboard afterwards (the browser caches the image; Ctrl+F5 also helps).
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "app"))
import bgfield  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--width", type=int, default=1920)
    ap.add_argument("--height", type=int, default=1080)
    ap.add_argument("--seed", type=int, default=67)
    a = ap.parse_args()
    out = bgfield.write(width=a.width, height=a.height, seed=a.seed)
    print(f"wrote {out} ({out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
