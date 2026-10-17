"""Copy the captured screenshots into docs/tutorial/img and write docs/tutorial/shots_data.js (picture sizes + the
on-screen rectangles of every control, which the tutorial uses to point at things).

    python scripts/tutorial/build_data.py <shots folder made by capture_shots.py>
(plain Python 3, no packages)
"""
import json
import shutil
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(__file__).resolve().parents[2] / "docs" / "tutorial"
(dst / "img").mkdir(parents=True, exist_ok=True)
# Hand-set rectangles where the automatic text match found the wrong element (checked by eye on the captured pictures, 2026-10-17).
# [x, y, w, h] in picture pixels; "tail" = the last `h` pixels of the page.
INTRO = [90, 97, 1120, 142]                      # the two "What this shows / How to use it" boxes sit right under the page bar on every page
OVERRIDES = {
    "pg_scans": {"masks": [660, 380, 550, 68], "roles": [90, 455, 1120, 95]},
    "pg_images": {"frame": [90, 405, 370, 90]},
    "pg_images_err": {"what": [90, 2270, 370, 179]},
    "pg_chrom": {"lib": [90, 285, 1120, 75], "sim": [90, 375, 1120, 70]},
}
TAILS = {"pg_probes": ("use", 140)}
data = {}
for j in sorted(src.glob("*.json")):
    d = json.loads(j.read_text(encoding="utf-8"))
    jpg = j.with_suffix(".jpg")
    if not jpg.exists():
        continue
    shutil.copy2(jpg, dst / "img" / jpg.name)
    d["src"] = f"img/{jpg.name}"
    if j.stem.startswith("pg_"):
        d["rects"]["intro"] = INTRO
    d["rects"].update(OVERRIDES.get(j.stem, {}))
    if j.stem in TAILS:
        k, hh = TAILS[j.stem]
        d["rects"][k] = [90, d["h"] - hh - 10, 1120, hh]
    data[j.stem] = d
(dst / "shots_data.js").write_text("window.SHOTS = " + json.dumps(data, ensure_ascii=False) + ";\n", encoding="utf-8")
print(len(data), "shots;", round(sum(f.stat().st_size for f in (dst / "img").glob("*")) / 1e6, 1), "MB of images")
