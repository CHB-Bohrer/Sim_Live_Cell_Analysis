"""Build vid/audio + vid/beats.json from the recorded Kokoro clips in docs/tutorial/audio (instead of tts.py's Microsoft voice).
Sentence list read from the player the same way tts.py does. Checks each clip's hash against the manifest (window.CLIPS) and the page's hashText.
Run from the video working folder (the one that gets vid/) with the shot venv python: python <repo>/scripts/tutorial/video/beats_from_kokoro.py"""
import json, pathlib, re, shutil, subprocess, sys
import imageio_ffmpeg
from playwright.sync_api import sync_playwright

CODE = pathlib.Path(__file__).resolve().parents[3]
TUT = CODE / "docs" / "tutorial"
BASE = (TUT / "index.html").as_uri()
OUT = pathlib.Path("vid/audio"); OUT.mkdir(parents=True, exist_ok=True)
FF = imageio_ffmpeg.get_ffmpeg_exe()

txt = (TUT / "audio" / "manifest.js").read_text(encoding="utf-8")
clips = json.loads(txt[txt.index("{"): txt.rindex(";")])
with sync_playwright() as p:
    b = p.chromium.launch(channel="msedge", headless=True); pg = b.new_page(); pg.goto(BASE)
    beats = pg.evaluate("T.scenes.flatMap(s => s.beats.map((b, i) => ({id: s.id + ':' + i, text: b.s})))")
    has_hash = pg.evaluate("typeof hashText === 'function'")
    hashes = pg.evaluate("(ts) => ts.map(t => hashText(t))", [bt["text"] for bt in beats]) if has_hash else None
    b.close()
print(len(beats), "beats; hash check", "on" if has_hash else "OFF (hashText not global)")
missing, mismatch, dur = [], [], {}
for i, bt in enumerate(beats):
    c = clips.get(bt["id"])
    if c is None:
        missing.append(bt["id"]); continue
    if hashes is not None and str(c.get("h")) != str(hashes[i]):
        mismatch.append(bt["id"])
    src = TUT / c["f"]
    dst = OUT / (bt["id"].replace(":", "_") + ".mp3")
    if not src.exists():
        missing.append(bt["id"]); continue
    shutil.copyfile(src, dst)
    r = subprocess.run([FF, "-i", str(dst)], capture_output=True, text=True)
    m = re.search(r"Duration: (\d+):(\d+):([\d.]+)", r.stderr)
    dur[bt["id"]] = int(m[1]) * 3600 + int(m[2]) * 60 + float(m[3])
print("missing clips:", missing or "none")
print("hash mismatch:", mismatch or "none")
if missing:
    sys.exit("stop: missing clips")
json.dump({"voice": "kokoro af_heart", "beats": beats, "dur": dur}, open("vid/beats.json", "w"), indent=0)
print("total audio minutes:", round(sum(dur.values()) / 60, 1))
