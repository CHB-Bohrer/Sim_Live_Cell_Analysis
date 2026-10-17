import asyncio, json, subprocess, re, pathlib, sys
import edge_tts, imageio_ffmpeg
from playwright.sync_api import sync_playwright
BASE = "file:///C:/Users/cbohr/OneDrive/Documents/Sim_Live_Cell_Analysis/docs/tutorial/index.html"
OUT = pathlib.Path("vid/audio"); OUT.mkdir(parents=True, exist_ok=True)
FF = imageio_ffmpeg.get_ffmpeg_exe()
VOICE = sys.argv[1] if len(sys.argv) > 1 else "en-US-JennyNeural"
with sync_playwright() as p:
    b = p.chromium.launch(channel="msedge", headless=True); pg = b.new_page(); pg.goto(BASE)
    beats = pg.evaluate("T.scenes.flatMap(s => s.beats.map((b, i) => ({id: s.id + ':' + i, text: b.s})))"); b.close()
print(len(beats), "beats")
sem = asyncio.Semaphore(4)
async def one(bt):
    f = OUT / (bt["id"].replace(":", "_") + ".mp3")
    if f.exists() and f.stat().st_size > 1000: return
    async with sem:
        for attempt in range(4):
            try:
                await edge_tts.Communicate(bt["text"], VOICE, rate="-4%").save(str(f)); return
            except Exception as e:
                print("retry", bt["id"], e); await asyncio.sleep(2)
        raise SystemExit("failed " + bt["id"])
async def main(): await asyncio.gather(*[one(b) for b in beats])
asyncio.run(main())
dur = {}
for bt in beats:
    f = OUT / (bt["id"].replace(":", "_") + ".mp3")
    r = subprocess.run([FF, "-i", str(f)], capture_output=True, text=True)
    m = re.search(r"Duration: (\d+):(\d+):([\d.]+)", r.stderr)
    dur[bt["id"]] = int(m[1]) * 3600 + int(m[2]) * 60 + float(m[3])
json.dump({"voice": VOICE, "beats": beats, "dur": dur}, open("vid/beats.json", "w"), indent=0)
print("total audio minutes:", round(sum(dur.values()) / 60, 1))
