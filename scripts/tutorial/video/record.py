import json, time, pathlib, shutil
from playwright.sync_api import sync_playwright
BASE = "file:///C:/Users/cbohr/OneDrive/Documents/Sim_Live_Cell_Analysis/docs/tutorial/index.html"
B = json.load(open("vid/beats.json"))
dur = {k: round(v + 0.12, 3) for k, v in B["dur"].items()}
vd = pathlib.Path("vid/rec"); shutil.rmtree(vd, ignore_errors=True); vd.mkdir()
with sync_playwright() as p:
    br = p.chromium.launch(channel="msedge", headless=True)
    ctx = br.new_context(viewport={"width": 1600, "height": 900}, record_video_dir=str(vd), record_video_size={"width": 1600, "height": 900})
    ctx.add_init_script("window.BEAT_DUR = " + json.dumps(dur) + ";")
    pg = ctx.new_page(); created = time.time()
    pg.goto(BASE + "#s=welcome&b=0&play=1")
    pg.wait_for_function("window.BEAT_LOG && window.BEAT_LOG.length > 0", timeout=60000)
    origin = pg.evaluate("performance.timeOrigin") / 1000
    print("recording started, offset", origin - created, flush=True)
    n = 0
    while True:
        time.sleep(30)
        n += 1
        st = pg.evaluate("[T.si, T.bi, T.playing]")
        if n % 10 == 0: print(round(n * 0.5), "min", st, flush=True)
        if not st[2]: break
    time.sleep(3)
    log = pg.evaluate("window.BEAT_LOG")
    json.dump({"origin_minus_created": origin - created, "log": log}, open("vid/beatlog.json", "w"))
    ctx.close(); br.close()
print("done", flush=True)
