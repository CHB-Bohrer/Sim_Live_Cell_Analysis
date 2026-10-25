"""Record the tutorial's narration as natural neural-voice clips that the player (docs/tutorial/index.html) plays itself.

    python make_clips.py [voice]        (default en-US-JennyNeural; others: en-US-AvaMultilingualNeural, en-US-AriaNeural, en-GB-SoniaNeural ...)

Run with the separate venv that has playwright + edge-tts (see docs/tutorial/README.md).  It reads every sentence from the player itself,
writes docs/tutorial/audio/<scene>_<n>.mp3 and docs/tutorial/audio/manifest.js (window.CLIPS = {"scene:n": {f, h}}), where h is a hash of the sentence
(the same function as engine.js hashText): the player only uses a clip whose hash matches the sentence, so an edited sentence falls back to the browser
voice until this script is run again. Only sentences that changed (or whose clip is missing) are synthesized, so re-running after an edit is quick.

The narration text is sent to Microsoft's online text-to-speech service (edge-tts); it is already public in this repository.
The clips are NOT committed (docs/tutorial/audio is in .gitignore, ~20 MB): publish them with a release if others need them.
"""
import asyncio
import json
import pathlib
import sys

import edge_tts
from playwright.sync_api import sync_playwright

TUT = pathlib.Path(__file__).resolve().parents[3] / "docs" / "tutorial"
OUT = TUT / "audio"
OUT.mkdir(exist_ok=True)
VOICE = sys.argv[1] if len(sys.argv) > 1 else "en-US-JennyNeural"

with sync_playwright() as p:
    b = p.chromium.launch(channel="msedge", headless=True)
    pg = b.new_page()
    pg.goto(TUT.joinpath("index.html").as_uri())
    beats = pg.evaluate("T.scenes.flatMap(s => s.beats.map((b, i) => ({key: s.id + ':' + i, text: b.s, h: hashText(b.s)})))")
    b.close()
print(len(beats), "sentences; voice", VOICE)

old = {}
mf = OUT / "manifest.js"
if mf.exists():
    try:
        txt = mf.read_text(encoding="utf-8").split("=", 1)[1]
        old = json.loads(txt[: txt.rindex(";")])
    except (ValueError, IndexError):
        old = {}
vf = OUT / "voice.txt"
if not vf.exists() or vf.read_text(encoding="utf-8").strip() != VOICE:     # a different voice: record everything again
    old = {}
sem = asyncio.Semaphore(4)


async def one(bt):
    f = OUT / (bt["key"].replace(":", "_") + ".mp3")
    o = old.get(bt["key"])
    if o and o["h"] == bt["h"] and f.exists() and f.stat().st_size > 1000:
        return
    async with sem:
        for attempt in range(5):
            try:
                await edge_tts.Communicate(bt["text"], VOICE).save(str(f))
                if f.stat().st_size > 1000:
                    return
            except Exception as e:  # noqa: BLE001  (network hiccups: retry)
                print("retry", bt["key"], e)
            await asyncio.sleep(2)
        raise SystemExit("failed " + bt["key"])


async def main():
    await asyncio.gather(*[one(b) for b in beats])

asyncio.run(main())
clips = {b["key"]: {"f": "audio/" + b["key"].replace(":", "_") + ".mp3", "h": b["h"]} for b in beats}
mf.write_text("window.CLIPS = " + json.dumps(clips) + ";\n", encoding="utf-8")
vf.write_text(VOICE, encoding="utf-8")
used = {b["key"].replace(":", "_") + ".mp3" for b in beats}
for f in OUT.glob("*.mp3"):                        # clips of sentences that no longer exist
    if f.name not in used:
        f.unlink()
print("clips:", len(beats), "total MB:", round(sum(f.stat().st_size for f in OUT.glob("*.mp3")) / 1e6, 1))
