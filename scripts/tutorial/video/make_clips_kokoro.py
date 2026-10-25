"""Record the tutorial's narration with Kokoro, a free text-to-speech model that runs on this PC (no text leaves the machine).

Why: the Microsoft online voices still sound synthetic; Kokoro (82M parameters, Apache 2.0) sounds noticeably more human.

    1.  <shot venv>\\python make_clips.py --export beats.json          (reads every sentence from the player; needs playwright)
    2.  <kokoro venv>\\python make_clips_kokoro.py beats.json [voice] [speed]

Setup, done once (2026-10-25; a venv of its own, NOT the simlive env):
    python -m venv C:\\Users\\<you>\\kokoro_tmp\\v ; pip install kokoro-onnx soundfile numpy
    model files from https://github.com/thewh1teagle/kokoro-onnx/releases/tag/model-files-v1.0 into C:\\Users\\<you>\\kokoro_tmp\\model :
    kokoro-v1.0.onnx (325 MB, sha256 7d5df8ecf7d4b1878015a32686053fd0eebe2bc377234608764cc0ef3636a6c5) and voices-v1.0.bin (28 MB, sha256 bca610b8308e8d99f32e6fe4197e7ec01679264efed0cac9140fe9c29f1fbf7d).

Writes docs/tutorial/audio/<scene>_<n>.mp3 + manifest.js, exactly like make_clips.py (same hash check: a clip is used only while its sentence is unchanged).
Only new / changed sentences are synthesized. Voice default af_heart (others: af_bella, am_michael, bf_emma ...). MP3 encoding uses the ffmpeg that imageio-ffmpeg ships in the
shot venv (FFMPEG below).  PRON respells names the model would otherwise misread; it only changes what is spoken, never the caption.
"""
import json
import pathlib
import re
import subprocess
import sys

import numpy as np
import soundfile as sf
from kokoro_onnx import Kokoro

REPO = pathlib.Path(__file__).resolve().parents[3]
OUT = REPO / "docs" / "tutorial" / "audio"
OUT.mkdir(exist_ok=True)
MODEL = pathlib.Path.home() / "kokoro_tmp" / "model"
FFMPEG = pathlib.Path.home() / "shotenv_tmp" / "v" / "Lib" / "site-packages" / "imageio_ffmpeg" / "binaries" / "ffmpeg-win-x86_64-v7.1.exe"
beats = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
VOICE = sys.argv[2] if len(sys.argv) > 2 else "af_heart"
SPEED = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
TAG = f"kokoro:{VOICE}:{SPEED}"

PRON = {                                      # spoken form of words the model would misread (whole words only)
    "Viterbi": "Vitterbee", "Cellpose": "Cell pose", "Trackastra": "Track astra", "Vutara": "Voo tara", "Gillespie": "Gillespie",
    "Mirny": "Meernee", "Levo": "Leevo", "Bruker": "Brooker", "CHOTA": "Chotta", "OpenMM": "Open M M", "YAML": "Yammel",
    "LNK": "L N K", "TRA": "T R A", "DET": "D E T", "VXL": "V X L", "MS2": "M S two", "polychrom": "poly chrome", "Pol": "Pol",
}
pron = re.compile(r"\b(" + "|".join(map(re.escape, PRON)) + r")\b")


def spoken(text: str) -> str:
    return pron.sub(lambda m: PRON[m.group(0)], text)


vf = OUT / "voice.txt"
old = {}
mf = OUT / "manifest.js"
if mf.exists() and vf.exists() and vf.read_text(encoding="utf-8").strip() == TAG:        # same voice: keep clips of unchanged sentences
    try:
        txt = mf.read_text(encoding="utf-8").split("=", 1)[1]
        old = json.loads(txt[: txt.rindex(";")])
    except (ValueError, IndexError):
        old = {}
k = Kokoro(str(MODEL / "kokoro-v1.0.onnx"), str(MODEL / "voices-v1.0.bin"))
done = 0
for i, b in enumerate(beats):
    f = OUT / (b["key"].replace(":", "_") + ".mp3")
    o = old.get(b["key"])
    if o and o["h"] == b["h"] and f.exists() and f.stat().st_size > 1000:
        continue
    samples, sr = k.create(spoken(b["text"]), voice=VOICE, speed=SPEED, lang="en-us")
    samples = np.concatenate([samples, np.zeros(int(0.12 * sr), dtype=samples.dtype)])          # a breath of silence at the end
    wav = f.with_suffix(".wav")
    sf.write(wav, samples, sr)
    subprocess.run([str(FFMPEG), "-y", "-v", "error", "-i", str(wav), "-ac", "1", "-c:a", "libmp3lame", "-b:a", "64k", str(f)], check=True)
    wav.unlink()
    done += 1
    if done % 10 == 0:
        print(f"{i + 1}/{len(beats)}", flush=True)
clips = {b["key"]: {"f": "audio/" + b["key"].replace(":", "_") + ".mp3", "h": b["h"]} for b in beats}
mf.write_text("window.CLIPS = " + json.dumps(clips) + ";\n", encoding="utf-8")
vf.write_text(TAG, encoding="utf-8")
used = {b["key"].replace(":", "_") + ".mp3" for b in beats}
for f in OUT.glob("*.mp3"):
    if f.name not in used:
        f.unlink()
print("sentences recorded:", done, "| clips:", len(beats), "| MB:", round(sum(f.stat().st_size for f in OUT.glob("*.mp3")) / 1e6, 1))
