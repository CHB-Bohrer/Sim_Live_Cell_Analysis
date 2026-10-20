"""Place each sentence's MP3 at the moment its caption appears in the recorded video, then encode H.264 + AAC.

The page logs when each sentence starts (vid/beatlog.json), but the recorder's frames reach the video 0.3-1.6 s late, by a
varying amount (measured 2026-10-20, median 0.7 s). Placing the audio at the logged time therefore makes the voice lead the
picture. So, unless --no-sync is given, the time at which the caption box changes is found in the recording for every
sentence (0.04 s steps) and the audio is placed there; a sentence never starts before the previous one has ended.
Usage (from the folder holding vid/): python mix.py <out.mp4> [--no-sync]
"""
import json, subprocess, numpy as np, glob, sys
from scipy.io import wavfile
import imageio_ffmpeg
FF = imageio_ffmpeg.get_ffmpeg_exe(); SR = 24000
CAPTION = (1100, 85, 388, 750)       # w, h, x, y of the caption box in the 1600 x 900 recording
FPS = 25

L = json.load(open("vid/beatlog.json")); off = L["origin_minus_created"]
starts = {k: off + ms / 1000 for k, ms in L["log"]}
vid = glob.glob("vid/rec/*.webm")[0]
out = sys.argv[1]


def caption_change(t):
    """Seconds after the logged start t at which the caption box first differs from how it looked 0.3 s before t (None if not found)."""
    w, h, x, y = CAPTION
    raw = subprocess.run([FF, "-v", "error", "-ss", f"{max(t - 0.3, 0):.3f}", "-i", vid, "-t", "3.3", "-vf", f"fps={FPS},crop={w}:{h}:{x}:{y},format=gray",
                          "-f", "rawvideo", "-"], capture_output=True).stdout
    fr = np.frombuffer(raw, dtype=np.uint8)
    if fr.size < w * h * 2:
        return None
    fr = fr[: fr.size // (w * h) * w * h].reshape(-1, h, w).astype(np.float32)
    d = np.abs(fr - fr[0]).mean(axis=(1, 2))
    hit = np.nonzero(d > 3)[0]
    return None if hit.size == 0 else hit[0] / FPS - 0.3


placed = dict(starts)
if "--no-sync" not in sys.argv:
    lags = {}
    for i, k in enumerate(sorted(starts, key=starts.get)):
        lags[k] = caption_change(starts[k]) if i > 0 else None
        if i % 20 == 0:
            print("sync", i, "of", len(starts), flush=True)
    found = [v for v in lags.values() if v is not None and v >= 0]
    med = float(np.median(found)) if found else 0.0
    for k, v in lags.items():
        placed[k] = starts[k] + (v if v is not None and v >= 0 else med)
    print(f"caption lag: median {med:.2f} s, range {min(found):.2f}-{max(found):.2f} s, not found {sum(v is None for v in lags.values()) - 1}")
    json.dump({k: round(v, 3) if v is not None else None for k, v in lags.items()}, open("vid/caption_lag.json", "w"), indent=0)

end = max(placed.values()) + 30
buf = np.zeros(int(end * SR), dtype=np.float32)
prev_end = 0
for k in sorted(placed, key=placed.get):
    f = f"vid/audio/{k.replace(':', '_')}.mp3"
    raw = subprocess.run([FF, "-v", "error", "-i", f, "-f", "f32le", "-ac", "1", "-ar", str(SR), "-"], capture_output=True).stdout
    a = np.frombuffer(raw, dtype=np.float32)
    i = max(int(placed[k] * SR), prev_end + int(0.1 * SR))     # never on top of the previous sentence
    buf[i:i + len(a)] += a[: len(buf) - i]
    prev_end = i + len(a)
buf = buf[: prev_end + 5 * SR]
wavfile.write("vid/narration.wav", SR, (np.clip(buf, -1, 1) * 32767).astype(np.int16))
print("audio seconds", len(buf) / SR)
cmd = [FF, "-y", "-v", "error", "-stats", "-i", vid, "-i", "vid/narration.wav", "-c:v", "libx264", "-preset", "medium", "-crf", "27", "-pix_fmt", "yuv420p", "-r", "25",
       "-c:a", "aac", "-b:a", "128k", "-movflags", "+faststart", "-shortest", out]
subprocess.run(cmd, check=True)
print("wrote", out)
