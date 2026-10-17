import json, subprocess, numpy as np, pathlib, glob, sys
from scipy.io import wavfile
import imageio_ffmpeg
FF = imageio_ffmpeg.get_ffmpeg_exe(); SR = 24000
L = json.load(open("vid/beatlog.json")); off = L["origin_minus_created"]
starts = {k: off + ms / 1000 for k, ms in L["log"]}
end = max(starts.values()) + 10
buf = np.zeros(int((end + 5) * SR), dtype=np.float32)
for k, t in starts.items():
    f = f"vid/audio/{k.replace(':', '_')}.mp3"
    raw = subprocess.run([FF, "-v", "error", "-i", f, "-f", "f32le", "-ac", "1", "-ar", str(SR), "-"], capture_output=True).stdout
    a = np.frombuffer(raw, dtype=np.float32); i = int(t * SR)
    buf[i:i + len(a)] += a[: len(buf) - i]
wavfile.write("vid/narration.wav", SR, (np.clip(buf, -1, 1) * 32767).astype(np.int16))
print("audio seconds", len(buf) / SR)
vid = glob.glob("vid/rec/*.webm")[0]
out = sys.argv[1]
cmd = [FF, "-y", "-v", "error", "-stats", "-i", vid, "-i", "vid/narration.wav", "-c:v", "libx264", "-preset", "medium", "-crf", "27", "-pix_fmt", "yuv420p", "-r", "25",
       "-c:a", "aac", "-b:a", "128k", "-movflags", "+faststart", "-shortest", out]
subprocess.run(cmd, check=True)
print("wrote", out)
