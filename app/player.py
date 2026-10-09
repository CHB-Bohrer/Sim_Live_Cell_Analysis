"""A movie player that runs entirely in the browser.

All frames are JPEG-encoded once and embedded in the page; play / pause / scrub / speed then run client-side with no
server round trip per frame, so playback is smooth (the old approach re-ran Python and matplotlib for every frame).
"""
import base64
import io
import json

import numpy as np
import streamlit.components.v1 as components
from PIL import Image


def _jpeg_b64(arr: np.ndarray, quality: int = 85) -> str:
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="JPEG", quality=quality)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()


def show_player(frames: list, labels: list | None = None, fps: int = 10, column_width_px: int = 900, key: str = "p"):
    """frames: list of uint8 RGB arrays (all the same shape). labels: one caption per frame (e.g. frame number)."""
    if not frames:
        return
    h, w = frames[0].shape[:2]
    uris = [_jpeg_b64(f) for f in frames]
    labels = labels or [f"frame {i}" for i in range(len(frames))]
    # The speed menu must contain `fps`: with fps=8 and no "8" option the select's value was "" -> NaN delay -> it never played.
    speeds = sorted({2, 5, 10, 20, 40, int(fps)})
    options = "".join(f'<option value="{s}">{s} fps</option>' for s in speeds)
    html = f"""
<div style="font-family:system-ui,sans-serif;color:#444">
  <img id="im{key}" style="width:100%;display:block;border-radius:4px;background:#000" />
  <div style="display:flex;gap:10px;align-items:center;margin-top:8px">
    <button id="pp{key}" style="padding:4px 14px;font-size:14px;cursor:pointer">&#9654; Play</button>
    <input id="sl{key}" type="range" min="0" max="{len(frames) - 1}" value="0" style="flex:1" />
    <span id="lb{key}" style="min-width:210px;font-size:13px"></span>
    <select id="fp{key}" style="font-size:13px">{options}</select>
  </div>
</div>
<script>
(function() {{
  const F = {json.dumps(uris)}, L = {json.dumps(labels)};
  const im = document.getElementById("im{key}"), sl = document.getElementById("sl{key}"),
        lb = document.getElementById("lb{key}"), pp = document.getElementById("pp{key}"),
        fp = document.getElementById("fp{key}");
  fp.value = "{fps}";
  let i = 0, playing = false, last = 0;
  F.forEach(u => {{ const p = new Image(); p.src = u; }});          // preload every frame
  function show(k) {{ i = k; im.src = F[k]; sl.value = k; lb.textContent = L[k]; }}
  function tick() {{          // a timer, not requestAnimationFrame: rAF is paused in hidden / embedded panes
    const ts = performance.now();
    if (playing && ts - last >= 1000 / (parseFloat(fp.value) || {int(fps)})) {{ last = ts; show((i + 1) % F.length); }}
  }}
  pp.onclick = () => {{ playing = !playing; pp.innerHTML = playing ? "&#10074;&#10074; Pause" : "&#9654; Play"; }};
  sl.oninput = () => {{ playing = false; pp.innerHTML = "&#9654; Play"; show(parseInt(sl.value)); }};
  document.addEventListener("keydown", e => {{
    if (e.key === "ArrowRight") show(Math.min(i + 1, F.length - 1));
    if (e.key === "ArrowLeft") show(Math.max(i - 1, 0));
    if (e.key === " ") {{ e.preventDefault(); pp.click(); }}
  }});
  show(0); setInterval(tick, 20);
}})();
</script>"""
    height = int(column_width_px * h / w) + 70
    components.html(html, height=height, scrolling=False)
