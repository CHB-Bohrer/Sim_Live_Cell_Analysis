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
<div style="font-family:system-ui,sans-serif;color:#444;max-width:{column_width_px}px;margin:0 auto">
  <div id="ld{key}" style="font-size:13px;margin-bottom:6px;color:#92400e">&#8987; Loading frames&hellip; <span id="ln{key}">0 / {len(frames)}</span>
    <div style="height:5px;background:#e5e7eb;border-radius:3px;margin-top:3px"><div id="lbar{key}" style="height:5px;width:0%;background:#0f766e;border-radius:3px"></div></div></div>
  <img id="im{key}" style="width:100%;display:block;border-radius:4px;background:#000" />
  <div style="display:flex;gap:10px;align-items:center;margin-top:8px;flex-wrap:wrap">
    <button id="pp{key}" disabled style="padding:4px 14px;font-size:14px;cursor:wait;min-width:150px">&#8987; Loading&hellip;</button>
    <input id="sl{key}" type="range" min="0" max="{len(frames) - 1}" value="0" style="flex:1;min-width:160px" />
    <span id="lb{key}" style="min-width:150px;font-size:13px"></span>
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
  const ld = document.getElementById("ld{key}"), ln = document.getElementById("ln{key}"), lbar = document.getElementById("lbar{key}");
  let loaded = 0;
  function ready() {{
    ld.innerHTML = "&#10003; Ready: press Play (or Space), drag the slider, or use the arrow keys.";
    ld.style.color = "#166534"; pp.disabled = false; pp.style.cursor = "pointer"; pp.innerHTML = "&#9654; Play";
  }}
  F.forEach(u => {{ const p = new Image();                             // preload every frame; Play is enabled when all are ready
    p.onload = p.onerror = () => {{ loaded++; ln.textContent = loaded + " / " + F.length; lbar.style.width = (100 * loaded / F.length) + "%";
                                    if (loaded === F.length) ready(); }};
    p.src = u; }});
  function show(k) {{ i = k; im.src = F[k]; sl.value = k; lb.textContent = L[k]; }}
  function tick() {{          // a timer, not requestAnimationFrame: rAF is paused in hidden / embedded panes
    const ts = performance.now();
    if (playing && ts - last >= 1000 / (parseFloat(fp.value) || {int(fps)})) {{ last = ts; show((i + 1) % F.length); }}
  }}
  pp.onclick = () => {{ if (pp.disabled) return; playing = !playing; pp.innerHTML = playing ? "&#10074;&#10074; Pause" : "&#9654; Play"; }};
  sl.oninput = () => {{ playing = false; if (!pp.disabled) pp.innerHTML = "&#9654; Play"; show(parseInt(sl.value)); }};
  document.addEventListener("keydown", e => {{
    if (e.key === "ArrowRight") show(Math.min(i + 1, F.length - 1));
    if (e.key === "ArrowLeft") show(Math.max(i - 1, 0));
    if (e.key === " " && !pp.disabled) {{ e.preventDefault(); pp.click(); }}
  }});
  show(0); setInterval(tick, 20);
}})();
</script>"""
    height = int(column_width_px * h / w) + 150       # image (width capped at column_width_px, so its height is known) + loading line + controls
    components.html(html, height=height, scrolling=False)
