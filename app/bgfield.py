"""The page background: a faint field of traced chromatin, the same one as on the Bohrer Lab website (chb-bohrer.github.io, its `.fieldbg` canvas).

Fibres are random walks of loci that wander back towards a centre; consecutive loci are joined by thin lines and each locus is a small cloud of dots whose
colour runs from grey through blue-grey to red along the fibre. The website draws it in the browser; here the same recipe is drawn once with PIL into a
transparent PNG (`app/static/bgfield.png`, served by Streamlit's static file serving, see .streamlit/config.toml) that `style.CSS` puts behind the whole page.
Regenerate with `scripts\\run.cmd python scripts\\make_bgfield.py` (deterministic for a given seed and size).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

STATIC = Path(__file__).resolve().parent / "static"
FILE = STATIC / "bgfield.png"

# the website's dark-theme values (--f1, --f2, --f3, --fline, --falpha)
STOPS = [(0x5E, 0x66, 0x74), (0x7A, 0x82, 0x90), (0xE8, 0x50, 0x6A)]
LINE = (230, 232, 238, round(0.09 * 255))
DOT_ALPHA = 0.5
DENSITY = 70000          # one fibre per this many pixels (the website's value)


def _colour(t: float) -> tuple:
    t = min(1.0, max(0.0, t)) * 2
    k = min(1, int(t))
    f = t - k
    a, b = STOPS[k], STOPS[k + 1]
    return tuple(round(x + (y - x) * f) for x, y in zip(a, b))


def render(width: int = 1920, height: int = 1080, seed: int = 67, supersample: int = 2) -> Image.Image:
    """Transparent RGBA image of the field (drawn at `supersample` x the size and reduced, so lines and dots are smooth)."""
    rng = np.random.default_rng(seed)
    s = supersample
    img = Image.new("RGBA", (width * s, height * s), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    for _ in range(max(4, round(width * height / DENSITY))):
        cx, cy, ang = rng.random() * width, rng.random() * height, rng.random() * 2 * np.pi
        loci = []
        for j in range(4 + int(rng.random() * 4)):                      # domains along the fibre
            if j:
                ang += rng.standard_normal() * 0.9
                step = 55 + rng.random() * 65
                cx += np.cos(ang) * step
                cy += np.sin(ang) * step * 0.75
            px, py = cx + rng.standard_normal() * 8, cy + rng.standard_normal() * 8
            for _ in range(5 + int(rng.random() * 6)):                  # loci of one domain: a walk pulled back to the domain centre
                px += rng.standard_normal() * 13
                py += rng.standard_normal() * 13
                px += 0.24 * (cx - px)
                py += 0.24 * (cy - py)
                loci.append((px, py))
        for a, b in zip(loci[:-1], loci[1:]):
            d.line([(a[0] * s, a[1] * s), (b[0] * s, b[1] * s)], fill=LINE, width=s)
        for i, (px, py) in enumerate(loci):
            col = _colour(i / (len(loci) - 1)) + (round(DOT_ALPHA * 255),)
            for _ in range(9 + int(rng.random() * 10)):
                x0, y0 = (px + rng.standard_normal() * 3.4 - 0.8) * s, (py + rng.standard_normal() * 3.4 - 0.8) * s
                d.rectangle([x0, y0, x0 + 1.6 * s, y0 + 1.6 * s], fill=col)
    if s > 1:                                                            # premultiplied, so transparent pixels do not darken the edges
        img = img.convert("RGBa").resize((width, height), Image.LANCZOS).convert("RGBA")
    return img


def write(path: Path = FILE, **kw) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    render(**kw).save(path, optimize=True)
    return path
