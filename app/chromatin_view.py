"""Fast renderer (numpy + PIL) that turns a saved chromatin simulation into movie frames for the in-browser player.

Two modes:
  * "snapshots": the WHOLE chain (every monomer) at each saved snapshot (every `full_every_blocks` production blocks),
    with the loop-extruding factors drawn as lines between their two legs.
  * "blocks": every production block (full time resolution) but only every `track_stride`-th monomer.
Beads are coloured by genomic position (blue -> yellow) and shaded by depth; the chain slowly rotates so the 3D shape is
visible. Two chosen loci are highlighted, joined by a line, with their distance printed.
"""
from pathlib import Path

import matplotlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont

try:
    FONT = ImageFont.load_default(size=14)
except TypeError:
    FONT = ImageFont.load_default()


def _lut(n: int) -> np.ndarray:
    return (matplotlib.colormaps["viridis"](np.linspace(0, 1, max(n, 2)))[:, :3] * 255).astype(np.float32)


def project(pos: np.ndarray, angle: float, radius: float, size: int):
    """Rotate about the vertical axis, orthographic projection. Returns pixel x, y and depth (larger = nearer)."""
    c, s = np.cos(angle), np.sin(angle)
    x = pos[:, 0] * c + pos[:, 2] * s
    z = -pos[:, 0] * s + pos[:, 2] * c
    k = 0.46 * size / radius
    return x * k + size / 2, pos[:, 1] * k + size / 2, z / radius


def render_frame(pos: np.ndarray, colors: np.ndarray, angle: float, radius: float, size: int = 620,
                 loci: dict | None = None, loops: np.ndarray | None = None, caption: str = "") -> np.ndarray:
    """pos (n,3), colors (n,3) 0-255. loci: {name: (index_in_pos, (r,g,b))}. loops: (m,2) indices into pos."""
    px, py, depth = project(pos, angle, radius, size)
    order = np.argsort(depth)                                   # far first, so near beads overwrite
    shade = (0.45 + 0.55 * (depth[order] + 1) / 2)[:, None]
    img = np.zeros((size, size, 3), np.float32)
    xi, yi = px[order].astype(int), py[order].astype(int)
    col = colors[order] * shade
    for dx in (0, 1):
        for dy in (0, 1):
            ok = (xi + dx >= 0) & (xi + dx < size) & (yi + dy >= 0) & (yi + dy < size)
            img[yi[ok] + dy, xi[ok] + dx] = col[ok]
    pil = Image.fromarray(img.astype(np.uint8))
    d = ImageDraw.Draw(pil, "RGBA")
    d.ellipse([size / 2 - 0.46 * size, size / 2 - 0.46 * size, size / 2 + 0.46 * size, size / 2 + 0.46 * size],
              outline=(120, 120, 120, 90))                       # the confining sphere (nucleus-like boundary)
    if loops is not None and len(loops):
        for a, b in loops:
            d.line([px[a], py[a], px[b], py[b]], fill=(255, 255, 255, 110), width=1)
    pts = []
    if loci:
        for name, (i, rgb) in loci.items():
            pts.append((px[i], py[i]))
            r = 7
            d.ellipse([px[i] - r, py[i] - r, px[i] + r, py[i] + r], fill=rgb + (255,), outline=(255, 255, 255, 255), width=2)
            d.text((px[i] + 10, py[i] - 8), name, fill=(255, 255, 255, 255), font=FONT, stroke_width=2, stroke_fill=(0, 0, 0, 255))
        if len(pts) == 2:
            d.line([pts[0][0], pts[0][1], pts[1][0], pts[1][1]], fill=(255, 255, 255, 200), width=2)
    d.rectangle([0, 0, size, 24], fill=(0, 0, 0, 180))
    d.text((8, 4), caption, fill=(240, 240, 240, 255), font=FONT)
    return np.asarray(pil)


def snapshot_frames(sim_dir: str, a_kb: int, b_kb: int, show_loops: bool, rotate_deg: float, max_frames: int = 120):
    """Whole chain at each saved snapshot."""
    import json
    import yaml
    from polychrom.hdf5_format import list_URIs, load_URI

    d = Path(sim_dir)
    meta = json.loads((d / "meta.json").read_text())
    cfg = yaml.safe_load((d / "params.yaml").read_text())
    full_every = int(cfg["dynamics"].get("full_every_blocks", 10))
    N, R = meta["n_monomers"], meta["confinement_radius"]
    uris = list_URIs(str(d / "blocks"))
    step = max(1, int(np.ceil(len(uris) / max_frames)))
    lef = np.load(d / "lef_positions.npy", mmap_mode="r") if (d / "lef_positions.npy").exists() else None
    colors = _lut(N)
    frames, labels = [], []
    for k in range(0, len(uris), step):
        pos = np.asarray(load_URI(uris[k])["pos"], np.float32)
        pos = pos - pos.mean(0)
        p = k * full_every
        loops = np.asarray(lef[p]) if (lef is not None and show_loops and p < len(lef)) else None
        dist = float(np.linalg.norm(pos[a_kb] - pos[b_kb]))
        cap = f"production block {p}   |   locus A ({a_kb} kb) to B ({b_kb} kb): {dist:.1f} monomer diameters"
        frames.append(render_frame(pos, colors, np.deg2rad(rotate_deg) * len(frames), R,
                                   loci={"A": (a_kb, (230, 60, 60)), "B": (b_kb, (60, 200, 255))}, loops=loops, caption=cap))
        labels.append(f"snapshot {len(frames)}  (block {p})   distance A-B {dist:.1f}")
    return frames, labels


def block_frames(sim_dir: str, a_kb: int, b_kb: int, show_loops: bool, rotate_deg: float, max_frames: int = 250):
    """Every production block (time resolution of the saved per-block data), every `track_stride`-th monomer."""
    import json

    d = Path(sim_dir)
    meta = json.loads((d / "meta.json").read_text())
    traj = np.load(d / "tracked_positions.npy", mmap_mode="r")
    beads = np.load(d / "tracked_beads.npy")
    lef = np.load(d / "lef_positions.npy", mmap_mode="r") if (d / "lef_positions.npy").exists() else None
    R, T = meta["confinement_radius"], len(traj)
    ia, ib = int(np.argmin(np.abs(beads - a_kb))), int(np.argmin(np.abs(beads - b_kb)))
    stride = int(beads[1] - beads[0]) if len(beads) > 1 else 1
    colors = _lut(len(beads))
    step = max(1, int(np.ceil(T / max_frames)))
    frames, labels = [], []
    for t in range(0, T, step):
        pos = np.asarray(traj[t], np.float32)
        pos = pos - pos.mean(0)
        loops = None
        if lef is not None and show_loops:
            legs = np.clip(np.rint(np.asarray(lef[t]) / stride).astype(int), 0, len(beads) - 1)   # snap legs to stored beads
            loops = legs[legs[:, 0] != legs[:, 1]]
        dist = float(np.linalg.norm(pos[ia] - pos[ib]))
        cap = f"block {t}/{T}   |   locus A ({beads[ia]} kb) to B ({beads[ib]} kb): {dist:.1f} monomer diameters"
        frames.append(render_frame(pos, colors, np.deg2rad(rotate_deg) * len(frames), R,
                                   loci={"A": (ia, (230, 60, 60)), "B": (ib, (60, 200, 255))}, loops=loops, caption=cap))
        labels.append(f"block {t}   distance A-B {dist:.1f}")
    return frames, labels
