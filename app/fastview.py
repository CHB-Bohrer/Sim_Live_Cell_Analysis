"""Fast frame renderers (numpy + PIL, ~10 ms per frame) used to pre-render movies for the in-browser player.

matplotlib is far too slow for this (hundreds of ms per figure); here overlays, outlines, ID text and circles are
drawn directly onto pixel arrays.
"""
import colorsys

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi
from skimage.segmentation import find_boundaries

try:
    FONT = ImageFont.load_default(size=13)
except TypeError:  # older Pillow
    FONT = ImageFont.load_default()
FONT_BIG = None


def color_for(i: int) -> np.ndarray:
    """Stable colour per ID (same rule as the dashboard), as 0-255 floats."""
    return 255 * np.array(colorsys.hsv_to_rgb((i * 0.61803398875) % 1.0, 0.65, 1.0))


def _lut(max_id: int, colour_of) -> np.ndarray:
    t = np.zeros((max_id + 1, 3))
    for k in range(1, max_id + 1):
        t[k] = colour_of(k)
    return t


def gray_to_rgb(img: np.ndarray, lo: float, hi: float) -> np.ndarray:
    g = np.clip((img.astype(float) - lo) / max(hi - lo, 1e-9), 0, 1) * 255
    return np.repeat(g[..., None], 3, axis=2)


def blend_labels(rgb: np.ndarray, lab: np.ndarray, lut: np.ndarray, alpha: float = 0.45) -> np.ndarray:
    rgb = rgb.astype(np.float32, copy=True)
    m = lab > 0                       # touch only labelled pixels (a small fraction of the image)
    rgb[m] = (1 - alpha) * rgb[m] + alpha * lut[lab[m]]
    return rgb


def outline(rgb: np.ndarray, mask: np.ndarray, colour) -> np.ndarray:
    if mask.any():
        rgb[find_boundaries(mask, mode="inner")] = colour
    return rgb


def draw_ids(pil: Image.Image, lab: np.ndarray, colour=(255, 255, 255)) -> None:
    ys, xs = np.nonzero(lab)
    if ys.size == 0:
        return
    ids, inv = np.unique(lab[ys, xs], return_inverse=True)
    n = np.bincount(inv)
    cy, cx = np.bincount(inv, ys) / n, np.bincount(inv, xs) / n   # centroids without scipy's slow per-label loop
    d = ImageDraw.Draw(pil)
    for i, y, x in zip(ids, cy, cx):
        d.text((x, y), str(i), fill=colour, font=FONT, anchor="mm", stroke_width=1, stroke_fill=(0, 0, 0))


def titled(pil: Image.Image, text: str, scale: float) -> Image.Image:
    """Resize a panel by `scale` and add a title strip on top."""
    w, h = max(1, int(pil.width * scale)), max(1, int(pil.height * scale))
    pil = pil.resize((w, h), Image.BILINEAR)
    out = Image.new("RGB", (w, h + 22), (24, 24, 24))
    out.paste(pil, (0, 22))
    ImageDraw.Draw(out).text((6, 4), text, fill=(235, 235, 235), font=FONT)
    return out


def hstack(panels: list, gap: int = 6) -> np.ndarray:
    h = max(p.height for p in panels)
    w = sum(p.width for p in panels) + gap * (len(panels) - 1)
    out = Image.new("RGB", (w, h), (0, 0, 0))
    x = 0
    for p in panels:
        out.paste(p, (x, 0)); x += p.width + gap
    return np.asarray(out)


# ----------------------------------------------------------------------------- full-field movie
def prepare_movie_state(R: dict, variant: str) -> dict:
    """Per-run lookups computed once (not per frame)."""
    V = R["v"][variant]
    m = V["matches"]
    ev = V["events"]
    cells = R["cells"].set_index(["cell_id", "t"]).sort_index()
    return {
        "lost": {t: g.gt_id.to_numpy() for t, g in m[m.pred_id == 0].groupby("t")},
        "cur": {t: dict(zip(g.pred_id, g.gt_id)) for t, g in m.groupby("t")},
        "events": {t: g for t, g in ev.groupby("t")},
        "cells": cells,
        "dom": V["dominant"],
        "parent_of": V["parent_of"],
        "lo": float(np.percentile(R["img"], 1)), "hi": float(np.percentile(R["img"], 99.8)),
        "gmax": int(R["gt"].max()),
    }


def movie_frame(R: dict, variant: str, t: int, st_: dict, show_image: bool = False, target_w: int = 560) -> np.ndarray:
    V = R["v"][variant]
    gt, trk, img = R["gt"][t], V["trk"][t], R["img"][t]
    px = R["cfg"]["optics"]["pixel_size_nm"] / 1000.0
    base = gray_to_rgb(img, st_["lo"], st_["hi"])
    scale = target_w / img.shape[1]

    # true cells: colour = true ID; yellow outline = missed by the segmentation
    g = blend_labels(base.copy(), gt, _lut(st_["gmax"], color_for))
    lost = st_["lost"].get(t, np.array([], int))
    outline(g, np.isin(gt, lost) & (gt > 0), (255, 230, 0))
    gp = Image.fromarray(g.astype(np.uint8)); draw_ids(gp, gt)

    # tracked cells: colour = the true cell the track mostly follows; purple = mask with no true cell under it
    dom = st_["dom"]
    pmax = int(max(trk.max(), 1))
    tr = blend_labels(base.copy(), trk, _lut(pmax, lambda p: color_for(dom[p]) if p in dom else (140, 140, 140)))
    cur = st_["cur"].get(t, {})
    nogt = [p for p in np.unique(trk) if p > 0 and p not in cur]
    outline(tr, np.isin(trk, nogt) & (trk > 0), (204, 51, 255))
    tp = Image.fromarray(tr.astype(np.uint8)); draw_ids(tp, trk)

    # identity switches at THIS frame: big red circle + label
    ev = st_["events"].get(t)
    if ev is not None:
        for r in ev.itertuples():
            try:
                row = st_["cells"].loc[(r.gt_id, t)]
            except KeyError:
                continue
            cx, cy, rad = row.x_um / px, row.y_um / px, 1.7 * row.bound_radius_um / px
            for pil in (gp, tp):
                ImageDraw.Draw(pil).ellipse([cx - rad, cy - rad, cx + rad, cy + rad], outline=(255, 30, 30),
                                            width=max(2, int(3 / scale)))
            par = st_["parent_of"].get(r.new_pred_id, 0)
            txt = f"track {r.old_pred_id}->{r.new_pred_id}" + (" (called a division)" if par > 0 else "")
            ImageDraw.Draw(tp).text((cx, cy - rad - 2), txt, fill=(255, 255, 255), font=FONT, anchor="md",
                                    stroke_width=2, stroke_fill=(220, 20, 20))
    panels = ([titled(Image.fromarray(base.astype(np.uint8)), "Image (what the tracker sees)", scale)] if show_image else [])
    panels += [titled(gp, "TRUE cells (ground-truth IDs)", scale), titled(tp, "TRACKED cells (tracker's IDs)", scale)]
    return hstack(panels)


def movie_frames(R: dict, variant: str, show_image: bool = False):
    st_ = prepare_movie_state(R, variant)
    T = len(R["gt"])
    frames = [movie_frame(R, variant, t, st_, show_image) for t in range(T)]
    labels = []
    for t in range(T):
        n_ev = len(st_["events"].get(t, []))
        n_lost = len(st_["lost"].get(t, []))
        labels.append(f"frame {t + 1}/{T}   " + (f"ERRORS: {n_ev} switch(es), {n_lost} missed" if n_ev or n_lost else "no errors"))
    return frames, labels


# ----------------------------------------------------------------------------- single cell and all-cells movies
def _norm(a, lo=None, hi=None):
    lo = np.percentile(a, 1) if lo is None else lo
    hi = np.percentile(a, 99.5) if hi is None else hi
    return np.clip((a.astype(float) - lo) / max(hi - lo, 1e-9), 0, 1)


def cell_composite(nuc: np.ndarray, loci: list, nuc_w: float = 0.35) -> np.ndarray:
    """(S,S) nucleus + list of (S,S) locus images -> RGB 0-255 (loci in red/green/blue)."""
    rgb = np.repeat(_norm(nuc)[..., None] * nuc_w, 3, axis=2)
    for k, im in enumerate(loci[:3]):
        im = im.astype(float)
        rgb[..., k] = np.clip(rgb[..., k] + np.clip((im - np.median(im)) / max(im.max() - np.median(im), 1e-9), 0, 1), 0, 1)
    return rgb * 255


def single_cell_frames(channels: dict, mask: np.ndarray, frames_df, loc_names: list, pos, truth, px_um: float,
                       target_w: int = 430):
    """Frames for one cell: [loci on nucleus (+ detected / true markers)] | [nucleus with mask outline]."""
    n, S = mask.shape[0], mask.shape[1]
    nuc = channels["nucleus"]
    nlo, nhi = float(np.percentile(nuc, 1)), float(np.percentile(nuc, 99.5))
    scale = target_w / S
    out, labels = [], []
    for i in range(n):
        t = int(frames_df.t.iloc[i])
        comp = cell_composite(nuc[i], [channels[c][i] for c in loc_names])
        a = Image.fromarray(comp.astype(np.uint8))
        d = ImageDraw.Draw(a)
        if pos is not None and len(pos):
            for q in pos[pos.t == t].itertuples():
                d.line([q.crop_x_px - 7, q.crop_y_px, q.crop_x_px + 7, q.crop_y_px], fill=(255, 255, 255), width=2)
                d.line([q.crop_x_px, q.crop_y_px - 7, q.crop_x_px, q.crop_y_px + 7], fill=(255, 255, 255), width=2)
        if truth is not None:
            fr = frames_df.iloc[i]
            for q in truth[truth.t == t].itertuples():
                x, y = q.x_um / px_um - fr.crop_x0, q.y_um / px_um - fr.crop_y0
                d.ellipse([x - 8, y - 8, x + 8, y + 8], outline=(0, 255, 255), width=1)
        b = outline(gray_to_rgb(nuc[i], nlo, nhi), mask[i], (255, 230, 0))
        out.append(hstack([titled(a, "Loci on nucleus  (+ found, o true)" if truth is not None else "Loci on nucleus  (+ found)", scale),
                           titled(Image.fromarray(b.astype(np.uint8)), "Nucleus + this cell's mask", scale)]))
        labels.append(f"frame {t}  ({i + 1}/{n} of this cell)")
    return out, labels


def all_cells_frames(cells_data: list, T: int, ncol: int = 4, tile: int = 150):
    """Montage of every cell side by side, one frame of the movie at a time. cells_data: [(cell_id, t->RGB tile dict)]"""
    n = len(cells_data)
    nrow = int(np.ceil(n / ncol))
    frames, labels = [], []
    for t in range(T):
        canvas = Image.new("RGB", (ncol * tile, nrow * (tile + 16)), (0, 0, 0))
        d = ImageDraw.Draw(canvas)
        for k, (cid, tiles) in enumerate(cells_data):
            x, y = (k % ncol) * tile, (k // ncol) * (tile + 16)
            if t in tiles:
                canvas.paste(Image.fromarray(tiles[t]), (x, y + 16))
            d.text((x + 4, y + 1), f"cell {cid}" + ("" if t in tiles else "  (absent)"), fill=(235, 235, 235), font=FONT)
        frames.append(np.asarray(canvas))
        labels.append(f"frame {t}")
    return frames, labels
