"""Interactive viewer for a simulated microscopy run (stage 3 output).

    scripts\\run.cmd python scripts\\view_run.py data\\runs\\<run_id>
    scripts\\run.cmd python scripts\\view_run.py data\\runs\\<run_id> --save movie.gif   # no window

Reads <run>/stage3_microscopy/{locus.tif, nucleus.tif, labels.tif}. Each is (T, Y, X) or (T, Z, Y, X);
z-stacks are max-projected (labels: the z-slice through the middle). Missing files are skipped.
Controls: slider = frame, keys left/right = step, space = play/pause, 'l' = toggle label overlay,
'i' = toggle cell-ID text.
"""
import argparse
from pathlib import Path

import matplotlib
import numpy as np
import tifffile


def load_stack(path: Path, is_label: bool = False):
    if not path.exists():
        return None
    a = tifffile.imread(path)
    if a.ndim == 4:  # T,Z,Y,X
        a = a[:, a.shape[1] // 2] if is_label else a.max(axis=1)
    if a.ndim != 3:
        raise ValueError(f"{path}: expected (T,Y,X) or (T,Z,Y,X), got {a.shape}")
    return a


def label_overlay(lab: np.ndarray) -> np.ndarray:
    """RGBA image, one stable color per cell ID (so persistent IDs keep their color)."""
    rgba = np.zeros(lab.shape + (4,), float)
    ids = np.unique(lab)
    ids = ids[ids > 0]
    cmap = matplotlib.colormaps["tab20"]
    for i in ids:
        rgba[lab == i] = (*cmap(int(i) % 20)[:3], 0.18)
    return rgba


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path, help="run directory, e.g. data/runs/demo")
    ap.add_argument("--save", type=Path, help="write .gif/.mp4/.png (first frame) instead of opening a window")
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--labels", type=Path,
                    help="label stack to overlay, relative to the run dir (default stage3_microscopy/labels.tif = "
                         "ground truth), e.g. stage4_segtrack/tracked_auto_masks.tif")
    args = ap.parse_args()
    if args.save:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from matplotlib.widgets import Slider

    d = args.run / "stage3_microscopy"
    loc, nuc, lab = (load_stack(d / "locus.tif"), load_stack(d / "nucleus.tif"),
                     load_stack(args.run / args.labels if args.labels else d / "labels.tif", True))
    panels = [(n, s) for n, s in (("locus", loc), ("nucleus", nuc)) if s is not None]
    show_ids_default = args.labels is not None
    if not panels:
        raise SystemExit(f"No locus.tif / nucleus.tif found in {d}")
    T = panels[0][1].shape[0]
    show = {"labels": lab is not None, "ids": show_ids_default}  # tracked labels: show IDs by default

    fig, axes = plt.subplots(1, len(panels), figsize=(5.5 * len(panels), 5.8), squeeze=False)
    axes = axes[0]
    ims, ovs, texts = [], [], [[] for _ in panels]
    for ax, (name, s) in zip(axes, panels):
        ims.append(ax.imshow(s[0], cmap="gray", vmin=np.percentile(s, 1), vmax=np.percentile(s, 99.8)))
        ovs.append(ax.imshow(label_overlay(lab[0]), interpolation="nearest") if lab is not None else None)
        ax.set_title(name)
        ax.axis("off")
    title = fig.suptitle("")

    def draw(t):
        t = int(t)
        for k, ((name, s), im) in enumerate(zip(panels, ims)):
            im.set_data(s[t])
            if ovs[k] is not None:
                ovs[k].set_data(label_overlay(lab[t]))
                ovs[k].set_visible(show["labels"])
            for tx in texts[k]:
                tx.remove()
            texts[k] = []
            if lab is not None and show["ids"]:
                for i in np.unique(lab[t]):
                    if i > 0:
                        yy, xx = np.nonzero(lab[t] == i)
                        texts[k].append(axes[k].text(xx.mean(), yy.mean(), str(i), color="yellow",
                                                     ha="center", va="center", fontsize=9))
        title.set_text(f"{args.run.name}   frame {t + 1}/{T}")

    if args.save:
        if args.save.suffix.lower() == ".png":
            draw(0)
            fig.savefig(args.save, dpi=150)
        else:
            anim = FuncAnimation(fig, lambda t: draw(t), frames=T, interval=1000 / args.fps)
            anim.save(args.save, fps=args.fps, dpi=100)
        print(f"Saved {args.save}")
        return

    fig.subplots_adjust(bottom=0.12, top=0.9)
    slider = Slider(fig.add_axes([0.2, 0.03, 0.6, 0.03]), "frame", 0, T - 1, valinit=0, valstep=1)
    slider.on_changed(lambda v: (draw(v), fig.canvas.draw_idle()))
    state = {"playing": False}
    timer = fig.canvas.new_timer(interval=1000 / args.fps)
    timer.add_callback(lambda: slider.set_val((slider.val + 1) % T) if state["playing"] else None)
    timer.start()

    def on_key(e):
        if e.key == "right":
            slider.set_val(min(slider.val + 1, T - 1))
        elif e.key == "left":
            slider.set_val(max(slider.val - 1, 0))
        elif e.key == " ":
            state["playing"] = not state["playing"]
        elif e.key == "l" and lab is not None:
            show["labels"] = not show["labels"]
            draw(slider.val); fig.canvas.draw_idle()
        elif e.key == "i" and lab is not None:
            show["ids"] = not show["ids"]
            draw(slider.val); fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", on_key)
    draw(0)
    plt.show()


if __name__ == "__main__":
    main()
