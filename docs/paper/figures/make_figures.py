r"""Draw the paper's figures from docs/paper/data/. Only stale figures are redrawn.

    scripts\run.cmd python docs\paper\figures\make_figures.py           # redraw what changed
    scripts\run.cmd python docs\paper\figures\make_figures.py --force   # redraw everything

A figure is stale when its SVG is missing or when the hash of its input file plus
this script differs from the one in .figure-cache.json. Nothing here runs a simulation.
Style: blue = simulation or analysis output, dashed grey = exact reference, orange = the
uncorrected alternative, pale bands = context. One y-axis per panel.
"""
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
CACHE = HERE / ".figure-cache.json"
BLUE, ORANGE, AQUA = "#257de0", "#f47230", "#00b894"
INK, MUTED, GRID, REF, BAND = "#191919", "#666666", "#e6e6e6", "#7a7a7a", "#ececec"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 9.5, "axes.titleweight": "bold",
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
                     "legend.frameon": False, "svg.fonttype": "none"})
FIGURES = {}


def figure(name, source=None):
    def wrap(fn):
        FIGURES[name] = (fn, source)
        return fn
    return wrap


def checks():
    return json.loads((DATA / "checks.json").read_text(encoding="utf-8"))


@figure("fig-pipeline")
def pipeline():
    fig, ax = plt.subplots(figsize=(7.4, 2.5))
    ax.axis("off"); ax.set_xlim(0, 7.4); ax.set_ylim(0, 2.5)
    top = [("1  Chromatin", "polymer library,\nprobes, transcription"), ("2  Cells", "moving, colliding,\ndividing nuclei"),
           ("3  Microscope", "optics, camera,\n28 error sources"), ("4  Segment\nand track", "Cellpose,\nTrackastra")]
    bottom = [("5  Isolate", "one small movie\nper cell"), ("6  Analyse", "spots, MS2,\nPol II starts"),
              ("7  Grade", "compare with\nthe known truth")]

    def box(x, y, head, body, colour):
        ax.add_patch(FancyBboxPatch((x, y), 1.55, 0.92, boxstyle="round,pad=0.02,rounding_size=0.08", fc="white", ec=colour, lw=1.4))
        ax.text(x + 0.775, y + 0.66, head, ha="center", va="center", fontsize=8.6, weight="bold", color=INK, linespacing=0.95)
        ax.text(x + 0.775, y + 0.27, body, ha="center", va="center", fontsize=7.2, color=MUTED, linespacing=1.05)

    for i, (h, b) in enumerate(top):
        x = 0.15 + i * 1.83
        box(x, 1.45, h, b, BLUE)
        if i < 3:
            ax.annotate("", (x + 1.81, 1.91), (x + 1.57, 1.91), arrowprops=dict(arrowstyle="->", color=MUTED))
    for i, (h, b) in enumerate(bottom):
        x = 5.64 - i * 1.83
        box(x, 0.12, h, b, AQUA)
        if i < 2:
            ax.annotate("", (x - 0.26, 0.58), (x - 0.02, 0.58), arrowprops=dict(arrowstyle="->", color=MUTED))
    ax.annotate("", (6.42, 1.06), (6.42, 1.43), arrowprops=dict(arrowstyle="->", color=MUTED))
    ax.text(0.15, 0.58, "The truth of every stage is saved,\nso stage 7 can grade each step.", fontsize=7.6, color=MUTED, va="center")
    return fig


@figure("val-probes", "checks.json")
def probes():
    r = checks()["probes"]["rows"]
    k = [x["k_off"] for x in r]
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    ax[0].plot(k, [x["mean_theory"] for x in r], "--", color=REF, label="exact theory")
    ax[0].plot(k, [x["mean_sim"] for x in r], "o", color=BLUE, ms=5, label="simulation")
    ax[0].set(xscale="log", xlabel="release rate k_off (per s)", ylabel="mean probes attached", title="Mean occupancy")
    ax[0].legend()
    ax[1].plot(k, [x["cv_theory"] for x in r], "--", color=REF)
    ax[1].plot(k, [x["cv_sim"] for x in r], "o", color=BLUE, ms=5)
    ax[1].set(xscale="log", xlabel="release rate k_off (per s)", ylabel="CV of probes attached", title="Fluctuations")
    for a in ax:
        a.set_xticks([0.03, 0.1, 0.3]); a.set_xticklabels(["0.03", "0.1", "0.3"]); a.minorticks_off()
    return fig


@figure("val-transcription", "checks.json")
def transcription():
    c = checks()["transcription"]
    t = [x / 60 for x in c["t_s"]]
    fig, ax = plt.subplots(2, 1, figsize=(7.0, 3.6), sharex=True)
    for a in ax:
        start = None
        for ti, on in zip(t + [t[-1]], c["contact"] + [False]):
            if on and start is None:
                start = ti
            if not on and start is not None:
                a.axvspan(start, ti, color=BAND, lw=0); start = None
    ax[0].plot(t, c["on_exact"], "--", color=REF, label="exact expectation")
    ax[0].plot(t, c["on_sim"], color=BLUE, lw=1.2, label="mean of 300 simulated cells")
    ax[0].set(ylabel="P(promoter ON)", title="Promoter state (grey bands: enhancer in contact)")
    ax[0].set_ylim(0, max(c["on_exact"]) * 1.45)
    ax[0].legend(ncol=2, loc="upper right", fontsize=7.5)
    ax[1].plot(t, c["ms2_exact"], "--", color=REF)
    ax[1].plot(t, c["ms2_sim"], color=BLUE, lw=1.2)
    ax[1].set(xlabel="time (min)", ylabel="MS2 loops on the gene", title="MS2 signal")
    return fig


@figure("val-ep-recovery", "checks.json")
def ep_recovery():
    c = checks()["ep_recovery"]
    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    ax.plot(c["d_nm"], c["true"], "--", color=REF, lw=1.6, label="curve put in")
    ax.plot(c["d_nm"], c["naive"], color=ORANGE, lw=1.4, label="fit ignoring measurement error")
    ax.plot(c["d_nm"], c["fit"], color=BLUE, lw=1.4, label="errors-in-variables fit")
    ax.set(xlabel="true enhancer-promoter distance (nm)", ylabel="P(burst ON)", title="Recovered curve")
    ax.legend(fontsize=7.5)
    return fig


def main():
    force = "--force" in sys.argv
    wanted = [a for a in sys.argv[1:] if not a.startswith("--")] or list(FIGURES)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    code = Path(__file__).read_bytes()
    for name in wanted:
        fn, source = FIGURES[name]
        h = hashlib.sha256(code + ((DATA / source).read_bytes() if source else b"")).hexdigest()[:16]
        out = HERE / f"{name}.svg"
        if not force and out.exists() and cache.get(name) == h:
            continue
        fig = fn()
        fig.tight_layout()
        fig.savefig(out, metadata={"Date": None})
        plt.close(fig)
        cache[name] = h
        print("drew", out.name)
    CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True) + "\n", newline="\n")


if __name__ == "__main__":
    main()
