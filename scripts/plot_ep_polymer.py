"""One-page picture of the E-P curve study with polymer distances (docs/EP_CURVE.md, update 2026-10-14 b).

    scripts\\run.cmd python scripts\\plot_ep_polymer.py            -> docs/img/ep_polymer_summary.png

A: true E-P distance of one cell for 4 minutes (the polymer, 55 nm per unit, 0.46 s per block) with the 400 nm working distance and the 10 s frames.
B: P(burst ON) against the distance NOW (nearly flat) and against a distance HELD (the curve the promoter really responds to), true distance and true ON state of 3000 cells.
C: the held-distance curve recovered from 600 cells through a 2D, a z-stack and a biplane microscope (95% bands), against the truth.
D: the contrast P(ON at 150 nm) / P(ON at 1000 nm) of 20 populations per microscope (scans of scripts/sweep_ep_polymer.py, experiment `modality`).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from simlive.io.runs import DATA_ROOT
from simlive.stage7_validation import ep_acquisition as aq
from simlive.stage7_validation import ep_curve as ep
from simlive.stage7_validation import ep_history as eh
from simlive.stage7_validation import ep_polymer as pp

OUT = Path(__file__).resolve().parents[1] / "docs" / "img" / "ep_polymer_summary.png"
COL = {"2d": "#d97706", "zstack": "#2563eb", "biplane": "#0d9488", "truth": "#111827"}
NAME = {"2d": "2D widefield", "zstack": "3D z-stack", "biplane": "3D biplane"}


def main() -> None:
    pool = pp.load_pool()
    cal = pp.calibrate(pool)
    setup = ep.EPSetup()
    ref_pop = pp.simulate_population(ep.EPSetup(n_cells=3000), pool, cal, 12345, with_trace=False, keep_blocks=False)
    K0 = aq.kernel_aniso(9.0, 9.0, 3, 0.0)
    ref = eh.analyse_history(ref_pop, ref_pop.d_true, K0)
    g = ep.D_GRID

    fig = plt.figure(figsize=(14, 9.2))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.22)

    # A ------------------------------------------------------------------ one cell's distance
    ax = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(4)
    v = pp.cell_vectors(pool, cal, int(240 / cal.block_s), rng)
    d = np.linalg.norm(v, axis=1)
    t = np.arange(len(d)) * cal.block_s
    ax.plot(t, d, color="#6b7280", lw=0.9, label="true E-P distance (polymer)")
    fr = np.arange(0, 240, 10.0)
    ax.plot(fr, d[np.minimum(np.round(fr / cal.block_s).astype(int), len(d) - 1)], "o", color=COL["biplane"], ms=4, label="one frame every 10 s")
    ax.axhline(400, color="#dc2626", ls="--", lw=1.2, label="400 nm working distance")
    ax.set_xlabel("time (s)"); ax.set_ylabel("E-P distance (nm)"); ax.set_ylim(0, 2200); ax.set_xlim(0, 240)
    ax.set_title("A  The distance changes within ~9 s", loc="left", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right"); ax.grid(alpha=0.25)

    # B ------------------------------------------------------------------ now vs held
    ax = fig.add_subplot(gs[0, 1])
    T = setup.n_frames
    dd = ref_pop.d_true.reshape(-1, T)
    on = ref_pop.on_true.reshape(-1, T).astype(float)
    edges = np.array([0, 150, 250, 350, 450, 600, 800, 1100, 1600, 3000.0])
    b = np.digitize(dd, edges) - 1
    centres = 0.5 * (edges[:-1] + edges[1:])
    p_now = np.array([on[b == i].mean() if (b == i).any() else np.nan for i in range(len(edges) - 1)])
    ax.plot(centres, p_now, "s-", color="#9ca3af", lw=1.8, label=f"distance NOW (contrast {np.interp(150, centres, p_now) / np.interp(1000, centres, p_now):.1f})")
    ax.plot(g, ref["fit"].g(g), color=COL["truth"], lw=2.4, label=f"distance HELD for ~{ref['fit'].tau_s:.0f} s (contrast {ref['contrast']:.1f})")
    ax.axvline(400, color="#dc2626", ls="--", lw=1.0)
    ax.set_xlim(0, 1500); ax.set_ylim(0, None)
    ax.set_xlabel("promoter-enhancer distance (nm)"); ax.set_ylabel("P(promoter in a burst)")
    ax.set_title("B  The promoter integrates the distance", loc="left", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9); ax.grid(alpha=0.25)

    # C ------------------------------------------------------------------ recovered curves
    ax = fig.add_subplot(gs[1, 0])
    pop = pp.simulate_population(setup, pool, cal, 21, with_trace=False)
    ax.plot(g, ref["fit"].g(g), color=COL["truth"], lw=2.6, label="truth (3000 cells, true distance and state)", zorder=5)
    for kind in ("2d", "zstack", "biplane"):
        m = aq.measure(pop, pop.blocks, aq.Acquisition(kind=kind), seed=21)
        K = aq.kernel_aniso(m.sigma_xy, m.sigma_z, m.dims, m.shift_nm)
        res = eh.analyse_history(pop, m.d, K, n_boot=80, seed=21)
        band = np.array([bb.g(g) for bb in res["boot"]])
        ax.fill_between(g, np.quantile(band, 0.025, 0), np.quantile(band, 0.975, 0), color=COL[kind], alpha=0.16, lw=0)
        ax.plot(g, res["fit"].g(g), color=COL[kind], lw=1.8, label=f"{NAME[kind]}: contrast {res['contrast']:.1f} (95% {res['contrast_ci'][0]:.1f}-{res['contrast_ci'][1]:.1f})")
    ax.set_xlim(0, 1500); ax.set_ylim(0, None)
    ax.set_xlabel("promoter-enhancer distance held at (nm)"); ax.set_ylabel("P(promoter in a burst)")
    ax.set_title("C  Recovered from 600 cells x 1 h (one population each)", loc="left", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8.5, loc="upper right"); ax.grid(alpha=0.25)

    # D ------------------------------------------------------------------ 20 populations
    ax = fig.add_subplot(gs[1, 1])
    f = DATA_ROOT / "sweeps" / "ep_polymer_modality" / "results.csv"
    if f.exists():
        df = pd.read_csv(f)
        order = [("truth (L0)", "truth", "true distance\nand ON state"), ("2d", "2d", "2D\nwidefield"), ("zstack", "zstack", "3D\nz-stack"), ("biplane", "biplane", "3D\nbiplane")]
        for i, (key, ck, label) in enumerate(order):
            x = df[df.modality == key].contrast.to_numpy()
            ax.scatter(np.full(len(x), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(x)), x, s=22, color=COL[ck], alpha=0.65, edgecolor="none")
            ax.errorbar(i + 0.32, x.mean(), yerr=x.std(), fmt="o", color="black", capsize=4, ms=5)
        ax.axhline(ref["contrast"], color="#111827", ls="--", lw=1.4)
        ax.text(-0.45, ref["contrast"] + 0.12, f"truth {ref['contrast']:.2f}", ha="left", fontsize=9)
        ax.axhline(1.0, color="#9ca3af", ls=":", lw=1)
        ax.set_xticks(range(4)); ax.set_xticklabels([o[2] for o in order])
        ax.set_ylabel("contrast: P(ON at 150 nm) / P(ON at 1000 nm)")
    ax.set_title("D  20 populations of 600 cells each (mean +- sd in black)", loc="left", fontsize=12, fontweight="bold")
    ax.grid(alpha=0.25, axis="y")
    fig.suptitle("E-P distance and P(burst ON): polymer distances (55 nm per unit, 0.46 s per block), contact fraction 0.10, three microscopes", fontsize=13, y=0.995)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130, bbox_inches="tight", facecolor="white")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
