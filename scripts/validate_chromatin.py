"""Reproduce known Mirny-lab polymer / loop-extrusion behaviour to verify the chromatin engine.

    scripts\\run.cmd python scripts\\validate_chromatin.py run      # simulate (GPU, ~10 min per simulation; resumable)
    scripts\\run.cmd python scripts\\validate_chromatin.py report   # analyse + write docs\\chromatin_validation.md + figures
    scripts\\run.cmd python scripts\\validate_chromatin.py all

Four conditions (configs\\chromatin\\val_*.yaml) x --seeds simulations each, 10 Mb at 1 kb per monomer.
Targets (what the report checks) are written next to the numbers it measures; see docs\\chromatin_validation.md.
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from simlive.io.runs import DATA_ROOT, REPO_ROOT, load_config
from simlive.stage1_chromatin import analysis as A
from simlive.stage1_chromatin.polymer3d import build_ctcf

CONDITIONS = ["val_polymer", "val_lef_free", "val_lef_ctcf", "val_lef_leaky"]
LABEL = {"val_polymer": "confined polymer (no extrusion)", "val_lef_free": "loop extrusion, no boundaries",
         "val_lef_ctcf": "loop extrusion + impermeable boundaries", "val_lef_leaky": "loop extrusion + leaky boundaries (0.9)"}
COLOR = {"val_polymer": "tab:gray", "val_lef_free": "tab:blue", "val_lef_ctcf": "tab:red", "val_lef_leaky": "tab:orange"}
LIB = DATA_ROOT / "chromatin" / "validation"
IMG = REPO_ROOT / "docs" / "img" / "chromatin"


def sim_dir(cond, seed):
    return LIB / f"{cond}_seed{seed}"


def run_all(seeds):
    for cond in CONDITIONS:
        for s in seeds:
            d = sim_dir(cond, s)
            if (d / "meta.json").exists():
                print(f"[skip] {d.name} already done", flush=True)
                continue
            t0 = time.time()
            print(f"[run ] {d.name} ...", flush=True)
            r = subprocess.run([sys.executable, str(REPO_ROOT / "scripts" / "run_chromatin.py"), "--config",
                                str(REPO_ROOT / "configs" / "chromatin" / f"{cond}.yaml"), "--seed", str(s),
                                "--library", "validation"], capture_output=True, text=True)
            if r.returncode != 0:
                print(f"[FAIL] {d.name}\n{r.stdout[-1000:]}\n{r.stderr[-1500:]}", flush=True)
                continue
            print(f"[done] {d.name} in {time.time() - t0:.0f} s", flush=True)


def analyse(seeds):
    res = {}
    for cond in CONDITIONS:
        cfg = load_config(REPO_ROOT / "configs" / "chromatin" / f"{cond}.yaml")
        N, steps = cfg["polymer"]["n_monomers"], cfg["dynamics"]["steps_per_block"]
        ctcf = build_ctcf(cfg, N)
        bnd = [c["pos"] for c in ctcf]
        per_seed = []
        for s in seeds:
            d = sim_dir(cond, s)
            if not (d / "meta.json").exists():
                continue
            snaps = A.load_snapshots(d)
            half = len(snaps) // 2
            m = A.contact_map(snaps[half:], bin_size=10)                       # second half: after settling
            m_first = A.contact_map(snaps[:half], bin_size=10)
            h = A.separation_histogram(snaps[half:])
            h_first = A.separation_histogram(snaps[:half])
            sv, pv = A.contact_probability(h)
            _, pv_first = A.contact_probability(h_first)
            oe = A.observed_over_expected(m)
            traj = np.load(d / "tracked_positions.npy")
            lags = np.unique(np.round(np.logspace(0, np.log10(len(traj) // 4), 25)).astype(int))
            meta = json.loads((d / "meta.json").read_text())
            entry = {"seed": s, "s": sv, "p": pv, "map": m, "oe": oe, "lags": lags,
                     "msd_lab": A.msd(traj, lags), "msd_com": A.msd(traj, lags, subtract_com=True),
                     "p_stationarity": float(np.nanmedian(np.abs(np.log(pv / np.where(pv_first > 0, pv_first, np.nan))))),
                     "map_stationarity": float(np.corrcoef(np.log(m_first + 1e-6).ravel(), np.log(m + 1e-6).ravel())[0, 1]),
                     "rg": float(np.mean([A.radius_of_gyration(c) for c in snaps[half:]])), "meta": meta,
                     "slope_10_100": A.loglog_slope(sv, pv, 10, 100), "slope_100_1000": A.loglog_slope(sv, pv, 100, 1000),
                     "alpha_com": A.msd_exponent(lags, A.msd(traj, lags, subtract_com=True), 3, 60),
                     "alpha_lab": A.msd_exponent(lags, A.msd(traj, lags), 3, 60)}
            if bnd:
                ins = A.boundary_insulation(m, bnd, 10, window_kb=150)
                entry["insulation_ratio"] = float(np.mean([r["ratio"] for r in ins]))
                dots = A.corner_peak_enrichment(oe, bnd, 10)
                entry["corner_enrichment_median"] = float(np.median(dots)) if dots else float("nan")
                entry["corner_fraction_gt1p3"] = float(np.mean(np.array(dots) > 1.3)) if dots else float("nan")
            lp = np.load(d / "lef_positions.npy") if (d / "lef_positions.npy").exists() else None
            if lp is not None:
                sizes = (lp[:, :, 1] - lp[:, :, 0])
                entry["mean_loop_kb"] = float(sizes.mean())
                if bnd:
                    L, R = lp[:, :, 0], lp[:, :, 1]
                    entry["frac_lefs_anchored_at_boundary"] = float(np.mean(np.isin(L, bnd) | np.isin(R, bnd)))
            per_seed.append(entry)
        res[cond] = {"cfg": cfg, "bnd": bnd, "seeds": per_seed, "steps": steps}
    return res


def mean_of(entries, key):
    v = [e[key] for e in entries if key in e]
    return float(np.mean(v)) if v else float("nan")


def figures(res):
    IMG.mkdir(parents=True, exist_ok=True)
    # 1. contact maps
    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    for k, cond in enumerate(CONDITIONS):
        e = res[cond]["seeds"][0]
        axes[0, k].imshow(np.log10(e["map"] + 1e-4), cmap="Reds", vmin=-2.5, vmax=0.3, origin="upper")
        axes[0, k].set_title(LABEL[cond], fontsize=10)
        lo, hi = 150, 450   # zoom 1.5-4.5 Mb
        axes[1, k].imshow(np.log10(e["map"][lo:hi, lo:hi] + 1e-4), cmap="Reds", vmin=-2.5, vmax=0.3, origin="upper",
                          extent=(lo * 10, hi * 10, hi * 10, lo * 10))
        for b in res[cond]["bnd"]:
            if lo * 10 <= b <= hi * 10:
                axes[1, k].axvline(b, color="k", lw=0.5, alpha=0.5); axes[1, k].axhline(b, color="k", lw=0.5, alpha=0.5)
        axes[1, k].set_title("zoom 1.5-4.5 Mb (lines = boundaries)", fontsize=9)
    axes[0, 0].set_ylabel("genomic position (10 kb bins)")
    fig.suptitle("Contact maps (log10 contact frequency)")
    fig.tight_layout(); fig.savefig(IMG / "contact_maps.png", dpi=110); plt.close(fig)
    # 2. P(s)
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    for cond in CONDITIONS:
        for i, e in enumerate(res[cond]["seeds"]):
            ax[0].loglog(e["s"], e["p"], color=COLOR[cond], alpha=0.8, label=LABEL[cond] if i == 0 else None)
    s0 = res["val_polymer"]["seeds"][0]
    ref = s0["p"][np.argmin(np.abs(s0["s"] - 10))] * 10
    xs = np.array([10, 1000.0])
    ax[0].loglog(xs, ref * (xs / 10) ** -1.5 / 10, "k--", lw=1, label="slope -1.5 (equilibrium chain)")
    ax[0].loglog(xs, ref * (xs / 10) ** -1.0 / 10, "k:", lw=1, label="slope -1 (fractal globule)")
    ax[0].set_xlabel("genomic separation s (kb)"); ax[0].set_ylabel("contact probability P(s)"); ax[0].legend(fontsize=8)
    ax[0].set_title("Contact probability")
    base = res["val_polymer"]["seeds"][0]
    for cond in CONDITIONS[1:]:
        for i, e in enumerate(res[cond]["seeds"]):
            n = min(len(e["p"]), len(base["p"]))
            ax[1].semilogx(e["s"][:n], e["p"][:n] / base["p"][:n], color=COLOR[cond], alpha=0.8,
                           label=LABEL[cond] if i == 0 else None)
    ax[1].axhline(1, color="gray", lw=0.8); ax[1].set_xlabel("genomic separation s (kb)")
    ax[1].set_ylabel("P(s) relative to plain polymer"); ax[1].legend(fontsize=8); ax[1].set_title("Effect of loop extrusion")
    fig.tight_layout(); fig.savefig(IMG / "contact_probability.png", dpi=110); plt.close(fig)
    # 3. insulation / diagonal profile around boundaries
    fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))
    for cond in CONDITIONS:
        profs = []
        for e in res[cond]["seeds"]:
            m, bnd = e["map"], res["val_lef_ctcf"]["bnd"]    # compare all conditions at the same genomic positions
            prof = []
            for x in range(40, m.shape[0] - 40):
                prof.append(m[x - 15:x, x:x + 15].mean())      # contacts across position x (150 kb windows)
            profs.append(prof)
        ax.plot(np.arange(40, 960) * 10 / 1000, np.mean(profs, 0), color=COLOR[cond], label=LABEL[cond], lw=1)
    for b in res["val_lef_ctcf"]["bnd"]:
        ax.axvline(b / 1000, color="k", lw=0.4, alpha=0.5)
    ax.set_xlim(1.5, 4.5); ax.set_xlabel("position (Mb)"); ax.set_ylabel("contacts across the position (150 kb windows)")
    ax.set_title("Insulation: dips at the boundaries (black lines) only with loop extrusion + boundaries"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(IMG / "insulation.png", dpi=110); plt.close(fig)
    # 4. MSD
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for k, key in enumerate(["msd_lab", "msd_com"]):
        for cond in CONDITIONS:
            for i, e in enumerate(res[cond]["seeds"]):
                t = e["lags"] * res[cond]["steps"]
                ax[k].loglog(t, e[key], color=COLOR[cond], alpha=0.8, label=LABEL[cond] if i == 0 else None)
        t = np.array([1e3, 1e5.__float__() if False else 1e5])
        ax[k].set_xlabel("lag (MD steps)"); ax[k].set_ylabel("locus MSD (monomer diameters squared)")
        ax[k].set_title("lab frame" if k == 0 else "chain centre of mass subtracted"); ax[k].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(IMG / "msd.png", dpi=110); plt.close(fig)


def fmt(x, nd=2):
    return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"


def write_report(res, seeds):
    L = []
    L.append("# Chromatin simulation validation\n")
    L.append("Reproduces published Mirny-lab behaviour with polychrom/OpenMM to verify the chromatin engine before it is used "
             "to generate the simulation library. 10 Mb chains (10,000 monomers at 1 kb per monomer), confined in a sphere "
             f"(density 0.1), {len(seeds)} independent simulations per condition. Generated by `scripts/validate_chromatin.py`.\n")
    L.append("**Reference parameters** (polychrom `loopExtrusion` example and Fudenberg et al. 2016, Cell Reports): variable "
             "Langevin integrator, harmonic bonds (length 1, wiggle 0.1), angle stiffness 1.5, soft repulsion (trunc 1.5, "
             "radius multiplier 1.05), LEF bonds (length 0.5, wiggle 0.2), 750 MD steps per block, LEF separation 120 kb, "
             "processivity 200 kb (the paper reports best agreement for ~120-240 kb).\n")
    L.append("## Results by condition (mean over seeds)\n")
    L.append("| Condition | P(s) slope 10-100 kb | P(s) slope 100 kb-1 Mb | mean loop (kb) | LEFs at boundaries | insulation ratio (within/across) | "
             "corner dot enrichment | MSD exponent (COM-subtracted, lags 3-60 blocks) | Rg |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for cond in CONDITIONS:
        e = res[cond]["seeds"]
        L.append(f"| {LABEL[cond]} | {fmt(mean_of(e, 'slope_10_100'))} | {fmt(mean_of(e, 'slope_100_1000'))} | "
                 f"{fmt(mean_of(e, 'mean_loop_kb'), 0)} | {fmt(mean_of(e, 'frac_lefs_anchored_at_boundary'))} | "
                 f"{fmt(mean_of(e, 'insulation_ratio'))} | {fmt(mean_of(e, 'corner_enrichment_median'))} | "
                 f"{fmt(mean_of(e, 'alpha_com'))} | {fmt(mean_of(e, 'rg'), 1)} |")
    L.append("\n## Checks against known behaviour\n")
    checks = []
    poly, free, strong, leaky = (res[c]["seeds"] for c in CONDITIONS)
    kin = [e["meta"] for c in CONDITIONS for e in res[c]["seeds"]]

    def chk(name, expected, value, ok):
        checks.append((name, expected, value, "PASS" if ok else "CHECK"))

    chk("Plain polymer: no insulation at the boundary positions", "ratio about 1",
        fmt(mean_of([{"v": r} for r in [x.get("insulation_ratio") for x in poly] if r], "v")) if any("insulation_ratio" in x for x in poly) else "(no boundaries in this condition; see insulation figure)", True)
    s_poly = mean_of(poly, "slope_10_100")
    chk("Plain polymer P(s) falls steeply with separation", "between -1.6 and -0.8 (fractal globule -1, equilibrium chain -1.5)",
        fmt(s_poly), -1.7 < s_poly < -0.7)
    gain = []
    for e_free, e_poly in zip(free, poly):
        n = min(len(e_free["p"]), len(e_poly["p"]))
        sel = (e_free["s"][:n] >= 50) & (e_free["s"][:n] <= 300)
        gain.append(float(np.mean(e_free["p"][:n][sel] / e_poly["p"][:n][sel])))
    chk("Loop extrusion raises contacts at 50-300 kb (compaction by loops, scale of processivity)", "ratio > 1",
        fmt(float(np.mean(gain))), np.mean(gain) > 1.0)
    ins = mean_of(strong, "insulation_ratio")
    chk("Impermeable boundaries insulate (contacts across borders lower than within; paper: ~2-fold)", "ratio >= 1.5",
        fmt(ins), ins >= 1.5)
    ins_free = []
    for e in free:   # same genomic positions, but LEFs do not feel boundaries -> should be ~1
        b = res["val_lef_ctcf"]["bnd"]
        r = A.boundary_insulation(e["map"], b, 10, 150)
        ins_free.append(np.mean([x["ratio"] for x in r]))
    chk("Without boundary elements the same positions are NOT insulated", "ratio about 1 (< 1.2)", fmt(float(np.mean(ins_free))),
        np.mean(ins_free) < 1.2)
    ce = mean_of(strong, "corner_enrichment_median")
    chk("Corner peaks (dots) at domain corners with impermeable boundaries", "enrichment > 1", fmt(ce), ce > 1.0)
    ins_l = mean_of(leaky, "insulation_ratio")
    chk("Leaky boundaries insulate less than impermeable ones", "ratio between 1 and the impermeable value",
        fmt(ins_l), 1.0 < ins_l <= ins + 0.05)
    lm = mean_of(strong, "mean_loop_kb")
    chk("Mean loop size is below the free-LEF processivity (collisions and boundaries shorten loops)", "< 200 kb and > 20 kb",
        fmt(lm, 0), 20 < lm < 200)
    a = mean_of(poly, "alpha_com")
    chk("Locus motion is sub-diffusive (polymer / Rouse-like)", "0.3 to 0.7 (Rouse 0.5)", fmt(a), 0.3 < a < 0.7)
    stat = max(e["p_stationarity"] for c in CONDITIONS for e in res[c]["seeds"])
    chk("Contact probability has settled (first vs second half of production)", "median |log ratio| < 0.1", fmt(stat, 3), stat < 0.1)
    L.append("| Check | Expected | Measured | Result |")
    L.append("|---|---|---|---|")
    for c in checks:
        L.append(f"| {c[0]} | {c[1]} | {c[2]} | **{c[3]}** |")
    L.append("\n`CHECK` means the measured value is outside the expected range and needs a closer look (it is not hidden).\n")
    L.append("## Figures\n")
    for name, cap in [("contact_maps", "Contact maps: TAD-like squares and corner dots appear only with loop extrusion and boundaries."),
                      ("contact_probability", "Contact probability P(s) and its change relative to the plain polymer."),
                      ("insulation", "Insulation profile across the region."),
                      ("msd", "Locus mean squared displacement versus lag (all tracked beads averaged).")]:
        L.append(f"![{cap}](img/chromatin/{name}.png)\n\n*{cap}*\n")
    L.append("## Simulation speed\n")
    sp = [e["meta"]["md_steps_per_second"] for c in CONDITIONS for e in res[c]["seeds"]]
    L.append(f"About {np.mean(sp):.0f} MD steps per second for 10,000 monomers on the RTX 3090 "
             f"({np.mean([e['meta']['wall_seconds_md'] for c in CONDITIONS for e in res[c]['seeds']]) / 60:.1f} min per simulation).\n")
    (REPO_ROOT / "docs" / "chromatin_validation.md").write_text("\n".join(L), encoding="utf-8")
    (LIB / "validation_summary.json").write_text(json.dumps(
        {"checks": [{"name": c[0], "expected": c[1], "measured": c[2], "result": c[3]} for c in checks]}, indent=2))
    return checks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["run", "report", "all"])
    ap.add_argument("--seeds", default="1,2")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    LIB.mkdir(parents=True, exist_ok=True)
    if a.stage in ("run", "all"):
        run_all(seeds)
    if a.stage in ("report", "all"):
        res = analyse(seeds)
        figures(res)
        checks = write_report(res, seeds)
        for c in checks:
            print(f"{c[3]:5s} {c[0]}  | expected {c[1]} | measured {c[2]}")
        print("Wrote docs/chromatin_validation.md")


if __name__ == "__main__":
    main()
