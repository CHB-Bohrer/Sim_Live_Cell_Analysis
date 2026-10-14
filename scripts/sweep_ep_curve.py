"""How well is the curve P(burst ON | promoter-enhancer distance) recovered, as a function of the measurement errors? (docs/EP_CURVE.md)

    scripts\\run.cmd python scripts\\sweep_ep_curve.py --exps ladder sigma shift outliers cells null --workers 12
    scripts\\run.cmd python scripts\\sweep_ep_curve.py --exps ms2 --workers 12          (deconvolution: slower)

Experiments (each writes $SIMLIVE_DATA/sweeps/ep_curve_<exp>/results.csv in the sweep layout, so the Scans page shows them):
  ladder    L0 true distance + true ON | L1 true distance + ON estimated from the MS2 trace | L2 measured distance + true ON | L3 both estimated
  sigma     localization precision x 2D / 3D, EIV fit vs naive fit (L2)
  shift     residual colour shift, corrected / not (L2)
  outliers  wrong picks and missed frames (L2)
  cells     number of cells: confidence interval width and power (L0, L2)
  null      coupling switched off: the estimated contrast must be 1 and its 95% interval must exclude 1 about 5% of the time (L2, L3)
  dynamics  how long the pair stays in contact (L0, L2): the reference curve changes with it
  contact   contact fraction 0.03 / 0.1 / 0.2 (L2)
  ms2       MS2 noise x frame interval (L1, L3), with / without the attenuation correction
  rules     the burst rule (window, starts needed) against the bias of the recovered contrast (L3)

Baseline measurement (PLACEHOLDERS, edit here): difference-vector precision 57 nm per axis (= 40 nm per locus), distances in the 2D image plane,
no residual colour shift, 5% missed frames, 2% wrong picks. Baseline cells: 600 x 1 h at 10 s frames; the enhancer acts over a working distance
of 400 nm (hill n = 4), fold 12 on k_on, contact fraction 0.1.
Metrics: contrast = P(ON at 150 nm) / P(ON at 1000 nm); covers = the 95% cell-bootstrap interval contains the reference value (the same
fit to the TRUE distance and ON state of 4000 cells); curve_err = RMS difference of the whole curve, weighted by how often each distance occurs.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):     # one BLAS thread per worker process (12 workers x 16 threads was 4x slower)
    os.environ.setdefault(_v, "1")

import argparse
import itertools
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace

import numpy as np
import pandas as pd

from simlive.io.runs import DATA_ROOT
from simlive.stage7_validation import ep_curve as ep

BASE_M = ep.Measurement(sigma_nm=57.0, dims=2, shift_nm=0.0, miss_rate=0.05, outlier_rate=0.02)
BASE_S = ep.EPSetup()
N_BOOT = 100
_REFS: dict = {}


def reference(setup: ep.EPSetup) -> dict:
    key = (setup.working_distance_nm, setup.hill_n, setup.contact_fraction, setup.loop_time_s, setup.relax_s, setup.fold, setup.frame_s,
           setup.block_s, tuple(sorted(setup.tx.items())))
    if key not in _REFS:
        _REFS[key] = ep.reference_curve(replace(setup, n_cells=4000, frame_s=min(setup.frame_s, 10.0)))
    return _REFS[key]


def metrics(a: dict, ref: dict, tag: str = "") -> dict:
    f, r = a["fit"], ref["fit"]
    row = {"contrast": a["contrast"], "contrast_ref": r.contrast(), "contrast_ratio": a["contrast"] / r.contrast(), "fold": f.fold, "d_half": f.d_half,
           "d_half_ref": r.d_half, "p_far": f.p_far, "p_far_ref": r.p_far, "curve_err": ep.curve_error(f, ref), "frames_used": a["frames_used"],
           "n_on": a["n_on"]}
    if "contrast_ci" in a:
        lo, hi = a["contrast_ci"]
        row.update({"contrast_lo": lo, "contrast_hi": hi, "contrast_ci_width": hi - lo, "contrast_covers": float(lo <= r.contrast() <= hi),
                    "excludes_one": float(lo > 1.0 or hi < 1.0), "d_half_lo": a["d_half_ci"][0], "d_half_hi": a["d_half_ci"][1],
                    "d_half_covers": float(a["d_half_ci"][0] <= r.d_half <= a["d_half_ci"][1])})
    return row


def task(spec: dict) -> dict:
    """One population + one analysis. spec: setup, meas, seed, coupled, level ('L0' | 'L2'), naive, ref, settings (the varied columns), name."""
    t0 = time.time()
    pop = ep.simulate_population(spec["setup"], spec["seed"], coupled=spec["coupled"], with_trace=False)
    m = spec["meas"] if spec["level"] == "L2" else ep.Measurement(sigma_nm=0.5, dims=3)
    a = ep.analyse(pop, m, naive=spec["naive"], n_boot=N_BOOT, seed=spec["seed"])
    row = {**spec["settings"], "seed": spec["seed"], "run": spec["name"], "seconds": time.time() - t0, **metrics(a, spec["ref"])}
    return row


def run_pool(specs: list[dict], workers: int) -> list[dict]:
    if workers > 1:
        with ProcessPoolExecutor(workers) as ex:
            return list(ex.map(task, specs, chunksize=2))
    return [task(s) for s in specs]


def write(exp: str, rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    lead = [c for c in df.columns if c not in ("seed", "run", "seconds") and not c.startswith(("contrast", "fold", "d_half", "p_far", "curve_err", "frames_used", "n_on",
                                                                                              "excludes", "sensitivity", "false_pos", "n_cells_used"))]
    cols = lead + ["seed", "run", "seconds"] + [c for c in df.columns if c not in lead + ["seed", "run", "seconds"]]
    out = DATA_ROOT / "sweeps" / f"ep_curve_{exp}"
    out.mkdir(parents=True, exist_ok=True)
    df[cols].to_csv(out / "results.csv", index=False)
    print(f"wrote {out / 'results.csv'} ({len(df)} rows)", flush=True)
    return df


def specs_for(exp: str, seeds: int) -> list[dict]:
    S, M = BASE_S, BASE_M
    ref = reference(S)
    out = []

    def add(settings, setup=S, meas=M, level="L2", naive=False, coupled=True, n=seeds, seed0=1000):
        r = reference(setup)
        for k in range(n):
            out.append({"setup": setup, "meas": meas, "level": level, "naive": naive, "coupled": coupled, "ref": r, "settings": settings,
                        "seed": seed0 + k, "name": f"{exp}_{len(out)}"})

    if exp == "sigma":
        for sg, dims, naive in itertools.product((10, 25, 57, 100, 200, 400), (2, 3), (False, True)):
            add({"sigma_nm": sg, "dims": dims, "fit": "naive" if naive else "eiv"}, meas=replace(M, sigma_nm=float(sg), dims=dims), naive=naive)
    elif exp == "shift":
        for sh, corrected in itertools.product((0, 25, 50, 100, 200), (False, True)):
            add({"shift_nm": sh, "corrected": corrected}, meas=replace(M, shift_nm=float(sh), assumed_shift_nm=float(sh) if corrected else 0.0))
    elif exp == "outliers":
        for o in (0.0, 0.02, 0.05, 0.1, 0.2, 0.4):
            add({"outlier_rate": o, "miss_rate": M.miss_rate}, meas=replace(M, outlier_rate=o))
        for mr in (0.0, 0.3, 0.6, 0.8):
            add({"outlier_rate": M.outlier_rate, "miss_rate": mr}, meas=replace(M, miss_rate=mr))
    elif exp == "cells":
        for n_cells, lvl in itertools.product((50, 100, 200, 400, 800, 1600, 3200), ("L0", "L2")):
            add({"n_cells": n_cells, "level": lvl}, setup=replace(S, n_cells=n_cells), level=lvl)
    elif exp == "null":
        add({"coupling": "off"}, coupled=False, n=max(seeds, 60), seed0=5000)
        add({"coupling": "on"}, coupled=True, n=max(seeds, 60), seed0=6000)
    elif exp == "dynamics":
        for lt, lvl in itertools.product((60.0, 240.0, 960.0, 3600.0), ("L0", "L2")):
            add({"loop_time_s": lt, "level": lvl}, setup=replace(S, loop_time_s=lt), level=lvl)
    elif exp == "contact":
        for cf in (0.03, 0.1, 0.2):
            add({"contact_fraction": cf}, setup=replace(S, contact_fraction=cf))
    elif exp == "frames":
        for fs in (10.0, 20.0, 40.0, 80.0):
            add({"frame_s": fs}, setup=replace(S, frame_s=fs))
    else:
        raise ValueError(exp)
    return out


def run_ladder_and_ms2(exp: str, seeds: int, workers: int) -> list[dict]:
    """L0-L3 (and the MS2 noise x frame interval scan): needs the deconvolution, so cells run in parallel inside each population."""
    rows = []
    grid = [{"ms2_sb": 0.05, "frame_s": 10.0}] if exp in ("ladder", "rules") else [{"ms2_sb": sb, "frame_s": fs} for sb in (0.05, 0.4, 1.0, 2.0) for fs in (10.0, 20.0)]
    for g in grid:
        setup = replace(BASE_S, n_cells=600 if exp == "ladder" else 300, frame_s=g["frame_s"])
        ref = reference(setup)
        for k in range(seeds):
            t0 = time.time()
            seed = 2000 + k
            pop = ep.simulate_population(setup, seed, noise_sb=g["ms2_sb"], noise_beta=0.03)
            est = ep.estimate_on(pop, workers=workers, seed=seed)
            if exp == "rules":                      # the burst rule against the bias of the recovered contrast (L3 only)
                for w, gb in ((5, 2), (5, 1), (3, 3), (8, 3), (12, 4), (2, 2), (8, 6)):
                    p_on = est.p_on(w, gb)
                    sc = est.scores(p_on)
                    for corrected in (False, True):
                        a = ep.analyse(pop, BASE_M, y=p_on, correct_calls=(sc["sensitivity"], sc["false_positive"]) if corrected else None, n_boot=N_BOOT, seed=seed)
                        rows.append({"rule_w": w, "rule_gb": gb, "corrected": corrected, "seed": seed, "run": f"{exp}_{len(rows)}", "seconds": time.time() - t0,
                                     "sensitivity": sc["sensitivity"], "false_positive": sc["false_positive"], **metrics(a, ref)})
                print(f"  rules seed {seed} ({time.time() - t0:.0f}s)", flush=True)
                continue
            w, gb = est.best_rule()
            p_on = est.p_on(w, gb)
            sc = est.scores(p_on)
            for level, y, meas in (("L0", None, ep.Measurement(sigma_nm=0.5, dims=3)), ("L1", p_on, ep.Measurement(sigma_nm=0.5, dims=3)),
                                   ("L2", None, BASE_M), ("L3", p_on, BASE_M)):
                if exp == "ms2" and level in ("L0", "L2"):
                    continue
                for corrected in ((False, True) if y is not None else (False,)):
                    a = ep.analyse(pop, meas, y=y, correct_calls=(sc["sensitivity"], sc["false_positive"]) if corrected else None, n_boot=N_BOOT, seed=seed)
                    rows.append({"level": level, "corrected": corrected, **g, "seed": seed, "run": f"{exp}_{len(rows)}", "seconds": time.time() - t0,
                                 "rule_w": w, "rule_gb": gb, "sensitivity": sc["sensitivity"], "false_positive": sc["false_positive"],
                                 "noise_sb_est": est.noise_sb, **metrics(a, ref)})
            print(f"  {exp} {g} seed {seed}: rule ({w:g}, {gb:g}) sens {sc['sensitivity']:.2f} fp {sc['false_positive']:.3f} ({time.time() - t0:.0f}s)", flush=True)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exps", nargs="+", default=["sigma", "shift", "outliers", "cells", "null", "dynamics", "contact", "frames"])
    ap.add_argument("--seeds", type=int, default=20, help="populations per setting (fast experiments)")
    ap.add_argument("--deconv-seeds", type=int, default=4, help="populations per setting for ladder / ms2 (deconvolution)")
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    for exp in a.exps:
        t0 = time.time()
        print(f"== {exp}", flush=True)
        rows = run_ladder_and_ms2(exp, a.deconv_seeds, a.workers) if exp in ("ladder", "ms2", "rules") else run_pool(specs_for(exp, a.seeds), a.workers)
        df = write(exp, rows)
        print(f"   {time.time() - t0:.0f} s", flush=True)
        key = [c for c in df.columns if c in ("level", "sigma_nm", "dims", "fit", "shift_nm", "corrected", "outlier_rate", "miss_rate", "n_cells", "coupling",
                                              "loop_time_s", "contact_fraction", "frame_s", "ms2_sb", "rule_w", "rule_gb")]
        agg = df.groupby(key).agg(contrast=("contrast", "mean"), ref=("contrast_ref", "mean"), sd=("contrast", "std"), curve_err=("curve_err", "mean"),
                                  cover=("contrast_covers", "mean"), excl1=("excludes_one", "mean"), ci_w=("contrast_ci_width", "mean"))
        print(agg.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
