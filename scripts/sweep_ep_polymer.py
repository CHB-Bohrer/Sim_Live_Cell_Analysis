"""E-P curve scans with E-P distances from the polymer simulations, three microscope modes, and a comparison with the invented distance process (docs/EP_CURVE.md).

    scripts\\run.cmd python scripts\\sweep_ep_polymer.py --exps modality null wrong cells contact frames colour photons focus timescale compare_sources --workers 12
    scripts\\run.cmd python scripts\\sweep_ep_polymer.py --exps ms2 --workers 12                 (deconvolution: slower)

Distances: 10 validation simulations, pairs 100 kb apart, 55.4 nm per polymer unit (the pair is within 400 nm 10% of the time) and 0.46 s per saved block (single-locus MSD at 60 s
= 6 D t with D = 0.002 um^2/s), see stage7_validation/ep_polymer.py. Acquisition: stage7_validation/ep_acquisition.py ("2d", "zstack", "biplane"; Vutara numbers: NA 1.3, n 1.406,
100 nm pixels, 575 / 665 nm labels, 0.15 s between colours, planes 0.8 um apart / 0.6 um for the biplane). Estimator: stage7_validation/ep_history.py (P(ON) from the distance
history, errors in variables, cell bootstrap) unless the experiment says "instant".
Results: $SIMLIVE_DATA/sweeps/ep_polymer_<exp>/results.csv (sweep layout: settings, seed, run, results).
Reference of every setting = the same estimator on the true distance and true ON state of 3000 cells.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import itertools
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace

import numpy as np
import pandas as pd

from simlive.io.runs import DATA_ROOT
from simlive.stage7_validation import ep_acquisition as aq
from simlive.stage7_validation import ep_curve as ep
from simlive.stage7_validation import ep_history as eh
from simlive.stage7_validation import ep_polymer as pp

N_BOOT = 100
BASE = ep.EPSetup()
CALIB0 = dict(working_nm=400.0, contact_fraction=0.10, D_um2_s=0.002, lag_s=60.0, applies_to="locus")
KINDS = ("2d", "zstack", "biplane")
_POOL = None
_REFS: dict = {}


def pool() -> pp.PolymerPool:
    global _POOL
    if _POOL is None:
        _POOL = pp.load_pool()
    return _POOL


def _ckey(c: dict) -> tuple:
    return tuple(sorted(c.items()))


def _skey(s: ep.EPSetup) -> tuple:
    return (s.working_distance_nm, s.hill_n, s.fold, s.frame_s, s.movie_s, s.exposure_s, tuple(sorted(s.tx.items())))


def make_ref(args) -> tuple:
    """Reference (true distance, true ON, 3000 cells) for one (calibration, setup, estimator)."""
    calib_args, setup, estimator = args
    cal = pp.calibrate(pool(), **calib_args)
    pop = pp.simulate_population(replace(setup, n_cells=3000), pool(), cal, 12345, with_trace=False, keep_blocks=False)
    K0 = aq.kernel_aniso(9.0, 9.0, 3, 0.0)
    a = eh.analyse_history(pop, pop.d_true, K0) if estimator == "history" else ep.analyse_measured(pop, pop.d_true, K0)
    hist, _ = np.histogram(pop.d_true, bins=np.r_[ep.D_GRID - 5.0, ep.D_GRID[-1] + 5.0])
    a["d_pdf"] = hist / hist.sum()
    a["contact_fraction_observed"] = float((pop.d_true < setup.working_distance_nm).mean())
    a["p_on_overall"] = float(pop.on_true.mean())
    a["nm_per_unit"], a["block_s"] = cal.nm_per_unit, cal.block_s
    return (_ckey(calib_args), _skey(setup), estimator), a


def metrics(a: dict, ref: dict) -> dict:
    f, r = a["fit"], ref["fit"]
    row = {"contrast": a["contrast"], "contrast_ref": r.contrast(), "contrast_ratio": a["contrast"] / r.contrast(), "d_half": f.d_half, "d_half_ref": r.d_half,
           "p_far": f.p_far, "p_far_ref": r.p_far, "curve_err": ep.curve_error(f, ref), "frames_used": a["frames_used"], "n_on": a["n_on"],
           "tau_s": getattr(f, "tau_s", np.nan), "tau_s_ref": getattr(r, "tau_s", np.nan)}
    if "contrast_ci" in a:
        lo, hi = a["contrast_ci"]
        row.update({"contrast_lo": lo, "contrast_hi": hi, "contrast_ci_width": hi - lo, "contrast_covers": float(lo <= r.contrast() <= hi), "excludes_one": float(lo > 1.0 or hi < 1.0),
                    "d_half_lo": a["d_half_ci"][0], "d_half_hi": a["d_half_ci"][1], "d_half_covers": float(a["d_half_ci"][0] <= r.d_half <= a["d_half_ci"][1])})
    return row


def task(spec: dict) -> dict:
    t0 = time.time()
    setup, source = spec["setup"], spec["source"]
    if source == "invented":
        pop = ep.simulate_population(setup, spec["seed"], coupled=spec["coupled"], with_trace=False)
    else:
        cal = pp.calibrate(pool(), **spec["calib_args"])
        pop = pp.simulate_population(setup, pool(), cal, spec["seed"], coupled=spec["coupled"], with_trace=False, keep_blocks=spec["level"] != "L0" and spec["acq"] is not None)
    lvl, est = spec["level"], spec["estimator"]
    if spec["acq"] is None:                                           # distance-only measurement (isotropic noise), no acquisition model
        if lvl == "L0":
            d, K = pop.d_true, aq.kernel_aniso(9.0, 9.0, 3, 0.0)
        else:
            m = spec["meas"]
            d = ep.measure_distance(pop.d_true, m, np.random.default_rng(spec["seed"]))
            K = ep.kernel_matrix(m.sigma_nm, m.dims, 0.0)
        extra = {}
    else:
        m = aq.measure(pop, pop.blocks, spec["acq"], seed=spec["seed"])
        d, K = m.d, aq.kernel_aniso(m.sigma_xy, m.sigma_z, m.dims, m.shift_nm)
        extra = {"sigma_xy_nm": m.sigma_xy, "sigma_z_nm": m.sigma_z, "detected": m.detected, **{k: v for k, v in m.parts.items() if k != "localization_z_nm" or m.dims == 3}}
    a = eh.analyse_history(pop, d, K, n_boot=N_BOOT, seed=spec["seed"]) if est == "history" else ep.analyse_measured(pop, d, K, n_boot=N_BOOT, seed=spec["seed"])
    ref = spec["ref"]
    return {**spec["settings"], "seed": spec["seed"], "run": spec["name"], "seconds": time.time() - t0, **metrics(a, ref), **extra}


def specs_for(exp: str, seeds: int) -> list[dict]:
    out = []

    def add(settings, setup=BASE, calib=None, acq=None, level="L2", estimator="history", coupled=True, source="polymer", meas=None, n=seeds, seed0=1000):
        calib = {**CALIB0, **(calib or {})}
        for k in range(n):
            out.append({"setup": setup, "calib_args": calib, "acq": acq, "meas": meas, "level": level, "estimator": estimator, "coupled": coupled, "source": source,
                        "settings": settings, "seed": seed0 + k, "name": f"{exp}_{len(out)}"})

    A = {k: aq.Acquisition(kind=k) for k in KINDS}
    if exp == "modality":
        add({"modality": "truth (L0)"}, level="L0")
        for k in KINDS:
            add({"modality": k}, acq=A[k], n=max(seeds, 20))
    elif exp == "null":
        for k in KINDS:
            add({"modality": k, "coupling": "off"}, acq=A[k], coupled=False, n=max(seeds, 60), seed0=5000)
    elif exp == "wrong":
        for k, w in itertools.product(KINDS, (0.0, 0.02, 0.05, 0.1, 0.2)):
            add({"modality": k, "wrong_pick_rate": w}, acq=replace(A[k], wrong_pick_rate=w))
    elif exp == "cells":
        for n_cells, k in itertools.product((100, 200, 400, 800, 1600, 3200), ("L0", "zstack", "biplane")):
            add({"modality": "truth (L0)" if k == "L0" else k, "n_cells": n_cells}, setup=replace(BASE, n_cells=n_cells), acq=None if k == "L0" else A[k],
                level="L0" if k == "L0" else "L2")
    elif exp == "contact":
        for cf, k in itertools.product((0.03, 0.1, 0.2), ("L0", "zstack", "biplane")):
            add({"modality": "truth (L0)" if k == "L0" else k, "contact_fraction": cf}, calib={"contact_fraction": cf}, acq=None if k == "L0" else A[k], level="L0" if k == "L0" else "L2")
    elif exp == "frames":
        for fs, k in itertools.product((10.0, 20.0, 40.0, 80.0), ("L0", "zstack", "biplane")):
            add({"modality": "truth (L0)" if k == "L0" else k, "frame_s": fs}, setup=replace(BASE, frame_s=fs), acq=None if k == "L0" else A[k], level="L0" if k == "L0" else "L2")
    elif exp == "colour":
        for dl, k in itertools.product((0.0, 0.05, 0.15, 0.3), KINDS):
            add({"modality": k, "colour_delay_s": dl}, acq=replace(A[k], colour_delay_s=dl, plane_time_s=2 * dl if k == "zstack" and dl > 0 else (0.15 if k == "zstack" else 0.0)))
    elif exp == "photons":
        for ps, k in itertools.product((0.03, 0.1, 0.3, 1.0, 3.0), ("zstack", "biplane")):
            add({"modality": k, "photon_scale": ps}, acq=replace(A[k], photon_scale=ps))
    elif exp == "focus":
        for fsd, k in itertools.product((0.0, 0.3, 0.6, 1.0, 1.5), ("2d", "biplane")):
            add({"modality": k, "focus_sd_um": fsd}, acq=replace(A[k], focus_sd_um=fsd))
    elif exp == "timescale":
        for (ap, lag), k in itertools.product((("locus", 10.0), ("locus", 30.0), ("locus", 60.0), ("locus", 120.0), ("pair", 60.0)), ("L0", "biplane")):
            add({"modality": "truth (L0)" if k == "L0" else k, "D_applies_to": ap, "D_lag_s": lag}, calib={"applies_to": ap, "lag_s": lag}, acq=None if k == "L0" else A[k],
                level="L0" if k == "L0" else "L2")
    elif exp == "compare_sources":
        iso = ep.Measurement(sigma_nm=57.0, dims=3, miss_rate=0.0, outlier_rate=0.0)
        for src, est, lvl in itertools.product(("invented", "polymer"), ("instant", "history"), ("L0", "L2")):
            add({"distance_source": src, "estimator": est, "level": lvl}, source=src, estimator=est, level=lvl, meas=iso, n=max(seeds, 16))
    else:
        raise ValueError(exp)
    return out


def _fill_refs(specs: list[dict], workers: int) -> None:
    keys = {}
    for sp in specs:
        k = (_ckey(sp["calib_args"]), _skey(sp["setup"]), sp["estimator"], sp["source"])
        keys.setdefault(k, (sp["calib_args"], sp["setup"], sp["estimator"], sp["source"]))
    todo = [v for v in keys.values()]
    print(f"  {len(todo)} references", flush=True)

    def ref_inputs(v):
        return (v[0], v[1], v[2])

    results = {}
    with ProcessPoolExecutor(min(workers, len(todo))) as ex:
        for v, (key, a) in zip(todo, ex.map(_ref_job, [(v[0], v[1], v[2], v[3]) for v in todo])):
            results[(key[0], key[1], key[2], v[3])] = a
    for sp in specs:
        sp["ref"] = results[(_ckey(sp["calib_args"]), _skey(sp["setup"]), sp["estimator"], sp["source"])]


def _ref_job(args):
    calib_args, setup, estimator, source = args
    if source == "invented":
        pop = ep.simulate_population(replace(setup, n_cells=3000), 12345, with_trace=False)
        K0 = aq.kernel_aniso(9.0, 9.0, 3, 0.0)
        a = eh.analyse_history(pop, pop.d_true, K0) if estimator == "history" else ep.analyse_measured(pop, pop.d_true, K0)
        hist, _ = np.histogram(pop.d_true, bins=np.r_[ep.D_GRID - 5.0, ep.D_GRID[-1] + 5.0])
        a["d_pdf"] = hist / hist.sum(); a["contact_fraction_observed"] = float((pop.d_true < setup.working_distance_nm).mean()); a["p_on_overall"] = float(pop.on_true.mean())
        return (_ckey(calib_args), _skey(setup), estimator), a
    return make_ref((calib_args, setup, estimator))


def write(exp: str, rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    res = ("contrast", "d_half", "p_far", "curve_err", "frames_used", "n_on", "tau_s", "excludes", "sigma_", "detected", "localization", "exposure_blur", "colour_delay", "piezo", "mean_",
           "sensitivity", "false_pos", "rule_", "noise_")
    lead = [c for c in df.columns if c not in ("seed", "run", "seconds") and not c.startswith(res)]
    cols = lead + ["seed", "run", "seconds"] + [c for c in df.columns if c not in lead + ["seed", "run", "seconds"]]
    out = DATA_ROOT / "sweeps" / f"ep_polymer_{exp}"
    out.mkdir(parents=True, exist_ok=True)
    df[cols].to_csv(out / "results.csv", index=False)
    print(f"wrote {out / 'results.csv'} ({len(df)} rows)", flush=True)
    return df


def run_ms2(seeds: int, workers: int) -> list[dict]:
    """L3: the ON state estimated from the MS2 trace (deconvolution, burst rule (3, 3) fixed in advance), measured distances, with / without the attenuation correction."""
    rows = []
    for kind, sb, fs in itertools.product(("biplane", "zstack"), (0.05, 0.4, 1.0, 2.0), (10.0, 20.0)):
        setup = replace(BASE, n_cells=300, frame_s=fs)
        calib_args = dict(CALIB0)
        ref = _ref_job((calib_args, setup, "history", "polymer"))[1]
        cal = pp.calibrate(pool(), **calib_args)
        for k in range(seeds):
            t0 = time.time()
            seed = 2000 + k
            pop = pp.simulate_population(setup, pool(), cal, seed, noise_sb=sb, noise_beta=0.03)
            est = ep.estimate_on(pop, workers=workers, seed=seed)
            p_on = est.p_on(3.0, 3.0)
            sc = est.scores(p_on)
            m = aq.measure(pop, pop.blocks, aq.Acquisition(kind=kind), seed=seed)
            K = aq.kernel_aniso(m.sigma_xy, m.sigma_z, m.dims, m.shift_nm)
            for corrected in (False, True):
                a = eh.analyse_history(pop, m.d, K, y=p_on, correct_calls=(sc["sensitivity"], sc["false_positive"]) if corrected else None, n_boot=N_BOOT, seed=seed)
                rows.append({"modality": kind, "ms2_sb": sb, "frame_s": fs, "corrected": corrected, "seed": seed, "run": f"ms2_{len(rows)}", "seconds": time.time() - t0,
                             "sensitivity": sc["sensitivity"], "false_positive": sc["false_positive"], "noise_sb_est": est.noise_sb, **metrics(a, ref)})
            print(f"  ms2 {kind} sb {sb} frame {fs}: seed {seed} sens {sc['sensitivity']:.2f} fp {sc['false_positive']:.3f} ({time.time() - t0:.0f}s)", flush=True)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exps", nargs="+", default=["modality", "null", "wrong", "cells", "contact", "frames", "colour", "photons", "focus", "timescale", "compare_sources"])
    ap.add_argument("--seeds", type=int, default=12)
    ap.add_argument("--ms2-seeds", type=int, default=4)
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    for exp in a.exps:
        t0 = time.time()
        print(f"== {exp}", flush=True)
        if exp == "ms2":
            rows = run_ms2(a.ms2_seeds, a.workers)
        else:
            specs = specs_for(exp, a.seeds)
            _fill_refs(specs, a.workers)
            with ProcessPoolExecutor(a.workers) as ex:
                rows = list(ex.map(task, specs, chunksize=2))
        df = write(exp, rows)
        print(f"   {time.time() - t0:.0f} s", flush=True)
        key = [c for c in df.columns if c in ("modality", "coupling", "wrong_pick_rate", "n_cells", "contact_fraction", "frame_s", "colour_delay_s", "photon_scale", "focus_sd_um", "D_applies_to",
                                              "D_lag_s", "distance_source", "estimator", "level", "ms2_sb", "corrected")]
        agg = df.groupby(key, dropna=False).agg(contrast=("contrast", "mean"), ref=("contrast_ref", "mean"), sd=("contrast", "std"), curve_err=("curve_err", "mean"),
                                                cover=("contrast_covers", "mean"), excl1=("excludes_one", "mean"), ci_w=("contrast_ci_width", "mean"))
        print(agg.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
