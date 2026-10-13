"""Convergence check of the MS2 -> Pol II samplers (docs/MCMC_CONVERGENCE.md).

  scripts\\run.cmd python scripts\\check_mcmc.py --run loci_demo_seed8_f5c88d            # every cell of a run, 3 samplers, 4 chains each
  scripts\\run.cmd python scripts\\check_mcmc.py --synthetic                              # constructed hard cases (known truth)

Writes <run>/stage6_analysis/<src>/mcmc_check.csv (one row per cell and sampler) and mcmc_burn_in.csv, or
$SIMLIVE_DATA/mcmc_check/synthetic.csv. Limits (R-hat < 1.05, ESS >= 400) are fixed in stage6_analysis/mcmc_check.py, not tuned here.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from simlive.io.runs import DATA_ROOT
from simlive.stage6_analysis.deconvolve import SAMPLERS, make_kernel, prepare_traces
from simlive.stage6_analysis.mcmc_check import burn_in_curve, check_convergence


def _trace(kern, starts_s, T, dt, sb, beta, rng):
    t = np.arange(T) * dt
    G = np.zeros(T)
    for s0 in starts_s:
        u = t - s0
        ok = (u > 0) & (u / kern.dt_fine_s < len(kern.values))
        G[ok] += kern.values[(u[ok] / kern.dt_fine_s).astype(int)]
    return G + rng.normal(size=T) * np.sqrt(sb ** 2 + beta * G)


def synthetic_cases(seed: int = 0) -> list[dict]:
    """Constructed traces with a known answer: sparse, dense bursts, a close pair (one start or two?), a start very near the movie start."""
    rng = np.random.default_rng(seed)
    k = make_kernel(24, 1.3, 5.0, 2.5, 116.0, 1.0)
    T, dt, sb, beta = 180, 10.0, 0.15, 0.02
    out = []
    for name, starts in (("sparse (5 starts)", [105.0, 140.0, 700.0, 1290.0, 1500.0]),
                         ("dense (35 starts, two bursts)", np.sort(np.r_[rng.uniform(100, 500, 20), rng.uniform(800, 1100, 15)])),
                         ("close pair (6 s apart)", [400.0, 406.0, 1200.0]),
                         ("one start at 40 s, two at 1700 s", [40.0, 1700.0, 1708.0]),
                         ("no start at all", [])):
        a = _trace(k, starts, T, dt, sb, beta, rng)
        rate = max(len(starts), 1) / (T * dt)
        out.append({"name": name, "a": a, "kern": k, "dt": dt, "prior": np.full(T, rate), "sb": sb, "beta": beta, "n_true": len(starts)})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run")
    ap.add_argument("--source", default="truth")
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--samplers", nargs="+", default=list(SAMPLERS), choices=SAMPLERS)
    ap.add_argument("--n-burn", type=int, default=500)
    ap.add_argument("--n-keep", type=int, default=1000)
    ap.add_argument("--cells", type=int, default=0, help="only the first N cells (0 = all)")
    a = ap.parse_args()
    rows, burn = [], []

    def go(label, cid, tr, kern, dt, prior, sb, beta, extra=None):
        for smp in a.samplers:
            cc = check_convergence(tr, kern, dt, prior, sb, beta, sampler=smp, n_burn=a.n_burn, n_keep=a.n_keep, seed=int(cid))
            r = {"case": label, "cell": cid, "n_starts_mean": float(np.mean([c["n_starts"].mean() for c in cc.chains])), **(extra or {}), **cc.row()}
            rows.append(r)
            print(f"{label:34s} {smp:12s} R-hat lp {cc.rhat_logpost:5.3f} n {cc.rhat_nstarts:5.3f} frames {cc.rhat_frames_max:5.3f} | "
                  f"ESS lp {cc.ess_logpost:6.0f} n {cc.ess_nstarts:6.0f} frames min {cc.ess_frames_min:6.0f} | diff {cc.max_chain_diff:5.2f}/min | "
                  f"{cc.seconds:5.1f}s | {'OK' if cc.converged else 'NOT CONVERGED: ' + '; '.join(cc.reasons)}", flush=True)
            if smp == "gibbs_shift":                                    # burn-in curve: chains started at n_burn = 0 on the same trace
                short = check_convergence(tr, kern, dt, prior, sb, beta, sampler=smp, n_burn=600, n_keep=20, seed=int(cid))
                for b in burn_in_curve(short.chains, sweeps=(0, 5, 10, 25, 50, 100, 250, 500), window=100):
                    burn.append({"case": label, "cell": cid, **b})

    if a.synthetic:
        for c in synthetic_cases():
            go(c["name"], 0, c["a"], c["kern"], c["dt"], c["prior"], c["sb"], c["beta"], {"n_true": c["n_true"]})
        out = DATA_ROOT / "mcmc_check"
    else:
        run = Path(a.run) if Path(a.run).exists() else DATA_ROOT / "runs" / a.run
        cfg = yaml.safe_load((run / "stage2_cells" / "params.yaml").read_text())
        d6 = run / "stage6_analysis" / a.source
        f = next(d6 / n for n in ("ms2_traces_cal.csv", "ms2_traces.csv") if (d6 / n).exists())
        traces, kern, dt, prior, sb, beta, _, _ = prepare_traces(pd.read_csv(f), cfg)
        for i, (cid, tr) in enumerate(traces.items()):
            if a.cells and i >= a.cells:
                break
            go(f"{run.name} cell {cid}", cid, tr, kern, dt, prior, sb, beta)
        out = d6
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / ("synthetic.csv" if a.synthetic else "mcmc_check.csv"), index=False)
    pd.DataFrame(burn).to_csv(out / ("synthetic_burn_in.csv" if a.synthetic else "mcmc_burn_in.csv"), index=False)
    df = pd.DataFrame(rows)
    print("\nconverged per sampler:\n", df.groupby("sampler").agg(cases=("converged", "size"), converged=("converged", "sum"),
          rhat_max=("rhat_frames_max", "max"), ess_min=("ess_nstarts", "min"), seconds=("seconds", "sum")).round(2))
    print("wrote", out)


if __name__ == "__main__":
    main()
