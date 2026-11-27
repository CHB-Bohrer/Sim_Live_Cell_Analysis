r"""Short validation checks whose numbers and figures appear in the paper.

    scripts\run.cmd python docs\paper\analysis\checks.py            # run all (CPU only; seconds once numba has compiled)
    scripts\run.cmd python docs\paper\analysis\checks.py probes     # run one

Each check feeds a model an input with a known answer and compares what comes out.
Limits are fixed here, in advance (they are the tolerances of the matching unit tests).
Results go to docs/paper/data/checks.json with a run manifest (commit, date, seed,
package versions); a check that is not re-run keeps its stored result.
"""
import datetime
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
OUT = HERE.parent / "data" / "checks.json"


def check_probes():
    """Probe binding: exact stochastic simulation against the exact stationary theory."""
    from simlive.stage1_chromatin import probes as P
    rows, ok = [], True
    for k_off in (0.03, 0.06, 0.1, 0.15, 0.25, 0.4):
        p = P.ProbeParams(n_probes=40, k_bind=0.1, k_scan=0.4, k_unbind=0.2, k_off=k_off)
        th = P.summary(p)
        d = P.sample_trace(P.simulate_trace(p, 60000.0, seed=3), np.arange(0, 60000.0, 1.0), 1e-3)
        a = d.n_attached.to_numpy()
        row = {"k_off": k_off, "mean_sim": float(a.mean()), "mean_theory": float(th["mean_attached"]),
               "cv_sim": float(a.std() / a.mean()), "cv_theory": float(th["cv_attached"])}
        ok &= abs(row["mean_sim"] / row["mean_theory"] - 1) < 0.03 and abs(row["cv_sim"] / row["cv_theory"] - 1) < 0.10
        rows.append(row)
    worst = max(abs(r["mean_sim"] / r["mean_theory"] - 1) for r in rows)
    return {"title": "Probe binding model against exact theory", "passed": bool(ok), "rows": rows,
            "limit": "mean within 3%, CV within 10%, at every release rate",
            "summary": f"largest error in the mean number of attached probes: {100 * worst:.1f}%",
            "params": {"n_probes": 40, "k_bind": 0.1, "k_scan": 0.4, "k_unbind": 0.2, "duration_s": 60000, "seed": 3}}


def check_transcription():
    """Promoter + Pol II + MS2 with a distance that alternates between contact and apart."""
    from simlive.stage1_chromatin import transcription as T
    fast = dict(k_on=0.004, k_off=0.05, k_init=0.2, coupling="contact", fold=15.0, f_const=0.0, elongation_kb_min=6.0,
                cassette_kb=0.6, gene_kb=2.4, n_loops=24, dwell_s=10.0)
    p = T.TranscriptionParams(**fast)
    d = np.where((np.arange(400) // 20) % 2 == 0, 80.0, 400.0)
    ft = np.arange(100, 4000.0, 25.0)
    n_cells = 300
    runs = [T.simulate_cell(p, d, 10.0, ft, 0.1, seed=s).frames for s in range(n_cells)]
    on = np.mean([r.promoter_on for r in runs], axis=0)
    ms2 = np.mean([r.ms2_loops for r in runs], axis=0)
    p_start, _ = T.expected_on(p, np.tile(d, 10), 10.0)
    exp_on = p_start[np.searchsorted(np.arange(len(p_start)) * 10.0, ft, side="right") - 1]
    exp_ms2 = T.expected_ms2(p, np.tile(d, 10), 10.0, ft)
    late = ft > 200
    e_on = float(abs(on.mean() - exp_on.mean()))
    e_ms2 = float(abs(ms2[late].mean() / exp_ms2[late].mean() - 1))
    ok = np.abs(on - exp_on).max() < 0.12 and e_on < 0.02 and e_ms2 < 0.10
    return {"title": "Transcription model against the exact expectation", "passed": bool(ok),
            "t_s": ft.tolist(), "on_sim": on.tolist(), "on_exact": exp_on.tolist(),
            "ms2_sim": ms2.tolist(), "ms2_exact": exp_ms2.tolist(), "contact": (d[(ft // 10).astype(int) % 400] < 150).tolist(),
            "limit": "mean P(ON) within 0.02; mean MS2 signal within 10%; no time point off by more than 0.12",
            "summary": f"mean P(ON) off by {e_on:.3f}; mean MS2 signal off by {100 * e_ms2:.1f}% ({n_cells} cells)",
            "params": {**fast, "n_cells": n_cells, "seeds": "0-299", "block_s": 10.0}}


def check_ep_recovery():
    """A known distance -> P(ON) curve, seen through 2D projection and localization error: is it recovered?"""
    from simlive.stage7_validation import ep_curve as ep
    rng = np.random.default_rng(2)
    n = 300000
    d = np.where(rng.random(n) < 0.12, rng.uniform(100, 380, n), 400 * 2 * np.exp(0.35 * rng.normal(size=n)))
    truth = (0.05, 8.0, 350.0)
    y = (rng.random(n) < ep.curve(d, *truth, 4.0)).astype(float)
    cells = rng.integers(0, 300, n)
    dm = ep.measure_distance(d, ep.Measurement(sigma_nm=70.0, dims=2), rng)
    agg = ep.Aggregates.build(cells, dm, y)
    K = ep.kernel_matrix(70.0, 2, 0.0)
    prior = ep.population_prior(agg.n.sum(0), K)
    fit = ep.fit_curve(agg, K, prior)
    naive = ep.fit_curve(agg, K, prior, naive=True)
    grid = ep.D_GRID
    sel = (grid > 100) & (grid < 1000)
    true_curve = ep.curve(grid, *truth, 4.0)
    err = lambda q: float(np.abs(q.g(grid) - true_curve)[sel].mean())     # noqa: E731
    contrast = lambda g: float(np.interp(150.0, grid, g) / np.interp(1000.0, grid, g))   # noqa: E731
    ok = abs(fit.fold / 8.0 - 1) < 0.2 and abs(fit.d_half / 350.0 - 1) < 0.15 and err(fit) < err(naive)
    return {"title": "Recovery of a known distance-to-activity curve through measurement error", "passed": bool(ok),
            "d_nm": grid[sel].tolist(), "true": true_curve[sel].tolist(), "fit": fit.g(grid)[sel].tolist(),
            "naive": naive.g(grid)[sel].tolist(),
            "fold_true": 8.0, "fold_fit": float(fit.fold), "d_half_true": 350.0, "d_half_fit": float(fit.d_half),
            "contrast_true": contrast(true_curve), "contrast_fit": contrast(fit.g(grid)), "contrast_naive": contrast(naive.g(grid)),
            "limit": "fold within 20%, half-distance within 15%, and closer to the truth than the uncorrected fit",
            "summary": (f"fold {fit.fold:.2f} (true 8), half-distance {fit.d_half:.0f} nm (true 350); "
                        f"contrast {contrast(fit.g(grid)):.2f} corrected, {contrast(naive.g(grid)):.2f} uncorrected, "
                        f"{contrast(true_curve):.2f} true"),
            "params": {"n_measurements": n, "n_cells": 300, "sigma_nm": 70.0, "dims": 2, "seed": 2}}


CHECKS = {"probes": check_probes, "transcription": check_transcription, "ep_recovery": check_ep_recovery}


def manifest():
    def git(*a):
        try:
            return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
        except OSError:
            return ""
    import numba
    import scipy
    return {"commit": git("rev-parse", "--short", "HEAD"), "dirty": bool(git("status", "--porcelain", "src")),
            "date": datetime.date.today().isoformat(), "python": platform.python_version(),
            "numpy": np.__version__, "scipy": scipy.__version__, "numba": numba.__version__, "platform": platform.system()}


def main():
    names = sys.argv[1:] or list(CHECKS)
    data = json.loads(OUT.read_text(encoding="utf-8")) if OUT.exists() else {}
    failed = False
    for name in names:
        res = CHECKS[name]()
        res["manifest"] = manifest()
        data[name] = res
        failed |= not res["passed"]
        print(f"{'PASS' if res['passed'] else 'FAIL'}  {name}: {res['summary']}")
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    lines = ["| Check | Compared against | Limit (set in advance) | Result | Run |", "| --- | --- | --- | --- | --- |"]
    kind = {"probes": "exact theory", "transcription": "exact expectation", "ep_recovery": "the curve put in"}
    for name in CHECKS:
        if name in data:
            r = data[name]; m = r["manifest"]
            lines.append(f"| {r['title']} | {kind[name]} | {r['limit']} | **{'pass' if r['passed'] else 'FAIL'}**: {r['summary']} | "
                         f"{m['date']}, commit `{m['commit']}` |")
    (OUT.parent / "checks-table.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
