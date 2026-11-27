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


def check_calibration():
    """Blind calibration: simulate dark / flat / bead acquisitions of one instrument, estimate it without looking."""
    from simlive.io.runs import REPO_ROOT, load_config
    from simlive.stage3_microscopy.calibration.run import run_calibration, settings_for
    from simlive.stage3_microscopy.render import entropy_for
    short = {"dark": {"n_frames": 40, "exposures_s": [0.01, 1.0, 5.0], "extra_frames": 8},
             "flat": {"levels_photons": [0, 300, 1000, 3000, 8000, 20000], "n_frames": 10},
             "flicker": {"level_photons": 3000, "n_frames": 80},
             "beads": {"n_beads": 30, "photons": 30000, "n_frames": 6, "n_frames_3d": 4, "z_um": 1.0}}
    cfg = load_config(REPO_ROOT / "configs" / "loci_demo.yaml", ["geometry.fov_um=[30, 30]"])
    _, table = run_calibration(cfg, entropy_for(3), {**settings_for(cfg), **short}, progress=lambda *_: None)
    # limits = the tolerances of tests/test_calibration.py (percent unless noted)
    limits = {"gain": 3, "read noise": 5, "row noise": 15, "offset pattern (DSNU)": 10, "dark current": 15,
              "pixel gain spread (PRNU)": 10, "PSF width (lateral sigma)": 3}
    rows, ok = [], True
    for _, r in table.iterrows():
        if r.truth != r.truth or r.estimate != r.estimate or float(r.truth) == 0:
            continue
        row = {"quantity": r.quantity, "channel": r.channel, "estimate": float(r.estimate), "truth": float(r.truth),
               "ratio": float(r.estimate) / float(r.truth)}
        if r.quantity in limits:
            ok &= abs(row["ratio"] - 1) * 100 < limits[r.quantity]
        elif r.quantity.startswith("chromatic shift"):
            ok &= abs(row["estimate"] - row["truth"]) < 0.03
        elif r.quantity == "hot pixels found":
            ok &= row["estimate"] == row["truth"]
        rows.append(row)
    within = sum(abs(x["ratio"] - 1) < 0.05 for x in rows)
    return {"title": "Blind calibration against the simulated instrument", "passed": bool(ok), "rows": rows,
            "limit": "gain, PSF width within 3%; read noise 5%; DSNU, PRNU 10%; dark current, row noise 15%; colour shift 0.03 px; every hot pixel",
            "summary": f"{within} of {len(rows)} estimates within 5% (short acquisitions; flicker and cosmic rays rough by design)",
            "params": {"config": "loci_demo", "fov_um": 30, "instrument_seed": 3, **{k: str(v) for k, v in short.items()}}}


def check_spot_precision():
    """Spot finder (DoG + Gaussian fit) on drawn spots: error against the physical limit (Mortensen et al. 2010)."""
    from scipy.special import erf
    from simlive.stage6_analysis import dog_viterbi as DV
    from simlive.stage6_analysis import spots as SP
    S, sigma, bg, read = 41, 0.85, 20.0, 1.5

    def spot(cy, cx, n):
        s2, i = np.sqrt(2) * sigma, np.arange(S)
        gy = 0.5 * (erf((i + 0.5 - cy) / s2) - erf((i - 0.5 - cy) / s2))
        gx = 0.5 * (erf((i + 0.5 - cx) / s2) - erf((i - 0.5 - cx) / s2))
        return n * np.outer(gy, gx)

    yy, xx = np.indices((S, S))
    region = np.hypot(yy - S // 2, xx - S // 2) <= 16
    rng = np.random.default_rng(1)
    sa2, b = sigma ** 2 + 1 / 12, bg + read ** 2
    rows = []
    for n in (150, 250, 400, 700, 1200, 2000, 3500, 6000, 10000):
        err, found = [], 0
        for _ in range(200):
            cy, cx = rng.uniform(17, 23, 2)
            img = rng.poisson(spot(cy, cx, n) + bg).astype(float) + rng.normal(0, read, (S, S)) + 100.0
            c = DV.dog_candidates(img, region)
            if c is None or not len(c["peaks"]):
                continue
            r = SP.refine_peak(c, 0, **SP.DEFAULTS)
            if np.hypot(r["pos"][0] - cy, r["pos"][1] - cx) < 3:
                found += 1
                err += [r["pos"][0] - cy, r["pos"][1] - cx]
        rms = float(np.sqrt(np.mean(np.square(err))))
        limit = float(np.sqrt(sa2 / n * (16 / 9 + 8 * np.pi * sa2 * b / n)))
        rows.append({"photons": n, "rms_px": rms, "limit_px": limit, "ratio": rms / limit, "found": found / 200})
    ok = all(r["ratio"] < 1.25 and r["found"] > 0.98 for r in rows)
    return {"title": "Spot position error against the physical limit", "passed": bool(ok), "rows": rows,
            "limit": "error within 1.25x the Mortensen limit and > 98% found, at every brightness",
            "summary": f"error {min(r['ratio'] for r in rows):.2f} to {max(r['ratio'] for r in rows):.2f}x the limit, 150 to 10,000 photons",
            "params": {"psf_sigma_px": sigma, "background_photons": bg, "read_noise": read, "spots_per_level": 200, "seed": 1}}


def check_deconvolution():
    """MS2 -> Pol II samplers on a problem small enough to sum all 2^14 start patterns exactly (tests/test_mcmc_check.py)."""
    import itertools
    from simlive.stage6_analysis.deconvolve import Kernel, deconvolve_trace
    kern = Kernel(1.0, np.array([0.0, 0.5, 1.0, 1.0, 0.0]))
    T, s, K = 5, 2, 5
    N = K + (T - 1) * s + 1
    idx = K + s * np.arange(T)

    def model(I):
        G = np.zeros(T)
        for n in range(T):
            u = idx[n] - np.flatnonzero(I)
            G[n] = kern.values[u[(u > 0) & (u < K)]].sum()
        return G

    rng = np.random.default_rng(5)
    truth = np.zeros(N, int); truth[[6, 11]] = 1
    sb, beta, p = 0.25, 0.05, 0.08
    G0 = model(truth)
    a = G0 + rng.normal(size=T) * np.sqrt(sb ** 2 + beta * G0)
    logw, states = [], []
    for bits in itertools.product((0, 1), repeat=N):
        I = np.array(bits); G = model(I); v = sb ** 2 + beta * np.clip(G, 0, None)
        logw.append((-0.5 * (a - G) ** 2 / v - 0.5 * np.log(v)).sum() + I.sum() * np.log(p / (1 - p))); states.append(I)
    w = np.exp(np.array(logw) - max(logw)); w /= w.sum()
    exact = (w[:, None] * np.array(states)).sum(0)
    samplers = {"gibbs_shift": {}, "gibbs_pt": {"n_replicas": 4, "beta_min": 0.1}, "chen_block": {"block_s": 4.0}}
    est = {name: deconvolve_trace(a, kern, 2.0, np.full(T, p), sb, beta, n_burn=300, n_keep=20000, thin=3, max_shift_s=3.0,
                                  seed=3, sampler=name, **kw)["p_fine"].tolist() for name, kw in samplers.items()}
    worst = {name: float(np.abs(np.array(e) - exact).max()) for name, e in est.items()}
    return {"title": "MS2 deconvolution samplers against the exact posterior", "passed": bool(max(worst.values()) < 0.03),
            "exact": exact.tolist(), "samplers": est, "worst": worst,
            "limit": "every start probability within 0.03 of the exact value, for every sampler",
            "summary": "largest error " + ", ".join(f"{k} {v:.3f}" for k, v in worst.items()),
            "params": {"fine_steps": int(N), "patterns": 2 ** int(N), "n_keep": 20000, "seed": 3}}


CHECKS = {"probes": check_probes, "transcription": check_transcription, "ep_recovery": check_ep_recovery,
          "calibration": check_calibration, "spot_precision": check_spot_precision, "deconvolution": check_deconvolution}


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
    kind = {"probes": "exact theory", "transcription": "exact expectation", "ep_recovery": "the curve put in",
            "calibration": "the simulated instrument", "spot_precision": "the physical limit", "deconvolution": "the exact posterior"}
    for name in CHECKS:
        if name in data:
            r = data[name]; m = r["manifest"]
            lines.append(f"| {r['title']} | {kind[name]} | {r['limit']} | **{'pass' if r['passed'] else 'FAIL'}**: {r['summary']} | "
                         f"{m['date']}, commit `{m['commit']}` |")
    (OUT.parent / "checks-table.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
