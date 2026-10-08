"""Stochastic model of fluorescent probes binding at a locus (the promoter and the enhancer each get their own).

Every probe is in one of three states and moves through the cycle

    FREE --(k_bind)--> BOUND --(k_scan)--> SCANNING --(k_off)--> FREE
                          \\--(k_unbind)--> FREE          (a probe that lands but leaves again without scanning)

  FREE      floating in the nucleoplasm near the locus; lands on the locus at rate k_bind per free probe.
  BOUND     just landed on the locus. Either starts scanning (k_scan) or lets go again (k_unbind).
  SCANNING  sliding along the locus looking for its site. Falls off (k_off).
Probes are fluorescent while BOUND or SCANNING ("attached"), so the locus brightness is proportional to
n_attached = n_bound + n_scanning, and this number fluctuates randomly in time. That shot-like binding noise is what
makes the locus signal imperfect in the simulated movies.

All rates are per second. The model is exact (Gillespie algorithm) and every probe is independent of the others, unless
`capacity` (the number of binding sites at the locus) is set: then landing is blocked in proportion to how full the
locus is, rate = k_bind * n_free * (1 - n_attached / capacity). Units: counts of probes; seconds.

Not modelled (kept simple on purpose): photobleaching of individual probes, probes' spatial position, the signal being
spread along the locus while scanning, and a finite pool being depleted by the other locus.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from numba import njit

LOCUS_NAMES = ("promoter", "enhancer")
# PLACEHOLDER values (the user will supply measured binding kinetics): a binder that stays ~10 s on average.
DEFAULTS = {
    "promoter": {"n_probes": 30, "k_bind": 0.05, "k_scan": 0.5, "k_unbind": 0.1, "k_off": 0.1, "capacity": None},
    "enhancer": {"n_probes": 30, "k_bind": 0.03, "k_scan": 0.3, "k_unbind": 0.1, "k_off": 0.05, "capacity": None},
}
PARAM_HELP = {
    "n_probes": "Probes in the pool that can reach this locus (free + attached). More probes = brighter and less noisy.",
    "k_bind": "Landing rate per free probe (1/s). Roughly: association rate x probe concentration.",
    "k_scan": "Rate at which a just-landed probe starts scanning (1/s).",
    "k_unbind": "Rate at which a just-landed probe lets go again without scanning (1/s). 0 = every landing leads to scanning.",
    "k_off": "Rate at which a scanning probe falls off (1/s). 1/k_off = how long the probe scans.",
    "capacity": "Binding sites at the locus. Empty = unlimited (probes never compete). Set it to make the locus saturate.",
}


@dataclass(frozen=True)
class ProbeParams:
    n_probes: int = 30
    k_bind: float = 0.05
    k_scan: float = 0.5
    k_unbind: float = 0.1
    k_off: float = 0.1
    capacity: int | None = None

    def __post_init__(self):
        if self.n_probes < 1 or min(self.k_bind, self.k_scan, self.k_off) <= 0 or self.k_unbind < 0:
            raise ValueError(f"need n_probes >= 1, k_bind/k_scan/k_off > 0, k_unbind >= 0, got {self}")
        if self.capacity is not None and self.capacity < 1:
            raise ValueError("capacity must be at least 1 (or empty for unlimited)")

    @classmethod
    def from_dict(cls, d: dict | None, locus: str = "promoter") -> "ProbeParams":
        return cls(**{**DEFAULTS[locus], **(d or {})})


def params_from_config(cfg: dict) -> dict[str, ProbeParams]:
    """{locus name: ProbeParams} from cfg['loci']['probes'] (missing entries take the DEFAULTS)."""
    pc = cfg["loci"].get("probes") or {}
    return {n: ProbeParams.from_dict(pc.get(n), n) for n in LOCUS_NAMES}


# ----------------------------------------------------------------------------- theory (exact, for the dashboard)
def generator(p: ProbeParams) -> np.ndarray:
    """Rate matrix Q of ONE probe over (FREE, BOUND, SCANNING): dP/dt = P Q."""
    return np.array([[-p.k_bind, p.k_bind, 0.0],
                     [p.k_unbind, -(p.k_scan + p.k_unbind), p.k_scan],
                     [p.k_off, 0.0, -p.k_off]])


def summary(p: ProbeParams) -> dict:
    """Exact stationary numbers for independent probes (capacity ignored; use simulate_trace for saturating loci)."""
    kb, ks, ku, ko = p.k_bind, p.k_scan, p.k_unbind, p.k_off
    pi = np.array([1.0, kb / (ks + ku), kb / (ks + ku) * ks / ko])
    pi /= pi.sum()
    p_att = pi[1] + pi[2]
    ev = np.sort(np.linalg.eigvals(generator(p)).real)          # one is 0 (stationary); the others are decay rates
    slow = -ev[ev < -1e-12].max()                               # slowest decay rate
    residence = 1 / (ks + ku) + ks / (ks + ku) / ko             # mean time a probe stays attached per visit
    return {"p_attached": p_att, "mean_free": p.n_probes * pi[0], "mean_bound": p.n_probes * pi[1],
            "mean_scanning": p.n_probes * pi[2], "mean_attached": p.n_probes * p_att,
            "sd_attached": float(np.sqrt(p.n_probes * p_att * (1 - p_att))),
            "cv_attached": float(np.sqrt((1 - p_att) / (p.n_probes * p_att))),
            "fraction_scanning": pi[2] / p_att, "residence_s": residence, "correlation_time_s": 1 / slow,
            "visits_per_s": p.n_probes * pi[0] * kb, "p_none_attached": (1 - p_att) ** p.n_probes}


# ----------------------------------------------------------------------------- exact simulation
@njit(cache=True)
def _gillespie(n, kb, ks, ku, ko, cap, t_end, seed):
    """All events of one locus from t=0 (everything free) to t_end. Returns times and (free, bound, scanning) after each."""
    np.random.seed(seed)
    nf, nb, ns = n, 0, 0
    times, F, B, S = [0.0], [n], [0], [0]
    t = 0.0
    while True:
        bind = kb * nf * (1.0 - (nb + ns) / cap) if cap > 0 else kb * nf
        if bind < 0.0:
            bind = 0.0
        r1, r2, r3, r4 = bind, ks * nb, ku * nb, ko * ns
        tot = r1 + r2 + r3 + r4
        t += -np.log(1.0 - np.random.random()) / tot
        if t >= t_end:
            break
        u = np.random.random() * tot
        if u < r1:
            nf -= 1; nb += 1
        elif u < r1 + r2:
            nb -= 1; ns += 1
        elif u < r1 + r2 + r3:
            nb -= 1; nf += 1
        else:
            ns -= 1; nf += 1
        times.append(t); F.append(nf); B.append(nb); S.append(ns)
    return np.array(times), np.array(F), np.array(B), np.array(S)


def burn_in_s(p: ProbeParams) -> float:
    """Time to forget the all-free start: 10 x the slowest relaxation time (exact for independent probes)."""
    return 10.0 * summary(p)["correlation_time_s"]


def simulate_trace(p: ProbeParams, duration_s: float, seed: int, burn_in: float | None = None):
    """Event trajectory over [0, duration_s] after discarding the burn-in. Returns (times, free, bound, scanning)."""
    b = burn_in_s(p) if burn_in is None else float(burn_in)
    t, F, B, S = _gillespie(p.n_probes, p.k_bind, p.k_scan, p.k_unbind, p.k_off, float(p.capacity or 0), b + duration_s,
                            int(seed) % (2**31))
    keep = np.searchsorted(t, b, side="right") - 1          # state in force at the start of the kept window
    return np.concatenate([[0.0], t[keep + 1:] - b]), F[keep:], B[keep:], S[keep:]


def sample_trace(trace, frame_times_s, exposure_s: float) -> pd.DataFrame:
    """Counts at the frame times, and the attached count averaged over each exposure window (what the camera sees)."""
    t, F, B, S = trace
    ft = np.asarray(frame_times_s, float)
    att = (B + S).astype(float)
    cum = np.concatenate([[0.0], np.cumsum(att[:-1] * np.diff(t))])

    def integral(u):
        i = np.clip(np.searchsorted(t, u, side="right") - 1, 0, len(t) - 1)
        return cum[i] + att[i] * (u - t[i])

    i = np.clip(np.searchsorted(t, ft, side="right") - 1, 0, len(t) - 1)
    mean = (integral(ft + exposure_s) - integral(ft)) / exposure_s
    return pd.DataFrame({"n_free": F[i], "n_bound": B[i], "n_scanning": S[i], "n_attached": B[i] + S[i],
                         "mean_attached": mean})


def simulate_probe_occupancy(cells: pd.DataFrame, cfg: dict, rng: np.random.Generator) -> pd.DataFrame:
    """Probe counts at every locus of every cell at every frame (independent stochastic models per cell and locus).

    Columns: t, cell_id, locus_id, locus_name, n_free, n_bound, n_scanning, n_attached, mean_attached (average over the
    exposure; this is what scales the rendered locus brightness). Locus 0 = promoter, locus 1 = enhancer.
    """
    dt, expo = float(cfg["acquisition"]["frame_interval_s"]), float(cfg["acquisition"]["exposure_s"])
    pars = params_from_config(cfg)
    rows = []
    for cid, g in cells.sort_values("t").groupby("cell_id"):
        ts = g.t.to_numpy()
        for k, name in enumerate(LOCUS_NAMES[: int(cfg["loci"]["n_loci"])]):
            tr = simulate_trace(pars[name], (ts.max() - ts.min() + 1) * dt + expo, int(rng.integers(2**31)))
            df = sample_trace(tr, (ts - ts.min()) * dt, expo)
            df.insert(0, "t", ts); df.insert(1, "cell_id", cid); df.insert(2, "locus_id", k); df.insert(3, "locus_name", name)
            rows.append(df)
    return pd.concat(rows, ignore_index=True)


# ----------------------------------------------------------------------------- parameter sweeps
def trace_stats(p: ProbeParams, duration_s: float, seed: int, dt_s: float = 1.0) -> dict:
    """Simulated numbers from one long trace sampled every dt_s."""
    tr = simulate_trace(p, duration_s, seed)
    d = sample_trace(tr, np.arange(0, duration_s, dt_s), 1e-3)
    a = d.n_attached.to_numpy().astype(float)
    return {"sim_mean_attached": a.mean(), "sim_cv": a.std() / a.mean() if a.mean() > 0 else np.nan,
            "sim_p_none": float((a == 0).mean()), "sim_events_per_s": (len(tr[0]) - 1) / duration_s}


def sweep(base: dict, grid: dict[str, list], locus: str = "promoter", reps: int = 3, duration_s: float = 3000.0,
          seed: int = 1) -> pd.DataFrame:
    """Every combination of `grid` ({parameter: values}) on top of `base` (parameter dict), `reps` independent simulations
    each. One row per combination x replicate, with the exact theory (`th_*`) next to the simulated numbers (`sim_*`)."""
    import itertools
    keys, rows, n = list(grid), [], 0
    for vals in itertools.product(*grid.values()):
        p = ProbeParams.from_dict({**base, **dict(zip(keys, vals))}, locus)
        th = {f"th_{k}": v for k, v in summary(p).items() if k in ("mean_attached", "cv_attached", "residence_s",
                                                                  "correlation_time_s", "p_none_attached")}
        for r in range(reps):
            n += 1
            rows.append({**dict(zip(keys, vals)), "rep": r, **th, **trace_stats(p, duration_s, seed * 100003 + n)})
    return pd.DataFrame(rows)
