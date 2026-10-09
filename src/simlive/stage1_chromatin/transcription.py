"""Stochastic transcription model whose rate can depend on the promoter-enhancer distance (stage 1c). Produces the MS2 signal.

Promoter (telegraph model) with the distance d(t) entering one chosen rate:

    OFF --k_on(d)--> ON --k_off(d)--> OFF           while ON, Pol II initiates at rate k_init(d)

Each initiation sends one Pol II along the MS2 cassette (5' end of the gene) and then the gene body at speed v. A Pol II carries
n_loops x (fraction of the cassette already transcribed) MS2 loops: its signal ramps up while it copies the cassette, stays at n_loops
while it travels through the gene body, and ends when it is released (after a random dwell at the end of the gene, mean `dwell_s`).
MS2 intensity = sum over Pol II, in units of "loops" (full signal of one Pol II = n_loops); the movie turns loops into photons.

Distance coupling (`coupling`): contact (sharp: 1 if d <= d_contact_nm), hill (smooth), exp (smooth), none (control: a constant f_const).
The coupling f(d) in [0, 1] multiplies the chosen rate by m = 1 + (fold - 1) f(d); k_off is DIVIDED by m (contact makes bursts last longer).
`coupled_rate` selects k_on (default), k_init or k_off.

All rates are per SECOND of real time and the distance series carries its own time step, so changing the movie's frame interval needs no
change here: the model runs in continuous time and is only SAMPLED at the frame times (with the exposure average the camera sees).
One-way coupling: the polymer does not feel transcription. Exact simulation (Gillespie, rates constant within a distance step); exact
expectation in `expected_*` (used to test the simulation and shown in the dashboard).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd
from numba import njit

COUPLINGS = ("contact", "hill", "exp", "none")
COUPLED_RATES = ("k_on", "k_init", "k_off")
# PLACEHOLDER values (the user's gene is not specified). Orders of magnitude from live-cell MS2 work on human genes, taken from
# search-result summaries of Larson-lab papers (full texts could not be opened from the build machine, so VERIFY before use):
# ON ~6 min and OFF ~45 min (MCF7 cells, Cell Reports 2021 MYC paper, Larson co-author); human Pol II elongation 1-6 kb/min across
# the literature (e.g. Darzacq 2007: ~1.9-4.3). The coupling strength, initiation rate and dwell are design choices, not literature values.
DEFAULTS = {
    "k_on": 1 / 3000.0, "k_off": 1 / 960.0, "k_init": 1 / 600.0,
    "coupling": "contact", "coupled_rate": "k_on", "fold": 12.0,
    "d_contact_nm": 150.0, "d_half_nm": 200.0, "hill_n": 4.0, "d_decay_nm": 150.0, "f_const": 0.1,
    "elongation_kb_min": 2.5, "cassette_kb": 1.3, "gene_kb": 5.0, "n_loops": 24, "dwell_s": 116.0,
}
HELP = {
    "k_on": "Rate (1/s) at which the OFF promoter turns ON, when the enhancer is far away. 1/k_on = mean OFF time.",
    "k_off": "Rate (1/s) at which the ON promoter turns OFF. 1/k_off = mean ON time (burst duration).",
    "k_init": "Rate (1/s) of Pol II initiation while the promoter is ON.",
    "fold": "How many times faster the coupled rate becomes when the enhancer is in full contact (k_off is divided by this).",
    "d_contact_nm": "Sharp contact: the enhancer 'touches' the promoter when their distance is at most this (nm).",
    "elongation_kb_min": "Pol II speed (kb per minute). Human genes: roughly 1-6 kb/min.",
    "cassette_kb": "Length of the MS2 cassette at the 5' end of the gene (kb). The signal ramps up while it is transcribed.",
    "gene_kb": "Length of the gene after the cassette (kb). The signal stays at full height while Pol II travels through it.",
    "n_loops": "Number of MS2 stem loops in the cassette (24 is common). Sets the signal of one Pol II.",
    "dwell_s": "Mean time (s) Pol II stays after reaching the end of the gene before it is released (random, exponential).",
}


@dataclass(frozen=True)
class TranscriptionParams:
    k_on: float = DEFAULTS["k_on"]
    k_off: float = DEFAULTS["k_off"]
    k_init: float = DEFAULTS["k_init"]
    coupling: str = DEFAULTS["coupling"]
    coupled_rate: str = DEFAULTS["coupled_rate"]
    fold: float = DEFAULTS["fold"]
    d_contact_nm: float = DEFAULTS["d_contact_nm"]
    d_half_nm: float = DEFAULTS["d_half_nm"]
    hill_n: float = DEFAULTS["hill_n"]
    d_decay_nm: float = DEFAULTS["d_decay_nm"]
    f_const: float = DEFAULTS["f_const"]
    elongation_kb_min: float = DEFAULTS["elongation_kb_min"]
    cassette_kb: float = DEFAULTS["cassette_kb"]
    gene_kb: float = DEFAULTS["gene_kb"]
    n_loops: int = DEFAULTS["n_loops"]
    dwell_s: float = DEFAULTS["dwell_s"]

    def __post_init__(self):
        if self.coupling not in COUPLINGS or self.coupled_rate not in COUPLED_RATES:
            raise ValueError(f"coupling must be one of {COUPLINGS}, coupled_rate one of {COUPLED_RATES}")
        if min(self.k_on, self.k_off, self.k_init, self.elongation_kb_min, self.dwell_s) <= 0 or self.fold < 1:
            raise ValueError("k_on, k_off, k_init, elongation_kb_min, dwell_s must be > 0 and fold >= 1")
        if min(self.d_contact_nm, self.d_half_nm, self.hill_n, self.d_decay_nm) <= 0 or not 0 <= self.f_const <= 1:
            raise ValueError("distance scales must be > 0 and f_const in [0, 1]")
        if self.cassette_kb < 0 or self.gene_kb < 0 or self.cassette_kb + self.gene_kb <= 0 or self.n_loops < 1:
            raise ValueError("gene lengths must be >= 0 (not both 0) and n_loops >= 1")

    @classmethod
    def from_dict(cls, d: dict | None) -> "TranscriptionParams":
        return cls(**{**DEFAULTS, **(d or {})})

    # derived times in seconds
    @property
    def t_cassette_s(self) -> float:
        return self.cassette_kb * 60.0 / self.elongation_kb_min

    @property
    def t_gene_s(self) -> float:                       # from initiation until Pol II reaches the end of the gene
        return (self.cassette_kb + self.gene_kb) * 60.0 / self.elongation_kb_min


def params_from_config(cfg: dict) -> TranscriptionParams:
    return TranscriptionParams.from_dict((cfg["loci"].get("transcription") or {}).get("params"))


# ----------------------------------------------------------------------------- coupling
def contact_fraction(p: TranscriptionParams, d_nm) -> np.ndarray:
    """f(d) in [0, 1]: how strongly the enhancer acts at distance d (nm)."""
    d = np.asarray(d_nm, float)
    if p.coupling == "contact":
        return (d <= p.d_contact_nm).astype(float)
    if p.coupling == "hill":
        return 1.0 / (1.0 + (d / p.d_half_nm) ** p.hill_n)
    if p.coupling == "exp":
        return np.exp(-d / p.d_decay_nm)
    return np.full(d.shape, p.f_const)


def rates_from_distance(p: TranscriptionParams, d_nm) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(k_on, k_off, k_init) in 1/s for each distance."""
    m = 1.0 + (p.fold - 1.0) * contact_fraction(p, d_nm)
    one = np.ones_like(m)
    kon = p.k_on * (m if p.coupled_rate == "k_on" else one)
    koff = p.k_off / (m if p.coupled_rate == "k_off" else one)
    kinit = p.k_init * (m if p.coupled_rate == "k_init" else one)
    return kon, koff, kinit


# ----------------------------------------------------------------------------- the profile of one Pol II
def _G_real(u, tc, Te, N):
    """Integral from 0 to u of ONE Pol II's signal: ramp over tc, plateau N, ends at Te (all seconds, arrays broadcast)."""
    w = np.clip(u, 0.0, Te)
    ramp = N * w * w / (2.0 * np.maximum(tc, 1e-12))
    return np.where(w <= tc, ramp, N * (w - tc / 2.0))


def _G_mean(u, p: TranscriptionParams):
    """Same, averaged over the random dwell (expected signal of one Pol II)."""
    tc, tg, tau, N = p.t_cassette_s, p.t_gene_s, p.dwell_s, p.n_loops
    u = np.maximum(u, 0.0)
    ramp = N * np.minimum(u, tc) ** 2 / (2.0 * max(tc, 1e-12))
    mid = np.where(u > tc, N * (np.minimum(u, tg) - tc), 0.0)
    tail = np.where(u > tg, N * tau * (1.0 - np.exp(-(u - tg) / tau)), 0.0)
    return ramp + mid + tail


# ----------------------------------------------------------------------------- exact simulation
@njit(cache=True)
def _telegraph(kon, koff, kinit, dt, seed, on0):
    """Switch times/states and initiation times of one promoter; rates are constant inside each step of length dt."""
    np.random.seed(seed)
    sw_t, sw_s, inits = [0.0], [on0], [0.0]
    inits.pop()
    state = on0
    for k in range(len(kon)):
        t, t1 = k * dt, (k + 1) * dt
        while True:
            r = kon[k] if state == 0 else koff[k] + kinit[k]
            t += -np.log(1.0 - np.random.random()) / r
            if t >= t1:
                break
            if state == 0:
                state = 1
                sw_t.append(t); sw_s.append(1)
            elif np.random.random() * r < koff[k]:
                state = 0
                sw_t.append(t); sw_s.append(0)
            else:
                inits.append(t)
    return np.array(sw_t), np.array(sw_s), np.array(inits)


@njit(cache=True)
def _G_one(u, tc, Te, N):
    w = min(max(u, 0.0), Te)
    return N * w * w / (2.0 * max(tc, 1e-12)) if w <= tc else N * (w - tc / 2.0)


@njit(cache=True)
def _sample_polymerases(ft, expo, inits, Te, tc, N):
    """Pol II count at each frame time and the MS2 signal averaged over [t, t + expo]. `inits` sorted ascending. Only Pol II that
    initiated within the longest transit before the frame can still be on the gene; earlier ones were finished before the window."""
    n = len(inits)
    n_pol, ms2 = np.zeros(len(ft), np.int64), np.zeros(len(ft))
    if n == 0:
        return n_pol, ms2
    tmax = Te.max()
    for i in range(len(ft)):
        a, b = ft[i], ft[i] + expo
        lo = np.searchsorted(inits, a - tmax, side="left")
        hi = np.searchsorted(inits, b, side="right")
        acc = 0.0
        for j in range(lo, hi):
            acc += _G_one(b - inits[j], tc, Te[j], N) - _G_one(a - inits[j], tc, Te[j], N)
            if inits[j] <= a and inits[j] + Te[j] > a:
                n_pol[i] += 1
        ms2[i] = acc / expo
    return n_pol, ms2


def pingpong(a: np.ndarray, n: int) -> tuple[np.ndarray, bool]:
    """Extend a to at least n points by appending alternately reversed and forward copies (no jump at the seams).
    Returns (array, True if it had to be extended)."""
    a = np.asarray(a)
    out, flip, looped = a, True, False
    while len(out) < n:
        out = np.concatenate([out, (a[::-1] if flip else a)[1:] if len(a) > 1 else a])
        flip, looped = not flip, True
    return out[:max(n, len(a))] if looped else out, looped


@dataclass
class CellTranscription:
    frames: pd.DataFrame                 # per frame: t, d_nm, in_contact, promoter_on, n_polII, ms2_loops (exposure average)
    events: pd.DataFrame                 # per Pol II: init_s, release_s (seconds from the cell's first frame; init may be < 0)
    switch_t: np.ndarray                 # promoter switch times (s from first frame, includes the pre-roll) ...
    switch_on: np.ndarray                # ... and the state after each switch
    d_nm: np.ndarray                     # distance series used, from t = 0 (one value per `dt_block_s`)
    dt_block_s: float
    looped: bool                         # True if the distance series was shorter than the movie and had to be repeated
    kon: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(0))


def simulate_cell(p: TranscriptionParams, d_nm, dt_block_s: float, frame_times_s, exposure_s: float,
                  seed: int) -> CellTranscription:
    """One cell: promoter switching, Pol II events and MS2 signal sampled at frame_times_s (seconds from the first frame).
    d_nm: promoter-enhancer distance, piecewise constant in steps of dt_block_s starting at t = 0."""
    rng = np.random.default_rng(seed)
    ft = np.asarray(frame_times_s, float)
    pre = p.t_gene_s + 5.0 * p.dwell_s                              # warm-up so Pol II already on the gene at t = 0 exist
    n_pre = int(np.ceil(pre / dt_block_s))
    n_need = int(np.ceil((ft.max() + exposure_s) / dt_block_s)) + 1
    d, looped = pingpong(np.asarray(d_nm, float), max(n_need, n_pre))
    d_ext = np.concatenate([d[:n_pre][::-1], d])                    # mirror image before t = 0 for the warm-up
    kon, koff, kinit = rates_from_distance(p, d_ext)
    on0 = int(rng.random() < kon[0] / (kon[0] + koff[0]))
    sw_t, sw_s, inits = _telegraph(kon, koff, kinit, dt_block_s, int(rng.integers(2**31)), on0)
    off = n_pre * dt_block_s
    sw_t, inits = sw_t - off, inits - off
    tc, tg, N = p.t_cassette_s, p.t_gene_s, p.n_loops
    Te = tg + rng.exponential(p.dwell_s, len(inits))                # time from initiation to release, per Pol II
    release = inits + Te
    # sample at the frames: instantaneous counts and the exposure-averaged MS2 signal
    on = sw_s[np.clip(np.searchsorted(sw_t, ft, side="right") - 1, 0, len(sw_t) - 1)]
    n_pol, ms2 = _sample_polymerases(ft, exposure_s, inits, Te, tc, float(N))
    dk = np.clip(np.floor(ft / dt_block_s).astype(int), 0, len(d) - 1)
    frames = pd.DataFrame({"t_s": ft, "d_nm": d[dk], "in_contact": contact_fraction(p, d[dk]) > 0.5 if p.coupling != "none" else False,
                           "promoter_on": on.astype(int), "n_polII": n_pol, "ms2_loops": ms2})
    keep = release > 0
    return CellTranscription(frames, pd.DataFrame({"init_s": inits[keep], "release_s": release[keep]}), sw_t, sw_s,
                             d[:n_need], dt_block_s, looped, kon[n_pre:n_pre + n_need])


# ----------------------------------------------------------------------------- exact expectations (theory)
def expected_on(p: TranscriptionParams, d_nm, dt_block_s: float, p0: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """P(promoter ON) at the start of every step and its mean over each step, for the GIVEN distance series (exact:
    rates are constant per step). p0 = starting probability (default: stationary for the first step)."""
    kon, koff, _ = rates_from_distance(p, d_nm)
    lam = kon + koff
    pinf = kon / lam
    p_start = np.empty(len(kon) + 1)
    p_start[0] = pinf[0] if p0 is None else p0
    decay = np.exp(-lam * dt_block_s)
    for k in range(len(kon)):
        p_start[k + 1] = pinf[k] + (p_start[k] - pinf[k]) * decay[k]
    mean = pinf + (p_start[:-1] - pinf) * (1.0 - decay) / (lam * dt_block_s)
    return p_start, mean


def expected_ms2(p: TranscriptionParams, d_nm, dt_block_s: float, times_s) -> np.ndarray:
    """Expected instantaneous MS2 signal (loops) at times_s for a distance series starting at t = 0 (promoter stationary at the start,
    no Pol II before t = 0, so only meaningful after about t_gene_s + dwell_s)."""
    _, pmean = expected_on(p, d_nm, dt_block_s)
    _, _, kinit = rates_from_distance(p, d_nm)
    r = kinit * pmean                                               # expected initiations per second in each step
    s = np.arange(len(r) + 1) * dt_block_s
    t = np.asarray(times_s, float)[:, None]
    return ((_G_mean(t - s[None, :-1], p) - _G_mean(t - s[None, 1:], p)) * r[None, :]).sum(1)


def summary(p: TranscriptionParams, d_nm: float) -> dict:
    """Exact stationary numbers if the distance stayed fixed at d_nm."""
    kon, koff, kinit = (float(x[0]) for x in rates_from_distance(p, [d_nm]))
    on = kon / (kon + koff)
    residence = p.t_gene_s + p.dwell_s
    ramp_area = p.n_loops * (residence - p.t_cassette_s / 2.0)       # integral of one Pol II's signal over its life
    return {"k_on": kon, "k_off": koff, "k_init": kinit, "on_fraction": on, "mean_off_min": 1 / kon / 60, "mean_on_min": 1 / koff / 60,
            "bursts_per_hour": 3600.0 * kon * koff / (kon + koff), "pol_per_burst": kinit / koff,
            "mean_polII": kinit * on * residence, "mean_ms2_loops": kinit * on * ramp_area,
            "transit_min": p.t_gene_s / 60.0, "residence_s": residence}


# ----------------------------------------------------------------------------- stand-in distance process (no library yet)
def standin_distance(n: int, dt_s: float, rng: np.random.Generator, contact_fraction_: float = 0.1,
                     loop_time_s: float = 240.0, open_nm: float = 350.0, contact_nm: float = 80.0,
                     relax_s: float = 30.0) -> np.ndarray:
    """INVENTED promoter-enhancer distance (nm): a hidden looped/open switch (looped `contact_fraction_` of the time, mean looped
    time `loop_time_s`) with confined fluctuations around `contact_nm` / `open_nm`. NOT from polymer simulations or data."""
    k_out = 1.0 / loop_time_s
    k_in = k_out * contact_fraction_ / (1.0 - contact_fraction_)
    state, a = rng.random() < contact_fraction_, np.exp(-dt_s / relax_s)
    x, out = rng.normal(), np.empty(n)
    for i in range(n):
        if rng.random() < (1 - np.exp(-(k_out if state else k_in) * dt_s)):
            state = not state
        x = a * x + np.sqrt(1 - a * a) * rng.normal()
        out[i] = abs(contact_nm * (1 + 0.3 * x)) if state else open_nm * np.exp(0.35 * x)
    return out


# ----------------------------------------------------------------------------- the movie pipeline
def simulate_transcription(cells: pd.DataFrame, traces: dict, cfg: dict, rng: np.random.Generator):
    """Run the model for every cell. `traces`: {cell_id: {"dt_block_s": float, "d_nm": array}} (distance series starting at the
    cell's first frame). Returns (frames table, events table); frame rows carry t (frame index), cell_id and the model state."""
    dt, expo = float(cfg["acquisition"]["frame_interval_s"]), float(cfg["acquisition"]["exposure_s"])
    p = params_from_config(cfg)
    fr, ev = [], []
    for cid, g in cells.sort_values("t").groupby("cell_id"):
        ts = g.t.to_numpy()
        tr = traces[int(cid)]
        c = simulate_cell(p, tr["d_nm"], tr["dt_block_s"], (ts - ts.min()) * dt, expo, int(rng.integers(2**31)))
        f = c.frames.copy()
        f.insert(0, "t", ts); f.insert(1, "cell_id", cid); f["looped_distance"] = c.looped
        e = c.events.copy(); e.insert(0, "cell_id", cid)
        fr.append(f); ev.append(e)
    return pd.concat(fr, ignore_index=True), pd.concat(ev, ignore_index=True)


def traces_from_truth(truth: pd.DataFrame, cfg: dict) -> dict:
    """Distance series for the toy (non-polymer) loci: 2D lab-frame distance between locus 0 and locus 1 at every frame (nm),
    piecewise constant over a frame. Used only when loci come from `toy_loci` (no polymer simulation to read from)."""
    dt, out = float(cfg["acquisition"]["frame_interval_s"]), {}
    for cid, g in truth.groupby("cell_id"):
        a = g[g.locus_id == 0].sort_values("t"); b = g[g.locus_id == 1].sort_values("t")
        out[int(cid)] = {"dt_block_s": dt, "looped": False,
                         "d_nm": np.hypot(a.y_um.to_numpy() - b.y_um.to_numpy(), a.x_um.to_numpy() - b.x_um.to_numpy()) * 1000.0}
    return out


def ms2_spots(tx: pd.DataFrame, truth: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Where and how bright the MS2 spots are: at the promoter (locus 0), photons = photons_per_loop_s x exposure x MS2 loops."""
    tc = cfg["loci"]["transcription"]
    prom = truth[truth.locus_id == 0][["t", "cell_id", "y_um", "x_um"]]
    m = prom.merge(tx[["t", "cell_id", "ms2_loops"]], on=["t", "cell_id"])
    m["photons"] = float(tc.get("photons_per_loop_s", 100.0)) * float(cfg["acquisition"]["exposure_s"]) * m.ms2_loops
    return m[["t", "cell_id", "y_um", "x_um", "photons"]]


def check_sampling(p: TranscriptionParams, frame_interval_s: float, exposure_s: float) -> list[tuple[str, str]]:
    """Plain-language warnings about whether the frame interval can resolve this gene's dynamics: [(level, message)]."""
    out, s = [], summary(p, 1e9)
    far_on, near_on = 1 / p.k_off, 1 / (p.k_off / (p.fold if p.coupled_rate == "k_off" else 1.0))
    ton = min(far_on, near_on)
    if frame_interval_s > ton / 2:
        out.append(("bad", f"Frame interval {frame_interval_s:g} s is longer than half the shortest ON time ({ton:.0f} s): bursts "
                           "will be missed or merged."))
    elif frame_interval_s > ton / 6:
        out.append(("warn", f"Only {ton / frame_interval_s:.0f} frames per ON period: burst durations will be coarse."))
    if frame_interval_s > p.t_cassette_s:
        out.append(("warn", f"Frame interval is longer than the time Pol II needs for the MS2 cassette ({p.t_cassette_s:.0f} s): "
                            "the ramp-up of the signal cannot be seen."))
    if exposure_s > frame_interval_s:
        out.append(("bad", "Exposure is longer than the frame interval."))
    if not out:
        out.append(("good", f"{frame_interval_s:g} s per frame resolves the dynamics: {ton / frame_interval_s:.0f} frames per shortest ON "
                            f"period, {p.t_gene_s / frame_interval_s:.0f} frames per Pol II transit."))
    return out
