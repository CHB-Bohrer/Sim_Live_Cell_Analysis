"""Stage 6: spot detection and tracking after Chen, Levo, Zoller & Gregor 2025 ("A conserved coupling of transcriptional ON and OFF periods
underlies bursting dynamics", Nat. Struct. Mol. Biol., doi 10.1038/s41594-025-01615-4; methods from the arXiv version 2304.08770v3, section
III.2 "Transcription spot tracking and signal integration"; no code was published, so this follows the text).

Their pipeline (MS2 transcription spots, one labelled allele per nucleus, nuclear MCP-GFP background, 220 nm pixels, 12 planes 1 um apart):
    1  split the movie into one 3D stack per nucleus, tracked over time            -> stage 5's single-cell crops
    2  3D median filter against salt-and-pepper noise                               -> `median_px` (x `median_z` in z-stacks)
    3  3D difference-of-Gaussians to find round objects; "tuned to high
       sensitivity", so several candidates per nucleus and time point                -> `dog_sigma_px` (inner, outer), `dog_snr`, `max_candidates`
    4  "at most one spot per nucleus at all times": a modified Viterbi algorithm
       picks the most likely trajectory through the candidates, scoring each link
       by pairwise distance with a diffusion-based potential "constructed
       iteratively from the random walk of all transcription spots"                 -> `viterbi_track`
    5  intensity = mean in an ellipsoid of diameters (7, 7, 5) px at the spot
       centroid minus the mean of a (19, 19, 7) px ellipsoid around it (shell
       only), removing free MCP-GFP                                                  -> `ellipsoid_intensity`, `intensity_inner_px/outer_px`

What the text leaves open, and the choices made here (each a setting):
- DoG widths: inner ~ the PSF sigma (1 px here), outer 3x (a band-pass that also removes the slowly varying background).
- The Viterbi score. Nodes = candidates; a frame may have no spot (OFF period, too dim): links may skip up to `viterbi_max_gap` frames at
  `viterbi_miss_cost` per skipped frame. Link cost = -log of the diffusion propagator, |d|^2 / (4 D g) + (ndim/2) log(4 pi D g) for a gap of g
  frames; node reward = log(SNR of the candidate) + log(region size / false candidates per frame) (the clutter term of a track
  likelihood: a candidate left off the path must be explained as a false one), which prefers real spots over noise peaks without letting a
  very bright one-frame event (cosmic ray) pull the path away. Positions are relative to the nucleus centroid (the crop follows the nucleus).
- The potential: D (px^2 per frame) is re-estimated from the chosen path's own steps (robust: median of |d|^2 / g) and the path recomputed,
  `viterbi_iterations` times. Here per cell and colour (the paper pools all spots of a movie); `viterbi_D_px2` fixes it instead.
- Position: the paper uses "the centroid coordinates" of the spot (`refine: centroid`); the default here is the Gaussian fit of spots.py,
  which measured better (docs/SPOT_FINDING.md).
- Intensity background: the paper averages the shell; the default here takes its median (`intensity_background`), so a cosmic ray or a
  neighbour's spot in the shell does not make the intensity negative (one such frame dropped the intensity-photon correlation 0.87 -> 0.75).
- Ellipsoid sizes are in pixels as in the paper. At 220 nm their inner 7 px diameter is ~4.6 PSF sigmas wide; at our 130 nm and ~0.85 px PSF
  sigma, 7 px is ~4.1 sigmas, so the same numbers fit; check them for other pixel sizes.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, median_filter

from simlive.stage6_analysis.spots import _reg, robust_bg_noise

MEDIAN_CHI2 = {1: 0.4549, 2: 1.3863, 3: 2.3660}       # median of chi-square with 1, 2, 3 degrees of freedom


def _pair(v, default):
    if v is None:
        return default
    v = list(v) if isinstance(v, (list, tuple)) else [v, 3 * float(v)]
    return float(v[0]), float(v[1])


def dog_candidates(img: np.ndarray, region: np.ndarray, *, median_px: int = 3, median_z=1, dog_sigma_px=(1.0, 3.0), dog_sigma_z=(0.8, 2.0),
                   dog_snr: float = 3.0, max_candidates: int = 5, min_sep_px: int = 4, steps: bool = False, **_) -> dict | None:
    """Steps 2-3 on one image (y, x) or stack (z, y, x) of one cell (`region`: its 2D mask). Returns None for a region too small, else
    {img, pre: {bg, sub, noise_smooth}, peaks (brightest first), values, snr} in the same form as spots.candidates (so spots.refine_peak works).
    `bg` = the outer Gaussian of the median-filtered image (the local background), `sub` = the DoG response."""
    from skimage.feature import peak_local_max

    img = np.asarray(img, np.float32)
    three_d = img.ndim == 3
    reg = _reg(region, img.shape)
    if reg.sum() < 20:
        return None
    k = max(int(median_px), 1)
    kz = int(median_z) if str(median_z).isdigit() else 1
    size = (max(kz, 1), k, k) if three_d else (k, k)
    med = median_filter(img, size=size, mode="nearest") if max(size) > 1 else img
    s1, s2 = _pair(dog_sigma_px, (1.0, 3.0))
    if three_d:
        z1, z2 = _pair(dog_sigma_z, (0.8, 2.0))
        inner, outer = gaussian_filter(med, (z1, s1, s1)), gaussian_filter(med, (z2, s2, s2))
    else:
        inner, outer = gaussian_filter(med, s1), gaussian_filter(med, s2)
    dog = (inner - outer).astype(np.float32)
    noise = robust_bg_noise(dog[reg])[1]
    work = np.where(reg, dog, 0.0)
    peaks = peak_local_max(work, min_distance=max(int(min_sep_px), 1), threshold_abs=float(dog_snr) * noise, exclude_border=False,
                           num_peaks=int(max_candidates))
    peaks = np.asarray(peaks, int).reshape(-1, img.ndim)
    peaks = peaks[reg[tuple(peaks.T)]] if len(peaks) else peaks
    vals = dog[tuple(peaks.T)] if len(peaks) else np.zeros(0, np.float32)
    order = np.argsort(vals)[::-1]
    out = {"img": img, "pre": {"bg": outer, "sub": dog, "noise_smooth": noise, "med": med}, "peaks": peaks[order], "values": vals[order],
           "snr": vals[order] / noise}
    if steps:
        out.update(thr=float(dog_snr) * noise)
    return out


def _link_cost(a: np.ndarray, b: np.ndarray, D: float, g: int, ndim: int) -> np.ndarray:
    """-log of the diffusion propagator for every pair (a_i -> b_j) over g frames (positions in px; D in px^2 per frame)."""
    d2 = ((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)
    return d2 / (4 * D * g) + 0.5 * ndim * np.log(4 * np.pi * D * g)


INT_Z2_CAP = 16.0           # brightness z^2 of one spot is capped (a real locus frame can be odd: an overlapping stray, a cosmic ray)


class IntensityModel:
    """The brightness of the spots a path follows, as a Gaussian AR(1) process in the DoG SNR (option `viterbi_intensity`; user 2026-10-26:
    "utilize the intensity of the individual loci ... for we explicitly simulate the intensity of the loci using Gillespie"):
        I_t = mu(t) + e_t,   e_t = rho e_{t-1} + innovation,   Var(e) = s^2,
    mu(t) a straight line (bleaching), fitted with s (robust) and rho (lag-1 correlation of consecutive spots) to the path's own spots, like D.
    Measurement noise is part of s and lowers rho, so the fit includes it. Every spot then pays (z^2 - 1) / 2 nats on top of the usual costs,
    z = (I - predicted) / predicted sd, predicted from the path's previous spot (gap g: mean mu + rho^g e_prev, sd s sqrt(1 - rho^2g)) or,
    at a start, from the level alone. A genuine spot costs nothing on average; a spot much dimmer (stray spot, noise peak) or much brighter
    (cosmic ray) than the path pays. The other copy of a two-copy cell is about as bright, so brightness neither favours nor penalises it
    unless rho is high and its brightness differs from the path's recent brightness; position decides as before.
    A first version scored log p_path(I) - log p_clutter(I) instead of log SNR, with the clutter fitted to the candidates the path leaves. In
    two-copy cells those include the other copy (locus-bright), so a real spot's brightness looked like clutter: on 12 sweep movies the share
    of loci found fell from 98% to 92-96% and other-copy picks rose 0 -> 0.8-1.5%. Dropped."""

    def __init__(self, t: np.ndarray, I: np.ndarray):
        t, I = np.asarray(t, float), np.asarray(I, float)
        if len(I) >= 6 and np.ptp(t) > 0:
            b, a = np.polyfit(t, I, 1)
        else:
            a, b = float(np.median(I)) if len(I) else 1.0, 0.0
        self.a, self.b = float(a), float(b)
        res = I - self.mu(t)
        mad = float(np.median(np.abs(res - np.median(res)))) * 1.4826 if len(res) else 0.0
        self.s = max(mad, 0.05 * max(abs(self.a), 1.0), 0.3)
        ok = np.diff(t) == 1
        if ok.sum() >= 5:
            r0, r1 = res[:-1][ok], res[1:][ok]
            self.rho = float(np.clip(np.sum(r0 * r1) / max(np.sqrt(np.sum(r0 ** 2) * np.sum(r1 ** 2)), 1e-9), 0.0, 0.95))
        else:
            self.rho = 0.0

    def mu(self, t):
        return self.a + self.b * np.asarray(t, float)

    def start_cost(self, t: int, I: np.ndarray) -> np.ndarray:
        z2 = ((np.asarray(I, float) - self.mu(t)) / self.s) ** 2
        return 0.5 * (np.minimum(z2, INT_Z2_CAP) - 1.0)

    def link_cost(self, s: int, Is: np.ndarray, t: int, It: np.ndarray) -> np.ndarray:
        """(z^2 - 1) / 2 of I_t = It[j] given I_s = Is[i], for every pair (gap t - s)."""
        r = self.rho ** (t - s)
        m = self.mu(t) + r * (np.asarray(Is, float)[:, None] - self.mu(s))
        sd = self.s * np.sqrt(max(1.0 - r * r, 1e-3))
        z2 = ((np.asarray(It, float)[None, :] - m) / sd) ** 2
        return 0.5 * (np.minimum(z2, INT_Z2_CAP) - 1.0)

    def as_dict(self) -> dict:
        return {"a": self.a, "b": self.b, "s": self.s, "rho": self.rho}


def _viterbi_once(rel: list, reward: list, D: float, max_gap: int, miss_cost: float, reacquire: float | None = None,
                  inten: tuple | None = None) -> list:
    """Lowest-cost path through the candidates (at most one per frame; frames may be skipped). Gaps of up to `max_gap` frames are bridged by
    the diffusion link; after a longer gap the path may pick the spot up again anywhere, at `reacquire` (-log of a uniform position over the
    cell, i.e. log of its size) plus the per-frame miss cost, so one path can hold several bursts separated by long OFF periods.
    reacquire=None: no re-acquisition (one connected stretch). `inten` = (IntensityModel, SNR per frame): the brightness of a spot is then
    scored given the path's previous spot (link) or the path's level (start / re-acquisition), see IntensityModel.
    Returns the index chosen per frame or None."""
    model, I = inten if inten is not None else (None, None)
    n = len(rel)
    ndim = next((r.shape[1] for r in rel if len(r)), 2)
    best = [np.full(len(r), np.inf) for r in rel]
    back = [[None] * len(r) for r in rel]
    run_min, run_arg = np.inf, None             # min over s <= t - max_gap - 1 of (best end cost at s - miss_cost * s)
    for t in range(n):
        s_new = t - max_gap - 1
        if reacquire is not None and s_new >= 0 and len(rel[s_new]) and np.isfinite(best[s_new]).any():
            j0 = int(np.argmin(best[s_new]))
            v = best[s_new][j0] - miss_cost * s_new
            if v < run_min:
                run_min, run_arg = v, (s_new, j0)
        if not len(rel[t]):
            continue
        start = model.start_cost(t, I[t]) if model is not None else 0.0
        best[t] = miss_cost * t - reward[t] + start                         # start here (skipping the frames before)
        if run_arg is not None:                                             # pick the spot up again after a long gap
            c = run_min + miss_cost * (t - 1) + reacquire - reward[t] + start
            better = c < best[t]
            best[t] = np.where(better, c, best[t])
            for j in np.nonzero(better)[0]:
                back[t][j] = run_arg
        for g in range(1, max_gap + 1):
            s = t - g
            if s < 0:
                break
            if not len(rel[s]) or not np.isfinite(best[s]).any():
                continue
            c = best[s][:, None] + _link_cost(rel[s], rel[t], D, g, ndim) + miss_cost * (g - 1) - reward[t][None, :]
            if model is not None:
                c = c + model.link_cost(s, I[s], t, I[t])
            i = np.argmin(c, axis=0)
            cj = c[i, np.arange(len(rel[t]))]
            better = cj < best[t]
            best[t] = np.where(better, cj, best[t])
            for j in np.nonzero(better)[0]:
                back[t][j] = (s, int(i[j]))
    ends = [(best[t][j] + miss_cost * (n - 1 - t), t, j) for t in range(n) for j in range(len(rel[t])) if np.isfinite(best[t][j])]
    choice = [None] * n
    if not ends:
        return choice
    _, t, j = min(ends)
    while True:
        choice[t] = j
        prev = back[t][j]
        if prev is None:
            break
        t, j = prev
    return choice


def estimate_D(rel: list, choice: list, D_min: float = 0.25, max_gap: int = 5) -> float:
    """Diffusion potential from the path's own steps: D = median(|d|^2 / g) / (2 median chi2_ndim), robust to a few wrong links. Only links
    over gaps of at most `max_gap` frames count (re-acquisitions after long gaps say nothing about the motion)."""
    pts = [(t, rel[t][j]) for t, j in enumerate(choice) if j is not None]
    steps = [((b - a) ** 2).sum() / (tb - ta) for (ta, a), (tb, b) in zip(pts[:-1], pts[1:]) if tb - ta <= max_gap]
    if len(steps) < 2:
        return max(1.0, D_min)
    ndim = len(pts[0][1])
    return max(float(np.median(steps)) / (2 * MEDIAN_CHI2[ndim]), D_min)


def viterbi_track(peaks: list, snr: list, centre: np.ndarray, *, viterbi_max_gap: int = 5, viterbi_miss_cost: float = 2.0,
                  viterbi_D_px2: float | None = None, viterbi_iterations: int = 3, z_scale: float = 1.0, search_volume: float | None = None,
                  reward_offset: list | None = None, viterbi_intensity: bool = False, model_out: dict | None = None, **_) -> tuple[list, float]:
    """Step 4. `peaks[i]`: candidate positions of frame i (crop px; (z, y, x) for stacks, z in planes), `snr[i]` their DoG signal-to-noise,
    `centre[i]` the nucleus centroid (y, x) in crop px; `z_scale` = plane spacing / pixel size, so distances are isotropic;
    `search_volume` = area (px^2) or volume (px^3, z scaled) of the cell region where candidates can appear; `reward_offset[i]` (optional,
    one value per candidate of frame i, nats) is added to the node rewards, e.g. the pair gate of multicopy.py (a candidate far from the
    other colour's spot of the same copy is penalised).
    A candidate left off the path must be explained as a false one, which appear uniformly over the region at (mean candidates per frame
    - 1) per frame: each node's reward therefore includes log(search_volume / false candidates per frame), the usual clutter term of a
    track likelihood. Without it, skipping a frame is cheaper than keeping a well-placed but faint spot.
    viterbi_intensity: after the usual passes, the brightness of the path's spots is modelled (`IntensityModel`) and the path found again
    with a brightness-consistency cost added to every spot, D and the model re-fitted each time
    (`viterbi_iterations` passes). `model_out` (a dict), if given, receives the final model's numbers.
    Returns (chosen candidate per frame or None, D in px^2 per frame)."""
    rel = []
    for i, p in enumerate(peaks):
        p = np.asarray(p, float).reshape(-1, np.asarray(p).shape[-1] if len(p) else 2)
        if not len(p):
            rel.append(np.zeros((0, 2)))
            continue
        q = p.copy()
        q[:, -2:] -= centre[i]
        if q.shape[1] == 3:
            q[:, 0] *= z_scale
        rel.append(q)
    counts = [len(r) for r in rel]
    false_per_frame = max(float(np.mean(counts)) - 1.0, 0.1) if counts else 1.0
    clutter = np.log(max(float(search_volume), 1.0) / false_per_frame) if search_volume else 0.0
    reward = [np.log(np.maximum(np.asarray(s, float), 1.0001)) + clutter for s in snr]
    if reward_offset is not None:
        reward = [r + np.asarray(o, float) for r, o in zip(reward, reward_offset)]
    reacq = float(np.log(max(float(search_volume), 1.0))) if search_volume else None
    gap, miss = int(viterbi_max_gap), float(viterbi_miss_cost)
    if viterbi_D_px2:
        D = float(viterbi_D_px2)
        choice = _viterbi_once(rel, reward, D, gap, miss, reacq)
    else:
        brightest = [0 if len(r) else None for r in rel]
        D = estimate_D(rel, brightest, max_gap=gap)
        choice = brightest
        for _ in range(max(int(viterbi_iterations), 1)):
            choice = _viterbi_once(rel, reward, D, gap, miss, reacq)
            D_new = estimate_D(rel, choice, max_gap=gap)
            if abs(D_new - D) < 0.02 * D:
                D = D_new
                break
            D = D_new
    if not viterbi_intensity or sum(c is not None for c in choice) < 6:
        return choice, D
    I = [np.asarray(v, float) for v in snr]
    for _ in range(max(int(viterbi_iterations), 1)):
        on = [(t, j) for t, j in enumerate(choice) if j is not None]
        model = IntensityModel(np.array([t for t, _ in on]), np.array([I[t][j] for t, j in on]))
        new = _viterbi_once(rel, reward, D, gap, miss, reacq, inten=(model, I))
        if not viterbi_D_px2:
            D = estimate_D(rel, new, max_gap=gap)
        if new == choice:
            break
        choice = new
    if model_out is not None:
        model_out.update(model.as_dict(), D_px2=D)
    return choice, D


def ellipsoid_intensity(img: np.ndarray, pos, inner=(7, 7, 5), outer=(19, 19, 7), background: str = "mean") -> tuple[float, int]:
    """Step 5: mean of the pixels inside an ellipsoid of diameters `inner` (x, y[, z] in px / planes) centred at `pos` ((y, x) or (z, y, x))
    minus the mean (Chen et al.) or median (`background="median"`: robust to a cosmic ray or another spot in the shell) of the pixels inside the
    larger `outer` ellipsoid but outside the inner one (`background="none"`: no subtraction). Returns (difference, number of inner pixels)."""
    img = np.asarray(img, np.float32)
    nd = img.ndim
    pos = np.asarray(pos, float)
    ri = np.array([inner[2] if nd == 3 else 0, inner[1], inner[0]][-nd:], float) / 2
    ro = np.array([outer[2] if nd == 3 else 0, outer[1], outer[0]][-nd:], float) / 2
    lo = np.maximum(np.floor(pos - ro).astype(int), 0)
    hi = np.minimum(np.ceil(pos + ro).astype(int) + 1, img.shape)
    sl = tuple(slice(a, b) for a, b in zip(lo, hi))
    grids = np.meshgrid(*[np.arange(a, b) for a, b in zip(lo, hi)], indexing="ij")
    ui = sum(((g - c) / r) ** 2 for g, c, r in zip(grids, pos, np.maximum(ri, 1e-6)))
    uo = sum(((g - c) / r) ** 2 for g, c, r in zip(grids, pos, np.maximum(ro, 1e-6)))
    win = img[sl]
    m_in, m_sh = ui <= 1.0, (uo <= 1.0) & (ui > 1.0)
    if not m_in.any() or not m_sh.any():
        return float("nan"), int(m_in.sum())
    if background == "none":            # no local background: the raw mean (only to show what the shell subtraction removes)
        return float(win[m_in].mean()), int(m_in.sum())
    bg = float(np.median(win[m_sh])) if background == "median" else float(win[m_sh].mean())
    return float(win[m_in].mean() - bg), int(m_in.sum())
