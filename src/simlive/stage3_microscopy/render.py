"""Stage 3: the microscope forward model. Ground truth (cells.csv, loci, MS2 spots, probe occupancy) -> recorded images.

Channels: `nucleus` (NLS-GFP; also carries the MS2 spots, same colour) and one channel per locus (`locus0`, `locus1`, ...).
Ground truth for the images: `labels` (integer mask per frame, value = persistent cell ID, 0 = background), the true locus positions,
and `imaging_truth` (what the sample position / focus / illumination actually were in each frame).

How a frame is made (see errors.py for every source and its switch): the objects are drawn at their true positions (shifted by stage
drift and, per colour, by chromatic shift) -> optical blur (diffraction, defocus, exposure motion, out-of-focus haze) -> photobleaching
-> uneven illumination and flicker -> colour cross-talk -> background (constant, autofluorescence) -> detector (quantum efficiency,
shot noise, dark current) -> camera electronics (hot pixels, cosmic rays, pixel gain, read noise, row noise, offset pattern, digitisation).

`Imager` renders any single frame on its own, with any subset of error sources switched on; the dashboard uses that to show how an
image is built up and which ground truth produced it. `render_nuclei` / `render_loci` are thin wrappers for one channel.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

from simlive.stage3_microscopy import errors as ERR
from simlive.stage3_microscopy.errors import ErrorSet


def psf_sigma_px(optics: dict) -> float:
    """Gaussian approximation to the widefield Airy PSF: sigma ~ 0.21 * lambda / NA."""
    return 0.21 * optics["wavelength_nm"] / optics["NA"] / optics["pixel_size_nm"]


def entropy_for(seed: int) -> int:
    """Imaging entropy of a run (the rendering stream of the pipeline's SeedSequence(seed).spawn(6))."""
    return int(np.random.SeedSequence(int(seed)).spawn(6)[1].generate_state(1, np.uint64)[0])


def _nucleus_mask(dy, dx, r, px_um, harm):
    """Pixels inside a deformed ellipse. dy, dx: pixel offsets from the cell centre; r: a cells.csv row."""
    from simlive.stage2_cells.motion import boundary_rho

    c, s = np.cos(r.angle_rad), np.sin(r.angle_rad)
    u, v = dy * s + dx * c, dy * c - dx * s        # rotate into the ellipse frame (u along the long axis)
    a_ax, b_ax = np.sqrt(r.aspect), 1 / np.sqrt(r.aspect)
    uu, vv = u / a_ax, v / b_ax
    amp = np.array([getattr(r, f"amp{int(k)}") for k in harm])
    phase = np.array([getattr(r, f"phase{int(k)}") for k in harm])
    r_px = r.radius_um / px_um
    inside = np.hypot(uu, vv) <= r_px * boundary_rho(np.arctan2(vv, uu), amp, phase, harm)
    return inside, uu / r_px, vv / r_px  # also the position in the nucleus' own frame (units of its radius)


_TEX_GRID, _TEX_EXTENT = 96, 1.5


def make_texture(tex: dict, radius_um: float, rng: np.random.Generator) -> np.ndarray:
    """Per-nucleus intensity texture on a grid covering +-1.5 nuclear radii, mean ~1 inside the nucleus.

    Fixed to the nucleus (rotates/deforms with it), so it is a feature a tracker can use. Components:
    spatially correlated chromatin-like variation (`contrast`, `correlation_um`) and dark nucleoli where nuclear
    NLS-GFP is excluded (`n_nucleoli` on average, `nucleolus_depth`).
    """
    g = _TEX_GRID
    cell = 2 * _TEX_EXTENT / g                          # grid spacing in nuclear radii
    field = gaussian_filter(rng.normal(size=(g, g)), (tex["correlation_um"] / radius_um) / cell, mode="wrap")
    field /= field.std()
    t = np.clip(1 + tex["contrast"] * field, 0.25, None)
    ax = np.linspace(-_TEX_EXTENT, _TEX_EXTENT, g)
    gx, gy = np.meshgrid(ax, ax)
    for _ in range(min(int(rng.poisson(tex.get("n_nucleoli", 0))), 6)):
        ang, rad, s = rng.uniform(0, 2 * np.pi), 0.65 * np.sqrt(rng.uniform()), rng.uniform(0.10, 0.18)
        d2 = (gx - rad * np.cos(ang)) ** 2 + (gy - rad * np.sin(ang)) ** 2
        t *= 1 - tex.get("nucleolus_depth", 0.5) * np.exp(-d2 / (2 * s**2))
    return t / t[np.hypot(gx, gy) <= 1].mean()


def _sample_texture(tex_map: np.ndarray, xn: np.ndarray, yn: np.ndarray) -> np.ndarray:
    from scipy.ndimage import map_coordinates

    to_grid = lambda c: (c + _TEX_EXTENT) / (2 * _TEX_EXTENT) * (_TEX_GRID - 1)
    return map_coordinates(tex_map, [to_grid(yn), to_grid(xn)], order=1, mode="nearest")


def _spots_from_arrays(shape, y, x, photons, sig: float) -> np.ndarray:
    """Expected photons per pixel from point sources at pixel coordinates (y, x): Gaussian PSF integrated over each pixel."""
    from scipy.special import erf

    ny, nx = shape
    out = np.zeros(shape, np.float32)
    sig = max(float(sig), 0.3)
    half, s2 = int(np.ceil(4 * sig)) + 1, np.sqrt(2) * sig
    for cy, cx, ph in zip(y, x, photons):
        y0, y1 = max(int(round(cy)) - half, 0), min(int(round(cy)) + half + 1, ny)
        x0, x1 = max(int(round(cx)) - half, 0), min(int(round(cx)) + half + 1, nx)
        if y1 <= y0 or x1 <= x0:
            continue
        ys, xs = np.arange(y0, y1), np.arange(x0, x1)
        py = 0.5 * (erf((ys + 0.5 - cy) / s2) - erf((ys - 0.5 - cy) / s2))
        px = 0.5 * (erf((xs + 0.5 - cx) / s2) - erf((xs - 0.5 - cx) / s2))
        out[y0:y1, x0:x1] += ph * np.outer(py, px)
    return out


def _spot_photons(shape, spots: pd.DataFrame, sig: float) -> np.ndarray:
    """Expected photons per pixel from point sources (columns y_px, x_px, photons)."""
    return _spots_from_arrays(shape, spots.y_px.to_numpy(), spots.x_px.to_numpy(), spots.photons.to_numpy(), sig)


class _NucleusModel:
    """The nuclei: brightness per cell, a fixed texture per nucleus, the masks and (optionally) MS2 spots."""

    def __init__(self, cells: pd.DataFrame, cfg: dict, entropy: int, spots: pd.DataFrame | None):
        self.cfg, self.entropy = cfg, entropy
        self.opt = cfg["optics"]
        self.px_um = self.opt["pixel_size_nm"] / 1000.0
        self.shape = tuple(int(round(v / self.px_um)) for v in cfg["geometry"]["fov_um"])
        self.yy, self.xx = np.mgrid[:self.shape[0], :self.shape[1]]
        self.harm = np.array([int(c[3:]) for c in cells.columns if c.startswith("amp")], float)
        self.by_t = {int(t): g for t, g in cells.groupby("t")}
        tex = cfg.get("nucleus_texture")
        if tex and tex.get("contrast", 0) == 0 and tex.get("n_nucleoli", 0) == 0:
            tex = None
        self.tex_cfg, self.textures, self.brightness = tex, {}, {}
        self.spots_by_t = {}
        if spots is not None and len(spots):
            for t, g in spots[spots.photons > 0].groupby("t"):
                self.spots_by_t[int(t)] = (g.y_um.to_numpy() / self.px_um, g.x_um.to_numpy() / self.px_um, g.photons.to_numpy())

    def _bright(self, cid) -> float:
        if cid not in self.brightness:
            z = np.random.default_rng(np.random.SeedSequence([self.entropy, 2_000_001, int(cid)])).normal()
            self.brightness[cid] = float(np.exp(self.opt.get("brightness_cv", 0.3) * z))
        return self.brightness[cid]

    def _texture(self, r):
        if r.cell_id not in self.textures:
            rng = np.random.default_rng(np.random.SeedSequence([self.entropy, 2_000_002, int(r.cell_id)]))
            self.textures[r.cell_id] = make_texture(self.tex_cfg, r.radius_um, rng)
        return self.textures[r.cell_id]

    def clean(self, t: int, shift) -> tuple[np.ndarray, np.ndarray]:
        """Relative nuclear intensity (brightness x texture) and the label mask of frame t, objects moved by `shift` = (dy, dx) px."""
        ny, nx = self.shape
        clean = np.zeros((ny, nx), np.float32)
        lab = np.zeros((ny, nx), np.uint16)
        grp = self.by_t.get(int(t))
        if grp is None:
            return clean, lab
        for r in grp.itertuples():
            cy, cx = r.y_um / self.px_um + shift[0], r.x_um / self.px_um + shift[1]
            half = int(np.ceil(r.bound_radius_um / self.px_um)) + 2
            y0, y1 = max(int(cy) - half, 0), min(int(cy) + half + 1, ny)
            x0, x1 = max(int(cx) - half, 0), min(int(cx) + half + 1, nx)
            if y1 <= y0 or x1 <= x0:
                continue
            inside, xn, yn = _nucleus_mask(self.yy[y0:y1, x0:x1] - cy, self.xx[y0:y1, x0:x1] - cx, r, self.px_um, self.harm)
            lab[y0:y1, x0:x1][inside] = r.cell_id
            val = self._bright(r.cell_id) * np.ones(inside.shape, np.float32)
            if self.tex_cfg:
                val = val * _sample_texture(self._texture(r), xn, yn)
            clean[y0:y1, x0:x1][inside] = val[inside]
        return clean, lab


class _LocusModel:
    """One locus colour channel: diffraction-limited point sources at the true locus positions."""

    def __init__(self, truth: pd.DataFrame, cfg: dict, locus_id: int, occupancy: pd.DataFrame | None):
        lc, acq, opt = cfg["loci"], cfg["acquisition"], cfg["optics"]
        self.px_um = opt["pixel_size_nm"] / 1000.0
        self.wavelength = lc["wavelengths_nm"][locus_id]
        pps = lc["photons_per_locus_s"] * acq["exposure_s"]
        attached, per_probe = None, 0.0
        if occupancy is not None:
            o = occupancy[occupancy.locus_id == locus_id]
            attached = dict(zip(zip(o.t, o.cell_id), o.mean_attached))
            per_probe = lc["probes"]["photons_per_probe_s"] * acq["exposure_s"]
        self.by_t = {}
        for t, g in truth[truth.locus_id == locus_id].groupby("t"):
            ph = np.array([pps if attached is None else per_probe * attached[(t, c)] for c in g.cell_id])
            self.by_t[int(t)] = (g.y_um.to_numpy() / self.px_um, g.x_um.to_numpy() / self.px_um, ph)


class Imager:
    """Renders frames of every channel from the ground truth. Deterministic: frame t of a given (config, truth, entropy) is always the
    same, whatever other frames or error sources are rendered."""

    def __init__(self, cfg: dict, cells: pd.DataFrame | None = None, truth: pd.DataFrame | None = None,
                 spots: pd.DataFrame | None = None, occupancy: pd.DataFrame | None = None, entropy: int = 0,
                 errors: ErrorSet | None = None):
        self.cfg, self.entropy = cfg, int(entropy)
        self.opt, self.acq = cfg["optics"], cfg["acquisition"]
        self.E = errors if errors is not None else ErrorSet(cfg)
        self.px_um = self.opt["pixel_size_nm"] / 1000.0
        self.shape = tuple(int(round(v / self.px_um)) for v in cfg["geometry"]["fov_um"])
        self.T = int(self.acq["n_frames"])
        self.exposure_s = float(self.acq["exposure_s"])
        self.nucleus = _NucleusModel(cells, cfg, self.entropy, spots) if cells is not None else None
        n_loci = int((cfg.get("loci") or {}).get("n_loci", 0)) if truth is not None else 0
        self.loci = [_LocusModel(truth, cfg, k, occupancy) for k in range(n_loci)]
        # channel index is fixed (nucleus 0, locus k -> k + 1) so a channel's noise does not depend on which others are rendered
        self.channel_index = {"nucleus": 0, **{f"locus{k}": k + 1 for k in range(n_loci)}}
        self.channels = ([] if self.nucleus is None else ["nucleus"]) + [f"locus{k}" for k in range(n_loci)]
        # per-frame states that every channel shares or that do not depend on which sources are switched on
        # unit-size random series, drawn whether or not the source is on and scaled by the parameters at use time, so switching a source
        # on/off or changing its size never changes any other random draw
        r = ERR.rng_fixed(self.entropy, "stage_drift")
        self._drift_walk = np.vstack([np.zeros(2), np.cumsum(r.normal(0, 1, (self.T - 1, 2)), axis=0)])
        self._drift_jit = r.normal(0, 1, (self.T, 2))
        self._focus_white = ERR.rng_fixed(self.entropy, "defocus").normal(size=self.T)
        self._flick_white = {c: ERR.rng_fixed(self.entropy, "flicker", i).normal(size=self.T) for c, i in self.channel_index.items()}

    # ------------------------------------------------------------------ per-frame ground truth of the imaging conditions
    @staticmethod
    def _ar_series(white: np.ndarray, memory: float, sd: float) -> np.ndarray:
        """Stationary AR(1) series with standard deviation sd built from unit white noise."""
        memory = float(np.clip(memory, 0.0, 0.999))
        x, s = np.zeros(len(white)), np.sqrt(1 - memory**2)
        x[0] = white[0]
        for i in range(1, len(white)):
            x[i] = memory * x[i - 1] + s * white[i]
        return sd * x

    def drift_px(self, t: int, E: ErrorSet | None = None) -> tuple[float, float]:
        E = E or self.E
        if not E.on("stage_drift"):
            return 0.0, 0.0
        d = self._drift_walk[int(t)] * E.par("stage_drift", "drift_px_per_frame") + self._drift_jit[int(t)] * E.par("stage_drift", "jitter_px")
        return float(d[0]), float(d[1])

    def defocus_um(self, t: int, E: ErrorSet | None = None) -> float:
        E = E or self.E
        if not E.on("defocus"):
            return 0.0
        key = ("focus", E.par("defocus", "rms_um"), E.par("defocus", "memory"))
        if key not in E.cache:
            E.cache[key] = self._ar_series(self._focus_white, key[2], key[1])
        return float(E.cache[key][int(t)])

    def channel_shift_px(self, channel: str, E: ErrorSet | None = None) -> tuple[float, float]:
        E = E or self.E
        if channel == "nucleus" or not E.on("chromatic_shift"):
            return 0.0, 0.0
        s = E.par("chromatic_shift", "shift_px").get(channel, [0.0, 0.0])
        return float(s[0]), float(s[1])

    def _defocus_sigma_px(self, t: int, E: ErrorSet) -> float:
        z = self.defocus_um(t, E)
        if z == 0.0:
            return 0.0
        n = float(E.par("defocus", "immersion_n"))
        na = min(float(self.opt["NA"]), 0.95 * n)
        rho = abs(z) * na / np.sqrt(n * n - na * na)                      # geometric blur-disc radius (um)
        return 0.5 * rho / self.px_um

    def sigma_px(self, channel: str, t: int, E: ErrorSet | None = None) -> float:
        """Total Gaussian blur of the channel in frame t (pixels): diffraction + defocus (+ locus motion), added in quadrature."""
        E = E or self.E
        wl = self.opt["wavelength_nm"] if channel == "nucleus" else self.cfg["loci"]["wavelengths_nm"][int(channel[5:])]
        s2 = psf_sigma_px({**self.opt, "wavelength_nm": wl}) ** 2 if E.on("psf_blur") else 0.0
        s2 += self._defocus_sigma_px(t, E) ** 2
        if channel != "nucleus" and E.on("exposure_motion"):
            s2 += (E.par("exposure_motion", "locus_sigma_nm") / 1000.0 / self.px_um) ** 2
        return float(np.sqrt(s2))

    def flicker(self, channel: str, t: int, E: ErrorSet | None = None) -> float:
        E = E or self.E
        if not E.on("flicker"):
            return 1.0
        cv, mem = E.par("flicker", "cv"), E.par("flicker", "memory")
        key = ("flick", channel, cv, mem)
        if key not in E.cache:
            E.cache[key] = self._ar_series(self._flick_white[channel], mem, cv)
        return float(1.0 + E.cache[key][int(t)])

    def bleach(self, channel: str, t: int, E: ErrorSet | None = None) -> float:
        E = E or self.E
        tau = self.opt.get("bleach_tau_s") if channel == "nucleus" else self.cfg["loci"].get("bleach_tau_s", self.opt.get("bleach_tau_s"))
        if not E.on("bleaching") or not tau:
            return 1.0
        return float(np.exp(-self.acq["frame_interval_s"] * t / tau))

    # ------------------------------------------------------------------ fixed fields
    def _flat(self, channel: str, E: ErrorSet) -> np.ndarray | None:
        if not E.on("flat_field"):
            return None
        key = ("flat", channel, E.par("flat_field", "vignette"), E.par("flat_field", "tilt"))
        if key not in E.cache:
            E.cache[key] = ERR.flat_field_map(self.shape, key[2], key[3], ERR.rng_fixed(self.entropy, "flat_field", self.channel_index[channel]))
        return E.cache[key]

    def _autofluor(self, channel: str, E: ErrorSet) -> np.ndarray | float:
        if not E.on("autofluorescence"):
            return 0.0
        p = E.params["autofluorescence"]
        key = ("af", channel, p["photons"], p["cv"], p["correlation_um"])
        if key not in E.cache:
            f = ERR.smooth_field(self.shape, p["correlation_um"] / self.px_um, ERR.rng_fixed(self.entropy, "autofluorescence", self.channel_index[channel]))
            E.cache[key] = np.clip(p["photons"] * (1 + p["cv"] * f), 0, None).astype(np.float32)
        return E.cache[key]

    def _camera(self, channel: str, E: ErrorSet) -> ERR.Camera:
        key = ("cam", channel)
        if key not in E.cache:
            E.cache[key] = ERR.Camera(self.shape, E, self.opt, self.exposure_s, self.entropy, self.channel_index[channel])
        return E.cache[key]

    # ------------------------------------------------------------------ one frame
    def signal(self, channel: str, t: int, E: ErrorSet | None = None):
        """Noise-free photons per pixel of one channel in frame t, after optics, bleaching, illumination (no background, no crosstalk).
        Returns (photons, labels or None)."""
        E = E or self.E
        dy, dx = self.drift_px(t, E)
        cy, cx = self.channel_shift_px(channel, E)
        sig = self.sigma_px(channel, t, E)
        bleach = self.bleach(channel, t, E)
        lab = None
        if channel == "nucleus":
            clean, lab = self.nucleus.clean(t, (dy, dx))
            ph = (gaussian_filter(clean, sig) if sig > 0 else clean) * (self.opt["photons_per_px_s"] * self.exposure_s) * bleach
            s = self.nucleus.spots_by_t.get(int(t))
            if s is not None:
                ph = ph + bleach * _spots_from_arrays(self.shape, s[0] + dy, s[1] + dx, s[2], sig)
        else:
            y, x, p = self.loci[int(channel[5:])].by_t.get(int(t), (np.zeros(0), np.zeros(0), np.zeros(0)))
            ph = _spots_from_arrays(self.shape, y + dy + cy, x + dx + cx, p * bleach, sig)
        if E.on("haze"):
            f = float(E.par("haze", "fraction"))
            ph = (1 - f) * ph + f * ERR.broad_blur(ph, E.par("haze", "sigma_um") / self.px_um)
        flat = self._flat(channel, E)
        if flat is not None:
            ph = ph * flat
        ph = ph * self.flicker(channel, t, E)
        return ph.astype(np.float32), lab

    def frame(self, t: int, E: ErrorSet | None = None, channels: list[str] | None = None) -> dict:
        """Render frame t. Returns {"images": {channel: ADU float32}, "labels": uint16 or None, "signals": {channel: photons}}."""
        E = E or self.E
        chans = [c for c in (channels or self.channels) if c in self.channels]
        need = set(chans)
        pairs = [(a, b, f) for a, b, f in E.par("crosstalk", "pairs")] if E.on("crosstalk") else []
        need |= {a for a, b, f in pairs if b in chans and a in self.channels}
        sig, labels = {}, None
        for c in need:
            sig[c], lab = self.signal(c, t, E)
            if c == "nucleus":
                labels = lab
        images = {}
        for c in chans:
            m = sig[c]
            for a, b, f in pairs:
                if b == c and a in sig:
                    m = m + f * sig[a]
            bg = (float(self.opt["background_photons"]) if E.on("background") else 0.0)
            af = self._autofluor(c, E)
            if not np.isscalar(af):
                fl = self._flat(c, E)
                af = af * (fl if fl is not None else 1.0) * self.flicker(c, t, E)
            images[c] = self._camera(c, E)(m + bg + af, t)
        return {"images": images, "labels": labels, "signals": sig}

    # ------------------------------------------------------------------ whole movie and its ground truth
    def render_movie(self, channels: list[str] | None = None, progress=None, dtype=np.float32) -> dict:
        """All frames. Returns {"images": {channel: (T, ny, nx) dtype}, "labels": (T, ny, nx) uint16 or None}. dtype=np.uint16 stores
        the camera counts as 16-bit integers (half the memory)."""
        chans = [c for c in (channels or self.channels) if c in self.channels]
        out = {c: np.zeros((self.T, *self.shape), dtype) for c in chans}
        lab = np.zeros((self.T, *self.shape), np.uint16) if "nucleus" in chans else None
        integer = np.issubdtype(np.dtype(dtype), np.integer)
        for t in range(self.T):
            f = self.frame(t, channels=chans)
            for c in chans:
                out[c][t] = np.rint(f["images"][c]).clip(0, 65535) if integer else f["images"][c]
            if lab is not None and f["labels"] is not None:
                lab[t] = f["labels"]
            if progress and t % 20 == 0:
                progress(f"  frame {t + 1}/{self.T}")
        return {"images": out, "labels": lab}

    def imaging_truth(self) -> pd.DataFrame:
        """What the imaging conditions really were, per frame: stage drift (px), focus error (um), power factor of every channel."""
        rows = []
        for t in range(self.T):
            dy, dx = self.drift_px(t)
            row = {"t": t, "drift_y_px": dy, "drift_x_px": dx, "defocus_um": self.defocus_um(t),
                   "defocus_sigma_px": self._defocus_sigma_px(t, self.E)}
            for c in self.channels:
                row[f"flicker_{c}"] = self.flicker(c, t)
                row[f"bleach_{c}"] = self.bleach(c, t)
                row[f"sigma_px_{c}"] = self.sigma_px(c, t)
            rows.append(row)
        return pd.DataFrame(rows)

    def static_truth(self) -> dict:
        """Constant imaging properties: per-channel shift and blur, pixel size, which sources were on and with what parameters."""
        return {"pixel_size_um": self.px_um, "exposure_s": self.exposure_s,
                "channel_shift_px": {c: self.channel_shift_px(c) for c in self.channels},
                "sources": self.E.summary(), "master_switch": self.E.master}


# ---------------------------------------------------------------------- ground truth in image coordinates
def truth_in_image(truth: pd.DataFrame, imaging_truth: pd.DataFrame, px_um: float, channel_shift_px: dict | None = None) -> pd.DataFrame:
    """Locus truth moved into image coordinates: + stage drift (y_um, x_um; the drift-corrected truth to score a localizer against)
    and, if `channel_shift_px` is given, extra columns y_drawn_um / x_drawn_um = where the spot is really drawn (drift + chromatic shift)."""
    d = imaging_truth.set_index("t")[["drift_y_px", "drift_x_px"]]
    out = truth.copy()
    dy, dx = d.loc[out.t, "drift_y_px"].to_numpy() * px_um, d.loc[out.t, "drift_x_px"].to_numpy() * px_um
    out["y_um"] = out.y_um + dy
    out["x_um"] = out.x_um + dx
    if channel_shift_px is not None:
        sh = np.array([channel_shift_px.get(f"locus{int(k)}", (0.0, 0.0)) for k in out.locus_id], float) * px_um
        out["y_drawn_um"], out["x_drawn_um"] = out.y_um + sh[:, 0], out.x_um + sh[:, 1]
    return out


def cells_in_image(cells: pd.DataFrame, imaging_truth: pd.DataFrame, px_um: float) -> pd.DataFrame:
    """cells.csv with positions moved by the stage drift, i.e. where the nuclei really are in the image."""
    d = imaging_truth.set_index("t")[["drift_y_px", "drift_x_px"]]
    out = cells.copy()
    out["y_um"] = out.y_um + d.loc[out.t, "drift_y_px"].to_numpy() * px_um
    out["x_um"] = out.x_um + d.loc[out.t, "drift_x_px"].to_numpy() * px_um
    return out


def imager_for_run(run, cfg: dict, errors: ErrorSet | None = None) -> Imager:
    """Rebuild the exact Imager of a finished run from its saved ground truth files (reproduces its images frame by frame)."""
    from pathlib import Path

    run = Path(run)
    cells = pd.read_csv(run / "stage2_cells" / "cells.csv")
    s1 = run / "stage1_chromatin"
    truth = pd.read_csv(s1 / "loci_truth.csv") if (s1 / "loci_truth.csv").exists() else None
    occ = pd.read_csv(s1 / "probe_occupancy.csv") if (s1 / "probe_occupancy.csv").exists() else None
    spots = None
    if (s1 / "transcription.csv").exists() and truth is not None:
        from simlive.stage1_chromatin.transcription import ms2_spots

        spots = ms2_spots(pd.read_csv(s1 / "transcription.csv"), truth, cfg)
    return Imager(cfg, cells, truth, spots, occ, entropy_for(cfg["seed"]), errors)


# ---------------------------------------------------------------------- one-channel wrappers (kept for tests and simple use)
def render_nuclei(cells: pd.DataFrame, cfg: dict, rng: np.random.Generator, spots: pd.DataFrame | None = None):
    """Nuclear channel (NLS-GFP) and its label masks. `spots` (columns t, y_um, x_um, photons) adds the MS2 spots in the same channel."""
    im = Imager(cfg, cells, spots=spots, entropy=int(rng.integers(2**63)))
    m = im.render_movie(["nucleus"])
    return m["images"]["nucleus"], m["labels"]


def render_loci(loci: pd.DataFrame, cfg: dict, rng: np.random.Generator, locus_id: int = 0, occupancy: pd.DataFrame | None = None):
    """One locus colour channel (diffraction-limited spots at the true positions); `occupancy` = stage 1b probe model."""
    im = Imager(cfg, truth=loci, occupancy=occupancy, entropy=int(rng.integers(2**63)))
    return im.render_movie([f"locus{locus_id}"])["images"][f"locus{locus_id}"]
