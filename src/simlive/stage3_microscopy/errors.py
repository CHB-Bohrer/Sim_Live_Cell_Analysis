"""Stage 3: every source of image error, in one place.

A real microscope image differs from the ground truth (where the nuclei and loci really are) in many independent ways. Each one is
a *source* with its own on/off switch and parameters, listed in `SOURCES` (the dashboard reads this table to explain them) and set
in the config block `imaging_errors:`. The physical order in which light is degraded is `ORDER`:

    sample position (stage drift, chromatic shift) -> optics (blur, defocus, haze, exposure motion) -> illumination (bleaching,
    flat-field, flicker) -> background (constant, autofluorescence, channel cross-talk) -> detection (QE, shot noise, dark current)
    -> camera electronics (hot pixels, cosmic rays, PRNU, read noise, row noise, offset pattern, ADC / saturation).

`baseline` sources are the original model (blur, bleaching, constant background, QE, shot noise, read noise) and stay on when
`imaging_errors.enabled` is false, so old configs still mean what they meant. With `enabled: true` every source follows its own flag.

Randomness: every source draws from its own stream, seeded by (entropy, frame, source, channel), so switching one source on or off
never changes the noise of another and any single frame can be re-rendered on its own (the dashboard's Images tab relies on this).
"""
from __future__ import annotations

import copy

import numpy as np
from scipy.ndimage import gaussian_filter

GROUPS = ["Sample position", "Optics", "Illumination", "Background", "Detection", "Camera electronics"]

# id -> spec. params: name -> (default, unit, plain-language meaning). `baseline`: part of the original model.
SOURCES: dict[str, dict] = {
    "stage_drift": dict(
        group="Sample position", title="Stage drift and vibration", baseline=False,
        what="The sample slowly slides relative to the camera (thermal drift, a stage that settles) and shakes a little between frames.",
        effect="Every object in a frame moves together by a fraction of a pixel. The true cell masks move with it, so tracking is still "
               "scored fairly; locus positions are scored against the drift-corrected truth.",
        params={"drift_px_per_frame": (0.15, "px", "step size of the slow random-walk drift, per frame"),
                "jitter_px": (0.05, "px", "frame-to-frame vibration (does not accumulate)")}),
    "chromatic_shift": dict(
        group="Sample position", title="Chromatic shift between colour channels", baseline=False,
        what="Different wavelengths focus at slightly different places (lateral chromatic aberration, imperfect filter-cube alignment), "
             "so each locus colour is imaged shifted relative to the nucleus channel.",
        effect="A constant offset per colour. It does not average away, and it biases the measured distance between two loci of "
               "different colours unless it is calibrated (e.g. with beads).",
        params={"shift_px": ({"locus0": [0.25, -0.15], "locus1": [-0.30, 0.20]}, "px (y, x)",
                             "offset of each locus channel relative to the nuclear channel")}),
    "psf_blur": dict(
        group="Optics", title="Diffraction blur (PSF)", baseline=True,
        what="A point of light is imaged as a small blob, not a point: sigma = 0.21 x wavelength / NA.",
        effect="Spots and nuclear edges are blurred; a spot's peak is lower the wider the blur. Sets the best possible localization.",
        params={}),
    "defocus": dict(
        group="Optics", title="Focus drift (defocus)", baseline=False,
        what="The focal plane wanders a little from frame to frame (autofocus residual, thermal drift). Out of focus the PSF widens.",
        effect="Spots get wider and dimmer in the same frames; the nuclear edge softens. Spot brightness fluctuates for a reason "
               "unrelated to biology.",
        params={"rms_um": (0.10, "um", "typical focus error (root mean square)"),
                "memory": (0.90, "", "how slowly the focus error changes from frame to frame (0 = new each frame, 1 = frozen)"),
                "immersion_n": (1.33, "", "refractive index of the immersion medium (sets how fast the blur grows with defocus)")}),
    "haze": dict(
        group="Optics", title="Out-of-focus haze (widefield)", baseline=False,
        what="A widefield microscope also collects light from above and below the focal plane; it arrives as a broad, faint glow.",
        effect="A fraction of every object's light is spread over several micrometres: lower contrast, brighter surroundings of "
               "bright nuclei and spots.",
        params={"fraction": (0.06, "", "fraction of the light that ends up in the broad glow"),
                "sigma_um": (2.5, "um", "width of the glow")}),
    "exposure_motion": dict(
        group="Optics", title="Locus motion during the exposure", baseline=False,
        what="A locus moves while the shutter is open, so its spot is smeared, not a perfect point.",
        effect="Spots are slightly wider than the PSF. Negligible at short exposures; large once the exposure is long compared with "
               "how fast chromatin moves (set it from your measured MSD).",
        params={"locus_sigma_nm": (25.0, "nm", "extra Gaussian smear of a locus spot")}),
    "bleaching": dict(
        group="Illumination", title="Photobleaching", baseline=True,
        what="Fluorophores are destroyed by the excitation light, so every channel dims exponentially during the movie.",
        effect="Signal falls over time (nuclear: optics.bleach_tau_s; loci: loci.bleach_tau_s); late frames are noisier.", params={}),
    "flat_field": dict(
        group="Illumination", title="Uneven illumination (flat-field)", baseline=False,
        what="The excitation light is brighter in the middle of the field than at the edges (vignetting) and may be tilted.",
        effect="The same nucleus looks dimmer near the edge of the field; thresholds and detection scores depend on position.",
        params={"vignette": (0.15, "", "fractional loss of brightness at the corners"),
                "tilt": (0.03, "", "fractional left-right/top-bottom brightness gradient")}),
    "flicker": dict(
        group="Illumination", title="Lamp / laser flicker", baseline=False,
        what="The excitation power fluctuates a little from frame to frame.",
        effect="Every object in the channel gets brighter or dimmer together in a frame.",
        params={"cv": (0.01, "", "relative size of the fluctuations (0.01 = 1%)"),
                "memory": (0.8, "", "how slowly the power changes between frames")}),
    "background": dict(
        group="Background", title="Constant background light", baseline=True,
        what="Stray light and camera-side background that adds the same number of photons to every pixel.",
        effect="Raises the floor under everything and adds shot noise of its own (optics.background_photons).", params={}),
    "autofluorescence": dict(
        group="Background", title="Autofluorescence", baseline=False,
        what="Cells, medium and plastic fluoresce a little on their own, unevenly across the field.",
        effect="A smooth, uneven background that is not the same at every position; hides faint spots.",
        params={"photons": (6.0, "photons/px", "mean extra background per pixel"),
                "cv": (0.3, "", "how much it varies across the field (relative)"),
                "correlation_um": (10.0, "um", "size of the patches")}),
    "crosstalk": dict(
        group="Background", title="Colour cross-talk (bleed-through)", baseline=False,
        what="A fluorophore's emission is not confined to its own channel: some of its light leaks into a neighbouring colour.",
        effect="Faint ghost copies of bright objects (the whole nucleus, or a locus) in other channels.",
        params={"pairs": ([["nucleus", "locus0", 0.02], ["locus0", "locus1", 0.02]], "[from, to, fraction]",
                          "fraction of the first channel's light that also lands in the second")}),
    "qe": dict(
        group="Detection", title="Quantum efficiency", baseline=True,
        what="Only a fraction of the photons that reach the sensor release an electron (optics.quantum_efficiency).",
        effect="Fewer detected events means more relative shot noise.", params={}),
    "shot_noise": dict(
        group="Detection", title="Photon shot noise", baseline=True,
        what="Light arrives in discrete photons, so the count in a pixel fluctuates (Poisson) even for a perfectly steady source.",
        effect="Grainy images: noise grows as the square root of the signal; the dominant error for dim spots.",
        params={"excess_noise_factor": (1.0, "", "extra multiplicative noise of EM-type cameras (1 = none, 1.41 = EMCCD)")}),
    "dark_current": dict(
        group="Detection", title="Dark current", baseline=False,
        what="Heat creates electrons in the sensor even in darkness.",
        effect="A small added signal with its own shot noise; matters for long exposures.",
        params={"e_per_s": (0.2, "e-/px/s", "dark electrons per pixel per second")}),
    "hot_pixels": dict(
        group="Camera electronics", title="Hot pixels", baseline=False,
        what="A few pixels always read too high (fixed positions).",
        effect="Isolated bright dots that look like tiny spots to a naive detector.",
        params={"fraction": (3e-4, "", "fraction of pixels that are hot"),
                "e_per_frame": (40.0, "e-", "extra electrons in a hot pixel")}),
    "cosmic_rays": dict(
        group="Camera electronics", title="Cosmic rays", baseline=False,
        what="High-energy particles occasionally hit the sensor and deposit charge in a few pixels.",
        effect="Rare, bright, tiny events in single frames; they can be mistaken for a transient locus signal.",
        params={"rate_per_frame": (0.05, "events", "average events per frame"),
                "e": (3000.0, "e-", "typical charge of one event")}),
    "prnu": dict(
        group="Camera electronics", title="Pixel gain differences (PRNU)", baseline=False,
        what="Every pixel converts light to charge with a slightly different gain.",
        effect="A fixed, speckle-like multiplicative pattern, proportional to the signal.",
        params={"cv": (0.01, "", "relative spread of pixel gains (0.01 = 1%)")}),
    "read_noise": dict(
        group="Camera electronics", title="Read noise", baseline=True,
        what="Reading the pixel adds electronic noise independent of the signal (optics.read_noise_e).",
        effect="Sets the noise floor in dark regions; dominates for very dim signals.",
        params={"pixel_cv": (0.25, "", "sCMOS: relative spread of the read noise between pixels (0 = identical pixels)")}),
    "row_noise": dict(
        group="Camera electronics", title="Row (banding) noise", baseline=False,
        what="CMOS readout adds a random offset to each row of pixels in each frame.",
        effect="Faint horizontal stripes that change from frame to frame.",
        params={"adu": (0.6, "ADU", "size of the per-row offset")}),
    "dsnu": dict(
        group="Camera electronics", title="Offset pattern (DSNU)", baseline=False,
        what="Each pixel has its own tiny fixed electronic offset.",
        effect="A fixed additive pattern, visible only in dim regions.",
        params={"adu": (1.0, "ADU", "spread of the per-pixel offset")}),
    "adc": dict(
        group="Camera electronics", title="Digitisation and saturation", baseline=False,
        what="The camera outputs whole numbers (ADU) with a limited range; very bright pixels clip at the full-well capacity.",
        effect="Rounding to integers (coarse in dim images) and a flat ceiling on bright objects.",
        params={"bit_depth": (16, "bits", "number of bits of the converter"),
                "full_well_e": (30000.0, "e-", "charge at which a pixel saturates")}),
}

# physical order in which a frame is degraded (used for the cumulative build-up in the dashboard)
ORDER = ["stage_drift", "chromatic_shift", "psf_blur", "defocus", "exposure_motion", "haze", "bleaching", "flat_field", "flicker",
         "background", "autofluorescence", "crosstalk", "qe", "shot_noise", "dark_current", "hot_pixels", "cosmic_rays", "prnu",
         "read_noise", "row_noise", "dsnu", "adc"]
assert set(ORDER) == set(SOURCES)
INDEX = {s: i + 1 for i, s in enumerate(ORDER)}
BASELINE = [s for s in ORDER if SOURCES[s]["baseline"]]
# values of baseline parameters that reproduce the original model when imaging_errors is off
LEGACY = {("read_noise", "pixel_cv"): 0.0, ("shot_noise", "excess_noise_factor"): 1.0}


class ErrorSet:
    """The resolved on/off state and parameters of every error source for one config (or a variant of it)."""

    def __init__(self, cfg: dict | None = None, enabled: set | None = None):
        blk = (cfg or {}).get("imaging_errors") or {}
        self.master = bool(blk.get("enabled", False))
        self.params: dict[str, dict] = {}
        self._on: dict[str, bool] = {}
        for sid, spec in SOURCES.items():
            user = blk.get(sid) or {} if self.master else {}
            self._on[sid] = bool(user.get("enabled", spec["baseline"])) if self.master else spec["baseline"]
            p = {k: copy.deepcopy(v[0]) for k, v in spec["params"].items()}
            for k in p:
                if (sid, k) in LEGACY and not self.master:
                    p[k] = LEGACY[(sid, k)]
                elif k in user:
                    p[k] = copy.deepcopy(user[k])
            self.params[sid] = p
        if enabled is not None:
            self._on = {sid: sid in enabled for sid in SOURCES}
        self.cache: dict = {}          # fixed patterns / cameras built for this exact set of switches (filled by the Imager)

    def on(self, sid: str) -> bool:
        return self._on[sid]

    def par(self, sid: str, name: str):
        return self.params[sid][name]

    def enabled_ids(self) -> list[str]:
        return [s for s in ORDER if self._on[s]]

    def with_enabled(self, ids) -> "ErrorSet":
        """A copy in which exactly `ids` are on (everything else off), with the same parameters."""
        out = ErrorSet.__new__(ErrorSet)
        out.master, out.params = self.master, copy.deepcopy(self.params)
        out._on = {sid: sid in set(ids) for sid in SOURCES}
        out.cache = {}
        return out

    def summary(self) -> dict:
        return {s: {"on": self._on[s], **self.params[s]} for s in ORDER}


def rng_for(entropy: int, t: int, sid: str, channel: int = 0) -> np.random.Generator:
    """Independent random stream for one source in one frame of one channel."""
    return np.random.default_rng(np.random.SeedSequence([int(entropy), int(t), INDEX[sid], int(channel)]))


def rng_fixed(entropy: int, sid: str, channel: int = 0, extra: int = 0) -> np.random.Generator:
    """Stream for a fixed (time-independent) pattern of a source, e.g. the pixel gain map of one channel."""
    return np.random.default_rng(np.random.SeedSequence([int(entropy), 1_000_000 + INDEX[sid], int(channel), int(extra)]))


def ar1(rng: np.random.Generator, n: int, memory: float, sd: float) -> np.ndarray:
    """Stationary AR(1) series with standard deviation sd and one-step memory in [0, 1)."""
    memory = float(np.clip(memory, 0.0, 0.999))
    x = np.zeros(n)
    x[0] = rng.normal(0, sd)
    s = sd * np.sqrt(1 - memory**2)
    for i in range(1, n):
        x[i] = memory * x[i - 1] + rng.normal(0, s)
    return x


def smooth_field(shape, corr_px: float, rng: np.random.Generator) -> np.ndarray:
    """Zero-mean, unit-variance random field with spatial correlation length corr_px (computed on a coarse grid, then upsampled)."""
    from scipy.ndimage import zoom

    k = max(1, int(corr_px // 4))
    small = (int(np.ceil(shape[0] / k)) + 2, int(np.ceil(shape[1] / k)) + 2)
    f = gaussian_filter(rng.normal(size=small), max(corr_px / k, 0.5), mode="reflect")
    f = zoom(f, k, order=1)[: shape[0], : shape[1]]
    return ((f - f.mean()) / (f.std() + 1e-12)).astype(np.float32)


def flat_field_map(shape, vignette: float, tilt: float, rng: np.random.Generator) -> np.ndarray:
    """Multiplicative illumination map (mean ~1): quadratic fall-off towards the corners plus a random-direction linear tilt."""
    ny, nx = shape
    yy, xx = np.mgrid[:ny, :nx]
    v, u = (yy - (ny - 1) / 2) / (ny / 2), (xx - (nx - 1) / 2) / (nx / 2)
    ang = rng.uniform(0, 2 * np.pi)
    m = 1 - vignette * (u**2 + v**2) / 2 + tilt * (np.cos(ang) * u + np.sin(ang) * v)
    return (m / m.mean()).astype(np.float32)


def broad_blur(img: np.ndarray, sigma_px: float) -> np.ndarray:
    """Gaussian blur with a large sigma, computed on a block-averaged copy (fast; the result is smooth anyway). Light is conserved
    except what leaves the field of view."""
    from scipy.ndimage import zoom

    k = int(max(1, sigma_px // 6))
    if k == 1:
        return gaussian_filter(img, sigma_px, mode="constant")
    ny, nx = img.shape
    py, px = (-ny) % k, (-nx) % k
    pad = np.pad(img, ((0, py), (0, px)))
    small = pad.reshape(pad.shape[0] // k, k, pad.shape[1] // k, k).sum((1, 3))             # photons per block
    small = gaussian_filter(small, sigma_px / k, mode="constant")
    big = zoom(small, k, order=1) / (k * k)
    return big[:ny, :nx].astype(np.float32)


class Camera:
    """Detection and camera electronics of one channel: photons per pixel in, ADU out. Fixed patterns (pixel gain, offsets, hot pixels,
    read-noise map) are created once from (entropy, channel); per-frame noise comes from per-source streams."""

    def __init__(self, shape, E: ErrorSet, opt: dict, exposure_s: float, entropy: int, channel: int):
        self.shape, self.E, self.opt, self.exposure_s, self.entropy, self.ch = tuple(shape), E, opt, exposure_s, entropy, channel
        ny, nx = shape
        self.qe = float(opt.get("quantum_efficiency", 0.8)) if E.on("qe") else 1.0
        self.gain = float(opt.get("gain_adu_per_e", 1.0))
        self.offset = float(opt.get("offset_adu", 100.0))
        self.read_e = float(opt["read_noise_e"]) if E.on("read_noise") else 0.0
        self.read_map = None
        if E.on("read_noise") and E.par("read_noise", "pixel_cv") > 0:
            self.read_map = np.clip(1 + E.par("read_noise", "pixel_cv") * rng_fixed(entropy, "read_noise", channel).normal(size=shape),
                                    0.3, None).astype(np.float32)
        self.prnu = (1 + E.par("prnu", "cv") * rng_fixed(entropy, "prnu", channel).normal(size=shape)).astype(np.float32) \
            if E.on("prnu") else None
        self.dsnu = (E.par("dsnu", "adu") * rng_fixed(entropy, "dsnu", channel).normal(size=shape)).astype(np.float32) \
            if E.on("dsnu") else None
        self.hot = None
        if E.on("hot_pixels"):
            r = rng_fixed(entropy, "hot_pixels", channel)
            self.hot = ((r.random(shape) < E.par("hot_pixels", "fraction")) * E.par("hot_pixels", "e_per_frame")
                        * r.uniform(0.5, 1.5, shape)).astype(np.float32)

    def __call__(self, photons: np.ndarray, t: int) -> np.ndarray:
        E, ny_nx = self.E, self.shape
        e_mean = np.clip(photons, 0, None).astype(np.float64) * self.qe
        e = rng_for(self.entropy, t, "shot_noise", self.ch).poisson(e_mean).astype(np.float64) if E.on("shot_noise") else e_mean
        f = float(E.par("shot_noise", "excess_noise_factor")) if E.on("shot_noise") else 1.0
        if f > 1.0:
            e += rng_for(self.entropy, t, "shot_noise", self.ch + 100).normal(size=ny_nx) * np.sqrt((f**2 - 1.0) * e)
        if E.on("dark_current"):
            lam = E.par("dark_current", "e_per_s") * self.exposure_s
            e += rng_for(self.entropy, t, "dark_current", self.ch).poisson(lam, ny_nx)
        if self.hot is not None:
            e += self.hot
        if E.on("cosmic_rays"):
            r = rng_for(self.entropy, t, "cosmic_rays", self.ch)
            for _ in range(int(r.poisson(E.par("cosmic_rays", "rate_per_frame")))):
                cy, cx = r.uniform(0, ny_nx[0]), r.uniform(0, ny_nx[1])
                y0, x0 = int(cy) - 3, int(cx) - 3
                ys, xs = np.arange(max(y0, 0), min(y0 + 7, ny_nx[0])), np.arange(max(x0, 0), min(x0 + 7, ny_nx[1]))
                if len(ys) and len(xs):
                    blob = np.exp(-((ys[:, None] - cy) ** 2 + (xs[None, :] - cx) ** 2) / (2 * 0.8**2))
                    e[ys[0]:ys[-1] + 1, xs[0]:xs[-1] + 1] += r.lognormal(np.log(E.par("cosmic_rays", "e")), 0.5) * blob
        if self.prnu is not None:
            e *= self.prnu
        if self.read_e > 0:
            sd = self.read_e * (self.read_map if self.read_map is not None else 1.0)
            e += rng_for(self.entropy, t, "read_noise", self.ch).normal(size=ny_nx) * sd
        adu = e * self.gain + self.offset
        if E.on("row_noise"):
            adu += rng_for(self.entropy, t, "row_noise", self.ch).normal(0, E.par("row_noise", "adu"), (ny_nx[0], 1))
        if self.dsnu is not None:
            adu += self.dsnu
        if E.on("adc"):
            top = min(2 ** int(E.par("adc", "bit_depth")) - 1, self.offset + E.par("adc", "full_well_e") * self.gain)
            adu = np.rint(adu).clip(0, top)
        else:
            adu = adu.clip(0, 65535)
        return adu.astype(np.float32)
