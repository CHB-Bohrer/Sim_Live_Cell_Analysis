# Sources of image error (stage 3)

Generated from `src/simlive/stage3_microscopy/errors.py` by `scripts/make_image_error_docs.py`; do not edit by hand.
All numbers below are **placeholders** until measured on the real microscope. Switch them in the config block `imaging_errors:`
(each source has `enabled: true|false` plus its parameters), or in the dashboard sidebar ("Image errors"). With
`imaging_errors.enabled: false` only the original model remains: blur, bleaching, constant background, quantum efficiency,
shot noise and read noise.

How a frame is made, in physical order: objects at their true positions -> stage drift and chromatic shift move them -> optical blur
(diffraction, defocus, exposure motion, out-of-focus haze) -> photobleaching -> uneven illumination and flicker -> colour cross-talk ->
background (constant, autofluorescence) -> detector (quantum efficiency, shot noise, dark current) -> camera electronics (hot pixels,
cosmic rays, pixel gain, read noise, row noise, offset pattern, digitisation and saturation).

Ground truth saved with every run (`stage3_microscopy/`): `labels.tif` (true masks, moved by the drift), `imaging_truth.csv` (per frame: drift,
focus error, blur per channel, light level per channel), `cells_image.csv` and `loci_truth_image.csv` (truth in image coordinates; the latter
also has where each colour's spot is really drawn), `imaging_errors.json` (which sources were on, with all parameters).
The dashboard's **Images** tab re-renders any frame from this truth, with any subset of sources, and checks it equals the saved image.

Reproducibility: every source draws from its own random stream seeded by (run, frame, source, channel), so switching one source never changes
another source's noise, and any frame can be rendered alone.


## Sample position

### Stage drift and vibration  (`stage_drift`)

The sample slowly slides relative to the camera (thermal drift, a stage that settles) and shakes a little between frames.

**Effect:** Every object in a frame moves together by a fraction of a pixel. The true cell masks move with it, so tracking is still scored fairly; locus positions are scored against the drift-corrected truth.

| parameter | default | unit | meaning |
|---|---|---|---|
| `drift_px_per_frame` | `0.15` | px | step size of the slow random-walk drift, per frame |
| `jitter_px` | `0.05` | px | frame-to-frame vibration (does not accumulate) |

### Chromatic shift between colour channels  (`chromatic_shift`)

Different wavelengths focus at slightly different places (lateral chromatic aberration, imperfect filter-cube alignment), so each locus colour is imaged shifted relative to the nucleus channel.

**Effect:** A constant offset per colour. It does not average away, and it biases the measured distance between two loci of different colours unless it is calibrated (e.g. with beads).

| parameter | default | unit | meaning |
|---|---|---|---|
| `shift_px` | `{'locus0': [0.25, -0.15], 'locus1': [-0.3, 0.2]}` | px (y, x) | offset of each locus channel relative to the nuclear channel |


## Optics

### Diffraction blur (PSF)  (`psf_blur`, original model)

A point of light is imaged as a small blob, not a point: sigma = 0.21 x wavelength / NA.

**Effect:** Spots and nuclear edges are blurred; a spot's peak is lower the wider the blur. Sets the best possible localization.

### Focus drift (defocus)  (`defocus`)

The focal plane wanders a little from frame to frame (autofocus residual, thermal drift). Out of focus the PSF widens.

**Effect:** Spots get wider and dimmer in the same frames; the nuclear edge softens. Spot brightness fluctuates for a reason unrelated to biology.

| parameter | default | unit | meaning |
|---|---|---|---|
| `rms_um` | `0.1` | um | typical focus error (root mean square) |
| `memory` | `0.9` |  | how slowly the focus error changes from frame to frame (0 = new each frame, 1 = frozen) |
| `immersion_n` | `1.33` |  | refractive index of the immersion medium (sets how fast the blur grows with defocus) |

### Locus motion during the exposure  (`exposure_motion`)

A locus moves while the shutter is open, so its spot is smeared, not a perfect point.

**Effect:** Spots are slightly wider than the PSF. Negligible at short exposures; large once the exposure is long compared with how fast chromatin moves (set it from your measured MSD).

| parameter | default | unit | meaning |
|---|---|---|---|
| `locus_sigma_nm` | `25.0` | nm | extra Gaussian smear of a locus spot |

### Out-of-focus haze (widefield)  (`haze`)

A widefield microscope also collects light from above and below the focal plane; it arrives as a broad, faint glow.

**Effect:** A fraction of every object's light is spread over several micrometres: lower contrast, brighter surroundings of bright nuclei and spots.

| parameter | default | unit | meaning |
|---|---|---|---|
| `fraction` | `0.06` |  | fraction of the light that ends up in the broad glow |
| `sigma_um` | `2.5` | um | width of the glow |


## Illumination

### Photobleaching  (`bleaching`, original model)

Fluorophores are destroyed by the excitation light, so every channel dims exponentially during the movie.

**Effect:** Signal falls over time (nuclear: optics.bleach_tau_s; loci: loci.bleach_tau_s); late frames are noisier.

### Uneven illumination (flat-field)  (`flat_field`)

The excitation light is brighter in the middle of the field than at the edges (vignetting) and may be tilted.

**Effect:** The same nucleus looks dimmer near the edge of the field; thresholds and detection scores depend on position.

| parameter | default | unit | meaning |
|---|---|---|---|
| `vignette` | `0.15` |  | fractional loss of brightness at the corners |
| `tilt` | `0.03` |  | fractional left-right/top-bottom brightness gradient |

### Lamp / laser flicker  (`flicker`)

The excitation power fluctuates a little from frame to frame.

**Effect:** Every object in the channel gets brighter or dimmer together in a frame.

| parameter | default | unit | meaning |
|---|---|---|---|
| `cv` | `0.01` |  | relative size of the fluctuations (0.01 = 1%) |
| `memory` | `0.8` |  | how slowly the power changes between frames |


## Background

### Constant background light  (`background`, original model)

Stray light and camera-side background that adds the same number of photons to every pixel.

**Effect:** Raises the floor under everything and adds shot noise of its own (optics.background_photons).

### Autofluorescence  (`autofluorescence`)

Cells, medium and plastic fluoresce a little on their own, unevenly across the field.

**Effect:** A smooth, uneven background that is not the same at every position; hides faint spots.

| parameter | default | unit | meaning |
|---|---|---|---|
| `photons` | `6.0` | photons/px | mean extra background per pixel |
| `cv` | `0.3` |  | how much it varies across the field (relative) |
| `correlation_um` | `10.0` | um | size of the patches |

### Colour cross-talk (bleed-through)  (`crosstalk`)

A fluorophore's emission is not confined to its own channel: some of its light leaks into a neighbouring colour.

**Effect:** Faint ghost copies of bright objects (the whole nucleus, or a locus) in other channels.

| parameter | default | unit | meaning |
|---|---|---|---|
| `pairs` | `[['nucleus', 'locus0', 0.02], ['locus0', 'locus1', 0.02]]` | [from, to, fraction] | fraction of the first channel's light that also lands in the second |


## Detection

### Quantum efficiency  (`qe`, original model)

Only a fraction of the photons that reach the sensor release an electron (optics.quantum_efficiency).

**Effect:** Fewer detected events means more relative shot noise.

### Photon shot noise  (`shot_noise`, original model)

Light arrives in discrete photons, so the count in a pixel fluctuates (Poisson) even for a perfectly steady source.

**Effect:** Grainy images: noise grows as the square root of the signal; the dominant error for dim spots.

| parameter | default | unit | meaning |
|---|---|---|---|
| `excess_noise_factor` | `1.0` |  | extra multiplicative noise of EM-type cameras (1 = none, 1.41 = EMCCD) |

### Dark current  (`dark_current`)

Heat creates electrons in the sensor even in darkness.

**Effect:** A small added signal with its own shot noise; matters for long exposures.

| parameter | default | unit | meaning |
|---|---|---|---|
| `e_per_s` | `0.2` | e-/px/s | dark electrons per pixel per second |


## Camera electronics

### Hot pixels  (`hot_pixels`)

A few pixels always read too high (fixed positions).

**Effect:** Isolated bright dots that look like tiny spots to a naive detector.

| parameter | default | unit | meaning |
|---|---|---|---|
| `fraction` | `0.0003` |  | fraction of pixels that are hot |
| `e_per_frame` | `40.0` | e- | extra electrons in a hot pixel |

### Cosmic rays  (`cosmic_rays`)

High-energy particles occasionally hit the sensor and deposit charge in a few pixels.

**Effect:** Rare, bright, tiny events in single frames; they can be mistaken for a transient locus signal.

| parameter | default | unit | meaning |
|---|---|---|---|
| `rate_per_frame` | `0.05` | events | average events per frame |
| `e` | `3000.0` | e- | typical charge of one event |

### Pixel gain differences (PRNU)  (`prnu`)

Every pixel converts light to charge with a slightly different gain.

**Effect:** A fixed, speckle-like multiplicative pattern, proportional to the signal.

| parameter | default | unit | meaning |
|---|---|---|---|
| `cv` | `0.01` |  | relative spread of pixel gains (0.01 = 1%) |

### Read noise  (`read_noise`, original model)

Reading the pixel adds electronic noise independent of the signal (optics.read_noise_e).

**Effect:** Sets the noise floor in dark regions; dominates for very dim signals.

| parameter | default | unit | meaning |
|---|---|---|---|
| `pixel_cv` | `0.25` |  | sCMOS: relative spread of the read noise between pixels (0 = identical pixels) |

### Row (banding) noise  (`row_noise`)

CMOS readout adds a random offset to each row of pixels in each frame.

**Effect:** Faint horizontal stripes that change from frame to frame.

| parameter | default | unit | meaning |
|---|---|---|---|
| `adu` | `0.6` | ADU | size of the per-row offset |

### Offset pattern (DSNU)  (`dsnu`)

Each pixel has its own tiny fixed electronic offset.

**Effect:** A fixed additive pattern, visible only in dim regions.

| parameter | default | unit | meaning |
|---|---|---|---|
| `adu` | `1.0` | ADU | spread of the per-pixel offset |

### Digitisation and saturation  (`adc`)

The camera outputs whole numbers (ADU) with a limited range; very bright pixels clip at the full-well capacity.

**Effect:** Rounding to integers (coarse in dim images) and a flat ceiling on bright objects.

| parameter | default | unit | meaning |
|---|---|---|---|
| `bit_depth` | `16` | bits | number of bits of the converter |
| `full_well_e` | `30000.0` | e- | charge at which a pixel saturates |
