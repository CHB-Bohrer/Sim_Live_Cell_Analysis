<!--
Section: Microscope (stage 3)
Update when: an image-error source is added or changed, or the z-stack, biplane, calibration or instrument presets change.
Code: src/simlive/stage3_microscopy, configs/instruments
Check: none
Sources: docs/IMAGE_ERRORS.md, docs/CALIBRATION_AND_3D.md, docs/INSTRUMENTS.md, docs/notes/imaging.md.
-->

# The microscope: from truth to camera counts

::: {.box}
**Short version.** A forward model turns positions and brightness into camera images. It contains 28 sources of image error, each of which can be switched on alone. Because the same instrument can also "image" dark frames, flat slides and beads, the pipeline can calibrate itself the way a real lab would and then check the calibration against the truth.
:::

## What a frame is made of

For every channel (the nuclear channel and one colour per locus) the `Imager` builds a frame in physical order: where things truly are, where the stage has drifted, blur by the point-spread function, illumination, background and autofluorescence, cross-talk between colours, then the camera: quantum efficiency, photon shot noise, dark current, hot pixels, cosmic rays, pixel-to-pixel gain, read noise, row noise, offset pattern and digitization. Electrons are drawn as $\text{Poisson}(\text{QE} \times \text{photons})$.

The 28 sources are one table (`errors.py`); `docs/IMAGE_ERRORS.md` is generated from it. Later additions include slice skew, the piezo stage, sequential colour acquisition, cytoplasm, stray spots and biplane registration.

Two properties make the model usable for experiments:

- **Each source has its own random stream.** Frame $t$ is identical whatever else is rendered, and switching one source never changes another source's noise. The dashboard rebuilds a saved frame from the settings and confirms it matches bit for bit.
- **Instrument and session are separate.** Fixed patterns (hot pixels, gain map, flat field) belong to the instrument. Drift, focus and flicker belong to the session. A calibration acquisition shares the instrument but has its own session, as in a real lab.

Truth is saved next to the images: label masks, drift, and each locus both where it truly is and where it is drawn (after chromatic shift). Analysis is scored against the true position, so the chromatic shift counts as error and the drift does not.

## Three dimensions

With z-stacks on, cells and loci have heights, the point-spread function is three-dimensional, and noise is drawn per plane. The piezo model gives each plane a true height that differs from the commanded one (random error, bow, hysteresis, settling lag). Beads on the slide appear in every image.

**Biplane mode** models the Bruker Vutara VXL: two planes 600 nm apart, each receiving half the light, recorded at the same instant. A z-stack is still needed to cover a nucleus about 10 µm tall.

## Instruments

A preset (`configs/instruments/*.yaml`) sets optics and camera for a named microscope. The Vutara presets use a 60x / 1.3 NA silicone objective, about 100 nm pixels and a camera with quantum efficiency 0.95. Each value in `docs/INSTRUMENTS.md` is marked as sourced or assumed. Bleed-through between colours is computed from published spectra: 0.036% with sequential colours, about 8.6% if both were excited at once.

## Calibration, and whether it works

The same instrument images dark frames, a flat slide at several light levels and a bead slide. Blind estimators then recover gain, read noise, hot pixels, flat field, flicker, chromatic shift, point-spread width and drift.

Recorded result (placeholder instrument, one seed): **53 of 53 compared quantities within tolerance**, for example gain within 0.3%, read noise within 0.1%, 40 of 40 hot pixels found with none false, vignetting 0.1496 against a true 0.15. Applying the calibration to a 2D movie reduced the median locus error from 54 to 27 nm and the 95th percentile from 4.7 µm to 76 nm.

A point the model makes concrete: beads on the slide cannot measure plane-height errors at the height of the nucleus. Only the piezo's own position log can.

All instrument values are placeholders until real dark, flat and bead images are supplied **[Assumption]**.
