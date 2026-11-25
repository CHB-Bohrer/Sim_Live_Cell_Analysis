# Cells and the microscope {#microscope}

::: summary
**In one sentence.** Nuclei move, bump into each other, change shape and divide. A forward model of the microscope then turns the true positions into camera images, with up to 28 kinds of imperfection that can each be switched on or off.
:::

## Cells

Each nucleus has its own ID, which it keeps for the whole movie. Nuclei drift and diffuse, collide, deform and rotate, and can divide. They can also leave and re-enter the field of view. The chromatin from chapter 1 is placed inside each nucleus, and the loci move with it. The speed is set by one number, `speed_scale`. The real cell speed, density and frame interval are still needed from experiments.

## The microscope

Each frame is drawn in the order the physics happens, from light to camera:

- **Optics:** stage drift, chromatic shift between colours, the point spread function (PSF), defocus and haze.
- **Sample and light:** exposure blur from motion, photobleaching (per cell, by its light dose), uneven illumination, flicker, background, autofluorescence, a cytoplasm glow and cross-talk between colours.
- **Camera:** quantum efficiency, shot noise, dark current, hot pixels, cosmic rays, per-pixel gain and offset patterns, read noise and row noise, and rounding to integers with saturation.
- **3D:** z-stacks with piezo errors, and the Bruker Vutara biplane mode (two planes 600 nm apart), with presets for that instrument.
- **Stray spots:** spots that look like a locus, but are not one.

Every error source has its own random-number stream. So switching one source off never changes the noise of another, and any frame can be redrawn bit for bit.

## How we know it works

Unit tests check each piece against theory. For example, the camera noise matches the photon-transfer formula, defocus conserves photons, and drift moves the spots, nuclei and truth together. A stronger end-to-end check is **blind calibration**. We simulate the calibration images a lab would take with the same virtual instrument: dark frames, flat-field slides and beads. Then we estimate the instrument from those images alone, without looking at the true settings.

![**Blind calibration recovers the instrument.** Each row is an estimated property, as estimate divided by true value; grey is ±5%. This figure uses *short* test acquisitions (40 dark frames, 30 beads), so the rare or slow quantities are rough. Cosmic-ray rate is based on a handful of events, and flicker on 80 frames. With the full default acquisitions, 53 of 53 quantities are within tolerance: gain within 0.3%, read noise within 0.1%, every hot pixel found with none false, and colour shift within 0.001 px. Made by `figure_scripts/fig_calibration.py`.](figures/calibration.png){width=85%}

Applying the calibration to a movie cut the locus error a lot. The median fell from 54 to 27 nm, and the 95th percentile from 4.7 µm to 76 nm.

## Limits

- **Placeholder sizes:** every error size is a placeholder until real dark, flat-field and bead images are available.
- **Optics:** colour shifts are the same everywhere in the field, and there is no lens distortion.
- **Motion blur:** nuclei are not motion-blurred (this is negligible at these exposures).
