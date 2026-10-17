# Instruments: what optics the simulation assumes

Until 2026-10-13 every config used a **generic placeholder** microscope: widefield, NA 1.2 water, 130 nm pixels, QE 0.8, read noise 1.5 e-,
15% vignetting, and 2% bleed-through of the nuclear GFP into the promoter's colour. None of these came from a real instrument.

An instrument preset now replaces them: `--set instrument=<name>` (or `instrument: <name>` in a config), files in `configs/instruments/`.
Explicit `--set` values still win.

## Bruker Vutara VXL (`vutara_vxl`), widefield (bypass) mode

| Setting | Value | Source |
|---|---|---|
| Objective | 60x, NA 1.3, silicone-oil immersion | [Bruker Vutara VXL page][b] |
| Immersion index | 1.406 | **ASSUMED** (silicone oil; not on the Vutara pages) |
| Camera | Orca Fusion BT sCMOS (optional; Orca Flash 4.0 V3 standard), QE 95% | [Bruker][b], [News-Medical 2025][nm] |
| Pixel size in the sample | ~100 nm | [News-Medical 2025][nm] |
| Read noise | 0.7 e- | **ASSUMED** (Fusion BT slow scan; check the camera sheet) |
| Illumination | top-hat (flat) from a square fibre | [Bruker][b] |
| Uneven illumination | vignette 3%, tilt 1% | **ASSUMED** small because of the top-hat |
| Lasers | 405, 488 (optional), 555, 638, 750 (optional) nm | [Bruker][b] |
| Colours | sequential, up to 5; simultaneous 2 only with an optional second camera | [Bruker][b] |
| Focus | active focus correction, drift < 30 nm over 10 min → defocus rms 30 nm | [Bruker][b] |
| Field of view | 200 x 200 µm widefield / bypass; 50 x 50 µm biplane | [Bruker][b] |
| Biplane | two planes 600 nm apart in the sample, 50:50 beam splitter, one camera (two halves) | [Bruker bi-plane page][bp]; built, see below |
| Dyes | GFP (nucleus, MCP) at 488 nm; promoter label in the 555 nm channel; enhancer label in the 638 nm channel | **ASSUMED** (which labels do you use?) |
| Time between colours | 0.15 s (0.1 s exposure + switching), order promoter → enhancer → GFP | **ASSUMED** |

### Bleed-through, from the spectra (not from a guess)

With sequential excitation, the 555 nm image is taken while the 488 nm laser is off. GFP therefore appears in it only if the 555 nm laser
excites it, and only the part of its emission that passes the red filter. From the [FPbase EGFP spectra][fp]:

- **Cross-excitation:** at 555 nm, 0.0042 of GFP's peak.
- **Emission share:** 6.2% of GFP's emission in a 570-620 nm band, vs 71.7% in 500-550 nm.
- **Bleed-through:** 0.0042 x 0.062 / 0.717 = **0.036%** of the GFP signal (equal laser power and exposure). That is 55x less than the
  old 2% placeholder.
- **If both colours were excited at once** with the same filters, it would be 0.062 / 0.717 = **8.6%**. That is the case for
  simultaneous two-camera detection.
- **GFP into the 638 nm channel:** 0.0023 x 0.006 / 0.717 ≈ 0.002%.
- **Promoter dye into the 638 nm channel:** assumed 0.2%; it depends on the red dye.

The price of sequential excitation is time. The colours are not taken at the same moment, and chromatin moves in between. This is the new
error source `sequential_colour`: it shifts each later colour by the chromatin motion over its delay (a Brownian bridge between frames;
the cell's slow drift is interpolated). No calibration can remove it. It can only be kept small:
- short exposures and fast switching;
- the two locus colours next to each other in the order.

The measured trade-off (sequential vs simultaneous on the same movie) is in [ERROR_EFFECTS.md](ERROR_EFFECTS.md).

### Not modelled yet
- In biplane mode: the true plane separation and the z calibration of the PSF (a bead z-scan) are taken as nominal; the camera halves' registration is a shift only (no rotation or scale).
- The Orca Flash 4.0 V3 alternative.
- Measured flat-field and read noise of your instrument.
- The real filter set and dyes. With your filter bands and dyes the bleed-through numbers can be recomputed the same way.

## Bruker Vutara VXL (`vutara_vxl_biplane`), biplane 3D mode (built 2026-10-17)
Same instrument and numbers as above, imaged the way Bruker describes the bi-plane mode: "simultaneous detection from two axially separated object
planes" on the two halves of one camera, planes **600 nm apart in the sample**, a 50:50 splitter, and the height of a spot from the intensity
difference between the two planes ([bp]). Bruker's page gives no axial range for one acquisition and says z-stacks are possible to more than 50 µm deep:
a nucleus (~10 µm) is much thicker than one pair of planes, so the stage still steps through it.

| Setting | Value | Source |
|---|---|---|
| Planes per stage position | 2, 0.6 µm apart, each with half of the light (`photon_fraction: 0.5`) | [Bruker bi-plane page][bp] |
| Stage step | 1.2 µm = 2 x the separation, so the planes tile the sample evenly every 0.6 µm (and the 3D localizer, which assumes even planes, is unchanged) | design choice |
| Stage positions | 11 (22 planes, z = 0.0 to 13.2 µm; the lowest pair straddles the fiducial beads at z = 0.3 µm) | **ASSUMED** |
| Time per stage position | 0.15 s | **ASSUMED** |
| Field of view | 40 x 40 µm by default (Bruker: 50 x 50 µm in biplane mode) | [Bruker][b] |
| Registration of the two camera halves | second plane shifted by (0.35, -0.25) px; its true height 30 nm above nominal | **ASSUMED** (new error source `biplane_registration`) |

Used by `--set instrument=vutara_vxl_biplane` or the sidebar's Microscope choice (it sets its own z-stack, beads and field). The piezo error and the
slice skew apply per stage position (both planes of a pair share them). Images per run are large: 22 planes x 3 colours x 400 x 400 pixels x 2 bytes.

**Calibration with the beads in the image** (`stage6_analysis/biplane.py`): a fiducial bead lies between the two planes of the lowest pair, so it
is bright in both; the difference of its two lateral centroids, median over beads and frames, is the registration offset. The locus finder then shifts
the second plane of every pair back before it looks for spots (`plane_shift_px`); the pipeline reports this as `truth_3d_registered` next to the
uncorrected `truth_3d`, and the estimate against the truth as `biplane_registration` in `scores.json`. Not measured by the beads: the plane separation
(one bead intensity ratio is one number for two unknowns) and the PSF z calibration, which needs a z-scan of a bead.
First measurement (seed 61, 2 cells, 12 frames, one run, so a demonstration, not evidence): offset estimated (0.321, -0.214) px for a true (0.35, -0.25);
the locus error barely changed (median 3D error 0.220 -> 0.213 µm) because it is dominated by the placeholder axial chromatic shift (median z error
-0.20 µm), which the dark / flat / bead-slide calibration (`calibration.simulate`, switched off in that run) is meant to remove. The same scene as an ordinary
22-plane z-stack gave the same 3D error (0.27 µm rms) and a larger lateral rms (0.10 vs 0.06 µm).

[b]: https://www.bruker.com/en/products-and-solutions/fluorescence-microscopy/super-resolution-microscopes/vutara-vxl.html
[bp]: https://www.bruker.com/en/products-and-solutions/fluorescence-microscopy/super-resolution-microscopes/3d-imaging-with-bi-plane.html
[nm]: https://www.news-medical.net/news/20250820/How-the-Vutara-VXL-simplifies-single-molecule-imaging.aspx
[fp]: https://www.fpbase.org/protein/egfp/
