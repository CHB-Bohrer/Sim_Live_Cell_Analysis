# Instruments: what optics the simulation assumes

Until 2026-10-13 every config used a **generic placeholder** microscope: widefield, NA 1.2 water, 130 nm pixels, QE 0.8, read noise 1.5 e-,
15% vignetting, and 2% bleed-through of the nuclear GFP into the promoter's colour. None of these came from a real instrument.

An instrument preset now replaces them: `--set instrument=<name>` (or `instrument: <name>` in a config), files in `configs/instruments/`.
Explicit `--set` values still win.

## Bruker Vutara VXL (`vutara_vxl`), 2D bypass mode

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
| Biplane (not built yet) | two planes 600 nm apart in the sample, 50:50 beam splitter, one camera | [Bruker bi-plane page][bp] |
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
- Biplane 3D: two simultaneous planes, half the photons each.
- The Orca Flash 4.0 V3 alternative.
- Measured flat-field and read noise of your instrument.
- The real filter set and dyes. With your filter bands and dyes the bleed-through numbers can be recomputed the same way.

[b]: https://www.bruker.com/en/products-and-solutions/fluorescence-microscopy/super-resolution-microscopes/vutara-vxl.html
[bp]: https://www.bruker.com/en/products-and-solutions/fluorescence-microscopy/super-resolution-microscopes/3d-imaging-with-bi-plane.html
[nm]: https://www.news-medical.net/news/20250820/How-the-Vutara-VXL-simplifies-single-molecule-imaging.aspx
[fp]: https://www.fpbase.org/protein/egfp/
