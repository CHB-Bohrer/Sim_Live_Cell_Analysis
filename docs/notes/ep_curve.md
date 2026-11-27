# E-P distance to P(burst ON) curve under measurement errors

Working notes moved verbatim from `CLAUDE.md` on 2026-11-27 (section order kept). Add new dated sections at the top; delete statements that are no longer true.

## E-P distance -> P(burst ON) curve under measurement errors (stage 7b, built 2026-10-14; docs/EP_CURVE.md)
- User's decisions (2026-10-14): readout = P(burst ON); fit the WHOLE curve (not one threshold) with much more data; enhancer with a longer working distance (hill, d_half 400 nm, n = 4,
  fold 12 on k_on; the old sharp contact was 150 nm); contact fraction 0.1 (exact: the stand-in distance process is solved so that P(d < 400 nm) = 0.100).
- Code: `stage7_validation/ep_curve.py` (`EPSetup`, `simulate_population`, `Measurement` / `measure_distance`, `kernel_matrix`, `estimate_on` (MS2 trace -> deconvolution -> burst rule),
  `fit_curve` / `analyse` / `reference_curve`), `scripts/sweep_ep_curve.py` (experiments ladder, sigma, shift, outliers, cells, null, dynamics, contact, frames, ms2, rules ->
  `$SIMLIVE_DATA/sweeps/ep_curve_<exp>/results.csv`, ~60 min for all), page 🔗 E–P curve (`app/ep_tab.py`), `tests/test_ep_curve.py` (6 tests).
- Method: errors-in-variables fit of p_far (1 + (F-1) / (1 + (d/d_half)^4)) integrating out the true distance (measurement kernel + Richardson-Lucy prior), cell bootstrap intervals,
  effect size = contrast P(ON at 150 nm)/P(ON at 1000 nm) (F alone is not identified at the null); reference = same fit to true distance + true ON of 4000 cells (contrast 4.07; 6.7 if the
  distance were fixed: the promoter's response delay flattens it).
- Findings (placeholders: invented distance dynamics, 57 nm per axis, 2D, 5% missed, 2% wrong picks): null false positive 3.3% (60 populations); 2D projection of the distance costs
  ~35% of the effect if ignored (2.7 vs 4.1) and is fully corrected by the EIV fit; localization precision matters only above ~200 nm; colour shift up to 100 nm harmless; missed frames
  harmless to 60%; WRONG PICKS are the bad one (5% -> -6%, 10% -> -13%, 20% -> -24% uncorrected; corrected since 2026-10-20 when their share is supplied, see below); cells: contrast sd 0.77 (50) -> 0.055 (3200), d_half scatter 11% -> 2%
  with a +3.6% offset; contact fraction 0.03 is biased low (2.7 vs 3.5). MS2 -> ON: the balanced-accuracy burst rule has 8% false positives next to real bursts, costs ~10% of the
  contrast, and the simple (sensitivity, false-positive) attenuation correction FAILS for it (9.2 vs 4.07); conservative rules (2,2), (3,3), (8,6) are nearly unbiased uncorrected
  and the correction is then acceptable (4.3-4.9). MS2 noise of 2 x one Pol II with 20 s frames halves the contrast. (An earlier note "intervals ~25% too narrow" was WITHDRAWN 2026-10-20: see below.)
- 2026-10-20 (docs/EP_CURVE.md "Update 2026-10-20"): (1) wrong-pick mixture `ep_curve.with_outliers` / `Measurement.assumed_outlier_rate` (kernel (1 - w) K + w q, q = a wrong
  pick's distance distribution; w must be SUPPLIED, it is not identifiable from distances): contrast 4.11-4.13 at 5-40% wrong picks with the true w (reference 4.07; uncorrected
  3.83 / 3.55 / 3.09 / 2.35); half the share corrects half, twice the share over-corrects and blows up at 40%. E-P page: checkbox "Analysis knows the wrong-pick share".
  (2) bootstrap replicates re-estimate the true-distance prior and refit on the full grid; the "25% too narrow" came from every scan setting reusing the SAME 20 populations:
  with 60 fresh populations coverage of the contrast is 0.95 (null: 6.7% exclude 1, nominal 5%); d_half intervals cover 88% (its +2-4% offset).
- NOT done: other coupling shapes; the real polymer-library distances (needs MSD calibration); linking sigma to photons / the movie's spot finder; a convergence
  check inside `estimate_on` (single chain); the real MS2 noise relative to one Pol II (user must supply).

### E-P curve with POLYMER distances and three microscopes (2026-10-14, second version; docs/EP_CURVE.md "Update 2026-10-14 (b)")
- User's decisions: 3D z-stacks AND biplane, compared (2D kept as a reference); D = 0.002 um^2/s is for a SINGLE LOCUS at a 1 MIN lag (assumed corrected for nucleus motion and drift; if it is for the E-P distance the block is 1.94 s instead of 0.46 s);
  cell type difference (paper = mouse ES cells) ignored for now. The earlier 2D baseline was MY choice (the repo's default configs and the Vutara preset are 2D "bypass"; z-stack config exists; biplane was "not built yet" and is now built at vector level).
- Paper checked (Jusuf et al. bioRxiv 2025, v1 main text only; Methods/supplement and the journal version not readable): mESC, Micro-C calibrated on BILD live imaging of 2 loops, "looping probability" = fraction of time looped ABOVE the random-polymer background; NO distances,
  no nm threshold, no diffusion numbers in what I could read. So nm per unit = 55.4 comes from the user's contact-fraction rule (P(d < 400 nm) = 0.100 for 100 kb pairs, pooled over the 10 validation sims), and block = 0.461 s from the MSD match (exponent 0.50 above one block, 0.57 below).
- Code: `stage7_validation/ep_polymer.py` (pool of 10 sims x 990 pairs, `calibrate`, `simulate_population`), `ep_acquisition.py` (`Acquisition` kinds 2d / zstack / biplane, `measure`, `kernel_aniso`: errors from Thompson / axial Fisher information, exposure blur, the 0.15 s colour delay, piezo,
  all with the extrapolated sub-block MSD), `ep_history.py` (P(ON) from the distance HISTORY), `scripts/sweep_ep_polymer.py` (experiments modality, null, wrong, cells, contact, frames, colour, photons, focus, timescale, compare_sources, ms2 -> `$SIMLIVE_DATA/sweeps/ep_polymer_<exp>`),
  page 🔗 E–P curve part 1 (live polymer + microscope), `tests/test_ep_polymer.py` (5 tests; skipped without the validation sims).
- FINDINGS: polymer E-P distance decorrelates in 9.2 s but the promoter integrates over ~45 s => P(ON | distance NOW) is nearly flat (contrast 1.36 vs 4.1 for the invented slow process); the curve to fit is P(ON) for a distance HELD, from the distance history (reference 4.19, tau 45 s found). Error budget per axis (nm): localization 6-9 in
  plane / 27 (z-stack) or 51 (biplane) in z, exposure motion 38-62, sequential colours 0.15 s apart = 120 (the biggest), piezo 12; total ~126-140. At 600 cells x 1 h: contrast 4.77 +- 1.00 (2D), 4.25 +- 0.72 (z-stack), 4.45 +- 0.80 (biplane), truth 3.80 +- 0.54; 95% interval width 5.1 / 3.9 / 3.1 / 2.3;
  null false positives 3-7%. Need ~10x more cells than with the invented process; frames every <= 10 s (reference falls 4.19 -> 1.78 at 80 s); biplane loses loci when focus is off by > 0.6 um or photons < 0.3x; simultaneous colours cut ~20% of the scatter; the ON state from the MS2 trace is biased high (calls cluster near the enhancer) and unstable.
- NOT done: wrong-pick filter from the locus's past positions (the spot finder already links with a diffusion model; the scan's wrong picks are random draws); nucleus "steadiness" (simulated nuclei: outline changes 0.29 um rms per 10 s frame, long/short radius 1.6, rotation 4.5 deg and centroid step 0.9 um per frame; awaiting which meaning the user intends);
  rendering real z-stack / biplane movies and localizing in them (the acquisition model is vector-level); loop-extrusion timing is too fast at 0.46 s per block (loop lifetime 46 s); the E-P polymer pool still uses the 10 validation simulations (the 100-simulation library exists since 2026-10-18 and could replace it).
