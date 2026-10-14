# P(burst ON) as a function of the promoter-enhancer distance, and what measurement errors do to it (2026-10-14)

Question (user): readout = P(promoter in a burst) vs E-P distance, the WHOLE curve (not one contact threshold); enhancer with a long working distance;
contact fraction 0.1; much more data; errors varied. Code: `src/simlive/stage7_validation/ep_curve.py`, scans `scripts/sweep_ep_curve.py`
(`$SIMLIVE_DATA/sweeps/ep_curve_<exp>/results.csv`, also readable on the 📈 Scans page), dashboard page 🔗 E–P curve (`app/ep_tab.py`), tests `tests/test_ep_curve.py`.

## The set-up (everything here is simulation; the numbers that stand in for the real experiment are PLACEHOLDERS)
- **Biology.** Coupling `hill`: f(d) = 1 / (1 + (d / 400 nm)^4) multiplies k_on by 1 + 11 f(d) (fold 12). The previous sharp contact was 150 nm; 400 nm is the
  "working distance" (f = 0.5). Gene: k_on 1/1200 s, k_off 1/45 s, k_init 0.15/s (bursts of ~7 Pol II), 24 loops, 2.5 kb/min, 1.3 + 2 kb.
- **Distance process (INVENTED, like toy_loci).** Hidden looped / open switch, mean 240 s per visit, fluctuations relaxing in 30 s; the looped fraction is solved so that
  the fraction of time with d < 400 nm is **exactly 0.100** (measured 0.0998 in 4000 cells). With the polymer library, d(t) will come from the simulations once nm_per_unit and
  the block duration are calibrated; at the moment the placeholder scale almost never reaches contact.
- **Data.** 600 cells x 1 h at 10 s frames (216,000 frames; ~9% of them ON). The reference is the same fit to the TRUE distance and TRUE ON state of 4000 cells.
- **What the curve contains.** The promoter responds to the distance with a delay (OFF -> ON takes ~100 s even in contact; bursts outlast the visit), so P(ON | d now) is
  flatter than the instantaneous k_on / (k_on + k_off): P(ON) at 150 nm / at 1000 nm is **4.1** in the simulated population and would be 6.7 if the distance were held fixed.
  The recovered quantity is the relation present in the population, dynamics included.

## Method
1. **Measurement** (`measure_distance`): 3D difference vector = true distance x random orientation + colour shift + Gaussian localization noise per axis; distance in 2D (image
   plane; shortens distances by pi/4 on average and blurs the curve) or 3D; missed frames; wrong picks (uniform in 0-2000 nm).
2. **ON state.** L2 uses the true state; L1 / L3 estimate it: observed MS2 trace (noise as the real pipeline would estimate it) -> Bayesian deconvolution (single gibbs_shift
   chain per cell) -> burst rule -> per-frame P(ON) = fraction of posterior samples called ON.
3. **Curve fit** (`fit_curve`): P(ON | d) = p_far (1 + (F - 1) / (1 + (d / d_half)^n)), n = 4. Errors in variables: the latent true distance is integrated out. A
   measured-distance bin b maps to P(true d | b) through the measurement kernel (non-central chi / Rice, orientation- and shift-averaged; `kernel_matrix`) and a
   population prior of true distances estimated from the measured ones (Richardson-Lucy). Grid search over (d_half, F, p_far) + Nelder-Mead polish. `naive=True` takes the
   measured distance as the true one.
4. **Intervals:** percentile bootstrap over CELLS (100 refits; the prior is kept fixed).
5. **Effect size = contrast** P(ON at 150 nm) / P(ON at 1000 nm). The fold F alone is not identified when there is no effect (any F goes with a tiny d_half); the contrast is.
   Also reported: d_half (working distance), p_far, and curve_err = RMS difference of the whole curve from the reference, weighted by how often each distance occurs.
6. **Checks (tests):** kernel columns sum to 1 and match a Monte-Carlo measurement; the contact fraction is 0.10; a known curve (fold 8, d_half 350 nm) is recovered through 70 nm 2D noise and the
   EIV fit beats the naive one; with coupling off the 95% interval contains 1; an attenuation correction undoes a known sensitivity / false-positive rate.

## Results (20 populations per setting; contrast reference = 4.07)
"covers" = the 95% interval contains the reference. Baseline measurement (placeholders): 57 nm per axis, 2D, 5% missed, 2% wrong picks, no colour shift.

**Null (coupling off, 60 populations):** contrast 1.015 +- 0.11; the interval excludes 1 in 3.3% of populations (nominal 5%). With coupling on (60 populations): the interval excludes 1 every time, coverage of the reference 93%.

**Localization precision and 2D vs 3D** (contrast; EIV / naive):

| sigma per axis (nm) | 2D EIV | 2D naive | 3D EIV | 3D naive |
|---|---|---|---|---|
| 10 | 3.99 | 2.72 | 3.95 | 3.96 |
| 57 | 4.00 | 2.75 | 3.95 | 3.99 |
| 100 | 4.00 | 2.66 | 3.95 | 3.94 |
| 200 | 4.01 | 2.17 | 3.81 | 3.05 |
| 400 | 3.15 | 1.42 | 3.25 | 1.50 |

- **Measuring the distance in 2D is the big systematic error** (contrast 2.7 instead of 4.1; curve error 10x larger); the errors-in-variables fit removes it completely up to
  ~200 nm precision. Localization noise alone matters only above ~200 nm: the distances are 100s of nm, the working distance is 400 nm.
- At 400 nm precision the effect is still detected but halved or worse and the intervals are wide.

**Colour shift left over** (2D): up to 100 nm harmless (contrast 3.99-4.02, whether or not the analysis knows it); at 200 nm uncorrected the contrast drops to 3.79 (coverage 70%), corrected 4.09.

**Wrong picks and missed frames:** missed frames are harmless up to 60% (contrast 4.0; interval width 0.64 -> 0.73). Wrong picks are not: 5% -> contrast 3.83 (coverage 55%), 10% -> 3.55 (10%),
20% -> 3.09 (0%), 40% -> 2.35. A wrong-pick mixture in the fit would remove most of this; not built.

**Number of cells** (L2; sd of the contrast over populations): 50 -> 0.77, 100 -> 0.55, 200 -> 0.36, 400 -> 0.26, 800 -> 0.17, 1600 -> 0.10, 3200 -> 0.055 (the 1/sqrt(n) of the cells).
The effect (contrast 4) is detected at every size from 50 cells; what many cells buy is the SHAPE: curve error 0.010 (50) -> 0.002 (3200); the scatter of the recovered working distance d_half falls from 11% (50 cells) to 4% (600) and 2% (3200). A small systematic offset remains: d_half 336 nm against 324 (+3.6%, the fitted hill shape is not exactly the shape the dynamics produce), so the d_half interval covers the reference less often as the cell number grows (95% at 50-100 cells, 80% at 600, 55% at 3200).

**Dynamics and design:** time per visit 60 / 240 / 960 / 3600 s: reference contrast 3.0 / 4.1 / 4.5 / 4.6 (longer visits let the promoter follow the distance; the recovery is unbiased in all).
Frame interval 10 / 20 / 40 / 80 s: contrast unchanged, interval width 0.64 / 0.69 / 0.77 / 0.91 (with the 1 h movie fixed). Contact fraction 0.03 / 0.1 / 0.2: reference 3.5 / 4.1 / 3.8, recovered
2.69 / 3.99 / 3.73: at 3% contact the estimate is biased low (coverage 45%): few contact frames and the fit shape cannot be pinned.

**From the MS2 trace to P(ON)** (L1 = true distance, L3 = measured; burst rule chosen on the simulated truth: best balanced accuracy):

| MS2 noise (sigma_b, units of one Pol II) | frame 10 s: contrast | frame 20 s: contrast | sensitivity / false positives |
|---|---|---|---|
| 0.05 | 3.68 | 3.63 | 10 s: 0.97 / 0.084; 20 s: 0.95 / 0.090 |
| 0.4 | 3.60 | 3.29 | 0.96 / 0.090; 0.94 / 0.109 |
| 1.0 | 2.82 | 2.58 | 0.95 / 0.138; 0.93 / 0.162 |
| 2.0 | 3.31 | 1.97 | 0.90 / 0.104; 0.89 / 0.236 |

- With a clean trace the burst call costs ~10% of the contrast (3.7 vs 4.07, curve error 0.08 vs 0.003): the call is slightly wrong around bursts (timing smear), and the false positives raise the far-distance level.
  The error that matters is the noise and the frame interval: at noise 2 x one Pol II and 20 s frames the contrast is halved; 95% intervals are then too narrow (coverage 0-50%) because they do not contain this bias.
- **The simple attenuation correction (observed = fp + (se - fp) P(ON)) FAILS for the balanced-accuracy rule:** contrast 9.2-10 instead of 4.07. The false positives are not independent of the
  distance (they sit next to real bursts, which are more frequent near the enhancer), so subtracting one constant overcorrects. With conservative rules (low false positives) it is close:

| burst rule (window frames, starts needed) | sensitivity | false positive | contrast uncorrected | corrected | curve error uncorrected / corrected |
|---|---|---|---|---|---|
| 2, 2 | 0.54 | 0.005 | 3.97 | 4.39 | 0.043 / 0.007 |
| 3, 3 | 0.53 | 0.004 | 3.96 | 4.28 | 0.045 / 0.006 |
| 8, 6 | 0.63 | 0.011 | 4.14 | 4.94 | 0.026 / 0.014 |
| 5, 1 (balanced accuracy) | 0.97 | 0.084 | 3.68 | 9.20 | 0.083 / 0.046 |
| 5, 2 | 0.89 | 0.037 | 4.00 | 6.27 | 0.030 / 0.025 |

  A conservative rule gives an almost unbiased contrast (uncorrected 3.96-4.14) with more scatter (sd of the contrast over populations 0.31-0.36 against 0.27 for the balanced rule, 4 populations each), and with the correction the absolute P(ON) level is right too (curve error 0.006-0.007).
  For real data the rule cannot be tuned on the truth: choose a conservative one from the start rate (page 🧮 MS2 → Pol II, step 6).

**Calibration of the intervals:** the cell-bootstrap intervals are ~25% too narrow (mean width 0.64 against 3.92 x sd over populations 0.19 = 0.75): coverage 85-90% at the nominal 95% in most scans.
Widen them by ~1.3 for honest 95%; not corrected in the code.

## What this means
1. With ~600 cells x 1 h and a 400 nm working distance the CURVE is recoverable: the contrast to 5% (sd over populations; interval +-8%), the working distance to 4% scatter (+3.6% offset), when the true ON state is known and the fit accounts for the measurement.
   Ignoring that the distance is measured in the image plane costs ~35% of the effect (2.7 instead of 4.1).
2. The measurement error that matters most is not localization precision (fine to ~200 nm): it is (i) 2D projection (fixable in the fit), (ii) wrong picks above ~5% (not fixed), (iii) the MS2-based ON
   estimate when the trace is noisy and frames are slow (not fixed; needs the user's real MS2 noise relative to one Pol II).
3. More cells narrow the intervals as 1/sqrt(n); the movie length and the time per visit set how many bursts occur in contact (contact fraction 0.03 is much harder than 0.1).

## Not done / caveats
- The distance process is invented; the real dynamics (library, calibrated) change the reference curve and the delays. The 2D projection assumes isotropic orientation.
- Hill shape with n = 4 assumed in the fit (the truth has it); a free n is implemented (`free_shape=True`) but not scanned. Other couplings (sharp contact, exp) not scanned.
- The localization noise is a parameter, not yet linked to the movie pipeline's spots (photons, PSF); the full-image pipeline (spot finding at the real positions) was not run for this scan; the measured
  quantities are simulated at the level of the distance vector. The earlier spot-finding scores (22-84 nm median error, 98% found) put the baseline in the right range.
- One gibbs_shift chain per cell for the ON estimate (no convergence check), 200 samples; a single gene (placeholder kinetics). MS2 noise, burst rule, frame interval grids are small (4 populations per setting).
- Homologous copies / sister chromatids, a shared probe pool and transcription feeding back on the polymer are not modelled.
