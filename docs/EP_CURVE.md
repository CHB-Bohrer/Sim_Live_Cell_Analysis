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
20% -> 3.09 (0%), 40% -> 2.35. **Fixed 2026-10-20** by a wrong-pick mixture in the fit when the share is supplied (see "Update 2026-10-20" below):
corrected 4.11-4.13 at 5-40%.

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

**Calibration of the intervals:** this said "~25% too narrow (coverage 85-90%)". **Withdrawn 2026-10-20:** every setting of a scan reuses the SAME 20 populations,
so the "consistent" shortfall was one sample of 20; with 60 fresh populations the contrast intervals cover the reference 95% of the time (see "Update 2026-10-20").

## What this means
1. With ~600 cells x 1 h and a 400 nm working distance the CURVE is recoverable: the contrast to 5% (sd over populations; interval +-8%), the working distance to 4% scatter (+3.6% offset), when the true ON state is known and the fit accounts for the measurement.
   Ignoring that the distance is measured in the image plane costs ~35% of the effect (2.7 instead of 4.1).
2. The measurement error that matters most is not localization precision (fine to ~200 nm): it is (i) 2D projection (fixable in the fit), (ii) wrong picks above ~5% (corrected since 2026-10-20 when their share is known), (iii) the MS2-based ON
   estimate when the trace is noisy and frames are slow (not fixed; needs the user's real MS2 noise relative to one Pol II).
3. More cells narrow the intervals as 1/sqrt(n); the movie length and the time per visit set how many bursts occur in contact (contact fraction 0.03 is much harder than 0.1).

## Not done / caveats
- The distance process is invented; the real dynamics (library, calibrated) change the reference curve and the delays. The 2D projection assumes isotropic orientation.
- Hill shape with n = 4 assumed in the fit (the truth has it); a free n is implemented (`free_shape=True`) but not scanned. Other couplings (sharp contact, exp) not scanned.
- The localization noise is a parameter, not yet linked to the movie pipeline's spots (photons, PSF); the full-image pipeline (spot finding at the real positions) was not run for this scan; the measured
  quantities are simulated at the level of the distance vector. The earlier spot-finding scores (22-84 nm median error, 98% found) put the baseline in the right range.
- One gibbs_shift chain per cell for the ON estimate (no convergence check), 200 samples; a single gene (placeholder kinetics). MS2 noise, burst rule, frame interval grids are small (4 populations per setting).
- Homologous copies / sister chromatids, a shared probe pool and transcription feeding back on the polymer are not modelled.

## Update 2026-10-20: wrong picks corrected in the fit; the interval check redone
- **Wrong-pick mixture** (`ep_curve.with_outliers`, `Measurement.assumed_outlier_rate`): the kernel becomes K' = (1 - w) K + w q, q = the measured-distance distribution of a
  wrong pick (uniform 0-2 um here, the simulation's model), used for the population prior and the EIV weights alike. A distance that may be a wrong pick is then trusted less
  automatically. w cannot be learned from the distances (a wider true-distance distribution explains them as well), so the analysis must supply it (spot finder statistics, or
  a simulation of the same movie). Scan `outliers` (20 populations, reference 4.07; column `assumed` = supplied share / true share):

  | wrong picks | not corrected | half the share | the true share | twice the share |
  |---|---|---|---|---|
  | 5% | 3.83 | 3.98 | 4.12 | 4.41 |
  | 10% | 3.55 | 3.83 | 4.11 | 4.74 |
  | 20% | 3.09 | 3.54 | 4.11 | 5.92 |
  | 40% | 2.35 | 2.97 | 4.13 | 52.9 (fails) |

  With the true share the bias is gone at every rate (coverage 0.80-0.90 with 20 populations). Overstating the share is worse than understating it: at high rates it blows up.
- **Bootstrap:** every replicate now re-estimates the true-distance distribution and refits on the full grid (`bootstrap_fits(refit_prior=True, coarse=False)`). Diagnosis
  on 30 fresh populations (`scripts` not kept; seeds 3000+): bootstrap SD / real SD of the contrast = 0.96 (true distances) and 0.99 (measured, new; 1.05 with the old
  coarse refits). Scan `null` (60 populations per arm, new bootstrap): coupling on: contrast 3.99 vs 4.07, interval-implied SD 0.155 vs real 0.143, **coverage 0.95**;
  coupling off: 6.7% of intervals exclude "no effect" (4 of 60; nominal 5%). The d_half intervals cover 88% (the +2-4% d_half offset above). Old results kept as
  `results_before_2026-10-20.csv` in `$SIMLIVE_DATA/sweeps/ep_curve_{outliers,null}`.

---

# Update 2026-10-14 (b): E-P distances from the polymer, three microscope modes

(The sections above describe the first version: an invented distance process, 2D noise only. They stay as they are; the comparison is in "What changed".)

## Where the numbers come from
- **Literature (checked, not used for a scale).** Jusuf, Grosse-Holz, Gabriele, Mach, Flyamer, Zechner, Giorgetti, Mirny, Hansen, "Genome-wide absolute quantification of chromatin looping" (bioRxiv 2025.01.13.632736, v1 of 15 Jan 2025).
  I read the main text; its Methods and supplement were not in the file, and the journal version (Nat Struct Mol Biol 2026) could not be opened. Mouse ES cells; deep Micro-C calibrated against earlier BILD live-imaging of two loops
  (Fbn2 505 kb, synthetic Npr3 164 kb); DNA FISH for four E-P pairs. "Looping probability" = the average fraction of alleles in which the loop exists, or the fraction of time looped at one allele, counting only interactions
  ABOVE the background of random polymer encounters. No nm threshold is given in the main text; the capture / interaction radii of its simulations are in the unread supplement. Reported: mean 2.3%, median 1.6%, max 26% over 36,804
  loops; CTCF-CTCF 3.2% vs cis-regulatory 1.1%; DNA FISH E-P: Sox2 19.8%, Mycn 12.7%, Klf4 5.8%, Nanog 5.3%; calibration slope k = 0.186 +- 0.047; CTCF loop lifetime ~20 min (cited from earlier work). **No distances in nm,
  no distance distribution at 100 kb, no diffusion coefficient.** So the length scale is NOT from the paper. (Search-summary numbers for the later version, 65,929 loops / mean 1.2%, are unverified.)
- **Length (nm per polymer unit) = 55.4**, chosen so that the E-P pair (100 kb apart) is within the 400 nm working distance 10% of the time (the user's contact fraction). With it the pooled E-P distance has median 792 nm, mean 827 nm, 5-95%:
  315-1458 nm. These distances cannot be checked against the paper. Its "looped" fraction (excess above background, 1-2% for a typical E-P pair) is a different quantity from "within 400 nm" (10%), which includes random encounters.
- **Time (seconds per saved block) = 0.461 s**: the single-locus MSD (3D, centre of mass removed) at 60 s is matched to 6 D t with D = 0.002 um^2/s (user: one locus, lag 1 min, assumed corrected for nucleus motion and drift). The simulated MSD exponent
  is 0.50 above one block (0.57 below it). If D referred to the E-P distance instead, a block would be 1.94 s. The choice of lag matters a lot (seconds per block at 10 / 30 / 60 / 120 s lag: 3.2 / 1.1 / 0.46 / 0.30 s).
- **Checks before the scans:** contact fraction of the pooled polymer pairs 0.1000 (population of 3000 cells 0.100); the distance autocorrelation falls to 1/e in 20 blocks = **9.2 s** (frame-to-frame correlation at 10 s: 0.41).
  Consequence of the 1-min lag: the loop extruder in the saved simulations steps once per block, i.e. a loop of 200 kb processivity lives ~100 blocks = 46 s, against ~20 min quoted in the paper for CTCF loops (from earlier work). The polymer
  dynamics at this block time are therefore too fast for loop extrusion; the plain-polymer part of the E-P distance does not depend on it.

## The main consequence: the instantaneous curve is almost flat
The promoter responds to the E-P distance over minutes (a burst starts at rate k_on(d) and lasts ~45 s), but with these dynamics the distance changes within ~9 s. P(burst ON | distance NOW) is therefore nearly flat (contrast P(ON at 150 nm) / P(ON at 1000 nm) =
**1.36**, against 4.1 with the invented process whose visits lasted 240 s). What can be recovered is P(ON) for a distance that is HELD, estimated from the distance HISTORY (`ep_history.py`: P(ON at t) = a + b H_t with H_t the exponentially weighted
history of h(d), tau and d_half fitted, errors in variables integrated out). With true distances and true ON states of 3000 cells it finds tau = 45 s (= 1/k_off, the burst length: correct), d_half = 528 nm and a contrast of **4.19** (4.1 for the instantaneous
fit of the invented process; 6.2 for the same history fit of the invented process, close to the 6.7 expected for a distance held fixed).

## The three microscopes (`ep_acquisition.py`): a vector-level model, NOT a rendered movie
Errors per axis of the measured E-P difference vector (medians; baseline settings: 6000 photons/s per locus, 0.1 s exposure, 0.15 s between the two colours, plane spacing 0.8 um / biplane 0.6 um, focus scatter 0.3 um):

| error (nm per axis) | 2D widefield | z-stack (17 planes, 0.3 s per plane) | biplane (2 planes at once, half the photons each) |
|---|---|---|---|
| localization, in plane | 8 | 6 | 9 |
| localization, z | not measured | 27 | 51 |
| motion during the exposure / plane window | 38 | 62 | 38 |
| the two colours 0.15 s apart (sequential excitation) | 120 | 125 | 120 |
| piezo repeatability (30 nm per plane) | | 12 | |
| mean delay between frame start and the measurement | 0.07 s | 2.5 s | 0.07 s |
| both loci seen | 97% | 100% | 99% |
| total per axis (xy / z) | 126 | 140 / 144 | 126 / 136 |

The error is dominated by MOTION, not photons: at D = 0.002 um^2/s at 1 min and an MSD exponent 0.5 (extrapolated below one block, 0.46 s), a locus moves ~110 nm per axis in the 0.15 s between the two colours. Removing that delay (both colours at once, e.g. a second camera) is the single biggest improvement.

## Results: 12-20 populations per setting, polymer distances. contrast = P(ON at 150 nm) / P(ON at 1000 nm) for a distance held; reference 4.19 (3000 cells, true distance and true ON state)

| setting | 2D | z-stack | biplane | truth (true distance, true ON) |
|---|---|---|---|---|
| 600 cells x 1 h, contrast (mean +- sd over 20 populations) | 4.77 +- 1.00 | 4.25 +- 0.72 | 4.45 +- 0.80 | 3.80 +- 0.54 |
| width of the 95% interval | 5.1 | 3.9 | 3.1 | 2.3 |
| interval contains the reference | 95% | 95% | 100% | 100% |
| no coupling: contrast; interval wrongly excludes 1 in (60 populations each) | 1.10; 6.7% | 1.33; 3.3% | 1.23; 3.3% | |

- **2D vs 3D:** 3D narrows the interval by 25-40% and removes the projection bias; biplane is slightly better than the z-stack (interval 3.1 vs 3.9), but the standard deviations of the contrast (0.80 vs 0.72) cannot be separated with 20 populations.
- **Wrong picks (contrast at 0 / 5 / 10 / 20%; reference 4.19):** 2D 4.73 / 4.67 / 4.37 / 3.68; z-stack 4.31 / 4.38 / 3.88 / 3.10; biplane 4.49 / 3.88 / 3.63 / 3.26. At 20% the contrast falls by 22-28%.
- **Number of cells (contrast sd over populations; biplane / z-stack / truth):** 100: 2.5 / 2.1 / 1.1; 200: 1.4 / 1.5 / 1.0; 400: 0.89 / 0.80 / 0.71; 800: 0.69 / 0.64 / 0.53; 1600: 0.73 / 0.61 / 0.49; 3200: 0.57 / 0.52 / 0.42. Interval width at 3200 cells: 1.6 / 1.9 / 1.3.
  Roughly 10x more cells than with the invented process are needed for the same precision (invented: sd 0.77 at 50 cells, 0.055 at 3200, an easier problem because the distance stayed in contact for minutes).
- **Contact fraction 0.03 / 0.10 / 0.20 (reference 4.64 / 4.19 / 3.18):** biplane 5.44 / 4.49 / 3.25; z-stack 6.52 / 4.31 / 3.62; truth 5.63 / 3.80 / 3.01. Rare contact gives a noisy, biased-high estimate (interval contains the reference 67-75% of the time).
- **Frame interval 10 / 20 / 40 / 80 s (reference contrast 4.19 / 2.56 / 2.38 / 1.78):** the reference itself falls: the distance changes within 10 s, so sparse frames cannot follow the history. Biplane 4.49 / 3.38 / 2.19 / 1.68. With these dynamics, frames every ~10 s or faster are needed.
- **Delay between the colours 0 / 0.05 / 0.15 / 0.30 s (contrast sd; biplane):** 0.63 / 0.71 / 0.77 / 0.93; z-stack 0.67 / 0.74 / 0.80 / 1.78; 2D 1.18 / 1.14 / 0.99 / 1.43. Simultaneous colours remove ~20% of the scatter of the 3D modes.
- **Photons (x 6000 per locus per s; detection of both loci / z precision):** biplane 0.3: 92% / 124 nm; 1: 99% / 51 nm; 3: 100% / 26 nm; below 0.3 the biplane loses the loci (0% at 0.1). Z-stack 0.1: 100% / 154 nm; 1: 27 nm; 3: 14 nm. The contrast barely changes with photons (4.3-4.8) while
  detection holds: localization is not the limit.
- **Focus scatter (um) 0 / 0.3 / 0.6 / 1.0 / 1.5 (both loci seen):** 2D 100 / 97 / 78 / 55 / 39%; biplane 100 / 99 / 86 / 64 / 46%. A single plane or a biplane has to be re-focused on the pair; the z-stack does not (it spans the nucleus).
- **Which lag the D refers to (reference contrast at 10 / 30 / 60 / 120 s: 5.48 / 5.16 / 4.19 / 2.90; D for the pair at 60 s: 5.32):** the biology you can see depends strongly on this assumption.
- **Invented vs polymer, same measurement (57 nm noise in 3D, no acquisition model):** invented: instantaneous 4.11, history 6.50 (reference 6.19); polymer: instantaneous **1.36**, history 3.73 (reference 4.19).
- **ON state from the MS2 trace (300 cells, burst rule (3 frames, 3 starts) fixed in advance; median contrast over 4 populations, reference 4.19 at 10 s frames, 2.56 at 20 s):** biplane 10 s: 4.6-5.5; z-stack 10 s: 7.8-7.9; biplane 20 s: 3.4-3.5; z-stack
  20 s: 3.6-3.9; the MS2 noise (0.05-2.0 x one Pol II) changes little. The call's sensitivity is 0.4-0.5 (10 s) and its false-positive rate 0.004; the calls are biased HIGH because sensitivity rises where bursts cluster (near the enhancer) and the far level collapses
  (a = 0.006 against 0.041). The simple attenuation correction makes it worse (6-10). 4 populations per setting, ratio statistics heavy-tailed: treat as a warning, not a measurement.

## What changed relative to the invented-process results above
- The effect of the enhancer on P(ON | distance now) shrinks from contrast 4.1 to 1.4 once the distance fluctuates as fast as the polymer (with D = 0.002 um^2/s at 1 min) does; the effect is still there in the history of the distance (4.2).
- Many more cells are needed (precision at 800 cells is what the invented process gave at ~50-100), and frames must be fast.
- The wrong-pick sensitivity is similar (20% wrong: -22 to -28%).
- The microscope matters through MOTION: the 0.15 s colour delay costs 120 nm per axis. The choice between 2D, z-stack and biplane matters less than that and than the cell number.
- Null control still good; interval coverage fine for the true-ON analyses; the MS2-based ON estimate is the weak link.

## Update 2026-10-26: measurement timing fixed (results above not yet re-run)
`ep_acquisition.measure` took the polymer block at round(t / block_s), while the distance that drives the promoter uses floor (block k holds for
[k, k + 1) x block_s; `transcription.simulate_cell`): on average the measurement was 0.23 s later than it claimed, and the E-P vector moves ~200 nm per axis
per 0.46 s block. The measured-vs-true distance error exceeded the model's own sigma by 32 nm (2D) and 72 nm (biplane); with floor by 6 and 21 nm (150
cells, 20 min, one seed). The EIV kernel did not know that extra error, so the modality comparison above (in particular biplane, which suffered most) should
be re-run with `scripts/sweep_ep_polymer.py`. Still open: a z-stack measures ~2.5 s after the frame time (plane scan), so against the frame-time truth its
error is ~2x what the kernel assumes (physical, but whether the promoter state should be taken at the measurement time instead is a modelling choice).

## Placeholders (not measured)
Photons 6000 / s per locus, background 5 photons per pixel, z spread 2 um, focus scatter 0.3 um, 17 planes 0.8 um / 0.3 s per plane, biplane planes 0.6 um apart, piezo 30 nm, 0.15 s colour delay (Vutara preset: assumed), the gene (k_on 1/1200 s,
k_off 1/45 s, k_init 0.15 /s), coupling = k_on x (1 + 11 h(d)) with h a Hill curve (d_half 400 nm, n = 4); the MSD below one block is extrapolated with the simulated exponent; the motion of a locus within an exposure is Gaussian; the errors of the two loci are independent.
