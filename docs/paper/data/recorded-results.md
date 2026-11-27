| Result | Value | Evidence | Source |
| --- | --- | --- | --- |
| Polymer engine reproduces known behaviour | 11 of 12 checks; insulation 1.36-fold against a target of 1.5 (miss) | 5 conditions, 2 seeds | `docs/chromatin_validation.md` |
| Loop-extrusion library | 100 of 100 simulations finished, 2026-10-18 | run log | `docs/notes/chromatin.md` |
| Blind calibration against truth | 53 of 53 quantities within tolerance | one seed | `docs/CALIBRATION_AND_3D.md` |
| Tracking quality against cell speed | CHOTA 0.99 at 0.7x, 0.91 at 1x, 0.85 at 1.25x | 5 seeds per setting | `docs/notes/status_detail.md` |
| Spot finder, 2D, calibrated | 98.0% found, 0 wrong picks, 26 nm median | library run | `docs/SPOT_FINDING.md` |
| Two-copy linker, copy number right | 100% / 98% / 90% at up to 0.5 / 1-2 / 4 strays per nucleus per minute | 25 movies | `docs/SPOT_FINDING.md` |
| Sampler convergence, dense gene | simple sampler 3 of 12; tempered 12 of 12 | re-run 2026-11-09 | `docs/MCMC_CONVERGENCE.md` |
| Polymerase starts recovered | 79% within 30 s; 0.68 false per true | showcase run | `docs/notes/spot_finding.md` |
| E-P contrast by imaging mode | 4.78 (2D), 4.34 (z-stack), 4.39 (biplane); reference 4.19 | 600 cells x 1 h, re-run 2026-11-09 | `docs/EP_CURVE.md` |
| E-P null test | effect reported in 3 to 7% of populations (nominal 5%) | 60 populations | `docs/EP_CURVE.md` |
| Logic audit | six defects found and fixed; earlier conclusions unchanged after re-running | 2026-10-26, re-run 2026-11-09 | `docs/notes/audits.md` |
