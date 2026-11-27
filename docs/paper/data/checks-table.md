| Check | Compared against | Limit (set in advance) | Result | Run |
| --- | --- | --- | --- | --- |
| Probe binding model against exact theory | exact theory | mean within 3%, CV within 10%, at every release rate | **pass**: largest error in the mean number of attached probes: 0.2% | 2026-11-27, commit `7cc8baa` |
| Transcription model against the exact expectation | exact expectation | mean P(ON) within 0.02; mean MS2 signal within 10%; no time point off by more than 0.12 | **pass**: mean P(ON) off by 0.001; mean MS2 signal off by 0.2% (300 cells) | 2026-11-27, commit `7cc8baa` |
| Recovery of a known distance-to-activity curve through measurement error | the curve put in | fold within 20%, half-distance within 15%, and closer to the truth than the uncorrected fit | **pass**: fold 8.01 (true 8), half-distance 347 nm (true 350); contrast 7.07 corrected, 4.39 uncorrected, 7.04 true | 2026-11-27, commit `7cc8baa` |
