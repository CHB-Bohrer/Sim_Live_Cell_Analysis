<!--
Section: Segmentation and tracking (stages 4 and 7a)
Update when: the segmenter, tracker adapter, tracking metrics or the evaluation protocol change.
Code: src/simlive/stage4_segtrack, src/simlive/stage7_validation/tracking_metrics.py, src/simlive/stage7_validation/identity.py, src/simlive/pipelines
Check: none (GPU: tests/test_segtrack_gpu.py)
Sources: docs/notes/status_detail.md "Findings so far".
-->

# Finding and following nuclei

::: {.box}
**Short version.** Cellpose finds nuclei and Trackastra links them through time. Tracking is scored twice, once on perfect masks and once on automatic ones, so linking errors and segmentation errors can be told apart. Errors come from crowding combined with speed.
:::

Segmentation uses Cellpose by default, with a threshold-and-watershed method as a simple baseline. Tracking uses Trackastra's pretrained model behind a small `Tracker` interface, so other trackers can be compared on identical data.

**Protocol.** The tracker runs on (a) the true masks, which isolates linking error, and (b) the automatic masks, which is the realistic case. Both are scored against true cell identities with the Cell Tracking Challenge metrics and with the project's own count of identity switches.

Recorded findings (simulated 2D nuclei, five seeds per setting):

| Cell speed (`speed_scale`) | Tracking quality (CHOTA) |
| --- | --- |
| up to 0.7 | at least 0.99 |
| 0.8 | 0.975 |
| 1.0 (default) | 0.91 |
| 1.25 | 0.85 |

- At half speed or slower, linking on perfect masks made no identity switch in 15 runs. About one switch per run remains from segmentation.
- Touching pairs are often called a false division.
- Texture inside nuclei did not help linking, and it broke the watershed baseline (0.89 to 0.78) while leaving Cellpose unaffected.
- Adding the full image-error model changed tracking little (0.909 with errors, 0.903 without, three seeds): bright nuclear images are robust.

**[Unverified]** The pretrained tracker has not been tested on a real annotated movie. That is the proper test, and it needs data from the owner.

Practical constraint: Cellpose runs in its own worker process, and the tracker must load before any heavy linear algebra, because of conflicting OpenMP libraries on Windows.
