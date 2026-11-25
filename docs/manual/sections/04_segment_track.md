# Finding and following cells {#tracking}

::: summary
**In one sentence.** Nuclei are segmented with Cellpose and linked through time with Trackastra. The result is scored with the Cell Tracking Challenge metrics against the true IDs. Then each tracked cell is cut out as its own small movie for locus analysis.
:::

## What it does

- **Segmentation:** Cellpose by default. A simple threshold-and-watershed method is kept as a baseline.
- **Tracking:** Trackastra (`general_2d` model) sits behind a swappable "tracker" interface, so Ultrack or TrackMate can be added later and run on identical data.
- **Two runs per movie:** the tracker runs (a) on the *true* masks, which isolates linking errors, and (b) on Cellpose masks, which is the realistic case.
- **Per-cell movies:** for each cell ID, `isolate_cells` writes a fixed-size crop of every channel plus the cell's own mask. This is done for the true IDs and for the tracked IDs, so the same locus analysis can run on both. The difference between the two is what tracking errors cost (chapter 7).

## How we know it works

Scores come from the standard Cell Tracking Challenge (CTC) metrics (DET, TRA, LNK, CHOTA) and our own count of identity switches, over at least 5 seeds per setting.

| Cell speed (× default) | ≤ 0.5 | 0.7 | 0.8 | 1.0 | 1.25 |
|---|---|---|---|---|---|
| CHOTA (1 = perfect) | ≈ 1 | ≥ 0.99 | 0.975 | 0.91 | 0.85 |

The errors come from crowding combined with speed. When two nuclei touch, Trackastra often calls them a false division. With true masks at half speed, linking was perfect (0 ID switches in 15 runs). Cellpose was not affected by texture inside the nuclei, but the watershed baseline was (0.89 → 0.78).

## Limits

- **Not checked on real data:** Trackastra's pretrained model has not been checked on real data. One annotated real movie would be the real test.
- **Broken tracks:** tracks that break or swap are not yet stitched back together (the rest of stage 5).
- **2D only:** tracking is 2D, on projections of z-stacks.
