# How this project works

**In one sentence:** we invent a microscope movie where we know exactly where every cell and every locus is, run
our analysis on it as if it were real data, and check how close the analysis gets to the known answer.

Diagrams are [Mermaid](https://mermaid.js.org/); GitHub draws them automatically.

**Colour key:** 🟩 green = built and tested · 🟨 yellow = a simple stand-in that works but will be replaced · ⬜ dashed grey = planned, not built yet

---

## 1. The big idea in four steps

```mermaid
flowchart LR
    A["<b>1. SIMULATE</b><br/>Invent a microscope movie<br/>where we know exactly where<br/>every cell and locus is"]:::built
    B["<b>2. ANALYZE</b><br/>Treat it like real data:<br/>find the cells, follow them,<br/>locate the loci in each cell"]:::built
    C["<b>3. GRADE</b><br/>Compare what the analysis found<br/>with the known truth"]:::built
    D1["<b>4a. LEARN: image errors</b><br/>What does each error of the microscope<br/>do to the locus and MS2 results,<br/>and how much does correcting it give back?"]:::built
    D2["<b>4b. LEARN: tracking</b><br/>How much do mistakes in following<br/>the cells change the locus results?"]:::planned
    A --> B --> C --> D1
    C --> D2

    classDef built fill:#d9f2d9,stroke:#2e7d32,color:#111;
    classDef planned fill:#ececec,stroke:#888,stroke-dasharray: 4 3,color:#111;
```

Why simulate? In a real movie nobody knows the right answer. In a simulated one we do, so we can measure
exactly how much each part of the analysis gets wrong.

---

## 2. What happens, step by step

Three zones. **Top:** we build the fake experiment (so we keep an "answer key"). **Middle:** we analyze the movie
while ignoring the answer key. **Bottom:** we grade the analysis against the answer key.

```mermaid
flowchart TD
    subgraph SIM["① SIMULATION: we invent the experiment, so we know the right answer"]
        direction TB
        S2["<b>Cells</b> (stage 2)<br/>Cells drift, bump into each other,<br/>change shape and sometimes divide.<br/>Each cell keeps its own ID number."]:::built
        LIB[("<b>CHROMATIN LIBRARY</b> (stage 1, made separately)<br/>many saved polymer simulations of a<br/>chromatin region (polychrom, GPU):<br/>confined polymer + loop extrusion<br/>with boundary elements")]:::standin
        S1["<b>Loci in each cell</b> (stage 1)<br/>One or more promoter + enhancer pairs per cell<br/>(copies, e.g. both alleles): toy motion, or one saved<br/>polymer simulation per copy (100-simulation library).<br/>nm / s come from the calibration above (assumed)."]:::standin
        S1B["<b>Probes and transcription</b> (stage 1b, 1c)<br/>Probe binding sets each locus' brightness;<br/>every copy's promoter bursts on its own, driven<br/>by its enhancer distance, Pol II carrying MS2 loops"]:::built
        CAL["<b>Calibration of the motion</b> (configs/calibration)<br/>nm and seconds per polymer unit from ASSUMED<br/>contact fraction and locus D;<br/>replace with your measured MSD"]:::standin
        S3["<b>Microscope</b> (stage 3)<br/>Turns cells and loci into recorded images (2D or 3D z-stacks)<br/>with 28 switchable image-error sources (drift, chromatic shift,<br/>blur, haze, per-cell bleaching, uneven light, cross-talk,<br/>stray decoy spots, camera noise...). Cells can walk in and out of the field.<br/>Nucleus + MS2 in one colour, each locus in its own."]:::built
        ICAL["<b>Calibration images</b> (stage 3)<br/>Dark frames, flat-field slide, bead slide<br/>of the same instrument"]:::built
        TRUTH[("<b>ANSWER KEY</b><br/>true cell outlines and ID numbers<br/>true position of every locus")]:::built
        MOVIE[("<b>THE MOVIE</b><br/>nucleus images +<br/>one image per locus colour")]:::built
    end

    subgraph ANA["② ANALYSIS: pretend we do NOT know the answer"]
        direction TB
        S4A["<b>Find the cells</b> (stage 4a)<br/>Draw an outline around every nucleus<br/>in every frame (Cellpose, GPU)"]:::built
        S4B["<b>Follow each cell through time</b> (stage 4b)<br/>Decide which outline in frame 2<br/>is the same cell as in frame 1<br/>(Trackastra, GPU)"]:::built
        S5["<b>Cut out each cell</b> (stage 5)<br/>Make one small movie per cell,<br/>centred on it, so cells can be<br/>analyzed side by side in parallel"]:::built
        SC["<b>Correct the images</b> (stage 6)<br/>with the calibration images, or<br/>from the movie alone: background and<br/>hot pixels, illumination, colour shift,<br/>per-cell bleaching, camera rows"]:::built
        S6["<b>Locate the loci</b> (stage 6)<br/>Chen et al. 2025 spot finder: difference of<br/>Gaussians + Viterbi paths; several copies per<br/>cell paired promoter-enhancer, stray spots rejected.<br/>MS2 trace read at each tracked promoter.<br/>A placeholder for YOUR analysis."]:::standin
        S6M["<b>MS2 → Pol II</b> (stage 6)<br/>Bayesian deconvolution of each trace<br/>into Pol II start times and bursts<br/>(Chen et al. 2025)"]:::built
    end

    subgraph GRADE["③ GRADING: compare with the answer key"]
        direction TB
        S7A["<b>Were the right cells followed?</b> (stage 7a)<br/>Standard tracking scores<br/>plus a count of identity switches"]:::built
        S7E["<b>What does each image error do?</b> (stage 7)<br/>Re-render the movie with one error at a time,<br/>analyse with and without its correction"]:::built
        S7B["<b>Did tracking mistakes change the loci results?</b> (stage 7b)<br/>Locus results scored through<br/>the tracker's IDs vs the true IDs"]:::built
    end

    CAL -.-> S1
    LIB -.->|"picked per cell"| S1
    S2 --> S1
    S1 --> S1B
    S1 --> S3
    S1B --> S3
    S2 --> S3
    S3 --> MOVIE
    S3 --> ICAL
    S2 --> TRUTH
    S1 --> TRUTH
    MOVIE --> S4A --> S4B --> S5 --> SC --> S6 --> S6M
    ICAL -.-> SC
    S4B --> S7A
    S6 --> S7B
    S6 --> S7E
    S6M --> S7E
    TRUTH ==> S7A
    TRUTH ==> S7E
    TRUTH ==> S7B

    classDef built fill:#d9f2d9,stroke:#2e7d32,color:#111;
    classDef standin fill:#fff3c4,stroke:#b8860b,color:#111;
    classDef planned fill:#ececec,stroke:#888,stroke-dasharray: 4 3,color:#111;
```

Thick arrows carry the answer key. The analysis zone never sees it; only the grading zone does.
(The chromatin library is simulated on its own, ahead of time, and checked against published Mirny-lab results
([validation report](chromatin_validation.md)). Each cell can already be given its own saved simulation; that path is
tested with a small fake library but has not yet been run with a real 100-simulation library. The loci need to know when
each cell exists, which is why stage 2 feeds stage 1.)

---

## 3. Why we care: what a tracking mistake does to a locus measurement

The tracker has to decide which cell in frame 2 is "the same cell" as in frame 1. When two cells touch, it can get
this wrong. Because every locus is analyzed *inside its cell's own movie*, a wrong decision puts the wrong cell's
locus into that movie.

```mermaid
flowchart LR
    subgraph OK["When the tracker gets it right"]
        direction TB
        a1["Frame 1<br/>Cell A, locus at position P1"] --> a2["Frame 2<br/>Cell A, locus at P2"] --> a3["Frame 3<br/>Cell A, locus at P3"]
        a3 --> ar["<b>Smooth locus path</b><br/>correct motion, correct<br/>distance between loci"]:::good
    end
    subgraph BAD["When the tracker swaps two touching cells"]
        direction TB
        b1["Frame 1<br/>Cell A, locus at P1"] --> b2["Frame 2<br/>'Cell A' label has moved to<br/>cell B: locus at Q2"] --> b3["Frame 3<br/>still cell B: locus at Q3"]
        b3 --> br["<b>Path jumps between two cells</b><br/>false motion, wrong distances"]:::bad
    end

    classDef good fill:#d9f2d9,stroke:#2e7d32,color:#111;
    classDef bad fill:#fbd5d5,stroke:#c62828,color:#111;
```

Stage 7b measures exactly this. We run the same locus analysis twice, once with the true cell IDs and once with the
tracker's IDs, and see how far the results drift apart.

```mermaid
flowchart LR
    MV["Simulated movie"]:::built
    MV --> T["Cells identified by<br/><b>TRUE IDs</b><br/>(from the answer key)"]:::built
    MV --> K["Find + follow the cells<br/>(Cellpose + Trackastra)"]:::built
    K --> TR["Cells identified by<br/><b>TRACKER IDs</b>"]:::built
    T --> A1["Locus analysis"]:::standin
    TR --> A2["Same locus analysis"]:::standin
    A1 --> R1["Result using true IDs"]:::built
    A2 --> R2["Result using tracker IDs"]:::built
    R1 --> CMP{{"<b>Difference = the damage<br/>caused by tracking mistakes</b><br/>(scores.json: propagation;<br/>table on the 📊 Metrics page)"}}:::built
    R2 --> CMP
    R1 --> ACC{{"<b>How accurate is the analysis itself?</b><br/>estimated vs true locus positions"}}:::built

    classDef built fill:#d9f2d9,stroke:#2e7d32,color:#111;
    classDef standin fill:#fff3c4,stroke:#b8860b,color:#111;
    classDef planned fill:#ececec,stroke:#888,stroke-dasharray: 4 3,color:#111;
```

The tracker is also tested twice: once on perfect cell outlines (this isolates its linking mistakes) and once on the
outlines Cellpose finds (the realistic case).

---

## 4. How each source of image error is handled

Every one of the 24 image-error sources is handled in one of four ways. The dashboard's 🖼 Images page (part 2) has a list to pick a
source: it explains it, shows the same frame without it, with it and corrected, and shows what it does to the results with and
without each correction. Those results are measured, not assumed: the movie is rendered again with only that error and analysed
both ways (`scripts/error_effects.py`).

```mermaid
flowchart LR
    ERRS["<b>An image-error source</b><br/>(one of 24)"]:::built
    ERRS --> FUND["<b>Fundamental: cannot be removed</b><br/>blur (PSF), photon shot noise,<br/>read noise, quantum efficiency.<br/>They set the best possible result."]:::built
    ERRS --> AN["<b>Handled by the analysis itself</b><br/>median filter (hot pixels, cosmic rays),<br/>shell background around each spot (background,<br/>haze, autofluorescence, cross-talk, offsets),<br/>Viterbi path (cosmic rays), per-cell crops (drift)"]:::built
    ERRS --> IMG["<b>Calibration images</b><br/>dark frames (offsets, hot pixels),<br/>flat-field slide (illumination, pixel gains),<br/>bead slide (colour shift, drift)"]:::built
    ERRS --> MOV["<b>From the movie alone</b> (no calibration images)<br/>background + hot-pixel map from frames free of nuclei,<br/>illumination from moving nuclei, colour shift from the<br/>mean enhancer-promoter vector, cross-talk<br/>unmixing, camera rows, per-cell bleaching<br/>from each cell's light dose"]:::built
    ERRS --> NOT["<b>Not corrected yet</b><br/>defocus in 2D, flicker of the locus colours,<br/>slice skew (3D)"]:::planned
    AN --> RES{{"<b>Measured effect on the results</b><br/>wrong picks, position error, distance error,<br/>intensity error and drift, MS2 trace"}}:::built
    IMG --> RES
    MOV --> RES

    classDef built fill:#d9f2d9,stroke:#2e7d32,color:#111;
    classDef planned fill:#ececec,stroke:#888,stroke-dasharray: 4 3,color:#111;
```

The movie-only estimates are checked against the truth and guarded against doing harm. The illumination map is used only if it
predicts held-out cells better than even light. The colour shift is applied only if it is larger than twice its standard error,
which needs many cells.

---

## 5. From the MS2 spot to the moments Pol II starts

```mermaid
flowchart LR
    SPOT["<b>MS2 spot intensity</b><br/>at the tracked promoter,<br/>every frame"]:::built
    CORR["<b>Corrected trace</b><br/>background removed,<br/>bleaching + illumination<br/>divided out"]:::built
    UNIT["<b>In units of one Pol II</b><br/>divided by the brightness<br/>of one fully loaded Pol II"]:::built
    KER["<b>What one Pol II looks like</b><br/>signal steps up while it copies<br/>the MS2 loops, stays while it<br/>travels, drops when released"]:::built
    NOISE["<b>Noise model</b><br/>from frame-to-frame changes"]:::built
    MCMC["<b>Bayesian deconvolution</b><br/>1,000 possible histories of<br/>Pol II starts that fit the trace"]:::built
    OUT["<b>Pol II starts per minute</b><br/>with uncertainty, and the<br/>probability of a burst"]:::built
    SPOT --> CORR --> UNIT --> MCMC
    KER --> MCMC
    NOISE --> MCMC
    MCMC --> OUT
    OUT -.->|"graded against the true starts<br/>and promoter states"| G["Answer key"]:::built

    classDef built fill:#d9f2d9,stroke:#2e7d32,color:#111;
```

The method is Chen, Levo, Zoller & Gregor 2025 (Nat. Struct. Mol. Biol.): same model, prior and burst rule. The sampler differs so
that it also works for genes that start Pol II only every few minutes. Dashboard page: 🧮 MS2 → Pol II.

---

## 6. What is built and what comes next

```mermaid
flowchart LR
    subgraph P1["Done"]
        direction TB
        d1["Simulated cells that move, deform,<br/>divide and walk in and out of view"]:::built
        d2["Microscope images with 28 error<br/>sources (incl. stray spots), 2D or 3D z-stacks"]:::built
        d3["Cellpose + Trackastra,<br/>scored against the truth"]:::built
        d4["One movie per cell, analyzed<br/>in parallel"]:::built
        d5["Chromatin engine (polychrom), validated;<br/>100-simulation library made"]:::built
        d6["Probe binding, distance-driven<br/>transcription and MS2"]:::built
        d7["Spot finder (Chen 2025), MS2 trace,<br/>MS2 → Pol II deconvolution"]:::built
        d8["Calibration from images AND<br/>from the movie; measured effect<br/>of every image error"]:::built
        d9["Dashboard (dark theme), players,<br/>sweeps, narrated tutorial linked"]:::built
        d10["Instrument preset: Bruker Vutara VXL (widefield mode),<br/>bleed-through from spectra"]:::built
        d11["Several promoter-enhancer copies per cell;<br/>multi-copy linker, copies paired"]:::built
    end
    subgraph P2["Next: make it real (needs your data)"]
        direction TB
        n2["Calibrate to your measured<br/>MSD (nm and seconds): the file exists,<br/>its inputs are still assumed"]:::standin
        n3["Measured values: MS2 brightness R,<br/>error sizes, real movie"]:::planned
        n4["Use MS2 as a third witness;<br/>sister chromatids"]:::planned
        n5["Vutara biplane 3D (two planes<br/>600 nm apart, 50:50 split)"]:::built
    end
    subgraph P3["Then: the real question"]
        direction TB
        m1["Plug in YOUR analysis"]:::planned
        m2["Measure how tracking mistakes<br/>change the locus results (7b): built,<br/>stitching broken tracks is not"]:::built
    end
    subgraph P4["Finally: trust it"]
        direction TB
        f1["Compare Ultrack and TrackMate<br/>on the same movies"]:::planned
        f2["Scan cell density and<br/>frame interval, many seeds"]:::planned
        f3["Check against a real,<br/>hand-tracked movie"]:::planned
    end
    P1 --> P2 --> P3 --> P4

    classDef built fill:#d9f2d9,stroke:#2e7d32,color:#111;
    classDef planned fill:#ececec,stroke:#888,stroke-dasharray: 4 3,color:#111;
```

---

## 7. Glossary

| Term | Plain meaning |
|---|---|
| **Locus** | A tagged spot on the DNA that shows up as a bright dot in the image. |
| **Segmentation** | Drawing an outline around every nucleus in an image. |
| **Tracking** | Deciding which outline in the next frame is the same cell. |
| **Identity switch** | The tracker gives a cell a different ID number than it had a frame ago (a swap, or a track that breaks into two). |
| **Ground truth / answer key** | What the simulation knows to be true: real cell outlines, IDs and locus positions. |
| **CHOTA, TRA, DET, LNK** | Standard tracking scores. 1.0 is perfect. CHOTA is the strictest of them. |
| **MSD** | Mean squared displacement: how far a locus typically moves over a given time. Used to calibrate the simulation. |
| **Stand-in** | A simple placeholder that makes the pipeline run end to end until the real piece is built. |

---

## 8. Files each stage reads and writes

Every stage writes its outputs into `data/runs/<run_id>/`, one folder per stage, so any stage can be re-run or
inspected on its own. Each folder also stores the settings used (`params.yaml`), the random seed (`seed.txt`) and the
software versions (`provenance.json`), so a run can always be reproduced.

| Stage | Folder | Main files |
|---|---|---|
| 1 loci (stand-in or saved polymer library) | `stage1_chromatin/` | `loci_truth.csv`: true locus positions |
| 1b probe binding | `stage1_chromatin/` | `probe_occupancy.csv`: probes attached per locus per frame |
| 1c transcription (MS2) | `stage1_chromatin/` | `transcription.csv` (distance, promoter state, Pol II, MS2 loops), `transcription_events.csv` (every Pol II) |
| 2 cells | `stage2_cells/` | `cells.csv`: position, size, shape and ID of every cell in every frame |
| 3 microscope | `stage3_microscopy/` | `nucleus.tif`, `locus0.tif`, `locus1.tif`, `labels.tif` (true cell IDs), `imaging_truth.csv`, `cells_image.csv`, `loci_truth_image.csv`, `bleach_cells.csv` (per-cell bleaching), `imaging_errors.json` (see docs/IMAGE_ERRORS.md), `calibration/` (estimated from simulated dark / flat / bead images) |
| 4 find + follow | `stage4_segtrack/` | `tracked_gt_masks.tif`, `tracked_auto_masks.tif`, `tracks_*.csv` |
| 5 cut out cells | `stage5_cells/<who identified them>/` | one folder per cell, plus `cell_summary.csv` |
| 6 locate loci | `stage6_analysis/<who identified them>/` | `loci_positions[_cal].csv`, `locus_distances.csv`, `ms2_traces[_cal].csv`, `ms2_deconvolution.csv` (+ samples, info), `nuclear_signals.csv`, `flatfield_est.npz`; `stage6_analysis/scores.json`, `bleach.csv` |
| 7 grading | `stage7_validation/` | `metrics.json`, `id_switch_events_*.csv`, `error_effects.csv` + `error_effects_selfcal/` (what each image error does, with and without correction) |

"Who identified them" is `truth` (the answer key's IDs), `tracked_gt_masks` or `tracked_auto_masks` (the tracker's IDs).

## 9. What runs where

| Work | Runs on |
|---|---|
| Polymer simulation (chromatin library) | GPU (OpenMM) |
| Finding cells, following cells | GPU (Cellpose, Trackastra), each in its own process |
| Cell motion, image rendering, scoring, MS2 deconvolution (numba) | CPU |
| Analyzing each cell | CPU, one cell per worker process, all cells at once |
| Dashboard | CPU; the movie players run in your browser |

The one RTX 3090 is shared, so GPU steps run one after another, never at the same time.
