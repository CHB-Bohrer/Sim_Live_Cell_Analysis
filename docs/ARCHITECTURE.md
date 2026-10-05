# Architecture and build plan

How the simulation-based validation pipeline fits together, what is built, and what comes next.
(Diagrams are [Mermaid](https://mermaid.js.org/); GitHub draws them automatically. Edit the text to change them.)

**Legend:** 🟩 built and tested · 🟨 stand-in (works, but a placeholder to be replaced) · ⬜ planned, not built

## 1. The pipeline

Every arrow is a set of files on disk (one folder per stage under `data/runs/<run_id>/`). A stage reads the previous
stage's files and writes its own, together with its resolved config, random seed and provenance.

```mermaid
flowchart TD
    CFG["Config YAML + random seed<br/>(configs/*.yaml, --set key=value)"]:::built

    CAL["Calibration<br/>polymer units to nm and s<br/>by matching locus MSD"]:::planned
    S1["Stage 1: chromatin dynamics<br/>NOW: toy loci, confined diffusion<br/>PLANNED: polychrom / OpenMM on the GPU"]:::standin
    S2["Stage 2: cells<br/>persistent IDs, motion, collisions,<br/>division, irregular deforming nuclei"]:::built
    S3["Stage 3: microscopy forward model<br/>PSF, noise, photon counts, bleaching<br/>nuclear channel + one colour per locus"]:::built

    GT[("Ground truth<br/>cells.csv, loci_truth.csv,<br/>true label masks")]:::built
    IMG[("Images<br/>nucleus.tif, locus0.tif, locus1.tif")]:::built

    S4A["Stage 4a: segmentation<br/>Cellpose (GPU) or watershed"]:::built
    S4B["Stage 4b: tracking, swappable adapter<br/>Trackastra (GPU)<br/>Ultrack and TrackMate planned"]:::built
    S5["Stage 5: isolate every cell<br/>one fixed-size movie per cell ID,<br/>processed in parallel"]:::built
    S6["Stage 6: per-cell locus analysis<br/>NOW: simple placeholder localizer<br/>LATER: your analysis"]:::standin
    S7A["Stage 7a: tracking accuracy<br/>Cell Tracking Challenge metrics<br/>+ identity switches"]:::built
    S7B["Stage 7b: error propagation<br/>how tracking errors change<br/>the locus results"]:::planned
    UI["Dashboard<br/>Movie, Metrics, Trajectories,<br/>Scans, Cells"]:::built

    CFG --> S1
    CFG --> S2
    CAL -.-> S1
    S2 -->|"toy loci need cell lifetimes;<br/>real polychrom runs independently"| S1
    S1 --> S3
    S2 --> S3
    S3 --> IMG
    S3 --> GT
    IMG --> S4A --> S4B
    S4B --> S5
    GT --> S5
    S5 --> S6
    S4B --> S7A
    GT --> S7A
    S6 --> S7B
    GT --> S7B
    S7A --> UI
    S7B -.-> UI
    S6 --> UI

    classDef built fill:#d9f2d9,stroke:#2e7d32,color:#111;
    classDef standin fill:#fff3c4,stroke:#b8860b,color:#111;
    classDef planned fill:#ececec,stroke:#888,stroke-dasharray: 4 3,color:#111;
```

## 2. How the validation works: two identity sources, one analysis

Stage 5 can isolate cells using either the **true** cell IDs or the **tracker's** IDs. The same stage 6 analysis runs
on both, so any difference between the two results is caused by tracking errors alone. Comparing both against the
simulation's true locus positions shows how accurate the analysis is.

```mermaid
flowchart LR
    MOVIE["Simulated movie<br/>(images + ground truth)"]:::built

    MOVIE --> TRUTHID["Cells identified by<br/>TRUE IDs"]:::built
    MOVIE --> TRACK["Segment + track<br/>(Cellpose + Trackastra)"]:::built
    TRACK --> TRACKID["Cells identified by<br/>TRACKED IDs"]:::built

    TRUTHID --> A1["Same stage 6 analysis"]:::standin
    TRACKID --> A2["Same stage 6 analysis"]:::standin

    A1 --> R1["Result with true IDs"]:::built
    A2 --> R2["Result with tracked IDs"]:::built

    R1 --> C1{{"Analysis accuracy<br/>estimate vs true loci"}}:::built
    MOVIE --> C1
    R1 --> C2{{"Effect of tracking errors<br/>tracked vs true IDs"}}:::planned
    R2 --> C2

    TRACK --> M["Tracking metrics<br/>TRA, DET, LNK, CHOTA,<br/>ID switches"]:::built
    MOVIE --> M

    classDef built fill:#d9f2d9,stroke:#2e7d32,color:#111;
    classDef standin fill:#fff3c4,stroke:#b8860b,color:#111;
    classDef planned fill:#ececec,stroke:#888,stroke-dasharray: 4 3,color:#111;
```

Tracking itself is also run twice: on **ground-truth masks** (isolates linking errors) and on **automatic
segmentation** (the realistic case).

## 3. Build roadmap

Boxes are in a sensible build order; arrows mean "needs". Items marked done are in the repository today.

```mermaid
flowchart TD
    D1["Environment: OpenMM CUDA, PyTorch CUDA,<br/>Trackastra, CuPy, Cellpose"]:::built
    D2["Cell motion, shapes, rendering,<br/>segmentation, tracking, scoring"]:::built
    D3["Sweep tool, dashboard,<br/>in-browser movie players"]:::built
    D4["Per-cell isolation +<br/>parallel per-cell analysis"]:::built

    N1["Real stage 1: polychrom chromatin on the GPU"]:::planned
    N2["Calibration: match locus MSD<br/>to your measured values"]:::planned
    N3["MS2 bursting, 3D z-stacks,<br/>optics presets: widefield, spinning disk, lattice"]:::planned
    N4["Your stage 6 analysis plugged in"]:::planned
    N5["Stage 7b: error propagation<br/>truth IDs vs tracked IDs"]:::planned
    N6["Ultrack and TrackMate adapters,<br/>compared on identical data"]:::planned
    N7["Scans: cell density, frame interval,<br/>tracker mode, many seeds"]:::planned
    N8["Check against a real, hand-tracked movie"]:::planned

    D1 --> D2 --> D3
    D2 --> D4
    N1 --> N2
    N2 --> N3
    D4 --> N5
    N4 --> N5
    N3 --> N5
    D3 --> N6
    D3 --> N7
    N5 --> N8
    N6 --> N8
    N7 --> N8

    classDef built fill:#d9f2d9,stroke:#2e7d32,color:#111;
    classDef planned fill:#ececec,stroke:#888,stroke-dasharray: 4 3,color:#111;
```

## 4. Files each stage reads and writes

| Stage | Folder in `data/runs/<run_id>/` | Main files |
|---|---|---|
| 1 chromatin (stand-in) | `stage1_chromatin/` | `loci_truth.csv` (true locus positions per cell and frame) |
| 2 cells | `stage2_cells/` | `cells.csv` (position, radius, shape, parent ID per cell and frame) |
| 3 microscopy | `stage3_microscopy/` | `nucleus.tif`, `locus0.tif`, `locus1.tif`, `labels.tif` (true cell IDs) |
| 4 segment + track | `stage4_segtrack/` | `auto_segmentation.tif`, `tracked_gt_masks.tif`, `tracked_auto_masks.tif`, `tracks_*.csv` |
| 5 isolate cells | `stage5_cells/<identity source>/` | `cells/cell_XXXX/{nucleus,locus0,locus1,mask}.tif`, `cell_summary.csv` |
| 6 analysis | `stage6_analysis/<identity source>/` | `loci_positions.csv`, `locus_distances.csv`, `scores.json` |
| 7 validation | `stage7_validation/` | `metrics.json`, `id_matches_*.csv`, `id_switch_events_*.csv` |

Every stage folder also holds `params.yaml` (the resolved config), `seed.txt` and `provenance.json` (git commit and
package versions), so any run can be reproduced from its config and seed. Identity sources: `truth`,
`tracked_gt_masks`, `tracked_auto_masks`.

## 5. Compute

| Work | Runs on |
|---|---|
| Chromatin simulation (OpenMM, planned) | GPU |
| Cellpose segmentation, Trackastra tracking | GPU (in separate processes) |
| Cell motion, image rendering, scoring | CPU |
| Per-cell isolation and analysis | CPU, one cell per worker process |
| Dashboard | CPU; the movie players run in your browser |

The single RTX 3090 is shared, so GPU stages run one after another, never at the same time.
