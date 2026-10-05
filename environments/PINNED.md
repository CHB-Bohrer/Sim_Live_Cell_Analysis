# Pinned software versions (saved 2026-10-05)

Nothing in this project updates itself: conda and pip only change packages when a command asks them to. These files
let you rebuild the *exact* working setup even if something is upgraded by accident or a package disappears.

## What is saved

| What | File / place | How to use it |
|---|---|---|
| Every conda + pip package with version | `simlive-env-2026-10-05.yml` | `conda env create -n simlive-restored -f environments/simlive-env-2026-10-05.yml` |
| Exact conda packages (URLs + checksums) | `simlive-conda-explicit-win64-2026-10-05.txt` | `conda create -n simlive-restored --file environments/simlive-conda-explicit-win64-2026-10-05.txt` (then `pip install -r` the pip file) |
| Every pip package with version | `simlive-pip-freeze-2026-10-05.txt` | `pip install -r environments/simlive-pip-freeze-2026-10-05.txt` |
| A complete copy of the working environment | conda env `simlive-frozen-2026-10-05` (~8 GB, never modified) | `scripts\run.cmd` uses `simlive`; to use the frozen one: `conda activate simlive-frozen-2026-10-05` |
| polychrom (installed from GitHub, not PyPI) | commit `11a870cac8a3b168a2e0e11995d2e899c3ceb657` ("Fix the sign of getLinkingNumber (#79)", 2026-10-04) | `powershell -File scripts\install_polychrom.ps1` installs exactly this commit |
| Code state | git tag `env-2026-10-05` | `git checkout env-2026-10-05` |

## Key versions

| Component | Version |
|---|---|
| Python | 3.12.14 |
| PyTorch | 2.14.1+cu126 (CUDA 12.6) |
| OpenMM | 8.4 |
| NumPy / SciPy | 2.5.3 / 1.17.1 |
| Cellpose | 4.2.1.1 |
| Trackastra | 0.5.6 |
| GPU / driver | NVIDIA GeForce RTX 3090 / 595.95 |

## Model weights (downloaded once, not in git)

Pretrained weights are separate files that could change upstream, so their fingerprints are recorded.

| Model | Location | Size | SHA-256 |
|---|---|---|---|
| Trackastra `general_2d` | `%LOCALAPPDATA%\trackastra\trackastra\models\general_2d\model.pt` | 41 MB | `5183FC65D759B09DAD4A772FFB94C0598A8687DADF871050D41E4F5DE1C37B6B` |
| Cellpose `cpsam` | `%USERPROFILE%\.cellpose\models\cpsam` | 1.2 GB | `E1440429EB384F95AFE32BCBA6510F90D518EAEDC917EDE549BED6804004ABE2` |

If a hash ever differs, the weights changed and results may not be comparable with earlier ones.

## Rules to keep this true

- Do not run `pip install -U`, `conda update`, or reinstall packages in `simlive` without saving new lock files first
  (copy this folder's commands, with a new date) and noting the change in CLAUDE.md.
- Install anything new into a *separate* test environment first; only move it into `simlive` once results still match.
- Simulation results record the git commit and package versions (`provenance.json`), so old runs stay interpretable.
