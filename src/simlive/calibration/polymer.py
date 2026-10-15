"""Polymer units -> nm and seconds: the calibration that every library movie and the E-P curve study share.

A saved polychrom simulation has no length or time scale of its own (1 monomer = 1 kb, positions in polymer length units, one saved
"block" = 750 MD steps). Two numbers turn it into a physical trajectory (`calibrate`):
  - nm per polymer unit: chosen so that a pair `sep_kb` apart is within the working distance `working_nm` for `contact_fraction`
    of the time (user, 2026-10-14: 100 kb pair, 400 nm, 0.1), unless a measured value is given;
  - seconds per saved block: the single-locus MSD (centre-of-mass removed, 3D) at lag `lag_s` is matched to 6 D lag_s with
    D = `D_um2_s` (user, 2026-10-14: a SINGLE locus, lag 1 min, D = 0.002 um^2/s). The simulated MSD exponent is ~0.5, so D measured at
    another lag would give a very different block time.
Both rest on assumptions (they are stored in the result); replace them with measured numbers as soon as there are any.

`Calibration` is saved as JSON in `configs/calibration/` (`scripts/calibrate_chromatin.py` writes it from the saved validation simulations).
`scale_from_config` is what `library_loci` calls: `loci.calibration: default` reads `configs/calibration/default.json`; an explicit
`loci.nm_per_unit` / `loci.block_duration_s` in the config (or a `--set`) still overrides that one number, and with no `calibration` key the
config's own numbers are used (the old PLACEHOLDER behaviour).
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from simlive.io.runs import DATA_ROOT, REPO_ROOT

CALIBRATION_DIR = REPO_ROOT / "configs" / "calibration"
VAL_DIR = DATA_ROOT / "chromatin" / "validation"
SIMS = ("val_polymer_seed1", "val_polymer_seed2", "val_lef_free_seed1", "val_lef_free_seed2", "val_lef_leaky_seed1", "val_lef_leaky_seed2",
        "val_lef_ctcf_seed1", "val_lef_ctcf_seed2", "val_lef_loops_seed1", "val_lef_loops_seed2")
LAGS = np.array([1, 2, 5, 10, 20, 50, 100, 200, 400], float)        # blocks


@dataclass
class PolymerPool:
    rel: np.ndarray               # (sims, blocks, pairs, 3) float32: vector between the two monomers, polymer units
    msd_single: np.ndarray        # pooled single-locus MSD at LAGS (units^2, 3D, centre of mass of the tracked monomers removed)
    msd_rel: np.ndarray           # pooled MSD of the relative vector (the E-P separation vector)
    names: tuple
    sep_kb: float

    @property
    def n_blocks(self) -> int:
        return self.rel.shape[1]

    def distances(self) -> np.ndarray:
        return np.linalg.norm(self.rel, axis=-1)

    def autocorr_time_blocks(self) -> float:
        """Blocks until the autocorrelation of the E-P distance falls to 1/e (averaged over pairs and simulations)."""
        d = self.distances()
        z = d - d.mean(1, keepdims=True)
        ac = [float((z[:, :-l] * z[:, l:]).mean() / (z * z).mean()) for l in range(1, 400)]
        ac = np.array(ac)
        i = int(np.argmax(ac < 1 / np.e))
        return float(i + 1) if ac.min() < 1 / np.e else float("nan")


def load_pool(sep_kb: float = 100.0, sims: tuple = SIMS, root: Path = VAL_DIR) -> PolymerPool:
    rels, ms, mr = [], [], []
    for s in sims:
        d = root / s
        tp = np.load(d / "tracked_positions.npy").astype(np.float32)                 # (T, beads, 3)
        beads = np.load(d / "tracked_beads.npy")
        step = int(round(sep_kb / float(beads[1] - beads[0])))
        x = tp - tp.mean(1, keepdims=True)
        rel = (tp[:, step:] - tp[:, :-step]).astype(np.float32)
        rels.append(rel)
        dv = x[:, step:] - x[:, :-step]
        ms.append([float(((x[l:] - x[:-l]) ** 2).sum(2).mean()) for l in LAGS.astype(int)])
        mr.append([float(((dv[l:] - dv[:-l]) ** 2).sum(2).mean()) for l in LAGS.astype(int)])
    return PolymerPool(np.stack(rels), np.mean(ms, 0), np.mean(mr, 0), tuple(sims), sep_kb)


@dataclass(frozen=True)
class Calibration:
    nm_per_unit: float
    block_s: float
    alpha: float                  # exponent of the single-locus MSD below one block (extrapolation)
    msd1_nm2: float               # single-locus MSD at one block (nm^2, 3D)
    msd_nm2: tuple                # single-locus MSD at LAGS (nm^2)
    msd_rel_nm2: tuple
    lag_s: float
    D_um2_s: float
    applies_to: str
    working_nm: float = float("nan")        # what nm_per_unit was derived from (nan: given, not derived)
    contact_fraction: float = float("nan")
    sep_kb: float = float("nan")
    sims: tuple = ()                        # the saved simulations the MSD and the distances came from

    def msd_single_s(self, tau_s):
        """Single-locus 3D MSD (nm^2) after tau_s seconds: the simulated table above one block, the power law with the small-lag exponent below it (the simulation
        stores nothing between blocks)."""
        tau = np.asarray(tau_s, float)
        L = tau / self.block_s
        tab = np.interp(np.log(np.maximum(L, 1e-9)), np.log(LAGS), np.log(self.msd_nm2))
        small = self.msd1_nm2 * np.maximum(L, 0.0) ** self.alpha
        return np.where(L < 1.0, small, np.exp(np.where(L > LAGS[-1], np.log(self.msd_nm2[-1]), tab)))

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2)

    @classmethod
    def from_json(cls, text: str) -> "Calibration":
        d = json.loads(text)
        for k in ("msd_nm2", "msd_rel_nm2", "sims"):
            d[k] = tuple(d.get(k, ()))
        return cls(**d)

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_json() + "\n", encoding="utf-8")
        return path


def calibrate(pool: PolymerPool, working_nm: float = 400.0, contact_fraction: float = 0.10, D_um2_s: float = 0.002, lag_s: float = 60.0,
              applies_to: str = "locus", nm_per_unit: float | None = None, block_s: float | None = None) -> Calibration:
    d = pool.distances()
    nm = float(nm_per_unit) if nm_per_unit else float(working_nm / np.quantile(d, contact_fraction))
    msd = np.array(pool.msd_single if applies_to == "locus" else pool.msd_rel, float)
    alpha = float(np.polyfit(np.log(LAGS[:2]), np.log(msd[:2]), 1)[0])
    target = 6.0 * D_um2_s * lag_s * 1e6 / nm ** 2                                  # units^2
    if block_s is None:
        if target < msd[0]:
            L = (target / msd[0]) ** (1.0 / alpha)
        else:
            L = float(np.exp(np.interp(np.log(target), np.log(msd), np.log(LAGS))))
        block_s = lag_s / L
    return Calibration(nm, float(block_s), alpha, float(np.array(pool.msd_single)[0]) * nm ** 2, tuple(np.array(pool.msd_single) * nm ** 2),
                       tuple(np.array(pool.msd_rel) * nm ** 2), lag_s, D_um2_s, applies_to,
                       working_nm=float(working_nm), contact_fraction=float(contact_fraction), sep_kb=float(pool.sep_kb), sims=tuple(pool.names))


def calibration_path(name: str | Path) -> Path:
    """`default` -> configs/calibration/default.json; a path (absolute, or with a suffix) is used as given."""
    p = Path(name)
    return p if p.suffix or p.is_absolute() else CALIBRATION_DIR / f"{name}.json"


def load_calibration(name: str | Path = "default") -> Calibration:
    p = calibration_path(name)
    if not p.exists():
        raise FileNotFoundError(f"no polymer calibration {p}: make it with scripts\\run.cmd python scripts\\calibrate_chromatin.py")
    return Calibration.from_json(p.read_text(encoding="utf-8"))


def scale_from_config(lc: dict) -> tuple[float, float, str]:
    """(nm per polymer unit, seconds per saved block, where they came from) for the `loci:` block of a config.
    `loci.calibration` names a saved calibration; `loci.nm_per_unit` / `loci.block_duration_s`, if present, override that one number."""
    cal = lc.get("calibration")
    nm, blk = lc.get("nm_per_unit"), lc.get("block_duration_s")
    if not cal:
        if nm is None or blk is None:
            raise KeyError("loci needs `calibration: <name>` or both `nm_per_unit` and `block_duration_s`")
        return float(nm), float(blk), "config numbers (placeholders unless measured)"
    c = load_calibration(cal)
    src = f"calibration '{cal}'"
    if nm is not None:
        src += f", nm_per_unit overridden by {nm}"
    if blk is not None:
        src += f", block_duration_s overridden by {blk}"
    return float(nm if nm is not None else c.nm_per_unit), float(blk if blk is not None else c.block_s), src


def default_scale() -> tuple[float, float, str]:
    """What a new library movie should start from (dashboard defaults): configs/calibration/default.json, else the old placeholders."""
    try:
        return scale_from_config({"calibration": "default"})
    except FileNotFoundError:
        return 50.0, 10.0, "placeholders (no configs/calibration/default.json)"
