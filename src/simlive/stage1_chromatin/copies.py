"""Several promoter-enhancer pairs per cell ("copies": the two homologous alleles of a diploid cell, or four after replication).

Config block `loci.copies` (absent = one copy per cell, exactly the old behaviour and the old random numbers):

    copies:
      n: [1, 2]                  # possible numbers of copies per cell ...
      weights: [0.5, 0.5]        # ... and the share of cells with each
      min_separation_um: 2.5     # territories of different copies at least this far apart in the image plane
      ms2: all                   # which copies carry the MS2 cassette: all | first (the other copies still transcribe, unlabelled)
      transcription: [{}, {k_on: 0.0016}]   # per-copy overrides of loci.transcription.params: separate bursting dynamics

Every copy is one promoter (locus 0, its colour) and one enhancer (locus 1, its colour), with its own chromatin trajectory, its own probe
binding and its own promoter / Pol II / MS2 process driven by ITS OWN promoter-enhancer distance. Copy 0 uses exactly the random numbers a
one-copy run uses; the extra copies draw from their own stream (`SeedSequence(seed).spawn(10)[8]` in the pipeline), so switching copies on
never changes copy 0.

`trace_id(cell_id, copy)` is one integer per (cell, copy) for tables that are keyed by a single id (MS2 traces, deconvolution): it equals
cell_id for copy 0, so one-copy runs keep their old keys.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

COPY_STRIDE = 100_000          # trace_id = cell_id + COPY_STRIDE * copy
MS2_CHOICES = ("all", "first")


@dataclass(frozen=True)
class CopySettings:
    n: tuple = (1,)
    weights: tuple = (1.0,)
    min_separation_um: float = 2.5
    ms2: str = "all"
    transcription: tuple = field(default_factory=tuple)      # per-copy dicts of transcription-parameter overrides

    def __post_init__(self):
        if len(self.n) != len(self.weights) or not self.n:
            raise ValueError("loci.copies: n and weights must be lists of the same (non-zero) length")
        if min(self.n) < 1 or max(self.n) > 8:
            raise ValueError("loci.copies.n: copies per cell must be between 1 and 8")
        if min(self.weights) < 0 or sum(self.weights) <= 0:
            raise ValueError("loci.copies.weights must be >= 0 and not all 0")
        if self.ms2 not in MS2_CHOICES:
            raise ValueError(f"loci.copies.ms2 must be one of {MS2_CHOICES}")

    @property
    def max_copies(self) -> int:
        return int(max(self.n))

    @property
    def active(self) -> bool:
        """True if any cell can have more than one copy (otherwise nothing changes)."""
        return any(k > 1 and w > 0 for k, w in zip(self.n, self.weights))

    def labelled(self, copy: int) -> bool:
        """Does this copy carry the MS2 cassette?"""
        return self.ms2 == "all" or int(copy) == 0

    def tx_overrides(self, copy: int) -> dict:
        return dict(self.transcription[copy]) if copy < len(self.transcription) and self.transcription[copy] else {}


def copy_settings(cfg: dict) -> CopySettings:
    c = (cfg.get("loci") or {}).get("copies") or {}
    n = c.get("n", [1])
    n = [int(n)] if np.isscalar(n) else [int(v) for v in n]
    w = c.get("weights") or [1.0] * len(n)
    w = [float(w)] if np.isscalar(w) else [float(v) for v in w]
    tx = c.get("transcription") or []
    return CopySettings(tuple(n), tuple(w), float(c.get("min_separation_um", 2.5)), str(c.get("ms2", "all")),
                        tuple(dict(d or {}) for d in tx))


def draw_copy_numbers(cell_ids, cfg: dict, rng: np.random.Generator | None) -> dict[int, int]:
    """Number of copies of every cell: {cell_id: n}. One copy everywhere (and no random draw) unless loci.copies allows more."""
    s = copy_settings(cfg)
    ids = [int(c) for c in cell_ids]
    if not s.active or rng is None:
        return {c: 1 for c in ids}
    p = np.asarray(s.weights, float) / sum(s.weights)
    k = rng.choice(np.asarray(s.n, int), size=len(ids), p=p)
    return {c: int(v) for c, v in zip(ids, k)}


def extra_centre(rng: np.random.Generator, others: list, min_sep: float, r_max: float, fits=None, tries: int = 2000) -> np.ndarray:
    """A territory centre (nucleus frame, nuclear radii, 3 components) for one more copy: inside r_max, at least `min_sep` (nuclear radii) from
    every centre in `others` in the image plane (x, y), and accepted by `fits(c)` (e.g. all of the copy's loci stay inside the nucleus).
    If no position satisfies the separation the best one found is returned (the copies then sit closer than asked)."""
    best, best_d = None, -1.0
    for _ in range(tries):
        c = rng.normal(size=3)
        c *= r_max * rng.random() ** (1 / 3) / np.linalg.norm(c)
        if fits is not None and not fits(c):
            continue
        d = min((float(np.hypot(*(c[:2] - np.asarray(o)[:2]))) for o in others), default=np.inf)
        if d >= min_sep:
            return c
        if d > best_d:
            best, best_d = c, d
    if best is None:
        raise ValueError("no territory position keeps this copy's loci inside the nucleus")
    return best


def trace_id(cell_id, copy) -> int:
    return int(cell_id) + COPY_STRIDE * int(copy)


def split_trace_id(tid) -> tuple[int, int]:
    tid = int(tid)
    return tid % COPY_STRIDE, tid // COPY_STRIDE


def trace_label(tid) -> str:
    c, k = split_trace_id(tid)
    return f"cell {c}" if k == 0 else f"cell {c}, copy {k}"


def with_copy(df, col: str = "copy"):
    """The table with a `copy` column (0 everywhere if it had none: runs made before copies existed)."""
    if df is None or col in df.columns:
        return df
    return df.assign(**{col: 0})


def keys(df, base=("cell_id", "t")) -> list[str]:
    """Join keys for pairing tables row by row: (cell_id, t) plus `copy` if the table has it."""
    return list(base) + (["copy"] if df is not None and "copy" in df.columns else [])
