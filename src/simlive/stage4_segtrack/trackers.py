"""Stage 4b: swappable tracker adapters.

A tracker takes (images, instance masks) and returns a `TrackResult`: masks relabeled so the value is the
track ID, plus a CTC-style track table (label, t1, t2, parent). Ultrack / TrackMate adapters implement the same
`Tracker` protocol and register in TRACKERS, so they run on identical synthetic data.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
import pandas as pd


@dataclass
class TrackResult:
    masks: np.ndarray       # (T,Y,X) integer; value = track ID, 0 = background
    table: pd.DataFrame     # columns: label, t1, t2, parent (0 = no parent)


class Tracker(Protocol):
    name: str

    def track(self, images: np.ndarray, masks: np.ndarray) -> TrackResult: ...


class TrackastraTracker:
    name = "trackastra"

    def __init__(self, model: str = "general_2d", mode: str = "greedy", device: str | None = None, **_):
        self.model_name, self.mode, self.device = model, mode, device

    def track(self, images: np.ndarray, masks: np.ndarray) -> TrackResult:
        from trackastra.model import Trackastra
        from trackastra.tracking import graph_to_ctc

        model = Trackastra.from_pretrained(self.model_name, device=self.device or "automatic")
        graph, _ = model.track(images, masks, mode=self.mode)
        table, tracked = graph_to_ctc(graph, masks, check=False)
        return TrackResult(np.asarray(tracked).astype(np.uint16), table)


TRACKERS = {"trackastra": TrackastraTracker}


def make_tracker(name: str, **params) -> Tracker:
    if name not in TRACKERS:
        raise ValueError(f"Unknown tracker {name!r}; available: {sorted(TRACKERS)}")
    return TRACKERS[name](**params)
