"""The result of a calibration (maps and numbers estimated from calibration images) and how to apply it to images and positions."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter


@dataclass
class Calibration:
    """Estimated properties of the imaging instrument, per channel, plus the bead-derived registration. Built by run.run_calibration."""
    channels: list
    shape: tuple
    run_exposure_s: float
    pixel_size_um: float
    ref: str = "nucleus"
    maps: dict = field(default_factory=dict)       # {channel: {"offset": (ny,nx), "flat": (ny,nx), "hot": (ny,nx) bool}}
    values: dict = field(default_factory=dict)     # {channel: {"gain_adu_per_e": ..., "read_noise_adu": ..., ...}}
    beads: dict = field(default_factory=dict)      # output of fit.fit_beads
    z_um: list = field(default_factory=list)

    # ------------------------------------------------------------------ registration
    def shift_px(self, channel: str) -> tuple[float, float]:
        s = (self.beads.get("shift_px") or {}).get(channel, [0.0, 0.0])
        return float(s[0]), float(s[1])

    def shift_z_um(self, channel: str) -> float:
        return float((self.beads.get("shift_z_um") or {}).get(channel, 0.0))

    def correct_positions(self, y_px, x_px, channel: str, z_um=None):
        """Remove the measured chromatic shift of `channel` from detected positions (pixels, and um along z)."""
        dy, dx = self.shift_px(channel)
        out = (np.asarray(y_px, float) - dy, np.asarray(x_px, float) - dx)
        return out if z_um is None else out + (np.asarray(z_um, float) - self.shift_z_um(channel),)

    # ------------------------------------------------------------------ images
    def correct(self, img: np.ndarray, channel: str, y0: int = 0, x0: int = 0, flat: bool = True, hot: bool = True) -> np.ndarray:
        """Dark-subtract, flat-field and hot-pixel-correct a (..., h, w) piece of an image whose top-left pixel is (y0, x0) in the full frame
        (parts outside the frame use the nearest edge pixel). Returns float32 counts: (counts - offset) / normalised flat."""
        m = self.maps[channel]
        ny, nx = self.shape
        h, w = img.shape[-2:]
        ix = np.ix_(np.clip(np.arange(y0, y0 + h), 0, ny - 1), np.clip(np.arange(x0, x0 + w), 0, nx - 1))
        out = img.astype(np.float32) - m["offset"][ix]
        if flat:
            out = out / np.maximum(m["flat"][ix], 0.2)
        if hot and m["hot"].any():
            hm = m["hot"][ix]
            if hm.any():
                size = (1,) * (out.ndim - 2) + (3, 3)
                out = np.where(hm, median_filter(out, size=size, mode="nearest"), out)
        return out.astype(np.float32)

    # ------------------------------------------------------------------ storage
    def save(self, directory) -> None:
        d = Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        arrays = {f"{c}__{k}": v for c, m in self.maps.items() for k, v in m.items()}
        np.savez_compressed(d / "calibration_maps.npz", **arrays)
        meta = {"channels": self.channels, "shape": list(self.shape), "run_exposure_s": self.run_exposure_s, "pixel_size_um": self.pixel_size_um,
                "ref": self.ref, "values": self.values, "beads": self.beads, "z_um": self.z_um}
        (d / "calibration.json").write_text(json.dumps(meta, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))

    @classmethod
    def load(cls, directory) -> "Calibration":
        d = Path(directory)
        meta = json.loads((d / "calibration.json").read_text())
        z = np.load(d / "calibration_maps.npz")
        maps: dict = {}
        for key in z.files:
            c, k = key.split("__")
            maps.setdefault(c, {})[k] = z[key]
        return cls(meta["channels"], tuple(meta["shape"]), meta["run_exposure_s"], meta["pixel_size_um"], meta["ref"], maps, meta["values"],
                   meta["beads"], meta.get("z_um", []))
