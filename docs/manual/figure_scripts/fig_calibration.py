"""Section 3 figure: blind calibration recovers the simulated instrument (CPU, a few minutes)."""
import numpy as np

from _style import plt, save, SIM, THEORY, GREY
from simlive.io.runs import REPO_ROOT, load_config
from simlive.stage3_microscopy.calibration.run import run_calibration, settings_for
from simlive.stage3_microscopy.render import entropy_for

FAST = {"dark": {"n_frames": 40, "exposures_s": [0.01, 1.0, 5.0], "extra_frames": 8},
        "flat": {"levels_photons": [0, 300, 1000, 3000, 8000, 20000], "n_frames": 10},
        "flicker": {"level_photons": 3000, "n_frames": 80},
        "beads": {"n_beads": 30, "photons": 30000, "n_frames": 6, "n_frames_3d": 4, "z_um": 1.0}}
cfg = load_config(REPO_ROOT / "configs" / "loci_demo.yaml", ["geometry.fov_um=[30, 30]"])
cal, table = run_calibration(cfg, entropy_for(3), {**settings_for(cfg), **FAST}, progress=lambda *_: None)
table.to_csv(__import__("_style").OUT / "calibration_table.csv", index=False)
t = table.dropna(subset=["estimate", "truth"]).copy()
t = t[np.isfinite(t.estimate.astype(float)) & np.isfinite(t.truth.astype(float))]
t["estimate"], t["truth"] = t.estimate.astype(float), t.truth.astype(float)
t = t[t.truth.abs() > 1e-9]
t["ratio"] = t.estimate / t.truth
g = t.groupby("quantity").ratio.agg(["min", "max", "median"]).sort_values("median")
fig, a = plt.subplots(figsize=(8, 0.28 * len(g) + 1.2))
y = np.arange(len(g))
a.hlines(y, g["min"], g["max"], color=SIM, lw=3, alpha=0.5)
a.plot(g["median"], y, "o", color=SIM)
a.axvline(1, color=THEORY, lw=1)
a.axvspan(0.95, 1.05, color=GREY, alpha=0.2, lw=0)
a.set_yticks(y, g.index, fontsize=8.5)
a.set(xlabel="estimate / true value (dot = median over colour channels, bar = range)", xlim=(0.3, 1.5),
      title="Blind calibration recovers the simulated microscope (short test acquisitions)")
fig.tight_layout()
save(fig, "calibration.png")
print(len(t), "quantities;", (abs(t.ratio - 1) < 0.05).mean())
