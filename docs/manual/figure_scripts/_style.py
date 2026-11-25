"""Shared look for the manual's figures (light background, readable when printed)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / "figures"
SIM, THEORY, ACCENT, GREY = "#2563eb", "#111827", "#dc2626", "#9ca3af"

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False, "figure.dpi": 100, "savefig.dpi": 200})


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / name)
