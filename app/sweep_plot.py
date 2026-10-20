"""Readable plots of parameter sweeps (📈 Scans and the 🔗 E–P curve scans).

A sweep varies several settings at once. Drawing one line for every combination of the settings other than the x axis gave a 16-entry legend on top of
the data. Here each setting gets ONE role: the x axis, the colour of the lines, the panel it is drawn in (small multiples sharing the y axis), or
"held at" one value. `default_roles` picks sensible roles (the most varied setting on x, then lines, then as many panels as stay readable); the dashboard
lets the user change them. The legend sits outside the plot.
"""
from __future__ import annotations

import itertools
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MAX_PANELS = 6
MAX_LINES = 8


def _order(values) -> list:
    return sorted(pd.unique(values), key=lambda v: (isinstance(v, str), v))


def default_roles(df: pd.DataFrame, settings: list[str], x: str | None = None) -> dict:
    """{"x": setting, "lines": [settings], "panels": [settings], "held": {setting: value}} for the sweep `df` (`x` given: that setting is the x axis)."""
    n = {c: df[c].nunique(dropna=False) for c in settings}
    numeric = {c: pd.api.types.is_numeric_dtype(df[c]) and not pd.api.types.is_bool_dtype(df[c]) for c in settings}
    ranked = sorted(settings, key=lambda c: (-n[c], settings.index(c)))
    # a numeric setting makes a real x axis; categories (methods, modes) are better as colours than as an x axis joined by lines
    x = x or sorted(settings, key=lambda c: (not (numeric[c] and n[c] > 1), -n[c], settings.index(c)))[0]
    rest = [c for c in ranked if c != x]
    lines = [c for c in rest if n[c] <= MAX_LINES][:1]
    rest = [c for c in rest if c not in lines]
    panels, count = [], 1
    for c in sorted(rest, key=lambda c: (n[c], settings.index(c))):      # fewest values first, while the grid stays small
        if count * n[c] <= MAX_PANELS:
            panels.append(c)
            count *= n[c]
    fixed = [c for c in rest if c not in panels]
    held = {}
    if fixed:                                  # a combination that exists in the sweep (sweeps are often not full grids): the most common one
        combo = df.groupby(fixed, dropna=False).size().idxmax()
        held = dict(zip(fixed, combo if isinstance(combo, tuple) else (combo,)))
    return {"x": x, "lines": lines, "panels": panels, "held": held}


def _label(cols, key) -> str:
    key = key if isinstance(key, tuple) else (key,)
    return ", ".join(f"{c}={v}" for c, v in zip(cols, key))


def draw(df: pd.DataFrame, xcol: str, ycol: str, ylabel: str, lines: list[str] | None = None, panels: list[str] | None = None,
         held: dict | None = None, hlines: list[tuple] | None = None, colours: list[str] | None = None) -> plt.Figure:
    """Mean ± sd over seeds of `ycol` against `xcol`, one coloured line per combination of `lines`, one panel per combination of `panels`, the settings in
    `held` fixed at the given values, individual seeds as faint dots. `hlines` = [(y, label, linestyle)] drawn in every panel. Legend outside, on the right."""
    lines, panels, held = lines or [], panels or [], held or {}
    d = df
    for c, v in held.items():
        d = d[d[c] == v]
    if d.empty:                                # the held values never occur together in this sweep
        fig, ax = plt.subplots(figsize=(8, 2.5))
        ax.axis("off")
        ax.text(0.5, 0.5, "No runs with these settings held at these values together:\nchange one of the 'held at' choices.", ha="center", va="center", fontsize=13)
        return fig
    colours = colours or plt.rcParams["axes.prop_cycle"].by_key()["color"]
    xs = _order(d[xcol])
    numeric = all(not isinstance(v, (str, bool, np.bool_)) for v in xs)
    xpos = {v: i for i, v in enumerate(xs)}
    panel_keys = list(itertools.product(*[_order(d[c]) for c in panels])) if panels else [()]
    line_keys = list(itertools.product(*[_order(d[c]) for c in lines])) if lines else [()]
    ncol = min(3, len(panel_keys))
    nrow = math.ceil(len(panel_keys) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.4 * ncol + 2.6, 3.9 * nrow + 0.4), sharey=True, squeeze=False)
    handles = {}
    for pi, pk in enumerate(panel_keys):
        ax = axes[pi // ncol][pi % ncol]
        dp = d
        for c, v in zip(panels, pk):
            dp = dp[dp[c] == v]
        for li, lk in enumerate(line_keys):
            g = dp
            for c, v in zip(lines, lk):
                g = g[g[c] == v]
            if g.empty:
                continue
            col = colours[li % len(colours)]
            m = g.groupby(xcol)[ycol].agg(["mean", "std"]).reindex(xs)
            # categories on x: points side by side, NOT joined by a line (there is no order between categories to draw)
            off = 0.0 if numeric else (li - (len(line_keys) - 1) / 2) * min(0.12, 0.6 / max(len(line_keys), 1))
            xv = list(xs) if numeric else [xpos[v] + off for v in xs]
            eb = ax.errorbar(xv, m["mean"], yerr=m["std"].fillna(0), marker="o", ms=6 if not numeric else 5, capsize=3, lw=1.8, color=col,
                             ls="-" if numeric else "none")
            ax.scatter([v if numeric else xpos[v] + off for v in g[xcol]], g[ycol], s=12, alpha=0.3, color=col, zorder=1)
            handles[_label(lines, lk) or ycol] = eb.lines[0]
        for y, lab, ls in hlines or []:
            h = ax.axhline(y, color="#e5e7eb", lw=1.1, ls=ls)
            handles[lab] = h
        ax.set_title(_label(panels, pk) or "", fontsize=13)
        if not numeric:
            ax.set_xticks(range(len(xs)))
            ax.set_xticklabels([str(v) for v in xs])
        elif xs[-1] > 0 and xs[0] > 0 and xs[-1] / xs[0] > 20:
            ax.set_xscale("log")
        ax.set_xlabel(xcol)
        if pi % ncol == 0:
            ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    for k in range(len(panel_keys), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    if len(handles) > 1:
        fig.legend(list(handles.values()), list(handles), loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=11, frameon=False,
                   title=", ".join(lines) if lines else None, alignment="left")
    if held:                                              # say which settings the plot is for (the y label already names the metric)
        fig.suptitle("held at " + ", ".join(f"{c}={v}" for c, v in held.items()), x=0.01, ha="left", fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


def controls(st, df: pd.DataFrame, settings: list[str], key: str) -> dict:
    """Streamlit widgets for the roles of the settings; returns the keyword arguments of `draw` (xcol, lines, panels, held)."""
    x = st.selectbox("X axis", settings, index=settings.index(default_roles(df, settings)["x"]), key=f"{key}_x")
    r = default_roles(df, settings, x)
    rest = [c for c in settings if c != x]
    if not rest:
        return {"xcol": x, "lines": [], "panels": [], "held": {}}
    c1, c2 = st.columns(2)
    lines = c1.multiselect("Colour lines by", rest, default=r["lines"], key=f"{key}_lines_{x}", help="One line per value of these settings.")
    left = [c for c in rest if c not in lines]
    panels = c2.multiselect("Separate panels by", left, default=[c for c in r["panels"] if c in left], key=f"{key}_panels_{x}",
                            help="One small plot per value of these settings (all share the y axis).")
    fixed = [c for c in left if c not in panels]
    held = {}
    if fixed:
        cols = st.columns(min(4, len(fixed)))
        for i, c in enumerate(fixed):
            vals = _order(df[c])
            dflt = r["held"].get(c, vals[0])
            held[c] = cols[i % len(cols)].selectbox(f"{c} held at", vals, index=vals.index(dflt) if dflt in vals else 0, key=f"{key}_held_{x}_{c}")
    n_p = int(np.prod([df[c].nunique() for c in panels])) if panels else 1
    if n_p > 12:
        st.warning(f"{n_p} panels is too many to read: hold one of the panel settings fixed instead.")
    return {"xcol": x, "lines": lines, "panels": panels, "held": held}
