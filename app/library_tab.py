"""The dashboard's 🗂 Library tab: watch a chromatin library being generated, and choose which simulations a movie uses.

Reads only files written by scripts/run_chromatin_library.py (<library>/status.json, progress.log) and
scripts/run_chromatin.py (<sim>/progress.json, meta.json). Nothing here imports torch or OpenMM.
"""
import json
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

import style
from simlive.io.runs import DATA_ROOT
from simlive.stage1_chromatin import library as L

CHROM = DATA_ROOT / "chromatin"
_COLS = {"seed": st.column_config.NumberColumn(format="%d"), "size (kb)": st.column_config.NumberColumn(format="%d"),
         "mean loop (kb)": st.column_config.NumberColumn(format="%.0f"), "run time (min)": st.column_config.NumberColumn(format="%.1f")}
STALE_S = 180          # a "running" library that has written nothing for this long is probably not running any more


def _json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def libraries() -> list[str]:
    return sorted(p.name for p in CHROM.glob("*") if p.is_dir() and not p.name.startswith("_")) if CHROM.exists() else []


def _gpu() -> dict | None:
    """GPU load from nvidia-smi (None if it is not available)."""
    exe = shutil.which("nvidia-smi")
    if not exe:
        return None
    try:
        out = subprocess.run([exe, "--query-gpu=name,utilization.gpu,memory.used,memory.total,temperature.gpu",
                              "--format=csv,noheader,nounits"], capture_output=True, text=True, timeout=5).stdout
        n, u, mu, mt, t = [x.strip() for x in out.strip().splitlines()[0].split(",")]
        return {"name": n, "util": float(u), "mem_used": float(mu), "mem_total": float(mt), "temp": float(t)}
    except (OSError, ValueError, IndexError, subprocess.TimeoutExpired):
        return None


def _hm(hours: float | None) -> str:
    if hours is None:
        return "not known yet"
    return f"{hours * 60:.0f} min" if hours < 1 else f"{hours:.1f} h"


def library_state(lib_dir: Path) -> dict:
    """Everything the monitor shows, from the files on disk."""
    st = _json(lib_dir / "status.json")
    finished = len(list(lib_dir.glob("*_seed*/meta.json")))
    cur = {}
    if st.get("current_sim"):
        cur = _json(lib_dir / st["current_sim"] / "progress.json")
    last = max([st.get("updated", 0), cur.get("updated", 0)])
    if not st:
        kind, text = "idle", "No generator has been started for this library"
    elif st["state"] == "running" and time.time() - last > STALE_S:
        kind, text = "warn", "Stalled: no sign of life for %d min (window closed, or the PC slept?)" % ((time.time() - last) / 60)
    elif st["state"] == "running":
        kind, text = "run", "Running"
    elif st["state"] == "finished":
        kind, text = "good", "Finished"
    else:
        kind, text = "idle", "Stopped (run the same command again to resume)"
    return {"status": st, "finished": finished, "current": cur, "kind": kind, "text": text}


@st.fragment(run_every=5)
def _monitor() -> None:
    libs = libraries()
    if not libs:
        st.info("No chromatin libraries yet. Start one from a terminal (it can run for many hours):\n\n"
                "`scripts\\run.cmd python scripts\\run_chromatin_library.py --config configs\\chromatin\\library_loop_extrusion.yaml "
                "--library loop_extrusion --n 100`\n\nThis page then shows its progress live.")
        return
    lib = st.selectbox("Library", libs, key="mon_lib")
    d = CHROM / lib
    S = library_state(d)
    s, cur = S["status"], S["current"]
    requested = int(s.get("n_requested") or _json(d / "library.json").get("n_requested") or S["finished"])

    c0, c1 = st.columns([3, 2])
    with c0:
        style.badge(st, S["text"], S["kind"])
    c1.caption(f"updated {time.strftime('%H:%M:%S')} · refreshes every 5 s")

    a, b, c, e = st.columns(4)
    a.metric("Simulations finished", f"{S['finished']} / {requested}")
    if cur and cur.get("phase") not in (None, "done"):
        b.metric("Now running", s.get("current_sim", "?").replace(f"{s.get('sim_name', '')}_", ""),
                 help="the seed of the simulation currently on the GPU")
        c.metric("Speed (MD steps/s)", f"{cur.get('steps_per_s', 0):,.0f}")
    else:
        b.metric("Now running", "none")
        c.metric("Average time", f"{s['avg_minutes_per_sim']:.1f} min each" if s.get("avg_minutes_per_sim") else "not known yet")
    e.metric("Library done in", _hm(s.get("eta_hours")) if s.get("state") == "running" else "—",
             help="estimated from the average of the finished simulations")

    st.progress(min(S["finished"] / max(requested, 1), 1.0), text=f"Library: {S['finished']} of {requested} simulations")
    if cur and cur.get("total") and cur.get("phase") != "done":
        frac = cur["block"] / cur["total"]
        left = (cur["total"] - cur["block"]) / cur["blocks_per_s"] / 60 if cur.get("blocks_per_s") else None
        st.progress(min(frac, 1.0), text=f"This simulation: {cur['phase']}, block {cur['block']:,} of {cur['total']:,}"
                    + (f"  ·  about {left:.0f} min left" if left is not None else ""))

    g = _gpu()
    if g:
        st.caption(f"GPU · {g['name']} · {g['util']:.0f}% busy · {g['mem_used'] / 1024:.1f} / {g['mem_total'] / 1024:.0f} GB · "
                   f"{g['temp']:.0f} °C")
    if s.get("failed_seeds"):
        st.error(f"{len(s['failed_seeds'])} simulation(s) failed (seeds {s['failed_seeds']}). The reason is in the log below. "
                 "Re-running the same command retries them.")
    log = d / "progress.log"
    if log.exists():
        with st.expander("Recent log"):
            st.code("\n".join(log.read_text(errors="replace").splitlines()[-12:]), language="text")


def _how_it_works() -> None:
    st.markdown("### How simulations reach a movie")
    c = st.columns(4)
    with c[0]:
        style.card(st, "Generate once", "The GPU simulates the chromatin region many times (different random seeds). "
                   "Each result is a <b>saved simulation</b>; together they are a <b>library</b>. Takes hours, done once.", 1)
    with c[1]:
        style.card(st, "Choose", "For a movie, each cell is given <b>its own</b> saved simulation. By default the choice is a "
                   "random draw; below you can limit the pool or <b>pin</b> a cell to a particular simulation.", 2)
    with c[2]:
        style.card(st, "Place in the nucleus", "The two loci are read from the simulation at the genomic positions you set, "
                   "converted to nanometres, randomly rotated, and placed inside the nucleus so they move with it.", 3)
    with c[3]:
        style.card(st, "Image and analyse", "Frames are rendered as a microscope would see them, then segmented, tracked and "
                   "analysed. The true locus positions are kept, so the analysis can be checked.", 4)


@st.fragment
def _picker() -> None:
    libs = libraries()
    if not libs:
        return
    st.markdown("### Choose which simulations a movie uses")
    lib = st.selectbox("Library", libs, key="pick_lib")
    df = L.list_library(lib)
    if df.empty:
        st.info("This library has no finished simulations yet; they appear here as soon as the first one completes.")
        return
    df = df.sort_values("seed").reset_index(drop=True)
    mode = st.radio("How should cells get their simulation?", ["random", "pool", "pin"], horizontal=True, key="pick_mode",
                    format_func={"random": "🎲 Random from the whole library (default)",
                                 "pool": "☑ Random from the ones I tick",
                                 "pin": "📌 Pin particular cells to particular simulations"}.get)
    with st.expander("What does each option do?", expanded=False):
        st.markdown(
            "- **Random from the whole library:** every cell gets a *different* simulation, drawn at random. The draw "
            "is controlled by the run's random seed, so the same seed always gives the same assignment.\n"
            "- **Random from the ones I tick:** the same, but only from the simulations you tick below (for example "
            "only those with long loops, or a hand-picked set of ten).\n"
            "- **Pin:** you decide which simulation a given cell gets (cell 1 → simulation X). Cells you do not pin "
            "are still drawn at random from the remaining ones.\n\n"
            "Every cell gets its own simulation as long as there are enough; if a movie has more cells than available "
            "simulations the dashboard turns on *reuse* and tells you, because cells then share a trajectory.")

    show = df[["sim_id", "seed", "n_monomers", "mean_loop_kb", "minutes"]].rename(columns={
        "sim_id": "simulation", "n_monomers": "size (kb)", "mean_loop_kb": "mean loop (kb)", "minutes": "run time (min)"})
    pool_ids = None
    if mode == "pool":
        show.insert(0, "use", True)
        ed = st.data_editor(show, hide_index=True, width="stretch", key=f"pool_{lib}", disabled=show.columns[1:].tolist(),
                            column_config={**_COLS, "use": st.column_config.CheckboxColumn("use", help="tick = may be drawn")})
        pool_ids = ed.loc[ed["use"], "simulation"].tolist()
        st.caption(f"{len(pool_ids)} of {len(show)} simulations ticked.")
    else:
        st.dataframe(show, hide_index=True, width="stretch", column_config=_COLS)

    n_cells = int(st.number_input("Number of cells in the movie", 1, 200, 8, key="pick_ncells"))
    seed = int(st.number_input("Run seed (the movie's random seed)", 0, 10_000, 1, key="pick_seed",
                               help="Seed 1 is the default. The preview uses the same random stream the movie will use."))
    pinned = {}
    if mode == "pin":
        st.caption("Choose a simulation for the cells you want to fix; leave '(random)' for the rest.")
        cols = st.columns(min(n_cells, 4))
        for cid in range(1, n_cells + 1):
            ch = cols[(cid - 1) % len(cols)].selectbox(f"Cell {cid}", ["(random)"] + df.sim_id.tolist(), key=f"pin_{lib}_{cid}")
            if ch != "(random)":
                pinned[cid] = ch
        if len(set(pinned.values())) < len(pinned):
            st.warning("Two cells are pinned to the same simulation, so they will move identically.")

    avail = len(pool_ids) if pool_ids is not None else len(df)
    reuse = n_cells - len(pinned) > avail - len(set(pinned.values()) & set(pool_ids or df.sim_id))
    if reuse:
        st.warning(f"{n_cells} cells but only {avail} simulations are available: cells will **share** simulations "
                   "(reuse is switched on in the generated settings).")
    if pool_ids is not None and not pool_ids:
        st.error("Tick at least one simulation.")
        return
    try:
        rng = np.random.default_rng(np.random.SeedSequence(seed).spawn(4)[2])
        asg = L.assign_simulations(lib, range(1, n_cells + 1), rng, replace=reuse, pool=pool_ids, pinned=pinned)
    except ValueError as e:
        st.error(str(e))
        return
    st.markdown("**Preview: what each cell will get** (exactly what a movie with this seed does, if no cell divides)")
    prev = pd.DataFrame({"cell": list(asg), "simulation": list(asg.values()),
                         "how": ["pinned" if c in pinned else "random draw" for c in asg]})
    st.dataframe(prev, hide_index=True, width="stretch")

    ov = [f"loci.source=library", f"loci.library={lib}", f"loci.replace={'true' if reuse else 'false'}"]
    if pool_ids is not None:
        ov.append(f"loci.sim_ids=[{', '.join(pool_ids)}]")
    if pinned:
        ov.append("loci.pin={" + ", ".join(f"{k}: {v}" for k, v in pinned.items()) + "}")
    st.session_state["library_choice"] = {"library": lib, "pool": pool_ids, "pinned": pinned}   # read by the sidebar
    with st.expander("Use this in a movie"):
        st.markdown("In the sidebar, open **New single-cell loci run** and choose **Saved chromatin library, with my "
                    "choice from the Library tab**. Or from a terminal / config file:")
        st.code("scripts\\run.cmd python scripts\\run_tracking_demo.py --config configs\\loci_library.yaml " +
                " ".join(f'--set "{o}"' for o in ov) + f" --set seed={seed} --set cells.n_cells={n_cells}", language="bash")


def overrides_for_movie(n_cells: int) -> list[str] | None:
    """`--set` overrides for a movie with n_cells cells, from the choice made in the picker (None if nothing chosen)."""
    c = st.session_state.get("library_choice")
    if not c:
        return None
    n_avail = len(c["pool"]) if c["pool"] is not None else len(L.list_library(c["library"]))
    ov = ["loci.source=library", f"loci.library={c['library']}",
          f"loci.replace={'true' if n_cells > n_avail else 'false'}"]
    if c["pool"] is not None:
        ov.append(f"loci.sim_ids=[{', '.join(c['pool'])}]")
    if c["pinned"]:
        ov.append("loci.pin={" + ", ".join(f"{k}: {v}" for k, v in c["pinned"].items()) + "}")
    return ov


def render() -> None:
    _how_it_works()
    st.markdown("### Library generation: live status")
    _monitor()
    _picker()
