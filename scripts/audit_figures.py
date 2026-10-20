"""Render dashboard pages headlessly (Streamlit AppTest) and save every matplotlib figure they draw.

    scripts\\run.cmd python scripts\\audit_figures.py <run name> <out dir> [page keys...]

Used to audit every dashboard figure (2026-10-19): look at each saved PNG and ask whether it answers its question. Pages without matplotlib
figures (Movie, Metrics, Library, Config) save nothing; sections behind toggles or buttons are not rendered.
Every figure is saved when the page closes it (plt.close is hooked), named <page>/<nn>_<function>_<line>.png, so each one can be traced
to the code that drew it. Exceptions on the page are written to <out>/<page>/errors.txt.
"""
import inspect
import sys
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
run_name, out = sys.argv[1], Path(sys.argv[2])
pages = sys.argv[3:] or ["movie", "metrics", "traj", "scan", "cells", "spots", "img", "cal", "stage", "tx", "ms2", "ep", "probe", "lib", "chrom", "cfg"]
STATE = {"page": None, "n": 0}
_close = plt.close
seen = set()


def close(fig=None):
    figs = [fig] if fig is not None and not isinstance(fig, str) else ([plt.gcf()] if fig is None else plt.get_fignums() and [plt.figure(n) for n in plt.get_fignums()])
    for f in figs or []:
        if f is None or id(f) in seen or not hasattr(f, "savefig"):
            continue
        seen.add(id(f))
        caller = next((fr for fr in inspect.stack()[1:] if "app" in fr.filename and "style.py" not in fr.filename), None)
        tag = f"{Path(caller.filename).stem}_{caller.function}_{caller.lineno}" if caller else "unknown"
        d = out / STATE["page"]
        d.mkdir(parents=True, exist_ok=True)
        STATE["n"] += 1
        try:
            f.savefig(d / f"{STATE['n']:02d}_{tag}.png", dpi=80, bbox_inches="tight")
        except Exception as e:  # noqa: BLE001
            (d / "errors.txt").open("a").write(f"save {tag}: {e}\n")
    return _close(fig) if fig is not None else _close()


plt.close = close
sys.path.insert(0, str(REPO / "app"))
from streamlit.testing.v1 import AppTest  # noqa: E402

import importlib  # noqa: E402

PAGES = None
for key in pages:
    STATE["page"], STATE["n"] = key, 0
    at = AppTest.from_file(str(REPO / "app" / "dashboard.py"), default_timeout=900)
    if PAGES is None:
        import ast
        src = (REPO / "app" / "dashboard.py").read_text(encoding="utf-8")
        line = next(l for l in src.splitlines() if l.startswith("PAGES = "))
        PAGES = ast.literal_eval(line.split("=", 1)[1].strip())
    at.session_state["page"] = PAGES[key]
    at.session_state["run"] = run_name
    try:
        at.run()
    except Exception:  # noqa: BLE001
        (out / key).mkdir(parents=True, exist_ok=True)
        (out / key / "errors.txt").open("a").write(traceback.format_exc())
    errs = [e.value for e in at.exception] if at.exception else []
    if errs:
        (out / key).mkdir(parents=True, exist_ok=True)
        (out / key / "errors.txt").open("a").write("\n\n".join(str(e) for e in errs))
    print(f"{key}: {STATE['n']} figures, {len(errs)} exceptions", flush=True)
