"""Every dashboard page opens without an exception (Streamlit AppTest, headless).

The pages need a finished run with tracking, which takes a GPU and minutes to make, so the test uses saved runs when they exist
(the showcase runs named below, or the run named in $SIMLIVE_TEST_RUNS, comma-separated) and is skipped otherwise. Sections behind
toggles and buttons are not opened. The "no runs yet" home page is checked with an empty data folder.
"""
import ast
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SRC = (REPO / "app" / "dashboard.py").read_text(encoding="utf-8")
PAGES = ast.literal_eval(next(line for line in SRC.splitlines() if line.startswith("PAGES = ")).split("=", 1)[1].strip())
SHOWCASE = ["loci_library_seed12_87a16d", "loci_3d_seed23_2163f4"]       # 2D library run with two copies + strays; 3D z-stack run with beads


def _runs() -> list[str]:
    import os

    from simlive.io.runs import DATA_ROOT

    names = [n for n in os.environ.get("SIMLIVE_TEST_RUNS", "").split(",") if n] or SHOWCASE
    return [n for n in names if (DATA_ROOT / "runs" / n / "stage7_validation" / "metrics.json").exists()]


@pytest.fixture(scope="module")
def apptest():
    pytest.importorskip("streamlit.testing.v1")
    sys.path.insert(0, str(REPO / "app"))
    from streamlit.testing.v1 import AppTest

    def run(run_name: str | None, page: str | None):
        at = AppTest.from_file(str(REPO / "app" / "dashboard.py"), default_timeout=600)
        if page:
            at.session_state["page"] = page
        if run_name:
            at.session_state["run"] = run_name
        at.run()
        return at
    yield run
    sys.path.remove(str(REPO / "app"))


@pytest.mark.parametrize("key", list(PAGES))
def test_every_page_opens(apptest, key):
    runs = _runs()
    if not runs:
        pytest.skip("no saved showcase run on this machine (set SIMLIVE_TEST_RUNS to a finished run)")
    for name in runs:
        at = apptest(name, PAGES[key])
        assert not at.exception, f"{key} on {name}: " + "\n".join(str(e.value) for e in at.exception)
        assert at.session_state["page"] == PAGES[key]


def test_home_without_runs(apptest, tmp_path, monkeypatch):
    import simlive.io.runs as runs_mod

    monkeypatch.setattr(runs_mod, "DATA_ROOT", tmp_path)
    at = apptest(None, None)
    assert not at.exception, "\n".join(str(e.value) for e in at.exception)
    assert any("tutorial" in m.value for m in at.markdown), "the tutorial link is missing on the empty home page"
