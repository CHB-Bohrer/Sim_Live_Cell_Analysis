"""The tutorial's bookkeeping stays consistent: every screenshot the scenes use exists, has a row in scripts/tutorial/shot_table.py, and that
row names capture keys and dashboard source files that exist (scripts/tutorial/update_tutorial.py relies on it to say what is stale)."""
import json
import re
import sys

from simlive.io.runs import REPO_ROOT

TUT = REPO_ROOT / "docs" / "tutorial"
sys.path.insert(0, str(REPO_ROOT / "scripts" / "tutorial"))


def _shots() -> dict:
    return json.loads((TUT / "shots_data.js").read_text(encoding="utf-8").split("=", 1)[1].rstrip().rstrip(";"))


def test_every_screenshot_has_a_row_and_a_picture():
    from shot_table import RUNS, SHOTS

    have = _shots()
    assert set(have) == set(SHOTS), f"shots_data.js and shot_table.py differ: {set(have) ^ set(SHOTS)}"
    capture = (REPO_ROOT / "scripts" / "tutorial" / "capture_shots.py").read_text(encoding="utf-8")
    wanted = set(re.findall(r'want\("(\w+)"\)', capture)) | set(re.findall(r'^\s+\("(sb_\w+)"', capture, re.M))
    for name, (keys, run, sources, _) in SHOTS.items():
        assert (TUT / have[name]["src"]).exists(), name
        assert run in RUNS and RUNS[run] in capture, f"{name}: run {run} is not one of capture_shots.py's"
        assert set(keys) <= wanted, f"{name}: capture_shots.py has no want() for {set(keys) - wanted}"
        for s in sources:
            assert (REPO_ROOT / s).exists(), f"{name}: {s} does not exist"


def test_scenes_only_use_screenshots_that_exist():
    have = _shots()
    text = (TUT / "scenes_a.js").read_text(encoding="utf-8") + (TUT / "scenes_b.js").read_text(encoding="utf-8")
    used = set(re.findall(r'shot:\s*"(\w+)"', text))
    assert used and used <= set(have), f"scenes use screenshots that are not there: {used - set(have)}"
