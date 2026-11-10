"""One command to keep the narrated tutorial (docs/tutorial) in step with the dashboard.

    scripts\\tutorial\\update.cmd                 = status: which screenshots are out of date, and a check of every sentence
    scripts\\tutorial\\update.cmd all             = re-take the stale screenshots, record the voice of edited sentences, check everything
    scripts\\tutorial\\update.cmd capture [shot ...]   re-take the named screenshots (no name: the stale ones; --all: every one)
    scripts\\tutorial\\update.cmd voice           record the sentences whose text changed (Kokoro, on this PC)
    scripts\\tutorial\\update.cmd check           open every sentence in the player: script errors, missing pictures / rectangles, long sentences, voice clips
    scripts\\tutorial\\update.cmd numbers         list every sentence that states a number (to compare with the current results)

(update.cmd runs this file with the screenshot venv's Python, which has playwright; the simlive env is only used to start the dashboard.)

How "stale" is decided: after every capture, docs/tutorial/shots_stamp.json records for each picture a hash of the dashboard source files that
draw it (scripts/tutorial/shot_table.py) and the newest result file of the run / sweeps it shows. A picture is stale when one of those changed.
The capture starts its own dashboard on port 8502 if none is running there and stops only that one afterwards.
"""
from __future__ import annotations

import datetime
import fnmatch
import hashlib
import json
import os
import pathlib
import re
import socket
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[1]
TUT = REPO / "docs" / "tutorial"
STAMP = TUT / "shots_stamp.json"
DATA = pathlib.Path(os.environ.get("SIMLIVE_DATA", pathlib.Path.home() / "SimLiveData"))
KOKORO_PY = pathlib.Path.home() / "kokoro_tmp" / "v" / "Scripts" / "python.exe"
PORT = 8502
MAX_CHARS = 250

sys.path.insert(0, str(HERE))
from shot_table import RUNS, SHOTS  # noqa: E402


# ----------------------------------------------------------------------------- what a picture depends on
def _hash(rel: str) -> str:
    p = REPO / rel
    return hashlib.sha1(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()[:12] if p.exists() else "missing"


def _newest(files) -> float:
    return max((f.stat().st_mtime for f in files if f.is_file()), default=0.0)


def run_signature(run: str) -> float:
    """Newest result file of a run (stage 6 and 7 tables and scores): re-analysing a run makes its pictures stale."""
    r = DATA / "runs" / run
    files = [f for d in ("stage6_analysis", "stage7_validation") for pat in ("*.csv", "*.json", "*/*.csv", "*/*.json") for f in (r / d).glob(pat)
             if "before_" not in str(f)]
    return _newest(files)


def sweep_signature(pattern: str) -> float:
    if not pattern:
        return 0.0
    return _newest(d / "results.csv" for d in (DATA / "sweeps").iterdir() if d.is_dir() and fnmatch.fnmatch(d.name, pattern)) if (DATA / "sweeps").exists() else 0.0


def fingerprint(shot: str) -> dict:
    keys, run, sources, sweeps = SHOTS[shot]
    return {"run": RUNS[run], "sources": {s: _hash(s) for s in sources}, "run_results": run_signature(RUNS[run]), "sweeps": sweep_signature(sweeps)}


def load_stamp() -> dict:
    return json.loads(STAMP.read_text(encoding="utf-8")) if STAMP.exists() else {}


def shots_in_tutorial() -> dict:
    txt = (TUT / "shots_data.js").read_text(encoding="utf-8")
    return json.loads(txt.split("=", 1)[1].rstrip().rstrip(";"))


def stale_shots(verbose: bool = True) -> list[str]:
    stamp, have = load_stamp(), shots_in_tutorial()
    out = []
    for name in sorted(set(have) - set(SHOTS)):
        print(f"  ?  {name}: in the tutorial but not in scripts/tutorial/shot_table.py (add a row)")
    for name in SHOTS:
        why = []
        if name not in have:
            why.append("no picture yet")
        elif name not in stamp:
            why.append("never stamped (captured before 2026-11-09)")
        else:
            old, new = stamp[name], fingerprint(name)
            why += [f"{s} changed" for s, h in new["sources"].items() if old["sources"].get(s) != h]
            if new["run"] != old["run"]:
                why.append(f"run is now {new['run']}")
            elif new["run_results"] > old["run_results"] + 1:
                why.append(f"results of {new['run']} are newer ({datetime.datetime.fromtimestamp(new['run_results']):%Y-%m-%d %H:%M})")
            if new["sweeps"] > old["sweeps"] + 1:
                why.append("its sweeps were re-run")
        if why:
            out.append(name)
        if verbose:
            when = stamp.get(name, {}).get("captured", "?")
            print(f"  {'STALE' if why else 'ok   '} {name:18s} captured {when:16s} {'; '.join(why)}")
    return out


# ----------------------------------------------------------------------------- capture
def _listening(port: int) -> bool:
    with socket.socket() as s:
        s.settimeout(0.5)
        return s.connect_ex(("127.0.0.1", port)) == 0


def capture(names: list[str]) -> None:
    if not names:
        print("no screenshot to take")
        return
    keys = sorted({k for n in names for k in SHOTS[n][0]})
    out = DATA / "tutorial" / f"shots_{datetime.datetime.now():%Y-%m-%d_%H%M}"
    proc = None
    if not _listening(PORT):
        print(f"starting a dashboard on port {PORT} for the screenshots ...", flush=True)
        log = open(DATA / "tutorial" / "capture_dashboard.log", "w", encoding="utf-8")
        proc = subprocess.Popen(["cmd", "/c", str(REPO / "scripts" / "run.cmd"), "python", "-m", "streamlit", "run", "app/dashboard.py", "--server.port", str(PORT),
                                 "--server.address", "127.0.0.1", "--server.headless", "true", "--browser.gatherUsageStats", "false"], cwd=REPO, stdout=log, stderr=log)
        t0 = time.time()
        while not _listening(PORT):
            if time.time() - t0 > 120 or proc.poll() is not None:
                raise SystemExit(f"the dashboard did not start: see {log.name}")
            time.sleep(1)
        time.sleep(5)
    try:
        print(f"capturing {len(names)} screenshot(s) into {out} (keys: {' '.join(keys)})", flush=True)
        env = {**os.environ, "TUT_URL": f"http://localhost:{PORT}"}
        subprocess.run([sys.executable, str(HERE / "capture_shots.py"), str(out), *keys], check=True, env=env)
    finally:
        if proc is not None:                       # stop only the dashboard started here (the whole process tree of that one PID)
            subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True)
    subprocess.run([sys.executable, str(HERE / "build_data.py"), str(out)], check=True)
    stamp = load_stamp()
    now = f"{datetime.datetime.now():%Y-%m-%d %H:%M}"
    got = {p.stem for p in out.glob("*.json")}
    for n in SHOTS:
        if n in got:
            stamp[n] = {"captured": now, **fingerprint(n)}
    STAMP.write_text(json.dumps(stamp, indent=1) + "\n", encoding="utf-8", newline="\n")
    missing = [n for n in names if n not in got]
    if missing:
        print("NOT captured (see the messages above):", " ".join(missing))


# ----------------------------------------------------------------------------- every sentence
def _player(pw):
    b = pw.chromium.launch(channel="msedge", headless=True)
    pg = b.new_page(viewport={"width": 1600, "height": 900})
    msgs = []
    pg.on("console", lambda m: msgs.append((m.type, m.text)) if m.type in ("error", "warning") else None)
    pg.on("pageerror", lambda e: msgs.append(("error", str(e))))
    pg.goto(TUT.joinpath("index.html").as_uri() + "#s=welcome&b=0")
    pg.wait_for_function("typeof T !== 'undefined' && T.scenes.length > 0 && T.si >= 0")
    return b, pg, msgs


def beats() -> list[dict]:
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        b, pg, _ = _player(pw)
        out = pg.evaluate("T.scenes.flatMap((s, i) => s.beats.map((b, j) => ({si: i, bi: j, scene: s.id, key: s.id + ':' + j, text: b.s, h: hashText(b.s)})))")
        b.close()
    return out


def check() -> int:
    """Show every sentence once in the player and collect what goes wrong. Returns the number of problems."""
    from playwright.sync_api import sync_playwright

    problems = []
    with sync_playwright() as pw:
        b, pg, msgs = _player(pw)
        bs = pg.evaluate("T.scenes.flatMap((s, i) => s.beats.map((b, j) => ({si: i, bi: j, key: s.id + ':' + j, text: b.s, h: hashText(b.s)})))")
        used = set(pg.evaluate("[...new Set(T.scenes.map(s => s.id))]"))
        for bt in bs:
            del msgs[:]
            pg.evaluate("([i, j]) => T.goto(i, j)", [bt["si"], bt["bi"]])
            pg.wait_for_timeout(60)
            for kind, text in msgs:
                if "favicon" not in text and "speechSynthesis" not in text:
                    problems.append(f"{bt['key']}: {kind}: {text[:160]}")
        b.close()
    long = [bt["key"] for bt in bs if len(bt["text"]) > MAX_CHARS]
    clips, no_clip = {}, []
    mf = TUT / "audio" / "manifest.js"
    if mf.exists():
        txt = mf.read_text(encoding="utf-8").split("=", 1)[1]
        clips = json.loads(txt[: txt.rindex(";")])
    for bt in bs:
        c = clips.get(bt["key"])
        if not c or c["h"] != bt["h"] or not (TUT / c["f"]).exists():
            no_clip.append(bt["key"])
    print(f"checked {len(bs)} sentences in {len(used)} scenes")
    for p in problems:
        print("  PROBLEM", p)
    if long:
        print(f"  note: {len(long)} sentence(s) longer than {MAX_CHARS} characters (fine with the recorded voice; some browser voices cut them off): {' '.join(long[:8])}{' ...' if len(long) > 8 else ''}")
    if no_clip:
        print(f"  {len(no_clip)} sentence(s) without a recorded voice clip (the browser voice reads them): {' '.join(no_clip[:12])}{' ...' if len(no_clip) > 12 else ''}"
              "\n  -> scripts\\tutorial\\update.cmd voice")
    if not problems:
        print("  no script errors, no missing pictures or rectangles")
    return len(problems)


def voice() -> None:
    if not KOKORO_PY.exists():
        print(f"voice: no Kokoro environment at {KOKORO_PY} (set-up: header of scripts/tutorial/video/make_clips_kokoro.py); the browser voice reads edited sentences")
        return
    f = pathlib.Path.home() / "kokoro_tmp" / "beats.json"
    f.write_text(json.dumps([{k: b[k] for k in ("key", "text", "h")} for b in beats()]), encoding="utf-8")
    subprocess.run([str(KOKORO_PY), str(HERE / "video" / "make_clips_kokoro.py"), str(f), "af_heart"], check=True)


_N = (r"(?:\d[\d.,]*|half|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|fifteen|twenty|twenty-\w+|thirty|forty|forty-\w+|fifty|fifty-\w+|sixty|"
      r"hundred|thousand|thousandths|nought|zero|point)")
_U = (r"(?:percent|nanomet\w+|microns?|pixels?|seconds?|minutes?|hours?|photons?|kilobases?|megabases?|monomers|fold|times|checks|sources|movies|strays|"
      r"cells|planes|frames|chains|stages|scenes|of (?:twelve|ten|\d+))")
QUANTITY = re.compile(rf"\b{_N}(?:[ -](?:{_N}|a|and|minus))*\s+{_U}\b|\bpercent\b|\bpoint (?:{_N})\b", re.I)


def numbers() -> None:
    """Sentences that state a quantity: the ones to re-read when results change. Written to docs/tutorial/NUMBERS.md."""
    rows = [b for b in beats() if QUANTITY.search(b["text"])]
    lines = ["# Sentences of the tutorial that state a number", "",
             "Made by `scripts\\tutorial\\update.cmd numbers`. When results or defaults change, re-read these against the dashboard and the docs;",
             "edit the sentence in `scenes_a.js` / `scenes_b.js` (search for its first words), then run `scripts\\tutorial\\update.cmd all`.", ""]
    scene = None
    for b in rows:
        if b["scene"] != scene:
            scene = b["scene"]
            lines += ["", f"## {scene}"]
        lines.append(f"- `{b['key']}` {b['text']}")
    (TUT / "NUMBERS.md").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(rows)} sentences with numbers -> {TUT / 'NUMBERS.md'}")


def main() -> None:
    args = sys.argv[1:]
    cmd = args[0] if args else "status"
    if cmd == "status":
        print("Screenshots:")
        st = stale_shots()
        print(f"{len(st)} stale" + (f": scripts\\tutorial\\update.cmd capture   (or: all)" if st else ""))
        print("\nSentences:")
        raise SystemExit(1 if check() else 0)
    if cmd == "capture":
        names = list(SHOTS) if "--all" in args else [a for a in args[1:] if a in SHOTS] or stale_shots(verbose=False)
        unknown = [a for a in args[1:] if a not in SHOTS and a != "--all"]
        if unknown:
            raise SystemExit(f"unknown shot(s) {unknown}; known: {' '.join(SHOTS)}")
        capture(names)
    elif cmd == "voice":
        voice()
    elif cmd == "check":
        raise SystemExit(1 if check() else 0)
    elif cmd == "numbers":
        numbers()
    elif cmd == "all":
        capture(list(SHOTS) if "--all" in args else stale_shots(verbose=False))
        voice()
        numbers()
        bad = check()
        print("\nScreenshots now:")
        stale_shots()
        raise SystemExit(1 if bad else 0)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
