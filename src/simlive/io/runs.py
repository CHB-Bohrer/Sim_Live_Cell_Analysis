"""Run-directory helpers, config loading and provenance. Every stage writes via these."""
from __future__ import annotations

import copy
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime
from functools import lru_cache
from importlib import metadata
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
# Where simulation outputs live (runs, sweeps, chromatin library). Defaults to <repo>/data, but scripts\run.cmd points
# it at C:\Users\<you>\SimLiveData so large files stay out of the OneDrive-synced repo folder.
DATA_ROOT = Path(os.environ.get("SIMLIVE_DATA", REPO_ROOT / "data"))


INSTRUMENTS = Path(__file__).resolve().parents[3] / "configs" / "instruments"


def _deep_merge(base: dict, top: dict) -> dict:
    out = dict(base)
    for k, v in top.items():
        out[k] = _deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def load_config(path: str | Path, overrides: list[str] | None = None) -> dict:
    """Load YAML and apply dotted overrides, e.g. ['motion.D_um2_s=0.05', 'seed=3']. If the config (or an override) names an
    `instrument` (a file in configs/instruments/, e.g. vutara_vxl), that instrument's settings are merged on top of the config, and the
    overrides are applied again so that explicit --set values still win."""
    cfg = _apply_overrides(yaml.safe_load(Path(path).read_text()), overrides)
    inst = cfg.get("instrument")
    if inst:
        f = INSTRUMENTS / f"{inst}.yaml"
        if not f.exists():
            raise FileNotFoundError(f"unknown instrument {inst!r}: no {f}")
        preset = yaml.safe_load(f.read_text()) or {}
        preset.pop("instrument_info", None)
        cfg = _apply_overrides(_deep_merge(cfg, preset), overrides)
    return cfg


def _apply_overrides(cfg: dict, overrides: list[str] | None) -> dict:
    for item in overrides or []:
        key, _, raw = item.partition("=")
        node = cfg
        parts = key.split(".")
        for p in parts[:-1]:
            node = node.setdefault(p, {})
        val = yaml.safe_load(raw)
        if isinstance(val, str):  # YAML reads "1e-05" as text; accept it as a number
            try:
                val = float(val)
            except ValueError:
                pass
        node[parts[-1]] = val
    return cfg


def run_id_for(cfg: dict) -> str:
    h = hashlib.sha1(json.dumps(cfg, sort_keys=True, default=str).encode()).hexdigest()[:6]
    return f"{cfg.get('name', 'run')}_seed{cfg['seed']}_{h}"


def stage_dir(run: Path, name: str) -> Path:
    d = run / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def _git(*args: str) -> str:
    try:
        # explicit UTF-8: with text=True alone Windows decodes with cp1252, and a diff holding other characters made stdout None
        return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, encoding="utf-8", errors="replace").stdout or ""
    except OSError:
        return ""


@lru_cache(maxsize=1)
def git_state() -> dict:
    """Commit of the code that is running, whether the working tree differed from it (tracked or untracked files) and a hash of that
    difference, so a run made from uncommitted code can be told apart from a clean one. Taken once per process (at the first call)."""
    commit = _git("rev-parse", "HEAD").strip()
    dirty = bool(_git("status", "--porcelain").strip())
    diff_hash = hashlib.sha256(_git("diff", "HEAD").encode("utf-8", "replace")).hexdigest() if dirty else ""
    return {"git_commit": commit, "git_dirty": dirty, "git_diff_sha256": diff_hash}


def write_provenance(d: Path, cfg: dict, seed: int) -> None:
    """params.yaml (resolved config), seed.txt, provenance.json (commit + dirty flag, versions, time)."""
    (d / "params.yaml").write_text(yaml.safe_dump(copy.deepcopy(cfg), sort_keys=False))
    (d / "seed.txt").write_text(str(seed))
    versions = {}
    for pkg in ("numpy", "scipy", "torch", "trackastra", "traccuracy", "scikit-image", "openmm"):
        try:
            versions[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            pass
    (d / "provenance.json").write_text(json.dumps({
        **git_state(), "time": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version.split()[0], "platform": platform.platform(), "versions": versions}, indent=2))
