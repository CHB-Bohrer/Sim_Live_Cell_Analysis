"""Run-directory helpers, config loading and provenance. Every stage writes via these."""
from __future__ import annotations

import copy
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]


def load_config(path: str | Path, overrides: list[str] | None = None) -> dict:
    """Load YAML and apply dotted overrides, e.g. ['motion.D_um2_s=0.05', 'seed=3']."""
    cfg = yaml.safe_load(Path(path).read_text())
    for item in overrides or []:
        key, _, raw = item.partition("=")
        node = cfg
        parts = key.split(".")
        for p in parts[:-1]:
            node = node.setdefault(p, {})
        node[parts[-1]] = yaml.safe_load(raw)
    return cfg


def run_id_for(cfg: dict) -> str:
    h = hashlib.sha1(json.dumps(cfg, sort_keys=True, default=str).encode()).hexdigest()[:6]
    return f"{cfg.get('name', 'run')}_seed{cfg['seed']}_{h}"


def stage_dir(run: Path, name: str) -> Path:
    d = run / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def write_provenance(d: Path, cfg: dict, seed: int) -> None:
    """params.yaml (resolved config), seed.txt, provenance.json (commit, versions, time)."""
    (d / "params.yaml").write_text(yaml.safe_dump(copy.deepcopy(cfg), sort_keys=False))
    (d / "seed.txt").write_text(str(seed))
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, capture_output=True,
                                text=True).stdout.strip()
    except OSError:
        commit = ""
    versions = {}
    for pkg in ("numpy", "scipy", "torch", "trackastra", "traccuracy", "scikit-image", "openmm"):
        try:
            versions[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            pass
    (d / "provenance.json").write_text(json.dumps({
        "git_commit": commit, "time": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version.split()[0], "platform": platform.platform(), "versions": versions}, indent=2))
