"""Configuration: roof/config/roof.yaml plus `--set KEY=VALUE` overrides."""
from __future__ import annotations

import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))       # roof/
CONFIG = os.path.join(ROOT, "config", "roof.yaml")
TEMPLATES = os.path.join(ROOT, "templates")


def runs_dir() -> str:
    """Output root: $ROOF_RUNS_DIR (e.g. on /scratch), else roof/runs."""
    return os.path.abspath(os.environ.get("ROOF_RUNS_DIR") or os.path.join(ROOT, "runs"))


def run_dir(name_or_path: str) -> str:
    """A run given by its name (under runs_dir()) or by a path."""
    if os.path.isdir(name_or_path):
        return os.path.abspath(name_or_path)
    return os.path.join(runs_dir(), name_or_path)


def _yaml():
    try:
        import yaml
    except ImportError:
        raise SystemExit("roofpipe needs PyYAML: `pip install --user pyyaml`, or load a "
                         "Python/Anaconda module that has it (slurm/env.sh: PYTHON_MODULE)")
    return yaml


def load(path: str | None = None, overrides=()) -> dict:
    with open(path or CONFIG, encoding="utf-8") as f:
        cfg = _yaml().safe_load(f)
    for item in overrides:
        apply_override(cfg, item)
    return cfg


def apply_override(cfg: dict, item: str):
    """`a.b.c=value`; the value is read as YAML. Only existing keys can be set."""
    if "=" not in item:
        raise SystemExit(f"--set expects KEY=VALUE, got {item!r}")
    key, raw = item.split("=", 1)
    parts = key.strip().split(".")
    node = cfg
    for p in parts[:-1]:
        if not isinstance(node.get(p), dict):
            raise SystemExit(f"--set {key}: the configuration has no section {p!r}")
        node = node[p]
    if parts[-1] not in node:
        raise SystemExit(f"--set {key}: unknown key {parts[-1]!r}")
    value = _yaml().safe_load(raw)
    old = node[parts[-1]]
    # YAML 1.1 reads 1e-8 as a string: keep numbers numbers
    if isinstance(old, (int, float)) and not isinstance(old, bool) and isinstance(value, str):
        try:
            value = float(value)
        except ValueError:
            raise SystemExit(f"--set {key}: {raw!r} is not a number")
    node[parts[-1]] = value


def dump(cfg: dict, path: str):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        _yaml().safe_dump(cfg, f, sort_keys=False)
