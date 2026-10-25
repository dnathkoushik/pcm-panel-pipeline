"""Configuration loading, study expansion and case naming."""
from __future__ import annotations

import copy
import itertools
import json
import os

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(ROOT, "config")
RUNS_DIR = os.environ.get("PCMPIPE_RUNS_DIR") or os.path.join(ROOT, "runs")


def _load(name):
    with open(os.path.join(CONFIG_DIR, name)) as f:
        return yaml.safe_load(f)


def set_dotted(d: dict, key: str, value):
    parts = key.split(".")
    for p in parts[:-1]:
        d = d.setdefault(p, {})
    d[parts[-1]] = value


def get_dotted(d: dict, key: str):
    for p in key.split("."):
        d = d[p]
    return d


def case_id(case: dict) -> str:
    p = case["panel"]
    tag = {"pyramid": "pyr", "tetrahedron": "tet", "hemisphere": "hem"}[p["shape"]]
    cid = (f"{tag}_V{p['volume_ml']:g}_r{round(100 * p['protrusion_ratio'])}"
           f"_t{p['shell_thickness_mm']:g}_g{p['gap_mm']:g}_dT{case['delta_T']:g}")
    if case["mesh"]["preset"] != "medium":
        cid += f"_{case['mesh']['preset']}"
    if case.get("id_suffix") and not cid.endswith("_" + case["id_suffix"]):  # override studies get own folders
        cid += f"_{case['id_suffix']}"
    return cid


def resolve(case: dict) -> dict:
    """Fill derived values (mesh preset, air temperature) in place."""
    case["mesh"]["active"] = copy.deepcopy(case["mesh"]["presets"][case["mesh"]["preset"]])
    ic = case["physics"]["initial_conditions"]
    ic["T_air"] = ic["T_pcm"] + float(case["delta_T"])
    air = case["physics"]["air"]
    air["T_ref"] = ic["T_air"]                     # Boussinesq reference, Eq. 17
    if air["beta"] == "auto":
        air["beta"] = 1.0 / air["T_ref"]
    case["id"] = case_id(case)
    case["dir"] = os.path.join(RUNS_DIR, case["id"])
    return case


def base_case() -> dict:
    case = _load("base_case.yaml")
    case["physics"] = _load("physics.yaml")
    return case


def study_names():
    return list(_load("studies.yaml").keys())


def expand_study(name: str) -> list[dict]:
    studies = _load("studies.yaml")
    if name not in studies:
        raise KeyError(f"unknown study '{name}', choose from {list(studies)}")
    spec = studies[name]
    keys = list(spec.get("vary", {}).keys())
    cases = []
    for combo in itertools.product(*(spec["vary"][k] for k in keys)):
        case = base_case()
        for k, v in spec.get("override", {}).items():
            set_dotted(case, k, v)
        if spec.get("override"):
            case["id_suffix"] = name
        for k, v in zip(keys, combo):
            set_dotted(case, k, v)
        case["study"] = name
        cases.append(resolve(case))
    return cases


def save_case(case: dict):
    os.makedirs(case["dir"], exist_ok=True)
    with open(os.path.join(case["dir"], "case.json"), "w") as f:
        json.dump(case, f, indent=2, default=str)


def load_case(case_dir: str) -> dict:
    with open(os.path.join(case_dir, "case.json")) as f:
        return json.load(f)
