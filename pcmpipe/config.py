"""Configuration loading, study expansion and case naming."""
from __future__ import annotations

import copy
import itertools
import json
import os

import numpy as np
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


def loads_active(case: dict) -> bool:
    return case.get("loads") not in (None, "none")


def case_id(case: dict) -> str:
    p = case["panel"]
    tag = {"pyramid": "pyr", "tetrahedron": "tet", "hemisphere": "hem"}[p["shape"]]
    cid = (f"{tag}_V{p['volume_ml']:g}_r{round(100 * p['protrusion_ratio'])}"
           f"_t{p['shell_thickness_mm']:g}_g{p['gap_mm']:g}")
    cid += f"_{case['loads_spec']['id_tag']}" if loads_active(case) else f"_dT{case['delta_T']:g}"
    if case.get("material_tag"):
        cid += f"_{case['material_tag']}"
    if case["mesh"]["preset"] != "medium":
        cid += f"_{case['mesh']['preset']}"
    if case.get("model", {}).get("shell") == "thin" and case["mesh"]["preset"] != "screen":
        cid += "_thin"
    if case.get("id_suffix") and not cid.endswith("_" + case["id_suffix"]):  # override studies get own folders
        cid += f"_{case['id_suffix']}"
    return cid


def apply_material(case: dict):
    """pcm_material (config/materials.yaml) -> physics.pcm and physics.solver."""
    name = case.get("pcm_material") or "paper_cla"
    mats = _load("materials.yaml")
    if name not in mats:
        raise KeyError(f"unknown pcm_material '{name}', choose from {list(mats)}")
    mat = copy.deepcopy(mats[name] or {})
    case["material_tag"] = mat.pop("id_tag", "")
    case["physics"]["solver"].update(mat.pop("solver", {}))
    case["physics"]["pcm"].update(mat)


def apply_loads(case: dict):
    """loads (config/loads.yaml): roof scenario initial temperatures; else the paper's."""
    ic = case["physics"]["initial_conditions"]
    air = case["physics"]["air"]
    if not loads_active(case):
        case["loads_spec"] = None
        ic["T_air"] = ic["T_pcm"] + float(case["delta_T"])
        air["T_ref"] = ic["T_air"]                 # Boussinesq reference, Eq. 17
        return
    from . import loads                            # numpy only
    specs = _load("loads.yaml")
    if case["loads"] not in specs:
        raise KeyError(f"unknown loads '{case['loads']}', choose from {list(specs)}")
    spec = case["loads_spec"] = copy.deepcopy(specs[case["loads"]])
    start = spec["start_clock_h"]
    T_room = loads.room_temperature_C(spec, start + np.arange(0.0, 24.0, 0.25), ROOT)
    ic["T_pcm"] = ic["T_shell"] = spec["initial"]["T_pcm_C"] + 273.15
    ic["T_air"] = float(T_room[0]) + 273.15
    air["T_ref"] = float(np.mean(T_room)) + 273.15  # middle of the day's air temperatures


def resolve(case: dict) -> dict:
    """Fill derived values (material, loads, mesh preset, air temperature) in place."""
    case["mesh"]["active"] = copy.deepcopy(case["mesh"]["presets"][case["mesh"]["preset"]])
    if case["mesh"].get("wall_layer_passes") is not None:   # one switch for every preset
        case["mesh"]["active"]["wall_layer_passes"] = case["mesh"]["wall_layer_passes"]
    apply_material(case)
    apply_loads(case)
    air = case["physics"]["air"]
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


def mode_names():
    return list(_load("base_case.yaml").get("modes", {"paper": {}}).keys())


def parse_sets(items) -> list:
    """--set key=value (value parsed as YAML: 0.5, [1, 2], {air: 0, pcm: 0}, none)."""
    out = []
    for item in items or []:
        key, sep, value = item.partition("=")
        if not sep or not key:
            raise ValueError(f"--set expects key=value, got '{item}'")
        out.append((key.strip(), yaml.safe_load(value)))
    return out


def expand_study(name: str, mode: str = "paper", sets=None) -> list[dict]:
    """Cases of a study. Overrides apply in the order: mode, study override, vary, --set."""
    studies = _load("studies.yaml")
    if name not in studies:
        raise KeyError(f"unknown study '{name}', choose from {list(studies)}")
    spec = studies[name]
    keys = list(spec.get("vary", {}).keys())
    cases = []
    for combo in itertools.product(*(spec["vary"][k] for k in keys)):
        case = base_case()
        if mode not in case.get("modes", {}):
            raise KeyError(f"unknown mode '{mode}', choose from {list(case.get('modes', {}))}")
        for k, v in case["modes"][mode].items():
            set_dotted(case, k, v)
        case["mode"] = mode
        for k, v in spec.get("override", {}).items():
            set_dotted(case, k, v)
        if spec.get("override") and spec.get("id_suffix", True):
            case["id_suffix"] = name
        for k, v in zip(keys, combo):
            set_dotted(case, k, v)
        for k, v in sets or []:
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
