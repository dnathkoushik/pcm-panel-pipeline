"""Step 1 -- geometry: resolve dimensions and build the CAD with FreeCAD."""
from __future__ import annotations

import json
import os
import shutil
import subprocess

from . import config, foam, geometry

CAD_SCRIPT = os.path.join(config.ROOT, "cad", "build_panel_freecad.py")


def write_cad_input(case: dict) -> str:
    geo = geometry.build_panel(case["panel"])
    cad_dir = os.path.join(case["dir"], "cad")
    os.makedirs(cad_dir, exist_ok=True)
    payload = {
        "geometry": geo.to_dict(),
        "out_dir": cad_dir,
        "stl_tolerance_mm": case["mesh"]["stl_tolerance_mm"],
        "ceiling_extension_mm": case["mesh"]["ceiling_extension_mm"],
        "rho_pcm": case["physics"]["pcm"]["rho"],
    }
    path = os.path.join(cad_dir, "panel.json")
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)
    print("  " + geometry.summary(geo))
    return path


def build(case: dict, dry: bool = False):
    path = write_cad_input(case)
    cad_dir = os.path.dirname(path)
    exe = foam.tool(case, "freecadcmd")
    print(f"  [cad] {' '.join(exe)} {CAD_SCRIPT}")
    if dry:
        return
    if shutil.which(exe[0]) is None:
        raise RuntimeError(f"FreeCAD command '{exe[0]}' not found: install FreeCAD or set "
                           "PCMPIPE_FREECADCMD / tools.freecadcmd in config/base_case.yaml")
    env = dict(os.environ, PCM_PANEL_JSON=path)
    with open(os.path.join(cad_dir, "log.freecad"), "w") as log:
        r = subprocess.run([*exe, CAD_SCRIPT], env=env, stdout=log, stderr=subprocess.STDOUT)
    report = os.path.join(cad_dir, "cad_report.json")
    if r.returncode != 0 or not os.path.exists(report):
        raise RuntimeError(f"FreeCAD build failed, see {cad_dir}/log.freecad")
    with open(report) as f:
        rep = json.load(f)
    print(f"  [cad] cavity {rep['cavity_volume_ml']:.1f} mL -> PCM {rep['pcm_mass_g']:.1f} g, "
          f"exposed area {rep['exposed_area_cm2']:.1f} cm2")


def cad_report(case: dict) -> dict:
    with open(os.path.join(case["dir"], "cad", "cad_report.json")) as f:
        return json.load(f)


def panel_geometry(case: dict) -> geometry.PanelGeometry:
    return geometry.build_panel(case["panel"])
