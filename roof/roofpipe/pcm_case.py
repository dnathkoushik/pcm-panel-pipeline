"""Phase 2: the concrete slab with one embedded PCM body (templates/concrete_slab with
templates/pcm_slab on top: the same slab, boundary conditions and solver settings as
the Phase 1 baseline, plus the PCM region, the snappyHexMesh set-up and PCM monitors).

    <runs>/<case_id>/                 e.g. pyr_N_h070
        cad/body.json                 geometry for FreeCAD (mm)
        cad/pcm_body.stl              FreeCAD (or roofpipe/stl.py where FreeCAD is missing)
        case/                         OpenFOAM case, meshed by scripts/run_pcm_case.sh
        boundary_conditions.csv, case_meta.json, config_used.yaml

`finalize` (after meshing) reads the meshed PCM volume from checkMesh and scales the
PCM density by V_exact / V_mesh. The geometry volumes are equal by construction; the
scaling makes the PCM mass, and so its latent store, equal as well despite the few
per cent that castellation and snapping add or remove.
"""
from __future__ import annotations

import copy
import json
import os
import re
import shutil

from . import case as base, config, geometry
from .case import num

LAYERS = ("concrete_slab", "pcm_slab")


def render_layers(layers, dst_dir: str, values: dict):
    """Render several template trees into one case; a later layer replaces files."""
    files = {}
    for layer in layers:
        src = os.path.join(config.TEMPLATES, layer)
        for dirpath, _, names in os.walk(src):
            for n in names:
                rel = os.path.normpath(os.path.join(os.path.relpath(dirpath, src), n))
                files[rel] = os.path.join(dirpath, n)
    missing = []
    for rel, src in sorted(files.items()):
        with open(src, encoding="utf-8") as f:
            text = f.read()

        def fill(m, where=rel):
            if m.group(1) not in values:
                missing.append(f"{where}: {m.group(0)}")
                return m.group(0)
            return str(values[m.group(1)])

        text = base._TOKEN.sub(fill, text)
        dst = os.path.join(dst_dir, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with open(dst, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
    if missing:
        raise SystemExit("template placeholders without a value:\n  " + "\n  ".join(missing))


def pcm_material(cfg: dict) -> dict:
    ps = cfg["pcm_study"]
    m = cfg["materials"][ps["pcm"]]
    lo, hi = (float(v) for v in ps["melting_range_C"])
    rho = float(m["rho_liquid"] if ps.get("rho_basis", "liquid") == "liquid" else m["rho_solid"])
    return {"name": ps["pcm"], "T_solidus": lo + 273.15, "T_liquidus": hi + 273.15,
            "L": float(m["latent_heat"]), "rho": rho,
            "cp_solid": float(m["cp_solid"]), "cp_liquid": float(m["cp_liquid"]),
            "k_solid": float(m["k_solid"]), "k_liquid": float(m["k_liquid"]),
            # reference values carried by the thermo package; the fvOption adds the rest
            "cp_ref": float(m["cp_liquid"]),
            "k_ref": 0.5 * (float(m["k_solid"]) + float(m["k_liquid"]))}


def background(cfg: dict, b: geometry.PcmBody):
    """Background block counts: odd nx, ny (the slab axis is a cell centre) and an nz
    for which no cell face of any refinement level lies on the base or apex plane."""
    s, mesh = cfg["slab"], cfg["pcm_study"]["mesh"]
    Lx, Ly, H = float(s["Lx"]), float(s["Ly"]), float(s["H"])
    cell = float(mesh["cell"])
    nx, ny = (max(3, round(L / cell)) for L in (Lx, Ly))
    nx += 1 - nx % 2
    ny += 1 - ny % 2
    nz = max(4, round(H / cell))
    top = int(mesh["feature_level"])

    def clear(n):
        for lvl in range(top + 1):
            dz = H / (n * 2 ** lvl)
            for z in (b.z_base, b.z_apex):
                r = z / dz
                if abs(r - round(r)) < 0.1:
                    return False
        return True

    for _ in range(20):
        if clear(nz):
            return nx, ny, nz
        nz += 1
    raise SystemExit(f"{b.case_id}: no background nz clears the PCM planes")


def _probes(b: geometry.PcmBody, H: float, n: int = 60, gap: float = 0.0015):
    """Points on the body axis every H/n, split into the concrete and the PCM sets;
    points within `gap` of the base or apex plane are left out (snapped surface)."""
    cx, cy = b.centre
    lo, hi = sorted((b.z_base, b.z_apex))
    conc, pcm = [], []
    for i in range(n):
        z = (i + 0.5) * H / n
        if lo + gap < z < hi - gap:
            pcm.append(z)
        elif z < lo - gap or z > hi + gap:
            conc.append(z)
    fmt = lambda zs: "\n".join(f"            ({num(cx + 1.37e-4)} {num(cy + 1.13e-4)} {num(z)})"  # noqa: E731
                               for z in zs)
    return fmt(conc), fmt(pcm), conc, pcm


def derive(cfg: dict, b: geometry.PcmBody, n_procs: int):
    cfg = copy.deepcopy(cfg)
    ps = cfg["pcm_study"]
    cfg["run"]["deltaT"] = float(ps["run"]["deltaT"])
    cfg["run"]["write_interval_h"] = float(ps["run"].get("write_interval_h",
                                                        cfg["run"]["write_interval_h"]))
    values, meta, dd, tab = base.derive(cfg, n_procs)
    s = cfg["slab"]
    Lx, Ly, H = float(s["Lx"]), float(s["Ly"]), float(s["H"])
    nx, ny, nz = background(cfg, b)
    pm = pcm_material(cfg)
    loc = geometry.locations(b, Lx, Ly, H)
    pc_txt, pp_txt, pc, pp = _probes(b, H)
    m = ps["mesh"]
    values.update({
        "NX": nx, "NY": ny, "NZ": nz, "CASE_ID": b.case_id,
        "N_BETWEEN": int(m["n_cells_between_levels"]), "PCM_LEVEL": int(m["pcm_level"]),
        "FEATURE_LEVEL": int(m["feature_level"]),
        "LOC_PCM": "(" + " ".join(num(v) for v in loc["pcm"]) + ")",
        "LOC_CONCRETE": "(" + " ".join(num(v) for v in loc["concrete"]) + ")",
        "PROBES_CONCRETE": pc_txt, "PROBES_PCM": pp_txt,
        "N_OUTER": int(ps["run"]["nOuterCorrectors"]),
        "PCM_NAME": pm["name"],
        "RHO_NOTE": f"{ps.get('rho_basis', 'liquid')} density, scaled by V_exact/V_mesh "
                    "after meshing (python3 -m roofpipe finalize)",
        "K_PCM_REF": num(pm["k_ref"]), "CP_PCM_REF": num(pm["cp_ref"]), "RHO_PCM": num(pm["rho"]),
        "T_SOLIDUS": num(pm["T_solidus"]), "T_LIQUIDUS": num(pm["T_liquidus"]),
        "LATENT": num(pm["L"]), "CP_SOLID": num(pm["cp_solid"]), "CP_LIQUID": num(pm["cp_liquid"]),
        "K_SOLID": num(pm["k_solid"]), "K_LIQUID": num(pm["k_liquid"]),
    })
    meta.update({
        "phase": 2, "kind": "pcm_embedded", "case_id": b.case_id,
        "body": b.to_dict(), "pcm": pm,
        "pcm_mass_exact_kg": pm["rho"] * b.volume,
        "latent_store_exact_kJ": pm["rho"] * b.volume * pm["L"] / 1000.0,
        "mesh": {"background": [nx, ny, nz], "cell": [Lx / nx, Ly / ny, H / nz],
                 "pcm_level": int(m["pcm_level"]), "feature_level": int(m["feature_level"])},
        "locations_in_mesh": loc,
        "probes": {"x": b.centre[0], "y": b.centre[1], "z_concrete": pc, "z_pcm": pp},
        "solver": {"nOuterCorrectors": int(ps["run"]["nOuterCorrectors"])},
    })
    return values, meta, cfg, tab


def build(cfg: dict, b: geometry.PcmBody, n_procs: int, force: bool = False) -> str:
    rdir = os.path.join(config.runs_dir(), b.case_id)
    cdir = os.path.join(rdir, "case")
    if " " in cdir:
        raise SystemExit(f"the case path contains spaces ({cdir}): OpenFOAM cannot handle that")
    if os.path.exists(cdir):
        if not force:
            raise SystemExit(f"{cdir} already exists (use --force to rewrite it)")
        shutil.rmtree(cdir)
        shutil.rmtree(os.path.join(rdir, "cad"), ignore_errors=True)
    os.makedirs(os.path.join(rdir, "cad"), exist_ok=True)

    values, meta, cfg_used, tab = derive(cfg, b, n_procs)
    render_layers(LAYERS, cdir, values)
    open(os.path.join(cdir, "case.foam"), "w").close()
    base._write_csv(os.path.join(rdir, "boundary_conditions.csv"), tab)

    mm = 1e3
    body_mm = {
        "case_id": b.case_id, "shape": b.shape, "orientation": b.orientation,
        "height": b.height * mm, "z_base": b.z_base * mm, "z_apex": b.z_apex * mm,
        "centre": [b.centre[0] * mm, b.centre[1] * mm],
        "footprint": [[x * mm, y * mm] for x, y in b.footprint], "radius": b.radius * mm,
        "volume_mm3": b.volume * 1e9,
    }
    s = cfg["slab"]
    with open(os.path.join(rdir, "cad", "body.json"), "w", encoding="utf-8") as f:
        json.dump({"body": body_mm, "out_dir": os.path.join(rdir, "cad"),
                   "slab": {"Lx": float(s["Lx"]) * mm, "Ly": float(s["Ly"]) * mm,
                            "H": float(s["H"]) * mm},
                   "stl_tolerance_mm": float(cfg["pcm_study"]["mesh"]["stl_tolerance_mm"])},
                  f, indent=2)
    with open(os.path.join(rdir, "case_meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    config.dump(cfg_used, os.path.join(rdir, "config_used.yaml"))
    print(f"  {b.case_id}: {geometry.SHAPES[b.shape]}, {b.orientation}, h {b.height * mm:.0f} mm, "
          f"V {b.volume * 1e6:.1f} mL, PCM {meta['pcm_mass_exact_kg'] * 1e3:.0f} g, "
          f"background {meta['mesh']['background']} -> {cdir}")
    return rdir


_VOL = re.compile(r"Total volume = ([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)")


def finalize(rdir: str) -> dict:
    """After meshing: PCM density so that the meshed mass equals the exact one."""
    cdir = os.path.join(rdir, "case")
    with open(os.path.join(rdir, "case_meta.json"), encoding="utf-8") as f:
        meta = json.load(f)
    vols = {}
    for r in ("concrete", "pcm"):
        with open(os.path.join(cdir, f"log.checkMesh.{r}"), encoding="utf-8",
                  errors="replace") as f:
            m = _VOL.findall(f.read())
        if not m:
            raise SystemExit(f"no 'Total volume' in {cdir}/log.checkMesh.{r}")
        vols[r] = float(m[-1])
    exact = meta["body"]["volume"]
    rho0 = meta["pcm"]["rho"]
    rho = rho0 * exact / vols["pcm"]
    tp = os.path.join(cdir, "constant", "pcm", "thermophysicalProperties")
    with open(tp, encoding="utf-8") as f:
        text = f.read()
    text, n = re.subn(r"(equationOfState\s*\{\s*rho\s+)[0-9.eE+-]+(;)",
                      lambda mo: f"{mo.group(1)}{num(rho)}{mo.group(2)}", text)
    if n != 1:
        raise SystemExit(f"could not set rho in {tp}")
    with open(tp, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    s = meta["slab_m"]
    meta["meshed"] = {
        "pcm_volume_m3": vols["pcm"], "concrete_volume_m3": vols["concrete"],
        "pcm_volume_error": vols["pcm"] / exact - 1.0,
        "total_volume_error": (vols["pcm"] + vols["concrete"]) / (s["Lx"] * s["Ly"] * s["H"]) - 1.0,
        "rho_pcm_used": rho, "pcm_mass_kg": rho * vols["pcm"],
    }
    with open(os.path.join(rdir, "case_meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    print(f"  {meta['case_id']}: meshed PCM {vols['pcm'] * 1e6:.1f} mL "
          f"({100 * meta['meshed']['pcm_volume_error']:+.2f} % vs exact), rho {rho0:g} -> "
          f"{rho:.2f} kg/m3, PCM mass {rho * vols['pcm'] * 1e3:.1f} g")
    return meta["meshed"]
