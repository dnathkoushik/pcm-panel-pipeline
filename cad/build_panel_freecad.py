"""FreeCAD (Part/OCC) builder for PCM ceiling panels.

Run headless with FreeCAD's interpreter; parameters come from a JSON file whose
path is given in the environment variable PCM_PANEL_JSON (freecadcmd treats extra
command-line arguments as documents to open, so env vars are the robust route):

    PCM_PANEL_JSON=runs/<case>/cad/panel.json freecadcmd cad/build_panel_freecad.py

The JSON is written by `pcmpipe.cad.write_cad_input` and holds the resolved
geometry from `pcmpipe.geometry` (mm). Outputs (in the JSON's out_dir):
    panel_outer.stl   envelope pushed `ceiling_extension` mm through the ceiling
                      (for snappyHexMesh), metres, solid name "panel_outer"
    pcm_cavity.stl    PCM cavity (inward offset of the envelope by the shell
                      thickness, sharp joins), metres, solid name "pcm_cavity"
    panel.step        shell + PCM solids, true geometry (mm), for manufacturing
    cad_report.json   OCC volumes and areas for QA

Construction (sharp-join offset done analytically, no OCC offset needed):
    outer = Slab(fp, -hb..0) U [Prism(fp) n U_i Cone_i]
    inner = [Slab(fp', -hb+t..-t) if flat bottom] U [Prism(fp') n U_i Cone'_i]
with fp' = fp inset by t, Cone_i the protrusion extended upwards, and Cone'_i the
same cone lifted by t / cos(face tilt) (planar faces) or the spheroid shrunk by t.
"""
import json
import math
import os
import sys

import FreeCAD as App  # noqa: F401  (provided by freecadcmd)
import Part
from FreeCAD import Vector, Matrix


def prism(poly, z0, z1):
    wire = Part.makePolygon([Vector(x, y, z0) for x, y in poly] + [Vector(poly[0][0], poly[0][1], z0)])
    return Part.Face(wire).extrude(Vector(0, 0, z1 - z0))


def inset(poly, d):
    n = len(poly)
    lines = []
    for i in range(n):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
        ex, ey = x2 - x1, y2 - y1
        L = math.hypot(ex, ey)
        lines.append(((x1 - d * ey / L, y1 + d * ex / L), (ex, ey)))
    out = []
    for i in range(n):
        (p, r), (q, s) = lines[i - 1], lines[i]
        den = r[0] * s[1] - r[1] * s[0]
        t = ((q[0] - p[0]) * s[1] - (q[1] - p[1]) * s[0]) / den
        out.append((p[0] + t * r[0], p[1] + t * r[1]))
    return out


def cone(apex, base_poly, z_base, z_top):
    """Pyramid with apex `apex` through `base_poly` at z_base, extended to z_top."""
    ax, ay, az = apex
    k = (z_top - az) / (z_base - az)
    top = [Vector(ax + k * (x - ax), ay + k * (y - ay), z_top) for x, y in base_poly]
    A = Vector(ax, ay, az)
    faces = [Part.Face(Part.makePolygon([A, top[i], top[(i + 1) % len(top)], A]))
             for i in range(len(top))]
    faces.append(Part.Face(Part.makePolygon(top + [top[0]])))
    solid = Part.Solid(Part.Shell(faces))
    if solid.Volume < 0:          # inconsistent face orientation -> flip
        solid.reverse()
    return solid


def spheroid(cx, cy, cz, a, c):
    s = Part.makeSphere(a)
    m = Matrix()
    m.scale(1.0, 1.0, c / a)
    s = s.transformGeometry(m)
    s.translate(Vector(cx, cy, cz))
    return s


def fuse_all(shapes):
    out = shapes[0]
    if len(shapes) > 1:
        out = out.fuse(shapes[1:])
    return out.removeSplitter()


def build(p, z_top, inner=False):
    t = p["shell_thickness"] if inner else 0.0
    hb, hp = p["slab_thickness"], p["protrusion_height"]
    fp = inset(p["footprint"], t) if inner else p["footprint"]
    z_cap = z_top - t if inner else z_top
    z_low = -(hb + hp) - 1.0
    lift = t / math.cos(math.radians(p["face_tilt_deg"])) if inner else 0.0
    prots = []
    for pr in p["protrusions"]:
        cx, cy = pr["centre"]
        if pr["kind"] == "spheroid":
            prots.append(spheroid(cx, cy, -hb, pr["radius"] - t, pr["height"] - t))
        else:
            apex = (cx, cy, -(hb + pr["height"]) + lift)
            base = pr["base_polygon"]
            # base polygon of the lifted cone at z = -hb + lift is the original one
            prots.append(cone(apex, base, -hb + lift, z_cap))
    parts = [prism(fp, z_low, z_cap).common(fuse_all(prots))]
    if not inner or p["has_flat_bottom"]:
        parts.append(prism(fp, -hb + t, z_cap))
    return fuse_all(parts)


def write_stl(shape, path, name, tol_mm):
    pts, tris = shape.tessellate(tol_mm)
    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for i, j, k in tris:
            a, b, c = pts[i], pts[j], pts[k]
            n = (b - a).cross(c - a)
            if n.Length < 1e-14:
                continue
            n.normalize()
            f.write(f" facet normal {n.x:.6e} {n.y:.6e} {n.z:.6e}\n  outer loop\n")
            for v in (a, b, c):
                f.write(f"   vertex {v.x * 1e-3:.9e} {v.y * 1e-3:.9e} {v.z * 1e-3:.9e}\n")
            f.write("  endloop\n endfacet\n")
        f.write(f"endsolid {name}\n")


def main():
    cfg_path = os.environ.get("PCM_PANEL_JSON")
    if not cfg_path:
        sys.exit("PCM_PANEL_JSON is not set")
    with open(cfg_path) as f:
        cfg = json.load(f)
    p, out = cfg["geometry"], cfg["out_dir"]
    os.makedirs(out, exist_ok=True)

    t, ext = p["shell_thickness"], cfg["ceiling_extension_mm"]
    outer_true = build(p, 0.0)
    outer_mesh = build(p, ext)
    cavity = build(p, 0.0, inner=True)
    shell = outer_true.cut(cavity)
    # Screening mode (shell as a thin thermal resistance): the PCM reaches the
    # ceiling, i.e. the cavity with its top raised by t (and through the ceiling
    # for meshing). Inside the room this is cavity + footprint' x t.
    cavity_thin = build(p, t, inner=True)
    cavity_thin_mesh = build(p, ext + t, inner=True)
    for name, s in (("outer", outer_true), ("cavity", cavity), ("shell", shell),
                    ("cavity_thin", cavity_thin)):
        if not s.isValid():
            sys.exit(f"invalid {name} solid")

    tol = cfg["stl_tolerance_mm"]
    write_stl(outer_mesh, os.path.join(out, "panel_outer.stl"), "panel_outer", tol)
    write_stl(cavity, os.path.join(out, "pcm_cavity.stl"), "pcm_cavity", tol)
    write_stl(cavity_thin_mesh, os.path.join(out, "pcm_cavity_thin.stl"), "pcm_cavity", tol)
    Part.Compound([shell, cavity]).exportStep(os.path.join(out, "panel.step"))

    top_area = abs(Part.Face(Part.makePolygon(
        [Vector(x, y, 0) for x, y in p["footprint"]] + [Vector(*p["footprint"][0], 0)])).Area)
    rho_pcm = cfg.get("rho_pcm", 840.0)
    report = {
        "envelope_volume_ml": outer_true.Volume / 1000.0,
        "cavity_volume_ml": cavity.Volume / 1000.0,
        "shell_volume_ml": shell.Volume / 1000.0,
        "pcm_mass_g": rho_pcm * cavity.Volume * 1e-9 * 1000.0,
        "pcm_mass_thin_g": rho_pcm * cavity_thin.Volume * 1e-9 * 1000.0,
        "cavity_thin_area_cm2": (cavity_thin.Area - abs(Part.Face(Part.makePolygon(
            [Vector(x, y, 0) for x, y in inset(p["footprint"], t)]
            + [Vector(*inset(p["footprint"], t)[0], 0)])).Area)) / 100.0,
        "exposed_area_cm2": (outer_true.Area - top_area) / 100.0,
        "cavity_area_cm2": cavity.Area / 100.0,
    }
    with open(os.path.join(out, "cad_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print("CAD OK", json.dumps(report))


main()
