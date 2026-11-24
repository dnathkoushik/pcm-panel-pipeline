"""FreeCAD (Part/OCC) builder for one PCM body embedded in the roof slab (Phase 2).

Run headless with FreeCAD's interpreter. The parameters come from a JSON file
named by the environment variable ROOF_PCM_JSON (freecadcmd treats command-line
arguments as documents to open):

    ROOF_PCM_JSON=runs/<case>/cad/body.json freecadcmd roof/cad/build_roof_pcm_freecad.py

The JSON is written by `roofpipe.pcm_case` from `roofpipe.geometry.body()`. Its
lengths are in mm. Outputs, in the JSON's out_dir:
    pcm_body.stl     closed surface of the PCM body, in METRES, solid "pcm_body"
                     (snappyHexMesh: castellated and snapped into the slab mesh)
    pcm_body.step    the body and the slab with the body cut out, in mm
    cad_report.json  OCC volume and area, and the volume enclosed by the STL, so
                     the equal-volume constraint can be checked on what is meshed

Shapes:
    tet, pyr  the triangular or square base polygon at z_base, joined to the apex
              by planar triangles (faces -> shell -> solid)
    hem       a sphere of radius r scaled to height h along z (a spheroid), cut
              at its equator; apex up (N) or down (I)
"""
import json
import math
import os
import sys

import FreeCAD as App  # noqa: F401  (provided by freecadcmd)
import Part
from FreeCAD import Matrix, Vector


def pyramid(poly, z_base, apex):
    base = [Vector(x, y, z_base) for x, y in poly]
    A = Vector(*apex)
    faces = [Part.Face(Part.makePolygon([base[i], base[(i + 1) % len(base)], A, base[i]]))
             for i in range(len(base))]
    faces.append(Part.Face(Part.makePolygon(base + [base[0]])))
    solid = Part.Solid(Part.Shell(faces))
    if solid.Volume < 0:            # inconsistent face orientation -> flip
        solid.reverse()
    return solid


def half_spheroid(cx, cy, z_base, r, h, up):
    s = Part.makeSphere(r)
    m = Matrix()
    m.scale(1.0, 1.0, h / r)
    s = s.transformGeometry(m)
    pad = 1.0
    z0 = 0.0 if up else -(h + pad)
    box = Part.makeBox(2 * (r + pad), 2 * (r + pad), h + pad, Vector(-(r + pad), -(r + pad), z0))
    half = s.common(box).removeSplitter()
    half.translate(Vector(cx, cy, z_base))
    return half


def stl_volume_and_write(shape, path, name, tol_mm):
    """Tessellate, write an ASCII STL in metres, return the enclosed volume [mm3]."""
    pts, tris = shape.tessellate(tol_mm)
    vol6 = 0.0
    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for i, j, k in tris:
            a, b, c = pts[i], pts[j], pts[k]
            n = (b - a).cross(c - a)
            if n.Length < 1e-14:
                continue
            vol6 += a.dot(b.cross(c))
            n.normalize()
            f.write(f" facet normal {n.x:.6e} {n.y:.6e} {n.z:.6e}\n  outer loop\n")
            for v in (a, b, c):
                f.write(f"   vertex {v.x * 1e-3:.9e} {v.y * 1e-3:.9e} {v.z * 1e-3:.9e}\n")
            f.write("  endloop\n endfacet\n")
        f.write(f"endsolid {name}\n")
    return abs(vol6) / 6.0, len(tris)


def main():
    path = os.environ.get("ROOF_PCM_JSON")
    if not path:
        sys.exit("ROOF_PCM_JSON is not set")
    with open(path) as f:
        cfg = json.load(f)
    b, out = cfg["body"], cfg["out_dir"]
    os.makedirs(out, exist_ok=True)

    cx, cy = b["centre"]
    if b["shape"] in ("tet", "pyr"):
        body = pyramid(b["footprint"], b["z_base"], (cx, cy, b["z_apex"]))
    elif b["shape"] == "hem":
        body = half_spheroid(cx, cy, b["z_base"], b["radius"], b["height"],
                             up=b["z_apex"] > b["z_base"])
    else:
        sys.exit(f"unknown shape {b['shape']!r}")
    if not body.isValid() or body.Volume <= 0:
        sys.exit("invalid PCM body solid")

    s = cfg["slab"]
    slab = Part.makeBox(s["Lx"], s["Ly"], s["H"], Vector(-s["Lx"] / 2, -s["Ly"] / 2, 0))
    concrete = slab.cut(body)

    stl_vol, n_tri = stl_volume_and_write(body, os.path.join(out, "pcm_body.stl"), "pcm_body",
                                          cfg["stl_tolerance_mm"])
    Part.Compound([concrete, body]).exportStep(os.path.join(out, "pcm_body.step"))

    analytic = b["volume_mm3"]
    report = {
        "case_id": b["case_id"],
        "volume_ml_occ": body.Volume / 1000.0,
        "volume_ml_stl": stl_vol / 1000.0,
        "volume_ml_analytic": analytic / 1000.0,
        "occ_vs_analytic": body.Volume / analytic - 1.0,
        "stl_vs_analytic": stl_vol / analytic - 1.0,
        "surface_area_cm2_occ": body.Area / 100.0,
        "stl_triangles": n_tri,
        "concrete_volume_ml": concrete.Volume / 1000.0,
        "bounding_box_mm": [body.BoundBox.XMin, body.BoundBox.XMax, body.BoundBox.YMin,
                            body.BoundBox.YMax, body.BoundBox.ZMin, body.BoundBox.ZMax],
    }
    with open(os.path.join(out, "cad_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    if abs(report["occ_vs_analytic"]) > 1e-3:
        sys.exit(f"OCC volume differs from the analytic one by "
                 f"{100 * report['occ_vs_analytic']:.3f} %: {json.dumps(report)}")
    print("CAD OK", json.dumps(report))


main()
