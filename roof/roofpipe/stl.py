"""PCM body surface without FreeCAD (clusters rarely have it).

The same bodies as cad/build_roof_pcm_freecad.py: planar triangles for the
tetrahedron and pyramid (exact), a latitude-longitude triangulation of the
half-spheroid. Writes cad/pcm_body.stl in metres and cad/cad_report.json with the
volume the STL encloses. Standard library only.
"""
from __future__ import annotations

import json
import math
import os


def _pyramid(body):
    cx, cy = body["centre"]
    zb, za = body["z_base"], body["z_apex"]
    base = [(x, y, zb) for x, y in body["footprint"]]
    apex = (cx, cy, za)
    centre = (cx, cy, zb)
    tris = []
    n = len(base)
    for i in range(n):
        a, b = base[i], base[(i + 1) % n]
        tris.append((a, b, apex))
        tris.append((b, a, centre))          # base fan
    return tris


def _half_spheroid(body, n_lon=144, n_lat=48):
    cx, cy = body["centre"]
    r, h, zb = body["radius"], body["height"], body["z_base"]
    sgn = 1.0 if body["z_apex"] > zb else -1.0

    def p(lat, lon):                         # lat 0 = rim, pi/2 = apex
        return (cx + r * math.cos(lat) * math.cos(lon), cy + r * math.cos(lat) * math.sin(lon),
                zb + sgn * h * math.sin(lat))

    tris = []
    for j in range(n_lat):
        a0, a1 = 0.5 * math.pi * j / n_lat, 0.5 * math.pi * (j + 1) / n_lat
        for i in range(n_lon):
            b0, b1 = 2 * math.pi * i / n_lon, 2 * math.pi * (i + 1) / n_lon
            if j == n_lat - 1:
                tris.append((p(a0, b0), p(a0, b1), p(a1, 0.0)))
            else:
                tris.append((p(a0, b0), p(a0, b1), p(a1, b1)))
                tris.append((p(a0, b0), p(a1, b1), p(a1, b0)))
        # rim -> base disc
    for i in range(n_lon):
        b0, b1 = 2 * math.pi * i / n_lon, 2 * math.pi * (i + 1) / n_lon
        tris.append((p(0.0, b1), p(0.0, b0), (cx, cy, zb)))
    return tris


def _vol6(t):
    (ax, ay, az), (bx, by, bz), (cx, cy, cz) = t
    return ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)


def write(rdir: str) -> dict:
    """Write cad/pcm_body.stl (m) from cad/body.json (mm)."""
    cad = os.path.join(rdir, "cad")
    with open(os.path.join(cad, "body.json"), encoding="utf-8") as f:
        body = json.load(f)["body"]
    tris = _half_spheroid(body) if body["shape"] == "hem" else _pyramid(body)
    vol = sum(_vol6(t) for t in tris) / 6.0
    if vol < 0:                              # outward normals
        tris = [(a, c, b) for a, b, c in tris]
        vol = -vol
    with open(os.path.join(cad, "pcm_body.stl"), "w", encoding="utf-8", newline="\n") as f:
        f.write("solid pcm_body\n")
        for a, b, c in tris:
            u = [b[i] - a[i] for i in range(3)]
            v = [c[i] - a[i] for i in range(3)]
            n = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
            ln = math.sqrt(sum(x * x for x in n)) or 1.0
            f.write(f" facet normal {n[0] / ln:.6e} {n[1] / ln:.6e} {n[2] / ln:.6e}\n  outer loop\n")
            for q in (a, b, c):
                f.write(f"   vertex {q[0] * 1e-3:.9e} {q[1] * 1e-3:.9e} {q[2] * 1e-3:.9e}\n")
            f.write("  endloop\n endfacet\n")
        f.write("endsolid pcm_body\n")
    report = {"case_id": body["case_id"], "writer": "roofpipe/stl.py",
              "volume_ml_stl": vol / 1000.0, "volume_ml_analytic": body["volume_mm3"] / 1000.0,
              "stl_vs_analytic": vol / body["volume_mm3"] - 1.0, "stl_triangles": len(tris)}
    with open(os.path.join(cad, "cad_report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"  {body['case_id']}: STL {len(tris)} triangles, {vol / 1000.0:.2f} mL "
          f"({100 * report['stl_vs_analytic']:+.3f} % vs exact)")
    return report
