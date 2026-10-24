"""Parametric dimension solver for PCM ceiling panels (pure Python, no deps).

Coordinate system: the ceiling is the plane z = 0, the panel hangs downwards.
The panel envelope is a flat "slab" (z in [-hb, 0]) over a convex footprint,
with an array of protrusions hanging below z = -hb.

Rules reproduced from Yang et al. (2023):
  * Reference design (protrusion_ratio 0.5, no gap) is fixed by the lateral edge
    angle: 43.0 deg for pyramids, 54.7356 deg (regular) for tetrahedra (Sec. 4.5).
  * Changing the protrusion ratio keeps the footprint and protrusion bases and
    changes the protrusion height (Sec. 4.5).  This reproduces the paper's
    29.2/43.0/52.6 deg and 40.3/54.7/63.2 deg edge angles exactly.
  * Changing the volume scales the panel self-similarly (Sec. 4.7).
  * The paper's 250 mL tetrahedral footprint (12.8 cm x 6.4 cm) and effective
    areas (232 / 282 cm^2, Table 3) are recovered (see tests/).
Assumptions (paper silent):
  * Gaps shrink each protrusion base by offsetting its edges inward by g/2 at a
    fixed footprint; protrusion height grows so volume and ratio are preserved.
  * Hemispheres: R_ref = fill * pitch / 2; other ratios/gaps give semi-spheroids.
All lengths are in millimetres, volumes in mm^3 unless stated otherwise.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, asdict

SQRT3 = math.sqrt(3.0)


@dataclass
class Protrusion:
    """One protrusion hanging below the slab."""
    kind: str                 # "pyramid" | "tetrahedron" | "spheroid"
    centre: tuple             # (x, y) of the apex / axis
    base_polygon: list        # [(x, y), ...] at z = -hb (empty for spheroids)
    radius: float = 0.0       # spheroid equatorial radius
    height: float = 0.0       # protrusion height below the slab


@dataclass
class PanelGeometry:
    shape: str
    volume_ml: float
    protrusion_ratio: float
    shell_thickness: float
    gap: float
    footprint: list                      # convex polygon [(x, y), ...], CCW
    slab_thickness: float                # hb
    protrusion_height: float             # hp (or spheroid semi-axis c)
    face_tilt_deg: float                 # protrusion face angle vs horizontal
    edge_angle_deg: float                # lateral edge angle vs base
    protrusions: list = field(default_factory=list)
    has_flat_bottom: bool = False        # slab underside exposed between protrusions
    exposed_area_mm2: float = 0.0        # envelope area excluding the ceiling face
    envelope_volume_mm3: float = 0.0
    z_tip: float = 0.0
    probes: dict = field(default_factory=dict)     # thermocouple positions (m)
    locations: dict = field(default_factory=dict)  # snappy locationsInMesh (m)
    slice_plane: dict = field(default_factory=dict)  # for contour plots (m)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["protrusions"] = [asdict(p) for p in self.protrusions]
        return d


# ----------------------------------------------------------------------------
# polygon helpers
# ----------------------------------------------------------------------------
def polygon_area(poly):
    a = 0.0
    for (x1, y1), (x2, y2) in zip(poly, poly[1:] + poly[:1]):
        a += x1 * y2 - x2 * y1
    return 0.5 * a


def polygon_perimeter(poly):
    return sum(math.dist(p, q) for p, q in zip(poly, poly[1:] + poly[:1]))


def inset_polygon(poly, d):
    """Offset every edge of a convex CCW polygon inward by d (sharp corners)."""
    n = len(poly)
    lines = []
    for i in range(n):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
        ex, ey = x2 - x1, y2 - y1
        L = math.hypot(ex, ey)
        nx, ny = -ey / L, ex / L          # inward normal for CCW polygons
        lines.append(((x1 + d * nx, y1 + d * ny), (ex, ey)))
    out = []
    for i in range(n):
        (p, r), (q, s) = lines[i - 1], lines[i]
        den = r[0] * s[1] - r[1] * s[0]
        t = ((q[0] - p[0]) * s[1] - (q[1] - p[1]) * s[0]) / den
        out.append((p[0] + t * r[0], p[1] + t * r[1]))
    return out


def _shift(poly, dx, dy):
    return [(x + dx, y + dy) for x, y in poly]


def half_spheroid_area(a, c):
    """Lateral area of half a spheroid with equatorial radius a, semi-axis c."""
    if abs(c - a) < 1e-12 * a:
        return 2.0 * math.pi * a * a
    if c < a:  # oblate
        e = math.sqrt(1.0 - (c / a) ** 2)
        full = 2 * math.pi * a * a * (1 + (1 - e * e) / e * math.atanh(e))
    else:      # prolate
        e = math.sqrt(1.0 - (a / c) ** 2)
        full = 2 * math.pi * a * a * (1 + c / (a * e) * math.asin(e))
    return 0.5 * full


# ----------------------------------------------------------------------------
# shape builders
# ----------------------------------------------------------------------------
def _pyramid(V, r, g, nx, ny, edge_deg):
    N = nx * ny
    # reference pitch: N * s^2 * hp / 3 = 0.5 V with hp = (s / sqrt2) tan(edge)
    s = (3.0 * 0.5 * V * math.sqrt(2.0) / (N * math.tan(math.radians(edge_deg)))) ** (1 / 3)
    Lx, Ly = nx * s, ny * s
    footprint = [(-Lx / 2, -Ly / 2), (Lx / 2, -Ly / 2), (Lx / 2, Ly / 2), (-Lx / 2, Ly / 2)]
    b = s - g
    if b <= 0:
        raise ValueError(f"gap {g} mm too large for pyramid pitch {s:.2f} mm")
    hp = 3.0 * r * V / (N * b * b)
    prots = []
    for i in range(nx):
        for j in range(ny):
            cx, cy = -Lx / 2 + (i + 0.5) * s, -Ly / 2 + (j + 0.5) * s
            sq = [(-b / 2, -b / 2), (b / 2, -b / 2), (b / 2, b / 2), (-b / 2, b / 2)]
            prots.append(Protrusion("pyramid", (cx, cy), _shift(sq, cx, cy), height=hp))
    slant = math.hypot(hp, b / 2)
    lateral = N * 4 * 0.5 * b * slant
    tilt = math.degrees(math.atan2(hp, b / 2))
    edge = math.degrees(math.atan2(hp, b / math.sqrt(2.0)))
    return footprint, prots, hp, tilt, edge, lateral, N * b * b, g > 0


def _tetrahedron(V, r, g, n_pairs, edge_deg):
    N = 2 * n_pairs
    # reference: each right triangular pyramid V_i = a^3 tan(edge) / 12
    a = (12.0 * 0.5 * V / (N * math.tan(math.radians(edge_deg)))) ** (1 / 3)
    h3 = a * SQRT3 / 2
    footprint = [(0, 0), (n_pairs * a, 0), (n_pairs * a + a / 2, h3), (a / 2, h3)]
    cx0 = sum(p[0] for p in footprint) / 4
    cy0 = sum(p[1] for p in footprint) / 4
    footprint = _shift(footprint, -cx0, -cy0)
    a2 = a - SQRT3 * g                     # edges offset inward by g/2
    if a2 <= 0:
        raise ValueError(f"gap {g} mm too large for tetrahedron edge {a:.2f} mm")
    hp = 4.0 * SQRT3 * (r * V / N) / (a2 * a2)
    r_in = a2 / (2 * SQRT3)                # inradius of the shrunk triangle
    prots = []
    for i in range(n_pairs):
        up = [(i * a, 0), ((i + 1) * a, 0), (i * a + a / 2, h3)]
        dn = [((i + 1) * a, 0), (i * a + 1.5 * a, h3), (i * a + a / 2, h3)]
        for tri in (up, dn):
            cx = sum(p[0] for p in tri) / 3
            cy = sum(p[1] for p in tri) / 3
            k = a2 / a
            shrunk = [(cx + k * (x - cx), cy + k * (y - cy)) for x, y in tri]
            prots.append(Protrusion("tetrahedron", (cx - cx0, cy - cy0),
                                    _shift(shrunk, -cx0, -cy0), height=hp))
    slant = math.hypot(hp, r_in)
    lateral = N * 3 * 0.5 * a2 * slant
    tilt = math.degrees(math.atan2(hp, r_in))
    edge = math.degrees(math.atan2(hp, a2 / SQRT3))
    return footprint, prots, hp, tilt, edge, lateral, N * SQRT3 / 4 * a2 * a2, g > 0


def _hemisphere(V, r, g, nx, ny, fill):
    N = nx * ny
    # reference: N * (2/3) pi R^3 = 0.5 V with R = fill * s / 2
    s = (0.5 * V / (N * (2.0 / 3.0) * math.pi * (fill / 2) ** 3)) ** (1 / 3)
    Lx, Ly = nx * s, ny * s
    footprint = [(-Lx / 2, -Ly / 2), (Lx / 2, -Ly / 2), (Lx / 2, Ly / 2), (-Lx / 2, Ly / 2)]
    R = fill * s / 2 - g / 2
    if R <= 0:
        raise ValueError(f"gap {g} mm too large for hemisphere pitch {s:.2f} mm")
    c = 3.0 * r * V / (2.0 * math.pi * N * R * R)
    prots = []
    for i in range(nx):
        for j in range(ny):
            cx, cy = -Lx / 2 + (i + 0.5) * s, -Ly / 2 + (j + 0.5) * s
            prots.append(Protrusion("spheroid", (cx, cy), [], radius=R, height=c))
    lateral = N * half_spheroid_area(R, c)
    tilt = 90.0                                  # vertical tangent at the equator
    edge = math.degrees(math.atan2(c, R))
    return footprint, prots, c, tilt, edge, lateral, N * math.pi * R * R, True


# ----------------------------------------------------------------------------
# public API
# ----------------------------------------------------------------------------
def build_panel(panel: dict) -> PanelGeometry:
    """Resolve all dimensions from the `panel` section of a case config."""
    shape = panel["shape"]
    V = float(panel["volume_ml"]) * 1000.0           # mL -> mm^3
    r = float(panel["protrusion_ratio"])
    t = float(panel["shell_thickness_mm"])
    g = float(panel.get("gap_mm", 0.0))
    if not 0.0 < r < 1.0:
        raise ValueError("protrusion_ratio must be in (0, 1)")

    if shape == "pyramid":
        fp, prots, hp, tilt, edge, lateral, covered, flat = _pyramid(
            V, r, g, int(panel["n_x"]), int(panel["n_y"]), float(panel["pyramid_edge_angle_deg"]))
    elif shape == "tetrahedron":
        fp, prots, hp, tilt, edge, lateral, covered, flat = _tetrahedron(
            V, r, g, max(1, int(panel["n_x"])), float(panel["tetrahedron_edge_angle_deg"]))
    elif shape == "hemisphere":
        fp, prots, hp, tilt, edge, lateral, covered, flat = _hemisphere(
            V, r, g, int(panel["n_x"]), int(panel["n_y"]), float(panel["hemisphere_fill"]))
    else:
        raise ValueError(f"unknown shape '{shape}'")

    A_fp = polygon_area(fp)
    hb = (1.0 - r) * V / A_fp
    if hb <= 2 * t:
        raise ValueError(f"slab thickness {hb:.2f} mm too thin for a {t} mm shell")
    if hp <= 2 * t:
        raise ValueError(f"protrusion height {hp:.2f} mm too small for a {t} mm shell")

    exposed = polygon_perimeter(fp) * hb + lateral + (A_fp - covered)
    geo = PanelGeometry(
        shape=shape, volume_ml=V / 1000.0, protrusion_ratio=r, shell_thickness=t, gap=g,
        footprint=fp, slab_thickness=hb, protrusion_height=hp, face_tilt_deg=tilt,
        edge_angle_deg=edge, protrusions=prots, has_flat_bottom=flat,
        exposed_area_mm2=exposed, envelope_volume_mm3=A_fp * hb + r * V,
        z_tip=-(hb + hp),
    )
    _add_points(geo)
    return geo


def inner_cone_lift(geo: PanelGeometry) -> float:
    """Upward shift of a planar-faced protrusion after an inward offset by t."""
    return geo.shell_thickness / math.cos(math.radians(geo.face_tilt_deg))


def _add_points(geo: PanelGeometry):
    """Thermocouples (Sec. 2.4), snappy seed points and the contour plane (m)."""
    mm = 1e-3
    p0 = geo.protrusions[0]
    cx, cy = p0.centre
    hb, hp, t = geo.slab_thickness, geo.protrusion_height, geo.shell_thickness
    z_tip = -(hb + hp)
    # top / mid / bottom at 1, 2/3, 1/3 of the protrusion height above its tip
    geo.probes = {name: [cx * mm, cy * mm, (z_tip + f * hp) * mm]
                  for name, f in (("top", 1.0), ("mid", 2 / 3), ("bot", 1 / 3))}
    eps = 1.234e-5   # keep seed points off cell faces
    xs, ys = (sum(p[0] for p in geo.footprint) / len(geo.footprint),
              sum(p[1] for p in geo.footprint) / len(geo.footprint))
    geo.locations = {
        "air": [xs * mm + eps, ys * mm + eps, -1.0 + eps],
        "shell": [xs * mm + eps, ys * mm + eps, -0.3 * t * mm + eps],
        "pcm": [cx * mm + eps, cy * mm + eps, (z_tip + 0.5 * hp) * mm + eps],
    }
    # vertical plane through the first protrusion axis; for square arrays it
    # cuts two apexes like the paper's Fig. 9a, for tetrahedra one apex (Fig. 9b)
    geo.slice_plane = {"origin": [cx * mm, cy * mm, -0.05], "normal": [0.0, 1.0, 0.0]}


def summary(geo: PanelGeometry) -> str:
    return (f"{geo.shape}: V={geo.volume_ml:.0f} mL  r={geo.protrusion_ratio:.2f}  "
            f"t={geo.shell_thickness:g} mm  gap={geo.gap:g} mm | footprint "
            f"{polygon_area(geo.footprint) / 100:.1f} cm2, slab hb={geo.slab_thickness:.2f} mm, "
            f"hp={geo.protrusion_height:.2f} mm, edge={geo.edge_angle_deg:.1f} deg, "
            f"exposed area={geo.exposed_area_mm2 / 100:.1f} cm2")
