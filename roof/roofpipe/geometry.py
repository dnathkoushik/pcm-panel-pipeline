"""Phase 2 geometry: one PCM body embedded in the slab unit cell (config: pcm_study).

Three shapes, two orientations, five heights; the slab (Lx x Ly x H) never changes.

  tet  tetrahedron: triangular pyramid on an equilateral base, apex above the
       base centroid
  pyr  pyramid: square base, apex above its centre
  hem  hemisphere: half of a spheroid on a circular base of radius r. It is a
       true hemisphere when the height equals r, a half-ellipsoid otherwise

Equal volume at every height. All three have V = A h / 3 when the tetrahedron and
pyramid bases share the area A = 2 pi r^2, twice the hemisphere's base circle
(half-spheroid: V = (2/3) pi r^2 h). Changing only the height changes the volume,
and the volume is the same for every shape at a given height.

Orientation
  N  normal:   flat base on the ceiling side (z = cover), apex towards the roof
  I  inverted: flat base on the roof side (z = H - cover), apex towards the ceiling

Lengths in metres. Standard library only.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field

SHAPES = {"tet": "tetrahedron", "pyr": "pyramid", "hem": "hemisphere"}
ORIENTATIONS = {"N": "normal: base on the ceiling side, apex towards the roof",
                "I": "inverted: base on the roof side, apex towards the ceiling"}


@dataclass
class PcmBody:
    case_id: str
    shape: str                      # tet | pyr | hem
    orientation: str                # N | I
    height: float                   # m, base plane to apex
    base_area: float                # m2
    volume: float                   # m3, analytic
    z_base: float                   # m, base plane
    z_apex: float                   # m
    centre: tuple                   # (x, y) of the base centroid = apex position
    footprint: list = field(default_factory=list)   # tet/pyr base polygon [(x, y), ...]
    radius: float = 0.0             # hem: base radius
    lateral_area: float = 0.0       # m2, sloped / curved surface
    surface_area: float = 0.0       # m2, lateral + base
    margin_xy: float = 0.0          # m, smallest concrete width between body and cell side

    def to_dict(self) -> dict:
        return asdict(self)


def case_id(shape: str, orientation: str, height_m: float) -> str:
    return f"{shape}_{orientation}_h{round(height_m * 1000):03d}"


def _spheroid_cap_area(a: float, c: float) -> float:
    """Curved area of half a spheroid with equatorial radius a and polar half-axis c."""
    if abs(a - c) < 1e-12:
        return 2.0 * math.pi * a * a
    if c < a:                                        # oblate
        e = math.sqrt(1.0 - c * c / (a * a))
        return math.pi * a * a + math.pi * c * c / e * math.atanh(e)
    e = math.sqrt(1.0 - a * a / (c * c))             # prolate
    return math.pi * a * a + math.pi * a * c / e * math.asin(e)


def _pyramid_lateral(poly, apex_xy, h) -> float:
    """Area of the triangular faces from the base polygon to the apex at height h."""
    ax, ay = apex_xy
    area = 0.0
    n = len(poly)
    for i in range(n):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
        u = (x2 - x1, y2 - y1, 0.0)
        v = (ax - x1, ay - y1, h)
        cx = u[1] * v[2] - u[2] * v[1]
        cy = u[2] * v[0] - u[0] * v[2]
        cz = u[0] * v[1] - u[1] * v[0]
        area += 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)
    return area


def body(cfg: dict, shape: str, orientation: str, height: float) -> PcmBody:
    """The PCM body of one Phase 2 case, checked against the slab and the cover."""
    if shape not in SHAPES:
        raise ValueError(f"shape {shape!r} not in {list(SHAPES)}")
    if orientation not in ORIENTATIONS:
        raise ValueError(f"orientation {orientation!r} not in {list(ORIENTATIONS)}")
    s, ps = cfg["slab"], cfg["pcm_study"]
    Lx, Ly, H = float(s["Lx"]), float(s["Ly"]), float(s["H"])
    r = float(ps["hemisphere_radius_mm"]) * 1e-3
    cover = float(ps["cover_mm"]) * 1e-3
    A = 2.0 * math.pi * r * r                       # tet / pyr base area

    if orientation == "N":
        z_base, z_apex = cover, cover + height
    else:
        z_base, z_apex = H - cover, H - cover - height
    if min(z_base, z_apex) < cover - 1e-9 or max(z_base, z_apex) > H - cover + 1e-9:
        raise ValueError(f"{shape} {orientation} h = {height * 1e3:g} mm does not fit between "
                         f"the {cover * 1e3:g} mm covers of the {H * 1e3:g} mm slab")

    footprint, radius = [], 0.0
    if shape == "pyr":
        a = math.sqrt(A)
        footprint = [(-a / 2, -a / 2), (a / 2, -a / 2), (a / 2, a / 2), (-a / 2, a / 2)]
        centre = (0.0, 0.0)
        lateral = _pyramid_lateral(footprint, centre, height)
        half_x, half_y = a / 2, a / 2
    elif shape == "tet":
        side = math.sqrt(4.0 * A / math.sqrt(3.0))
        tri_h = side * math.sqrt(3.0) / 2.0
        # bounding box centred on the cell axis; apex above the triangle's centroid
        footprint = [(-side / 2, -tri_h / 2), (side / 2, -tri_h / 2), (0.0, tri_h / 2)]
        centre = (0.0, -tri_h / 2 + tri_h / 3)
        lateral = _pyramid_lateral(footprint, centre, height)
        half_x, half_y = side / 2, tri_h / 2
    else:
        radius = r
        centre = (0.0, 0.0)
        lateral = _spheroid_cap_area(r, height) - math.pi * r * r
        half_x = half_y = r
    margin = min(Lx / 2 - half_x, Ly / 2 - half_y)
    min_margin = float(ps.get("min_side_margin_mm", 10.0)) * 1e-3
    if margin < min_margin:
        raise ValueError(f"{shape}: the base leaves only {margin * 1e3:.1f} mm of concrete to "
                         f"the cell side (pcm_study.min_side_margin_mm = {min_margin * 1e3:g})")
    return PcmBody(
        case_id=case_id(shape, orientation, height), shape=shape, orientation=orientation,
        height=height, base_area=A if shape != "hem" else math.pi * r * r,
        volume=A * height / 3.0, z_base=z_base, z_apex=z_apex, centre=centre,
        footprint=footprint, radius=radius, lateral_area=lateral,
        surface_area=lateral + (A if shape != "hem" else math.pi * r * r), margin_xy=margin)


def study_cases(cfg: dict) -> list:
    """Every Phase 2 case: shapes x orientations x heights (config order)."""
    ps = cfg["pcm_study"]
    return [body(cfg, sh, o, float(h) * 1e-3)
            for sh in ps["shapes"] for o in ps["orientations"] for h in ps["heights_mm"]]


def locations(b: PcmBody, Lx: float, Ly: float, H: float) -> dict:
    """A point inside each region, for snappyHexMesh locationsInMesh. The PCM point
    sits on the axis a quarter of the height from the base; the concrete point near
    a cell corner at mid-depth, clear of every body (min_side_margin_mm >= 10).
    The small odd offsets keep both off the faces of the background mesh."""
    cx, cy = b.centre
    z_in = b.z_base + 0.25 * (b.z_apex - b.z_base)
    return {"pcm": (cx + 1.37e-4, cy + 1.13e-4, z_in),
            "concrete": (0.5 * Lx - 0.0129, 0.5 * Ly - 0.0127, 0.5 * H + 1.1e-4)}
