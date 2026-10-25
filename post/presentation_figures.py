"""Presentation figures that need no CFD run (mid-evaluation material).

    python post/presentation_figures.py        ->  runs/_figures/presentation/

  design_catalogue.png        3D views of candidate panels (shape x design variant)
  geometry_validation.png/md  our geometry vs. the numbers published in the paper
  latent_heat_verification.png
                              enthalpy-porosity source (as implemented in OpenFOAM)
                              vs. a fully converged enthalpy solution, 1D slab
"""
from __future__ import annotations

import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))
from pcmpipe import config  # noqa: E402
from pcmpipe.geometry import build_panel, polygon_area  # noqa: E402

OUT = os.path.join(config.RUNS_DIR, "_figures", "presentation")
BASE = config.base_case()["panel"]
FACE = {"pyramid": "#8fa8e8", "tetrahedron": "#ef9a9a", "hemisphere": "#a5d6a7"}


def panel(**kw):
    d = dict(BASE)
    d.update(kw)
    return build_panel(d)


def draw_panel(ax, geo, half_span, depth):
    """Envelope of the panel (mm), hanging below the ceiling z = 0, on a fixed scale."""
    hb, col = geo.slab_thickness, FACE[geo.shape]
    fp = geo.footprint
    # matplotlib does not depth-sort separate 3D artists: an opaque slab underside
    # would hide the (surface-plotted) hemispheres, so draw it translucent there
    under_alpha = 0.25 if geo.shape == "hemisphere" else 0.95
    ax.add_collection3d(Poly3DCollection([[(x, y, -hb) for x, y in fp]], facecolor=col,
                                         edgecolor="k", linewidths=0.3, alpha=under_alpha))
    polys = []
    for (x1, y1), (x2, y2) in zip(fp, fp[1:] + fp[:1]):            # slab sides
        polys.append([(x1, y1, 0), (x2, y2, 0), (x2, y2, -hb), (x1, y1, -hb)])
    for p in geo.protrusions:
        cx, cy = p.centre
        if p.kind == "spheroid":
            u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi / 2, 15))
            ax.plot_surface(cx + p.radius * np.cos(u) * np.sin(v),
                            cy + p.radius * np.sin(u) * np.sin(v),
                            -hb - p.height * np.cos(v), color=col, linewidth=0,
                            antialiased=True, shade=True)
        else:
            apex = (cx, cy, -hb - p.height)
            b = p.base_polygon
            for (x1, y1), (x2, y2) in zip(b, b[1:] + b[:1]):
                polys.append([(x1, y1, -hb), (x2, y2, -hb), apex])
    ax.add_collection3d(Poly3DCollection(polys, facecolor=col, edgecolor="k",
                                         linewidths=0.3, alpha=0.95))
    ax.set_xlim(-half_span, half_span)
    ax.set_ylim(-half_span, half_span)
    ax.set_zlim(-depth, 0)
    ax.set_box_aspect((2 * half_span, 2 * half_span, depth))     # true proportions
    ax.view_init(elev=-22, azim=-58)
    ax.set_axis_off()


def extent(geo):
    xs = [abs(x) for x, _ in geo.footprint] + [abs(y) for _, y in geo.footprint]
    return max(xs), -geo.z_tip


def design_catalogue():
    variants = [("Reference\n250 mL, 50 %", {}),
                ("Protrusion ratio 70 %", {"protrusion_ratio": 0.7}),
                ("Gap 5 mm", {"gap_mm": 5.0}),
                ("Volume 750 mL", {"volume_ml": 750.0})]
    shapes = ["pyramid", "tetrahedron", "hemisphere"]
    fig = plt.figure(figsize=(4.0 * len(variants), 3.2 * len(shapes)))
    for i, s in enumerate(shapes):
        geos = [panel(shape=s, **kw) for _, kw in variants]
        half_span = max(extent(g)[0] for g in geos) * 1.02          # same scale per row
        depth = max(extent(g)[1] for g in geos) * 1.02
        for j, ((title, kw), geo) in enumerate(zip(variants, geos)):
            ax = fig.add_subplot(len(shapes), len(variants), i * len(variants) + j + 1,
                                 projection="3d")
            draw_panel(ax, geo, half_span, depth)
            ax.set_title(f"{s.capitalize()} | {title.splitlines()[0]}"
                         + (f"\n{title.splitlines()[1]}" if "\n" in title else "")
                         + f"\nA = {geo.exposed_area_mm2 / 100:.0f} cm$^2$, "
                         f"A/V = {geo.exposed_area_mm2 / 100 / geo.volume_ml:.2f} cm$^2$/mL",
                         fontsize=9)
    fig.suptitle("Candidate PCM ceiling panels (envelope, 1 mm ABS shell not drawn)", fontsize=13)
    fig.tight_layout()
    save(fig, "design_catalogue")


def geometry_validation():
    rows = []
    pyr, tet = panel(shape="pyramid"), panel(shape="tetrahedron")
    rows.append(("Effective area, pyramid (Table 3)", "232 cm²", f"{pyr.exposed_area_mm2 / 100:.0f} cm²"))
    rows.append(("Effective area, tetrahedron (Table 3)", "282 cm²", f"{tet.exposed_area_mm2 / 100:.0f} cm²"))
    a = math.dist(tet.footprint[0], tet.footprint[1]) / 10
    b = math.dist(tet.footprint[0], tet.footprint[3]) / 10
    rows.append(("Tetrahedral footprint (Fig. 4)", "12.8 × 6.4 cm", f"{a:.1f} × {b:.1f} cm"))
    for r, p_ref, t_ref in ((0.3, 29.2, 40.3), (0.5, 43.0, 54.7), (0.7, 52.6, 63.2)):
        rows.append((f"Edge angle at {r:.0%} protrusion (Sec. 4.5)",
                     f"{p_ref}° / {t_ref}°",
                     f"{panel(protrusion_ratio=r).edge_angle_deg:.1f}° / "
                     f"{panel(shape='tetrahedron', protrusion_ratio=r).edge_angle_deg:.1f}°"))
    for V, ref in ((250, 1.13), (500, 0.895), (750, 0.782)):
        g = panel(shape="tetrahedron", volume_ml=V)
        rows.append((f"Tetrahedron area/volume, {V} mL (Sec. 4.7)", f"{ref} cm²/mL",
                     f"{g.exposed_area_mm2 / 100 / V:.3f} cm²/mL"))
    rows.append(("Surface/volume, pyramid / tetra (Table 3)", "1.33 / 1.41 cm⁻¹",
                 f"{(pyr.exposed_area_mm2 + polygon_area(pyr.footprint)) / 100 / 250:.2f} / "
                 f"{(tet.exposed_area_mm2 + polygon_area(tet.footprint)) / 100 / 250:.2f} cm⁻¹"))

    fig, ax = plt.subplots(figsize=(11, 0.32 * len(rows) + 0.4))
    ax.axis("off")
    tab = ax.table(cellText=[list(r) for r in rows], colLabels=["Quantity", "Paper", "This work"],
                   loc="center", cellLoc="left", colWidths=[0.5, 0.22, 0.28])
    tab.auto_set_font_size(False)
    tab.set_fontsize(10)
    tab.scale(1, 1.4)
    for (r, c), cell in tab.get_celld().items():
        if r == 0:
            cell.set_facecolor("#dfe6f5")
            cell.set_text_props(weight="bold")
    ax.set_title("Parametric geometry reproduces the paper's panels", fontsize=12)
    save(fig, "geometry_validation")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "geometry_validation.md"), "w", encoding="utf-8") as f:
        f.write("| Quantity | Paper | This work |\n|---|---|---|\n")
        f.writelines(f"| {a} | {b} | {c} |\n" for a, b, c in rows)


def latent_heat_verification():
    from test_latent_source_1d import run_fvoption, run_reference
    ref = run_reference()
    runs = {n: run_fvoption(n_outer=n) for n in (1, 2)}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.plot(ref[:, 0] / 60, ref[:, 1], "k-", lw=3, alpha=0.35, label="converged enthalpy method")
    for n, st in ((1, "--"), (2, ":")):
        ax1.plot(runs[n][:, 0] / 60, runs[n][:, 1], st, lw=1.6,
                 label=f"OpenFOAM source term, {n} outer corrector{'s' if n > 1 else ''}")
    ax1.set_xlabel("Time (min)")
    ax1.set_ylabel("Liquid fraction")
    ax1.set_title("1D PCM slab heated from one side (Eqs. 1, 15, 18-21)", fontsize=10)
    ax1.legend(fontsize=8, frameon=False)
    for n in (1, 2):
        ax2.semilogy(runs[n][:, 0] / 60, np.abs(runs[n][:, 1] - ref[:, 1]) + 1e-12,
                     label=f"{n} outer corrector{'s' if n > 1 else ''}")
    ax2.set_xlabel("Time (min)")
    ax2.set_ylabel("|liquid fraction error|")
    ax2.set_title("Error vs. converged solution", fontsize=10)
    ax2.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    save(fig, "latent_heat_verification")


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {os.path.join(OUT, name)}.png")


if __name__ == "__main__":
    design_catalogue()
    geometry_validation()
    latent_heat_verification()
