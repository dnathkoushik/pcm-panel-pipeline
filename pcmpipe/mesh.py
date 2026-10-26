"""Step 2 -- meshing with snappyHexMesh (multi-region, conformal interfaces).

Workflow inside runs/<case>/case:
  blockMesh                 room background mesh (paper Sec. 2.4 dimensions)
  surfaceFeatureExtract     sharp edges of the panel and cavity STLs
  snappyHexMesh -overwrite  castellate + snap; locationsInMesh -> cellZones
                            air | shell | pcm; snappy creates the interface faceZones
  splitMeshRegions -cellZones -overwrite
                            constant/{air,shell,pcm}/polyMesh with conformal
                            mappedWall patches air_to_shell, shell_to_pcm, ...
  refineWallLayer           optional wall-normal "inflation" layers on the fluid
                            side of both interfaces (mesh.wall_layer_passes, off
                            by default: see base_case.yaml). Interface faces are
                            not split, so the conformal coupling survives.
  checkMesh -region <r>
The shell is resolved with >= `shell_cells_across` (at least 2) cells through its
thickness. With a single cell, one cell-to-cell step can cross both surfaces and
snappy's locationsInMesh zone walk leaks (all cells end up in one zone).
"""
from __future__ import annotations

import glob
import math
import os
import re
import shutil

from . import foam, geometry
from .foam import vec, num

def thin_shell(case: dict) -> bool:
    """Screening mode: shell as a thermal resistance, not meshed (model.shell: thin)."""
    return case.get("model", {}).get("shell", "resolved") == "thin"


def regions(case: dict) -> tuple:
    return ("air", "pcm") if thin_shell(case) else ("air", "shell", "pcm")


def fluid_interfaces(case: dict) -> dict:
    if thin_shell(case):
        return {"air": "air_to_pcm", "pcm": "pcm_to_air"}
    return {"air": "air_to_shell", "pcm": "pcm_to_shell"}


def surfaces(case: dict) -> list:
    """(name in snappyHexMeshDict, STL written by FreeCAD); first = air-facing surface."""
    if thin_shell(case):
        return [("pcm_cavity", "pcm_cavity_thin.stl")]
    return [("panel_outer", "panel_outer.stl"), ("pcm_cavity", "pcm_cavity.stl")]


def levels(case: dict) -> dict:
    m = case["mesh"]["active"]
    t = case["panel"]["shell_thickness_mm"] * 1e-3
    base = m["base_cell"]
    n_across = max(2, m["shell_cells_across"])     # 1 breaks the zone walk, see above
    shell = m["surface_level_min"] if thin_shell(case) else max(
        m["surface_level_min"], math.ceil(math.log2(base * n_across / t) - 1e-9))
    return {"base": base, "shell": shell, "pcm": min(m["pcm_level"], shell),
            "air": m["air_distance_levels"], "plume": m["plume_level"],
            "shell_cell_mm": base / 2 ** shell * 1e3}


def write_dicts(case: dict, geo: geometry.PanelGeometry, cdir: str):
    room = case["physics"]["room"]
    lv = levels(case)
    sysd = os.path.join(cdir, "system")
    nx, ny, nz = (max(1, round(room[k] / lv["base"])) for k in ("Lx", "Ly", "Lz"))
    x0, x1, y0, y1 = -room["Lx"] / 2, room["Lx"] / 2, -room["Ly"] / 2, room["Ly"] / 2
    z0, z1 = -room["Lz"], 0.0

    foam.write(os.path.join(sysd, "blockMeshDict"), f"""
scale 1;
vertices
(
    ({num(x0)} {num(y0)} {num(z0)}) ({num(x1)} {num(y0)} {num(z0)})
    ({num(x1)} {num(y1)} {num(z0)}) ({num(x0)} {num(y1)} {num(z0)})
    ({num(x0)} {num(y0)} {num(z1)}) ({num(x1)} {num(y0)} {num(z1)})
    ({num(x1)} {num(y1)} {num(z1)}) ({num(x0)} {num(y1)} {num(z1)})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );
boundary
(
    floor   {{ type wall; faces ((0 3 2 1)); }}
    ceiling {{ type wall; faces ((4 5 6 7)); }}
    walls   {{ type wall; faces ((0 1 5 4) (1 2 6 5) (2 3 7 6) (3 0 4 7)); }}
);
""")

    feat = case["mesh"]["feature_angle"]
    foam.write(os.path.join(sysd, "surfaceFeatureExtractDict"), "\n".join(f"""
{s}.stl
{{
    extractionMethod    extractFromSurface;
    includedAngle       {feat};
    extractFromSurfaceCoeffs {{ includedAngle {feat}; }}
    subsetFeatures {{ nonManifoldEdges no; openEdges yes; }}
    writeObj            yes;
}}""" for s, _ in surfaces(case)))

    # plume refinement box below the panel (inverted plume, Sec. 3.1)
    xs = [p[0] * 1e-3 for p in geo.footprint]
    ys = [p[1] * 1e-3 for p in geo.footprint]
    pad = 0.05
    plume_min = (min(xs) - pad, min(ys) - pad, z0)
    plume_max = (max(xs) + pad, max(ys) + pad, geo.z_tip * 1e-3)
    air_levels = " ".join(f"({num(d)} {l})" for d, l in lv["air"])
    loc = geo.locations
    S = lv["shell"]
    names = [s for s, _ in surfaces(case)]
    # air-side distance refinement needs its own geometry entry; in thin mode the
    # only surface is the cavity, so a copy of it (pcm_cavity_air.stl) is used
    near = "pcm_cavity_air" if thin_shell(case) else "panel_outer"
    geom = "\n".join(f"    {s}.stl {{ type triSurfaceMesh; name {s}; }}"
                     for s in names + ([near] if near not in names else []))
    feats = "\n".join(f'        {{ file "{s}.eMesh"; level {S}; }}' for s in names)
    surfs = "\n".join(f"        {s} {{ level ({S} {S}); }}" for s in names)
    locs = "\n".join(f"        ({vec(loc[r])} {r})" for r in regions(case))

    foam.write(os.path.join(sysd, "snappyHexMeshDict"), f"""
castellatedMesh true;
snap            true;
addLayers       false;   // layers are added after splitMeshRegions (refineWallLayer)

geometry
{{
{geom}
    plumeBox {{ type searchableBox; min {vec(plume_min)}; max {vec(plume_max)}; }}
}}

castellatedMeshControls
{{
    maxLocalCells       5000000;
    maxGlobalCells      40000000;
    minRefinementCells  10;
    maxLoadUnbalance    0.10;
    nCellsBetweenLevels {case["mesh"]["n_cells_between_levels"]};

    features
    (
{feats}
    );

    refinementSurfaces
    {{
        // no faceZone here: unnamed surfaces stop the locationsInMesh walk;
        // named (faceZone) surfaces would be crossed by it (see ESI tutorial
        // heatTransfer/chtMultiRegionFoam/snappyMultiRegionHeater)
{surfs}
    }}

    resolveFeatureAngle 30;

    refinementRegions
    {{
        pcm_cavity  {{ mode inside;   levels ((1e15 {lv["pcm"]})); }}
        {near} {{ mode distance; levels ({air_levels}); }}
        plumeBox    {{ mode inside;   levels ((1e15 {lv["plume"]})); }}
    }}

    // flood fill from each point, bounded by the surfaces -> cellZones;
    // faceZones between different cellZones are created automatically
    locationsInMesh
    (
{locs}
    );
    allowFreeStandingZoneFaces false;
}}

snapControls
{{
    nSmoothPatch            3;
    tolerance               2.0;
    nSolveIter              100;
    nRelaxIter              5;
    nFeatureSnapIter        15;
    implicitFeatureSnap     false;
    explicitFeatureSnap     true;
    multiRegionFeatureSnap  true;
}}

addLayersControls
{{
    relativeSizes       true;
    layers              {{ }};
    expansionRatio      1.2;
    finalLayerThickness 0.5;
    minThickness        0.1;
    nGrow               0;
    featureAngle        60;
    nRelaxIter          3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals      3;
    nSmoothThickness    10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle  90;
    nBufferCellsNoExtrude 0;
    nLayerIter          50;
}}

meshQualityControls
{{
    #includeEtc "caseDicts/meshQualityDict"
    nSmoothScale    4;
    errorReduction  0.75;
}}

mergeTolerance 1e-6;
""")
    return lv


def write_mesh_system(cdir: str, n_procs: int):
    sysd = os.path.join(cdir, "system")
    foam.write(os.path.join(sysd, "controlDict"), """
application     snappyHexMesh;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1;
deltaT          1;
writeControl    timeStep;
writeInterval   1;
writeFormat     binary;
writePrecision  12;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;
""")
    foam.write(os.path.join(sysd, "fvSchemes"), """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    foam.write(os.path.join(sysd, "fvSolution"), "")
    foam.write(os.path.join(sysd, "decomposeParDict"),
               f"numberOfSubdomains {n_procs};\nmethod scotch;\n")


def _refine_wall_layers(cdir: str, region: str, patch: str, passes: int, frac: float, dry: bool):
    if passes <= 0:
        return
    tmp = os.path.join(cdir, f"_wallLayer_{region}")
    if not dry:
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(os.path.join(tmp, "constant"))
        shutil.copytree(os.path.join(cdir, "system"), os.path.join(tmp, "system"))
        shutil.copytree(os.path.join(cdir, "constant", region, "polyMesh"),
                        os.path.join(tmp, "constant", "polyMesh"))
    for i in range(passes):
        foam.run(["refineWallLayer", "-overwrite", f"({patch})", str(frac)], tmp,
                 log=f"refineWallLayer.{i}", dry=dry)
    if not dry:
        dst = os.path.join(cdir, "constant", region, "polyMesh")
        shutil.rmtree(dst)
        shutil.copytree(os.path.join(tmp, "constant", "polyMesh"), dst)
        for f in glob.glob(os.path.join(tmp, "log.*")):
            shutil.copy(f, os.path.join(cdir, f"log.{region}.{os.path.basename(f)[4:]}"))
        shutil.rmtree(tmp)


def cell_counts(cdir: str, region_names=("air", "shell", "pcm")) -> dict:
    counts = {}
    for r in region_names:
        log = os.path.join(cdir, f"log.checkMesh.{r}")
        if os.path.exists(log):
            m = re.search(r"^\s*cells:\s+(\d+)", open(log).read(), re.M)
            counts[r] = int(m.group(1)) if m else None
    return counts


def build(case: dict, n_procs: int = 1, dry: bool = False) -> dict:
    cdir = os.path.join(case["dir"], "case")
    geo = geometry.build_panel(case["panel"])
    if not dry:
        foam.check_openfoam()
        shutil.rmtree(cdir, ignore_errors=True)
    tri = os.path.join(cdir, "constant", "triSurface")
    os.makedirs(tri, exist_ok=True)
    if not dry:
        for name, src in surfaces(case):
            shutil.copy(os.path.join(case["dir"], "cad", src), os.path.join(tri, f"{name}.stl"))
        if thin_shell(case):
            shutil.copy(os.path.join(tri, "pcm_cavity.stl"), os.path.join(tri, "pcm_cavity_air.stl"))
    write_mesh_system(cdir, n_procs)
    lv = write_dicts(case, geo, cdir)
    if thin_shell(case):
        print(f"  thin shell ({case['panel']['shell_thickness_mm']:g} mm as a thermal resistance), "
              f"surface cells {lv['shell_cell_mm']:.3f} mm")
    else:
        print(f"  shell level {lv['shell']} -> {lv['shell_cell_mm']:.3f} mm cells "
              f"({case['panel']['shell_thickness_mm'] / lv['shell_cell_mm']:.1f} across the shell)")

    foam.run(["blockMesh"], cdir, dry=dry)
    foam.run(["surfaceFeatureExtract"], cdir, dry=dry)
    if n_procs > 1:
        foam.run(["decomposePar", "-force"], cdir, log="decomposePar.mesh", dry=dry)
        foam.run([*foam.tool(case, "mpirun"), "-np", str(n_procs), "snappyHexMesh",
                  "-parallel", "-overwrite"], cdir, log="snappyHexMesh", dry=dry)
        foam.run(["reconstructParMesh", "-constant"], cdir, dry=dry)
        if not dry:
            for d in glob.glob(os.path.join(cdir, "processor*")):
                shutil.rmtree(d)
    else:
        foam.run(["snappyHexMesh", "-overwrite"], cdir, dry=dry)
    foam.run(["splitMeshRegions", "-cellZones", "-overwrite"], cdir, dry=dry)

    if not dry:
        # hexRef8 bookkeeping is invalid after refineWallLayer; not needed downstream
        for r in regions(case):
            for f in ("cellLevel", "pointLevel", "level0Edge", "surfaceIndex", "refinementHistory"):
                p = os.path.join(cdir, "constant", r, "polyMesh", f)
                if os.path.exists(p):
                    os.remove(p)
        shutil.rmtree(os.path.join(cdir, "constant", "polyMesh"), ignore_errors=True)

    m = case["mesh"]
    passes = m["active"].get("wall_layer_passes") or {}
    for region, patch in fluid_interfaces(case).items():
        _refine_wall_layers(cdir, region, patch, passes.get(region, 0),
                            m["wall_layer_edge_fraction"], dry)
    for r in regions(case):
        foam.run(["checkMesh", "-region", r, "-constant"], cdir, log=f"checkMesh.{r}", dry=dry)

    counts = cell_counts(cdir, regions(case))
    if counts and all(counts.values()):
        print(f"  cells: {counts}  total {sum(counts.values()):,}")
    return counts
