"""Phase 3, ParaView part: what needs the 3D fields, from the few snapshots of the
analysed day (written every write_interval_h), into small files.

    pvpython roof/post/pv_extract.py <run_dir> [--purge-fields]
    (headless servers: xvfb-run -a pvpython ...)

Writes <run_dir>/results/pv/:
    slice_<hh>h.csv          vertical mid-plane through the PCM axis (x-z): T [C] and, in
                             the PCM, liquidFraction, at cell centres
    contour_T_<hh>h.png      T on that plane, the PCM outline and its melting front
                             (liquid fraction 0.5); fixed 25-60 C scale for every case
    surface_metrics.csv      per snapshot: area-averaged T of the roof and ceiling faces and
                             their heat flux from the written wallHeatFlux field; a check on
                             the in-solver time series (every 300 s), which
                             `python3 -m roofpipe check` already turned into
                             results/timeseries_lastday.csv
--purge-fields then deletes the field folders (keeps 0 and the latest, for a restart), so
a finished case is a few MB.

Needs only the standard library and paraview.simple. Reads case_meta.json.
"""
import argparse
import csv
import json
import os
import shutil

from paraview import servermanager as sm
from paraview.simple import (Calculator, CellDatatoPointData, ColorBy, Contour, CreateView, GetColorTransferFunction,
                             GetScalarBar, IntegrateVariables, MergeBlocks, OpenFOAMReader,
                             Render, SaveData, SaveScreenshot, Show, Slice, Text, UpdatePipeline)

DAY = 86400.0
T_RANGE = (25.0, 55.0)      # C, the same colour scale for every case and time


def blocks(reader, region, kind):
    """Reader block names of one region: kind 'internalMesh' or a patch name."""
    avail = list(reader.MeshRegions.Available)
    out = [r for r in avail if region in r.split("/") and r.split("/")[-1] == kind]
    if not out:
        raise RuntimeError(f"{region}/{kind} not in {avail}")
    return out


def reader(cfoam, names, arrays):
    rd = OpenFOAMReader(FileName=cfoam)
    rd.MeshRegions = names
    rd.CellArrays = arrays
    rd.UpdatePipelineInformation()
    return rd


def plane_slice(src, origin):
    sl = Slice(Input=MergeBlocks(Input=src))
    sl.SliceType = "Plane"
    sl.SliceType.Origin = origin
    sl.SliceType.Normal = [0.0, 1.0, 0.0]
    return sl


def scalar(proxy, name):
    data = sm.Fetch(proxy)
    arr = data.GetCellData().GetArray(name) or data.GetPointData().GetArray(name)
    return float(arr.GetValue(0)) if arr is not None else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--purge-fields", action="store_true")
    a = ap.parse_args()
    run = os.path.abspath(a.run_dir)
    with open(os.path.join(run, "case_meta.json")) as f:
        meta = json.load(f)
    pcm = meta.get("kind") == "pcm_embedded"
    cfoam = os.path.join(run, "case", "case.foam")
    out = os.path.join(run, "results", "pv")
    os.makedirs(out, exist_ok=True)
    s = meta["slab_m"]
    Lx, H = s["Lx"], s["H"]
    area = meta["face_area_m2"]
    t_end = meta["time"]["end_time_s"]
    clock0 = meta["time"].get("start_clock_h", 0.0)
    cy = meta["body"]["centre"][1] if pcm else 0.0
    origin = [0.0, cy + 1.1e-4, 0.5 * H]

    probe = OpenFOAMReader(FileName=cfoam)
    probe.UpdatePipelineInformation()
    times = [t for t in probe.TimestepValues if t >= t_end - DAY - 1.0]
    if not times:
        raise SystemExit(f"no field snapshots on the analysed day in {cfoam}")

    regions = ["concrete", "pcm"] if pcm else ["concrete"]
    rd_c = reader(cfoam, blocks(probe, "concrete", "internalMesh"), ["T"])
    sl_c = Calculator(Input=plane_slice(rd_c, origin), ResultArrayName="T_C", Function="T-273.15",
                      AttributeType="Cell Data")
    if pcm:
        rd_p = reader(cfoam, blocks(probe, "pcm", "internalMesh"), ["T", "liquidFraction"])
        sl_p0 = plane_slice(rd_p, origin)
        sl_p = Calculator(Input=sl_p0, ResultArrayName="T_C", Function="T-273.15",
                          AttributeType="Cell Data")
        front = Contour(Input=CellDatatoPointData(Input=plane_slice(
            reader(cfoam, blocks(probe, "pcm", "internalMesh"), ["liquidFraction"]), origin)))
        front.ContourBy = ["POINTS", "liquidFraction"]
        front.Isosurfaces = [0.5]
    # surface metrics: patch area averages of T and wallHeatFlux
    faces = {}
    for patch in ("roof", "ceiling"):
        rp = reader(cfoam, blocks(probe, "concrete", patch), ["T", "wallHeatFlux"])
        faces[patch] = IntegrateVariables(Input=MergeBlocks(Input=rp), DivideCellDataByVolume=1)

    rows = []
    for t in times:
        hh = (clock0 + (t - (t_end - DAY)) / 3600.0)
        tag = f"{hh:04.1f}".replace(".", "p")
        row = {"time_s": t, "hour": hh}
        for patch, integ in faces.items():
            UpdatePipeline(time=t, proxy=integ)
            row[f"T_{patch}_C"] = scalar(integ, "T") - 273.15
            q = scalar(integ, "wallHeatFlux")       # into the concrete
            row[f"q_{patch}_into_slab_W_m2"] = q
        rows.append(row)

        UpdatePipeline(time=t, proxy=sl_c)
        SaveData(os.path.join(out, f"slice_concrete_{tag}h.csv"), proxy=sl_c,
                 FieldAssociation="Cell Data")
        if pcm:
            UpdatePipeline(time=t, proxy=sl_p)
            SaveData(os.path.join(out, f"slice_pcm_{tag}h.csv"), proxy=sl_p,
                     FieldAssociation="Cell Data")

        rv = CreateView("RenderView")
        rv.ViewSize = [1400, 760]
        rv.ViewTime = t
        rv.OrientationAxesVisibility = 0
        rv.Background = [1, 1, 1]
        if hasattr(rv, "UseColorPaletteForBackground"):
            rv.UseColorPaletteForBackground = 0
        lut = GetColorTransferFunction("T_C")
        lut.AutomaticRescaleRangeMode = "Never"
        for src in ([sl_c, sl_p] if pcm else [sl_c]):
            d = Show(src, rv)
            ColorBy(d, ("CELLS", "T_C"))
            d.SetScalarBarVisibility(rv, True)
        # one fixed scale for every case and time (after ColorBy, which rescales)
        lut.ApplyPreset("Inferno (matplotlib)", True)
        lut.RescaleTransferFunction(T_RANGE[0], T_RANGE[1])
        bar = GetScalarBar(lut, rv)
        bar.AutomaticLabelFormat = 0
        bar.LabelFormat = "%-.0f"
        bar.RangeLabelFormat = "%-.0f"
        bar.Title = "T (C)"
        bar.ComponentTitle = ""
        bar.TitleColor = [0, 0, 0]
        bar.LabelColor = [0, 0, 0]
        if pcm:
            fd = Show(front, rv)
            fd.DiffuseColor = [1, 1, 1]
            fd.AmbientColor = [1, 1, 1]
            fd.LineWidth = 3.0
            fd.ColorArrayName = ["POINTS", ""]      # solid colour (not ColorBy(rep, None))
            outline = Show(sl_p0, rv)
            outline.SetRepresentationType("Feature Edges")
            outline.ColorArrayName = ["CELLS", ""]
            outline.DiffuseColor = [0.2, 0.9, 1.0]
            outline.AmbientColor = [0.2, 0.9, 1.0]
            outline.LineWidth = 2.0
        name = meta.get("case_id", "pure concrete baseline")
        txt = Show(Text(Text=f"{name}   {int(hh) % 24:02d}:{int(round(60 * (hh % 1))) % 60:02d} "
                             f"on the analysed day   (white: melting front, cyan: PCM outline)"), rv)
        txt.FontSize = 14
        txt.Color = [0, 0, 0]
        txt.WindowLocation = "Upper Center"
        rv.CameraPosition = [0.0, origin[1] - 1.0, 0.5 * H]
        rv.CameraFocalPoint = [0.0, origin[1], 0.5 * H]
        rv.CameraViewUp = [0, 0, 1]
        rv.CameraParallelProjection = 1
        rv.CameraParallelScale = 0.62 * max(H, Lx * 760 / 1400.0)
        Render(rv)
        SaveScreenshot(os.path.join(out, f"contour_T_{tag}h.png"), rv, ImageResolution=[1400, 760])
        print(f"  {name}: {hh:5.1f} h  roof {row['T_roof_C']:.2f} C  ceiling {row['T_ceiling_C']:.2f} C")

    with open(os.path.join(out, "surface_metrics.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: f"{v:.6g}" for k, v in r.items()})

    if a.purge_fields:
        case = os.path.join(run, "case")
        tdirs = sorted((float(d), d) for d in os.listdir(case) if _is_time(d))
        keep = {"0", tdirs[-1][1]} if tdirs else {"0"}
        for _, d in tdirs:
            if d not in keep:
                shutil.rmtree(os.path.join(case, d))
        print(f"  purged {len(tdirs) - len(keep & {d for _, d in tdirs})} field folders")
    print(f"  -> {out}")


def _is_time(d):
    try:
        float(d)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    main()
