"""Temperature and velocity contours on the panel symmetry plane (paper Fig. 9).

Run with ParaView's interpreter:
    pvpython post/pv_contours.py runs/<case> --times 50 150 250 300
Writes runs/<case>/post/contour_{T,U}_<min>min.png and the slice data as CSV.
Needs only the standard library + paraview.simple (no pcmpipe import).
"""
import argparse
import json
import os

from paraview.simple import (Calculator, ColorBy, CreateView, GetColorTransferFunction,
                             GetScalarBar, MergeBlocks, OpenFOAMReader, Render, SaveData, SaveScreenshot,
                             Show, Slice, Text, UpdatePipeline)


def internal_regions(reader, wanted):
    avail = list(reader.MeshRegions.Available)
    out = [r for r in avail if r.endswith("internalMesh") and any(w in r.split("/") for w in wanted)]
    if not out:
        raise RuntimeError(f"regions {wanted} not found in {avail}")
    return out


def nearest_time(reader, t_s):
    times = list(reader.TimestepValues)
    return min(times, key=lambda x: abs(x - t_s)) if times else 0.0


def make_slice(cfoam, regions, arrays, plane):
    rd = OpenFOAMReader(FileName=cfoam)
    rd.MeshRegions = internal_regions(rd, regions)
    rd.CellArrays = arrays
    rd.UpdatePipelineInformation()
    merged = MergeBlocks(Input=rd)
    sl = Slice(Input=merged)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = plane["origin"]
    sl.SliceType.Normal = plane["normal"]
    return rd, sl


def render(source, array, rng, label, title, path, view_box, size, t):
    rv = CreateView("RenderView")
    rv.ViewSize = size
    rv.ViewTime = t
    rv.OrientationAxesVisibility = 0
    rv.Background = [1, 1, 1]
    if hasattr(rv, "UseColorPaletteForBackground"):   # ParaView >= 5.10 ignores
        rv.UseColorPaletteForBackground = 0            # Background otherwise
    d = Show(source, rv)
    ColorBy(d, ("POINTS", array))
    lut = GetColorTransferFunction(array)
    lut.ApplyPreset("Jet", True)
    lut.RescaleTransferFunction(rng[0], rng[1])
    lut.AutomaticRescaleRangeMode = "Never"
    d.SetScalarBarVisibility(rv, True)
    bar = GetScalarBar(lut, rv)
    bar.Title = label
    bar.ComponentTitle = ""
    bar.TitleColor = [0, 0, 0]          # the default white is invisible on white
    bar.LabelColor = [0, 0, 0]
    bar.TitleFontSize = 18
    bar.LabelFontSize = 16
    txt = Text(Text=title)
    td = Show(txt, rv)
    td.FontSize = 14
    td.Color = [0, 0, 0]
    td.WindowLocation = "Upper Center"
    (x0, x1), (z0, z1), y = view_box
    cx, cz = 0.5 * (x0 + x1), 0.5 * (z0 + z1)
    rv.CameraPosition = [cx, y - 1.0, cz]
    rv.CameraFocalPoint = [cx, y, cz]
    rv.CameraViewUp = [0, 0, 1]
    rv.CameraParallelProjection = 1
    aspect = size[0] / float(size[1])
    rv.CameraParallelScale = max(0.5 * (z1 - z0), 0.5 * (x1 - x0) / aspect)
    Render(rv)
    SaveScreenshot(path, rv, ImageResolution=size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--times", nargs="+", type=float, default=[50, 150, 250, 300])
    a = ap.parse_args()

    with open(os.path.join(a.case_dir, "case_meta.json")) as f:
        meta = json.load(f)
    with open(os.path.join(a.case_dir, "case.json")) as f:
        cj = json.load(f)
    post = cj["post"]
    # model.air: coefficient -- the air is meshed but not solved (no fields after 0/)
    air = ["air"] if cj.get("model", {}).get("air", "resolved") == "resolved" else []
    cfoam = os.path.join(a.case_dir, "case", "case.foam")
    out = os.path.join(a.case_dir, "post")
    os.makedirs(out, exist_ok=True)

    plane = meta["slice_plane"]
    b = meta["panel_bounds_m"]
    cx = 0.5 * (b["x"][0] + b["x"][1])
    W, D = post["view_half_width_m"], post["view_depth_m"]
    if not air:      # nothing to show around the panel: frame the panel itself
        W = 0.6 * (b["x"][1] - b["x"][0])
        D = 1.25 * (b["z"][1] - b["z"][0])
    # a white band above the ceiling (no mesh there) keeps the title off the panel
    view_box = ((cx - W, cx + W), (-D, 0.08 * D), plane["origin"][1])
    size = post["image_size"]

    rT, sT = make_slice(cfoam, [*air, "shell", "pcm"], ["T"], plane)
    cT = Calculator(Input=sT, ResultArrayName="T_C", Function="T-273.15")
    rU, sU = make_slice(cfoam, [*air, "pcm"], ["U"], plane)
    cU = Calculator(Input=sU, ResultArrayName="Umag_mm_s", Function="mag(U)*1000")

    last = max(list(rT.TimestepValues), default=0.0)
    for t_min in a.times:
        if 60.0 * t_min > last + 1.0:    # run ended earlier (melted, or a wall-clock cap)
            print(f"  skip {t_min:g} min: last written time is {last / 60:.0f} min")
            continue
        for rd, calc, arr, rng, label in (
                (rT, cT, "T_C", post["T_range_C"], "T (C)"),
                (rU, cU, "Umag_mm_s", post["U_range_mm_s"], "|U| (mm/s)")):
            t = nearest_time(rd, 60.0 * t_min)
            UpdatePipeline(time=t, proxy=calc)
            tag = "T" if arr == "T_C" else "U"
            when = f"{t:.1f} s" if t < 60 else f"{t / 60:.0f} min"
            title = f"{meta['id']}  {label}  t = {when}"
            render(calc, arr, rng, label, title,
                   os.path.join(out, f"contour_{tag}_{t_min:g}min.png"), view_box, size, t)
            SaveData(os.path.join(out, f"slice_{tag}_{t_min:g}min.csv"), proxy=calc)
            print(f"  wrote contour_{tag}_{t_min:g}min.png (t = {t:.0f} s)")


if __name__ == "__main__":
    main()
