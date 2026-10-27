"""3D melting front: PCM liquid fraction in a cut-away view of the panel.

Run with ParaView's interpreter:
    pvpython post/pv_liquid_fraction_3d.py runs/<case> --times 60 120 240
Writes runs/<case>/post/liquid_fraction_3d_<min>min.png. The PCM is cut on the
panel symmetry plane (front half removed) and coloured by liquid fraction
(0 solid, 1 liquid); the gamma = 0.5 iso-surface, i.e. the melting front, is
drawn in dark grey; the uncut PCM is a faint silhouette. The camera looks up at
the panel from below the ceiling. Needs only paraview.simple (no pcmpipe import).
"""
import argparse
import json
import math
import os

from paraview.simple import (CellDatatoPointData, Clip, ColorBy, Contour, CreateView,
                             ExtractSurface, GetAnimationScene, GetColorTransferFunction,
                             GetScalarBar, MergeBlocks, OpenFOAMReader, Render, ResetCamera,
                             SaveScreenshot, Show, Text)


def nearest_time(reader, t_s):
    times = list(reader.TimestepValues)
    return min(times, key=lambda x: abs(x - t_s)) if times else 0.0


def plane_clip(source, plane):
    c = Clip(Input=source)
    c.ClipType = "Plane"
    c.ClipType.Origin = plane["origin"]
    c.ClipType.Normal = plane["normal"]
    return c


def camera(rv, meta):
    """Below the ceiling, on the cut side of the symmetry plane, slightly to the side."""
    b = meta["panel_bounds_m"]
    c = [0.5 * (b["x"][0] + b["x"][1]), 0.5 * (b["y"][0] + b["y"][1]), 0.5 * b["z"][0]]
    n = meta["slice_plane"]["normal"]
    side = [n[1], -n[0], 0.0]                      # horizontal, in the cut plane
    d = [n[i] + 0.35 * side[i] for i in range(3)]
    d[2] = -0.55
    norm = math.sqrt(sum(v * v for v in d))
    rv.CameraFocalPoint = c
    rv.CameraPosition = [c[i] + d[i] / norm for i in range(3)]
    rv.CameraViewUp = [0.0, 0.0, 1.0]
    ResetCamera(rv)
    rv.GetActiveCamera().Dolly(1.15)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--times", nargs="+", type=float, default=[50, 150, 250, 300])
    a = ap.parse_args()

    with open(os.path.join(a.case_dir, "case_meta.json")) as f:
        meta = json.load(f)
    with open(os.path.join(a.case_dir, "case.json")) as f:
        post = json.load(f)["post"]
    out = os.path.join(a.case_dir, "post")
    os.makedirs(out, exist_ok=True)

    rd = OpenFOAMReader(FileName=os.path.join(a.case_dir, "case", "case.foam"))
    pcm = [r for r in rd.MeshRegions.Available
           if r.endswith("internalMesh") and "pcm" in r.split("/")]
    if not pcm:
        raise RuntimeError(f"pcm region not found in {list(rd.MeshRegions.Available)}")
    rd.MeshRegions = pcm
    rd.CellArrays = ["liquidFraction"]
    rd.UpdatePipelineInformation()
    merged = MergeBlocks(Input=rd)
    points = CellDatatoPointData(Input=merged)
    plane = meta["slice_plane"]
    cut = plane_clip(points, plane)
    front = Contour(Input=points)
    front.ContourBy = ["POINTS", "liquidFraction"]
    front.Isosurfaces = [0.5]
    front_cut = plane_clip(front, plane)
    silhouette = ExtractSurface(Input=merged)

    size = post["image_size"]
    rv = CreateView("RenderView")
    rv.ViewSize = size
    rv.OrientationAxesVisibility = 0
    rv.Background = [1, 1, 1]
    if hasattr(rv, "UseColorPaletteForBackground"):   # ParaView >= 5.10 ignores
        rv.UseColorPaletteForBackground = 0            # Background otherwise

    d = Show(cut, rv)
    ColorBy(d, ("POINTS", "liquidFraction"))
    lut = GetColorTransferFunction("liquidFraction")
    lut.ApplyPreset("Cool to Warm", True)
    lut.RescaleTransferFunction(0.0, 1.0)
    lut.AutomaticRescaleRangeMode = "Never"
    d.SetScalarBarVisibility(rv, True)
    bar = GetScalarBar(lut, rv)
    bar.Title = "Liquid fraction"
    bar.ComponentTitle = ""
    bar.TitleColor = [0, 0, 0]
    bar.LabelColor = [0, 0, 0]

    # solid colours; ColorBy(rep, None) fails on ParaView 5.11 for an uncoloured rep
    df = Show(front_cut, rv)
    df.ColorArrayName = ["POINTS", ""]
    df.DiffuseColor = [0.25, 0.25, 0.25]
    ds = Show(silhouette, rv)
    ds.ColorArrayName = ["POINTS", ""]
    ds.DiffuseColor = [0.6, 0.6, 0.6]
    ds.Opacity = 0.12

    title = Text()
    dt = Show(title, rv)
    dt.FontSize = 14
    dt.Color = [0, 0, 0]
    dt.WindowLocation = "Upper Center"

    camera(rv, meta)
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    clock = meta.get("start_clock_h")
    last = max(list(rd.TimestepValues), default=0.0)
    for t_min in a.times:
        if 60.0 * t_min > last + 1.0:    # run ended earlier (melted, or a wall-clock cap)
            print(f"  skip {t_min:g} min: last written time is {last / 60:.0f} min")
            continue
        t = nearest_time(rd, 60.0 * t_min)
        scene.AnimationTime = t
        rv.ViewTime = t
        when = f"t = {t:.1f} s" if t < 60 else f"t = {t / 60:.0f} min"
        if clock is not None:
            h = clock + t / 3600.0
            when += f" ({int(h) % 24:02d}:{int(round(60 * (h % 1))) % 60:02d})"
        title.Text = f"{meta['id']}   liquid fraction, melting front (grey)   {when}"
        Render(rv)
        path = os.path.join(out, f"liquid_fraction_3d_{t_min:g}min.png")
        SaveScreenshot(path, rv, ImageResolution=size)
        print(f"  wrote {os.path.basename(path)} (t = {t:.0f} s)")


if __name__ == "__main__":
    main()
