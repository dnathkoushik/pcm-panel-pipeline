"""Volume-averaged PCM liquid fraction and temperature for every written time.

    pvpython post/pv_extract_transient.py runs/<case>
Writes runs/<case>/post/pv_transient.csv. This cross-checks the in-run
`pcmAverages` function object (which samples every 30 s) on the written fields.
"""
import argparse
import csv
import os

from paraview import servermanager
from paraview.simple import IntegrateVariables, OpenFOAMReader, UpdatePipeline


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    a = ap.parse_args()

    rd = OpenFOAMReader(FileName=os.path.join(a.case_dir, "case", "case.foam"))
    pcm = [r for r in rd.MeshRegions.Available if r.endswith("internalMesh") and "pcm" in r.split("/")]
    rd.MeshRegions = pcm
    rd.CellArrays = ["liquidFraction", "T"]
    rd.UpdatePipelineInformation()
    integ = IntegrateVariables(Input=rd)
    integ.DivideCellDataByVolume = 1

    out = os.path.join(a.case_dir, "post")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "pv_transient.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_min", "liquid_fraction", "T_mean_K", "pcm_volume_m3"])
        for t in rd.TimestepValues:
            UpdatePipeline(time=t, proxy=integ)
            data = servermanager.Fetch(integ)
            cd = data.GetCellData()
            lf = cd.GetArray("liquidFraction")
            T = cd.GetArray("T")
            vol = cd.GetArray("Volume")
            w.writerow([t / 60.0,
                        lf.GetValue(0) if lf else "",
                        T.GetValue(0) if T else "",
                        vol.GetValue(0) if vol else ""])
    print(f"  wrote {path}")


if __name__ == "__main__":
    main()
