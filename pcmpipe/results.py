"""Step 4b -- read function-object output and compute the paper's metrics.

Metrics (Sec. 4):
  t80, t100        first time the PCM volume-averaged liquid fraction reaches 0.8
                   and `melt_threshold_100` (default 0.999) [min]
  melting rate     m_PCM * d(gamma)/dt  [g/min] (Fig. 10, right axis)
  average rate     m_PCM / t100 [g/min] (Sec. 4.1: 0.582 / 0.482 g/min)
  normalised rate  average rate / panel volume [g/min/L] (Fig. 18)
Paper convention: Fig. 18 values equal 170 g / t100 / V for every volume, i.e. the
250 mL PCM mass is used for all volumes. `normalised_rate_paper` reproduces it;
`normalised_rate` uses the PCM mass actually in the panel.
"""
from __future__ import annotations

import csv
import glob
import json
import os

import numpy as np

PAPER_REFERENCE_MASS_G = 170.0


def _read_table(paths):
    """Concatenate whitespace tables (with '#' headers) from restart folders."""
    rows, header = [], None
    for p in sorted(paths, key=lambda s: float(os.path.basename(os.path.dirname(s)) or 0)):
        with open(p) as f:
            for line in f:
                if line.startswith("#"):
                    if "Time" in line:
                        header = line[1:].split()
                    continue
                parts = line.replace("(", " ").replace(")", " ").split()
                if parts:
                    rows.append([float(x) for x in parts])
    if not rows:
        return header, np.empty((0, 0))
    data = np.array(rows)
    _, idx = np.unique(data[:, 0], return_index=True)     # drop restart overlaps
    return header, data[np.sort(idx)]


def _find(cdir, fo_name, filename):
    return [p for p in glob.glob(os.path.join(cdir, "postProcessing", "**", filename), recursive=True)
            if f"{os.sep}{fo_name}{os.sep}" in p]


def liquid_fraction_series(cdir):
    header, data = _read_table(_find(cdir, "pcmAverages", "volFieldValue.dat"))
    if data.size == 0:
        raise FileNotFoundError(f"no pcmAverages output in {cdir}/postProcessing")
    cols = {name: i for i, name in enumerate(header)}
    g = data[:, cols.get("volAverage(liquidFraction)", 1)]
    T = data[:, cols.get("volAverage(T)", 2)]
    return data[:, 0], g, T


def probe_series(cdir):
    paths = _find(cdir, "pcmProbes", "T")
    _, data = _read_table(paths)
    return data   # columns: time, top, mid, bot


def first_crossing(t, y, level):
    i = np.argmax(y >= level)
    if y[i] < level:
        return None
    if i == 0:
        return t[0]
    return t[i - 1] + (level - y[i - 1]) * (t[i] - t[i - 1]) / (y[i] - y[i - 1])


def metrics(case_dir: str, melt_threshold_100: float = 0.999) -> dict:
    with open(os.path.join(case_dir, "case_meta.json")) as f:
        meta = json.load(f)
    cdir = os.path.join(case_dir, "case")
    t, g, T = liquid_fraction_series(cdir)
    t_min = t / 60.0
    m = meta.get("pcm_mass_g") or PAPER_REFERENCE_MASS_G
    t80 = first_crossing(t_min, g, 0.8)
    t100 = first_crossing(t_min, g, melt_threshold_100)
    rate = m * np.gradient(g, t_min) if len(t_min) > 2 else np.zeros_like(t_min)
    V_L = meta["volume_ml"] / 1000.0
    out = dict(meta)
    out.update(t80_min=t80, t100_min=t100, final_liquid_fraction=float(g[-1]),
               end_time_min=float(t_min[-1]))
    if t100:
        out["average_melting_rate_g_min"] = m / t100
        out["normalised_rate_g_min_L"] = m / t100 / V_L
        out["normalised_rate_paper_g_min_L"] = PAPER_REFERENCE_MASS_G / t100 / V_L
        out["average_melting_rate_per_area_g_min_cm2"] = (
            m / t100 / meta.get("exposed_area_cm2", meta["exposed_area_cm2_analytic"]))
    with open(os.path.join(case_dir, "series.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_min", "liquid_fraction", "T_mean_K", "melting_rate_g_min"])
        w.writerows(zip(t_min, g, T, rate))
    with open(os.path.join(case_dir, "metrics.json"), "w") as f:
        json.dump(out, f, indent=2)
    return out
