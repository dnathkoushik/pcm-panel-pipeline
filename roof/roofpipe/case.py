"""Phase 1: the pure concrete roof slab as an OpenFOAM case (templates/concrete_slab).

    <runs>/<name>/
        case/                    OpenFOAM case: blockMesh, then chtMultiRegionFoam on
                                 the region `concrete` (scripts/run_slab_case.sh)
        boundary_conditions.csv  outdoor air, GHI, sol-air and indoor temperatures on
                                 the run's time axis (the tables in 0/concrete/T)
        reference_1d.csv         the same problem solved in 1D (slab1d.py)
        case_meta.json           every derived number: mesh, coefficients, initial state
        config_used.yaml         the configuration after --set overrides
"""
from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil

from . import __version__, config, slab1d, weather

TEMPLATE = "concrete_slab"
DAY_S = weather.DAY_S
_TOKEN = re.compile(r"@([A-Z][A-Z0-9_]*)@")


def num(x) -> str:
    return f"{float(x):.9g}"


def render(src_dir: str, dst_dir: str, values: dict):
    """Copy the template tree, filling every @TOKEN@; fail on a token without a value."""
    missing = []
    for dirpath, _, files in os.walk(src_dir):
        rel = os.path.relpath(dirpath, src_dir)
        for name in sorted(files):
            with open(os.path.join(dirpath, name), encoding="utf-8") as f:
                text = f.read()

            def fill(m, where=os.path.join(rel, name)):
                if m.group(1) not in values:
                    missing.append(f"{where}: {m.group(0)}")
                    return m.group(0)
                return str(values[m.group(1)])

            text = _TOKEN.sub(fill, text)
            dst = os.path.normpath(os.path.join(dst_dir, rel, name))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with open(dst, "w", encoding="utf-8", newline="\n") as f:
                f.write(text)
    if missing:
        raise SystemExit("template placeholders without a value:\n  " + "\n  ".join(missing))


def _split2(n: int):
    """n = px * py with px >= py as close as possible (geometric decompositions)."""
    py = math.isqrt(n)
    while n % py:
        py -= 1
    return n // py, py


def _whole_steps(name: str, value: float, dt: float):
    steps = value / dt
    if abs(steps - round(steps)) > 1e-9 * max(1.0, steps):
        raise SystemExit(f"{name} = {value:g} s is not a whole number of time steps "
                         f"(run.deltaT = {dt:g} s)")


def derive(cfg: dict, n_procs: int):
    """Template values and metadata of the baseline case."""
    s, mat, m, run = cfg["slab"], cfg["materials"]["concrete"], cfg["mesh"], cfg["run"]
    Lx, Ly, H = (float(s[k]) for k in ("Lx", "Ly", "H"))
    k, rho, cp = (float(mat[x]) for x in ("k", "rho", "cp"))
    nx = max(1, round(Lx / float(m["cell_xy"])))
    ny = max(1, round(Ly / float(m["cell_xy"])))
    nz = max(2, round(H / float(m["cell_z"])))
    dx, dy, dz = Lx / nx, Ly / ny, H / nz

    dt = float(run["deltaT"])
    n_days = int(run["n_days"])
    if n_days < 2:
        raise SystemExit("run.n_days must be >= 2 (spin-up, then the analysed day)")
    end = n_days * DAY_S
    sample = float(run["sample_interval_s"])
    profile = float(run["profile_interval_s"])
    write = 3600.0 * float(run["write_interval_h"])
    for name, value in (("run.sample_interval_s", sample), ("run.profile_interval_s", profile),
                        ("run.write_interval_h", write), ("run.n_days", end),
                        ("one day", DAY_S)):
        _whole_steps(name, value, dt)

    dd, tab = weather.tables(cfg)

    # Initial state: mid-depth of the steady profile under the daily-mean conditions
    first_day = [i for i, t in enumerate(tab["time_s"]) if t < DAY_S]
    Tsa_mean = sum(tab["T_solair_C"][i] for i in first_day) / len(first_day)
    Tin_mean = sum(tab["T_indoor_C"][i] for i in first_day) / len(first_day)
    R_out, R_slab, R_in = 1.0 / dd.h_out, H / k, 1.0 / dd.h_in_eff
    q_mean = (Tsa_mean - Tin_mean) / (R_out + R_slab + R_in)
    T0_C = Tsa_mean - q_mean * (R_out + 0.5 * R_slab)

    # T(z) probes at the cell centres of the column at (or next to) the slab centre
    xc = 0.5 * dx if nx % 2 == 0 else 0.0
    yc = 0.5 * dy if ny % 2 == 0 else 0.0
    z_probe = [(i + 0.5) * dz for i in range(nz)]
    probes = "\n".join(f"            ({num(xc)} {num(yc)} {num(z)})" for z in z_probe)

    method = str(run["decomposition"])
    if n_procs > 1 and method not in ("scotch", "ptscotch", "hierarchical", "simple", "kahip",
                                      "metis"):
        raise SystemExit(f"run.decomposition: unsupported method {method!r}")
    px, py = _split2(n_procs)

    values = {
        "X0": num(-0.5 * Lx), "X1": num(0.5 * Lx), "Y0": num(-0.5 * Ly), "Y1": num(0.5 * Ly),
        "H": num(H), "NX": nx, "NY": ny, "NZ": nz,
        "K_CONCRETE": num(k), "RHO_CONCRETE": num(rho), "CP_CONCRETE": num(cp),
        "T0": num(T0_C + 273.15),
        "TA_ROOF": weather.foam_rows(tab["time_s"], tab["T_solair_C"], " " * 12),
        "TA_CEILING": weather.foam_rows(tab["time_s"], tab["T_indoor_C"], " " * 12),
        "H_OUT": num(dd.h_out), "H_IN": num(dd.h_in_eff), "INDOOR_MODE": dd.indoor_mode,
        "END_TIME": num(end), "DELTA_T": num(dt), "N_DAYS": n_days,
        "WRITE_INTERVAL": num(write), "SAMPLE_INTERVAL": num(sample),
        "PROFILE_INTERVAL": num(profile), "PROFILE_PROBES": probes,
        "H_TOL": num(run["energy_tolerance"]),
        "NP": n_procs, "DECOMP_METHOD": method, "PX": px, "PY": py,
    }
    meta = {
        "roofpipe_version": __version__,
        "phase": 1, "kind": "pure_concrete_baseline",
        "slab_m": {"Lx": Lx, "Ly": Ly, "H": H},
        "face_area_m2": Lx * Ly,
        "concrete": {"k": k, "rho": rho, "cp": cp,
                     "heat_capacity_J_m2K": rho * cp * H,
                     "diffusivity_m2_s": k / (rho * cp)},
        "mesh": {"nx": nx, "ny": ny, "nz": nz, "cells": nx * ny * nz,
                 "dx": dx, "dy": dy, "dz": dz},
        "time": {"deltaT_s": dt, "n_days": n_days, "end_time_s": end,
                 "analysed_day_s": [end - DAY_S, end],
                 "start_clock_h": float(run.get("start_clock_h", 0.0)),
                 "sample_interval_s": sample, "write_interval_s": write},
        "boundary": {
            "h_out": dd.h_out, "h_in": dd.h_in, "h_in_effective": dd.h_in_eff,
            "indoor_mode": dd.indoor_mode, "absorptance": dd.absorptance,
            "longwave_correction_K": dd.longwave_K,
            "T_solair_C": [min(tab["T_solair_C"]), max(tab["T_solair_C"])],
            "T_out_C": [min(tab["T_out_C"]), max(tab["T_out_C"])],
            "T_indoor_C": [min(tab["T_indoor_C"]), max(tab["T_indoor_C"])],
            "GHI_peak_W_m2": max(tab["GHI_W_m2"]),
        },
        "sun_clock_h": {"sunrise": dd.sun.sunrise_h, "solar_noon": dd.sun.noon_h,
                        "sunset": dd.sun.sunset_h},
        "initial": {"T0_C": T0_C, "T_solair_mean_C": Tsa_mean, "T_indoor_mean_C": Tin_mean,
                    "q_mean_steady_W_m2": q_mean,
                    "U_value_W_m2K": 1.0 / (R_out + R_slab + R_in)},
        "probes": {"x": xc, "y": yc, "z": z_probe},
        "n_procs": n_procs, "decomposition": method,
    }
    return values, meta, dd, tab


def _write_csv(path: str, cols: dict):
    keys = list(cols)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(keys)
        for row in zip(*(cols[k] for k in keys)):
            w.writerow([f"{v:.8g}" for v in row])


def build_baseline(cfg: dict, name: str, n_procs: int, force: bool = False) -> str:
    rdir = os.path.abspath(name) if os.sep in name else os.path.join(config.runs_dir(), name)
    cdir = os.path.join(rdir, "case")
    if " " in cdir:
        raise SystemExit(f"the case path contains spaces ({cdir}): OpenFOAM cannot handle "
                         "that. Clone the repository to a path without spaces, or set "
                         "ROOF_RUNS_DIR.")
    if os.path.exists(cdir):
        if not force:
            raise SystemExit(f"{cdir} already exists. Use --force to rewrite it (this deletes "
                             "its results); the job script continues an existing case instead.")
        shutil.rmtree(cdir)
    os.makedirs(rdir, exist_ok=True)

    values, meta, dd, tab = derive(cfg, n_procs)
    render(os.path.join(config.TEMPLATES, TEMPLATE), cdir, values)
    open(os.path.join(cdir, "case.foam"), "w").close()          # ParaView reader handle
    _write_csv(os.path.join(rdir, "boundary_conditions.csv"), tab)

    s, mat, run = cfg["slab"], meta["concrete"], meta["time"]
    ref = slab1d.solve(
        H=float(s["H"]), nz=meta["mesh"]["nz"], k=mat["k"], rho=mat["rho"], cp=mat["cp"],
        h_top=dd.h_out, Ta_top=slab1d.table_function(tab["time_s"], tab["T_solair_C"]),
        h_bot=dd.h_in_eff, Ta_bot=slab1d.table_function(tab["time_s"], tab["T_indoor_C"]),
        T0=meta["initial"]["T0_C"], dt=run["deltaT_s"], t_end=run["end_time_s"],
        sample_dt=run["sample_interval_s"])
    _write_csv(os.path.join(rdir, "reference_1d.csv"), ref)

    with open(os.path.join(rdir, "case_meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    config.dump(cfg, os.path.join(rdir, "config_used.yaml"))

    last = [i for i, t in enumerate(ref["time_s"]) if t >= run["end_time_s"] - DAY_S]
    gain = [ref["q_ceiling_gain_W_m2"][i] for i in last]
    b = meta["boundary"]
    print(f"  case written : {cdir}")
    print(f"  mesh         : {meta['mesh']['nx']} x {meta['mesh']['ny']} x "
          f"{meta['mesh']['nz']} = {meta['mesh']['cells']} cells, region concrete")
    print(f"  run          : {run['n_days']} days, dt {run['deltaT_s']:g} s, "
          f"{n_procs} MPI rank(s)")
    print(f"  sun (IST)    : rise {dd.sun.sunrise_h:.2f} h, noon {dd.sun.noon_h:.2f} h, "
          f"set {dd.sun.sunset_h:.2f} h; GHI peak {b['GHI_peak_W_m2']:.0f} W/m2")
    print(f"  sol-air      : {b['T_solair_C'][0]:.1f} .. {b['T_solair_C'][1]:.1f} C; "
          f"indoor ({dd.indoor_mode}) {b['T_indoor_C'][0]:.1f} .. {b['T_indoor_C'][1]:.1f} C")
    print(f"  1D reference : ceiling gain {min(gain):.1f} .. {max(gain):.1f} W/m2 on the "
          f"last day (mean {sum(gain) / len(gain):.1f}) -> {rdir}/reference_1d.csv")
    return rdir
