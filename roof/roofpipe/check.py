"""Phase 1 gates: is the pure concrete baseline right?

Reads only the function-object time series (case/postProcessing/**), never the
fields, and writes <run>/results/:
    timeseries_lastday.csv   the analysed (last) day, every sample_interval_s
    summary.json             gates, daily metrics, deviations
    quicklook.png            temperatures and heat fluxes (if matplotlib is available)

Gates (all must pass before Phase 2):
    complete    the solver reached endTime
    periodic    surface temperatures of the last day repeat the day before within
                PERIODIC_K: the spin-up has washed out the initial state
    reference   surface temperatures equal the 1D solution (same discretisation,
                reference_1d.csv) within REFERENCE_K: boundary conditions, tables,
                properties and the 3D mesh are what the configuration says
    energy      last day: heat in through the roof - heat out through the ceiling -
                change of stored heat, within ENERGY_FRAC of the gross roof heat flow

Sign convention: q_roof_in > 0 heats the slab from outside; q_ceiling_gain > 0 is
heat gained by the room (cooling load), < 0 heat lost by it (heating load).
"""
from __future__ import annotations

import csv
import glob
import json
import os

from . import slab1d

DAY_S = 86400.0
PERIODIC_K = 0.05
REFERENCE_K = 0.05
ENERGY_FRAC = 0.01


# ----------------------------------------------------------------- reading
def _read_csv(path: str) -> dict:
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    head, data = rows[0], rows[1:]
    return {h: [float(r[i]) for r in data] for i, h in enumerate(head)}


def _files(case_dir: str, fo_name: str, pattern: str) -> list:
    """Output files of one function object; one per start time (restarts), oldest first."""
    found = [p for p in glob.glob(os.path.join(case_dir, "postProcessing", "**", pattern),
                                  recursive=True)
             if fo_name in os.path.normpath(p).split(os.sep)]

    def start(p):
        try:
            return float(os.path.basename(os.path.dirname(p)))
        except ValueError:
            return 0.0

    return sorted(found, key=lambda p: (start(p), p))


def _rows(path: str):
    """(header names, data token lists) of an OpenFOAM .dat file."""
    header, data = None, []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                names = s.lstrip("#").split()
                if names and names[0] == "Time":
                    header = names
                continue
            data.append(s.split())
    return header, data


def scalar_series(case_dir: str, fo_name: str, pattern: str) -> dict:
    """{time: value} of a single-value function object (later restarts win)."""
    out = {}
    for p in _files(case_dir, fo_name, pattern):
        for tok in _rows(p)[1]:
            try:
                out[round(float(tok[0]), 3)] = float(tok[1])
            except (ValueError, IndexError):
                continue
    return out


def wall_heat_flux(case_dir: str) -> dict:
    """{(time, patch): integral heat flow [W] into the region}, from wallHeatFlux."""
    out = {}
    for p in _files(case_dir, "surfaceHeatFlux", "wallHeatFlux*.dat"):
        header, data = _rows(p)
        col = 4
        if header:
            for name in ("integral", "Q"):
                if name in header:
                    col = header.index(name)
                    break
        for tok in data:
            try:
                out[(round(float(tok[0]), 3), tok[1])] = float(tok[col])
            except (ValueError, IndexError):
                continue
    return out


# ----------------------------------------------------------------- helpers
def _trapz(t, y) -> float:
    return sum(0.5 * (y[i] + y[i + 1]) * (t[i + 1] - t[i]) for i in range(len(t) - 1))


def _argmax(t, y):
    i = max(range(len(y)), key=y.__getitem__)
    return t[i], y[i]


def _argmin(t, y):
    i = min(range(len(y)), key=y.__getitem__)
    return t[i], y[i]


def _gate(ok, value, limit, unit, note=""):
    return {"pass": bool(ok) if ok is not None else None, "value": value, "limit": limit,
            "unit": unit, "note": note}


# -------------------------------------------------------------------- main
def evaluate(run_dir: str) -> dict:
    case_dir = os.path.join(run_dir, "case")
    with open(os.path.join(run_dir, "case_meta.json"), encoding="utf-8") as f:
        meta = json.load(f)
    bc = _read_csv(os.path.join(run_dir, "boundary_conditions.csv"))
    ref = _read_csv(os.path.join(run_dir, "reference_1d.csv"))
    f_sa = slab1d.table_function(bc["time_s"], bc["T_solair_C"])
    f_in = slab1d.table_function(bc["time_s"], bc["T_indoor_C"])
    f_out = slab1d.table_function(bc["time_s"], bc["T_out_C"])
    f_ghi = slab1d.table_function(bc["time_s"], bc["GHI_W_m2"])

    T_roof = scalar_series(case_dir, "roofSurface", "surfaceFieldValue*.dat")
    T_ceil = scalar_series(case_dir, "ceilingSurface", "surfaceFieldValue*.dat")
    T_mean = scalar_series(case_dir, "slabAverage", "volFieldValue*.dat")
    whf = wall_heat_flux(case_dir)
    if not T_roof or not T_ceil:
        raise SystemExit(f"no roofSurface/ceilingSurface output under {case_dir}/postProcessing "
                         "(has the solver run? see case/log.chtMultiRegionFoam)")

    tm, b, c = meta["time"], meta["boundary"], meta["concrete"]
    h_out, h_in = b["h_out"], b["h_in_effective"]
    area = meta["face_area_m2"]
    H = meta["slab_m"]["H"]
    end, dt = tm["end_time_s"], tm["deltaT_s"]
    times = sorted(set(T_roof) & set(T_ceil))
    t_last = times[-1]
    complete = t_last >= end - 0.5 * dt

    # analysed day: the last full day available
    t0 = t_last - DAY_S
    day = [t for t in times if t >= t0 - 1e-6]
    gates, metrics = {}, {}
    gates["complete"] = _gate(complete, t_last, end, "s", "last sampled time vs endTime")

    # ---- periodicity
    prev = [t for t in day if round(t - DAY_S, 3) in T_roof and round(t - DAY_S, 3) in T_ceil]
    if t0 >= DAY_S - 1e-6 and len(prev) == len(day):
        dev = max(max(abs(T_roof[t] - T_roof[round(t - DAY_S, 3)]),
                      abs(T_ceil[t] - T_ceil[round(t - DAY_S, 3)])) for t in day)
        gates["periodic"] = _gate(dev <= PERIODIC_K, dev, PERIODIC_K, "K",
                                  "max |T_s(t) - T_s(t - 24 h)|, roof and ceiling")
    else:
        gates["periodic"] = _gate(None, None, PERIODIC_K, "K", "less than two days computed")

    # ---- 1D reference (every common sample of the run)
    rT = {round(t, 3): (a, b_) for t, a, b_ in zip(ref["time_s"], ref["T_roof_C"],
                                                   ref["T_ceiling_C"])}
    common = [t for t in times if t in rT]
    if common:
        dev = max(max(abs(T_roof[t] - 273.15 - rT[t][0]), abs(T_ceil[t] - 273.15 - rT[t][1]))
                  for t in common)
        gates["reference"] = _gate(dev <= REFERENCE_K, dev, REFERENCE_K, "K",
                                   f"max |T_s - T_s,1D| over {len(common)} samples")
    else:
        gates["reference"] = _gate(None, None, REFERENCE_K, "K", "no common sample times")

    # ---- last-day series
    clock0 = tm.get("start_clock_h", 0.0)
    rows = []
    rq = {round(t, 3): (q1, q2) for t, q1, q2 in zip(ref["time_s"], ref["q_roof_in_W_m2"],
                                                     ref["q_ceiling_gain_W_m2"])}
    k_eff = 1.0 / (1.0 / b["h_out"] + H / c["k"] + 1.0 / b["h_in"])   # paper Eq. 4, plain slab
    for t in day:
        Tr, Tc = T_roof[t] - 273.15, T_ceil[t] - 273.15
        row = {
            "time_s": t,
            "hour": clock0 + (t - t0) / 3600.0,
            "clock_h": (clock0 + t / 3600.0) % 24.0,
            "T_out_C": f_out(t), "GHI_W_m2": f_ghi(t), "T_solair_C": f_sa(t),
            "T_indoor_C": f_in(t),
            "T_roof_C": Tr, "T_ceiling_C": Tc,
            "T_mean_C": T_mean[t] - 273.15 if t in T_mean else float("nan"),
            "q_roof_in_W_m2": h_out * (f_sa(t) - Tr),
            "q_ceiling_gain_W_m2": h_in * (Tc - f_in(t)),
            "q_paper_eq3_W_m2": k_eff * (Tr - Tc),
            "ref_T_roof_C": rT[t][0] if t in rT else float("nan"),
            "ref_T_ceiling_C": rT[t][1] if t in rT else float("nan"),
            "ref_q_ceiling_gain_W_m2": rq[t][1] if t in rq else float("nan"),
            "whf_roof_in_W_m2": whf[(t, "roof")] / area if (t, "roof") in whf else float("nan"),
            "whf_ceiling_gain_W_m2": (-whf[(t, "ceiling")] / area if (t, "ceiling") in whf
                                      else float("nan")),
        }
        rows.append(row)

    ts = [r["time_s"] for r in rows]
    hr = [r["hour"] for r in rows]
    q_in = [r["q_roof_in_W_m2"] for r in rows]
    q_g = [r["q_ceiling_gain_W_m2"] for r in rows]
    Wh = 1.0 / 3600.0
    roof_in = _trapz(ts, q_in) * Wh
    gross = _trapz(ts, [abs(q) for q in q_in]) * Wh
    gain = _trapz(ts, q_g) * Wh
    cooling = _trapz(ts, [max(q, 0.0) for q in q_g]) * Wh
    heating = -_trapz(ts, [min(q, 0.0) for q in q_g]) * Wh

    # ---- energy balance of the analysed day
    if day[0] in T_mean and day[-1] in T_mean:
        stored = c["rho"] * c["cp"] * H * (T_mean[day[-1]] - T_mean[day[0]]) * Wh
        resid = (roof_in - gain - stored) / gross if gross > 0 else float("nan")
        gates["energy"] = _gate(abs(resid) <= ENERGY_FRAC, resid, ENERGY_FRAC, "-",
                                "(roof in - ceiling out - stored) / gross roof heat flow")
    else:
        stored = None
        gates["energy"] = _gate(None, None, ENERGY_FRAC, "-", "no slabAverage output")

    # ---- daily metrics (paper Sec. 3.1-3.5 quantities for the plain slab)
    Tsa = [r["T_solair_C"] for r in rows]
    Tr = [r["T_roof_C"] for r in rows]
    Tc = [r["T_ceiling_C"] for r in rows]
    t_sa_max, sa_max = _argmax(hr, Tsa)
    t_c_max, c_max = _argmax(hr, Tc)
    t_c_min, c_min = _argmin(hr, Tc)
    t_r_max, r_max = _argmax(hr, Tr)
    t_r_min, r_min = _argmin(hr, Tr)
    t_g_max, g_max = _argmax(hr, q_g)
    t_g_min, g_min = _argmin(hr, q_g)
    metrics.update({
        "analysed_window_s": [day[0], day[-1]],
        "T_roof_C": {"max": r_max, "at_h": t_r_max % 24, "min": r_min, "at_h_min": t_r_min % 24},
        "T_ceiling_C": {"max": c_max, "at_h": t_c_max % 24, "min": c_min,
                        "at_h_min": t_c_min % 24, "mean": sum(Tc) / len(Tc)},
        "T_solair_C": {"max": sa_max, "at_h": t_sa_max % 24, "min": min(Tsa)},
        "decrement_factor": (c_max - c_min) / (sa_max - min(Tsa)),
        "time_lag_h": (t_c_max - t_sa_max) % 24.0,
        "q_ceiling_gain_W_m2": {"max": g_max, "at_h": t_g_max % 24, "min": g_min,
                                "at_h_min": t_g_min % 24, "mean": gain / 24.0},
        "daily_Wh_m2": {"roof_in": roof_in, "ceiling_gain_net": gain,
                        "cooling_load": cooling, "heating_load": heating,
                        "total_load": cooling + heating, "stored_change": stored},
    })
    wr = [abs(r["whf_roof_in_W_m2"] - r["q_roof_in_W_m2"]) for r in rows
          if r["whf_roof_in_W_m2"] == r["whf_roof_in_W_m2"]]
    wc = [abs(r["whf_ceiling_gain_W_m2"] - r["q_ceiling_gain_W_m2"]) for r in rows
          if r["whf_ceiling_gain_W_m2"] == r["whf_ceiling_gain_W_m2"]]
    metrics["wallHeatFlux_vs_bc_W_m2"] = {"roof_max_abs": max(wr) if wr else None,
                                          "ceiling_max_abs": max(wc) if wc else None,
                                          "note": "diagnostic, not a gate"}
    passed = all(g["pass"] for g in gates.values())
    return {"run": run_dir, "passed": passed, "gates": gates, "metrics": metrics,
            "rows": rows}


def _plot(res: dict, path: str) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    rows = res["rows"]
    h = [r["hour"] for r in rows]
    col = lambda k: [r[k] for r in rows]          # noqa: E731
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    a1.plot(h, col("T_solair_C"), color="#e6550d", lw=1.2, label="sol-air")
    a1.plot(h, col("T_out_C"), color="#636363", lw=1.0, ls="--", label="outdoor air")
    a1.plot(h, col("T_roof_C"), color="#b30000", lw=2.0, label="roof surface (CFD)")
    a1.plot(h, col("T_ceiling_C"), color="#08519c", lw=2.0, label="ceiling surface (CFD)")
    a1.plot(h, col("ref_T_ceiling_C"), color="k", lw=0.8, ls=":", label="ceiling, 1D reference")
    a1.plot(h, col("T_indoor_C"), color="#6baed6", lw=1.0, label="indoor ambient")
    a1.set_ylabel("Temperature (°C)")
    a1.legend(fontsize=8, ncol=2)
    a2.plot(h, col("q_roof_in_W_m2"), color="#b30000", lw=1.5, label="into the roof face")
    a2.plot(h, col("q_ceiling_gain_W_m2"), color="#08519c", lw=2.0, label="ceiling -> room (gain)")
    a2.axhline(0.0, color="k", lw=0.5)
    a2.set_ylabel("Heat flux (W/m²)")
    a2.set_xlabel("Clock time on the analysed day (h, IST)")
    a2.set_xticks(range(int(h[0]), int(h[0]) + 25, 3))
    a2.legend(fontsize=8)
    ok = "all gates pass" if res["passed"] else "GATES FAILED"
    a1.set_title(f"Phase 1 pure concrete slab, last simulated day ({ok})", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return True


def main(run_dir: str, plot: bool = True) -> int:
    res = evaluate(run_dir)
    out = os.path.join(run_dir, "results")
    os.makedirs(out, exist_ok=True)
    rows = res.pop("rows")
    with open(os.path.join(out, "timeseries_lastday.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(list(rows[0]))
        for r in rows:
            w.writerow([f"{v:.6g}" for v in r.values()])
    with open(os.path.join(out, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2)
    res["rows"] = rows
    plotted = plot and _plot(res, os.path.join(out, "quicklook.png"))

    print(f"Phase 1 check: {run_dir}")
    for name, g in res["gates"].items():
        state = {True: "PASS", False: "FAIL", None: "n/a "}[g["pass"]]
        val = "-" if g["value"] is None else f"{g['value']:.4g}"
        print(f"  [{state}] {name:<10} {val:>10} (limit {g['limit']:g} {g['unit']})  {g['note']}")
    m = res["metrics"]
    d = m["daily_Wh_m2"]
    print(f"  ceiling T {m['T_ceiling_C']['min']:.2f} .. {m['T_ceiling_C']['max']:.2f} C "
          f"(peak at {m['T_ceiling_C']['at_h']:.2f} h), roof T {m['T_roof_C']['min']:.2f} .. "
          f"{m['T_roof_C']['max']:.2f} C")
    print(f"  decrement factor {m['decrement_factor']:.3f}, time lag {m['time_lag_h']:.2f} h")
    print(f"  ceiling gain {m['q_ceiling_gain_W_m2']['min']:.1f} .. "
          f"{m['q_ceiling_gain_W_m2']['max']:.1f} W/m2; daily cooling load "
          f"{d['cooling_load']:.0f} Wh/m2, heating load {d['heating_load']:.0f} Wh/m2")
    print(f"  -> {out}/summary.json, timeseries_lastday.csv"
          + (", quicklook.png" if plotted else " (no matplotlib: no plot)"))
    print("PHASE 1 PASSED" if res["passed"] else "PHASE 1 NOT PASSED: see the gates above")
    return 0 if res["passed"] else 1
