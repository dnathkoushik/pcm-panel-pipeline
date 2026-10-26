"""Thermal loads for a roof installation (config/loads.yaml).

  outdoor air T_out(t) and global horizontal irradiance I(t)
    -> sol-air temperature  T_sa = T_out + a I / h_out - eps dR / h_out   (ASHRAE)
    -> 1D transient conduction through the roof layers, h_out to T_sa outside,
       h_in to T_room inside, repeated until the day is periodic
    -> ceiling (slab underside) temperature T_ceiling(t)

The slab is solved here, not in OpenFOAM: its thermal mass (125 mm concrete,
time constant ~5 h) sets the phase and amplitude of the ceiling temperature,
and a 0.25-0.75 L panel is far too small to change it.

Weather is either `synthetic` (daily max/min + daily irradiation: piecewise
cosine temperature with the minimum at sunrise, half-sine irradiance between
sunrise and sunset) or the hottest day of a month in an EPW file.
"""
from __future__ import annotations

import csv
import math
import os

import numpy as np

DAY_S = 86400.0


def solar_times(doy: int, lat_deg: float, lon_deg: float, utc_offset_h: float):
    """Sunrise, sunset and solar noon as local clock hours."""
    decl = math.radians(23.45) * math.sin(2 * math.pi * (284 + doy) / 365.0)   # Cooper
    B = 2 * math.pi * (doy - 81) / 364.0
    eot_min = 9.87 * math.sin(2 * B) - 7.53 * math.cos(B) - 1.5 * math.sin(B)
    noon = 12.0 - (lon_deg - 15.0 * utc_offset_h) / 15.0 - eot_min / 60.0
    ws = math.degrees(math.acos(-math.tan(math.radians(lat_deg)) * math.tan(decl)))
    return noon - ws / 15.0, noon + ws / 15.0, noon


def _synthetic(spec: dict, clock_h: np.ndarray):
    w, loc = spec["weather"], spec["location"]
    rise, sset, _ = solar_times(spec["day_of_year"], loc["lat_deg"], loc["lon_deg"],
                                loc["utc_offset_h"])
    mean = 0.5 * (w["T_max_C"] + w["T_min_C"])
    amp = 0.5 * (w["T_max_C"] - w["T_min_C"])
    t_hi = w["t_Tmax_h"]
    h = np.mod(clock_h - rise, 24.0)              # hours since sunrise
    rising = h <= t_hi - rise
    T = np.where(rising, mean - amp * np.cos(np.pi * h / (t_hi - rise)),
                 mean + amp * np.cos(np.pi * (h - (t_hi - rise)) / (24.0 - (t_hi - rise))))
    N = sset - rise
    peak = math.pi * 1000.0 * w["daily_global_horizontal_kWh_m2"] / (2.0 * N)   # W/m2
    s = np.mod(clock_h, 24.0) - rise
    I = np.where((s > 0) & (s < N), peak * np.sin(np.pi * np.clip(s, 0, N) / N), 0.0)
    return T, I


def read_epw_day(path: str, month: int):
    """Hourly dry-bulb [C] and global horizontal irradiance [W/m2] of the hottest
    day of `month`. EPW hour k covers (k-1, k]: its value is placed at k - 0.5."""
    days = {}
    with open(path, newline="", encoding="latin-1") as f:
        for i, row in enumerate(csv.reader(f)):
            if i < 8 or len(row) < 14:
                continue
            if int(row[1]) == month:
                days.setdefault(int(row[2]), []).append((int(row[3]), float(row[6]), float(row[13])))
    if not days:
        raise ValueError(f"no data for month {month} in {path}")
    day = max(days, key=lambda d: max(r[1] for r in days[d]))
    rows = sorted(days[day])
    return day, np.array([r[0] - 0.5 for r in rows]), np.array([r[1] for r in rows]), \
        np.array([r[2] for r in rows])


def weather(spec: dict, clock_h: np.ndarray, root: str = "."):
    w = spec["weather"]
    if w.get("source", "synthetic") == "synthetic":
        return _synthetic(spec, clock_h)
    _, hrs, T, I = read_epw_day(os.path.join(root, w["epw_file"]), w.get("epw_month", 5))
    x = np.mod(clock_h, 24.0)
    return (np.interp(x, hrs, T, period=24.0), np.interp(x, hrs, I, period=24.0))


def room_temperature_C(spec: dict, clock_h, root: str = "."):
    T_out, _ = weather(spec, np.asarray(clock_h, dtype=float), root)
    return T_out + spec["room"].get("T_offset_K", 0.0)


def slab_periodic(clock_h, T_sa, T_room, layers, h_out, h_in, dt=60.0, cells_per_layer=20,
                  max_days=20, tol=1e-3):
    """Implicit 1D conduction, repeated over the (periodic) day until two days agree.
    Returns the inner surface temperature and the heat flux into the room [W/m2]."""
    k, rc, dx = [], [], []
    for L in layers:
        n = cells_per_layer
        k += [L["k"]] * n
        rc += [L["rho"] * L["cp"]] * n
        dx += [L["thickness"] / n] * n
    k, rc, dx = map(np.array, (k, rc, dx))
    n = len(k)
    R = np.r_[0.5 * dx[0] / k[0] + 1.0 / h_out,                      # T_sa -> cell 0
              0.5 * dx[:-1] / k[:-1] + 0.5 * dx[1:] / k[1:],        # between cells
              0.5 * dx[-1] / k[-1] + 1.0 / h_in]                     # cell n-1 -> T_room
    G = 1.0 / R
    C = rc * dx / dt
    M = np.diag(C + G[:-1] + G[1:]) - np.diag(G[1:-1], 1) - np.diag(G[1:-1], -1)
    Minv = np.linalg.inv(M)
    steps = int(round(DAY_S / dt))
    t_day = np.arange(steps) * dt / 3600.0
    sa = np.interp(t_day, clock_h, T_sa, period=24.0)
    rm = np.interp(t_day, clock_h, T_room, period=24.0)
    T = np.full(n, float(np.mean(T_sa)))
    prev = None
    for _ in range(max_days):
        Tin = np.empty(steps)
        for i in range(steps):
            b = C * T
            b[0] += G[0] * sa[i]
            b[-1] += G[-1] * rm[i]
            T = Minv @ b
            Tin[i] = T[-1]
        if prev is not None and np.abs(Tin - prev).max() < tol:
            break
        prev = Tin
    # inner surface from the last half cell and the film
    k_half = k[-1] / (0.5 * dx[-1])
    Ts = (k_half * Tin + h_in * rm) / (k_half + h_in)
    q = h_in * (Ts - rm)
    clock = np.asarray(clock_h, dtype=float)
    return (np.interp(np.mod(clock, 24.0), t_day, Ts, period=24.0),
            np.interp(np.mod(clock, 24.0), t_day, q, period=24.0))


def profile(spec: dict, root: str = ".") -> dict:
    """Boundary histories on the CFD time axis (t = 0 at start_clock_h)."""
    step = 60.0 * spec.get("table_step_min", 10)
    t = np.arange(0.0, 3600.0 * spec.get("duration_h", 24.0) + 0.5 * step, step)
    clock = spec["start_clock_h"] + t / 3600.0
    T_out, I = weather(spec, clock, root)
    r = spec["roof"]
    T_sa = T_out + r["absorptance"] * I / r["h_out"] - r["longwave_correction_K"]
    T_room = T_out + spec["room"].get("T_offset_K", 0.0)
    # the slab sees a fine, full day (not only the CFD window)
    fine = np.arange(0.0, 24.0, 0.25)
    To_f, I_f = weather(spec, fine, root)
    sa_f = To_f + r["absorptance"] * I_f / r["h_out"] - r["longwave_correction_K"]
    rm_f = To_f + spec["room"].get("T_offset_K", 0.0)
    T_c, q = slab_periodic(fine, sa_f, rm_f, r["layers"], r["h_out"], r["h_in"])
    T_ceiling = np.interp(np.mod(clock, 24.0), fine, T_c, period=24.0)
    q_ceiling = np.interp(np.mod(clock, 24.0), fine, q, period=24.0)
    return {"time_s": t, "clock_h": clock, "T_out_C": T_out, "I_global_W_m2": I,
            "T_solair_C": T_sa, "T_ceiling_C": T_ceiling, "T_room_C": T_room,
            "q_ceiling_W_m2": q_ceiling}


def foam_table(time_s, values_K) -> str:
    """Inline Function1 table for uniformFixedValue / externalWallHeatFluxTemperature."""
    rows = " ".join(f"({t:g} {v:.3f})" for t, v in zip(time_s, values_K))
    return f"table ({rows})"


def recharge_hours(p: dict, T_solidus_K: float) -> float:
    """Hours per day the room air is below the PCM solidus (night-time recharge)."""
    fine = np.arange(0.0, 24.0, 0.1)
    T_room = np.interp(fine, np.mod(p["clock_h"], 24.0), p["T_room_C"], period=24.0)
    return float(0.1 * np.sum(T_room + 273.15 < T_solidus_K))


def write_profile(case: dict, out_dir: str, root: str = ".") -> dict:
    spec = case["loads_spec"]
    p = profile(spec, root)
    os.makedirs(out_dir, exist_ok=True)
    keys = list(p)
    with open(os.path.join(out_dir, "loads_profile.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for row in zip(*(p[k] for k in keys)):
            w.writerow([f"{v:.4g}" for v in row])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pcm = case["physics"]["pcm"]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    h = p["clock_h"]
    offset = spec["room"].get("T_offset_K", 0.0)
    ax.plot(h, p["T_solair_C"], color="#e6550d", lw=1.2, label="sol-air (roof top)")
    if offset:
        ax.plot(h, p["T_out_C"], color="#636363", lw=1.2, ls="--", label="outdoor air")
    ax.plot(h, p["T_ceiling_C"], color="#b30000", lw=2.0, label="ceiling = slab underside (CFD)")
    ax.plot(h, p["T_room_C"], color="#3182bd", lw=2.0,
            label="room walls/floor (CFD)" + ("" if offset else " = outdoor air"))
    ax.axhspan(pcm["T_solidus"] - 273.15, pcm["T_liquidus"] - 273.15, color="#31a354", alpha=0.25,
               label=f"PCM melting range ({pcm['T_solidus'] - 273.15:.1f}-{pcm['T_liquidus'] - 273.15:.1f} °C)")
    ax.set_xlabel("Clock time (h, IST)")
    ax.set_ylabel("Temperature (°C)")
    ticks = np.arange(math.floor(h[0]), math.ceil(h[-1]) + 1, 3)
    ax.set_xticks(ticks, [f"{int(x) % 24:02d}:00" for x in ticks])
    ax.set_xlim(h[0], h[-1])
    ax2 = ax.twinx()
    ax2.fill_between(h, p["I_global_W_m2"], color="#fdae6b", alpha=0.25, lw=0)
    ax2.set_ylabel("Global horizontal irradiance (W/m²)")
    ax2.set_ylim(0, 3.0 * max(1.0, p["I_global_W_m2"].max()))
    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title(f"{spec['location']['name']}, day {spec['day_of_year']}: roof boundary conditions",
                 fontsize=10)
    fig.savefig(os.path.join(out_dir, "loads_profile.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    hrs = recharge_hours(p, pcm["T_solidus"])
    print(f"  loads: ceiling {p['T_ceiling_C'].min():.1f}-{p['T_ceiling_C'].max():.1f} C, "
          f"room {p['T_room_C'].min():.1f}-{p['T_room_C'].max():.1f} C, "
          f"sol-air max {p['T_solair_C'].max():.1f} C -> "
          f"{os.path.join(out_dir, 'loads_profile.png')}")
    if hrs < 6.0:
        print(f"  WARNING: room air is below the PCM solidus ({pcm['T_solidus'] - 273.15:.1f} C) "
              f"for only {hrs:.1f} h/day: the PCM will not re-solidify overnight without "
              f"night ventilation or cooling.")
    return p
