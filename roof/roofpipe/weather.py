"""Design-day boundary conditions: outdoor air, solar irradiance, sol-air temperature
and the indoor ambient seen by the ceiling (config: site, weather, roof_surface, indoor).

One design day, repeated for every simulated day (a function of the clock time only):

  T_out(t)  piecewise cosine: T_min at sunrise, T_max at t_Tmax_h, back to T_min at
            the next sunrise
  I(t)      global horizontal irradiance. Shape: ASHRAE (1985) clear-sky model from
            the solar altitude, or a half sine between sunrise and sunset. Scaled so
            that the day integrates to daily_GHI_kWh_m2.
  T_sa(t) = T_out + a I / h_out - eps dR / h_out          (ASHRAE sol-air temperature)

The roof flux  q = a I - eps dR + h_out (T_out - T_s) = h_out (T_sa - T_s)  is then a
single convective condition (externalWallHeatFluxTemperature, mode coefficient).
Standard library only.
"""
from __future__ import annotations

import math

DAY_S = 86400.0

# ASHRAE (1985) clear-sky model, 21st of each month: apparent extraterrestrial
# irradiance A [W/m2], optical depth B [-], diffuse factor C [-]
_ASHRAE_ABC = (
    (1230.0, 0.142, 0.058), (1215.0, 0.144, 0.060), (1186.0, 0.156, 0.071),
    (1136.0, 0.180, 0.097), (1104.0, 0.196, 0.121), (1088.0, 0.205, 0.134),
    (1085.0, 0.207, 0.136), (1107.0, 0.201, 0.122), (1151.0, 0.177, 0.092),
    (1192.0, 0.160, 0.073), (1221.0, 0.149, 0.063), (1233.0, 0.142, 0.057),
)
_MONTH_END = (31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
INDOOR_MODES = ("conditioned", "ventilated", "specimen")


def month_index(day_of_year: int) -> int:
    """0 = January (non-leap year)."""
    for m, end in enumerate(_MONTH_END):
        if day_of_year <= end:
            return m
    return 11


class Sun:
    """Solar geometry of one day: Cooper declination, equation of time after Duffie &
    Beckman, local clock time (IST) of solar noon, sunrise and sunset."""

    def __init__(self, day_of_year: int, lat_deg: float, lon_deg: float, utc_offset_h: float):
        self.lat = math.radians(lat_deg)
        self.decl = math.radians(23.45) * math.sin(2.0 * math.pi * (284 + day_of_year) / 365.0)
        b = 2.0 * math.pi * (day_of_year - 81) / 364.0
        eot_min = 9.87 * math.sin(2.0 * b) - 7.53 * math.cos(b) - 1.5 * math.sin(b)
        self.noon_h = 12.0 - (lon_deg - 15.0 * utc_offset_h) / 15.0 - eot_min / 60.0
        x = -math.tan(self.lat) * math.tan(self.decl)
        half_day_h = math.degrees(math.acos(max(-1.0, min(1.0, x)))) / 15.0
        self.sunrise_h = self.noon_h - half_day_h
        self.sunset_h = self.noon_h + half_day_h

    def sin_altitude(self, clock_h: float) -> float:
        w = math.radians(15.0 * (clock_h - self.noon_h))
        return (math.sin(self.lat) * math.sin(self.decl)
                + math.cos(self.lat) * math.cos(self.decl) * math.cos(w))


class DesignDay:
    """Every boundary history as a function of the clock hour (any real number; the
    day repeats)."""

    def __init__(self, cfg: dict):
        site, w = cfg["site"], cfg["weather"]
        doy = int(site["day_of_year"])
        self.sun = Sun(doy, float(site["lat_deg"]), float(site["lon_deg"]),
                       float(site["utc_offset_h"]))
        self.T_max = float(w["T_max_C"])
        self.T_min = float(w["T_min_C"])
        self.t_max_h = float(w["t_Tmax_h"])
        if not self.sun.sunrise_h < self.t_max_h < self.sun.sunrise_h + 24.0:
            raise ValueError(f"weather.t_Tmax_h = {self.t_max_h} must be after sunrise "
                             f"({self.sun.sunrise_h:.2f} h)")
        self.shape = w.get("irradiance_shape", "ashrae_clear_sky")
        if self.shape not in ("ashrae_clear_sky", "half_sine"):
            raise ValueError(f"weather.irradiance_shape: unknown shape {self.shape!r}")
        self.abc = _ASHRAE_ABC[month_index(doy)]
        n = 1440                                    # 1-minute midpoint rule over the day
        daily_Wh = sum(self._shape(24.0 * (i + 0.5) / n) for i in range(n)) * 24.0 / n
        target_Wh = 1000.0 * float(w["daily_GHI_kWh_m2"])
        self.ghi_scale = target_Wh / daily_Wh if daily_Wh > 0.0 else 0.0

        r = cfg["roof_surface"]
        self.absorptance = float(r["absorptance"])
        self.h_out = float(r["h_out"])
        self.longwave_K = float(r["longwave_correction_K"])

        ind = cfg["indoor"]
        self.indoor_mode = ind["mode"]
        if self.indoor_mode not in INDOOR_MODES:
            raise ValueError(f"indoor.mode must be one of {INDOOR_MODES}, not {self.indoor_mode!r}")
        self.h_in = float(ind["h_in"])
        self.T_set = float(ind["T_set_C"])
        self.T_offset = float(ind.get("T_offset_K", 0.0))
        if self.indoor_mode == "specimen":
            ins = ind["insulation"]
            self.h_in_eff = 1.0 / (1.0 / self.h_in + float(ins["thickness"]) / float(ins["k"]))
        else:
            self.h_in_eff = self.h_in

    # ------------------------------------------------------------------ outdoor
    def _shape(self, clock_h: float) -> float:
        h = clock_h % 24.0
        if self.shape == "half_sine":
            day = self.sun.sunset_h - self.sun.sunrise_h
            x = h - self.sun.sunrise_h
            return math.sin(math.pi * x / day) if 0.0 < x < day else 0.0
        s = self.sun.sin_altitude(h)
        if s <= 0.0:
            return 0.0
        A, B, C = self.abc
        dni = A * math.exp(-B / s)
        return dni * (s + C)                        # beam on the horizontal + diffuse

    def ghi(self, clock_h: float) -> float:
        return self.ghi_scale * self._shape(clock_h)

    def T_out(self, clock_h: float) -> float:
        h = clock_h % 24.0
        rise, hi = self.sun.sunrise_h, self.t_max_h
        mean = 0.5 * (self.T_max + self.T_min)
        amp = 0.5 * (self.T_max - self.T_min)
        if rise <= h <= hi:
            return mean - amp * math.cos(math.pi * (h - rise) / (hi - rise))
        since_max = (h - hi) % 24.0
        return mean + amp * math.cos(math.pi * since_max / (24.0 - (hi - rise)))

    def T_solair(self, clock_h: float) -> float:
        return (self.T_out(clock_h) + self.absorptance * self.ghi(clock_h) / self.h_out
                - self.longwave_K)

    # ------------------------------------------------------------------- indoor
    def T_indoor(self, clock_h: float) -> float:
        """Ambient behind the ceiling film: room air, or shaded outdoor air (specimen)."""
        if self.indoor_mode == "conditioned":
            return self.T_set
        if self.indoor_mode == "ventilated":
            return self.T_out(clock_h) + self.T_offset
        return self.T_out(clock_h)


def tables(cfg: dict):
    """Boundary histories on the run's time axis (t = 0 at run.start_clock_h), every
    run.table_step_min, end point included. Returns (DesignDay, dict of lists)."""
    dd = DesignDay(cfg)
    run = cfg["run"]
    step = 60.0 * float(run["table_step_min"])
    end = DAY_S * float(run["n_days"])
    n = int(round(end / step))
    if abs(n * step - end) > 1e-6:
        raise ValueError("run.n_days x 86400 s must be a multiple of run.table_step_min")
    start = float(run.get("start_clock_h", 0.0))
    out = {k: [] for k in ("time_s", "clock_h", "T_out_C", "GHI_W_m2", "T_solair_C",
                           "T_indoor_C")}
    for i in range(n + 1):
        t = i * step
        c = (start + t / 3600.0) % 24.0
        out["time_s"].append(t)
        out["clock_h"].append(c)
        out["T_out_C"].append(dd.T_out(c))
        out["GHI_W_m2"].append(dd.ghi(c))
        out["T_solair_C"].append(dd.T_solair(c))
        out["T_indoor_C"].append(dd.T_indoor(c))
    return dd, out


def foam_rows(time_s, values_C, indent: str = "") -> str:
    """Rows of an inline OpenFOAM Function1 table, temperatures in K."""
    return "\n".join(f"{indent}({t:.10g} {v + 273.15:.4f})" for t, v in zip(time_s, values_C))
