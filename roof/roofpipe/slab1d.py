"""1D finite-volume reference for the concrete slab.

The same discretisation as the OpenFOAM baseline: nz uniform cells through H,
implicit Euler with the same time step, the same boundary tables (linear
interpolation, as OpenFOAM's Function1 table) and the same film coefficients.
At a convective face the surface temperature is, as in OpenFOAM's mixed condition,

    T_s = (h T_a + (k/d) T_c) / (h + k/d),      d = dz/2 (cell centre to face)

so the 3D case, which is one-dimensional by symmetry, must reproduce this solution
to linear-solver tolerance. Pure Python: a tridiagonal solve per step, factorised once.

Sign convention (used throughout roofpipe):
    q_roof_in        heat flux into the slab through the roof face        [W/m2]
    q_ceiling_gain   heat flux out of the slab into the room (heat gain)  [W/m2]
"""
from __future__ import annotations

import bisect


def table_function(times, values):
    """Piecewise-linear interpolation, clamped at both ends."""
    times, values = list(times), list(values)

    def f(t: float) -> float:
        if t <= times[0]:
            return values[0]
        if t >= times[-1]:
            return values[-1]
        i = bisect.bisect_right(times, t) - 1
        w = (t - times[i]) / (times[i + 1] - times[i])
        return values[i] + w * (values[i + 1] - values[i])

    return f


def solve(H, nz, k, rho, cp, h_top, Ta_top, h_bot, Ta_bot, T0, dt, t_end, sample_dt):
    """Ta_top, Ta_bot: functions of time [s]. Temperatures in C (the problem is linear).
    Index 0 is the roof (top) cell, nz-1 the ceiling (bottom) cell."""
    if nz < 2:
        raise ValueError("need at least 2 cells through the slab")
    dz = H / nz
    d = 0.5 * dz
    C = rho * cp * dz / dt                      # per unit area
    G = k / dz                                  # cell to cell
    G_top = 1.0 / (d / k + 1.0 / h_top)         # top cell to the outdoor ambient
    G_bot = 1.0 / (d / k + 1.0 / h_bot)         # bottom cell to the indoor ambient
    k_d = k / d

    # tridiagonal matrix a T[i-1] + b T[i] + c T[i+1] = rhs, constant in time
    a = [-G] * nz
    c = [-G] * nz
    b = [C + 2.0 * G] * nz
    a[0] = 0.0
    c[-1] = 0.0
    b[0] = C + G + G_top
    b[-1] = C + G + G_bot
    cp_ = [0.0] * nz                             # Thomas: modified upper diagonal
    den = [0.0] * nz
    den[0] = b[0]
    cp_[0] = c[0] / den[0]
    for i in range(1, nz):
        den[i] = b[i] - a[i] * cp_[i - 1]
        cp_[i] = c[i] / den[i]

    n_steps = int(round(t_end / dt))
    every = int(round(sample_dt / dt))
    if every < 1 or abs(every * dt - sample_dt) > 1e-9 * sample_dt:
        raise ValueError("sample_dt must be a whole number of time steps")
    T = [float(T0)] * nz
    out = {k_: [] for k_ in ("time_s", "T_roof_C", "T_ceiling_C", "T_mean_C",
                             "q_roof_in_W_m2", "q_ceiling_gain_W_m2")}

    def record(t, Tt, Tb):
        Ts_top = (h_top * Tt + k_d * T[0]) / (h_top + k_d)
        Ts_bot = (h_bot * Tb + k_d * T[-1]) / (h_bot + k_d)
        out["time_s"].append(t)
        out["T_roof_C"].append(Ts_top)
        out["T_ceiling_C"].append(Ts_bot)
        out["T_mean_C"].append(sum(T) / nz)
        out["q_roof_in_W_m2"].append(h_top * (Tt - Ts_top))
        out["q_ceiling_gain_W_m2"].append(h_bot * (Ts_bot - Tb))

    record(0.0, Ta_top(0.0), Ta_bot(0.0))
    dp = [0.0] * nz
    for n in range(1, n_steps + 1):
        t = n * dt
        Tt, Tb = Ta_top(t), Ta_bot(t)           # implicit Euler: boundary values at t_n+1
        rhs0 = C * T[0] + G_top * Tt
        dp[0] = rhs0 / den[0]
        for i in range(1, nz):
            r = C * T[i]
            if i == nz - 1:
                r += G_bot * Tb
            dp[i] = (r - a[i] * dp[i - 1]) / den[i]
        T[-1] = dp[-1]
        for i in range(nz - 2, -1, -1):
            T[i] = dp[i] - cp_[i] * T[i + 1]
        if n % every == 0:
            record(t, Tt, Tb)
    return out
