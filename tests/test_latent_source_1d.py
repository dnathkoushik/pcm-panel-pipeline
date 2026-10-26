"""1D verification of the pcmEnthalpyPorositySource energy formulation.

Mirrors the C++ source: the solver advances h = CpRef (T - Tstd) with constant
CpRef/kappaRef and the fvOption adds
    -(rho E^k - rho E0)/dt + Sp h^k - Sp h + div((k(T) - kappaRef) grad T),
    Sp = rho dE/dT / (CpRef dt),
iterated nOuter times per step (PIMPLE outer correctors). Reference: the
fully converged enthalpy method rho (H(T) - H0)/dt = div(k(T) grad T),
H = int cp dT + gamma L (paper Eqs. 1, 15, 18-20), solved by Newton to 1e-10.
Slab of PCM heated from one side at 298.15 K, adiabatic on the other.

The fvOption also has a kink limiter (correct() in the C++ source): an outer
iteration may not carry a cell across Tsolidus or Tliquidus; it stops just past
the kink and is re-linearised. Without it Newton oscillates for pure PCMs with
a ~1 K melting range (n-octadecane), see the Octadecane tests.
"""
import contextlib
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

TS, TL, LAT = 286.0, 293.0, 115000.0
CPS, CPL, KS, KL, RHO = 1830.0, 2210.0, 0.183, 0.169, 840.0
CPREF, KREF, TSTD = CPL, KL, 298.15


def gamma(T):
    return np.clip((T - TS) / (TL - TS), 0.0, 1.0)


def cp_true(T):
    return CPS + (CPL - CPS) * gamma(T)


def k_true(T):
    return KS + (KL - KS) * gamma(T)


def E(T):
    """Same piecewise formula as Ecell() in the C++ source."""
    dT = TL - TS
    g = gamma(T)
    liq = (CPL - CPREF) * (T - TL) + LAT
    mush = (CPREF - CPS) * (TL - T) - 0.5 * (CPL - CPS) * dT * (1 - g**2) + g * LAT
    sol = CPREF * (TL - T) - CPS * (TS - T) - 0.5 * (CPS + CPL) * dT
    return np.where(T >= TL, liq, np.where(T > TS, mush, sol))


def dEdT(T):
    return cp_true(T) - CPREF + np.where((T > TS) & (T < TL), LAT / (TL - TS), 0.0)


def H(T):
    """Total enthalpy relative to Tstd: CpRef (T - Tstd) + E(T) (+ const)."""
    return CPREF * (T - TSTD) + E(T)


def kink_limit(T0, Tn):
    """Same rule as pcmEnthalpyPorositySource::correct()."""
    tol = 1e-3 * (TL - TS)
    for kink in (TS, TL):
        up = (T0 < kink - tol) & (Tn > kink)
        down = (T0 > kink + tol) & (Tn < kink)
        Tn = np.where(up, kink + tol, np.where(down, kink - tol, Tn))
    return Tn


@contextlib.contextmanager
def material(**props):
    """Temporarily swap the module-level PCM properties (TS, TL, LAT, ...)."""
    saved = {k: globals()[k] for k in props}
    globals().update(props)
    try:
        yield
    finally:
        globals().update(saved)


def diffusion_matrix(k_face, dx, n, T_wall):
    """FV operator div(k grad T) with Dirichlet at x=0 (half cell) and adiabatic end."""
    A = np.zeros((n, n))
    b = np.zeros(n)
    for i in range(n - 1):
        c = k_face[i + 1] / dx**2
        A[i, i] -= c
        A[i, i + 1] += c
        A[i + 1, i + 1] -= c
        A[i + 1, i] += c
    c0 = KREF / (0.5 * dx) / dx        # wall flux with kappaRef, as in the coupled BC
    A[0, 0] -= c0
    b[0] += c0 * T_wall
    return A, b


def run_fvoption(n=40, L=0.01, dt=0.5, t_end=3600.0, n_outer=3, T0=278.15, Tw=298.15,
                 limiter=True):
    dx = L / n
    T = np.full(n, T0)
    E0 = E(T)
    hist = []
    for step in range(int(t_end / dt)):
        T_old = T.copy()
        h_old = CPREF * (T_old - TSTD)
        E0 = E(T_old)                                   # stored at the first call of a step
        for _ in range(n_outer):
            Ek, sp = E(T), RHO * dEdT(T) / (CPREF * dt)
            hk = CPREF * (T - TSTD)
            kex = k_true(T) - KREF
            kex_face = np.r_[0.0, 0.5 * (kex[1:] + kex[:-1])]
            Aex, _ = diffusion_matrix(kex_face, dx, n, 0.0)
            Aex[0, 0] += KREF / (0.5 * dx) / dx            # no extra wall flux (kEx = 0 there)
            Aref, bref = diffusion_matrix(np.full(n, KREF), dx, n, Tw)
            # matrix for h: rho/dt h - div(kRef/CpRef grad h) + Sp h = rhs
            M = np.diag(RHO / dt + sp) - Aref / CPREF
            rhs = (RHO / dt * h_old
                   + bref - Aref @ np.full(n, TSTD) * 0.0    # Dirichlet term in T units
                   - RHO * (Ek - E0) / dt + sp * hk
                   + Aex @ T)
            # Aref acts on T = Tstd + h/CpRef: the Tstd part has zero row sums except the wall row
            rhs += Aref @ np.full(n, TSTD)
            h = np.linalg.solve(M, rhs)
            Tn = TSTD + h / CPREF
            T = kink_limit(T, Tn) if limiter else Tn
        hist.append((dt * (step + 1), gamma(T).mean(), RHO * dx * (H(T) - H(np.full(n, T0))).sum()))
    return np.array(hist)


def run_reference(n=40, L=0.01, dt=0.5, t_end=3600.0, T0=278.15, Tw=298.15):
    dx = L / n
    T = np.full(n, T0)
    hist = []
    for step in range(int(t_end / dt)):
        H_old = H(T)
        for _ in range(100):                            # Newton on rho (H(T) - H0)/dt = div(k grad T)
            kf = k_true(T)
            k_face = np.r_[0.0, 0.5 * (kf[1:] + kf[:-1])]
            A, b = diffusion_matrix(k_face, dx, n, Tw)
            R = RHO * (H(T) - H_old) / dt - (A @ T + b)
            J = np.diag(RHO * (CPREF + dEdT(T)) / dt) - A
            dTn = np.linalg.solve(J, -R)
            T = T + dTn
            if np.abs(dTn).max() < 1e-10:
                break
        hist.append((dt * (step + 1), gamma(T).mean(), RHO * dx * (H(T) - H(np.full(n, T0))).sum()))
    return np.array(hist)


def run_explicit_reference(n=40, L=0.01, dt=0.05, every=0.5, t_end=3600.0, T0=278.15, Tw=298.15):
    """Explicit enthalpy method (conservative by construction, no iteration),
    for melting ranges too narrow for Newton without a limiter."""
    dx = L / n
    Tg = np.linspace(min(T0, Tw) - 30.0, max(T0, Tw) + 30.0, 200001)
    Hg = H(Tg)
    T = np.full(n, T0)
    Hc = H(T)
    hist, nsub = [], int(round(every / dt))
    for step in range(int(t_end / every)):
        for _ in range(nsub):
            kf = k_true(T)
            q = np.zeros(n + 1)
            q[1:-1] = -0.5 * (kf[1:] + kf[:-1]) * (T[1:] - T[:-1]) / dx
            q[0] = -KREF * (T[0] - Tw) / (0.5 * dx)
            Hc = Hc - dt * (q[1:] - q[:-1]) / dx / RHO
            T = np.interp(Hc, Hg, Tg)
        hist.append((every * (step + 1), gamma(T).mean(), RHO * dx * (Hc - H(np.full(n, T0))).sum()))
    return np.array(hist)


# n-octadecane, Faden et al. (2019) Materials 12, 2974 (config/materials.yaml)
OCTADECANE = dict(TS=300.63, TL=301.63, LAT=236980.0, CPS=1921.0, CPL=2230.0,
                  KS=0.334, KL=0.152, RHO=776.9, CPREF=2230.0, KREF=0.152)


class LatentSource(unittest.TestCase):
    def test_E_continuous_and_derivative(self):
        for Tk in (TS, TL):
            self.assertAlmostEqual(float(E(Tk - 1e-9)), float(E(Tk + 1e-9)), places=3)
        T = np.linspace(275, 300, 2001)
        T = T[(np.abs(T - TS) > 1e-3) & (np.abs(T - TL) > 1e-3)]
        fd = (E(T + 1e-6) - E(T - 1e-6)) / 2e-6
        np.testing.assert_allclose(fd, dEdT(T), rtol=1e-5, atol=1e-3)

    def test_total_enthalpy_matches_paper_definition(self):
        # H(T2) - H(T1) = int cp dT + L across the full melt (Eqs. 18-20)
        T1, T2 = 278.15, 298.15
        exact = CPS * (TS - T1) + 0.5 * (CPS + CPL) * (TL - TS) + CPL * (T2 - TL) + LAT
        self.assertAlmostEqual(float(H(T2) - H(T1)), exact, places=6)

    def test_fvoption_matches_converged_enthalpy_method(self):
        ref = run_reference()
        fvo = run_fvoption(n_outer=3)
        self.assertLess(np.abs(fvo[:, 1] - ref[:, 1]).max(), 0.01)        # liquid fraction
        self.assertLess(abs(fvo[-1, 2] - ref[-1, 2]) / ref[-1, 2], 0.005)  # stored energy

    def test_single_outer_corrector_still_close(self):
        ref = run_reference(t_end=1800.0)
        fvo = run_fvoption(t_end=1800.0, n_outer=1)
        self.assertLess(np.abs(fvo[:, 1] - ref[:, 1]).max(), 0.03)

    def test_limiter_inactive_in_the_band(self):
        T0 = np.array([TS + 1.0, TS + 3.0, TL + 2.0, TS - 2.0])
        Tn = np.array([TS + 2.0, TS + 4.0, TL + 5.0, TS - 1.0])
        np.testing.assert_array_equal(kink_limit(T0, Tn), Tn)


class Octadecane(unittest.TestCase):
    """Pure PCM with a 1 K melting range, heated at 40 C from 24 C."""

    @classmethod
    def setUpClass(cls):
        with material(**OCTADECANE):
            cls.ref = run_explicit_reference(T0=297.15, Tw=313.15)

    def _run(self, **kw):
        with material(**OCTADECANE):
            return run_fvoption(T0=297.15, Tw=313.15, **kw)

    def test_limiter_matches_reference(self):
        fvo = self._run(dt=0.5, n_outer=3)
        self.assertLess(np.abs(fvo[:, 1] - self.ref[:, 1]).max(), 0.05)
        self.assertLess(abs(fvo[-1, 2] - self.ref[-1, 2]) / self.ref[-1, 2], 0.005)

    def test_plain_newton_fails(self):
        # documents why the limiter exists: without it the stored energy is badly wrong
        fvo = self._run(dt=0.5, n_outer=3, limiter=False)
        self.assertGreater(abs(fvo[-1, 2] - self.ref[-1, 2]) / self.ref[-1, 2], 0.05)


if __name__ == "__main__":
    unittest.main(verbosity=2)
