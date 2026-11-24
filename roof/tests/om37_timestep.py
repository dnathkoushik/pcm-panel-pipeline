"""Time-step check for OM37 in the roof slab (pcm_study.run.deltaT).

Runs the repository's 1D mirror of pcmEnthalpyPorositySource (tests/test_latent_source_1d.py:
the same Newton linearisation and kink limiter as the C++ fvOption) with the OM37 properties
of roof/config/roof.yaml, and compares it with the fully converged enthalpy method at a
5 s step. Case: 50 mm of solid OM37 at 30 C, one face raised to 50 C (the roof-side
extreme), 6 h, 2.5 mm cells as in the Phase 2 mesh.

    python3 roof/tests/om37_timestep.py          (numpy; e.g. inside the Docker image)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tests"))
import test_latent_source_1d as t1d  # noqa: E402

OM37 = dict(TS=36.4 + 273.15, TL=44.06 + 273.15, LAT=179000.0, CPS=2550.0, CPL=2630.0,
            KS=0.16, KL=0.13, RHO=860.0, CPREF=2630.0, KREF=0.145)
CASE = dict(n=20, L=0.05, t_end=6 * 3600.0, T0=303.15, Tw=323.15)


def main():
    rows = []
    with t1d.material(**OM37):
        ref = t1d.run_reference(dt=5.0, **CASE)
        for dt in (10.0, 30.0, 60.0, 120.0):
            for n_outer in (2, 3):
                fvo = t1d.run_fvoption(dt=dt, n_outer=n_outer, **CASE)
                # compare on the coarse run's times
                g_ref = np.interp(fvo[:, 0], ref[:, 0], ref[:, 1])
                e_ref = np.interp(fvo[:, 0], ref[:, 0], ref[:, 2])
                rows.append((dt, n_outer, np.abs(fvo[:, 1] - g_ref).max(),
                             abs(fvo[-1, 2] - e_ref[-1]) / e_ref[-1], fvo[-1, 1]))
    print(f"reference (converged, dt 5 s): liquid fraction after 6 h {ref[-1, 1]:.4f}, "
          f"stored {ref[-1, 2] / 1e3:.1f} kJ/m2")
    print("  dt [s]  outer  max |d gamma|  stored-energy error")
    for dt, n, dg, de, g in rows:
        print(f"  {dt:6.0f}  {n:5d}  {dg:13.4f}  {100 * de:+18.3f} %")


if __name__ == "__main__":
    main()
