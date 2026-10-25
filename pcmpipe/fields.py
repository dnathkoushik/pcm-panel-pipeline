"""0/<region>/ fields: initial (Eqs. 25-27) and boundary conditions (Sec. 3.4)."""
from __future__ import annotations

import os

from . import foam
from .foam import num

ROOM_WALLS = '"(floor|ceiling|walls)"'


def _field(path, cls, dims, internal, boundary):
    body = f"dimensions      {dims};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for patch, bc in boundary:
        if any(c in patch for c in ".*|()[]") and not patch.startswith('"'):
            patch = f'"{patch}"'                 # regex patch names must be quoted
        body += f"    {patch}\n    {{\n"
        for line in bc.strip().splitlines():
            body += f"        {line.strip()}\n"
        body += "    }\n"
    body += "}\n"
    foam.write(path, body, cls=cls, location=os.path.basename(os.path.dirname(path)))


def _coupled(T, method):
    return f"""
type            compressible::turbulentTemperatureCoupledBaffleMixed;
Tnbr            T;
kappaMethod     {method};
value           uniform {num(T)};
"""


def write_fields(case: dict, cdir: str):
    ph = case["physics"]
    ic = ph["initial_conditions"]
    wf = ph["air"]["wall_functions"]
    p0 = num(ic["p"])
    zero_dir = os.path.join(cdir, "0")

    # ------------------------------------------------------------------ air
    d = os.path.join(zero_dir, "air")
    Ta = ic["T_air"]
    _field(os.path.join(d, "T"), "volScalarField", "[0 0 0 1 0 0 0]", f"uniform {num(Ta)}", [
        (ROOM_WALLS, "type zeroGradient;"),               # adiabatic room walls
        ("air_to_shell", _coupled(Ta, "fluidThermo")),
    ])
    _field(os.path.join(d, "U"), "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", [
        (ROOM_WALLS, "type noSlip;"),
        ("air_to_shell", "type noSlip;"),
    ])
    _field(os.path.join(d, "p"), "volScalarField", "[1 -1 -2 0 0 0 0]", f"uniform {p0}", [
        (".*", f"type calculated; value uniform {p0};"),
    ])
    _field(os.path.join(d, "p_rgh"), "volScalarField", "[1 -1 -2 0 0 0 0]", f"uniform {p0}", [
        (".*", f"type fixedFluxPressure; value uniform {p0};"),
    ])
    k0, e0 = num(ic["k_air"]), num(ic["epsilon_air"])
    _field(os.path.join(d, "k"), "volScalarField", "[0 2 -2 0 0 0 0]", f"uniform {k0}", [
        (".*", f"type kqRWallFunction; value uniform {k0};"),
    ])
    _field(os.path.join(d, "epsilon"), "volScalarField", "[0 2 -3 0 0 0 0]", f"uniform {e0}", [
        (".*", f"type epsilonWallFunction; value uniform {e0};"),
    ])
    _field(os.path.join(d, "nut"), "volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", [
        (".*", f"""type nutkWallFunction; Cmu {num(wf['Cmu'])}; kappa {num(wf['kappa'])};
                   E {num(wf['E'])}; value uniform 0;"""),
    ])
    _field(os.path.join(d, "alphat"), "volScalarField", "[1 -1 -1 0 0 0 0]", "uniform 0", [
        (".*", f"""type compressible::alphatJayatillekeWallFunction;
                   Prt {num(ph['air']['turbulence']['Prt'])}; value uniform 0;"""),
    ])

    # ------------------------------------------------------------------ pcm
    d = os.path.join(zero_dir, "pcm")
    Tp = ic["T_pcm"]
    _field(os.path.join(d, "T"), "volScalarField", "[0 0 0 1 0 0 0]", f"uniform {num(Tp)}", [
        ("pcm_to_shell", _coupled(Tp, "fluidThermo")),
    ])
    _field(os.path.join(d, "U"), "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", [
        (".*", "type noSlip;"),                           # Eq. 24: u = 0
    ])
    _field(os.path.join(d, "p"), "volScalarField", "[1 -1 -2 0 0 0 0]", f"uniform {p0}", [
        (".*", f"type calculated; value uniform {p0};"),
    ])
    _field(os.path.join(d, "p_rgh"), "volScalarField", "[1 -1 -2 0 0 0 0]", f"uniform {p0}", [
        (".*", f"type fixedFluxPressure; value uniform {p0};"),
    ])

    # ---------------------------------------------------------------- shell
    d = os.path.join(zero_dir, "shell")
    Ts = ic["T_shell"]
    _field(os.path.join(d, "T"), "volScalarField", "[0 0 0 1 0 0 0]", f"uniform {num(Ts)}", [
        ("ceiling", "type zeroGradient;"),                # adiabatic panel base, assumption (3)
        ("shell_to_air", _coupled(Ts, "solidThermo")),
        ("shell_to_pcm", _coupled(Ts, "solidThermo")),   # Eq. 24 flux continuity
    ])
    _field(os.path.join(d, "p"), "volScalarField", "[1 -1 -2 0 0 0 0]", f"uniform {p0}", [
        (".*", f"type calculated; value uniform {p0};"),
    ])
