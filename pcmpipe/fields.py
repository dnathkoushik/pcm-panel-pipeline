"""0/<region>/ fields: initial (Eqs. 25-27) and boundary conditions (Sec. 3.4).

With a loads scenario (config/loads.yaml) the ceiling and the room walls/floor
follow the tables from pcmpipe/loads.py instead of being adiabatic.
"""
from __future__ import annotations

import os

from . import foam, loads
from .foam import num


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


def _coupled(T, method, layer=None):
    """Coupled interface; `layer` = (thickness, conductivity) of a thin resistive layer."""
    extra = (f"thicknessLayers ({num(layer[0])});\nkappaLayers     ({num(layer[1])});\n"
             if layer else "")
    return f"""
type            compressible::turbulentTemperatureCoupledBaffleMixed;
Tnbr            T;
kappaMethod     {method};
{extra}value           uniform {num(T)};
"""


def _table(time_s, T_C, T0):
    """Prescribed temperature history (Function1 table, held after the last entry)."""
    return f"""
type            uniformFixedValue;
uniformValue    {loads.foam_table(time_s, [v + 273.15 for v in T_C])};
value           uniform {num(T0)};
"""


def _table_through_layer(time_s, T_C, T0, layer):
    """Thin-shell mode, PCM on the ceiling: slab temperature behind the ABS top wall."""
    return f"""
type            externalWallHeatFluxTemperature;
mode            coefficient;
Ta              {loads.foam_table(time_s, [v + 273.15 for v in T_C])};
h               uniform 1e5;       // contact with the slab; the resistance is the layer
thicknessLayers ({num(layer[0])});
kappaLayers     ({num(layer[1])});
kappaMethod     fluidThermo;
value           uniform {num(T0)};
"""


def _room_air(case, T0, profile, layer=None):
    """model.air: coefficient -- the room air as a natural-convection coefficient
    (physics.air_side). Paper set-up: closed adiabatic room, its air gives up the
    heat the panel takes; loads scenario: room air follows T_room(t)."""
    ph = case["physics"]
    hs = ph["air_side"]
    C = hs["C"][case["panel"]["shape"]]
    if profile is None:
        a, room = ph["air"], ph["room"]
        Ta = f"constant {num(ph['initial_conditions']['T_air'])}"
        room_cap = f"roomHeatCapacity {num(a['rho'] * a['cp'] * room['Lx'] * room['Ly'] * room['Lz'])};\n"
    else:
        Ta = loads.foam_table(profile["time_s"], [v + 273.15 for v in profile["T_room_C"]])
        room_cap = ""
    extra = (f"thicknessLayers ({num(layer[0])});\nkappaLayers     ({num(layer[1])});\n"
             if layer else "")
    return f"""
type            naturalConvectionTemperature;
Ta              {Ta};
C               {num(C)};
Cstable         {num(hs['stable_factor'] * C)};
n               {num(hs['n'])};
{extra}{room_cap}value           uniform {num(T0)};
"""


def write_fields(case: dict, cdir: str, profile: dict | None = None):
    """`profile`: loads.profile() of the case's scenario, or None (paper set-up)."""
    ph = case["physics"]
    ic = ph["initial_conditions"]
    wf = ph["air"]["wall_functions"]
    p0 = num(ic["p"])
    zero_dir = os.path.join(cdir, "0")
    thin = case.get("model", {}).get("shell", "resolved") == "thin"
    air = case.get("model", {}).get("air", "resolved") == "resolved"
    # thin: ABS shell as a 1D resistance t/k between air and PCM (no shell region)
    layer = (case["panel"]["shell_thickness_mm"] * 1e-3, ph["shell"]["k"]) if thin else None
    air_if, pcm_if = ("air_to_pcm", "pcm_to_air") if thin else ("air_to_shell", "pcm_to_shell")
    # the air is meshed and its 0/ fields are written either way; with model.air:
    # coefficient it is not solved and the panel surface carries _room_air()

    # ------------------------------------------------------------------ air
    d = os.path.join(zero_dir, "air")
    Ta = ic["T_air"]
    # ".*" comes first: later regexes and exact names override it. It catches the
    # few stray faces snappy can leave on the panel_outer / pcm_cavity patches.
    roof = []
    if profile is not None:      # slab underside and room surfaces follow the loads tables
        t = profile["time_s"]
        roof = [("ceiling", _table(t, profile["T_ceiling_C"], Ta)),
                ("(walls|floor)", _table(t, profile["T_room_C"], Ta))]
    _field(os.path.join(d, "T"), "volScalarField", "[0 0 0 1 0 0 0]", f"uniform {num(Ta)}", [
        (".*", "type zeroGradient;"),                     # adiabatic room walls
        *roof,
        (air_if, _coupled(Ta, "fluidThermo", layer)),
    ])
    _field(os.path.join(d, "U"), "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", [
        (".*", "type noSlip;"),
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
    roof = []
    if profile is not None and thin:     # PCM reaches the ceiling through the ABS top wall
        roof = [("ceiling", _table_through_layer(profile["time_s"], profile["T_ceiling_C"],
                                                 Tp, layer))]
    _field(os.path.join(d, "T"), "volScalarField", "[0 0 0 1 0 0 0]", f"uniform {num(Tp)}", [
        (".*", "type zeroGradient;"),
        *roof,
        (pcm_if, _coupled(Tp, "fluidThermo", layer) if air or not thin
         else _room_air(case, Tp, profile, layer)),
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

    if thin:
        return

    # ---------------------------------------------------------------- shell
    d = os.path.join(zero_dir, "shell")
    Ts = ic["T_shell"]
    roof = []
    if profile is not None:
        roof = [("ceiling", _table(profile["time_s"], profile["T_ceiling_C"], Ts))]
    _field(os.path.join(d, "T"), "volScalarField", "[0 0 0 1 0 0 0]", f"uniform {num(Ts)}", [
        (".*", "type zeroGradient;"),                     # incl. ceiling: adiabatic panel base
        *roof,                                            #   unless a loads scenario is set
        ("shell_to_air", _coupled(Ts, "solidThermo") if air else _room_air(case, Ts, profile)),
        ("shell_to_pcm", _coupled(Ts, "solidThermo")),   # Eq. 24 flux continuity
    ])
    _field(os.path.join(d, "p"), "volScalarField", "[1 -1 -2 0 0 0 0]", f"uniform {p0}", [
        (".*", f"type calculated; value uniform {p0};"),
    ])
