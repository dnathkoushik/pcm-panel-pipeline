"""Step 3 -- chtMultiRegionFoam case: constant/, system/ and 0/ for air|shell|pcm.

Paper -> OpenFOAM mapping
  Eqs. 2-11 (air, Realizable k-eps)   air region, heRhoThermo + rhoConst, Boussinesq
                                      buoyancy fvOption (Eq. 17), RAS realizableKE,
                                      standard wall functions (Eq. 23)
  Eq. 12 (ABS shell)                  shell region, heSolidThermo (constant properties)
  Eqs. 13-21 (PCM, laminar)           pcm region, heRhoThermo + rhoConst, Boussinesq
                                      buoyancy (Eq. 16) and pcmEnthalpyPorositySource fvOptions
  Eq. 24 (interfaces)                 turbulentTemperatureCoupledBaffleMixed + noSlip
  SIMPLE / 2nd-order upwind / 0.5 s   PIMPLE (transient CHT) / linearUpwind / maxDeltaT 0.5
"""
from __future__ import annotations

import json
import os

from . import cad, config, fields, foam, loads, mesh
from .foam import num, vec


def _thermo_fluid(f, rho0, cp, mu, Pr, molw):
    return f"""
thermoType
{{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}}

// Eqs. 4/15 have no dp/dt term. With an incompressible EoS the
// pressure level floats in the closed room, and dp/dt made T unbounded (CI run).
dpdt            no;

// Boussinesq approximation as in the paper (Fluent): constant density in every
// equation, buoyancy from the boussinesqBuoyancySource fvOption. The Boussinesq
// EoS would put rho(T) into continuity and energy as well; in the sealed PCM
// capsule and the closed room that broke mass conservation and drove single
// cells tens of kelvin out of bounds (CI diagnostic: gone with beta = 0).
mixture
{{
    specie          {{ molWeight {num(molw)}; }}
    equationOfState {{ rho {num(rho0)}; }}
    thermodynamics  {{ Cp {num(cp)}; Hf 0; }}
    transport       {{ mu {num(mu)}; Pr {num(Pr)}; }}
}}
"""


def _buoyancy(beta, T0, g):
    """Boussinesq buoyancy rho0 beta (T0 - T) g (Eqs. 16/17) for a fluid region's fvOptions."""
    return f"""
boussinesq
{{
    type            boussinesqBuoyancySource;
    active          yes;
    selectionMode   all;
    beta            {num(beta)};
    TRef            {num(T0)};
    g               {vec(g)};
}}
"""


def write_constant(case: dict, cdir: str):
    ph = case["physics"]
    c = os.path.join(cdir, "constant")
    thin = mesh.thin_shell(case)
    solids = "" if thin else "shell"     # thin: shell is a resistance in the BC, no region
    fluids = "air pcm" if mesh.air_resolved(case) else "pcm"   # coefficient: air not solved
    foam.write(os.path.join(c, "regionProperties"),
               f"regions\n(\n    fluid ({fluids})\n    solid ({solids})\n);")
    # openfoam.com chtMultiRegionFoam reads one case-level g (as in its tutorials)
    foam.write(os.path.join(c, "g"), f"dimensions [0 1 -2 0 0 0 0];\nvalue {vec(ph['gravity'])};",
               cls="uniformDimensionedVectorField", location="constant")
    for r in mesh.regions(case):
        foam.write(os.path.join(c, r, "radiationProperties"), "radiation off;\nradiationModel none;")

    a = ph["air"]
    foam.write(os.path.join(c, "air", "thermophysicalProperties"), _thermo_fluid(
        "air", a["rho"], a["cp"], a["mu"], a["mu"] * a["cp"] / a["k"], a["molWeight"]))
    foam.write(os.path.join(c, "air", "fvOptions"), _buoyancy(a["beta"], a["T_ref"], ph["gravity"]))
    tu = a["turbulence"]
    foam.write(os.path.join(c, "air", "turbulenceProperties"), f"""
simulationType  RAS;

RAS
{{
    RASModel        {tu['model']};
    turbulence      on;
    printCoeffs     on;
    Prt             {num(tu['Prt'])};
    {tu['model']}Coeffs
    {{
        A0          {num(tu['A0'])};
        C2          {num(tu['C2'])};
        sigmak      {num(tu['sigmak'])};
        sigmaEps    {num(tu['sigmaEps'])};
    }}
}}
""")

    # PCM: the thermo package carries the liquid cp and k as reference values;
    # the fvOption restores the solid/mushy values (see pcmEnthalpyPorositySource.H)
    p = ph["pcm"]
    foam.write(os.path.join(c, "pcm", "thermophysicalProperties"), _thermo_fluid(
        "pcm", p["rho"], p["cp_liquid"], p["mu"],
        p["mu"] * p["cp_liquid"] / p["k_liquid"], p["molWeight"]))
    foam.write(os.path.join(c, "pcm", "turbulenceProperties"), "simulationType  laminar;")
    foam.write(os.path.join(c, "pcm", "fvOptions"), _buoyancy(p["beta"], p["T_ref"], ph["gravity"]) + f"""
pcmPhaseChange
{{
    type            pcmEnthalpyPorositySource;
    active          yes;
    selectionMode   all;
    energyField     h;

    Tsolidus        {num(p['T_solidus'])};    // Eq. 1
    Tliquidus       {num(p['T_liquidus'])};
    L               {num(p['latent_heat'])};   // Eq. 20
    Amush           {num(p['A_mush'])};        // Eq. 21
    epsilon         {num(p['eps_mush'])};
    CpSolid         {num(p['cp_solid'])};      // Table 2, linear in the mushy zone
    CpLiquid        {num(p['cp_liquid'])};
    kappaSolid      {num(p['k_solid'])};
    kappaLiquid     {num(p['k_liquid'])};
    latentConvection true;                     // div(rho u gamma L) in Eq. 15
    kinkLimiter     true;                      // outer iterations stop at Tsol/Tliq
}}
""")

    if thin:
        return
    s = ph["shell"]
    foam.write(os.path.join(c, "shell", "thermophysicalProperties"), f"""
thermoType
{{
    type            heSolidThermo;
    mixture         pureMixture;
    transport       constIso;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}}

mixture
{{
    specie          {{ molWeight {num(s['molWeight'])}; }}
    transport       {{ kappa {num(s['k'])}; }}
    thermodynamics  {{ Hf 0; Cp {num(s['cp'])}; }}
    equationOfState {{ rho {num(s['rho'])}; }}
}}
""")


def _fluid_schemes(div_scheme, turbulent, pcm):
    lu = div_scheme
    extra = ""
    if turbulent:
        extra += f"    div(phi,k)      Gauss {lu} grad(k);\n    div(phi,epsilon) Gauss {lu} grad(epsilon);\n"
    if pcm:
        extra += f"    div(phi,pcmExtraEnthalpy) Gauss {lu} grad(pcmExtraEnthalpy);\n"
    return f"""
ddtSchemes      {{ default Euler; }}

gradSchemes
{{
    default         Gauss linear;
    limited         cellLimited Gauss linear 1;
    grad(U)         $limited;
    grad(h)         $limited;
    grad(k)         $limited;
    grad(epsilon)   $limited;
    grad(pcmExtraEnthalpy) $limited;
}}

divSchemes
{{
    default         none;
    div(phi,U)      Gauss {lu} grad(U);
    div(phi,h)      Gauss {lu} grad(h);
    div(phi,K)      Gauss linear;
{extra}    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}}

laplacianSchemes {{ default Gauss linear limited corrected 0.5; }}
interpolationSchemes {{ default linear; }}
snGradSchemes   {{ default limited corrected 0.5; }}
wallDist        {{ method meshWave; }}
"""


def _fluid_solution(sv, turbulent):
    tol_e = num(sv["residuals"]["energy"])
    tol_o = num(sv["residuals"]["others"])
    turb = "|k|epsilon" if turbulent else ""
    return f"""
solvers
{{
    "rho.*"         {{ solver PCG; preconditioner DIC; tolerance 0; relTol 0; }}
    p_rgh           {{ solver GAMG; smoother GaussSeidel; tolerance 1e-6; relTol 0.01; }}
    // maxIter: in the room air GAMG can stall just above 1e-6 and then ran its
    // default 1000 sweeps every step (4 s/step, long screening runs 27 Oct)
    p_rghFinal      {{ $p_rgh; relTol 0; maxIter 50; }}
    "(U{turb})"     {{ solver PBiCGStab; preconditioner DILU; tolerance {tol_o}; relTol 0.1; }}
    "(U{turb})Final" {{ solver PBiCGStab; preconditioner DILU; tolerance {tol_o}; relTol 0; }}
    h               {{ solver PBiCGStab; preconditioner DILU; tolerance {tol_e}; relTol 0.1; }}
    hFinal          {{ $h; relTol 0; }}
}}

PIMPLE
{{
    momentumPredictor        on;
    nCorrectors              {sv['nCorrectors']};
    nNonOrthogonalCorrectors {sv['nNonOrthogonalCorrectors']};
    pRefCell                 0;      // closed room / sealed capsule
    pRefValue                1e5;
}}

// Momentum/pressure relaxed on the non-final outer corrector only; h stays at 1
// so the latent-heat Newton iteration converges as in tests/test_latent_source_1d.py
relaxationFactors
{{
    fields      {{ p_rgh 0.7; p_rghFinal 1; }}
    equations   {{ "(U|k|epsilon)" 0.7; "(U|k|epsilon)Final" 1; "h.*" 1; }}
}}
"""


def write_system(case: dict, cdir: str, meta: dict):
    ph, run = case["physics"], case["run"]
    sv = ph["solver"]
    sysd = os.path.join(cdir, "system")
    if run["end_time_min"] == "auto":
        end_s = 60.0 * run["max_end_time_min"]
    else:
        end_s = 60.0 * float(run["end_time_min"])
    probes = " ".join(vec(v) for v in meta["probes"].values())
    dt_sample = num(run["sample_interval_s"])
    thin = mesh.thin_shell(case)
    pcm_if = mesh.fluid_interfaces(case)["pcm"]
    shell_minmax = "" if thin else """
    shellTMinMax
    {
        type            fieldMinMax;
        libs            (fieldFunctionObjects);
        region          shell;
        fields          (T);
        mode            magnitude;
        location        true;
        writeControl    timeStep;
        writeInterval   5;
    }
"""
    air_functions = "" if not mesh.air_resolved(case) else """
    airTMinMax
    {
        type            fieldMinMax;
        libs            (fieldFunctionObjects);
        region          air;
        fields          (T U);
        mode            magnitude;
        location        true;
        writeControl    timeStep;
        writeInterval   5;
    }

    airResiduals
    {
        type            solverInfo;
        libs            (utilityFunctionObjects);
        region          air;
        fields          (U p_rgh h k epsilon);
        writeControl    timeStep;
        writeInterval   20;
    }
"""

    foam.write(os.path.join(sysd, "controlDict"), f"""
application     chtMultiRegionFoam;
libs            ("libpcmEnthalpyPorosity.so");

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {num(end_s)};
deltaT          {num(sv.get('deltaT_initial', sv['deltaT']) if sv['adjustTimeStep'] else sv['deltaT'])};

writeControl    adjustableRunTime;
writeInterval   {num(60 * run['write_interval_min'])};
purgeWrite      0;
writeFormat     binary;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   10;
runTimeModifiable true;

adjustTimeStep  {'yes' if sv['adjustTimeStep'] else 'no'};
maxCo           {num(sv['maxCo'])};
maxDi           {num(sv['maxDi'])};
maxDeltaT       {num(sv['maxDeltaT'])};

functions
{{
    pcmAverages        // liquid fraction (Figs. 10, 11) and mean PCM temperature (Figs. 13-17)
    {{
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        region          pcm;
        writeControl    adjustableRunTime;
        writeInterval   {dt_sample};
        log             false;
        writeFields     false;
        regionType      all;
        operation       volAverage;
        fields          (liquidFraction T);
    }}

    pcmProbes          // thermocouples top / mid / bot (Sec. 2.4, Fig. 8)
    {{
        type            probes;
        libs            (sampling);
        region          pcm;
        writeControl    adjustableRunTime;
        writeInterval   {dt_sample};
        fields          (T liquidFraction);
        probeLocations  ({probes});
    }}

    pcmWallHeatFlux    // heat flow into the PCM, for the energy balance check
    {{
        type            wallHeatFlux;
        libs            (fieldFunctionObjects);
        region          pcm;
        patches         ({pcm_if});
        writeControl    adjustableRunTime;
        writeInterval   {dt_sample};
        writeFields     false;
    }}

    pcmTMinMax         // where are the extreme temperatures? (boundedness check)
    {{
        type            fieldMinMax;
        libs            (fieldFunctionObjects);
        region          pcm;
        fields          (T U);
        mode            magnitude;
        location        true;
        writeControl    timeStep;
        writeInterval   5;
    }}
{shell_minmax}{air_functions}
    pcmResiduals
    {{
        type            solverInfo;
        libs            (utilityFunctionObjects);
        region          pcm;
        fields          (U p_rgh h);
        writeControl    timeStep;
        writeInterval   20;
    }}
{_stop_block(run) if run['end_time_min'] == 'auto' else ''}}}
""")
    foam.write(os.path.join(sysd, "fvSchemes"), "")
    foam.write(os.path.join(sysd, "fvSolution"), f"PIMPLE\n{{\n    nOuterCorrectors {sv['nOuterCorrectors']};\n}}")
    write_decomposition(cdir, run["n_procs"], run["decomposition"], mesh.regions(case))

    div = "linearUpwind" if sv["div_scheme"] == "linearUpwind" else sv["div_scheme"]
    for r, turb, pcm in (("air", True, False), ("pcm", False, True)):
        foam.write(os.path.join(sysd, r, "fvSchemes"), _fluid_schemes(div, turb, pcm))
        foam.write(os.path.join(sysd, r, "fvSolution"), _fluid_solution(sv, turb))
    if thin:
        return
    foam.write(os.path.join(sysd, "shell", "fvSchemes"), """
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes      { default none; }
laplacianSchemes { default Gauss linear limited corrected 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited corrected 0.5; }
""")
    foam.write(os.path.join(sysd, "shell", "fvSolution"), f"""
solvers
{{
    h               {{ solver PCG; preconditioner DIC; tolerance {num(sv['residuals']['energy'])}; relTol 0.1; }}
    hFinal          {{ $h; relTol 0; }}
}}

PIMPLE
{{
    nNonOrthogonalCorrectors {sv['nNonOrthogonalCorrectors']};
}}
""")


def write_decomposition(cdir: str, n_procs: int, method: str = "scotch",
                        region_names=("air", "shell", "pcm")):
    """decomposeParDict for the case and every region (must match the MPI rank count)."""
    body = f"numberOfSubdomains {n_procs};\nmethod {method};"
    for sub in ("", *region_names):
        foam.write(os.path.join(cdir, "system", sub, "decomposeParDict"), body)


def _stop_block(run):
    return f"""
    meltComplete       // stop once the PCM is melted (run.end_time_min: auto)
    {{
        type            runTimeControl;
        libs            (utilityFunctionObjects);
        conditions
        {{
            melted
            {{
                type            minMax;
                functionObject  pcmAverages;
                fields          (volAverage(liquidFraction));
                mode            maximum;
                value           {num(run['stop_liquid_fraction'])};
            }}
        }}
    }}
"""


def case_meta(case: dict) -> dict:
    geo = cad.panel_geometry(case)
    meta = {
        "id": case["id"], "study": case.get("study"), "shape": geo.shape,
        "volume_ml": geo.volume_ml, "protrusion_ratio": geo.protrusion_ratio,
        "shell_thickness_mm": geo.shell_thickness, "gap_mm": geo.gap,
        "delta_T": None if config.loads_active(case) else case["delta_T"],
        "mesh_preset": case["mesh"]["preset"],
        "pcm_material": case.get("pcm_material") or "paper_cla",
        "loads": case.get("loads") if config.loads_active(case) else None,
        "start_clock_h": (case["loads_spec"]["start_clock_h"] if config.loads_active(case)
                          else None),
        "T_melt": case["physics"]["pcm"]["T_pc"],
        "T_solidus": case["physics"]["pcm"]["T_solidus"],
        "T_liquidus": case["physics"]["pcm"]["T_liquidus"],
        "T_pcm0": case["physics"]["initial_conditions"]["T_pcm"],
        "T_air": case["physics"]["initial_conditions"]["T_air"],
        "probes": geo.probes, "slice_plane": geo.slice_plane,
        "panel_bounds_m": {
            "x": [min(p[0] for p in geo.footprint) * 1e-3, max(p[0] for p in geo.footprint) * 1e-3],
            "y": [min(p[1] for p in geo.footprint) * 1e-3, max(p[1] for p in geo.footprint) * 1e-3],
            "z": [geo.z_tip * 1e-3, 0.0]},
        "exposed_area_cm2_analytic": geo.exposed_area_mm2 / 100.0,
        "latent_heat": case["physics"]["pcm"]["latent_heat"],
        "mode": case.get("mode", "paper"),
        "shell_model": "thin" if mesh.thin_shell(case) else "resolved",
        "air_model": "resolved" if mesh.air_resolved(case) else "coefficient",
    }
    try:
        rep = cad.cad_report(case)
        if mesh.thin_shell(case):   # PCM fills the cavity up to the ceiling
            meta.update(pcm_mass_g=rep["pcm_mass_thin_g"],
                        exposed_area_cm2=rep["cavity_thin_area_cm2"])
        else:
            meta.update(pcm_mass_g=rep["pcm_mass_g"], exposed_area_cm2=rep["exposed_area_cm2"])
        meta["cavity_volume_ml"] = rep["cavity_volume_ml"]
    except (FileNotFoundError, KeyError):   # dry run without CAD / old CAD report
        meta["pcm_mass_g"] = None
    return meta


def build(case: dict):
    cdir = os.path.join(case["dir"], "case")
    meta = case_meta(case)
    profile = (loads.write_profile(case, case["dir"], config.ROOT)
               if config.loads_active(case) else None)
    write_constant(case, cdir)
    write_system(case, cdir, meta)
    fields.write_fields(case, cdir, profile)
    with open(os.path.join(case["dir"], "case_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    open(os.path.join(cdir, "case.foam"), "w").close()   # ParaView reader handle
    print(f"  case written: {cdir}")
