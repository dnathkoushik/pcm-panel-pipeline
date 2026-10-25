"""Step 3 -- chtMultiRegionFoam case: constant/, system/ and 0/ for air|shell|pcm.

Paper -> OpenFOAM mapping
  Eqs. 2-11 (air, Realizable k-eps)   air region, heRhoThermo + Boussinesq EoS (Eq. 17),
                                      RAS realizableKE, standard wall functions (Eq. 23)
  Eq. 12 (ABS shell)                  shell region, heSolidThermo (constant properties)
  Eqs. 13-21 (PCM, laminar)           pcm region, heRhoThermo + Boussinesq EoS (Eq. 16),
                                      laminar, pcmEnthalpyPorositySource fvOption
  Eq. 24 (interfaces)                 turbulentTemperatureCoupledBaffleMixed + noSlip
  SIMPLE / 2nd-order upwind / 0.5 s   PIMPLE (transient CHT) / linearUpwind / maxDeltaT 0.5
"""
from __future__ import annotations

import json
import os

from . import cad, fields, foam
from .foam import num, vec


def _thermo_fluid(f, rho0, T0, beta, cp, mu, Pr, molw):
    return f"""
thermoType
{{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState Boussinesq;
    specie          specie;
    energy          sensibleEnthalpy;
}}

mixture
{{
    specie          {{ molWeight {num(molw)}; }}
    equationOfState {{ rho0 {num(rho0)}; T0 {num(T0)}; beta {num(beta)}; }}
    thermodynamics  {{ Cp {num(cp)}; Hf 0; }}
    transport       {{ mu {num(mu)}; Pr {num(Pr)}; }}
}}
"""


def write_constant(case: dict, cdir: str):
    ph = case["physics"]
    c = os.path.join(cdir, "constant")
    foam.write(os.path.join(c, "regionProperties"), "regions\n(\n    fluid (air pcm)\n    solid (shell)\n);")
    for r in ("air", "pcm"):
        foam.write(os.path.join(c, r, "g"), f"dimensions [0 1 -2 0 0 0 0];\nvalue {vec(ph['gravity'])};",
                   cls="uniformDimensionedVectorField", location=f"constant/{r}")
    for r in ("air", "pcm", "shell"):
        foam.write(os.path.join(c, r, "radiationProperties"), "radiation off;\nradiationModel none;")

    a = ph["air"]
    foam.write(os.path.join(c, "air", "thermophysicalProperties"), _thermo_fluid(
        "air", a["rho"], a["T_ref"], a["beta"], a["cp"], a["mu"], a["mu"] * a["cp"] / a["k"], a["molWeight"]))
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
        "pcm", p["rho"], p["T_ref"], p["beta"], p["cp_liquid"], p["mu"],
        p["mu"] * p["cp_liquid"] / p["k_liquid"], p["molWeight"]))
    foam.write(os.path.join(c, "pcm", "turbulenceProperties"), "simulationType  laminar;")
    foam.write(os.path.join(c, "pcm", "fvOptions"), f"""
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
}}
""")

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
    p_rghFinal      {{ $p_rgh; relTol 0; }}
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

relaxationFactors {{ equations {{ ".*" 1; }} }}
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

    foam.write(os.path.join(sysd, "controlDict"), f"""
application     chtMultiRegionFoam;
libs            ("libpcmEnthalpyPorosity.so");

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {num(end_s)};
deltaT          {num(sv['deltaT'])};

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
        patches         (pcm_to_shell);
        writeControl    adjustableRunTime;
        writeInterval   {dt_sample};
        writeFields     false;
    }}

    airResiduals
    {{
        type            solverInfo;
        libs            (utilityFunctionObjects);
        region          air;
        fields          (U p_rgh h k epsilon);
        writeControl    timeStep;
        writeInterval   20;
    }}

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
    write_decomposition(cdir, run["n_procs"], run["decomposition"])

    div = "linearUpwind" if sv["div_scheme"] == "linearUpwind" else sv["div_scheme"]
    for r, turb, pcm in (("air", True, False), ("pcm", False, True)):
        foam.write(os.path.join(sysd, r, "fvSchemes"), _fluid_schemes(div, turb, pcm))
        foam.write(os.path.join(sysd, r, "fvSolution"), _fluid_solution(sv, turb))
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


def write_decomposition(cdir: str, n_procs: int, method: str = "scotch"):
    """decomposeParDict for the case and every region (must match the MPI rank count)."""
    body = f"numberOfSubdomains {n_procs};\nmethod {method};"
    for sub in ("", "air", "pcm", "shell"):
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
        "delta_T": case["delta_T"], "mesh_preset": case["mesh"]["preset"],
        "T_pcm0": case["physics"]["initial_conditions"]["T_pcm"],
        "T_air": case["physics"]["initial_conditions"]["T_air"],
        "probes": geo.probes, "slice_plane": geo.slice_plane,
        "panel_bounds_m": {
            "x": [min(p[0] for p in geo.footprint) * 1e-3, max(p[0] for p in geo.footprint) * 1e-3],
            "y": [min(p[1] for p in geo.footprint) * 1e-3, max(p[1] for p in geo.footprint) * 1e-3],
            "z": [geo.z_tip * 1e-3, 0.0]},
        "exposed_area_cm2_analytic": geo.exposed_area_mm2 / 100.0,
        "latent_heat": case["physics"]["pcm"]["latent_heat"],
    }
    try:
        rep = cad.cad_report(case)
        meta.update(pcm_mass_g=rep["pcm_mass_g"], cavity_volume_ml=rep["cavity_volume_ml"],
                    exposed_area_cm2=rep["exposed_area_cm2"])
    except FileNotFoundError:   # dry run without CAD
        meta["pcm_mass_g"] = None
    return meta


def build(case: dict):
    cdir = os.path.join(case["dir"], "case")
    meta = case_meta(case)
    write_constant(case, cdir)
    write_system(case, cdir, meta)
    fields.write_fields(case, cdir)
    with open(os.path.join(case["dir"], "case_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    open(os.path.join(cdir, "case.foam"), "w").close()   # ParaView reader handle
    print(f"  case written: {cdir}")
