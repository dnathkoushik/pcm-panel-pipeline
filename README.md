# PCM ceiling panel CFD pipeline

[![CI](https://github.com/dnathkoushik/pcm-panel-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/dnathkoushik/pcm-panel-pipeline/actions/workflows/ci.yml)

This is an automated, open-source pipeline to design, mesh, simulate and evaluate
phase change material (PCM) ceiling panels. All physics follows
**Yang et al., "PCM-based ceiling panels for passive cooling in buildings: A CFD
modelling", Energy & Buildings 285 (2023) 112898** (Sections 2.2, 2.3 and 3).

```
config/*.yaml ──► FreeCAD (CAD, STL, STEP) ──► snappyHexMesh (3 regions) ──► chtMultiRegionFoam
                                                                        + pcmEnthalpyPorositySource
      runs/_figures/<study>/*.png ◄── matplotlib ◄── metrics.json ◄── function objects / pvpython
```

| Stage | Tool | Code |
|---|---|---|
| 1. Geometry | FreeCAD (Part/OCC), headless `freecadcmd` | `pcmpipe/geometry.py`, `cad/build_panel_freecad.py` |
| 2. Meshing | blockMesh, snappyHexMesh, splitMeshRegions, refineWallLayer | `pcmpipe/mesh.py` |
| 3. Solver | chtMultiRegionFoam + custom fvOption (C++) | `pcmpipe/case_setup.py`, `pcmpipe/fields.py`, `openfoam/src/` |
| 4. Execution / post | MPI, function objects, `pvpython` | `pcmpipe/run.py`, `pcmpipe/results.py`, `post/pv_*.py` |
| 5. Figures | matplotlib | `post/plot_results.py` |

## Directory structure

```
pcm_panel_pipeline/
├── README.md
├── requirements.txt
├── config/
│   ├── physics.yaml          # paper ground truth; every value tagged [PAPER]/[ASSUMPTION]/[INTERPRET]
│   ├── base_case.yaml        # panel defaults, mesh presets, run and post settings, tool paths
│   ├── materials.yaml        # PCMs: paper_cla (the paper's), n_octadecane (Faden et al. 2019)
│   ├── loads.yaml            # thermal loads: none (paper set-up) | kharagpur_roof
│   └── studies.yaml          # parametric studies (cartesian products of overrides)
├── pcmpipe/                  # Python driver:  python3 -m pcmpipe <step> --study <name>
│   ├── __main__.py           # CLI
│   ├── config.py             # YAML loading, study expansion, case ids
│   ├── geometry.py           # parametric dimension solver (pure Python)
│   ├── cad.py                # calls FreeCAD
│   ├── mesh.py               # blockMesh/snappyHexMesh dictionaries + meshing sequence
│   ├── case_setup.py         # constant/ and system/ for air | shell | pcm
│   ├── fields.py             # 0/ initial and boundary conditions
│   ├── loads.py              # sol-air temperature + 1D roof slab -> ceiling/room T(t) tables
│   ├── run.py                # decomposePar, mpirun chtMultiRegionFoam, reconstructPar
│   ├── results.py            # melting times, rates, normalised rates
│   └── foam.py               # dictionary writer, command runner
├── cad/build_panel_freecad.py
├── openfoam/
│   ├── Allwmake
│   └── src/pcmEnthalpyPorositySource/   # enthalpy-porosity fvOption (Eqs. 1, 15, 18-21)
├── post/
│   ├── pv_contours.py        # pvpython: T and |U| contours at 50/150/250/300 min (Fig. 9)
│   ├── pv_liquid_fraction_3d.py # pvpython: 3D cut-away liquid fraction + melting front
│   ├── pv_extract_transient.py  # pvpython: liquid fraction / mean T for every written time
│   └── plot_results.py       # paper-style figures (Figs. 8, 10-13, 15, 17, 18)
├── validation/paper_reference.json   # numbers quoted in the paper, overlaid on plots
├── tests/
│   ├── test_geometry.py      # geometry solver vs. the paper's Table 3, Sec. 4.5, Sec. 4.7
│   ├── test_latent_source_1d.py  # latent-heat formulation (incl. n-octadecane) vs. reference
│   ├── test_loads.py         # roof loads model, materials, --set overrides
│   └── test_pipeline_offline.py  # every study: dictionaries, results, plots (no OpenFOAM)
├── scripts/
│   ├── docker_run.sh         # run anything inside the container (Linux/macOS/Windows)
│   ├── install_ubuntu.sh     # native Ubuntu/Debian/WSL2
│   ├── install_macos.sh      # native macOS (Homebrew, OpenFOAM.app)
│   └── run_all_studies.sh
├── Dockerfile, docker/entrypoint.sh   # OpenFOAM + FreeCAD + ParaView image
├── Makefile                  # make test | smoke | study | docker-*
├── .github/workflows/ci.yml  # tests on Linux/macOS/Windows + Docker smoke run
└── runs/                     # created at run time: one folder per case + _figures/
```

Each case folder `runs/<case_id>/` contains:
- `cad/`: STL, STEP and `cad_report.json`
- `case/`: the OpenFOAM case
- `post/`: contour PNGs and CSVs
- `case.json`, `case_meta.json`, `metrics.json`, `series.csv`

Case ids encode the parameters, for example `pyr_V250_r50_t1_g0_dT20`: pyramid, 250 mL,
50 % protrusion ratio, 1 mm shell, no gap, ΔT = 20 K. Identical cases shared by
several studies are only run once.

## 1. Installation

Clone the repository to a path **without spaces**; OpenFOAM and `wmake` cannot handle them.

```bash
git clone https://github.com/dnathkoushik/pcm-panel-pipeline.git
cd pcm-panel-pipeline
```

Choose one of the three setups:

| | Linux | macOS | Windows |
|---|---|---|---|
| **A. Docker** (recommended, same environment as CI) | ✓ | ✓ Intel and Apple silicon | ✓ via WSL2 |
| **B. Native** | Ubuntu 22.04/24.04, Debian 12 | Homebrew + OpenFOAM.app | WSL2 Ubuntu |

You need Python 3.10+ with `numpy`, `matplotlib` and `pyyaml`. OpenFOAM must be
**v2306 or newer from openfoam.com** (ESI). The Foundation releases (openfoam.org) use a
different fvOption API and will not compile the library.

### A. Docker (any OS)

Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) (macOS/Windows)
or Docker Engine (Linux), then:

```bash
make docker-build      # ~10-15 min once: OpenFOAM v2406, FreeCAD, ParaView, Python
make docker-test       # Python tests inside the container
make docker-smoke      # real run: CAD -> mesh -> chtMultiRegionFoam -> pvpython -> plots
scripts/docker_run.sh  # interactive shell; the repo is mounted at /work
```

Any pipeline command can be run as `scripts/docker_run.sh python3 -m pcmpipe ...`.
Results land in `runs/` on your machine.
- **Apple silicon:** the image is amd64 and runs under emulation, so it is several times
  slower than native. Give Docker at least 8 GB of memory (Settings → Resources).
- **Windows:** run the commands from a WSL2 terminal.
- **Without make:** use `docker build --platform linux/amd64 -t pcm-panel-pipeline .` and
  `scripts/docker_run.sh make smoke`.

### B1. Native Linux (Ubuntu / Debian / WSL2)

```bash
bash scripts/install_ubuntu.sh   # OpenFOAM v2406 (apt repo), FreeCAD, ParaView, Python, compiles the fvOption
source ~/.bashrc
make test && make smoke
```

On a headless server, ParaView needs a virtual display:
`export PCMPIPE_PVPYTHON="xvfb-run -a pvpython"`.

On Windows, first run `wsl --install -d Ubuntu-24.04` in an administrator PowerShell,
then clone *inside* the Linux home folder. `/mnt/c` is very slow for OpenFOAM.

### B2. Native macOS

```bash
bash scripts/install_macos.sh    # brew: gerlero/openfoam/openfoam, FreeCAD, ParaView; Python venv
openfoam                         # shell with the OpenFOAM environment
source .venv/bin/activate
make lib test smoke
```

[OpenFOAM.app](https://github.com/gerlero/openfoam-app) provides native openfoam.com
builds: the latest version on Apple silicon, up to v2506 on Intel Macs. It compiles with
Apple Clang. If `make lib` fails on macOS, use Docker. This route is not covered by CI.

### Tool locations

FreeCAD, ParaView and MPI are looked up in this order:
1. the environment variables `PCMPIPE_FREECADCMD`, `PCMPIPE_PVPYTHON` and `PCMPIPE_MPIRUN`
2. `tools:` in `config/base_case.yaml`
3. `PATH`
4. `/Applications/*.app` on macOS

Values may include arguments, e.g. `PCMPIPE_MPIRUN="mpirun --oversubscribe"`.
`PCMPIPE_RUNS_DIR` moves the output folder.

## 2. Verify the installation

```bash
make test        # = python3 -m unittest discover -s tests -v   (~40 s, no OpenFOAM needed)
make smoke NP=4  # CAD, coarse smoke mesh, 2 simulated minutes, post, plots
```

`make test` covers:
- the geometry against the paper's numbers
- the latent-heat discretisation against a converged enthalpy solution
- writing all 39 cases, then running `results` and every plot on synthetic solver output

The smoke run checks the OpenFOAM plumbing; it is not physics. Inspect it before launching
long runs:
- `runs/pyr_V250_r50_t1_g0_dT20_smoke/case/log.*`: every application should end with `End`.
- `log.checkMesh.{air,shell,pcm}` should show three regions.
- The `chtMultiRegionFoam` log should print `pcmEnthalpyPorositySource: ... initial liquid fraction 0`.

Every push runs the same checks on GitHub Actions (`.github/workflows/ci.yml`). The logs,
STLs, contour images and figures of the CI smoke run can be downloaded as the
`smoke-results` artifact.

## 3. Running the pipeline step by step

Every step takes `--study <name>` (from `config/studies.yaml`), plus `--case <text>` to
select one case and `--np <N>` for the number of MPI ranks.

```bash
python3 -m pcmpipe cad     --study validation              # FreeCAD: STL/STEP + volumes/areas
python3 -m pcmpipe mesh    --study validation --np-mesh 8  # snappyHexMesh, 3 regions, wall layers
python3 -m pcmpipe setup   --study validation              # 0/, constant/, system/
python3 -m pcmpipe run     --study validation --np 8       # chtMultiRegionFoam (stops when melted)
python3 -m pcmpipe results --study validation              # t80, t100, rates -> metrics.json
python3 -m pcmpipe post    --study validation              # pvpython contours + transient CSV
python3 -m pcmpipe plots   --study validation              # runs/_figures/validation/
# or all of the above:
python3 -m pcmpipe all --study validation --np 8 --np-mesh 8
```

Run `setup` after `mesh`, because `mesh` recreates the case folder. Use `--dry-run` to
write the dictionaries and print the commands without running OpenFOAM.

### Two modes: screening and paper

The 1 mm ABS shell needs ~0.5 mm cells to resolve. Those cells cap the time step at
about 0.015 s (Courant limit), against the paper's 0.5 s. In the CI smoke run that meant
about 14 s per step on 2 cores, so a full 5-hour melt in paper mode takes weeks on a
laptop. Every study can therefore run in one of two modes (`--mode`, or `MODE=` with make):

| | `--mode screening` | `--mode paper` (default) |
|---|---|---|
| Purpose | rank many designs | confirm the best 2–3 designs, validate against the paper |
| ABS shell | thin thermal resistance t/k in the air–PCM coupling (`thicknessLayers`, `kappaLayers`); no shell cells | solid region, Eq. 12 |
| PCM | fills the cavity up to the (adiabatic) ceiling: cavity + footprint × t | the cavity |
| Mesh preset | `screen`: 1 mm surface cells | `medium` |
| Results folder | `runs/<case>_screen/`, figures in `runs/_figures/<study>_screening/` | `runs/<case>/`, `runs/_figures/<study>/` |

The only physics screening mode drops is the heat stored in the shell itself. For the
250 mL pyramid that is about 70 J/K (35 g of ABS), or 1.0–1.4 kJ over a 15–20 K rise:
5–7 % of the PCM's 21 kJ latent heat. The bias is similar for every design, so screening can
rank designs but not give absolute melting times. Its geometry approximation (PCM up to the
ceiling, about 5 % more PCM for 1 mm) is reported in `metrics.json` through the actual PCM
mass. Rank designs on screening results, then re-run the shortlist in paper mode:

```bash
make study STUDY=shapes MODE=screening NP=8      # or: python3 -m pcmpipe all --study shapes --mode screening
make study STUDY=volume MODE=screening NP=8
make study STUDY=validation NP=8                 # paper mode: compare with the paper's 292 / 353 min
```

### Presentation figures (no CFD needed)

`python3 post/presentation_figures.py` writes to `runs/_figures/presentation/`:
- `design_catalogue.png`: 3D views of the candidate panels
- `geometry_validation.png` and `.md`: our geometry against the paper's Table 3 and Sec. 4.5/4.7 numbers
- `latent_heat_verification.png`: the OpenFOAM phase-change source term against a converged enthalpy solution

### Studies

| Study | Varies | Paper | Main figure |
|---|---|---|---|
| `validation` | pyramid, tetrahedron | Sec. 3.7, 4.1 | probe temperatures (Fig. 8), liquid fraction + melting rate (Fig. 10) |
| `shapes` | + hemisphere | new | liquid fraction curves |
| `delta_T` | ΔT = 15, 20, 25, 30 K | Sec. 4.2 | melting time vs. ΔT (Fig. 11) |
| `shell_thickness` | 1, 3, 5 mm | Sec. 4.3 | Fig. 12 |
| `protrusion_ratio` | 30, 50, 70 % | Sec. 4.5 | Fig. 15 |
| `gap` | 0, 3, 5 mm | Sec. 4.6 | Fig. 17 |
| `volume` | 250, 500, 750 mL | Sec. 4.7 | melting time + normalised melting rate (Fig. 18) |
| `mesh_independence` | coarse / medium / fine | Sec. 3.5 | Fig. 7 |
| `smoke` | smoke mesh, 2 min | | plumbing test |
| `octadecane_roof` | pyramid, tetrahedron, hemisphere | new | n-octadecane under a Kharagpur roof (section 3b) |
| `smoke_roof` | smoke mesh | | plumbing test of the roof set-up |

To add a design, add a study to `studies.yaml`. Any key in `base_case.yaml` can be
varied with a dotted path, for example `panel.n_x`, `panel.hemisphere_fill` or `mesh.preset`.

### Outputs

- `runs/_figures/<study>/`:
  - `melting_time_bars.png`: 80 % solid bars, 100 % hatched bars, as in the paper
  - `liquid_fraction_curves.png`
  - `mean_pcm_temperature.png`
  - `normalised_rate_{physical,paper}.png` (volume study)
  - `liquid_fraction_and_rate.png` and `probe_temperatures.png` (validation and shapes studies)
  - `summary.csv`
  - `liquid_fraction_vs_time.csv`, `mean_pcm_temperature_vs_time.csv`: every case's curve
    on one time axis (one column per case, temperature in °C)
  - `liquid_fraction_3d_montage.png`: the 3D views below, one row per case
- `runs/<case>/post/liquid_fraction_3d_<t>min.png`: PCM cut on the symmetry plane,
  coloured by liquid fraction (blue solid, red liquid), melting front (γ = 0.5) in grey.
- `runs/<case>/post/contour_{T,U}_<t>min.png`: symmetry-plane contours (paper Fig. 9),
  on the paper's colour ranges (14–25 °C, 0–22.8 mm/s), plus `slice_*.csv`.
- `runs/<case>/case/postProcessing/`: liquid fraction and mean PCM temperature every 30 s,
  thermocouple probes (top/mid/bottom at 1, 2/3, 1/3 of the protrusion height), wall heat
  flow into the PCM, and solver residuals.

**Normalised melting rate.** The paper's Fig. 18 values equal 170 g / t100 / V for
*every* volume: the 250 mL PCM mass was used for the 500 and 750 mL panels too.
For example, 170 g / 430 min / 0.5 L = 0.79, which matches the published 0.786.
`normalised_rate_physical.png` uses the PCM mass actually in each panel.
`normalised_rate_paper.png` reproduces the paper's convention for comparison.

## 3b. n-octadecane panels under a Kharagpur roof (PMAY retrofit)

Same panels and solver, two switches (both in `config/`):

```yaml
pcm_material: n_octadecane    # materials.yaml (paper_cla = the paper's PCM)
loads: kharagpur_roof         # loads.yaml (none = the paper's ΔT set-up)
```

The `octadecane_roof` study sets both for the three shapes. It starts at 06:00 and runs
until the PCM is melted or 18:00, whichever comes first:

```bash
# boundary conditions only (no OpenFOAM): runs/<case>/loads_profile.{csv,png}
python3 -m pcmpipe loads --study octadecane_roof
# 16 GB laptop, 4 cores: screening mode (~85-200k cells, thin-shell resistance)
make docker-study STUDY=octadecane_roof MODE=screening NP=4
# workstation: paper mode (resolved 1 mm ABS shell, ~1 M cells for the smoke mesh, more for medium/fine)
make docker-study STUDY=octadecane_roof MODE=paper NP=16
# any setting for one call, e.g.:
python3 -m pcmpipe all --study octadecane_roof --mode screening --np 4 \
    --set "mesh.wall_layer_passes={air: 0, pcm: 0}" --set run.max_end_time_min=480
```

**Material** (`materials.yaml`). Values are from Faden et al. 2019, *Materials* 12, 2974:
the recommended correlations, evaluated at T_m or 305 K.

| | value | | value |
|---|---|---|---|
| T_m | 301.13 K (28.0 °C) | latent heat | 236.98 kJ/kg |
| solidus / liquidus | 300.63 / 301.63 K (T_m ∓ 0.5 K, model band) | ρ liquid (Boussinesq ρ0) | 776.9 kg/m³ |
| cp solid / liquid | 1921 / 2230 J/(kg K) | k solid / liquid | 0.334 / 0.152 W/(m K) |
| μ (305 K) | 3.63 mPa s | β | 8.7e-4 1/K |

These go into the same dictionaries as the paper's PCM:

```
// constant/pcm/thermophysicalProperties
mixture
{
    specie          { molWeight 254.5; }
    equationOfState { rho 776.9; }                  // rhoConst: Boussinesq rho0
    thermodynamics  { Cp 2230; Hf 0; }              // liquid cp; the fvOption restores the rest
    transport       { mu 0.00363; Pr 53.26; }
}
// constant/pcm/fvOptions
boussinesq
{
    type boussinesqBuoyancySource;  selectionMode all;
    beta 0.00087;  TRef 301.13;  g (0 0 -9.81);    // rho0 beta (TRef - T) g
}
pcmPhaseChange
{
    type pcmEnthalpyPorositySource;  selectionMode all;  energyField h;
    Tsolidus 300.63;  Tliquidus 301.63;  L 236980;  Amush 1e+06;  epsilon 0.001;
    CpSolid 1921;  CpLiquid 2230;  kappaSolid 0.334;  kappaLiquid 0.152;
    latentConvection true;  kinkLimiter true;
}
```

The solid is 12 % denser than the liquid (867.9 kg/m³, the least certain value in the
review). The Boussinesq model has one density, so shrinkage voids are not modelled.

**Loads** (`loads.yaml` → `pcmpipe/loads.py`). The panel is fixed under an uninsulated
125 mm RCC roof slab:
- **Day.** 15 May, from 06:00 IST.
- **Outdoor air.** 27–38 °C (Kharagpur May mean low/high, weatherapi.com), lowest at
  sunrise and highest at 14:30.
- **Solar.** Clear-sky global horizontal irradiation 6.80 kWh/m²/day (NASA POWER
  climatology; all-sky 5.98), spread as a half-sine from sunrise to sunset.
- **Sol-air temperature.** T_sa = T_out + αI/h_o − 3.9 K, with α = 0.7 and h_o = 17 W/(m² K)
  (ASHRAE).
- **Slab.** 1D transient conduction through the slab (k 1.40, ρ 2100, cp 840), with
  h_i = 6.13 W/(m² K) to the room, repeated until the day is periodic.

The result for this design day:

| | range |
|---|---|
| roof top (sol-air) | 23–66 °C |
| ceiling = slab underside | 26.8 °C at 06:00 → 47.7 °C at ~15:45 (3.5 h after the sol-air peak) |
| room air, walls, floor | 27–38 °C (ventilated room: T_room = T_out) |

In the CFD these become tables on the top and bottom of the component:

```
// 0/shell/T (paper mode) and 0/air/T: slab underside
ceiling        { type uniformFixedValue;  uniformValue table ((0 299.98) (600 299.99) ... (86400 299.98)); }
// 0/air/T: the room around the panel
"(walls|floor)" { type uniformFixedValue;  uniformValue table ((0 300.41) ...); }
// 0/pcm/T (screening mode: PCM up to the ceiling behind the 1 mm ABS top wall)
ceiling        { type externalWallHeatFluxTemperature;  mode coefficient;  Ta table (...);
                 h uniform 1e5;  thicknessLayers (0.001);  kappaLayers (0.163);  kappaMethod fluidThermo; }
```

The PCM and shell start solid at 24 °C; the room air starts at T_room(06:00). For a
constant load instead, set `weather.T_max_C` = `T_min_C` and
`daily_global_horizontal_kWh_m2: 0`. For measured weather, give an EPW file
(`weather.source: epw`, e.g. the ISHRAE/TMYx Kolkata file from climate.onebuilding.org).

**Caveats:**
- **Night-time recharge.** Kharagpur May nights stay at or above about 27 °C, so the room
  is below n-octadecane's solidus (27.5 °C) for only about 3 h a day. Without night
  ventilation or cooling the PCM would not re-solidify, and a second hot day would find
  it already liquid. The 06:00 "charged" start is therefore an assumption, and
  `pcmpipe loads` warns about it.
- **Contact.** The panel is assumed to be in perfect thermal contact with the slab.
- **Slab properties.** The slab values are generic dense concrete; replace them with the
  house's construction (add a weathering course or tiles as extra `layers`).

**Outputs**, in `runs/_figures/octadecane_roof[_screening]/`:
- `liquid_fraction_vs_time.csv` and `liquid_fraction_curves.png`
- `mean_pcm_temperature_vs_time.csv` and `mean_pcm_temperature.png`, with the melting range marked
- `liquid_fraction_3d_montage.png`
- per case: `runs/<case>/post/liquid_fraction_3d_<t>min.png`

## 4. How the paper's model maps to OpenFOAM

| Paper | OpenFOAM implementation |
|---|---|
| Air: Eqs. 2–11, incompressible turbulent flow, Realizable k-ε | `air` fluid region: `heRhoThermo` + `rhoConst`, buoyancy ρ0 β (T_ref − T) g from the `boussinesqBuoyancySource` fvOption (Eq. 17, T_ref = ambient), `RAS realizableKE` (C2 = 1.9, σk = 1). No buoyancy production in k-ε, as in Eqs. 7–8 |
| Standard wall functions, Eq. 23 (κ = 0.4187, Cμ = 0.09) | `nutkWallFunction` (Cmu 0.09, kappa 0.4187, E 9.793), `epsilonWallFunction`, `kqRWallFunction`, `alphatJayatillekeWallFunction` (Prt 0.85) |
| ABS shell: Eq. 12 | `shell` solid region, `heSolidThermo` (ρ 1070, k 0.163, cp 1990) |
| PCM: Eqs. 13–15, laminar | `pcm` fluid region: `heRhoThermo` + `rhoConst`, `boussinesqBuoyancySource` (Eq. 16), `simulationType laminar` |
| Boussinesq approximation (density constant except in the buoyancy term, as in Fluent) | `rhoConst` + the buoyancy fvOption, not OpenFOAM's `Boussinesq` EoS. That EoS also puts ρ(T) into continuity and energy; in the sealed capsule and the closed room total mass is then not conserved, and the continuity error showed up as single PCM cells cooling at ~40 K/s (CI diagnostic: of six variants only β = 0 removed it) |
| Liquid fraction, Eq. 1; H = h_sen + γL, Eqs. 18–20; properties linear in the mushy zone | `pcmEnthalpyPorositySource` fvOption (see below) |
| Mushy-zone sink, Eq. 21: A_mush = 1e6, ε = 0.001 | same fvOption, implicit `fvm::Sp` sink |
| Interfaces, Eq. 24: flux continuity, u = 0 | `compressible::turbulentTemperatureCoupledBaffleMixed` on conformal `mappedWall` patches, `noSlip` |
| Adiabatic room walls, adiabatic panel base | `zeroGradient` T on `floor/ceiling/walls`; shell `ceiling` patch adiabatic (`loads: none`; section 3b has the roof set-up) |
| Initial conditions, Eqs. 25–27 | 0/: T_PCM = T_shell = 278.15 K, T_air = 278.15 + ΔT, k = 0.3, ε = 0.5 |
| SIMPLE, second-order upwind, Δt = 0.5 s | PIMPLE (transient CHT), `Gauss linearUpwind`, `deltaT 0.5`, `maxDeltaT 0.5`, `adjustTimeStep` (maxCo 4, maxDi 10) |
| Residuals 1e-6 (energy), 1e-3 (others) | linear-solver absolute tolerances: h 1e-6; U, k, ε 1e-3 (p_rgh 1e-6, see below) |
| ~1.1e6 tetra cells, mesh study at 8.6e5 / 1.1e6 / 1.7e6 | hex-dominant snappy mesh, presets `coarse/medium/fine`; check `log.checkMesh.*` and tune `config/base_case.yaml` |
| 3.51 × 3.16 × 2.25 m room, panel at ceiling centre | blockMesh room, panel at (0, 0, 0) hanging into −z |

### The enthalpy-porosity fvOption

OpenFOAM's thermo package solves for h = cp_ref (T − T_std), with constant cp_ref and
k_ref (the liquid values). `pcmEnthalpyPorositySource` adds the remaining terms of
Eq. 15 exactly:

```
S_h = − ∂(ρE)/∂t − ∇·(φE) + ∇·((k(T) − k_ref)∇T)
E(T) = h_sen(T) − cp_ref·T + γ(T)·L          with h_sen = ∫cp(T)dT, cp and k linear in γ
S_U = − A_mush (1 − γ)² / (γ³ + ε) · U          (Eq. 21, applied as a sink)
```

∂(ρE)/∂t is Newton-linearised in T, so the PIMPLE outer correctors converge the latent
heat. `tests/test_latent_source_1d.py` runs the same scheme in 1D against a fully
converged enthalpy solution:
- with 2 outer correctors the liquid fraction matches to 4 decimals, with 0.000 % energy error
- with 1 corrector the liquid fraction is within 0.7 %

**Kink limiter.** Newton on the piecewise-linear E(T) overshoots when a cell crosses the
solidus or liquidus in one iteration. With the paper's 7 K band this is harmless. For a
pure PCM such as n-octadecane (1 K band, twice the latent heat) plain Newton oscillates:
in the 1D test the stored energy comes out 22.5 % too high. After each energy solve the
fvOption (`correct()`) therefore stops an update that would jump across Tsolidus or
Tliquidus just past the kink, and the next outer corrector re-linearises in the new phase.
With 3 outer correctors and Δt = 0.5 s the n-octadecane error is 0.0003 in liquid fraction
and 0.003 % in energy (the `Octadecane` tests); the paper case is unchanged.
`kinkLimiter false;` in `constant/pcm/fvOptions` switches it off.

The excess conductivity is zero on the interface faces, so the coupled heat flux uses
k_ref = k_liquid. The PCM next to the shell melts first, so this is exact after the first
minutes.

## 5. Assumptions and interpretations (all in `config/physics.yaml`)

The paper does not give these values. Each one is a single line to change.

| Item | Value used | Why |
|---|---|---|
| Solidus / liquidus | 286 K / 293 K | Table 2 gives only T_pc = 286 K. Fig. 8 CFD curves show the latent stage starting at ~13 °C and ending at ~20 °C |
| PCM viscosity | 0.008 Pa·s | not given. Typical for liquid capric/lauric acid mixtures near 20 °C. **Sensitive**: it controls convection in the melt |
| PCM β | 8.5e-4 1/K | not given. Typical for liquid fatty acids |
| Air properties | Fluent defaults: ρ 1.225, μ 1.7894e-5, cp 1006.43, k 0.0242; β = 1/T_ref | the paper used Fluent 2020 R2 |
| σε (Eq. 8) | 1.2 | Eq. 8 prints σk in the ε equation. 1.2 is the Fluent/standard value; set `sigmaEps: 1.0` for a literal reading |
| Panel layout | 2×2 pyramids; 4 tetrahedra on a 60° parallelogram | Fig. 4 / Fig. 16. Reproduces Table 3 areas (228 vs 232 cm², 282 vs 282 cm²) and Sec. 4.5 edge angles exactly |
| Gaps | protrusion bases shrink by g/2 per edge at a fixed footprint; height grows to keep volume | the paper does not say how volume is kept. Equal PCM mass keeps melting times comparable |
| Hemispheres (new) | R_ref = 0.95 × pitch/2; other ratios/gaps give semi-spheroids | not in the paper |
| Shell offset | sharp-join inward offset (exact face distance t; slightly thicker at valleys) | standard CAD "thickness" |
| PCM mass | ρ × CAD cavity volume (≈ 182 g for the 250 mL pyramid) | the paper states 170 g, which is not consistent with 250 mL, a 1 mm shell and ρ = 840 kg/m³. Rates use the actual mass |
| Residuals | p_rgh linear tolerance 1e-6, tighter than the paper's 1e-3 | Fluent's scaled continuity residual is not OpenFOAM's p_rgh residual. A loose p_rgh breaks mass conservation in a closed room |

## 6. Status, known issues and limitations

**Verified** (CI on every push, plus local tests):
- Python side on Linux, macOS and Windows (Python 3.10/3.12): geometry against the paper's
  numbers, the latent-heat scheme against a converged enthalpy solution, case generation for
  every study, and results/plots on synthetic solver output.
- Docker image (OpenFOAM v2406, FreeCAD 0.21, ParaView 5.11).
- The `pcmEnthalpyPorositySource` and `boussinesqBuoyancySource` fvOptions compile against
  OpenFOAM v2406.
- FreeCAD builds valid, watertight panels. For the 250 mL pyramid the cavity is 217.2 mL,
  giving 182.4 g of PCM.
- snappyHexMesh produces the three regions with conformal interfaces (smoke preset: air 323k,
  shell 266k, PCM 360k cells; screening smoke: air 42k, PCM 29k).
- All three smoke runs go end to end (CAD, mesh, 0.6 s of `chtMultiRegionFoam`, results,
  pvpython figures, plots): paper mode (16 min on 2 cores), screening mode (1.5 min) and the
  n-octadecane Kharagpur-roof screening run (`uniformFixedValue` tables and
  `externalWallHeatFluxTemperature` accepted). In all three the PCM minimum stays exactly at
  its initial temperature (`postProcessing/pcm/pcmTMinMax`).

**Fixed (kept here because earlier results are affected):** single PCM cells used to drift
below every temperature around them. In the screening smoke run (paper PCM, ΔT = 20 K) one
cell lost about 40 K/s, reaching 260.9 K at 0.5 s from 278.15 K; in paper mode a run
crashed. A CI diagnostic ran six variants: switching off the latent source, latent
convection or the Darcy sink, or moving the pressure reference cell, changed nothing;
β = 0 removed it exactly. The cause was OpenFOAM's `Boussinesq` equation of state (density
varying in continuity and energy in a sealed region). The cases now use `rhoConst` plus the
`boussinesqBuoyancySource` fvOption (section 4). The same change took the paper-mode air
spike at the pyramid grooves from 124 m/s to 0.67 m/s. **Results computed before commit
69ac77c should be discarded.**

**Known issues** (open):
1. **Air velocity at the groove roots, paper mode.** The air cells at the roots of the
   V-grooves between touching pyramids (gap 0, just below the slab underside) still reach
   0.67 m/s at 0.1 s, decaying to 0.44 m/s at 0.5 s. That is about 10 times the free-fall
   buoyant velocity, so it is a local artefact of the snapped groove-root cells. It holds
   the air Courant limit at Δt ≈ 5 ms. The screening mesh does not show it (4 mm/s).
2. **Cost of paper mode.** Even without that, the 0.5 mm shell cells cap the time step at
   about 0.015 s, against the paper's 0.5 s. Use screening mode for design studies
   (section 3).
3. The medium/fine presets are larger than the paper's 1.1 M cells, because the smoke preset
   alone is already 0.95 M (1.2 M with wall layers). Check `log.checkMesh.*` and tune `config/base_case.yaml`.

**Limitations:**
- Inflation layers are **off in every preset**. They are made with `refineWallLayer` after
  `splitMeshRegions`, which keeps the interfaces conformal, but on these geometries the
  split cells come out crushed. In CI the paper-mode smoke run (1 pass) had 2 highly
  skewed PCM faces and stopped after two time steps with "Negative initial temperature"
  from the Boussinesq thermo. Re-enable with `--set "mesh.wall_layer_passes={air: 1, pcm: 1}"`.
- The solver uses PIMPLE with 2 outer correctors for the paper's PCM and 3 for n-octadecane
  (`solver` block in `materials.yaml`), as the verification tests require. `residualControl` is not used because its syntax differs
  between releases; `solverInfo` logs the residuals.
- The paper's hanging-distance study (Sec. 4.4) is not included.
- Native macOS (OpenFOAM.app, Apple Clang) is not covered by CI; Docker is.
