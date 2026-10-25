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
│   └── studies.yaml          # parametric studies (cartesian products of overrides)
├── pcmpipe/                  # Python driver:  python3 -m pcmpipe <step> --study <name>
│   ├── __main__.py           # CLI
│   ├── config.py             # YAML loading, study expansion, case ids
│   ├── geometry.py           # parametric dimension solver (pure Python)
│   ├── cad.py                # calls FreeCAD
│   ├── mesh.py               # blockMesh/snappyHexMesh dictionaries + meshing sequence
│   ├── case_setup.py         # constant/ and system/ for air | shell | pcm
│   ├── fields.py             # 0/ initial and boundary conditions
│   ├── run.py                # decomposePar, mpirun chtMultiRegionFoam, reconstructPar
│   ├── results.py            # melting times, rates, normalised rates
│   └── foam.py               # dictionary writer, command runner
├── cad/build_panel_freecad.py
├── openfoam/
│   ├── Allwmake
│   └── src/pcmEnthalpyPorositySource/   # enthalpy-porosity fvOption (Eqs. 1, 15, 18-21)
├── post/
│   ├── pv_contours.py        # pvpython: T and |U| contours at 50/150/250/300 min (Fig. 9)
│   ├── pv_extract_transient.py  # pvpython: liquid fraction / mean T for every written time
│   └── plot_results.py       # paper-style figures (Figs. 8, 10-13, 15, 17, 18)
├── validation/paper_reference.json   # numbers quoted in the paper, overlaid on plots
├── tests/
│   ├── test_geometry.py      # geometry solver vs. the paper's Table 3, Sec. 4.5, Sec. 4.7
│   ├── test_latent_source_1d.py  # latent-heat formulation vs. converged enthalpy method
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

## 4. How the paper's model maps to OpenFOAM

| Paper | OpenFOAM implementation |
|---|---|
| Air: Eqs. 2–11, incompressible turbulent flow, Realizable k-ε | `air` fluid region: `heRhoThermo` + `Boussinesq` EoS (Eq. 17, T_ref = ambient), `RAS realizableKE` (C2 = 1.9, σk = 1). No buoyancy production in k-ε, as in Eqs. 7–8 |
| Standard wall functions, Eq. 23 (κ = 0.4187, Cμ = 0.09) | `nutkWallFunction` (Cmu 0.09, kappa 0.4187, E 9.793), `epsilonWallFunction`, `kqRWallFunction`, `alphatJayatillekeWallFunction` (Prt 0.85) |
| ABS shell: Eq. 12 | `shell` solid region, `heSolidThermo` (ρ 1070, k 0.163, cp 1990) |
| PCM: Eqs. 13–15, laminar | `pcm` fluid region: `heRhoThermo` + `Boussinesq` EoS (Eq. 16), `simulationType laminar` |
| Liquid fraction, Eq. 1; H = h_sen + γL, Eqs. 18–20; properties linear in the mushy zone | `pcmEnthalpyPorositySource` fvOption (see below) |
| Mushy-zone sink, Eq. 21: A_mush = 1e6, ε = 0.001 | same fvOption, implicit `fvm::Sp` sink |
| Interfaces, Eq. 24: flux continuity, u = 0 | `compressible::turbulentTemperatureCoupledBaffleMixed` on conformal `mappedWall` patches, `noSlip` |
| Adiabatic room walls, adiabatic panel base | `zeroGradient` T on `floor/ceiling/walls`; shell `ceiling` patch adiabatic |
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

## 6. Limitations and what was not verified

- **None of the OpenFOAM, FreeCAD or ParaView steps have been executed yet.** The machine
  where this was written had no OpenFOAM, FreeCAD or ParaView. What was tested:
  - the geometry solver (`tests/test_geometry.py`, against the paper's numbers)
  - the latent-heat discretisation (`tests/test_latent_source_1d.py`)
  - dictionary generation for all 39 cases (`--dry-run` / `setup`)
  - results parsing and every plot, on synthetic data

  Expect to fix small version-specific issues on the first `smoke` run. Most likely:
  1. **Compiling the fvOption** (`openfoam/Allwmake`). Written against the openfoam.com
     v2306+ `fv::cellSetOption` API. Check `wmake` output first.
  2. **Multi-region snappyHexMesh.** Zones come from `locationsInMesh` (the surfaces carry no faceZone; snappy creates them). If
     `splitMeshRegions` does not produce exactly `air`, `shell`, `pcm`, open the snappy
     log and `checkMesh` output. Raise `shell_cells_across` if the 1 mm shell
     breaks up.
  3. **`refineWallLayer`** (the inflation layers) skips cells it cannot cut and reports them
     in `log.<region>.refineWallLayer.*`. Set `wall_layer_passes` to 0 to disable it.
- Inflation layers are made by splitting the wall-adjacent cells (`refineWallLayer`) after
  `splitMeshRegions`. This keeps the shell/fluid interfaces conformal. snappyHexMesh's own
  layer addition does not reliably apply to internal faceZones between regions.
- Compute cost: at the medium preset, a 250 mL case is ~1.1 M cells and up to 48 000 time
  steps at Δt = 0.5 s. Expect roughly one to two days on 8 cores per case. The full study
  set is ~39 unique cases; run `validation` first and compare with
  `validation/paper_reference.json`.
- The paper's hanging-distance study (Sec. 4.4) is not included.
- `residualControl` for chtMultiRegionFoam outer loops is not used, because its syntax
  differs between releases. Convergence comes from 2 outer correctors plus the tolerances
  above. `solverInfo` logs the residuals for checking.
