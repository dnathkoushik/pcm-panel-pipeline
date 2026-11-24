# PCM embedded in a concrete roof slab

This directory holds the roof-embedded continuation of the project. PCM panels sit
**inside a solid concrete roof slab** (pure concrete, no air gap), and the slab is
simulated over the 24-hour thermal cycle of a design day. It replicates and extends

> P. Jaffar Abass, S. Muthulingam (2024), *Comprehensive assessment of PCM integrated
> roof for passive building design: A study in energo-economics*, Energy & Buildings
> 317, 114387, doi:10.1016/j.enbuild.2024.114387

The ceiling-panel pipeline in the repository root (`pcmpipe/`, after Yang et al. 2023)
is unchanged. This directory only reuses its solver choice (openfoam.com
`chtMultiRegionFoam`) and the boundary-condition syntax its CI verified on v2406.

| Phase | Content | Status |
|---|---|---|
| 1 | Pure concrete slab (the paper's Ref-SU): geometry, mesh, OpenFOAM case, SLURM job, verification gates | **written, waiting for the PARAM Shakti run** |
| 2 | 3 shapes × 2 orientations × 5 volumes = 30 PCM cases, one SLURM job each | after Phase 1 passes |
| 3 | `pvpython` extraction to small CSVs; paper-style temperature, heat-flux and thermal-load figures | after Phase 2 |

## Layout

```
roof/
├── config/roof.yaml            every input, tagged [PAPER] / [SOURCE] / [ASSUMPTION]
├── roofpipe/                   case generator: python3 -m roofpipe ... (stdlib + PyYAML)
│   ├── weather.py              design day: solar geometry, GHI, T_out, sol-air, indoor ambient
│   ├── slab1d.py               1D finite-volume reference of the slab (same discretisation)
│   ├── case.py                 renders templates/concrete_slab -> runs/<name>/case
│   └── check.py                Phase 1 gates, last-day time series, daily metrics
├── templates/concrete_slab/    the OpenFOAM case; placeholders are filled from the config
│   ├── 0/concrete/{T,p}
│   ├── constant/{regionProperties,g}, constant/concrete/{thermophysical,radiation}Properties
│   └── system/{blockMeshDict,controlDict,fvSchemes,fvSolution,decomposeParDict}, system/concrete/*
├── scripts/run_slab_case.sh    blockMesh -> decomposePar -> chtMultiRegionFoam -> reconstructPar
└── slurm/
    ├── env.sh                  OpenFOAM + Python environment (edit once)
    └── phase1_baseline.sbatch  Phase 1 job: generate, run, check
```

## Phase 1: the pure concrete baseline

### Model

| | |
|---|---|
| Domain | the paper's unit cell, 300 × 300 × 150 mm (Sec. 2.4, Fig. 4). z = 0 is the ceiling face, z = 150 mm the roof face |
| Mesh | blockMesh, 30 × 30 × 60 hexahedra (10 × 10 × 2.5 mm), 54 000 cells, copied to the region `concrete` |
| Solver | `chtMultiRegionFoam`, one solid region (`heSolidThermo`, constant properties). Phase 2 adds the PCM regions to this set-up |
| Concrete | M20 RCC with the rebar smeared in: k 1.58 W/m K, ρ 2288 kg/m³, cp 880 J/kg K (ECBC 2017 / SP 41; the paper gives no values) |
| Roof face | `externalWallHeatFluxTemperature`, `mode coefficient`: q = h_out (T_sa − T_s), h_out = 17 W/m² K. T_sa = T_out + α I / h_out − 3.9 K (ASHRAE sol-air), α = 0.7 |
| Weather | Rupnagar (30.97 N, 76.52 E), 5 June. T_out: 27 °C at sunrise to 43 °C at 15:00 (the paper's summer average, Sec. 2.1). GHI: ASHRAE clear-sky shape scaled to 6.67 kWh/m²/day (NASA POWER June mean, all-sky; its annual mean, 5.22, is the paper's figure) |
| Ceiling face | q = h_in (T_s − T_in), h_in = 6.13 W/m² K (ASHRAE, heat flow downward). `indoor.mode`: **conditioned** room at 26 °C (default) · **ventilated** room at T_out · **specimen**, the paper's rig: 25 mm XPS under the slab, shaded outdoor air below |
| Sides | adiabatic (insulated test unit, or symmetry planes of a repeated slab) |
| Time | 3 days of the same design day, Δt = 30 s fixed, implicit Euler. Days 1-2 wash out the uniform initial state; day 3 is analysed |
| Output | function objects only, every 300 s: surface T of both faces, wall heat flux, mean slab T, T(z) at the 60 cell centres hourly, min/max, residuals. Full fields every 6 h, for restarts and Phase 3 contours |

### Run it on PARAM Shakti

1. **Clone** to a path without spaces (OpenFOAM cannot handle spaces), and check out
   this branch:
   ```bash
   git clone -b roof-embedded https://github.com/dnathkoushik/pcm-panel-pipeline.git
   cd pcm-panel-pipeline/roof
   ```
2. **Environment.** Edit the settings block of `slurm/env.sh`. You need:
   - an **openfoam.com** release, v2306 or newer (`module avail openfoam`; set
     `OPENFOAM_MODULE` or `OPENFOAM_BASHRC`);
   - `python3` with **PyYAML** (`PYTHON_MODULE=...`, or `pip install --user pyyaml`).

   Check it on the login node: `bash -c 'source slurm/env.sh && roof_check_env'`.
3. **Partition.** Uncomment `#SBATCH --partition=` (and `--account=` if needed) in
   `slurm/phase1_baseline.sbatch`; `sinfo -s` lists the partitions.
4. **Submit** from `roof/`:
   ```bash
   sbatch slurm/phase1_baseline.sbatch
   # the paper's insulated test rig, as a second baseline:
   sbatch --export=ALL,CASE_NAME=baseline_specimen,ROOF_SET="indoor.mode=specimen" slurm/phase1_baseline.sbatch
   ```
   The job writes the case, runs `scripts/run_slab_case.sh` and then
   `python3 -m roofpipe check`. Its log is `roof-p1-baseline-<jobid>.out` in `roof/`.

To generate a case without running it (e.g. to inspect the dictionaries), run
`python3 -m roofpipe baseline --name test --np 8` from `roof/`.

### Phase 1 is successful when

The job log ends with `PHASE 1 PASSED`. The four gates are:

| Gate | Test | Limit |
|---|---|---|
| complete | the solver reached endTime (`End` in `log.chtMultiRegionFoam`) | |
| periodic | surface temperatures of day 3 repeat those of day 2 | 0.05 K |
| reference | roof and ceiling surface temperatures match `reference_1d.csv`, a 1D finite-volume solution with the same cells, time step, tables and film coefficients. The 3D slab is one-dimensional by symmetry, so a deviation means a wrong boundary condition, table, property or mesh | 0.05 K |
| energy | day 3: roof heat in − ceiling heat out − change in stored heat | 1 % of the gross roof heat flow |

**Sanity anchor (conditioned mode).** For a linear, periodic problem the daily-mean
flux equals the steady flux under daily-mean conditions:
U = 1 / (1/17 + 0.15/1.58 + 1/6.13) = 3.16 W/m² K, mean T_sa = 35.0 + 0.7 × 277.9 / 17 − 3.9 = 42.5 °C,
so the mean ceiling gain is about 3.16 × (42.5 − 26) ≈ 52 W/m², or ≈ 1250 Wh/m² per day,
and the mean ceiling surface temperature is about 26 + 52 / 6.13 ≈ 34.5 °C.

Please send back `runs/<name>/results/summary.json`, `results/quicklook.png` (if
matplotlib is available on the cluster) and, if a gate fails, the tail of
`case/log.chtMultiRegionFoam`.

### Files of a run

```
runs/<name>/
├── case/                      OpenFOAM case; case.foam opens it in ParaView
│   └── postProcessing/        roofSurface, ceilingSurface, surfaceHeatFlux, slabAverage,
│                              thicknessProfile, concreteTMinMax, concreteResiduals
├── boundary_conditions.csv    T_out, GHI, T_sa, T_in on the run's time axis (the 0/ tables)
├── reference_1d.csv           1D solution, every 300 s
├── case_meta.json             mesh, coefficients, sun times, initial state, U-value
├── config_used.yaml           configuration after the --set overrides
└── results/                   from `roofpipe check`
    ├── summary.json           gates, daily metrics (peak T and times, decrement factor,
    │                          time lag, cooling/heating loads Wh/m², stored-heat change)
    ├── timeseries_lastday.csv the analysed day (temperatures, heat fluxes, 1D reference,
    │                          the paper's Eq. 3 flux K_eff (T_ext − T_int))
    └── quicklook.png
```

Sign convention: `q_roof_in` > 0 heats the slab from outside; `q_ceiling_gain` > 0 is
heat gained by the room (cooling load), < 0 heat it loses (heating load). The thermal
load is the time integral of `q_ceiling_gain` (paper Eq. 9), split by sign.

### Container instead of a module

If PARAM Shakti has no openfoam.com module, use the repository's Docker image through
Apptainer (single node, MPI inside the container). On a machine with Docker, from the
repository root:
```bash
make docker-build
docker save pcm-panel-pipeline:latest -o pcm.tar
```
Copy `pcm.tar` to the cluster, then `apptainer build pcm.sif docker-archive://pcm.tar`,
and set `PCM_SIF=/path/to/pcm.sif` in `slurm/env.sh`.

### Assumptions to confirm

Every value is one line in `config/roof.yaml`.

| Item | Value | Note |
|---|---|---|
| Concrete k, ρ, cp | 1.58, 2288, 880 | not in the paper. Measured M20 values would be better |
| h_ext, h_int | 17 and 6.13 W/m² K | Eq. 4 uses them but the paper gives no values (ASHRAE) |
| Roof absorptance | 0.70 | bare weathered concrete |
| Design day | 43/27 °C, 6.67 kWh/m² | the paper's Fig. 7 measurements are not tabulated |
| Indoor side | conditioned room at 26 °C | the thermal load is then the HVAC load. The paper's specimen was insulated below: `indoor.mode=specimen` |
| Rebar | smeared into the concrete | two 10 mm meshes, Sec. 2.4 |
