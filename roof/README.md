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
| 1 | Pure concrete slab (the paper's Ref-SU): geometry, mesh, OpenFOAM case, SLURM job, verification gates | **passed** (tag `roof-phase1`); results in `results/phase1_baseline_concrete/` |
| 2 | 3 shapes × 2 orientations × 5 heights = 30 PCM cases, one SLURM job each | **passed**, all 30 cases (tag `roof-phase2`); results in `results/phase2/` |
| 3 | `pvpython` extraction to small files; paper-style temperature, heat-flux, HGR and thermal-load figures | **done** (tag `roof-phase3`); `results/phase3/`, summary in [Results](#results) |

Every phase was run end to end in the repository's Docker image (OpenFOAM v2406,
FreeCAD 0.21, ParaView 5.11) on a 12-core workstation. The same scripts run on PARAM
Shakti through SLURM.

## Results

Design day, room held at 26 °C, every case against the plain slab at the same time step
(`baseline_dt60`: peak heat gain 101.5 W/m² at 17:15, cooling load 1,253 Wh/m² a day,
peak ceiling 42.56 °C). All 30 cases pass every gate.

| | HGR_max (Eq. 8) | Cooling load cut | Peak ceiling T |
|---|---|---|---|
| Best: pyramid, normal, 1.03 L (`pyr_N_h100`) | 16.2 % | 5.6 % (1,183 Wh/m²) | −2.68 K |
| Tetrahedron, normal, 1.03 L | 16.2 % | 5.3 % | −2.69 K |
| Hemisphere, normal, 1.03 L | 12.3 % | 3.8 % | −2.04 K |
| Pyramid, normal, 0.72 L (the paper's 600 g) | 14.2 % | 4.8 % | −2.35 K |
| Pyramid at 1.03 L **without latent heat** (control) | 7.4 % | 5.7 % | −1.22 K |

- Volume matters most: about +1 point of HGR_max per 150 mL, still rising at 1.03 L.
- At equal volume the tetrahedron and pyramid (base 30,788 mm²) beat the hemisphere
  (15,394 mm²) by 3–4 points of HGR_max; they are within 0.3 points of each other.
- Normal beats inverted by 1.2 points at 40 mm; the gap closes by 100 mm. Small inverted
  bodies melt through (liquid fraction 1.0 from 17:25 to 18:55) at the room's peak.
- The control run shows that the daily load cut comes from OM37's low conductivity
  (0.145 against 1.58 W/m·K), while melting more than doubles the peak cut.
- Every case refreezes completely each night; the normal 600 g pyramid uses 54 % of
  its latent store, the inverted one 86 %.

Full table: `results/phase3/summary_all.csv`; figures `results/phase3/fig1`–`fig7`; the
complete report, with the comparison against the paper's measurements, is
`results/Roof_PCM_Simulation_Report.pdf`.

## Layout

```
roof/
├── config/roof.yaml            every input, tagged [PAPER] / [SOURCE] / [ASSUMPTION]
├── roofpipe/                   case generator: python3 -m roofpipe ... (stdlib + PyYAML)
│   ├── weather.py              design day: solar geometry, GHI, T_out, sol-air, indoor ambient
│   ├── slab1d.py               1D finite-volume reference of the slab (same discretisation)
│   ├── case.py                 Phase 1: renders templates/concrete_slab -> runs/<name>/case
│   ├── geometry.py             Phase 2: the 30 PCM bodies, equal volumes, fit checks
│   ├── pcm_case.py             Phase 2: concrete_slab + pcm_slab -> runs/<id>/case; finalize
│   ├── stl.py                  Phase 2: PCM body STL where FreeCAD is missing
│   └── check.py                gates, last-day time series, daily metrics (both phases)
├── cad/build_roof_pcm_freecad.py   FreeCAD: PCM body solid, STL (m), STEP, volume report
├── templates/
│   ├── concrete_slab/          the OpenFOAM case of the slab (Phase 1, and the base of Phase 2)
│   └── pcm_slab/               overlay: PCM region, its fvOptions, snappyHexMesh, interface BCs
├── scripts/
│   ├── run_slab_case.sh        blockMesh -> decomposePar -> chtMultiRegionFoam -> reconstructPar
│   ├── run_pcm_case.sh         CAD -> snappyHexMesh -> splitMeshRegions -> mass -> solve -> prune
│   └── run_all_local.sh        every Phase 2 case in Docker, N at a time (no SLURM)
├── post/
│   ├── pv_extract.py           pvpython: slices, T contours + melting front, field checks
│   ├── pvpython_headless.sh    pvpython with a private Xvfb on headless machines
│   └── collect.py              every case vs the baseline: tables, time series, figures
├── slurm/
│   ├── env.sh                  OpenFOAM + Python environment (edit once)
│   ├── phase1_baseline.sbatch  Phase 1 job: generate, run, check
│   ├── phase2_case.sbatch      one PCM case (CASE_ID)
│   ├── phase2_submit_all.sh    builds the PCM library, submits every case not yet passed
│   └── phase3_post.sbatch      pv_extract on every case, then collect
├── tests/om37_timestep.py      time-step check of the PCM model against the converged method
└── results/                    kept in git: phase1_*/, phase2_baseline_dt60/, phase2/<case>/, phase3/
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

**Result** (Docker, OpenFOAM v2406, 4 ranks, 87 s; `results/phase1_baseline_concrete/`):
all four gates pass. Periodic to 0.015 K, 1.4 × 10⁻⁵ K from the 1D reference, energy
balance closed to 8 × 10⁻⁸. Ceiling 27.2–42.6 °C, peaking at 17:15 (time lag 4.1 h,
decrement factor 0.33); room heat gain 7.6–101.6 W/m², mean 52.2; daily cooling load
1,253 Wh/m², no heating load (the ceiling never drops below 26 °C).

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

## Phase 2: 30 PCM cases

### Design

| | |
|---|---|
| Shapes | **tetrahedron**: equilateral-triangle base, apex above its centroid · **pyramid**: square base · **hemisphere**: half-spheroid on a circle of r = 70 mm, a true hemisphere at h = 70 mm |
| Orientations | **N**ormal: base 25 mm above the ceiling face, apex up · **I**nverted: base 25 mm below the roof face, apex down |
| Heights | 40, 55, 70, 85, 100 mm. The slab stays 150 mm, with 25 mm of concrete cover above and below the body |
| Equal volumes | V = A·h/3 for all three shapes when the tetrahedron and pyramid bases share A = 2πr² = 30,788 mm² (pyramid 175.5 mm square, tetrahedron 266.6 mm triangle). V = 410.5, 564.4, 718.4, 872.3, 1026.3 mL at every shape; 718.4 mL × 860 kg/m³ = 618 g is the paper's 600 g per void former |
| PCM | OM37 (Table 1): melting band 36.4–44.06 °C (DSC onset/end, Sec. 2.2), L = 179 kJ/kg, cp 2550/2630 J/kg·K, k 0.16/0.13 W/m·K (solid/liquid), ρ = 860 kg/m³ (liquid: a capsule filled with melt) |
| PCM model | a solid region with conduction and latent heat: the repository's enthalpy fvOption `pcmEnthalpyPorositySource` with `latentConvection false`. Natural convection in the melt is neglected: heated from above, the melt lies on top of the solid (stable). No aluminium shell: the PCM touches the concrete |
| Mesh | blockMesh 61 × 61 × 31 (≈ 4.9 mm), snappyHexMesh level 1 (2.5 mm) inside and on the PCM, level 2 (1.2 mm) on its edges; `splitMeshRegions -cellZones` makes the regions `concrete` and `pcm` with a conformal interface. The background counts keep every cell face off the base and apex planes. About 216,000 cells |
| Equal mass | after meshing, `roofpipe finalize` scales the PCM density by V_exact / V_mesh (snapping changes the volume by up to ±0.2 %), so every shape at a height holds the same mass |
| Time | 3 days, Δt = 60 s, 2 outer correctors. `tests/om37_timestep.py`: within 0.002 in liquid fraction and 0.06 % in stored energy of the converged enthalpy method (30 s: 0.001 / 0.03 %). The cases are compared with a pure slab run at the same step, `baseline_dt60`, which differs from the 30 s baseline by 0.01 K in peak ceiling temperature and 0.1 W/m² in peak gain |
| Gates | complete; periodic (both surfaces and the mean PCM temperature, 0.05 K); energy balance including the PCM's sensible and latent heat (1 %) |

The Phase 1 baseline is the reference for every case. Its 2.5 mm vertical cells against
the ≈ 4.9 mm background of Phase 2 change the surface temperatures by 0.002 K and the
cooling load by 0.00 % (1D check), so the comparison is not affected by the mesh.

### Run it

On PARAM Shakti, from `roof/`, after Phase 1 (same `slurm/env.sh`; set the partition in
`slurm/phase2_case.sbatch`):

```bash
bash slurm/phase2_submit_all.sh          # compiles the PCM library once, then 30 jobs
squeue -u $USER                          # logs in roof/logs/
DRY=1 bash slurm/phase2_submit_all.sh    # only print the sbatch lines
bash slurm/phase2_submit_all.sh pyr_     # only some cases (ids containing "pyr_")
```

A case is about 216,000 cells and takes about 6 min on 4 cores of a workstation at
Δt = 60 s (2 h are requested). Cases that passed are skipped on resubmission; a stopped job continues.

On a workstation with Docker, from the repository root:

```bash
JOBS=3 NP=3 bash roof/scripts/run_all_local.sh      # all cases, 3 at a time
```

One case by hand (from `roof/`):

```bash
python3 -m roofpipe cases                             # the 30 ids and their volumes
python3 -m roofpipe pcm --case pyr_N_h070 --np 4
bash scripts/run_pcm_case.sh runs/pyr_N_h070 4        # MESH_ONLY=1: stop after the mesh
python3 -m roofpipe check pyr_N_h070
```

Each case keeps `cad/` (STL, STEP, `cad_report.json`), `case/` (the fields of the
analysed day only), `case_meta.json` (with the meshed volume and the PCM mass) and
`results/` (`summary.json`, `timeseries_lastday.csv`, `quicklook.png`).

## Phase 3: extraction and figures

```bash
sbatch slurm/phase3_post.sbatch                 # PARAM Shakti (PURGE=1: delete fields after)
# or, per case and then all cases (from roof/):
bash post/pvpython_headless.sh post/pv_extract.py runs/pyr_N_h070
python3 post/collect.py
```

The time series are not extracted from stored fields: the solver's function objects
write them every 300 s while it runs (a few hundred kB per case). `pvpython` works only
on the analysed day's field snapshots (every 3 h), for what needs the 3D fields:

- `runs/<id>/results/pv/slice_{concrete,pcm}_<hh>h.csv`: the vertical mid-plane through
  the PCM axis, T and liquid fraction at cell centres
- `runs/<id>/results/pv/contour_T_<hh>h.png`: T on one fixed 25–55 °C scale for every case
  and time, with the PCM outline and its melting front (liquid fraction 0.5)
- `runs/<id>/results/pv/surface_metrics.csv`: roof/ceiling T and heat flux from the fields,
  a check on the function objects

`post/collect.py` writes `results/phase3/`:

| File | Content |
|---|---|
| `summary_all.csv` | one row per case and the baseline: PCM volume and mass, peak heat gain and its time, HGR_max (Eq. 8), daily cooling/heating/total load (Eq. 9) and its reduction, peak ceiling T and its reduction, time lag, decrement factor, liquid-fraction range, latent heat cycled |
| `timeseries_T_roof_C.csv`, `timeseries_T_ceiling_C.csv` | roof and ceiling surface temperature, one column per case |
| `timeseries_q_roof_in_W_m2.csv`, `timeseries_q_ceiling_gain_W_m2.csv` | outdoor heat flux (into the slab) and indoor heat flux (into the room) |
| `timeseries_load_cumulative_Wh_m2.csv`, `timeseries_liquid_fraction.csv` | cumulative cooling load; PCM liquid fraction |
| `fig1_temperature_h070.png` | temperature variation vs time, T_amb, T_ext, T_int,Ref, T_int,Exp (paper Figs. 7–8) |
| `fig2_heatflux_h070.png` | indoor heat flux vs time (Fig. 9) |
| `fig3_thermal_load.png` | daily cooling load per case against the plain slab (Fig. 13), and the load each case saves |
| `fig4_hgr_vs_volume.png` | maximum heat-gain reduction vs PCM volume (Eq. 8, Fig. 12) |
| `fig5_liquid_fraction_h070.png` | melting and solidification over the day |
| `fig6_ceiling_T_all.png` | ceiling temperature of all 30 cases |
| `fig7_melting_front_pyr_h070.png` | mid-plane T contours with the melting front, normal and inverted pyramid at 70 mm, 12:00–21:00 (composed from the case contours) |
| `contours/<case>/` | every case's `pv_extract` contours (every 3 h) and `surface_metrics.csv`; the slice CSVs stay in `runs/` |
| `control_pyr_N_h100_no_latent/` | the best design with OM37's latent heat set to 100 J/kg (`--set materials.OM37.latent_heat=100`) |
| `specimen_*` | the paper's test-unit set-up (`--set indoor.mode=specimen`): plain slab and the 600 g pyramid |
