# Figures

Figures used in the mid-evaluation report and slides, kept in the repository so they do
not depend on a local run or on CI artifacts (GitHub deletes those after 90 days). Each
section says how to regenerate its figures.

## design/: no CFD needed

Regenerate with `python3 post/presentation_figures.py` (writes `runs/_figures/presentation/`).

| File | Shows |
|---|---|
| `design_catalogue.png` | the three shapes × four design variants (reference 250 mL, 70 % protrusion, 5 mm gap, 750 mL), common scale per row |
| `geometry_validation.png`, `.md` | our parametric geometry against every geometric number the paper publishes |
| `latent_heat_verification.png` | 1D PCM slab: the OpenFOAM phase-change source term against a converged enthalpy solution |

## loads/: Kharagpur roof boundary conditions

Regenerate with `python3 -m pcmpipe loads --study octadecane_roof --mode screening --case pyr_`
(writes `runs/<case>/loads_profile.{png,csv}`; no OpenFOAM needed).

| File | Shows |
|---|---|
| `kharagpur_roof_loads.png`, `.csv` | 24 h from 06:00 IST on a hot May day: sol-air temperature, ceiling (slab underside) and room temperatures applied in the CFD, global horizontal irradiance, and the n-octadecane melting range |

## pipeline_check/: CI smoke runs, 0.6 s of flow

From the CI run of commit d353e6f (GitHub Actions run 40457762599, 27 Oct 2026), smoke
meshes on 2 cores. They show that the whole chain runs: three coupled regions, the first
conduction into the PCM and the start of the air flow. **0.6 s is far too short for
melting: the liquid fraction is 0 everywhere, so these are not results.** In all three runs
the lowest PCM temperature stays at its starting value (278.15 K, 297.15 K). Regenerate with
`make docker-smoke MODE=paper`, `MODE=screening`, or `MODE=screening SMOKE=smoke_roof`.

| File | Shows |
|---|---|
| `paper_mode_T_0.6s.png`, `paper_mode_U_0.6s.png` | paper mode (resolved 1 mm shell, 0.95 M cells): temperature and air speed on the symmetry plane |
| `paper_mode_liquid_fraction_3d_0.6s.png` | the same run: 3D cut-away of the PCM coloured by liquid fraction (all solid) |
| `screening_T_0.6s.png`, `screening_liquid_fraction_3d_0.6s.png` | screening mode (shell as a thin resistance, 71 k cells), paper set-up |
| `kharagpur_roof_T_0.6s.png`, `kharagpur_roof_liquid_fraction_3d_0.6s.png` | screening mode, n-octadecane under the Kharagpur roof at 06:00 (colour range 24–40 °C) |

The 2D contour titles in these images sit over the panel; `post/pv_contours.py` now leaves a
white band above the ceiling for them.
