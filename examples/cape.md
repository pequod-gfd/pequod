# Wind-driven rectangular cape

The cape is a land rectangle attached to the north wall and extending south. It leaves one connected coastline. The three-layer prototype uses the qg3 stratification, CMM modal inversion and baroclinic mass constraints, with one outer-wall constant per layer. The cape has axial edges and creates wet-cell fractions of 1/4, 1/2 and 3/4 along the coast.

The prototype starts from rest and evolves relative PV with beta-plane advection, top-layer wind curl and lateral vorticity diffusion. It supports first-order upwind and experimental CABARET advection, a cut-cell outgoing CFL limit and implicit viscosity. A wet cell is the water-covered part of the control volume around a grid node. A wet face is the open portion of an edge shared by two such volumes. All walls are stationary and impermeable. Viscous runs apply finite Navier slip on the outer walls and the cape.

## Build and run

From the repository root, with GNU Fortran, MPI, BLAS and LAPACK installed:

```sh
make -j4 cape
make run-cape CAPE_RANKS=4 CAPE_RUN_DIR=run-cape-01
```

The first command builds `build/cape-native/qg3_cape` in an isolated native build directory. The second runs the supplied `examples/cape.nml`. To supply a different namelist, use a fresh run directory:

```sh
cp examples/cape.nml my-cape.nml
# Edit my-cape.nml.
make run-cape CAPE_CONFIG=my-cape.nml CAPE_RANKS=4 CAPE_RUN_DIR=run-cape-02
```

To launch the executable directly from the repository root:

```sh
mkdir run-cape-03
(cd run-cape-03 && mpirun -np 4 ../build/cape-native/qg3_cape ../examples/cape.nml)
```

## Namelist

| Parameter | Meaning |
| --- | --- |
| `nx`, `ny` | Embedding node counts east–west and north–south |
| `cape_start` | Number of grid intervals east of the northwest corner to the cape's west edge |
| `cape_width` | Number of grid cells occupied eastward along the north wall |
| `cape_height` | Number of grid cells occupied southward from the north wall |
| `basin_km` | North–south span in km. The east–west grid uses the same spacing |
| `wind_profile` | `'cape_single'` for the localized gyre east of the cape or `'double_gyre'` for the basin-wide profile |
| `wind_stress` | Zonal stress amplitude per unit density in cm²/s² |
| `wind_shift`, `wind_tilt`, `wind_asymmetry` | Additional parameters for `'double_gyre'` |
| `beta_gradient` | Planetary vorticity gradient in cm⁻¹ s⁻¹ |
| `advection` | `'upwind'` by default or `'cabaret'` for staggered face PV |
| `run_days`, `max_dt_seconds` | Run length in days and maximum step size in seconds. The outgoing CFL may reduce the step |
| `viscosity_m2_s` | Lateral viscosity in m²/s. Default 100. Set it to zero for the inviscid experiment |
| `slip_length_km` | Navier slip length in km. Default 120. Viscous runs require a finite positive value |
| `snapshot_days` | Time between qg3 snapshots in model days. Default 1. The initial and final states are also saved |

The cape must stay inside both side walls and leave a strict-fluid row south of its tip. Grid spacing is `basin_km/(ny-1)`. The supplied 129×129 configuration has 30 km spacing, attaches 720 km east of the northwest corner, and removes a 480 km × 480 km cape.

## CABARET

Add `advection = 'cabaret'` inside `&cape` in a copy of the namelist, then pass it with `CAPE_CONFIG`. Start in a fresh run directory. The same snapshot files work with `plots/movies.py`.

Full cells extrapolate from the opposing-face history and the half-step cell PV. Cut cells use the same half-step state with a fit to actual wet centroids and incoming face history. The limiter bounds departures from a locally fitted affine characteristic. Each face owner integrates one shared PV flux. If velocity reverses along the face, the integral uses the appropriate upstream trace on each portion. The driver retries steps that exceed an outgoing CFL of 0.45.

CABARET advances the cell state in two half steps, using the old face flux in the predictor and the new face flux in the corrector. Both half steps apply wind and use Crank–Nicolson for viscosity. The face extrapolation adjusts for the first half-step source, including the wall PV projection, and estimates the second half-step source. This standalone cape integrator diagnoses velocity and wall PV at the half step and at the accepted final state. It uses the two most recent half-step velocities to extrapolate the new face velocity on unequal time steps.

The limited scheme reproduced an affine transport trajectory to roundoff. On a smooth manufactured trajectory with spatial and time refinement, the finest observed RMS orders were about 1.86 over the basin and 1.80 in cut cells and their immediate neighbors. In the coupled viscous run, time refinement gave about second order in streamfunction and variable order in PV at very small CFL. The measured stability covers the tested flows and step sizes.

## Wind forcing

The default `wind_profile = 'cape_single'` applies steady top-layer wind from the cape's east face to the east wall and from the latitude of the cape tip to the north wall. The patch follows the cape automatically when its dimensions or attachment change. In the supplied configuration it spans 1200–3840 km east and 3360–3840 km north, measured from the southwest corner. The model prints these bounds at startup.

```fortran
  wind_profile = 'cape_single'
  wind_stress = 0.8
```

Within the patch, let `xi` and `eta` run from zero to one from west to east and from south to north. The zonal stress is `-wind_stress * sin(pi*xi)**2 * sin(pi*eta/2)**2`, with zero meridional stress. Its curl has one sign inside the patch and vanishes outside. Stress tapers to zero at the west, east and south edges and is nonzero along the northern wall, supplying the net wind input in the circulation budget. A negative `wind_stress` reverses the curl and zero switches off the wind. The resulting circulation can extend beyond the forcing patch.

To change the zonal shape within the same patch, edit `cape_gyre_stress` in `qg3p_irregular_forcing.F90` and rebuild with `make cape`. The helper converts a point to patch coordinates and evaluates the profile with the model's acceleration scaling. `cape_gyre_line` clips the integration to the patch, and `accumulate_wind_curl` uses the same stress integral on shared faces and the coast. If a new formula has a breakpoint inside the patch, split the segment integral there.

The earlier basin-wide forcing is available with `wind_profile = 'double_gyre'`. For that profile, `wind_shift` and `wind_tilt` set the divide at `y = wind_shift + wind_tilt*x`, where eastward `x` and northward `y` each range from −0.5 to +0.5. Require `abs(wind_shift) + abs(wind_tilt)/2 < 0.5`. Increasing the positive `wind_asymmetry` multiplies the southern branch and divides the northern branch by the same factor. A value of 1 gives equal amplitude multipliers.

## Viscosity and slip

The viscous term diffuses relative vorticity. It uses the same symmetric, length-weighted open-face difference as the discrete Poisson balance. The Navier law relates integrated wall circulation to `−slip_length × wall integral of vorticity`. The cape scheme approximates that integral as the wet-cell vorticity value times the sum of its two wall half-segment lengths. It uses each cell's own wet area, including at the cape's convex and reentrant corners. This pairing gives a discrete energy dissipation identity. The upwind option uses Heun advection with an implicit predictor and a Crank–Nicolson diffusion corrector.

The wall PV input printed by the model comes from the cut-cell budget. The wall integral is a low-order quadrature. Pointwise convergence near reentrant corners remains open. The prototype starts from rest.

## Outputs

- `cape_geometry.dat` holds the flagged, checksummed geometry. Validate it with `python3 scripts/qg_irregular --check run-cape-01/cape_geometry.dat`.
- `qg3_<ny>_s.<index>` and `qg3_<ny>_s.<index>.d` are standard single-level qg3 primary snapshot headers and data. MPI runs append the rank number. They contain the two timestamps, three streamfunction fields and three relative PV fields in native units. Relative PV is a wet-cell mean, including along the coast. Land slots contain zero and the movie script reads the geometry to shade them grey.
- `cape_state.csv` holds the final wet-node streamfunction in m²/s and wet-cell mean relative PV in s⁻¹. It includes global indices, layer, position in km, node kind and wet fraction.
- Standard output reports wet area, elapsed time, peak outgoing CFL, transport divergence, baroclinic mass residual and integrated PV budget residual. Viscous runs also report the wall-law residual and diagnosed wall PV input by layer. CABARET reports CFL rejections and the largest per-step count of cut traces that reverted to the half-level mean because their fit was rank-deficient.

Render the primary snapshots with Python, NumPy, Matplotlib and FFmpeg.

```sh
python3 plots/movies.py run-cape-01 --fps 3
```

This creates `fields.png` and `fields.mp4` in the run directory. Use a fresh run directory for each experiment. The model refuses to overwrite existing output files. The primary snapshots are visualization files for this experiment. Rectangular qg3 pickup and time-average outputs use additional files.
