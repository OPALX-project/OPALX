## Integrator comparison (2026-09-07)

Historical experiment: the Boris new-search columns use the unrestricted secant version. The final implementation restores Boris bisection; RK4/DOP853 retain secant. See README commit-scope notes and commit-boris-fallback results.

Stable hard-edge DBA; current Release executable, one MPI rank and four OpenMP threads.
All 15 COF runs verified the zero closed orbit with unchanged solver tolerances.
Tables report principal fractional phases; integer and oriented branches are not determined.

### Matched COF setting: DT=5e-12 s, h=3e-5

| Integrator | Qx | Qy | max absolute phase error | max transverse matrix error | Classification |
|---|---:|---:|---:|---:|---|
| BORIS | 0.229891459894 | 0.310462384101 | 1.219e-06 | 1.644e-05 | STABLE |
| RK4 | 0.229892431323 | 0.310463603175 | 4.555e-10 | 4.687e-09 | STABLE |
| DOP853 | 0.229892431493 | 0.310463603175 | 6.263e-10 | 6.312e-09 | STABLE |

### Timestep check, fixed h=3e-5

| Integrator | DT [s] | max absolute phase error | max transverse matrix error | Classification |
|---|---:|---:|---:|---|
| BORIS | 2e-11 | 1.912e-05 | 2.562e-04 | UNSTABLE |
| BORIS | 1e-11 | 4.842e-06 | 6.529e-05 | STABLE |
| BORIS | 5e-12 | 1.219e-06 | 1.644e-05 | STABLE |
| DOP853 | 2e-11 | 3.637e-10 | 1.097e-08 | STABLE |
| DOP853 | 1e-11 | 4.229e-10 | 2.972e-08 | STABLE |
| DOP853 | 5e-12 | 6.263e-10 | 6.312e-09 | STABLE |
| RK4 | 2e-11 | 2.939e-10 | 4.270e-09 | STABLE |
| RK4 | 1e-11 | 4.198e-10 | 4.651e-09 | STABLE |
| RK4 | 5e-12 | 4.555e-10 | 4.687e-09 | STABLE |

### Perturbation check, fixed DT=5e-12 s

| Integrator | h | max absolute phase error | max modulus error | Classification |
|---|---:|---:|---:|---|
| BORIS | 1e-05 | 1.219e-06 | 7.367e-09 | STABLE |
| BORIS | 3e-05 | 1.219e-06 | 2.604e-09 | STABLE |
| BORIS | 0.0001 | 1.217e-06 | 1.684e-09 | STABLE |
| DOP853 | 1e-05 | 1.619e-09 | 1.931e-08 | NON_UNIT_CIRCLE |
| DOP853 | 3e-05 | 6.263e-10 | 5.126e-10 | STABLE |
| DOP853 | 0.0001 | 4.644e-09 | 6.141e-09 | STABLE |
| RK4 | 1e-05 | 1.463e-10 | 6.784e-09 | STABLE |
| RK4 | 3e-05 | 4.555e-10 | 9.969e-10 | STABLE |
| RK4 | 0.0001 | 4.266e-09 | 2.578e-09 | STABLE |

### Matched old/new exit search

DOP853, RK4 and Boris each launch their COF-verified zero orbit at the same DT=5e-12 s and h=3e-5, Richardson level 0.
The direct COF comparison uses the corresponding integrator and actual beam momentum for the slope conversion.

| Integrator | max old/new 6D difference | old/new exit trials | old/new exit seconds | old/new wall seconds |
|---|---:|---:|---:|---:|
| BORIS | 1.651e-05 | 14400 / 2341 | 1.369 / 0.085 | 2.64 / 1.38 |
| RK4 | 2.429e-08 | 14400 / 1000 | 7.170 / 0.433 | 11.25 / 4.53 |
| DOP853 | 6.890e-09 | 14400 / 1100 | 20.285 / 1.388 | 30.97 / 11.91 |

Boris remains sensitive to exit localization: the old/new 6D difference is 1.65e-5. The new transverse map agrees more closely with the direct Boris COF Jacobian (3.13e-7 versus 1.10e-5), but the old map has a smaller transverse analytic error (7.10e-6 versus 1.68e-5). This does not establish an accuracy winner; it motivates isolating boundary/integration subdivision sensitivity before accepting the optimization for Boris.
Boris exhibits approximately second-order phase-error reduction over the three timesteps. Its coarse-step UNSTABLE label corresponds to only about 1.9e-8 modulus excess; it is a strict numerical diagnostic, not evidence that the analytic lattice is unstable.
RK4 and DOP853 are close to the finite-difference/roundoff floor for this case; these data do not rank their asymptotic orders. Larger h adds differentiation error, and smaller h may add noise. No tolerance was relaxed.
Matrix errors refer to slope coordinates, with position scaled by 1 m. Error ranges are observations against the analytic ideal-lattice reference, not rigorous uncertainty bounds. Timings are individual sequential runs. These checks do not validate the ISIS lattice or its coarse Boris step.

Reproduce grids with run_benchmark.py --integrator METHOD --threads 4 --dt DT --h H --name FRESH_NAME; map comparisons with trusted_map_check.py --integrator METHOD --prefix FRESH_PREFIX (the named COF baseline must exist). summarize_integrators.py regenerates this report and the CSV tables from retained runs.
