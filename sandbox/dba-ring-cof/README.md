# Stable non-achromatic DBA-ring COF benchmark

Separate from the historical achromatic map-2 and unmatched finite-fringe
map-3 fixtures. Current status: COF agrees with the analytic hard-edge matrix
and principal fractional tunes after correcting opposite-arm aperture selection.
Historical failure evidence is retained below.

## Reproduce the analytic reference

Run with `/Users/adelmann/.venv-h6/bin/python`:

```
python analytic_scan.py
python -m unittest discover -s . -p 'test_*.py'
python build_input.py
python run_benchmark.py --name my-validation --h 3e-5
```

The single-family scan samples K1 from -15 to +15 m^-2 in steps of 0.005.
822 sampled settings are stable in both planes. The selected K1=+2.84 m^-2
keeps the original geometry, with no second family needed. The working point
minimises squared distance to principal ring phases (0.23,0.31), subject to
cell half-trace magnitudes below 0.95. This is not a chromaticity optimisation.

Predictions: circumference 25.766370614359175 m; principal fractional phases
0.2298924308671466 and 0.31046360297978914; periodic horizontal dispersion
9.337460094228828 m at the entrance section. Conjugate phases are equally
represented by the spectrum; no integer tune assignment is claimed.

`analytic-results/` contains the complete scan CSV, summary, cell and ring
6D slope matrices, and 4D mechanical-momentum matrix. At the on-axis design
orbit the conversion is T M T^-1 with T=diag(1,p0,1,p0), p0=p/(m_e c).
The inherited map-2 normalized momentum differs slightly from the current
BEAM PC-to-mass conversion (490.23677597553325 versus 490.23677663749817).
Use the actual run momentum for the final coordinate conversion comparison.

`opalx/cof.in` is the new hard-edge COF input. `opalx/run.log` preserves the
first single-rank DOP853 run at DT=1e-11 s. Geometry closure passes, but the
run fails with `Aperture/material loss in C4_D2`, evaluations=1. This historical
failure is diagnosed below; no apertures were enlarged to hide it.
Use a fresh result directory/prefix for subsequent runs.

Diagnosis: `diagnose_aperture.py` reproduces a false positive geometrically.
At the exact design midpoint of C1_D1, local coordinates are (0,0,0.5) m
in C1_D1, but (-7.81051177665,0,0.5) m in opposite drift C4_D2.
Both satisfy 0 <= z < L. CofCmd checks every element's native loss predicate;
ElementBase then rejects the point outside C4_D2's aperture despite the ray
belonging to C1_D1. This proves the current all-element aperture guard can
reject the ideal orbit. It does not establish the absence of other tracking
issues. No production fix or aperture change was made during diagnosis.

## Tracked comparison after the fix

COF now selects the nearest finite conventional design centreline independently
of aperture size. Circular SBEND arcs and straight/RBEND axes are used, including
endpoints; roundoff-level ties are all checked. Native predicates still detect
real aperture losses. Cyclotron domains and field selection are unchanged.
This local-orbit assignment is not a model for intersecting pipe solids.

All physics runs use one MPI rank and one OpenMP thread, DOP853, K1=2.84 m^-2.
`results/nearest-body-dt1e-11-h3e-5/` is the selected benchmark, with DT=1e-11 s,
position FD=3e-5 m and normalized-momentum FD=p0*3e-5 (slope-equivalent).
Newton residual and eigenvalue tolerances are unchanged. It reports STABLE:

| Plane | Analytic principal phase | Tracked minus analytic |
|---|---:|---:|
| Horizontal | 0.2298924308671466 | +4.229470e-10 |
| Vertical | 0.3104636029797891 | +1.950489e-10 |

Maximum matrix-entry error is 2.972264e-8 after conversion to slope coordinates
(position entries scaled by 1 m for a dimensionless comparison). The coordinate
conversion uses the actual BEAM normalized momentum from stdout. Maximum
eigenvalue-modulus error is 1.89e-9, within the unchanged 1e-8 band.
The on-axis seed remains (0,0,0,0), with return residual x=4.70e-15 m and
px=-5.78e-12 in p/(mc); relative energy drift is -2.56e-15.

`results/nearest-body-displaced/` starts at x=y=0.1 mm and px=py=0.001.
It converges to x=-4.21e-13 m, px=4.72e-11, with transverse coordinates
consistent with the design orbit and STABLE eigenanalysis.

The saved scans include DT=1e-11 and 5e-12 s and FD amplitudes 3e-6, 1e-5,
3e-5 and 1e-4. At the smaller FD values the horizontal modulus is noisier and
the strict classification is NON_UNIT_CIRCLE (roughly 2e-8 to 6e-8 modulus
error); at 3e-5 and 1e-4 it is STABLE. Thus tiny FD steps are not intrinsically
more accurate; no classification threshold was relaxed to select this result.
Each run preserves input, stdout, binary/input hashes, matrix difference,
tune CSV and comparison JSON. The runner reports eigenphases even in unstable
cases; only STABLE classifications justify calling them stable mode tunes.

The timestep-halved check `results/nearest-body-dt5e-12-h3e-5/` is also STABLE:
max matrix error 6.312143e-9; horizontal tune 0.2298924314934086, vertical tune
0.3104636031748715. Maximum modulus error is 5.13e-10. This confirms agreement
at another timestep without requiring monotonically decreasing phase errors
at the finite-difference noise floor.

Three analytic tests and seven CTest suites pass, including dedicated two-rank
COF consistency/failure checks and new opposite-arm/local-loss tests. The source
and user manual document the spatial aperture assignment and its limitation.

Next: off-momentum dispersion/chromaticity, spectral cross-checks and finite
fringes remain separate future stages, not claims of this hard-edge benchmark.


## Trusted baseline: current four-thread study (2026-09-07)

The new `results/trusted-20260907-*` runs retain the same analytic hard-edge
lattice and use one MPI rank with four OpenMP threads. This is a trusted reference
for this idealized lattice, not a validation of ISIS or finite-fringe optics.
The analytic formulas reproduce the retained matrix, and all three analytic
sanity tests pass. No physics source or numerical acceptance tolerance changed.

Selected COF baseline: DOP853, DT=5e-12 s, h=3e-5 m/slope-equivalent, mechanical
momentum perturbations p0*h using the actual BEAM p0. The zero design orbit is
verified; maximum closure residuals are 2.16e-13 m and 3.28e-11 p/(mc).
Principal fractional tunes are Qx=0.2298924314934086 and Qy=0.3104636031748715,
with analytic errors +6.26e-10 and +1.95e-10. Conjugate phases are 1-Q; integer
and oriented tune branches are not determined. Maximum modulus error is 5.13e-10.

| DT [s] | h [m/slope-equivalent] | max transverse matrix error | max tune error | classification |
|---|---|---|---|---|
| 2e-11 | 3e-5 | 1.10e-8 | 3.64e-10 | STABLE |
| 1e-11 | 3e-5 | 2.97e-8 | 4.23e-10 | STABLE |
| 5e-12 | 3e-5 | 6.31e-9 | 6.26e-10 | STABLE |
| 5e-12 | 1e-5 | 1.49e-8 | 1.62e-9 | NON_UNIT_CIRCLE |
| 5e-12 | 1e-4 | 5.28e-8 | 4.64e-9 | STABLE |

These observed errors support a tested envelope below 5e-9 in principal tune and
6e-8 in the scaled transverse matrix entries. They are not rigorous uncertainty
bounds. Nonmonotonic timestep errors indicate no resolved asymptotic integration
order here. The smaller perturbation fails the unchanged 1e-8 modulus band;
we do not hide this failure or assume that smaller perturbations are better.
Matrix comparisons use slope coordinates, with position scaled by 1 m.

### Old versus new exit localization

`trusted_map_check.py` generates identical TRACK inputs from the stable COF
lattice and launches the verified zero orbit with the actual p0. Both use
DOP853, DT=5e-12 s, h=3e-5 in all six map coordinates and Richardson level 0.
The comparison uses the saved pre-exit-search executable and current executable;
SHA256 hashes are retained in each result directory. No rebuild is required to
reproduce the comparison while those binaries remain available.

| Measurement | Old bisection | Safeguarded secant |
|---|---:|---:|
| max 6D analytic matrix error | 6.24e-8 | 6.93e-8 |
| max transverse analytic matrix error | 6.70e-9 | 8.55e-9 |
| max transverse difference from direct COF | 1.30e-8 | 1.49e-8 |
| exit integrations | 14400 | 1100 |
| exit time [s] | 20.28 | 1.39 |
| process wall time [s] | 30.97 | 11.91 |

The two 6D matrices differ by at most 6.89e-9. All 9084 design-path rows and
printed reference return diagnostics agree exactly. Both transverse eigenpairs
lie within 6.4e-10 of the unit circle, and principal phase errors are below
1.9e-10. Timings are individual sequential runs, not statistical speedup estimates.
The secant result is comparable in accuracy, not systematically more accurate.

The direct COF derivative and the product of segment derivatives are different
finite-difference approximations. COF probes mechanical momentum, whereas map
rays probe slopes; the linear conversion is exact at the zero fixed point, but
finite probe amplitudes differ at higher order. Agreement at the observed scale
is the relevant check, not bitwise equality between those methods. The analytic
6D reference retains its historical gamma in the longitudinal term; its p0
conversion difference is negligible for these error levels and does not enter
the transverse tune reference.

Reproduce the COF grid with `run_benchmark.py --threads 4 --dt ... --h ... --name ...`
and fresh directory names. `summarize_trusted_study.py` recomputes the analytic
reference check and collects the named study results without rerunning tracking.
`trusted_map_check.py` runs both map binaries sequentially and verifies input,
particle, trajectory and reference-return agreement. The initial map input was
rejected for combining FROMFILE with BEAM PC; that rejected attempt is retained
in `trusted-20260907-map-old-input-rejected`. The corrected input takes its
particle momentum only from the distribution file, as required by TRACK.

Raw summaries: [COF CSV](results/trusted-20260907-cof-summary.csv),
[COF JSON](results/trusted-20260907-cof-summary.json),
[map comparison](results/trusted-20260907-map-comparison.json).

The additional `trusted-20260907-displaced` check starts at x=y=0.1 mm and
px=py=0.001 p/(mc), using the selected timestep and perturbation. It converges
back to x=8.29e-14 m, px=3.98e-13 p/(mc), y=-2.37e-18 m,
py=9.96e-16 p/(mc), and reports STABLE. Its largest transverse matrix error
is 4.42e-9. This checks Newton recovery from a displaced seed as well as closure
of the on-axis design seed.

The subsequent [Boris/RK4/DOP853 comparison](INTEGRATOR_COMPARISON.md) contains the full timestep and perturbation grids, matched-setting summary, and old/new exit-search comparison. Boris retains material search sensitivity on this trusted lattice; see the qualification in that report.


## Commit scope and retained artifacts

The implementation being committed retains the original exit-plane bisection
for Boris (including LF2). RK4 and DOP853 use the safeguarded secant search.
The `trusted-boris` study records the pre-fallback implementation and demonstrates
why this restriction is necessary; it is not the behavior of the final Boris code.
`results/commit-boris-fallback-map-comparison.json` verifies that the final Boris
matrix is identical at printed precision to the saved original-bisection result.
Both runs use 14400 exit iterations and identical reference trajectories.

To reproduce the historical unrestricted comparison, supply its executable with
`trusted_map_check.py --new-binary PATH`; `--old-binary PATH` selects the baseline.
The unrestricted executable is not bundled; supply a build of the pre-fallback
implementation to repeat that historical experiment.
Executable hashes in the retained provenance identify the actual measured binaries.
Executables are not committed. Reproduction requires compatible built binaries;
benchmark scripts default to the local `omp-build/src/opalx` executable.

Versioned artifacts include scripts and their source dependencies, the input,
analytic reference, compact trusted-study results and fallback verification.
Full trajectories, stdout logs, executable files and earlier historical run
folders are local artifacts. References to those historical runs above record
provenance; they are not required by the trusted-study summary scripts.
Python analysis requires numpy and pandas; tests use unittest.
