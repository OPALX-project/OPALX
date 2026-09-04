# Linear transfer-map validation inputs

These cases exercise `OPTION, ENABLELINEARTRANSFERMAPS = TRUE` against systems with closed-form
first-order maps. They use one on-axis 250 MeV electron, no collective field, hard-edge elements,
and explicit `ELEMEDGE` placement.

- `drift/map-2-drift.in`: one 1 m drift; validates the complete 6x6 map.
- `quadrupole/map-2-quadrupole.in`: one 0.4 m quadrupole with `K1=1 m^-2`; validates the complete
  6x6 map.
- `fodo/map-2-fodo.in`: a symmetric thick-lens FODO cell; validates all entries of the product of
  five exact maps.
- `dba/map-2-dba.in`: two 30 degree sector bends around a thick quadrupole. Its strength is
  the analytic solution for zero exit dispersion in the hard-edge model; all entries of the
  resulting 6x6 map are checked.

For a drift of length L, the nontrivial entries in the coordinates printed by OPALX are

    R12 = R34 = L,    R56 = L/gamma^2.

In the OPALX normal-gradient convention, positive `K1` focuses vertically and defocuses
horizontally. For a horizontally focusing quadrupole, with u=sqrt(-K1), the horizontal block is

    [ cos(uL)       sin(uL)/u ]
    [ -u sin(uL)    cos(uL)   ]

and the vertical block replaces trigonometric functions by their hyperbolic counterparts. Positive
`K1` exchanges the two planes. The FODO reference is the ordered product of these drift and
quadrupole matrices.

The DBA checker uses the complete hard-edge sector-bend matrix. Its horizontal-dispersion block is

    [ cos(theta)       rho sin(theta)   rho(1-cos(theta)) ]
    [ -sin(theta)/rho  cos(theta)       sin(theta)        ]
    [ 0                 0                1                 ]

on `(x,x',delta)`. The remaining nontrivial entries are `R34=rho*theta`, `R51=-sin(theta)`,
`R52=-rho*(1-cos(theta))`, and
`R56=rho*(sin(theta)-theta)+rho*theta/gamma^2`. The bend case retains a slightly looser
interpretation because the transported reference frame and mechanical-momentum coordinates may
differ from canonical textbook conventions.

The reference matrices are continuous-s solutions. The inputs use `DT=2 ps` for the straight
cases and `DT=0.2 ps` for the DBA. Map-enabled reference and shadow-ray integration now subdivide
steps at field-support transitions; no exact element maps are used in OPALX. The tightened
full-matrix tolerances are `1e-8` (drift), `1e-6` (quadrupole/FODO), and `1e-5` (DBA), covering
the measured discretization and finite-difference errors without changing the analytic model.

Each directory is self-contained apart from the shared reference-particle distribution in its
parent. For example, run one case with

    cd drift
    ../../../omp-build/src/opalx --info 2 map-2-drift.in > map-2-drift.out 2>&1

From `sandbox/map-2`, run all cases with `--info 2` and write the complete output to a matching
`.out` file inside each case directory with

    ~/.venv-h6/bin/python check_maps.py ../../omp-build/src/opalx

The script also writes the calculated closed-form matrix to `analytic-map.txt` in every case
directory. It exits nonzero if any entry exceeds its documented numerical tolerance.
`--input-root <directory>` uses an alternate self-contained case tree while still writing the
`.out` and analytic files here, allowing validation without changing an edited input deck.

## DBA time-step convergence

The DBA input deliberately uses `MAXSTEPS=1`: the transfer map is constructed by OrbitThreader
before the single production-tracking step. From `sandbox/map-2/dba`, sweep the four decade values
from `DT=1e-10 s` through `DT=1e-13 s` without modifying the input file:

    ~/.venv-h6/bin/python convergence_dt.py ../../../omp-build/src/opalx

The script saves each complete `--info 2` output as `map-2-dba-dt-<dt>.out`, all 36 matrix
entries as `map-2-dba-dt-<dt>-matrix.txt`, the numerical table as `convergence-dt.csv`, and a
log-log summary as `convergence-dt.png`. `--input <file>` and `--output-dir <directory>` allow
isolated comparisons. The input template is never modified.

Results with boundary subdivision and unchanged finite-difference amplitudes (`1e-3`) are:

| DT [s] | Full 6x6 error | R16 [m] | R26 | abs(det(M)-1) | Canonical-J error | Difference from finest map |
|---:|---:|---:|---:|---:|---:|---:|
| 1e-10 | invalid: too long for QACH | -- | -- | -- | -- | -- |
| 1e-11 | 4.795159e-5 | -1.886340e-5 | -5.631439e-6 | 1.765608e-7 | 2.434575e-7 | 5.334683e-5 |
| 1e-12 | 4.864807e-6 | -1.594085e-6 | -4.164608e-7 | 1.764936e-7 | 2.295652e-7 | 5.304300e-7 |
| 1e-13 | 5.395237e-6 | -1.416116e-6 | -3.617474e-7 | 1.765487e-7 | 2.296295e-7 | 0 |

Compared with the finest numerical map, reducing DT from `1e-11` to `1e-12` reduces the error
by 100.6 (observed order 2.00). This removes the previously observed first-order boundary error.
The comparison with the exact *linear* map reaches a finite-difference floor, so it need not
decrease monotonically at the finest steps. R16 and R26 approach small nonzero finite-difference
residuals; the analytic values remain zero. Matrix stdout now retains twelve decimal places
in scientific notation so printing does not mask the time-discretization trend.

The study was repeated with two MPI ranks; all printed matrix entries and diagnostic residuals
matched the single-rank values. For the local Open MPI launcher, run with

    ~/.venv-h6/bin/python convergence_dt.py ../../../omp-build/src/opalx \
      --ranks 2 --mpi-args '--map-by slot:OVERSUBSCRIBE --bind-to none' \
      --output-dir /tmp/map-2-dba-mpi2

Boundary subdivision uses each element's field-support query and also limits steps relative
to the shortest support extent. It does not estimate field-table interpolation error or
guarantee resolution of grazing crossings or arbitrarily small transverse features. The
production-particle integrator and the existing input time-step safety check are unchanged.
