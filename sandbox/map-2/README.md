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

The reference matrices are continuous-s solutions. The OPALX maps use time integration with
`DT=2 ps`; endpoint and integration errors therefore scale with approximately `c*DT = 0.6 mm`.
The checker tolerances are numerical convergence tolerances, not changes to the analytic maps.

Each directory is self-contained apart from the shared reference-particle distribution in its
parent. For example, run one case with

    cd drift
    ../../../omp-build/src/opalx --info 2 map-2-drift.in > map-2-drift.out 2>&1

From `sandbox/map-2`, run all cases with `--info 2` and write the complete output to a matching
`.out` file inside each case directory with

    ~/.venv-h6/bin/python check_maps.py ../../omp-build/src/opalx

The script also writes the calculated closed-form matrix to `analytic-map.txt` in every case
directory. It exits nonzero if any entry exceeds its documented numerical tolerance.

## DBA time-step convergence

The DBA input deliberately uses `MAXSTEPS=1`: the transfer map is constructed by OrbitThreader
before the single production-tracking step. From `sandbox/map-2/dba`, sweep the four decade values
from `DT=1e-10 s` through `DT=1e-13 s` without modifying the input file:

    ~/.venv-h6/bin/python convergence_dt.py ../../../omp-build/src/opalx

The script saves each complete `--info 2` output as `map-2-dba-dt-<dt>.out`, the numerical table
as `convergence-dt.csv`, and a log-log summary as `convergence-dt.png`.

The current results are:

| DT [s] | Status | Full 6x6 error | max(abs(R16), abs(R26)) | abs(det(M)-1) | Canonical-J error |
|---:|:---|---:|---:|---:|---:|
| 1e-10 | invalid: too long for QACH | -- | -- | -- | -- |
| 1e-11 | OK | 1.143394e-3 | 1.143394e-3 | 3.271809e-7 | 1.113492e-3 |
| 1e-12 | OK | 1.036400e-4 | 2.950285e-5 | 1.679018e-7 | 8.112797e-5 |
| 1e-13 | OK | 8.104655e-6 | 8.104655e-6 | 1.764114e-7 | 3.112249e-6 |

The observed orders of the full-matrix error are 1.04 and 1.11. Thus this range shows roughly
first-order convergence with `DT`. The determinant residual has already reached an approximately
`2e-7` floor, while the full matrix and canonical-form residual continue to improve.
