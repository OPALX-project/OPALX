# map-3: DBA with finite bend fringes

First stage: **analytic fields on the nominal design line**, with `HGAP=0.1 m`
and `FINT=0.1`. Bend supports are separated from the explicitly field-free
drifts. This is not yet a matched achromat, an analytic transfer matrix, or an
OPALX tracking result. A subsequent three-case OPALX reproduction is recorded
below; no new convergence sweep has been started.

The authoritative equations, units and limitations are in the physics manual:
[map-3 fringe DBA](https://github.com/OPALX-project/opalx-manual/blob/main/physics/map-3-fringe-dba.qmd).

From the OPALX repository root, reproduce the figures and run the checks:

```sh
MPLCONFIGDIR=/tmp/opalx-map3-mpl ~/.venv-h6/bin/python sandbox/map-3/dba/analytic_fields.py
MPLCONFIGDIR=/tmp/opalx-map3-mpl ~/.venv-h6/bin/python -m unittest discover -s sandbox/map-3/dba -p 'test_*.py' -v
ctest --test-dir omp-build -R '^TestBendFieldModel$' --output-on-failure
```

`dba/parameters.json` is the single study-parameter source. Generated files in
`dba/analytic/` include centreline/off-axis PNG and SVG figures, all three local
field components in `fields-vs-s.csv`, the labelled layout in `layout.csv`, and
integrals and source hashes in `summary.json`. The layout is calculated, not
read from an `ElementPositions.sdds` produced by tracking.

The quadrupole uses the old map-2 K1 only as an **unmatched seed**. Matching
remains future work; an independent reference comparison is now available below.

## OPALX reproduction: one archived-best case per integrator

`dba/map-3-dba.in` reproduces the analytic model without changing its dipole
normalization or quadrupole strength. `MAXSTEPS=1` applies to production tracking;
OrbitThreader computes the full map over the requested 8.294395102393196 m.
The production `.stat` file is not a full field history. The field comparison
uses each map integrator's own `data/map-3-dba_DesignPath.dat` instead.

The original runner selected the minimum full-map error among **all saved valid
single-rank observations** in the stopped map-2 study. The default runner now
replays those published settings and audits the committed map-3 artifacts;
it needs no private campaign directories. Optional `--reselect-map2` repeats
the original selection when the complete map-2 archive is available.
These are observed minima, not guaranteed portable optima;
in particular, the RK4 point may benefit from error cancellation.

| Map integrator | DT [s] | Starting epsilon | Richardson levels |
|:--|--:|--:|--:|
| BORIS | 1e-13 | 1e-5 | 1 |
| RK4 | 3e-11 | 3e-5 | 0 |
| DOP853 | 1e-11 | 0.003 | 1 |

Each epsilon value is used in all six coordinates (metres for positions,
dimensionless for slopes and relative momentum). The three cases ran
**sequentially**, with one MPI rank and one OpenMP thread, and all completed.
Every case retains `map-3-dba.in`, `map-3-dba.out` (`--info 2`), `timing.dat`,
`matrix.txt`, the design path and recorded element supports.
The original runner, field validator and comparison script are retained in
`opalx-best-settings/source/` to verify their historical source hashes (before
the portable replay/audit changes). Executables are not committed; field validation audits
source and result hashes, with optional `--require-current-build` for a strict
same-machine executable check. New runs always record their own executable hash.

```sh
cmake --build omp-build --target opalx_exe -j 4
# Existing runs are protected: reproducing requires a new output directory.
~/.venv-h6/bin/python sandbox/map-3/dba/run_best_cases.py --output /tmp/opalx-map3-repeat
~/.venv-h6/bin/python sandbox/map-3/dba/validate_opalx_fields.py --output /tmp/opalx-map3-repeat
# Inspect saved results without running OPALX again:
~/.venv-h6/bin/python sandbox/map-3/dba/run_best_cases.py --report-only
~/.venv-h6/bin/python sandbox/map-3/dba/validate_opalx_fields.py
```

Results are in `dba/opalx-best-settings/results.csv`; differences between all 36
matrix entries are summarized in `pairwise-differences.csv`. These differences
are **not errors against an analytic map**. All cases visit all five elements
without recorded support overlaps, and logged D1/D2 fields are exactly zero.
The analytic field evaluated at the logged global positions agrees within the
eight-significant-digit log precision. Boundary samples whose rounded position
is ambiguous at a hard quadrupole face are explicitly counted and excluded.
Detailed checks and timings are in `field-validation-summary.csv` and each
case's `field-validation.json`/`.csv`.

The actual orbit differs from the nominal design line: compare fields at the
same global position, not merely at the same value of `s`. The final logged
trial midpoint can lie beyond ZSTOP before endpoint localization; it is not the
map's final reference plane. Recorded B2 support ends before ZSTOP here, leaving
a short final field-free interval included in the combined map.

The unchanged 1e-6 determinant/canonical-J diagnostic threshold passes for
BORIS and DOP853; RK4 is slightly above it. No tolerances were loosened and no
additional tuning runs were made. Quantitative results and physical caveats
are documented in the physics manual's map-3 reproduction section.

## Independent reference and analytic thin-edge comparison

`variational_reference.py` integrates analytic field sensitivities without
OPALX or shadow-ray finite differences. It retains actual particle momentum,
hard quadrupole-face crossing sensitivity and the reference-normal final plane.
This is a numerical evaluation of the distributed analytic fields, not a
closed-form map. `compare_reference.py` also evaluates the separate CERN
thin-edge analytic approximation with full-gap FINT convention. Differences
from that approximation include model and reference-orbit differences, not
only numerical integration errors. See the manual for the equations and the
identified full-gap/half-gap discrepancy in the current OPALX edge coefficient.

```sh
# Requires a fresh output directory if reproduced again.
~/.venv-h6/bin/python sandbox/map-3/dba/variational_reference.py --output /tmp/opalx-map3-reference
~/.venv-h6/bin/python sandbox/map-3/dba/compare_reference.py --reference /tmp/opalx-map3-reference --output /tmp/opalx-map3-comparison
```

Saved reference refinements are in `dba/variational-reference-validated/`;
all 36 signed differences for each method, a heatmap, the separate CAS matrix,
summary table and provenance are in `dba/map-comparison/`. No new OPALX runs
or production-code changes were made for this comparison. The earlier
`dba/variational-reference/` contains the initial numerical cross-check, retained
only as intermediate analysis; use the validated directory for reported results.
