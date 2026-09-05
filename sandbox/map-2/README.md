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

## DBA perturbation-amplitude convergence

At fixed `DT=1e-13 s` and `MAXSTEPS=1`, the six centered-difference amplitudes were varied
together. The tabulated epsilon is in metres for `(x,y,zeta)` and dimensionless for
`(x',y',delta)`; it is not a common physical unit for all six coordinates. The same hard-edge
DBA and analytic full matrix were used throughout, with analytic `R16=R26=0`.

| Epsilon | R16 [m] | R26 | Full 6x6 error | abs(det(M)-1) | Canonical-J error |
|---:|---:|---:|---:|---:|---:|
| 1e-3 | -1.416116e-6 | -3.617474e-7 | 5.395237e-6 | 1.765487e-7 | 2.296295e-7 |
| 1e-4 | -8.487746e-9 | +1.875595e-8 | 9.610314e-8 | 1.864153e-9 | 6.034381e-8 |
| 1e-5 | -1.081298e-7 | +1.318702e-7 | 4.252702e-7 | 6.645471e-9 | 7.375886e-7 |
| 1e-6 | -7.989188e-6 | -7.976136e-7 | 7.989188e-6 | 1.983597e-8 | 1.747449e-5 |

The full-matrix error is the largest absolute entrywise difference in the printed coordinate
convention (entries have mixed units). `1e-4` is the best **tested** amplitude for this case,
reducing that error by about 56 relative to `1e-3`. Smaller amplitudes worsen the result.
This is consistent with numerical ray errors being amplified by finite differencing; the sweep
does not separate trajectory integration, boundary resolution and roundoff contributions, nor
establish a universal optimal amplitude for other lattices or field maps. The canonical-J
residual remains a coordinate-dependent diagnostic, not a general symplecticity certificate.

Every amplitude was run with one and two MPI ranks. All 36 printed entries and both reported
diagnostics matched exactly between rank counts. Neither the differentiation order nor the
integrator was changed. The original `1e-3` source defaults and executable were restored afterward.

`dba/convergence_perturbations.py` now uses the runtime map options, so no temporary source
edits are needed. It rebuilds the configured target and records the experiment. Run from the
repository root, changing `--epsilon` for each amplitude:

    ~/.venv-h6/bin/python sandbox/map-2/dba/convergence_perturbations.py \
      --input sandbox/map-2/dba/map-2-dba.in --epsilon 1e-3

Repeat for `1e-3`, `1e-4`, `1e-5`, and `1e-6`. Do not run competing builds/source edits during
the experiment. To check selectable Richardson refinement, add `--richardson-levels 1`
(or 2, 3, 4) and use a separate `--output-dir` for each level. For example:

    ~/.venv-h6/bin/python sandbox/map-2/dba/convergence_perturbations.py \
      --epsilon 1e-3 --richardson-levels 1 --output-dir sandbox/map-2/dba/richardson/levels-1

Repeat at `5e-4` and `2.5e-4` to test moderate starting amplitudes without shrinking the
perturbations directly into the small-amplitude noise regime. The input options and their
mathematics are documented in the OPALX manual (Runtime Options / Numerical linear transfer maps).
The input template is never modified; temporary run copies use `DT=1e-13 s`. These recorded
runs used the committed input template (different diagnostic output frequencies from the working
deck, identical physics). Each `result.json` includes that complete template and source/binary
hashes for provenance.

Results live in `dba/perturbations/eps-<epsilon>/ranks-<ranks>/`: full `--info 2` stdout in
`map-2-dba-dt-1e-13.out`, all 36 entries in `map-2-dba-dt-1e-13-matrix.txt`, and diagnostics
in `result.json`. The aggregate [CSV](dba/perturbations/convergence-perturbations.csv) and
[plot](dba/perturbations/convergence-perturbations.png) can be regenerated without building or
running OPALX:

    ~/.venv-h6/bin/python sandbox/map-2/dba/convergence_perturbations.py --summarize-only

### Richardson verification

With the same DBA, `DT=1e-13 s`, and the runtime settings above, the complete matrices give:

| Starting epsilon | Richardson levels | R16 [m] | R26 | Full 6x6 error |
|---:|---:|---:|---:|---:|
| 1e-3 | 0 | -1.416116e-6 | -3.617474e-7 | 5.395237e-6 |
| 5e-4 | 0 | -3.398227e-7 | -8.725519e-8 | 1.351277e-6 |
| 2.5e-4 | 0 | -9.084693e-8 | -3.263752e-8 | 3.402266e-7 |
| 1e-3 | 1 | +1.894120e-8 | +4.242136e-9 | 2.015959e-8 |
| 5e-4 | 1 | -7.855008e-9 | -1.443163e-8 | 2.519152e-8 |
| 2.5e-4 | 1 | -8.741044e-9 | -4.866523e-9 | 1.094175e-7 |
| 1e-3 | 2 | -9.641421e-9 | -1.567655e-8 | 2.616852e-8 |

All seven configurations were repeated on two MPI ranks, with identical printed matrices
and diagnostics. Level zero at `1e-3` reproduces the earlier default matrix exactly.
One refinement at starting epsilon `1e-3` reduces the full-matrix error by about 268,
but further refinement does not improve it monotonically. This study does not demonstrate
`1e-10` map accuracy. The separate relativistic-drift unit tests verify the formal second-,
fourth-, and sixth-order differentiation regimes and all supported tableau levels.

Full `--info 2` output, matrices and provenance are in `dba/richardson/levels-<L>/eps-<epsilon>/`.
Each level also contains its own [level-0 CSV](dba/richardson/levels-0/convergence-perturbations.csv),
[level-1 CSV](dba/richardson/levels-1/convergence-perturbations.csv), or
[level-2 CSV](dba/richardson/levels-2/convergence-perturbations.csv), plus a dispersion plot.

## Joint DBA time-step, amplitude and Richardson study

`dba/convergence_grid.py` extends those one-dimensional studies to a joint grid:

- `DT = {1e-10, 3e-11, 1e-11, 3e-12, 1e-12, 3e-13, 1e-13}` seconds.
- Starting `epsilon = {1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6}`.
- Richardson levels `L = 0, 1, 2, 3, 4`: centered differences and formal perturbation
  orders 4, 6, 8, 10. These are **not** different time-integration orders.

All cases use the same boundary-split Boris/LF2 trajectory integrator, the same hard-edge
DBA, all six starting amplitudes varied together, and `MAXSTEPS=1`. No exact element maps
are used by OPALX; the analytic product in `check_maps.py` is only the external reference.
Neither the production source nor the user's input template is modified by the runner.
The 385 single-rank configurations include the intentionally too-coarse time step; the
existing `OrbitThreader::checkElementLengths` rejection for QACH is retained, not bypassed.
The two largest amplitudes extend the initial 315-case design to better expose the
truncation-dominated Richardson regime; they are study points, not default settings.

For an explicit, dimensionally consistent full-map comparison, define

$$
S=\operatorname{diag}(1\,\mathrm{m},1,1\,\mathrm{m},1,1\,\mathrm{m},1),\qquad
A=S^{-1}MS.
$$

The plotted overall error is the maximum absolute entry difference

$$
E_\infty=\max_{i,j}|A_{ij}-A_{ij}^{\mathrm{exact}}|.
$$

This is dimensionless with the stated scales, and numerically equals the earlier raw SI
entrywise error. It is not a relative error or the matrix infinity norm based on row sums.
The CSV also reports the relative Frobenius error, the worst entry, all 36 signed entries,
`abs(det(M)-1)`, the coordinate-dependent canonical-J diagnostic, and wall times. Timings
come from concurrent independent processes and are not isolated performance benchmarks.

Reproduce from the repository root:

```sh
cmake --build omp-build --target opalx_exe -j8
~/.venv-h6/bin/python sandbox/map-2/dba/convergence_grid.py --workers 4
~/.venv-h6/bin/python sandbox/map-2/dba/validate_convergence_grid.py \
  sandbox/map-2/dba/convergence-grid
```

The runner fixes `OMP_NUM_THREADS=1`, runs independent cases in temporary directories,
and saves every configured deck, complete `--info 2` `.out`, matrix and checksum under
`dba/convergence-grid/L<L>/dt-<dt>/eps-<epsilon>/ranks-<ranks>/`. Re-running the same command
resumes verified results; changed source, binary, distribution or input requires a new
`--output-dir`. The top-level provenance records the executable/source hashes and input,
and `source-changes.patch` captures the production-source diff against the recorded Git revision.
Archived decks retain the template's relative distribution path; use the runner to replay
them in its self-contained temporary layout.

Aggregate data and figures (PNG and editable vector SVG) can be regenerated without OPALX:

```sh
~/.venv-h6/bin/python sandbox/map-2/dba/convergence_grid.py --summarize-only
```

- [All cases and all 36 entries](dba/convergence-grid/results.csv).
- [Full map, R16 and R26 versus DT](dba/convergence-grid/errors-vs-dt.png), at two fixed amplitudes.
- [Full map, R16 and R26 versus epsilon](dba/convergence-grid/errors-vs-epsilon.png), at two fixed time steps.
- [Complete three-dimensional grid](dba/convergence-grid/complete-grid.png), including invalid cases.
- [Signed R16 and R26](dba/convergence-grid/signed-dispersion.png), showing sign changes hidden by log-magnitude plots.
- [Best tested amplitude at each DT](dba/convergence-grid/best-vs-dt.png); selection minimizes the full-map
  error, and both dispersion entries are taken from that same map.
- [Difference from the finest numerical map](dba/convergence-grid/time-differences.png), at fixed
  differentiation settings. The finest map is not an exact reference.
- [Adjacent-point observed orders](dba/convergence-grid/observed-orders.csv); orders near an error
  floor or a sign cancellation are not reliable asymptotic-order estimates.

### Joint-study results

All 385 primary configurations completed: 330 valid maps and 55 expected QACH guard
rejections at `DT=1e-10 s`. The best **tested full-matrix** result for each level is:

| Richardson L | DT [s] | Starting epsilon | Full-map error E_inf | R16 [m] | R26 |
|---:|---:|---:|---:|---:|---:|
| 0 | 3e-13 | 1e-4 | 6.165533e-8 | -1.594706e-8 | -2.065981e-8 |
| 1 | 1e-13 | 3e-3 | 7.478657e-9 | -7.478657e-9 | +7.614869e-10 |
| 2 | 1e-13 | 3e-2 | 1.141348e-8 | -8.072299e-9 | -1.570483e-9 |
| 3 | 1e-13 | 3e-2 | 9.118535e-9 | -5.625266e-9 | -2.407758e-9 |
| 4 | 1e-13 | 1e-1 | 9.361484e-9 | -7.103862e-9 | -2.120369e-9 |

The best configuration is L=1, epsilon=3e-3, DT=1e-13 s. Its relative Frobenius
error is `9.865179e-10`, but its largest scaled entry error remains `7.478657e-9`:
these are different accuracy criteria. Its largest error is R16. A nearly zero
dispersion entry alone does not establish accuracy of the other 35 entries.

The convergence regimes are distinct:

1. **Time discretization:** at L=1 and epsilon=1e-3, differences from the finest
   numerical map have observed DT orders 2.001 and 2.004 over the intervals
   `1e-11 -> 3e-12 -> 1e-12 s`. Richardson refinement in epsilon does not change
   the Boris time-integration order.
2. **Finite-difference truncation:** at DT=1e-13 s, L=0 gives epsilon orders
   2.000 and 2.000 over `1e-2 -> 3e-3 -> 1e-3`. L=1 gives orders 4.015 and
   4.167 over `1e-1 -> 3e-2 -> 1e-2`, before reaching the numerical floor.
   The DBA grid does not resolve clean sixth-, eighth- or tenth-order regimes;
   formal order is not a promise of accuracy at arbitrary amplitudes.
3. **Small-amplitude degradation:** at DT=1e-13 s, L=4 and starting epsilon=1e-6
   give full-map error `2.737e-4`, despite the higher formal order. The finest rays
   then use epsilon/16. The behavior is consistent with numerical ray errors
   amplified by differentiation; this study does not separately identify boundary,
   trajectory roundoff and integration contributions.

More refinement and a smaller DT are therefore not universally better. L=1 offers
the best accuracy among these tested points with only 24 rays per segment; L=4 uses
60. The best results for L>=1 are around `1e-8` in the maximum-entry criterion, and
this study does **not** demonstrate `1e-10` full-map accuracy. Optimal settings are
benchmark-specific; error cancellation can also improve individual grid points.
No runtime defaults were changed on the basis of this study. This is a hard-edge
DBA benchmark, not a test of field-table interpolation or a closed-orbit solver.

Thirteen selected configurations were repeated with two MPI ranks: all levels at
DT=3e-12/epsilon=1e-3; L=0 at DT=3e-13/epsilon=1e-4; L=1..4 at
DT=1e-13/epsilon=3e-3; and L=2..4 at DT=1e-13/epsilon=1e-1. All 36 printed
entries and both reported diagnostics agree exactly with their single-rank peers.
The earlier epsilon=1e-3, DT=1e-13 matrices for L=0,1,2 were also reproduced exactly.
Individual repeats can be made with `--ranks 2 --dt <dt> --epsilon <epsilon> --levels <L>`.

The independent audit verifies case coverage, every archived artifact checksum, every
matrix entry against stdout, both error norms, and MPI agreement. All eight fast
runner tests pass; they cover cache provenance, corruption detection, unexpected
failures, and expected-rejection handling:

```sh
~/.venv-h6/bin/python sandbox/map-2/dba/test_convergence_grid.py
```

The figures were regenerated from the complete dataset and visually checked, including
the equation labels and signed-axis zero region. `analysis-provenance.json` records the
postprocessing-script hashes. The original input deck and all earlier studies are preserved.

### Bookkeeping accuracy comparison

The follow-up study is stored separately in `dba/convergence-grid-bookkeeping`;
`dba/convergence-grid` remains the historical baseline. It uses the **same 385
input configurations**, analytic oracle, coordinate scales, Boris pusher,
Richardson levels and field-support tolerances. The change is limited to physical
half drifts, compensated position/time/path accumulation, speed-integrated path
length and tracked start/stop localization. The derivation and limitations are in
the physics manual's *Numerical linear transfer maps: Position, clock and path
bookkeeping* section; the equations are also in `ExternalFieldRayTracker` Doxygen.

Reproduce without overwriting the before-study:

```sh
~/.venv-h6/bin/python sandbox/map-2/dba/convergence_grid.py \
  --output-dir sandbox/map-2/dba/convergence-grid-bookkeeping
~/.venv-h6/bin/python sandbox/map-2/dba/repeat_grid_mpi.py \
  sandbox/map-2/dba/convergence-grid sandbox/map-2/dba/convergence-grid-bookkeeping
~/.venv-h6/bin/python sandbox/map-2/dba/convergence_grid.py --summarize-only \
  --output-dir sandbox/map-2/dba/convergence-grid-bookkeeping
~/.venv-h6/bin/python sandbox/map-2/dba/validate_convergence_grid.py \
  sandbox/map-2/dba/convergence-grid-bookkeeping
~/.venv-h6/bin/python sandbox/map-2/dba/compare_convergence_grids.py \
  sandbox/map-2/dba/convergence-grid sandbox/map-2/dba/convergence-grid-bookkeeping \
  --output sandbox/map-2/dba/convergence-grid-bookkeeping/comparison
```

The runner refuses to mix different executable/source/input hashes in one output
directory. After future source changes, select a new output directory. Archived
CSV comparison and `--summarize-only` do not require the original executable.
To audit the historical baseline against its earlier matrices, add
`--check-legacy` to its validator invocation; that check is intentionally not
applied to the changed numerical implementation.

The comparison checks identical template/distribution/analytic-reference inputs,
matching acceptance/rejection and every configured input-file hash. All paired
errors, signed entries and gains are saved in
[paired-results.csv](dba/convergence-grid-bookkeeping/comparison/paired-results.csv).
A gain is `error_before / error_after`; greater than one means improvement, not a
statistical significance claim. Each optimized row selects its amplitude/time
step by the **full-map max-entry error**, not by a dispersion entry separately.

- [Matched L=1 DT and epsilon cuts](dba/convergence-grid-bookkeeping/comparison/before-after.png).
- [All five levels versus epsilon at DT=1e-13](dba/convergence-grid-bookkeeping/comparison/all-levels-vs-epsilon.png).
- [Signed R16 and R26 at matched settings](dba/convergence-grid-bookkeeping/comparison/signed-before-after.png).
- [Full-map gain across every valid pair](dba/convergence-grid-bookkeeping/comparison/paired-full-map-gains.png).
- [Best full-map-selected results before/after](dba/convergence-grid-bookkeeping/comparison/best-before-after.csv).

PNG and vector SVG versions are generated from the archived data. Fast runner and
comparison guards can be run together with:

```sh
~/.venv-h6/bin/python -m unittest discover -s sandbox/map-2/dba -p 'test_*grid*.py'
```

#### Completed comparison: robustness improves, peak accuracy does not

All 385 configurations completed (330 valid maps, the same 55 QACH rejections).
The full sweep took 862.9 s with four independent workers, reusing two validated
smoke-run points. The 13 previous MPI validation points were repeated: every
printed matrix entry and both diagnostics match their new single-rank peers
exactly. The independent audits passed for both archived studies; the new source
snapshots also match the recorded hashes.

Best **full-map max-entry errors** over the same grid:

| L | Before | After | After epsilon | Before/after gain |
|---:|---:|---:|---:|---:|
| 0 | 6.165533e-8 | 9.296484e-9 | 3e-6 | 6.63 |
| 1 | 7.478657e-9 | 8.366484e-9 | 1e-5 | 0.894 |
| 2 | 1.141348e-8 | 9.152485e-9 | 1e-4 | 1.247 |
| 3 | 9.118535e-9 | 9.393484e-9 | 3e-5 | 0.971 |
| 4 | 9.361484e-9 | 9.283484e-9 | 3e-4 | 1.008 |

Every after optimum uses DT=1e-13 s; the earlier optimum settings are recorded
above and in the comparison CSV. The new global optimum is about **11.9% worse**
than the previous global optimum. At that new optimum, R16=-5.338758e-9 m and
R26=-2.974798e-10, but R12 limits the full-map error to 8.366484e-9. These small
dispersion entries therefore do not establish 1e-10 accuracy for the full matrix.

The gain is much stronger at matched small amplitudes. At DT=1e-13 s:

| L | Same epsilon | Before full-map error | After full-map error | Gain |
|---:|---:|---:|---:|---:|
| 1 | 3e-3 | 7.478657e-9 | 9.762484e-9 | 0.766 |
| 1 | 1e-3 | 2.015959e-8 | 9.824485e-9 | 2.05 |
| 1 | 1e-5 | 4.235746e-6 | 8.366484e-9 | 506 |
| 1 | 1e-6 | 1.721532e-5 | 2.139090e-8 | 805 |
| 4 | 1e-6 | 2.736751e-4 | 2.937909e-7 | 932 |

Across all 330 pairs, 132 errors decrease, 197 increase and one is unchanged at
stored precision; the median before/after gain is 0.964. Thus this is **not** a
uniform improvement across the whole grid. It strongly reduces small-amplitude
noise and broadens the near-1e-8 plateau, while the time-discretization-dominated
region remains similar. L=1/epsilon=1e-3 differences from the finest numerical map
give observed DT orders 1.997 and 2.008 on 1e-11 -> 3e-12 -> 1e-12 s, consistent
with the unchanged second-order integrator.

No default perturbations, integration method, Richardson formula or existing
test tolerance was changed. This increment is useful bookkeeping infrastructure,
but does not achieve a lower full-map accuracy floor on this benchmark. Higher
integration order remains a next candidate, not a result established by this
study. Momentum accumulation, boundary localization and field interpolation
have not been separately isolated; no field-table or closed-orbit accuracy claim
is made. The combined package was compared here, not an ablation of each change.

Validation: 87/87 CTests, 20/20 focused tests on both one and two MPI ranks, 12/12
Python study/comparison tests, and all four existing analytic optics regressions
pass. Doxygen has no warnings; its formulas and the physics-manual formulas were
rendered and visually checked. All ten new-study/comparison figures were reviewed.
