# Part 1 completion audit {#tps_part1_completion_audit}

Current implementation inventory, 2026-11-24, for the issue-692 branch in
opalx-tps. This replaces the earlier proposal. **The local Part 1 implementation
scope is complete with the explicitly retained numerical backends below. Final
regression acceptance remains pending the push/PR and CI results.** No CI result
or external machine convergence claim is implied by local unit-test success.

## Implemented linear migration

| Area | Current implementation |
| --- | --- |
| Segment construction/composition | Existing 6D finite-difference rays and Richardson refinement; first-order TPS composition in unique reference-path order. |
| Map storage | Immutable TPS segments and combined polynomial in a complete beamline-owned TransferMapPass, including ownerless field tails. Matrices are derived Jacobians. |
| Element ownership | Shared immutable segment handles with separate per-owner attachment ordinals. Owner queries pin records; the complete pass is independent of element ownership. |
| Identity and invalidation | Process-local generation/segment IDs; tested clone, clear, recalculation, invalid-publication, swap, merge and geometry-addition semantics. |
| Runtime teardown | Explicit release before IPPL finalization plus a finalize-hook fallback; retained handles remain safe after payload revocation. |
| 6D input solve | TPS/Kokkos pivoted LU on D_in and transpose solve of D_in^T M^T=D_out^T. No handwritten builder elimination or inverse-product map calculation. |
| Construction diagnostics | TPS/Kokkos determinant and canonical-form products; singular determinant remains a nonfatal diagnostic. |
| COF map/eigenanalysis | Checked first-order 4D TPS snapshot with reference/residual; physical Jacobian passed to the retained LAPACK analysis. |

The coordinate contracts remain distinct: 6D transport uses mechanical slopes
(x,xprime,y,yprime,zeta,delta); the fixed-energy 4D COF section uses mechanical
(x,px,y,py). No implicit conversion or canonical-coordinate claim is introduced.
Compensated reference metadata remains separate from the rounded TPS chart.
See [persistent storage](persistent-storage.md),
[6D numerical contracts](linear-map-algebra.md) and the Doxygen mathematical pages.

## Accepted numerical backends

Accepted policy, 2026-11-24: retain the current COF Newton QR for Part 1.
Its replacement is outside the completed Part 1 implementation scope.
For closure residual r=T(u)-u and M=dT/du, ClosedOrbitSolver solves

\f[ S^{-1}(M-I)S v=-S^{-1}r,\qquad \Delta u=S v. \f]

The existing column-pivoted Householder QR, recomputed trailing column norms,
relative rankTolerance (default 1e-12), reported pivot ratio, coordinate scaling,
damping and failure classification remain unchanged. The pivot ratio is not a
condition number. TPS's 6D LU floor/condition policy is not a replacement for this
rank-revealing Newton policy. A future QR backend change requires separate tests
of rank thresholds, pivoting, scaling, damping and rejected/accepted corrections.
Existing ClosedOrbitSolver tests cover coupled/scaled systems, column permutation,
rank rejection, nonlinear/lost-trial damping and failure statuses.

LAPACK spectrum/SVD, full-tune trajectory winding and the independent spectral
sampling path are also intentionally retained. They perform different tasks from
TPS polynomial storage/composition. Dense arrays at construction, Newton, LAPACK
and reporting boundaries are legitimate buffers, not duplicate persistent maps.

The 6D solve retains strict pivot rejection below 1e-18 and condition rejection
above 1e12, using the original unscaled infinity-norm measure. LU rounding near
thresholds may differ from Gauss–Jordan. One formerly exact dense-vs-TPS product
assertion is now checked against independent long-double dot products with a
calculated roundoff bound after a one-ULP difference was observed; physical and
finite-difference tolerances were not widened. Details are in the numerical report.

## Final unused-code cleanup

Repository source/test/build references were checked before removal:

- Removed Algorithms/LinearTransferMap.h, an unused compatibility forwarding header,
  and its Algorithms install-list entry. Include Structure/LinearTransferMap.h.
- Removed the unused global vector_t alias and prod_matrix_vector template from
  Algorithms/Matrix.h. No in-repository callers remain to migrate.
- Removed the commented getUnit6x6 function from Distribution.cpp and the unused
  direct Matrix.h include in Distribution.h.

The public compatibility header/helper names are intentionally retired; old API
compatibility is not a requirement of this migration. Active 3D coordinate-frame
and device geometry operations remain. The matrix6x6_t construction/reporting
buffer and test-reference products remain in use. The standalone historical TPS
baseline is untouched. No tests are added merely to mirror these deletions;
compilation and the existing numerical suite validate remaining consumers.

## Validation and CI acceptance

After cleanup, the OpenMP Release build and all 176 local CTest cases pass,
including 38 TPS-labelled cases and one-/two-rank adapters/storage regressions.
Coverage includes drift/fringe/Richardson, the tabulated solenoid field-map fixture,
solver status/damping and spectrum/full-tune tests. The full test run took 65.45s.
Focused Doxygen builds with zero warnings; the generated HTML contains the retained
QR policy and current inventory. Local Markdown links, removed-symbol/include
references and git diff --check are clean. No numerical tolerances changed during
cleanup. Build environment and log paths are recorded in HANDOFF-TPS.md.

The existing workflows are unchanged:

- CPU Serial builds/runs unit tests on non-draft pull requests targeting `master`.
- CPU OpenMP additionally requires the `smoke-compile-checks` label.
- CPU Serial Smoke Tests also requires that label and runs the configured
  regression-tests-x cases: Drift-1-fromfile, Drift-4-multi-emit-open and
  FodoCell-fromfile.

A branch push alone does not start these PR-only workflows without a qualifying
PR. CI acceptance must be recorded after the configured checks actually finish;
local success does not substitute for those results. The opt-in PSI COF machine
benchmark was not run locally. No new full-machine convergence/performance or
GPU map-tracking validation is claimed by this cleanup or by the smoke inventory.

## Outside Part 1

Nonlinear map construction/analysis, canonical normal forms, device batch tracking,
space-charge split operators and performance benchmarking remain later work.
Host TPS algebra remains serialized. Fixed-run global-order enforcement and the
previously deferred TPSDependencies.cmake consolidation remain follow-ups, not
Part 1 closure requirements. Direct field/geometry edits still require caller
map invalidation/recalculation as documented in the storage contract.
