# 6D linear-map numerical backend {#tps_linear_map_algebra}

The builder now delegates its input-difference right solve and map diagnostics to
`Algorithms/TpsLinearMapAlgebra`. Finite-difference rays, field selection, frames,
exit localization and Richardson extrapolation remain unchanged. Dense matrices
are local construction/reporting buffers; persistent TPS storage is unchanged.
The working implementation removes the old Gauss–Jordan inverse, explicit
inverse-product map calculation, hand-written determinant elimination and legacy
6D canonical-form matrix products from the builder. The historical baseline is
untouched. Active 3D geometry helpers in Algorithms/Matrix.h remain in use.

For M D_in=D_out, factor D_in once and solve D_in^T M^T=D_out^T with TPS LU's new
transpose solve. Factoring D_in^T instead would change the pivots tested against
the absolute floor. Pivot magnitudes strictly below 1e-18 are rejected, as are
unscaled original infinity-norm condition measures strictly above 1e12. Equality
is accepted.
Identity solves still calculate an inverse only to retain the existing condition
measure, never to calculate M. LU and the old Gauss–Jordan algorithm can round
pivots/conditions differently near thresholds; bitwise equivalence is not promised.

TPS LU additionally exposes the determinant (diagonal product with backend swap
parity), a free determinant helper returning zero on exact singularity, and the
maximum absolute row-sum norm. The determinant helper applies no absolute pivot
floor. Ordinary scalar underflow remains possible; OPALX explicitly rejects
nonfinite inputs, solve outputs and diagnostic results. The canonical-form
residual remains max_ij abs((M^T J M-J)_ij), not a Frobenius norm. Mechanical slopes
remain noncanonical; these diagnostics do not alter or symplectify the map.

## Regression coverage

LU scalar support follows the selected Kokkos Getrf/Getrs backend; the 6D path
and these LU regressions use double. The standalone TPS regression checks pivoted
real transpose solves, multiple RHS columns, determinant parity with several 6D permutations, exact
singularity, invalid shapes/floors/nonfinite values and strict floor boundaries.
The OPALX tests check a known nonsymmetric right solve with coupled random maps,
row permutations and mixed coordinate scales; analytic original-vs-transpose
condition norms; well-conditioned tiny matrices; independent pivot/condition
limits; an original-vs-transpose pivot-floor counterexample; analytic block
triangular determinants; singular diagnostic behavior; an independent dense
canonical-form oracle; and explicit nonfinite/overflow failures. These tests run
on one and two MPI ranks, alongside existing drift, overlap, fringe, tabulated
field-map and Richardson regressions.

The fringe-map composition check originally required exact equality between a
legacy dense product and TPS composition. The new solve exposed a one-ULP
coefficient difference (0.14932478836759086 versus 0.14932478836759089). Both
products are now checked against an independent long-double six-term dot product
with the forward-error bound gamma_12 times sum(abs(a_k*b_k)), where
gamma_12=12*epsilon/(1-12*epsilon). This accounts for six products and additions;
there is no added absolute floor. Existing physical/finite-difference tolerances
are unchanged.

The COF Newton pivoted QR is intentionally retained for Part 1. Its relative
rank tolerance, coordinate scaling and damping remain unchanged; it is not
replaced by the 6D construction solve's absolute-pivot LU policy. This numerical-policy boundary was accepted for Part 1 on 2026-11-24. Any QR backend replacement is a
separate future numerical change, with rank and failure-classification regressions.
LAPACK and the full-tune winding calculation remain unchanged. Device tracking
and nonlinear construction are outside this migration.
