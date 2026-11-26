// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
#ifndef OPAL_TPS_LINEAR_MAP_ALGEBRA_H
#define OPAL_TPS_LINEAR_MAP_ALGEBRA_H
#include "Algorithms/Matrix.h"

namespace opalx::maps {
    /**
     * @brief Solve M D_in=D_out using pivoted TPS/Kokkos LU on D_in.
     *
     * Solves D_in^T M^T=D_out^T with the same factorization. The reported condition
     * is ||D_in||_infinity ||D_in^-1||_infinity, computed by identity solves. It is
     * unscaled and mixes the six coordinate units, as in the original builder.
     * Reject exact singularity, LU pivot magnitudes <1e-18 and nonfinite conditions
     * or conditions >1e12. Equality is accepted. Backend rounding near these boundaries
     * can differ from Gauss–Jordan; the floor applies to the original matrix's LU,
     * not a separately factored transpose. Inputs/output must be finite.
     * Serialized host use with an active Kokkos runtime; input matrices are unchanged.
     * @throws OpalException On invalid input or a rejected/failed numerical solve.
     */
    matrix6x6_t solveInputDifferences(
            const matrix6x6_t& input, const matrix6x6_t& output, double& condition);

    /// Diagnostics in the existing mixed-unit 6D mechanical slope coordinates.
    struct LinearMapDiagnostics {
        double determinantResidual;  ///< |det(M)-1|; exactly singular maps give one.
        double symplecticResidual;   ///< max_ij |(M^T J M-J)_ij|, not a Frobenius norm.
    };
    /**
     * @brief Kokkos LU determinant and matrix products for the canonical-form residual.
     *
     * J=diag(J2,J2,J2), J2=[0,1;-1,0]. Slopes/mechanical momenta are not globally
     * canonical: these are diagnostics, not universally preserved physical invariants.
     * No symplectification is applied. Singular determinant calculation is nonfatal;
     * nonfinite inputs/results and other backend failures throw OpalException.
     */
    LinearMapDiagnostics diagnoseLinearMap(const matrix6x6_t& map);
}  // namespace opalx::maps
#endif
