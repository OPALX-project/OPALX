/** @file CanonicalModes6D.h
 * @brief Host normal modes and Williamson normalization for canonical 6D RMS covariances.
 */
#ifndef OPALX_CANONICAL_MODES_6D_H
#define OPALX_CANONICAL_MODES_6D_H

#include <array>
#include <complex>
#include <vector>
#include "Algorithms/SpaceChargeMap.h"

namespace opalx::maps::modes6d {
    using Matrix     = SpaceChargeMatrix;
    using Emittances = std::array<double, 3>;

    /** @brief Numerical acceptance criteria, not physical stability margins.
     * scales equilibrate z=D*z_scaled; they are positive and need not be symplectic.
     * The scaled canonical form is D^T*J*D, never bare J. Phase distances use turns.
     * Residuals use Frobenius norms; covariance residuals use diagonal equilibration.
     * The symplectic residual is ||N^T G N-G||_F/||G||_F for N=D^-1 M D and
     * G=D^T J D/max_k(D_qk D_pk). Basis normalization uses ||S^T J S-J||_F/sqrt(6).
     * Eigenpair residuals divide ||N v-lambda v||_2 by max(||N||_F,|lambda|)*||v||_2.
     * The reciprocal condition uses singular values of LAPACK's real packed eigenbasis.
     */
    struct Settings {
        std::array<double, 6> scales{1, 1, 1, 1, 1, 1};
        double modulusTolerance              = 1e-8;
        double eigenResidualTolerance        = 1e-10;
        double symplecticTolerance           = 1e-7;
        double covarianceTolerance           = 1e-7;
        double minimumPhaseSeparation        = 1e-6;
        double minimumBasisRcond             = 1e-10;
        double minimumRelativeSymplecticNorm = 1e-10;
        double minimumSubspaceOverlap        = 0.5;  ///< Mode-tracking policy, not a physics error.
        double minimumAssignmentGap = 1e-3;  ///< Best-minus-next total overlap of three modes.
    };

    enum class Status {
        Stable,         ///< Three separated, resolved elliptic canonical pairs.
        NonSymplectic,  ///< Canonical-form residual exceeds tolerance.
        Unstable,       ///< At least one growing eigenvalue.
        NonUnitCircle,  ///< Nonconservative spectrum without a growing eigenvalue.
        Marginal,       ///< Real unit eigenvalue or phase close to zero/half a turn.
        Degenerate,     ///< Oscillatory phases are unresolved; no unique mode assignment.
        IllConditioned  ///< Eigenbasis or symplectic normalization is unreliable.
    };

    struct Mode {
        std::complex<double> eigenvalue;  ///< Positive-imaginary member, not a full tune.
        double phase  = 0;                ///< Principal phase in (0,1/2), in turns.
        int kreinSign = 0;                ///< sign(Re(v)^T J Im(v)), either +1 or -1.
    };

    /** @brief Spectrum and accepted real canonical basis S, with S^T J S=J.
     * J=diag([0,1;-1,0],...). For mode k, the eigenvector is S[:,2k]+i*sign*S[:,2k+1].
     * Only Stable results have modes/basis; other statuses retain spectral diagnostics.
     * Initial order is increasing principal phase, with no x/y/z or full-tune assignment.
     * For an eigenvector v=a+i b, normalize by sqrt(|a^T J b|), and multiply the
     * second real basis column by sign(a^T J b). Both Krein signs are accepted.
     * Repeated phases are conservatively unsupported even when the map is semisimple.
     */
    struct Analysis {
        Status status = Status::Marginal;
        Matrix map{}, basis{};
        std::array<std::complex<double>, 6> eigenvalues{};
        std::vector<Mode> modes;
        double maximumEigenResidual     = 0;
        double basisReciprocalCondition = 0;
        double symplecticResidual       = 0;
        double normalizationResidual    = 0;
    };

    /** @brief Analyze a physical canonical matrix in (q1,p1,q2,p2,q3,p3).
     * q uses metres; p uses species p/(mc), hence action/emittance uses metres, without pi.
     * Real LAPACK on host values; no Kokkos runtime, MPI collectives or device work.
     * Spectral stability of this frozen particle map is NOT collective envelope stability.
     * Degenerate/marginal cases are unsupported, not proof of physical instability.
     * No eigenvalue snapping, symplectification or longitudinal constraint is imposed.
     * NonSymplectic takes precedence over spectral statuses. NonStable results never
     * expose a partial normalized basis. Finite inputs may still fail if chosen scales
     * make intermediate values unrepresentable; use appropriate numerical scales.
     * @throws OpalException Invalid inputs/settings, LAPACK failure or inaccurate eigenpairs.
     * @test TestCanonicalModes6D checks coupled canonical similarities, both Krein signs,
     * extreme scales, instability, marginal/degenerate spectra and one/two-rank agreement.
     */
    Analysis analyze(const Matrix& map, const Settings& settings = {});

    /** @brief First-order, origin-preserving, same-chart canonical TPS adapter.
     * Requires exact boundary metadata equality, interleaved canonical pairs, q [m], p [p/(mc)].
     * Removes local coordinate scales exactly once. Requires active Kokkos/serialized TPS use.
     * The caller must establish periodic orbit/time/frame/gauge identification explicitly;
     * an ordinary transit map with different boundary reference times is rejected.
     */
    Analysis analyze(const tps::Map<double>& map, const Settings& settings = {});

    /** @brief Sigma=S diag(eps1,eps1,eps2,eps2,eps3,eps3) S^T.
     * Three strictly positive normalized RMS mode emittances [m], without pi; ordered by the
     * supplied analysis, not by laboratory planes. Checks SPD and M Sigma M^T=Sigma.
     * A negative Krein sign still contributes a positive covariance. Neutral longitudinal
     * motion cannot supply a unique stable 6D covariance and is rejected.
     */
    Matrix invariantCovariance(const Analysis&, const Emittances&, const Settings& settings = {});

    /** @brief Williamson decomposition Sigma=S diag(eps_k,eps_k) S^T, S^T J S=J.
     * eps_k [m] are sorted increasingly. Equal emittances are supported but their basis is
     * nonunique. No physical mode labels or continuation identity are implied.
     */
    struct Williamson {
        Matrix basis{};
        Emittances emittances{};
        double normalizationResidual  = 0;
        double reconstructionResidual = 0;
    };

    /** @brief Cholesky + real skew-Schur decomposition of a finite symmetric SPD covariance.
     * Carries D^T J D through numerical equilibration; no mixed-unit eigenvalue floor or
     * clipping of nonpositive emittances. Unresolvable/nonrepresentable results throw.
     * With Sigma=D c L L^T D and G=D^T J D/w, form A=L^T G L. Real Schur vectors
     * Q are paired/oriented so Q^T A Q=diag(kappa_k J2), kappa_k>0. Then
     * S=D L Q diag(1/sqrt(w kappa_k)) and eps_k=c w kappa_k. The repeated factor
     * for each pair ensures S^T J S=J and reconstructs Sigma. c and w are positive
     * numerical normalization factors; no physical action is discarded.
     * Covariance symmetry/SPD checks use diagonal equilibration before Cholesky;
     * reconstruction checks use the same diagonal metric. Equal eps do not require
     * a unique Q. This factorization is independent of the particle-map spectrum.
     * @test TestCanonicalModes6D checks reconstruction, trace invariants of (J Sigma)^2p,
     * prescribed-action restoration, equal emittances and covariance mixing.
     */
    Williamson williamson(const Matrix& covariance, const Settings& settings = {});
    /** @brief Restore prescribed eigenemittances in the supplied Williamson basis.
     * Targets correspond to its three basis pairs. This preserves SPD and canonical actions;
     * it is not a periodic fixed-point solver or an automatic physical mode assignment.
     */
    Matrix withEmittances(const Williamson&, const Emittances&, const Settings& settings = {});

    /** @brief Reorder current modes/basis into reference identity using 2D subspace overlap.
     * Overlap is phase independent in the fixed numerical scales. Enumerates all six
     * assignments, requires adequate overlap, a unique best assignment and matching Krein
     * signs. Ambiguous crossings throw rather than silently exchange prescribed emittances.
     * Only permutes the current mode pairs. A prescribed emittance array subsequently
     * follows reference identity, so it must not be re-sorted by principal phase.
     */
    Analysis trackModes(
            const Analysis& current, const Analysis& reference, const Settings& settings = {});
}  // namespace opalx::maps::modes6d
#endif
