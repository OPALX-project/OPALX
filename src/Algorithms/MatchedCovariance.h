#ifndef OPAL_MATCHED_COVARIANCE_H
#define OPAL_MATCHED_COVARIANCE_H

#include "Algorithms/ClosedOrbitOptics.h"

/** Host-only zero-current transverse matching of a static magnetic RING.
 * Requires active Kokkos and serialized host TPS algebra. Uses the retained COF
 * linearization through TPS and the existing LAPACK eigensolver; no ray retracing,
 * map projection, eigenvalue snapping, particle generation or 6D equilibrium.
 *
 * On the fixed-energy section dz=dH=0, physical canonical momentum is p+Q*A,
 * with SI charge Q=q_e*e and u=p/(mc).
 * In a local symmetric-gauge differential chart, eta=C*xi with
 * C(1,2)=-b/2, C(3,0)=b/2 and unit diagonal, b=q_e*c*B_n/E_rest [1/m].
 * q_e is in elementary-charge units, E_rest in eV, B_n in tesla. Thus the
 * mechanical symplectic form is Omega=C^T*J*C, not generally J. Gauge gradients
 * change the canonical chart, not the returned mechanical covariance.
 * This is a tangent chart at the solved orbit, not a global vector potential.
 *
 * For each positive-imaginary eigenvector v=a+i*d of C*M*C^-1, undo the
 * eigensolver similarity scales and normalize |a^T J d|=1. Then
 * Sigma_can=sum eps_k*(a_k*a_k^T+d_k*d_k^T) and Sigma=C^-1*Sigma_can*C^-T.
 * eps_k are normalized RMS mode eigenemittances [m], without a pi factor.
 * Modes are ordered by ascending positive-imaginary fractional phase in (0,1/2),
 * not by plane or oriented/full tune. Both Krein signs are permitted.
 *
 * Distinct stable modes, a well-conditioned symplectic eigenbasis and SPD
 * covariance are required. The canonical form is interleaved:
 * \f[ J=\operatorname{diag}\left(\begin{pmatrix}0&1\\-1&0\end{pmatrix},
 *                                  \begin{pmatrix}0&1\\-1&0\end{pmatrix}\right). \f]
 * With D=diag(COF.scales), N=D^-1*(C*M*C^-1)*D and K=D^T*J*D,
 * the symplectic residual is ||N^T*K*N-K||_F/||K||_F. The periodic residual
 * is ||L*S*L^T-S||_F/||S||_F, L=D^-1*M*D, S=D^-1*Sigma*D^-T.
 * The mode normalization residual is ||B^T*J*B-J_sign||_F/||J_sign||_F,
 * with J_sign containing the measured mode signatures. Diagonal equilibration
 * precedes the LAPACK Cholesky positivity check; no mixed-unit pivot cutoff.
 * @test TestMatchedCovariance checks analytic Twiss, coupling, scale undoing,
 * magnetic charts and rejection paths. TestMatchCmd checks command attributes,
 * invalid named references, clones and atomic JSON output using C++ fixtures.
 *
 * Residuals are finite-precision diagnostics, not
 * certification of physical integration accuracy; refine DT and FDSTEP.
 * Relative COF kinetic-energy drift must be <=1e-8, the existing eigenanalysis
 * unit-circle band, before a fixed-energy covariance is accepted.
 */
class MatchedCovariance {
public:
    using Matrix = OneTurnMap::Matrix;
    struct Settings {
        /// Relative Frobenius ceilings in the COF similarity-scaled coordinates.
        double symplecticTolerance = 1e-7;
        double periodicTolerance   = 1e-7;
        /// Fractional phase separation from each other and from 0 and 1/2.
        double minimumModeSeparation = 1e-6;
    };
    struct Result {
        Matrix covariance{};  ///< Mechanical section covariance; units are products of m,p/(mc).
        std::array<double, 2> emittances{}, tunes{};
        double symplecticResidual = 0, periodicResidual = 0, normalizationResidual = 0;
    };
    /// Local mechanical-to-canonical tangent matrix; validates finite species/field data.
    static Matrix canonicalTransform(double normalFieldTesla, double chargeE, double massEV);
    /// Throws OpalException for invalid inputs or unsupported/uncertified matching.
    static Result compute(const ClosedOrbitOptics&, const std::array<double, 2>&, const Settings&);
    static Result compute(const ClosedOrbitOptics& optics, const std::array<double, 2>& eps) {
        return compute(optics, eps, Settings{});
    }
};
#endif
