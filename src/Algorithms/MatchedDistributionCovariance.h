#ifndef OPAL_MATCHED_DISTRIBUTION_COVARIANCE_H
#define OPAL_MATCHED_DISTRIBUTION_COVARIANCE_H

#include "Algorithms/ClosedOrbitOptics.h"

/** Tangent conversion from a fixed-section MATCH to common-time sampling.
 * The sampler coordinates are orbit-local (x,ux,y,uy,z,uz), in metres and p/(mc),
 * with +z along the solved momentum. SIGMAZ and SIGMAPZ prescribe the population
 * marginal widths of those local longitudinal coordinates. All input correlations
 * are replaced by the MATCH model and its required coordinate conversion.
 *
 * Let Q rotate orbit-local perturbations into the COF section and u be the
 * reference mechanical momentum there. Intersect the section using the local
 * magnetic flow: delta_t=-delta_r_z/v_z,
 * delta_r_section=delta_r-v*delta_r_z/v_z,
 * delta_u_section=delta_u-f*delta_r_z/v_z. Here
 * f/v_z=(q_e*c/E_rest)*(u cross B)/u_z; B is in tesla and E_rest in eV.
 * The 6x6 matrix G maps sampler deviations to
 * (section x,ux,y,uy, sampler z,uz). Thus, with zero prescribed model cross
 * correlations, Sigma_sample=G^-1*diag(Sigma_MATCH,sigma_z^2,sigma_uz^2)*G^-T.
 * This includes arrival-time/force correlations even at a non-field-free section.
 * It preserves the two specified longitudinal marginal widths. It does not
 * assert longitudinal equilibrium, or change the fixed-energy 4D map model.
 *
 * Host-only TPS/Kokkos algebra; requires an active Kokkos runtime. Value results
 * own no runtime storage. Both longitudinal widths must be positive and finite.
 */
class MatchedDistributionCovariance {
public:
    using Matrix = std::array<std::array<double, 6>, 6>;
    /// Local section-intersection Jacobian, including the reference magnetic force.
    static Matrix sectionProjection(const ClosedOrbitOptics&);
    /// Validate SPD and return the common-time covariance; no Gaussian cutoff bias here.
    static Matrix compute(
            const ClosedOrbitOptics&, const OneTurnMap::Matrix&, double sigmaZ, double sigmaPZ);
};
#endif
