/**
 * @file SpaceChargeMap.h
 * @brief Host 6D tangent charts and TPS kicks for analytic RMS space charge.
 */
#ifndef OPALX_SPACE_CHARGE_MAP_H
#define OPALX_SPACE_CHARGE_MAP_H

#include <array>
#include "Algorithms/TpsLinearTransferMap.h"
#include "SpaceCharge/Analytic3D/FieldTensor.h"

namespace opalx::maps {

    using SpaceChargeMatrix = Kokkos::Array<Kokkos::Array<double, 6>, 6>;

    /** @brief Transport centered covariance by a 6D first-order origin-preserving TPS map.
     * Returns J*Sigma*J^T with equilibrated SPD checks. Input covariance is in the map's input
     * chart; output is in its output chart. Reject nonlinear/order-zero maps, affine offsets,
     * nonfinite Jacobians and numerically non-SPD covariances. Active Kokkos/serialized host use.
     */
    SpaceChargeMatrix transportSpaceChargeCovariance(
            const tps::Map<double>& map, const SpaceChargeMatrix& covariance);

    /** @brief Reference event and EXTERNAL fields sampled there; fields use laboratory axes. */
    struct SpaceChargeReference {
        LinearTransferMapReference orbit;
        opalx::spacecharge::analytic3d::Vector
                electricFieldLab{};  ///< V/m; no self-field included.
        opalx::spacecharge::analytic3d::Vector magneticFieldLab{};  ///< T; no self-field included.
        double massEV  = 0;  ///< Positive species rest energy [eV], not macroparticle rest energy.
        double chargeE = 0;  ///< Signed species charge in elementary-charge units; zero is allowed.
    };

    /** @brief Six perturbations around a reference event; all position coordinates use metres. */
    enum class SpaceChargeChart {
        SectionSlopes,         ///< (x,x',y,y',zeta,delta), existing external-map chart.
        CommonTimeMechanical,  ///< (rx,ux,ry,uy,rz,uz), u=p/(mc), in the stored section axes.
        CommonTimeCanonical  ///< (rx,pi_x,ry,pi_y,rz,pi_z), pi in units mc, local symmetric gauge.
    };

    /**
     * @brief Invertible 6D tangent transformations at one positively crossed section.
     *
     * SectionSlopes uses zeta=-beta*c*(t-t0), delta=|u|/U-1, a=u_x/u_s and b=u_y/u_s.
     * Let u0 be the reference momentum in section axes, U=|u0|, gamma=sqrt(1+U^2),
     * beta=U/gamma, n=u0/U and kappa=chargeE*c/massEV. Fields E,B below are rotated
     * from lab to section axes. The external reference force gives
     * \f[g=\dot u_0/(\beta c)=\kappa[E/(\beta c)+n\mathbin{\times} B].\f]
     * A particle at the section at t=t0-zeta/(beta*c) is projected to time t0:
     * \f[r=(x,y,0)+n\zeta,\qquad \Delta u=\Delta u_{section}+g\zeta,\f]
     * \f[\Delta u_{section}=u_s(e_x-n n_x)\Delta x'
     *                         +u_s(e_y-n n_y)\Delta y'+u_0\delta.\f]
     * The inverse uses zeta=r_z/n_z, x=r_x-n_x*zeta, y=r_y-n_y*zeta, removes
     * g*zeta, and recovers slopes and relative total momentum from that section momentum.
     * This includes transverse reference slopes and longitudinal energy/arrival-time coupling.
     *
     * At common time choose the local differential symmetric gauge A=(B cross r)/2,
     * A(reference)=0. The normalized canonical momentum is
     * \f[\pi=\Delta u+\frac{\kappa}{2}(B\mathbin{\times}r).\f]
     * All three B components enter this transform. Its three kappa*B parameters [1/m]
     * participate in exact TPS chart matching. Both common-time charts also record species rest
     * energy, making the mc momentum normalization explicit even when B=0. The canonical form is
     * interleaved
     * J=diag(J2,J2,J2), J2=[0,1;-1,0]. A mechanical chart instead has C^T*J*C.
     * This is a local tangent gauge, not a global vector potential or a finite-amplitude
     * coordinate transformation. Reference translations/rotations remain explicit metadata.
     *
     * Requires a finite, right-handed orthonormal frame, finite external fields, positive
     * reference section momentum and positive finite species rest energy. Require n_s>1e-12,
     * the frame-orthonormality accuracy of TPS metadata, to resolve the crossing direction.
     * Reversed/grazing sections are rejected, not projected. Arrays are value-owned.
     * TPS algebra and covariance products require active Kokkos and serialized host use.
     * @test TestSpaceChargeMap checks independent event/force linearizations, inverse and
     * covariance round trips, canonical forms, tilted sections and external TPS composition.
     */
    class SpaceChargeCoordinates {
    public:
        using Matrix = SpaceChargeMatrix;
        explicit SpaceChargeCoordinates(const SpaceChargeReference& reference);

        /// Exact metadata for one of the three supported charts; invalid enum values throw.
        const tps::PhaseSpaceChart& chart(SpaceChargeChart kind) const;
        /// Linear Jacobian from one chart to another; no frame/reference shift is implicit.
        Matrix jacobian(SpaceChargeChart from, SpaceChargeChart to) const;
        /// First-order TPS coordinate map with explicit, nonempty traversal provenance.
        tps::Map<double> map(
                SpaceChargeChart from, SpaceChargeChart to, const std::string& occurrenceId) const;
        /** @brief Sigma_to=T*Sigma_from*T^T, validated symmetric/SPD in both charts.
         * Diagonal equilibration precedes Cholesky; no mixed-unit eigenvalue floor.
         * Asymmetry <=64*epsilon*sqrt(Sigma_ii)*sqrt(Sigma_jj) is averaged.
         * @throws OpalException Invalid/non-SPD covariance or nonfinite numerical result.
         */
        Matrix convertCovariance(const Matrix&, SpaceChargeChart from, SpaceChargeChart to) const;

        double gamma() const { return gamma_m; }
        double beta() const { return beta_m; }
        double massEV() const { return massEV_m; }
        double chargeE() const { return chargeE_m; }
        const opalx::spacecharge::analytic3d::Vector& direction() const { return direction_m; }

    private:
        std::array<tps::PhaseSpaceChart, 3> charts_m;
        Kokkos::Array<Matrix, 3> toMechanical_m{}, fromMechanical_m{};
        opalx::spacecharge::analytic3d::Vector direction_m{};
        double gamma_m, beta_m, massEV_m, chargeE_m;
    };

    /** @brief Thin map and field diagnostics, frozen at the supplied midpoint covariance. */
    struct SpaceChargeKick {
        tps::Map<double> polynomial;  ///< First-order map in the requested chart, zero affine term.
        opalx::spacecharge::analytic3d::TensorResult restTensor;
        opalx::spacecharge::analytic3d::Matrix restSpatialCovariance;  ///< m^2, section axes.
        opalx::spacecharge::analytic3d::Matrix forceGradient;  ///< F/(species charge) per r, V/m^2.
        opalx::spacecharge::analytic3d::Matrix forceAbsoluteError;  ///< Quadrature estimate, V/m^2.
    };

    /**
     * @brief Build a frozen-density analytic SC kick from a covariance in the specified chart.
     *
     * Convert the FULL input covariance to common-time mechanical coordinates before extracting
     * its spatial block C. For a single narrow-spread quasistatic reference boost, let
     * D=I+(gamma-1)*n*n^T. Evaluate the existing analytic tensor A_rest at C_rest=D*C*D.
     * At reference velocity the symmetric electric-equivalent Lorentz-force gradient is
     * \f[H=D A_{rest}D/\gamma.\f]
     * Using signed reference-path increment h [m], apply
     * \f[r^+=r^-,\qquad \Delta u^+=\Delta u^-+
     *             \frac{h\,chargeE}{\beta\,massEV}H r^-.\f]
     * massEV is in eV; no extra c or elementary-charge factor enters this path-length formula.
     * Negative h is permitted for inverse-map checks; h=0 or zero force gives exact identity.
     * A position-only kick also leaves the canonical gauge term unchanged. Symmetry of H makes
     * the canonical shear symplectic and preserves canonical covariance eigenemittances.
     * Maps in SectionSlopes are obtained by explicit chart conjugation, not relabelled canonical.
     *
     * The centroid stays on the reference, reference fields are external only, and the reference
     * velocity/boost and RMS density are frozen during this kick. UniformEllipsoid is its
     * internal linear field; GaussianRms is the RMS linear closure, not its nonlinear field.
     * Error diagnostics propagate quadrature estimates only; no relativistic E/B cancellation
     * is used to compute H. Existing LorentzBoost supplies an independent regression oracle.
     * No external half-transports, NS scheduling, overlaps, particle pusher or matching iteration
     * are implemented here. Requires active Kokkos and serialized host TPS algebra.
     * @throws OpalException Invalid covariance/inputs, analytic tensor failure, overflow or
     * invalid traversal provenance. The caller's global TPS truncation setting is preserved.
     */
    SpaceChargeKick makeSpaceChargeKick(
            const SpaceChargeCoordinates& coordinates, const SpaceChargeMatrix& covariance,
            SpaceChargeChart chart, opalx::spacecharge::analytic3d::Profile profile,
            double totalCharge, double pathLength, const std::string& occurrenceId,
            const opalx::spacecharge::analytic3d::IntegrationSettings& settings = {});
}  // namespace opalx::maps
#endif
