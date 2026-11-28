/**
 * @file FieldTensor.h
 * @brief Host analytic free-space 3D self-field tensors for RMS moment transport.
 */
#ifndef OPALX_SPACE_CHARGE_ANALYTIC3D_FIELD_TENSOR_H
#define OPALX_SPACE_CHARGE_ANALYTIC3D_FIELD_TENSOR_H

#include <Kokkos_Array.hpp>
#include <cstddef>

namespace opalx::spacecharge::analytic3d {

    using Vector = Kokkos::Array<double, 3>;
    using Matrix = Kokkos::Array<Vector, 3>;

    /** @brief Spatial profiles; neither value defines a 6D particle sampler. */
    enum class Profile {
        UniformEllipsoid,  ///< Uniform spatial density, semiaxes sqrt(5) times the RMS sizes.
        GaussianRms        ///< Least-squares linear force for an unbounded 3D Gaussian.
    };

    /** @brief Controls the dimensionless ellipsoid quadrature, independent of SI scales. */
    struct IntegrationSettings {
        /// Relative error target for EACH principal coefficient, in [1e-13, 1e-2).
        double relativeTolerance = 1e-10;
        /// Maximum number of integration intervals, including the initial aspect-ratio partition.
        std::size_t maximumIntervals = 1024;
    };

    /** @brief Rest-frame linear electric field E'=gradient*r, r measured from the centroid. */
    struct TensorResult {
        Matrix gradient{};  ///< Symmetric electric-field gradient [V/m^2], including signed charge.
        /** Estimated componentwise absolute quadrature error [V/m^2].
         * Rotated using absolute eigenvector products. Excludes covariance/eigensolver roundoff,
         * physical model error and subsequent frame transformations; not a rigorous error bound.
         */
        Matrix absoluteError{};
        std::size_t integrandEvaluations =
                0;  ///< Three components evaluated together at each node.
        std::size_t integrationIntervals = 0;  ///< Final adaptive partition size.
    };

    /**
     * @brief Compute the rest-frame RMS electric-field gradient from spatial covariance.
     *
     * @param profile Uniform spatial ellipsoid or Gaussian RMS closure.
     * @param covariance Centered, symmetric positive-definite REST-FRAME covariance [m^2].
     * @param totalCharge Signed total physical bunch charge [C], not macroparticle count.
     * @param settings Dimensionless quadrature tolerance and work limit.
     *
     * With covariance eigenvalues sigma_i^2 and orthogonal eigenvectors R, define
     * \f[
     * N_i=\frac{\sigma_1\sigma_2\sigma_3}{2}
     * \int_0^\infty\frac{ds}{(\sigma_i^2+s)\sqrt{\prod_j(\sigma_j^2+s)}},
     * \qquad \sum_i N_i=1.
     * \f]
     * Writing k=1/(4*pi*epsilon_0), the principal coefficients are
     * \f[
     * A_i=\frac{kQ N_i}{\sigma_1\sigma_2\sigma_3}
     * \begin{cases}3/(5\sqrt{5})&\text{uniform ellipsoid},\\
     *                 1/(2\sqrt{\pi})&\text{Gaussian RMS}.\end{cases}
     * \f]
     * The result is R*diag(A_i)*R^T. The uniform ellipsoid has semiaxes
     * sqrt(5)*sigma_i and its INTERNAL field is exactly linear. It does not specify
     * a classical transverse KV phase-space distribution or the exterior field.
     *
     * For Gaussian matching A_RMS=<E' r^T>*covariance^-1=<dE'/dr>. Gaussian
     * convolution gives A_RMS(C)=A_core(2*C)=A_core(C)/(2*sqrt(2)). This is an RMS
     * closure, NOT the core gradient or a nonlinear field evaluator. It assumes
     * controlled emittances; it does not describe halo or nonlinear emittance growth.
     * See Sacherer, RMS Envelope Equations with Space Charge, PAC 1971, p.1105.
     *
     * Normalized eigenvalues and a scale-aware partition remove the SI size/charge
     * from adaptive Gauss-Kronrod integration. Equal axes need no special branch;
     * the returned tensor is independent of eigenvectors within a degenerate subspace.
     * Covariance asymmetry up to 64*epsilon*sqrt(C_ii)*sqrt(C_jj) is averaged.
     * Zero charge validates all inputs and returns exact zeros without quadrature.
     *
     * @note Host-only, no Kokkos initialization or MPI collective required. Free space,
     * electrostatic rest-frame model; no pipe, image charge, boost, species-force
     * factor, path length, TPS kick or periodic matching is applied here. A caller
     * must transform a common-time lab covariance to the rest frame consistently.
     * @throws OpalException Invalid profile/settings, nonfinite or non-SPD covariance,
     * nonfinite charge, quadrature exhaustion or unrepresentable numerical result.
     * @test TestAnalytic3DFieldTensor covers independent sphere/spheroid limits,
     * Gaussian force averaging, rotations, scaling, relativistic conversion and failures.
     */
    [[nodiscard]] TensorResult evaluateRestFrameTensor(
            Profile profile, const Matrix& covariance, double totalCharge,
            const IntegrationSettings& settings = {});

}  // namespace opalx::spacecharge::analytic3d
#endif
