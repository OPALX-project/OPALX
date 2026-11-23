/**
 * @file BeamMoments.h
 * @brief Moments of one or several particle containers, and the function that computes them.
 *
 * Positions R are in metres and momenta P in beta*gamma, both in the frame the particles are in
 * when the moments are computed. Phase-space indices follow (x, px, y, py, z, pz).
 */

#ifndef OPALX_MOMENTS_BEAM_MOMENTS_H
#define OPALX_MOMENTS_BEAM_MOMENTS_H

#include "Moments/MomentRequest.h"
#include "PartBunch/PartBunchFwd.h"
#include "Types/Vector.h"

#include <array>
#include <cstddef>
#include <span>

namespace opalx::moments {

    /**
     * @brief Moments of a particle set at the time they were computed.
     *
     * A value holds the groups listed by available(). Reading a group that is not available throws
     * an OpalException, so a caller cannot use a quantity that was never computed. Averages
     * divide by max(N, 1): an empty set has zero means, spreads and bounds, and every derived
     * quantity whose denominator vanishes (correlation without spread, geometric emittance at rest)
     * is zero.
     */
    class BeamMoments {
    public:
        using Vector3 = ippl::Vector<double, 3>;

        [[nodiscard]] Request available() const { return available_m; }
        /// Global number of particles.
        [[nodiscard]] std::size_t count() const { return count_m; }

        [[nodiscard]] const Vector3& minR() const;
        /// Largest signed position per axis.
        [[nodiscard]] const Vector3& maxR() const;

        [[nodiscard]] const Vector3& meanR() const;

        [[nodiscard]] const Vector3& meanP() const;
        [[nodiscard]] double meanGamma() const;
        /// @f$\sqrt{\langle p_z\rangle^2 + 1}@f$ from the mean longitudinal momentum.
        [[nodiscard]] double meanGammaZ() const;
        /// Mean kinetic energy in MeV, with each container's reference mass when it has one.
        [[nodiscard]] double meanKineticEnergy() const;

        /// Central second moment of phase-space coordinates @p i and @p j (each below 6).
        [[nodiscard]] double covariance(unsigned i, unsigned j) const;
        [[nodiscard]] Vector3 rmsR() const;
        [[nodiscard]] Vector3 rmsP() const;
        /// Correlation @f$\langle x p_x\rangle / (\sigma_x \sigma_{p_x})@f$ per axis.
        [[nodiscard]] Vector3 correlationRP() const;
        /// rms emittance @f$\sqrt{\sigma_x^2\sigma_{p_x}^2 - \langle x p_x\rangle^2}@f$ in m rad.
        [[nodiscard]] Vector3 normalizedEmittance() const;
        /// Normalized emittance over @f$\sqrt{\langle\gamma\rangle^2 - 1}@f$.
        [[nodiscard]] Vector3 geometricEmittance() const;
        /// rms kinetic energy spread in MeV.
        [[nodiscard]] double stdKineticEnergy() const;
        /**
         * @brief Temperature per axis in the centre-of-momentum frame, in K.
         *
         * Each momentum u = beta*gamma is boosted with beta_c = <u>/<gamma> into the frame where
         * the mean momentum vanishes, and k_B T_i = m c^2 <u'_i^2 / gamma'>. This equals k_B T for
         * a Maxwell-Juettner distribution and m c^2 rms_p_i^2 in the non-relativistic limit. The
         * mass is the one used for the kinetic energy.
         * @note Particles are sampled at equal lab time, not at equal rest-frame time.
         * @note Over several containers beta_c is not mass-weighted, so the value describes one
         *       species, for example a beam and its decay daughters of the same kind.
         */
        [[nodiscard]] const Vector3& restFrameTemperature() const;

        [[nodiscard]] const Vector3& meanPolarization() const;
        /// @f$\sqrt{\langle P_i^2\rangle - \langle P_i\rangle^2}@f$ per component.
        [[nodiscard]] const Vector3& rmsPolarization() const;
        /// @f$|\langle \vec P\rangle|@f$, the depolarization observable.
        [[nodiscard]] double meanPolarizationMagnitude() const;

    private:
        friend BeamMoments compute(std::span<const ParticleContainer_t* const>, Request);

        void require(Request groups) const;

        Request available_m;
        std::size_t count_m = 0;

        Vector3 minR_m{0.0};
        Vector3 maxR_m{0.0};
        Vector3 meanR_m{0.0};
        Vector3 meanP_m{0.0};
        double meanGamma_m         = 0.0;
        double meanKineticEnergy_m = 0.0;

        std::array<std::array<double, 6>, 6> covariance_m{};
        double stdKineticEnergy_m = 0.0;
        Vector3 restFrameTemperature_m{0.0};

        Vector3 meanPolarization_m{0.0};
        Vector3 rmsPolarization_m{0.0};
    };

    /**
     * @brief Compute the requested groups over the union of @p containers.
     *
     * Each pass reduces every container into one set of sums with one allreduce per MPI operator,
     * so several beams cost no extra communication. Empty containers contribute nothing. The
     * groups a requested group depends on are computed too (Spread needs both means).
     * Polarization is available only if every container tracks spin.
     *
     * Collective on ippl::Comm: every rank calls it with the same containers and request, also
     * ranks without local particles. It blocks until the device reductions have finished.
     * @pre The containers' particles are in one coordinate frame.
     */
    [[nodiscard]] BeamMoments compute(
            std::span<const ParticleContainer_t* const> containers, Request request);

    /** @brief compute() for one container. */
    [[nodiscard]] BeamMoments compute(const ParticleContainer_t& particles, Request request);

}  // namespace opalx::moments

#endif  // OPALX_MOMENTS_BEAM_MOMENTS_H
