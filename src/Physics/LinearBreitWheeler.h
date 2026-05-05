#ifndef OPAL_LINEAR_BREIT_WHEELER_HH
#define OPAL_LINEAR_BREIT_WHEELER_HH

#include "OPALTypes.h"

#include <cstdint>
#include <random>

namespace Physics {
    namespace LinearBreitWheeler {

        /**
         * @brief Host-only cached data for repeated linear Breit-Wheeler event sampling.
         *
         * This structure stores the fixed two-photon geometry used in the first OPALX
         * Breit-Wheeler validation path. It is intentionally scalar and host-side,
         * matching the staged approach already used for linear Compton benchmarking.
         *
         * The incoming state consists of two photons:
         *
         * - a high-energy photon with energy @f$E_\gamma@f$ and direction
         *   @f$\hat n_\gamma@f$,
         * - a laser photon with energy @f$\omega_L@f$ and direction @f$\hat n_L@f$.
         *
         * The corresponding two-photon invariant is
         * @f[
         *   s = 2 E_\gamma \omega_L (1 - \hat n_\gamma \cdot \hat n_L).
         * @f]
         *
         * Pair creation is kinematically allowed only for @f$s \ge 4 m_e^2@f$.
         */
        struct SamplingKernel {
            double highEnergyPhotonEnergyGeV =
                    0.0;  ///< Incoming high-energy photon energy in the lab frame [GeV].
            double laserPhotonEnergyGeV =
                    0.0;  ///< Incoming laser-photon energy in the lab frame [GeV].
            Vector_t<double, 3> highEnergyDirection = Vector_t<double, 3>(
                    0.0);  ///< Normalized incoming high-energy photon direction.
            Vector_t<double, 3> laserDirection =
                    Vector_t<double, 3>(0.0);  ///< Normalized incoming laser-photon direction.
            double invariantSGeV2 = 0.0;       ///< Two-photon invariant @f$s@f$ in GeV^2.
            double pairBetaCM = 0.0;  ///< Outgoing lepton speed @f$\beta = \sqrt{1 - 4m_e^2/s}@f$
                                      ///< in the pair CM frame.
            double cmLorentzGamma =
                    0.0;  ///< Lorentz factor of the CM frame relative to the laboratory frame.
            Vector_t<double, 3> cmVelocity =
                    Vector_t<double, 3>(0.0);  ///< Dimensionless CM boost velocity.
            double rejectionUpperBound = 0.0;  ///< Conservative rejection-sampling envelope for the
                                               ///< unpolarized CAIN-aligned angular kernel.
        };

        /**
         * @brief One sampled linear Breit-Wheeler event.
         *
         * The first benchmark path keeps only the unpolarized kinematics. The event is
         * represented by the center-of-momentum scattering variables and by the final
         * laboratory-frame electron and positron four-momenta.
         */
        struct SampledEvent {
            double scatteringCosineCM =
                    0.0;             ///< Center-of-momentum scattering cosine @f$\cos\theta@f$.
            double azimuthCM = 0.0;  ///< Center-of-momentum azimuth @f$\phi@f$ in rad.
            double electronEnergyLabGeV =
                    0.0;  ///< Outgoing electron energy in the lab frame [GeV].
            double positronEnergyLabGeV =
                    0.0;  ///< Outgoing positron energy in the lab frame [GeV].
            Vector_t<double, 3> electronMomentumLabGeV = Vector_t<double, 3>(
                    0.0);  ///< Outgoing electron three-momentum in the lab frame [GeV/c].
            Vector_t<double, 3> positronMomentumLabGeV = Vector_t<double, 3>(
                    0.0);  ///< Outgoing positron three-momentum in the lab frame [GeV/c].
        };

        /**
         * @brief Convert a laser wavelength to a single-photon energy.
         *
         * The conversion uses
         * @f[
         *   \omega_L = \frac{2 \pi \hbar c}{\lambda_L}.
         * @f]
         *
         * @param wavelength_m Laser wavelength in m.
         * @return Laser photon energy in GeV.
         */
        double photonEnergyFromWavelengthGeV(double wavelength_m);

        /**
         * @brief Two-photon invariant @f$s@f$ for Breit-Wheeler pair production.
         *
         * For two incoming photons with energies @f$E_1@f$ and @f$E_2@f$ and unit
         * directions @f$\hat n_1@f$ and @f$\hat n_2@f$,
         * @f[
         *   s = 2 E_1 E_2 (1 - \hat n_1 \cdot \hat n_2).
         * @f]
         *
         * @param photon1EnergyGeV First photon energy in GeV.
         * @param photon2EnergyGeV Second photon energy in GeV.
         * @param photon1Direction First photon direction.
         * @param photon2Direction Second photon direction.
         * @return Invariant @f$s@f$ in GeV^2.
         */
        double invariantSGeV2(
                double photon1EnergyGeV, double photon2EnergyGeV,
                const Vector_t<double, 3>& photon1Direction,
                const Vector_t<double, 3>& photon2Direction);

        /**
         * @brief Threshold invariant for @f$\gamma + \gamma \to e^- + e^+@f$.
         *
         * The exact kinematic threshold is
         * @f[
         *   s_\text{thr} = 4 m_e^2.
         * @f]
         */
        double thresholdInvariantSGeV2();

        /**
         * @brief Check whether pair creation is kinematically allowed.
         * @param invariantSGeV2 Two-photon invariant @f$s@f$ in GeV^2.
         * @return True if @f$s \ge 4 m_e^2@f$.
         */
        bool isAboveThreshold(double invariantSGeV2);

        /**
         * @brief Outgoing lepton speed in the pair center-of-momentum frame.
         *
         * Above threshold,
         * @f[
         *   \beta = \sqrt{1 - \frac{4m_e^2}{s}}.
         * @f]
         *
         * @param invariantSGeV2 Two-photon invariant @f$s@f$ in GeV^2.
         * @return Dimensionless center-of-momentum lepton speed.
         */
        double pairBetaCM(double invariantSGeV2);

        /**
         * @brief Total unpolarized linear Breit-Wheeler cross section.
         *
         * The first OPALX implementation uses the standard unpolarized Breit-Wheeler
         * total cross section in terms of the center-of-momentum lepton speed
         * @f$\beta = \sqrt{1 - 4m_e^2/s}@f$:
         * @f[
         *   \sigma_{\gamma\gamma\to e^+e^-} = \frac{\pi r_e^2}{2}(1-\beta^2)
         *   \left[
         *     (3-\beta^4) \ln\!\left(\frac{1+\beta}{1-\beta}\right)
         *     - 2\beta(2-\beta^2)
         *   \right].
         * @f]
         *
         * @param invariantSGeV2 Two-photon invariant @f$s@f$ in GeV^2.
         * @return Total cross section in m^2.
         */
        double totalCrossSection(double invariantSGeV2);

        /**
         * @brief Map CAIN's proposal parameter @f$z \in [-1,1]@f$ to the CM scattering cosine.
         *
         * The linear CAIN generator samples a proposal variable @f$z@f$, transforms it to
         * an auxiliary variable @f$Y(z)@f$, and then sets
         * @f[
         *   \cos\theta = 1 - 2Y.
         * @f]
         * This helper exposes that mapping directly so the local angular kernel can be
         * tested pointwise in the same parameterization used by the rejection sampler.
         *
         * @param invariantSGeV2 Two-photon invariant @f$s@f$ in GeV^2.
         * @param proposalZ CAIN proposal variable @f$z@f$ in [-1,1].
         * @return Center-of-momentum scattering cosine corresponding to @p proposalZ.
         */
        double proposalZToScatteringCosineCM(double invariantSGeV2, double proposalZ);

        /**
         * @brief Evaluate the unpolarized CAIN-aligned angular kernel in proposal coordinates.
         *
         * This is the dimensionless angular weight used by the current OPALX rejection
         * sampler after the CAIN change of variables from the proposal parameter @f$z@f$
         * to the center-of-momentum scattering angle. It is not the physical cross section
         * in SI units; it is the local kernel shape that controls accepted event sampling.
         *
         * @param invariantSGeV2 Two-photon invariant @f$s@f$ in GeV^2.
         * @param proposalZ CAIN proposal variable @f$z@f$ in [-1,1].
         * @return Non-negative unpolarized angular weight.
         */
        double unpolarizedAngularWeight(double invariantSGeV2, double proposalZ);

        /**
         * @brief Build a cached sampling kernel for repeated sampled events.
         *
         * The initial OPALX sampled-event path mirrors the linear CAIN event generator:
         * it keeps the fixed two-photon geometry, computes the corresponding CM boost,
         * and precomputes a rejection-sampling envelope for the unpolarized CAIN-style
         * angular weight.
         */
        SamplingKernel makeSamplingKernel(
                double highEnergyPhotonEnergyGeV, double laserPhotonEnergyGeV,
                const Vector_t<double, 3>& highEnergyDirection,
                const Vector_t<double, 3>& laserDirection);

        /**
         * @brief Build a deterministic host RNG from `Options::seed`.
         * @param streamIndex Optional stream index for reproducible stream splitting.
         * @return Seeded host-side random engine.
         */
        std::mt19937_64 makeHostRandomEngine(std::uint64_t streamIndex = 0);

        /**
         * @brief Sample one unpolarized linear Breit-Wheeler event.
         *
         * The angular proposal and acceptance rule are aligned with the linear CAIN
         * event generator `LNBWGN`, but restricted here to the unpolarized case. The
         * sampled outgoing electron and positron momenta are returned in the lab frame.
         *
         * @param kernel Cached sampling kernel created by @ref makeSamplingKernel.
         * @param engine Host-side random engine.
         * @return Sampled event.
         */
        SampledEvent sampleEvent(const SamplingKernel& kernel, std::mt19937_64& engine);

    }  // namespace LinearBreitWheeler
}  // namespace Physics

#endif  // OPAL_LINEAR_BREIT_WHEELER_HH
