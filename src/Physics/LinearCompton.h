#ifndef OPAL_LINEAR_COMPTON_HH
#define OPAL_LINEAR_COMPTON_HH

#include "OPALTypes.h"

#include <cstdint>
#include <random>

namespace Physics {
namespace LinearCompton {

/**
 * @brief Host-only cached data for repeated linear-Compton event sampling.
 *
 * The first Monte Carlo validation path in OPALX is intentionally scalar and
 * host-side. This structure caches the fixed beam / laser geometry, the
 * incoming photon energy in the electron rest frame, and the rejection-sampler
 * envelope used to sample @f$\cos\Theta^*@f$ from the unpolarized
 * Klein-Nishina kernel.
 *
 * The kernel is designed for deterministic testing and validation:
 * create it once for a fixed benchmark geometry, seed a host RNG from
 * `Options::seed`, and then draw sampled events via @ref sampleEvent.
 */
struct SamplingKernel {
    double electronTotalEnergyGeV = 0.0;  ///< Incoming electron total energy in the lab frame [GeV].
    double laserPhotonEnergyGeV = 0.0;  ///< Incoming laser-photon energy in the lab frame [GeV].
    Vector_t<double, 3> beamDirection = Vector_t<double, 3>(0.0);  ///< Normalized incoming electron direction in the lab frame.
    Vector_t<double, 3> laserDirection = Vector_t<double, 3>(0.0);  ///< Normalized incoming laser direction in the lab frame.
    double incomingPhotonEnergyERFGeV = 0.0;  ///< Laser-photon energy after Lorentz transforming into the electron rest frame [GeV].
    double rejectionUpperBoundSolidAngleERF = 0.0;  ///< Conservative upper envelope for @f$d\sigma/d\Omega^*@f$ used by the rejection sampler.
};

/**
 * @brief Sampled single-photon Compton event.
 *
 * The event is parameterized by the electron-rest-frame scattering variables
 * @f$(\cos\Theta^*, \phi^*)@f$ and the corresponding scattered photon energy
 * and direction after boosting back to the laboratory frame.
 */
struct SampledEvent {
    double scatteringCosineERF = 0.0;  ///< Rest-frame polar scattering cosine @f$\cos\Theta^*@f$.
    double azimuthERF = 0.0;  ///< Rest-frame azimuth @f$\phi^*@f$ around the incoming-photon axis [rad].
    double scatteredPhotonEnergyERFGeV = 0.0;  ///< Outgoing photon energy in the electron rest frame [GeV].
    double scatteredPhotonEnergyLabGeV = 0.0;  ///< Outgoing photon energy in the laboratory frame [GeV].
    Vector_t<double, 3> scatteredPhotonDirectionLab = Vector_t<double, 3>(0.0);  ///< Normalized outgoing photon direction in the laboratory frame.
};

/**
 * @brief Convert a laser wavelength to a single-photon energy.
 *
 * The conversion uses
 * @f[
 *   \\omega_L = \\frac{2 \\pi \\hbar c}{\\lambda_L}.
 * @f]
 *
 * @param wavelength_m Laser wavelength in m.
 * @return Laser photon energy in GeV.
 * @throws OpalException If @p wavelength_m is not strictly positive.
 */
double photonEnergyFromWavelengthGeV(double wavelength_m);

/**
 * @brief Compute the electron Lorentz factor from total energy.
 * @param electronTotalEnergyGeV Electron total energy in GeV.
 * @return Dimensionless Lorentz factor @f$\\gamma = E_e / m_e@f$.
 * @throws OpalException If the energy is below the electron rest energy.
 */
double electronGamma(double electronTotalEnergyGeV);

/**
 * @brief Compute the electron speed in units of c from total energy.
 * @param electronTotalEnergyGeV Electron total energy in GeV.
 * @return Dimensionless speed @f$\\beta = \\sqrt{1 - 1/\\gamma^2}@f$.
 * @throws OpalException If the energy is below the electron rest energy.
 */
double electronBeta(double electronTotalEnergyGeV);

/**
 * @brief Compute the incoming photon energy in the electron rest frame.
 *
 * For a laboratory-frame beam direction @f$\\hat n_e@f$ and laser direction
 * @f$\\hat n_L@f$, the electron-rest-frame incoming photon energy is
 * @f[
 *   \\omega_1 = \\gamma \\omega_L (1 - \\beta \\cos \\alpha),
 *   \\qquad \\cos \\alpha = \\hat n_e \\cdot \\hat n_L.
 * @f]
 *
 * @param electronTotalEnergyGeV Electron total energy in GeV.
 * @param laserPhotonEnergyGeV Laser photon energy in GeV.
 * @param beamDirection Laboratory-frame incoming electron direction.
 * @param laserDirection Laboratory-frame incoming laser direction.
 * @return Incoming photon energy in the electron rest frame in GeV.
 * @throws OpalException If the energy is unphysical or either direction is zero.
 */
double restFrameIncomingPhotonEnergyGeV(double electronTotalEnergyGeV,
                                        double laserPhotonEnergyGeV,
                                        const Vector_t<double, 3>& beamDirection,
                                        const Vector_t<double, 3>& laserDirection);

/**
 * @brief Compute the dimensionless Klein-Nishina parameter @f$\\kappa = \\omega_1/m_e@f$.
 * @param incomingPhotonEnergyERFGeV Incoming photon energy in the electron rest frame in GeV.
 * @return Dimensionless rest-frame photon energy parameter.
 * @throws OpalException If the energy is not strictly positive.
 */
double invariantKappa(double incomingPhotonEnergyERFGeV);

/**
 * @brief Compute CAIN's linear-Compton invariant @f$x = 2 \\omega_1 / m_e@f$.
 * @param incomingPhotonEnergyERFGeV Incoming photon energy in the electron rest frame in GeV.
 * @return Dimensionless invariant @f$x@f$.
 * @throws OpalException If the energy is not strictly positive.
 */
double invariantX(double incomingPhotonEnergyERFGeV);

/**
 * @brief Compute CAIN's linear-Compton invariant directly from laboratory quantities.
 * @param electronTotalEnergyGeV Electron total energy in GeV.
 * @param laserPhotonEnergyGeV Laser photon energy in GeV.
 * @param beamDirection Laboratory-frame incoming electron direction.
 * @param laserDirection Laboratory-frame incoming laser direction.
 * @return Dimensionless invariant @f$x = 2 \\omega_1 / m_e@f$.
 */
double invariantX(double electronTotalEnergyGeV,
                  double laserPhotonEnergyGeV,
                  const Vector_t<double, 3>& beamDirection,
                  const Vector_t<double, 3>& laserDirection);

/**
 * @brief Minimum scattered photon energy in the electron rest frame.
 *
 * This is the back-scattered endpoint at @f$\\cos\\Theta^* = -1@f$:
 * @f[
 *   \\omega_{2,\\min} = \\frac{\\omega_1}{1 + 2 \\omega_1/m_e}.
 * @f]
 */
double scatteredPhotonEnergyMinERFGeV(double incomingPhotonEnergyERFGeV);

/**
 * @brief Maximum scattered photon energy in the electron rest frame.
 *
 * This is the forward-scattered endpoint at @f$\\cos\\Theta^* = 1@f$:
 * @f[
 *   \\omega_{2,\\max} = \\omega_1.
 * @f]
 */
double scatteredPhotonEnergyMaxERFGeV(double incomingPhotonEnergyERFGeV);

/**
 * @brief Compute the scattered photon energy in the electron rest frame.
 *
 * The Compton recoil formula is
 * @f[
 *   \\omega_2(\\Theta^*) =
 *   \\frac{\\omega_1}{1 + (\\omega_1/m_e)(1 - \\cos\\Theta^*)}.
 * @f]
 *
 * @param incomingPhotonEnergyERFGeV Incoming photon energy in the electron rest frame in GeV.
 * @param scatteringCosineERF Rest-frame scattering cosine @f$\\cos\\Theta^*@f$.
 * @return Scattered photon energy in GeV.
 * @throws OpalException If the energy is unphysical or the cosine is out of range.
 */
double scatteredPhotonEnergyERFGeV(double incomingPhotonEnergyERFGeV,
                                   double scatteringCosineERF);

/**
 * @brief Unpolarized Klein-Nishina differential cross section in solid angle.
 *
 * The exact electron-rest-frame kernel is
 * @f[
 *   \frac{d\sigma}{d\Omega^*}
 *   = \frac{r_e^2}{2}
 *     \left(\frac{\omega_2}{\omega_1}\right)^2
 *     \left(
 *       \frac{\omega_1}{\omega_2}
 *       + \frac{\omega_2}{\omega_1}
 *       - \sin^2\Theta^*
 *     \right),
 * @f]
 * where @f$\omega_2@f$ is the Compton-recoiled photon energy at the same
 * scattering cosine.
 *
 * @param incomingPhotonEnergyERFGeV Incoming photon energy in the electron rest frame in GeV.
 * @param scatteringCosineERF Rest-frame scattering cosine @f$\cos\Theta^*@f$.
 * @return Differential cross section in units of m^2 / sr.
 */
double differentialCrossSectionSolidAngleERF(double incomingPhotonEnergyERFGeV,
                                             double scatteringCosineERF);

/**
 * @brief Incoming laser-photon direction in the electron rest frame.
 *
 * The laboratory-frame photon direction is Lorentz transformed into the rest
 * frame of an electron moving along @p beamDirection. The returned direction is
 * normalized and can be used as the axis of the Klein-Nishina scattering cone.
 *
 * @param electronTotalEnergyGeV Electron total energy in GeV.
 * @param beamDirection Laboratory-frame incoming electron direction.
 * @param laserDirection Laboratory-frame incoming laser direction.
 * @return Unit vector for the incoming photon in the electron rest frame.
 */
Vector_t<double, 3> restFrameIncomingPhotonDirection(double electronTotalEnergyGeV,
                                                     const Vector_t<double, 3>& beamDirection,
                                                     const Vector_t<double, 3>& laserDirection);

/**
 * @brief Laboratory-frame photon direction from electron-rest-frame scattering angles.
 *
 * The outgoing photon direction is defined in the electron rest frame by
 * @f$\cos\Theta^*@f$ and @f$\phi^*@f$ around the incoming photon axis. The
 * corresponding direction is then Lorentz boosted back to the laboratory
 * frame along the incoming electron beam direction.
 *
 * @param electronTotalEnergyGeV Electron total energy in GeV.
 * @param beamDirection Laboratory-frame incoming electron direction.
 * @param laserDirection Laboratory-frame incoming laser direction.
 * @param scatteringCosineERF Rest-frame scattering cosine @f$\cos\Theta^*@f$.
 * @param azimuthERF Rest-frame azimuth angle @f$\phi^*@f$ in rad.
 * @return Unit vector for the scattered photon direction in the laboratory frame.
 */
Vector_t<double, 3> labPhotonDirection(double electronTotalEnergyGeV,
                                       const Vector_t<double, 3>& beamDirection,
                                       const Vector_t<double, 3>& laserDirection,
                                       double scatteringCosineERF,
                                       double azimuthERF);

/**
 * @brief Laboratory-frame photon energy from electron-rest-frame scattering angles.
 *
 * The outgoing photon is defined by a rest-frame scattering cosine
 * @f$\cos\Theta^*@f$ relative to the incoming photon direction and an azimuth
 * angle @f$\phi^*@f$ around that axis. The photon is then boosted back to the
 * laboratory frame along the electron beam direction.
 *
 * @param electronTotalEnergyGeV Electron total energy in GeV.
 * @param laserPhotonEnergyGeV Laser photon energy in GeV.
 * @param beamDirection Laboratory-frame incoming electron direction.
 * @param laserDirection Laboratory-frame incoming laser direction.
 * @param scatteringCosineERF Rest-frame scattering cosine @f$\cos\Theta^*@f$.
 * @param azimuthERF Rest-frame azimuth angle @f$\phi^*@f$ in rad.
 * @return Scattered photon energy in the laboratory frame in GeV.
 */
double labPhotonEnergyGeV(double electronTotalEnergyGeV,
                          double laserPhotonEnergyGeV,
                          const Vector_t<double, 3>& beamDirection,
                          const Vector_t<double, 3>& laserDirection,
                          double scatteringCosineERF,
                          double azimuthERF);

/**
 * @brief Unpolarized Klein-Nishina differential cross section in the scattered-photon energy.
 *
 * The underlying angular kernel is
 * @f[
 *   \\frac{d\\sigma}{d\\Omega^*}
 *   = \\frac{r_e^2}{2}
 *     \\left(\\frac{\\omega_2}{\\omega_1}\\right)^2
 *     \\left(
 *       \\frac{\\omega_1}{\\omega_2}
 *       + \\frac{\\omega_2}{\\omega_1}
 *       - \\sin^2\\Theta^*
 *     \\right),
 * @f]
 * and the one-dimensional energy spectrum integrated over azimuth is
 * @f[
 *   \\frac{d\\sigma}{d\\omega_2}
 *   = 2\\pi \\frac{d\\sigma}{d\\Omega^*} \\frac{m_e}{\\omega_2^2}.
 * @f]
 *
 * The valid support is
 * @f$\\omega_2 \\in [\\omega_{2,\\min}, \\omega_{2,\\max}]@f$.
 *
 * @param incomingPhotonEnergyERFGeV Incoming photon energy in the electron rest frame in GeV.
 * @param scatteredPhotonEnergyERFGeV Scattered photon energy in the electron rest frame in GeV.
 * @return Differential cross section in units of m^2 / GeV.
 * @throws OpalException If the input energy is outside the physical support.
 */
double differentialCrossSectionOmegaERF(double incomingPhotonEnergyERFGeV,
                                        double scatteredPhotonEnergyERFGeV);

/**
 * @brief Total unpolarized Klein-Nishina cross section.
 *
 * With @f$\\kappa = \\omega_1/m_e@f$, the exact total cross section is
 * @f[
 *   \\sigma_{\\mathrm{KN}} = 2\\pi r_e^2
 *   \\left[
 *     \\frac{1+\\kappa}{\\kappa^3}
 *     \\left(
 *       \\frac{2\\kappa(1+\\kappa)}{1+2\\kappa} - \\ln(1+2\\kappa)
 *     \\right)
 *     + \\frac{\\ln(1+2\\kappa)}{2\\kappa}
 *     - \\frac{1+3\\kappa}{(1+2\\kappa)^2}
 *   \\right].
 * @f]
 *
 * @param incomingPhotonEnergyERFGeV Incoming photon energy in the electron rest frame in GeV.
 * @return Total cross section in m^2.
 * @throws OpalException If the energy is not strictly positive.
 */
double totalCrossSection(double incomingPhotonEnergyERFGeV);

/**
 * @brief Thomson limit cross section.
 * @return The classical Thomson cross section @f$8\\pi r_e^2 / 3@f$ in m^2.
 */
double thomsonCrossSection();

/**
 * @brief Rest-frame scattering cosine corresponding to a lab-forward photon.
 *
 * For a photon emitted parallel to the laboratory beam direction, the scattered
 * photon points along the boost axis in the electron rest frame. The scattering
 * cosine is therefore the aberrated incoming-photon cosine,
 * @f[
 *   \\cos\\Theta^*_{\\mathrm{lab\\,fwd}} =
 *   \\frac{\\cos\\alpha - \\beta}{1 - \\beta \\cos\\alpha}.
 * @f]
 */
double restFrameScatteringCosineForLabForwardPhoton(double electronTotalEnergyGeV,
                                                    const Vector_t<double, 3>& beamDirection,
                                                    const Vector_t<double, 3>& laserDirection);

/**
 * @brief Forward scattered photon energy in the laboratory frame.
 *
 * This helper reproduces the first OPALX benchmark path: compute the incoming
 * rest-frame photon energy, apply the Compton recoil at the scattering cosine
 * corresponding to a lab-forward photon, and boost the result back to the lab.
 *
 * @param electronTotalEnergyGeV Electron total energy in GeV.
 * @param laserPhotonEnergyGeV Laser photon energy in GeV.
 * @param beamDirection Laboratory-frame incoming electron direction.
 * @param laserDirection Laboratory-frame incoming laser direction.
 * @return Forward scattered photon energy in GeV.
 */
double labForwardPhotonEnergyGeV(double electronTotalEnergyGeV,
                                 double laserPhotonEnergyGeV,
                                 const Vector_t<double, 3>& beamDirection,
                                 const Vector_t<double, 3>& laserDirection);

/**
 * @brief Create a deterministic host-side random engine for Monte Carlo validation.
 *
 * The engine is seeded from `Options::seed`. For the validation path, a
 * negative seed falls back to a fixed constant instead of a time-dependent
 * seed so that unit tests remain reproducible.
 *
 * @param streamIndex Optional stream offset for independent deterministic streams.
 * @return Host-side Mersenne-Twister engine.
 */
std::mt19937_64 makeHostRandomEngine(std::uint64_t streamIndex = 0);

/**
 * @brief Build a cached sampling kernel for repeated event generation.
 *
 * This precomputes the incoming photon energy in the electron rest frame and a
 * safe rejection-sampling envelope for the unpolarized
 * @f$d\sigma/d\Omega^*@f$ kernel.
 *
 * The intended usage is:
 *
 * 1. fix one benchmark geometry or one electron state,
 * 2. call @ref makeSamplingKernel once,
 * 3. reuse the returned kernel for many calls to @ref sampleEvent.
 *
 * If the incoming electron energy or direction varies from particle to
 * particle, as in the finite-beam benchmark, a separate kernel must be built
 * for each sampled electron state before drawing the event.
 *
 * @param electronTotalEnergyGeV Electron total energy in GeV.
 * @param laserPhotonEnergyGeV Laser photon energy in GeV.
 * @param beamDirection Laboratory-frame incoming electron direction.
 * @param laserDirection Laboratory-frame incoming laser direction.
 * @return Cached host-side sampling kernel.
 */
SamplingKernel makeSamplingKernel(double electronTotalEnergyGeV,
                                  double laserPhotonEnergyGeV,
                                  const Vector_t<double, 3>& beamDirection,
                                  const Vector_t<double, 3>& laserDirection);

/**
 * @brief Sample one unpolarized linear-Compton event on the host.
 *
 * The sampler draws @f$\cos\Theta^*@f$ by rejection from the exact
 * Klein-Nishina angular kernel and samples @f$\phi^*@f$ uniformly on
 * @f$[0, 2\pi)@f$. The returned event contains both rest-frame and
 * laboratory-frame observables.
 *
 * This is a validation-oriented Monte Carlo path: it is scalar, host-only, and
 * deterministic when the engine is created from @ref makeHostRandomEngine with
 * a fixed `Options::seed`. It is not yet intended to define production photon
 * creation semantics for tracking.
 *
 * @param kernel Cached sampling kernel created by @ref makeSamplingKernel.
 * @param engine Host-side random engine, typically created by
 *        @ref makeHostRandomEngine.
 * @return Sampled linear-Compton event.
 */
SampledEvent sampleEvent(const SamplingKernel& kernel, std::mt19937_64& engine);

}  // namespace LinearCompton
}  // namespace Physics

#endif  // OPAL_LINEAR_COMPTON_HH
