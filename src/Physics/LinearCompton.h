#ifndef OPAL_LINEAR_COMPTON_HH
#define OPAL_LINEAR_COMPTON_HH

#include "OPALTypes.h"

namespace Physics {
namespace LinearCompton {

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

}  // namespace LinearCompton
}  // namespace Physics

#endif  // OPAL_LINEAR_COMPTON_HH
