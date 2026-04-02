#include "Physics/LinearCompton.h"

#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace {
    Vector_t<double, 3> normalizeDirection(const Vector_t<double, 3>& direction,
                                           const char* where,
                                           const char* argument) {
        const double norm = std::sqrt(dot(direction, direction));
        if (norm <= 0.0) {
            throw OpalException(where, std::string("\"") + argument + "\" must be a non-zero vector.");
        }

        return direction / norm;
    }

    double clampCosine(double value) {
        return std::max(-1.0, std::min(1.0, value));
    }

    void ensureStrictlyPositive(double value, const char* where, const char* argument) {
        if (value <= 0.0) {
            throw OpalException(where, std::string("\"") + argument + "\" must be greater than 0.");
        }
    }

    void ensureCosineRange(double cosine, const char* where, const char* argument) {
        if (cosine < -1.0 - 1.0e-12 || cosine > 1.0 + 1.0e-12) {
            throw OpalException(where, std::string("\"") + argument + "\" must lie in [-1, 1].");
        }
    }
}

namespace Physics {
namespace LinearCompton {

double photonEnergyFromWavelengthGeV(double wavelength_m) {
    constexpr const char* where = "Physics::LinearCompton::photonEnergyFromWavelengthGeV()";
    ensureStrictlyPositive(wavelength_m, where, "wavelength_m");
    return Physics::two_pi * Physics::h_bar * Physics::c / wavelength_m;
}

double electronGamma(double electronTotalEnergyGeV) {
    constexpr const char* where = "Physics::LinearCompton::electronGamma()";
    if (electronTotalEnergyGeV <= Physics::m_e) {
        throw OpalException(where, "Electron total energy must be larger than the electron rest energy.");
    }

    return electronTotalEnergyGeV / Physics::m_e;
}

double electronBeta(double electronTotalEnergyGeV) {
    const double gamma = electronGamma(electronTotalEnergyGeV);
    return std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

double restFrameIncomingPhotonEnergyGeV(double electronTotalEnergyGeV,
                                        double laserPhotonEnergyGeV,
                                        const Vector_t<double, 3>& beamDirection,
                                        const Vector_t<double, 3>& laserDirection) {
    constexpr const char* where = "Physics::LinearCompton::restFrameIncomingPhotonEnergyGeV()";
    ensureStrictlyPositive(laserPhotonEnergyGeV, where, "laserPhotonEnergyGeV");

    const double gamma = electronGamma(electronTotalEnergyGeV);
    const double beta = electronBeta(electronTotalEnergyGeV);
    const Vector_t<double, 3> normalizedBeamDirection =
        normalizeDirection(beamDirection, where, "beamDirection");
    const Vector_t<double, 3> normalizedLaserDirection =
        normalizeDirection(laserDirection, where, "laserDirection");
    const double cosAlpha = clampCosine(dot(normalizedBeamDirection, normalizedLaserDirection));

    return gamma * laserPhotonEnergyGeV * (1.0 - beta * cosAlpha);
}

double invariantKappa(double incomingPhotonEnergyERFGeV) {
    constexpr const char* where = "Physics::LinearCompton::invariantKappa()";
    ensureStrictlyPositive(incomingPhotonEnergyERFGeV, where, "incomingPhotonEnergyERFGeV");
    return incomingPhotonEnergyERFGeV / Physics::m_e;
}

double invariantX(double incomingPhotonEnergyERFGeV) {
    return 2.0 * invariantKappa(incomingPhotonEnergyERFGeV);
}

double invariantX(double electronTotalEnergyGeV,
                  double laserPhotonEnergyGeV,
                  const Vector_t<double, 3>& beamDirection,
                  const Vector_t<double, 3>& laserDirection) {
    return invariantX(restFrameIncomingPhotonEnergyGeV(electronTotalEnergyGeV,
                                                       laserPhotonEnergyGeV,
                                                       beamDirection,
                                                       laserDirection));
}

double scatteredPhotonEnergyMinERFGeV(double incomingPhotonEnergyERFGeV) {
    const double kappa = invariantKappa(incomingPhotonEnergyERFGeV);
    return incomingPhotonEnergyERFGeV / (1.0 + 2.0 * kappa);
}

double scatteredPhotonEnergyMaxERFGeV(double incomingPhotonEnergyERFGeV) {
    invariantKappa(incomingPhotonEnergyERFGeV);
    return incomingPhotonEnergyERFGeV;
}

double scatteredPhotonEnergyERFGeV(double incomingPhotonEnergyERFGeV,
                                   double scatteringCosineERF) {
    constexpr const char* where = "Physics::LinearCompton::scatteredPhotonEnergyERFGeV()";
    const double kappa = invariantKappa(incomingPhotonEnergyERFGeV);
    ensureCosineRange(scatteringCosineERF, where, "scatteringCosineERF");
    const double clampedCosine = clampCosine(scatteringCosineERF);
    return incomingPhotonEnergyERFGeV / (1.0 + kappa * (1.0 - clampedCosine));
}

double differentialCrossSectionOmegaERF(double incomingPhotonEnergyERFGeV,
                                        double scatteredPhotonEnergyERFGeV) {
    constexpr const char* where = "Physics::LinearCompton::differentialCrossSectionOmegaERF()";
    const double kappa = invariantKappa(incomingPhotonEnergyERFGeV);
    const double minEnergy = scatteredPhotonEnergyMinERFGeV(incomingPhotonEnergyERFGeV);
    const double maxEnergy = scatteredPhotonEnergyMaxERFGeV(incomingPhotonEnergyERFGeV);
    if (scatteredPhotonEnergyERFGeV < minEnergy - 1.0e-15 ||
        scatteredPhotonEnergyERFGeV > maxEnergy + 1.0e-15) {
        throw OpalException(where, "\"scatteredPhotonEnergyERFGeV\" is outside the physical Compton support.");
    }

    const double clampedEnergy = std::max(minEnergy, std::min(maxEnergy, scatteredPhotonEnergyERFGeV));
    const double energyRatio = clampedEnergy / incomingPhotonEnergyERFGeV;
    const double cosine = 1.0 - Physics::m_e * (1.0 / clampedEnergy - 1.0 / incomingPhotonEnergyERFGeV);
    const double clampedCosine = clampCosine(cosine);
    const double sin2 = std::max(0.0, 1.0 - clampedCosine * clampedCosine);
    const double dsigma_dOmega = 0.5 * Physics::r_e * Physics::r_e * energyRatio * energyRatio
        * (1.0 / energyRatio + energyRatio - sin2);
    const double jacobian = Physics::m_e / (clampedEnergy * clampedEnergy);
    (void)kappa;
    return Physics::two_pi * dsigma_dOmega * jacobian;
}

double totalCrossSection(double incomingPhotonEnergyERFGeV) {
    const double kappa = invariantKappa(incomingPhotonEnergyERFGeV);

    // The exact closed form suffers from cancellation for kappa << 1.
    // Use the low-energy expansion around the Thomson limit there.
    if (kappa < 1.0e-3) {
        return thomsonCrossSection() * (1.0 - 2.0 * kappa + 5.2 * kappa * kappa);
    }

    const double logTerm = std::log(1.0 + 2.0 * kappa);
    const double prefactor = 2.0 * Physics::pi * Physics::r_e * Physics::r_e;
    const double bracket = ((1.0 + kappa) / (kappa * kappa))
            * ((2.0 * (1.0 + kappa)) / (1.0 + 2.0 * kappa) - logTerm / kappa)
        + logTerm / (2.0 * kappa)
        - (1.0 + 3.0 * kappa) / ((1.0 + 2.0 * kappa) * (1.0 + 2.0 * kappa));
    return prefactor * bracket;
}

double thomsonCrossSection() {
    return (8.0 / 3.0) * Physics::pi * Physics::r_e * Physics::r_e;
}

double restFrameScatteringCosineForLabForwardPhoton(double electronTotalEnergyGeV,
                                                    const Vector_t<double, 3>& beamDirection,
                                                    const Vector_t<double, 3>& laserDirection) {
    constexpr const char* where = "Physics::LinearCompton::restFrameScatteringCosineForLabForwardPhoton()";
    const double beta = electronBeta(electronTotalEnergyGeV);
    const Vector_t<double, 3> normalizedBeamDirection =
        normalizeDirection(beamDirection, where, "beamDirection");
    const Vector_t<double, 3> normalizedLaserDirection =
        normalizeDirection(laserDirection, where, "laserDirection");
    const double cosAlpha = clampCosine(dot(normalizedBeamDirection, normalizedLaserDirection));
    return clampCosine((cosAlpha - beta) / (1.0 - beta * cosAlpha));
}

double labForwardPhotonEnergyGeV(double electronTotalEnergyGeV,
                                 double laserPhotonEnergyGeV,
                                 const Vector_t<double, 3>& beamDirection,
                                 const Vector_t<double, 3>& laserDirection) {
    const double gamma = electronGamma(electronTotalEnergyGeV);
    const double beta = electronBeta(electronTotalEnergyGeV);
    const double omega1 = restFrameIncomingPhotonEnergyGeV(electronTotalEnergyGeV,
                                                           laserPhotonEnergyGeV,
                                                           beamDirection,
                                                           laserDirection);
    const double scatteringCosine = restFrameScatteringCosineForLabForwardPhoton(electronTotalEnergyGeV,
                                                                                  beamDirection,
                                                                                  laserDirection);
    const double omega2 = scatteredPhotonEnergyERFGeV(omega1, scatteringCosine);
    return gamma * (1.0 + beta) * omega2;
}

}  // namespace LinearCompton
}  // namespace Physics
