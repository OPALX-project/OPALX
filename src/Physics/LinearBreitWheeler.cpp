#include "Physics/LinearBreitWheeler.h"

#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
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

    void ensureStrictlyPositive(double value, const char* where, const char* argument) {
        if (value <= 0.0) {
            throw OpalException(where, std::string("\"") + argument + "\" must be greater than 0.");
        }
    }

    double clampCosine(double value) {
        return std::max(-1.0, std::min(1.0, value));
    }

    std::uint64_t deterministicHostSeed() {
        return Options::seed >= 0 ? static_cast<std::uint64_t>(Options::seed) : 123456789ULL;
    }

    std::uint64_t mixSeed(std::uint64_t seed, std::uint64_t streamIndex) {
        const std::uint64_t golden = 0x9E3779B97F4A7C15ULL;
        std::uint64_t mixed = seed + golden * (streamIndex + 1ULL);
        mixed ^= (mixed >> 30);
        mixed *= 0xBF58476D1CE4E5B9ULL;
        mixed ^= (mixed >> 27);
        mixed *= 0x94D049BB133111EBULL;
        mixed ^= (mixed >> 31);
        return mixed;
    }

    Vector_t<double, 3> chooseReferenceDirection(const Vector_t<double, 3>& direction) {
        Vector_t<double, 3> reference(0.0);
        reference(0) = 1.0;
        if (std::abs(direction(0)) > 0.9) {
            reference(0) = 0.0;
            reference(1) = 1.0;
        }
        return reference;
    }

    void orthonormalFrame(const Vector_t<double, 3>& axis,
                          Vector_t<double, 3>& e1,
                          Vector_t<double, 3>& e2,
                          const char* where) {
        e1 = cross(axis, chooseReferenceDirection(axis));
        e1 = normalizeDirection(e1, where, "e1");
        e2 = normalizeDirection(cross(axis, e1), where, "e2");
    }

    double cainProposalY(double invariantSGeV2, double z) {
        const double cz = 2.0 * std::log(invariantSGeV2 / (4.0 * Physics::m_e * Physics::m_e));
        const double expCZ = std::exp(cz);
        const double expCZZ = std::exp(cz * z);
        return (1.0 + (1.0 - expCZZ) / (expCZ - 1.0)) / (1.0 + expCZZ);
    }

    double cainProposalJacobian(double invariantSGeV2, double z, double y) {
        const double cz = 2.0 * std::log(invariantSGeV2 / (4.0 * Physics::m_e * Physics::m_e));
        const double expCZ = std::exp(cz);
        const double expCZZ = std::exp(cz * z);
        return -cz * expCZZ / ((expCZ - 1.0) * (1.0 + expCZZ))
            - y * cz * expCZZ / (1.0 + expCZZ);
    }

    double cainUnpolarizedWeight(double invariantSGeV2, double z) {
        const double threshold = 4.0 * Physics::m_e * Physics::m_e;
        if (invariantSGeV2 <= threshold) {
            return 0.0;
        }

        const double w = 0.5 * std::sqrt(invariantSGeV2);
        const double beta = std::sqrt(std::max(0.0, 1.0 - threshold / invariantSGeV2));
        const double root = std::sqrt(std::max(0.0, w * w - Physics::m_e * Physics::m_e));
        const double x1 = Physics::m_e * Physics::m_e / (w * w + w * root);
        const double y = cainProposalY(invariantSGeV2, z);
        const double sm = 2.0 * Physics::m_e * Physics::m_e / ((1.0 + beta) * x1)
            * (beta * (1.0 - 2.0 * y) - 1.0);
        const double um = -2.0 * Physics::m_e * Physics::m_e / ((1.0 + beta) * x1)
            * (beta * (1.0 - 2.0 * y) + 1.0);
        const double costh1 = 1.0 + 2.0 * Physics::m_e * Physics::m_e * (1.0 / um + 1.0 / sm);
        const double d = -(sm / um + um / sm);
        const double f0 = d - (1.0 - costh1 * costh1);
        const double jacobian = cainProposalJacobian(invariantSGeV2, z, y);
        const double weight = f0 * beta * (1.0 + beta) * x1 * jacobian / 4.0;
        return std::max(0.0, weight);
    }

    double rejectionUpperBound(double invariantSGeV2) {
        constexpr int samples = 4096;
        double upper = 0.0;
        for (int i = 0; i <= samples; ++i) {
            const double z = -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(samples);
            upper = std::max(upper, cainUnpolarizedWeight(invariantSGeV2, z));
        }
        return std::max(upper * 1.0001, std::numeric_limits<double>::min());
    }

    Vector_t<double, 3> centerOfMomentumVelocity(double energy1GeV,
                                                 double energy2GeV,
                                                 const Vector_t<double, 3>& dir1,
                                                 const Vector_t<double, 3>& dir2) {
        return (energy1GeV * dir1 + energy2GeV * dir2) / (energy1GeV + energy2GeV);
    }

    Vector_t<double, 3> boostMomentum(const Vector_t<double, 3>& momentumCM,
                                      double energyCMGeV,
                                      const Vector_t<double, 3>& betaCM,
                                      double gammaCM,
                                      const char* where) {
        const double beta2 = dot(betaCM, betaCM);
        if (beta2 <= 0.0) {
            return momentumCM;
        }
        if (beta2 >= 1.0) {
            throw OpalException(where, "CM boost speed must stay below c.");
        }

        const double betaDotP = dot(betaCM, momentumCM);
        const double factor = ((gammaCM - 1.0) * betaDotP / beta2) + gammaCM * energyCMGeV;
        return momentumCM + factor * betaCM;
    }
}

namespace Physics {
namespace LinearBreitWheeler {

double photonEnergyFromWavelengthGeV(double wavelength_m) {
    constexpr const char* where = "Physics::LinearBreitWheeler::photonEnergyFromWavelengthGeV()";
    ensureStrictlyPositive(wavelength_m, where, "wavelength_m");
    return Physics::two_pi * Physics::h_bar * Physics::c / wavelength_m;
}

double invariantSGeV2(double photon1EnergyGeV,
                      double photon2EnergyGeV,
                      const Vector_t<double, 3>& photon1Direction,
                      const Vector_t<double, 3>& photon2Direction) {
    constexpr const char* where = "Physics::LinearBreitWheeler::invariantSGeV2()";
    ensureStrictlyPositive(photon1EnergyGeV, where, "photon1EnergyGeV");
    ensureStrictlyPositive(photon2EnergyGeV, where, "photon2EnergyGeV");
    const auto dir1 = normalizeDirection(photon1Direction, where, "photon1Direction");
    const auto dir2 = normalizeDirection(photon2Direction, where, "photon2Direction");
    return 2.0 * photon1EnergyGeV * photon2EnergyGeV * (1.0 - clampCosine(dot(dir1, dir2)));
}

double thresholdInvariantSGeV2() {
    return 4.0 * Physics::m_e * Physics::m_e;
}

bool isAboveThreshold(double sGeV2) {
    return sGeV2 >= thresholdInvariantSGeV2();
}

double pairBetaCM(double sGeV2) {
    constexpr const char* where = "Physics::LinearBreitWheeler::pairBetaCM()";
    ensureStrictlyPositive(sGeV2, where, "sGeV2");
    if (!isAboveThreshold(sGeV2)) {
        return 0.0;
    }
    return std::sqrt(std::max(0.0, 1.0 - thresholdInvariantSGeV2() / sGeV2));
}

double totalCrossSection(double sGeV2) {
    constexpr const char* where = "Physics::LinearBreitWheeler::totalCrossSection()";
    ensureStrictlyPositive(sGeV2, where, "sGeV2");
    if (!isAboveThreshold(sGeV2)) {
        return 0.0;
    }

    const double beta = pairBetaCM(sGeV2);
    const double oneMinusBeta2 = std::max(0.0, 1.0 - beta * beta);
    const double logTerm = std::log((1.0 + beta) / (1.0 - beta));
    return 0.5 * Physics::pi * Physics::r_e * Physics::r_e * oneMinusBeta2
        * ((3.0 - std::pow(beta, 4)) * logTerm - 2.0 * beta * (2.0 - beta * beta));
}

SamplingKernel makeSamplingKernel(double highEnergyPhotonEnergyGeV,
                                  double laserPhotonEnergyGeV,
                                  const Vector_t<double, 3>& highEnergyDirection,
                                  const Vector_t<double, 3>& laserDirection) {
    constexpr const char* where = "Physics::LinearBreitWheeler::makeSamplingKernel()";
    ensureStrictlyPositive(highEnergyPhotonEnergyGeV, where, "highEnergyPhotonEnergyGeV");
    ensureStrictlyPositive(laserPhotonEnergyGeV, where, "laserPhotonEnergyGeV");

    SamplingKernel kernel;
    kernel.highEnergyPhotonEnergyGeV = highEnergyPhotonEnergyGeV;
    kernel.laserPhotonEnergyGeV = laserPhotonEnergyGeV;
    kernel.highEnergyDirection = normalizeDirection(highEnergyDirection, where, "highEnergyDirection");
    kernel.laserDirection = normalizeDirection(laserDirection, where, "laserDirection");
    kernel.invariantSGeV2 = invariantSGeV2(highEnergyPhotonEnergyGeV,
                                           laserPhotonEnergyGeV,
                                           kernel.highEnergyDirection,
                                           kernel.laserDirection);
    if (!isAboveThreshold(kernel.invariantSGeV2)) {
        throw OpalException(where, "Photon-photon invariant is below the Breit-Wheeler threshold.");
    }

    kernel.pairBetaCM = pairBetaCM(kernel.invariantSGeV2);
    kernel.cmVelocity = centerOfMomentumVelocity(highEnergyPhotonEnergyGeV,
                                                 laserPhotonEnergyGeV,
                                                 kernel.highEnergyDirection,
                                                 kernel.laserDirection);
    const double beta2 = dot(kernel.cmVelocity, kernel.cmVelocity);
    if (beta2 >= 1.0) {
        throw OpalException(where, "Center-of-momentum boost speed must stay below c.");
    }
    kernel.cmLorentzGamma = 1.0 / std::sqrt(std::max(1.0e-30, 1.0 - beta2));
    kernel.rejectionUpperBound = rejectionUpperBound(kernel.invariantSGeV2);
    return kernel;
}

std::mt19937_64 makeHostRandomEngine(std::uint64_t streamIndex) {
    return std::mt19937_64(mixSeed(deterministicHostSeed(), streamIndex));
}

SampledEvent sampleEvent(const SamplingKernel& kernel, std::mt19937_64& engine) {
    constexpr const char* where = "Physics::LinearBreitWheeler::sampleEvent()";
    if (!isAboveThreshold(kernel.invariantSGeV2)) {
        throw OpalException(where, "Sampling kernel is below the Breit-Wheeler threshold.");
    }

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const double sqrtS = std::sqrt(kernel.invariantSGeV2);
    const double energyCMGeV = 0.5 * sqrtS;
    const double momentumCM = std::sqrt(std::max(0.0, energyCMGeV * energyCMGeV - Physics::m_e * Physics::m_e));

    double z = 0.0;
    double y = 0.0;
    double scatteringCosine = 0.0;
    while (true) {
        z = 1.0 - 2.0 * unit(engine);
        y = cainProposalY(kernel.invariantSGeV2, z);
        scatteringCosine = 1.0 - 2.0 * y;
        const double accept = cainUnpolarizedWeight(kernel.invariantSGeV2, z) / kernel.rejectionUpperBound;
        if (unit(engine) <= accept) {
            break;
        }
    }

    const double azimuth = Physics::two_pi * unit(engine);
    Vector_t<double, 3> cmAxis1(0.0);
    Vector_t<double, 3> cmAxis2(0.0);
    const Vector_t<double, 3> highEnergyCM = normalizeDirection(
        kernel.highEnergyDirection - kernel.cmVelocity,
        where,
        "highEnergyCM");
    orthonormalFrame(highEnergyCM, cmAxis1, cmAxis2, where);
    const double sinTheta = std::sqrt(std::max(0.0, 1.0 - scatteringCosine * scatteringCosine));
    const Vector_t<double, 3> electronDirectionCM = normalizeDirection(
        scatteringCosine * highEnergyCM
            + sinTheta * (std::cos(azimuth) * cmAxis1 + std::sin(azimuth) * cmAxis2),
        where,
        "electronDirectionCM");

    const Vector_t<double, 3> electronMomentumCM = momentumCM * electronDirectionCM;
    const Vector_t<double, 3> positronMomentumCM = -electronMomentumCM;
    const Vector_t<double, 3> electronMomentumLab = boostMomentum(electronMomentumCM,
                                                                  energyCMGeV,
                                                                  kernel.cmVelocity,
                                                                  kernel.cmLorentzGamma,
                                                                  where);
    const Vector_t<double, 3> positronMomentumLab = boostMomentum(positronMomentumCM,
                                                                  energyCMGeV,
                                                                  kernel.cmVelocity,
                                                                  kernel.cmLorentzGamma,
                                                                  where);

    SampledEvent event;
    event.scatteringCosineCM = scatteringCosine;
    event.azimuthCM = azimuth;
    event.electronMomentumLabGeV = electronMomentumLab;
    event.positronMomentumLabGeV = positronMomentumLab;
    event.electronEnergyLabGeV = std::sqrt(dot(electronMomentumLab, electronMomentumLab) + Physics::m_e * Physics::m_e);
    event.positronEnergyLabGeV = std::sqrt(dot(positronMomentumLab, positronMomentumLab) + Physics::m_e * Physics::m_e);
    return event;
}

}  // namespace LinearBreitWheeler
}  // namespace Physics
