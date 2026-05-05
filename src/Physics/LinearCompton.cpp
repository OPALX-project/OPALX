#include "Physics/LinearCompton.h"

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
    Vector_t<double, 3> normalizeDirection(
            const Vector_t<double, 3>& direction, const char* where, const char* argument) {
        const double norm = std::sqrt(dot(direction, direction));
        if (norm <= 0.0) {
            throw OpalException(
                    where, std::string("\"") + argument + "\" must be a non-zero vector.");
        }

        return direction / norm;
    }

    double clampCosine(double value) { return std::max(-1.0, std::min(1.0, value)); }

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

    Vector_t<double, 3> chooseReferenceDirection(const Vector_t<double, 3>& direction) {
        Vector_t<double, 3> reference(0.0);
        reference(0) = 1.0;
        if (std::abs(direction(0)) > 0.9) {
            reference(0) = 0.0;
            reference(1) = 1.0;
        }
        return reference;
    }

    Vector_t<double, 3> scatteringPlaneAxis(
            const Vector_t<double, 3>& incomingDirectionERF,
            const Vector_t<double, 3>& boostDirectionLab, const char* where) {
        Vector_t<double, 3> axis =
                boostDirectionLab
                - dot(boostDirectionLab, incomingDirectionERF) * incomingDirectionERF;
        const double axisNorm = std::sqrt(dot(axis, axis));
        if (axisNorm > 1.0e-14) {
            return axis / axisNorm;
        }

        axis = cross(incomingDirectionERF, chooseReferenceDirection(incomingDirectionERF));
        return normalizeDirection(axis, where, "scatteringPlaneAxis");
    }

    Vector_t<double, 3> scatteredDirectionERF(
            const Vector_t<double, 3>& incomingDirectionERF,
            const Vector_t<double, 3>& boostDirectionLab, double scatteringCosineERF,
            double azimuthERF, const char* where) {
        ensureCosineRange(scatteringCosineERF, where, "scatteringCosineERF");
        const double clampedCosine = clampCosine(scatteringCosineERF);
        const double sinTheta      = std::sqrt(std::max(0.0, 1.0 - clampedCosine * clampedCosine));
        const Vector_t<double, 3> axis1 =
                scatteringPlaneAxis(incomingDirectionERF, boostDirectionLab, where);
        const Vector_t<double, 3> axis2 =
                normalizeDirection(cross(incomingDirectionERF, axis1), where, "axis2");
        return normalizeDirection(
                clampedCosine * incomingDirectionERF
                        + sinTheta * (std::cos(azimuthERF) * axis1 + std::sin(azimuthERF) * axis2),
                where, "scatteredDirectionERF");
    }

    std::uint64_t deterministicHostSeed() {
        return Options::seed >= 0 ? static_cast<std::uint64_t>(Options::seed) : 123456789ULL;
    }

    std::uint64_t mixSeed(std::uint64_t seed, std::uint64_t streamIndex) {
        const std::uint64_t golden = 0x9E3779B97F4A7C15ULL;
        std::uint64_t mixed        = seed + golden * (streamIndex + 1ULL);
        mixed ^= (mixed >> 30);
        mixed *= 0xBF58476D1CE4E5B9ULL;
        mixed ^= (mixed >> 27);
        mixed *= 0x94D049BB133111EBULL;
        mixed ^= (mixed >> 31);
        return mixed;
    }

    double rejectionUpperBoundSolidAngleERF(double incomingPhotonEnergyERFGeV) {
        constexpr int samples = 4096;
        double upperBound     = 0.0;
        for (int i = 0; i <= samples; ++i) {
            const double mu = -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(samples);
            upperBound      = std::max(
                    upperBound, Physics::LinearCompton::differentialCrossSectionSolidAngleERF(
                                        incomingPhotonEnergyERFGeV, mu));
        }

        return std::max(upperBound * 1.0001, std::numeric_limits<double>::min());
    }
}  // namespace

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
                throw OpalException(
                        where,
                        "Electron total energy must be larger than the electron rest energy.");
            }

            return electronTotalEnergyGeV / Physics::m_e;
        }

        double electronBeta(double electronTotalEnergyGeV) {
            const double gamma = electronGamma(electronTotalEnergyGeV);
            return std::sqrt(1.0 - 1.0 / (gamma * gamma));
        }

        double restFrameIncomingPhotonEnergyGeV(
                double electronTotalEnergyGeV, double laserPhotonEnergyGeV,
                const Vector_t<double, 3>& beamDirection,
                const Vector_t<double, 3>& laserDirection) {
            constexpr const char* where =
                    "Physics::LinearCompton::restFrameIncomingPhotonEnergyGeV()";
            ensureStrictlyPositive(laserPhotonEnergyGeV, where, "laserPhotonEnergyGeV");

            const double gamma = electronGamma(electronTotalEnergyGeV);
            const double beta  = electronBeta(electronTotalEnergyGeV);
            const Vector_t<double, 3> normalizedBeamDirection =
                    normalizeDirection(beamDirection, where, "beamDirection");
            const Vector_t<double, 3> normalizedLaserDirection =
                    normalizeDirection(laserDirection, where, "laserDirection");
            const double cosAlpha =
                    clampCosine(dot(normalizedBeamDirection, normalizedLaserDirection));

            return gamma * laserPhotonEnergyGeV * (1.0 - beta * cosAlpha);
        }

        Vector_t<double, 3> restFrameIncomingPhotonDirection(
                double electronTotalEnergyGeV, const Vector_t<double, 3>& beamDirection,
                const Vector_t<double, 3>& laserDirection) {
            constexpr const char* where =
                    "Physics::LinearCompton::restFrameIncomingPhotonDirection()";
            const double gamma = electronGamma(electronTotalEnergyGeV);
            const double beta  = electronBeta(electronTotalEnergyGeV);
            const Vector_t<double, 3> normalizedBeamDirection =
                    normalizeDirection(beamDirection, where, "beamDirection");
            const Vector_t<double, 3> normalizedLaserDirection =
                    normalizeDirection(laserDirection, where, "laserDirection");
            const double cosAlpha =
                    clampCosine(dot(normalizedBeamDirection, normalizedLaserDirection));
            const double denominator = 1.0 - beta * cosAlpha;
            const Vector_t<double, 3> transverseComponent =
                    normalizedLaserDirection - cosAlpha * normalizedBeamDirection;
            const Vector_t<double, 3> directionERF =
                    transverseComponent / (gamma * denominator)
                    + ((cosAlpha - beta) / denominator) * normalizedBeamDirection;
            return normalizeDirection(directionERF, where, "directionERF");
        }

        Vector_t<double, 3> labPhotonDirection(
                double electronTotalEnergyGeV, const Vector_t<double, 3>& beamDirection,
                const Vector_t<double, 3>& laserDirection, double scatteringCosineERF,
                double azimuthERF) {
            constexpr const char* where = "Physics::LinearCompton::labPhotonDirection()";
            const double gamma          = electronGamma(electronTotalEnergyGeV);
            const double beta           = electronBeta(electronTotalEnergyGeV);
            const Vector_t<double, 3> normalizedBeamDirection =
                    normalizeDirection(beamDirection, where, "beamDirection");
            const Vector_t<double, 3> incomingDirectionERF = restFrameIncomingPhotonDirection(
                    electronTotalEnergyGeV, beamDirection, laserDirection);
            const Vector_t<double, 3> outgoingDirectionERF = scatteredDirectionERF(
                    incomingDirectionERF, normalizedBeamDirection, scatteringCosineERF, azimuthERF,
                    where);
            const double parallelCosine = dot(outgoingDirectionERF, normalizedBeamDirection);
            const Vector_t<double, 3> transverseDirection =
                    outgoingDirectionERF - parallelCosine * normalizedBeamDirection;
            const Vector_t<double, 3> boostedDirection =
                    transverseDirection + gamma * (parallelCosine + beta) * normalizedBeamDirection;
            return normalizeDirection(boostedDirection, where, "boostedDirection");
        }

        std::mt19937_64 makeHostRandomEngine(std::uint64_t streamIndex) {
            return std::mt19937_64(mixSeed(deterministicHostSeed(), streamIndex));
        }

        SamplingKernel makeSamplingKernel(
                double electronTotalEnergyGeV, double laserPhotonEnergyGeV,
                const Vector_t<double, 3>& beamDirection,
                const Vector_t<double, 3>& laserDirection) {
            constexpr const char* where = "Physics::LinearCompton::makeSamplingKernel()";
            ensureStrictlyPositive(laserPhotonEnergyGeV, where, "laserPhotonEnergyGeV");

            SamplingKernel kernel;
            kernel.electronTotalEnergyGeV = electronTotalEnergyGeV;
            kernel.laserPhotonEnergyGeV   = laserPhotonEnergyGeV;
            kernel.beamDirection  = normalizeDirection(beamDirection, where, "beamDirection");
            kernel.laserDirection = normalizeDirection(laserDirection, where, "laserDirection");
            kernel.incomingPhotonEnergyERFGeV = restFrameIncomingPhotonEnergyGeV(
                    electronTotalEnergyGeV, laserPhotonEnergyGeV, kernel.beamDirection,
                    kernel.laserDirection);
            kernel.rejectionUpperBoundSolidAngleERF =
                    rejectionUpperBoundSolidAngleERF(kernel.incomingPhotonEnergyERFGeV);
            return kernel;
        }

        SampledEvent sampleEvent(const SamplingKernel& kernel, std::mt19937_64& engine) {
            constexpr const char* where = "Physics::LinearCompton::sampleEvent()";
            ensureStrictlyPositive(
                    kernel.incomingPhotonEnergyERFGeV, where, "incomingPhotonEnergyERFGeV");
            ensureStrictlyPositive(
                    kernel.rejectionUpperBoundSolidAngleERF, where,
                    "rejectionUpperBoundSolidAngleERF");

            std::uniform_real_distribution<double> cosineDistribution(-1.0, 1.0);
            std::uniform_real_distribution<double> azimuthDistribution(0.0, Physics::two_pi);
            std::uniform_real_distribution<double> acceptDistribution(
                    0.0, kernel.rejectionUpperBoundSolidAngleERF);

            for (int attempt = 0; attempt < 1000000; ++attempt) {
                const double scatteringCosineERF      = cosineDistribution(engine);
                const double differentialCrossSection = differentialCrossSectionSolidAngleERF(
                        kernel.incomingPhotonEnergyERFGeV, scatteringCosineERF);
                if (acceptDistribution(engine) > differentialCrossSection) {
                    continue;
                }

                SampledEvent event;
                event.scatteringCosineERF         = scatteringCosineERF;
                event.azimuthERF                  = azimuthDistribution(engine);
                event.scatteredPhotonEnergyERFGeV = scatteredPhotonEnergyERFGeV(
                        kernel.incomingPhotonEnergyERFGeV, event.scatteringCosineERF);
                event.scatteredPhotonEnergyLabGeV = labPhotonEnergyGeV(
                        kernel.electronTotalEnergyGeV, kernel.laserPhotonEnergyGeV,
                        kernel.beamDirection, kernel.laserDirection, event.scatteringCosineERF,
                        event.azimuthERF);
                event.scatteredPhotonDirectionLab = labPhotonDirection(
                        kernel.electronTotalEnergyGeV, kernel.beamDirection, kernel.laserDirection,
                        event.scatteringCosineERF, event.azimuthERF);
                return event;
            }

            throw OpalException(
                    where, "Rejection sampler failed to accept an event within the attempt limit.");
        }

        double invariantKappa(double incomingPhotonEnergyERFGeV) {
            constexpr const char* where = "Physics::LinearCompton::invariantKappa()";
            ensureStrictlyPositive(incomingPhotonEnergyERFGeV, where, "incomingPhotonEnergyERFGeV");
            return incomingPhotonEnergyERFGeV / Physics::m_e;
        }

        double invariantX(double incomingPhotonEnergyERFGeV) {
            return 2.0 * invariantKappa(incomingPhotonEnergyERFGeV);
        }

        double invariantX(
                double electronTotalEnergyGeV, double laserPhotonEnergyGeV,
                const Vector_t<double, 3>& beamDirection,
                const Vector_t<double, 3>& laserDirection) {
            return invariantX(restFrameIncomingPhotonEnergyGeV(
                    electronTotalEnergyGeV, laserPhotonEnergyGeV, beamDirection, laserDirection));
        }

        double scatteredPhotonEnergyMinERFGeV(double incomingPhotonEnergyERFGeV) {
            const double kappa = invariantKappa(incomingPhotonEnergyERFGeV);
            return incomingPhotonEnergyERFGeV / (1.0 + 2.0 * kappa);
        }

        double scatteredPhotonEnergyMaxERFGeV(double incomingPhotonEnergyERFGeV) {
            invariantKappa(incomingPhotonEnergyERFGeV);
            return incomingPhotonEnergyERFGeV;
        }

        double scatteredPhotonEnergyERFGeV(
                double incomingPhotonEnergyERFGeV, double scatteringCosineERF) {
            constexpr const char* where = "Physics::LinearCompton::scatteredPhotonEnergyERFGeV()";
            const double kappa          = invariantKappa(incomingPhotonEnergyERFGeV);
            ensureCosineRange(scatteringCosineERF, where, "scatteringCosineERF");
            const double clampedCosine = clampCosine(scatteringCosineERF);
            return incomingPhotonEnergyERFGeV / (1.0 + kappa * (1.0 - clampedCosine));
        }

        double differentialCrossSectionSolidAngleERF(
                double incomingPhotonEnergyERFGeV, double scatteringCosineERF) {
            constexpr const char* where =
                    "Physics::LinearCompton::differentialCrossSectionSolidAngleERF()";
            invariantKappa(incomingPhotonEnergyERFGeV);
            ensureCosineRange(scatteringCosineERF, where, "scatteringCosineERF");
            const double clampedCosine = clampCosine(scatteringCosineERF);
            const double omega2 =
                    scatteredPhotonEnergyERFGeV(incomingPhotonEnergyERFGeV, clampedCosine);
            const double energyRatio = omega2 / incomingPhotonEnergyERFGeV;
            const double sin2        = std::max(0.0, 1.0 - clampedCosine * clampedCosine);
            return 0.5 * Physics::r_e * Physics::r_e * energyRatio * energyRatio
                   * (1.0 / energyRatio + energyRatio - sin2);
        }

        double differentialCrossSectionOmegaERF(
                double incomingPhotonEnergyERFGeV, double scatteredPhotonEnergyERFGeV) {
            constexpr const char* where =
                    "Physics::LinearCompton::differentialCrossSectionOmegaERF()";
            const double minEnergy = scatteredPhotonEnergyMinERFGeV(incomingPhotonEnergyERFGeV);
            const double maxEnergy = scatteredPhotonEnergyMaxERFGeV(incomingPhotonEnergyERFGeV);
            if (scatteredPhotonEnergyERFGeV < minEnergy - 1.0e-15
                || scatteredPhotonEnergyERFGeV > maxEnergy + 1.0e-15) {
                throw OpalException(
                        where,
                        "\"scatteredPhotonEnergyERFGeV\" is outside the physical Compton support.");
            }

            const double clampedEnergy =
                    std::max(minEnergy, std::min(maxEnergy, scatteredPhotonEnergyERFGeV));
            const double cosine =
                    1.0 - Physics::m_e * (1.0 / clampedEnergy - 1.0 / incomingPhotonEnergyERFGeV);
            const double clampedCosine = clampCosine(cosine);
            const double dsigma_dOmega = differentialCrossSectionSolidAngleERF(
                    incomingPhotonEnergyERFGeV, clampedCosine);
            const double jacobian = Physics::m_e / (clampedEnergy * clampedEnergy);
            return Physics::two_pi * dsigma_dOmega * jacobian;
        }

        double totalCrossSection(double incomingPhotonEnergyERFGeV) {
            const double kappa = invariantKappa(incomingPhotonEnergyERFGeV);

            // The exact closed form suffers from cancellation for kappa << 1.
            // Use the low-energy expansion around the Thomson limit there.
            if (kappa < 1.0e-3) {
                return thomsonCrossSection() * (1.0 - 2.0 * kappa + 5.2 * kappa * kappa);
            }

            const double logTerm   = std::log(1.0 + 2.0 * kappa);
            const double prefactor = 2.0 * Physics::pi * Physics::r_e * Physics::r_e;
            const double bracket =
                    ((1.0 + kappa) / (kappa * kappa))
                            * ((2.0 * (1.0 + kappa)) / (1.0 + 2.0 * kappa) - logTerm / kappa)
                    + logTerm / (2.0 * kappa)
                    - (1.0 + 3.0 * kappa) / ((1.0 + 2.0 * kappa) * (1.0 + 2.0 * kappa));
            return prefactor * bracket;
        }

        double thomsonCrossSection() {
            return (8.0 / 3.0) * Physics::pi * Physics::r_e * Physics::r_e;
        }

        double restFrameScatteringCosineForLabForwardPhoton(
                double electronTotalEnergyGeV, const Vector_t<double, 3>& beamDirection,
                const Vector_t<double, 3>& laserDirection) {
            constexpr const char* where =
                    "Physics::LinearCompton::restFrameScatteringCosineForLabForwardPhoton()";
            const double beta = electronBeta(electronTotalEnergyGeV);
            const Vector_t<double, 3> normalizedBeamDirection =
                    normalizeDirection(beamDirection, where, "beamDirection");
            const Vector_t<double, 3> normalizedLaserDirection =
                    normalizeDirection(laserDirection, where, "laserDirection");
            const double cosAlpha =
                    clampCosine(dot(normalizedBeamDirection, normalizedLaserDirection));
            return clampCosine((cosAlpha - beta) / (1.0 - beta * cosAlpha));
        }

        double labPhotonEnergyGeV(
                double electronTotalEnergyGeV, double laserPhotonEnergyGeV,
                const Vector_t<double, 3>& beamDirection, const Vector_t<double, 3>& laserDirection,
                double scatteringCosineERF, double azimuthERF) {
            constexpr const char* where = "Physics::LinearCompton::labPhotonEnergyGeV()";
            ensureStrictlyPositive(laserPhotonEnergyGeV, where, "laserPhotonEnergyGeV");
            ensureCosineRange(scatteringCosineERF, where, "scatteringCosineERF");

            const double gamma = electronGamma(electronTotalEnergyGeV);
            const double beta  = electronBeta(electronTotalEnergyGeV);
            const Vector_t<double, 3> normalizedBeamDirection =
                    normalizeDirection(beamDirection, where, "beamDirection");
            const Vector_t<double, 3> incomingDirectionERF = restFrameIncomingPhotonDirection(
                    electronTotalEnergyGeV, beamDirection, laserDirection);
            const double omega1 = restFrameIncomingPhotonEnergyGeV(
                    electronTotalEnergyGeV, laserPhotonEnergyGeV, beamDirection, laserDirection);
            const double clampedCosine                     = clampCosine(scatteringCosineERF);
            const Vector_t<double, 3> outgoingDirectionERF = scatteredDirectionERF(
                    incomingDirectionERF, normalizedBeamDirection, clampedCosine, azimuthERF,
                    where);
            const double omega2 = scatteredPhotonEnergyERFGeV(omega1, clampedCosine);
            return gamma * omega2
                   * (1.0 + beta * dot(outgoingDirectionERF, normalizedBeamDirection));
        }

        double labForwardPhotonEnergyGeV(
                double electronTotalEnergyGeV, double laserPhotonEnergyGeV,
                const Vector_t<double, 3>& beamDirection,
                const Vector_t<double, 3>& laserDirection) {
            const double scatteringCosine = restFrameScatteringCosineForLabForwardPhoton(
                    electronTotalEnergyGeV, beamDirection, laserDirection);
            return labPhotonEnergyGeV(
                    electronTotalEnergyGeV, laserPhotonEnergyGeV, beamDirection, laserDirection,
                    scatteringCosine, 0.0);
        }

    }  // namespace LinearCompton
}  // namespace Physics
