/**
 * @file BeamMoments.cpp
 * @brief Checked accessors and derived quantities of BeamMoments.
 */

#include "Moments/BeamMoments.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

namespace opalx::moments {

    namespace {
        std::string describe(Request groups) {
            static constexpr std::pair<Moment, const char*> names[] = {
                    {Moment::Bounds, "Bounds"},
                    {Moment::MeanPosition, "MeanPosition"},
                    {Moment::MeanMomentum, "MeanMomentum"},
                    {Moment::Spread, "Spread"},
                    {Moment::Polarization, "Polarization"}};
            std::string text;
            for (const auto& [moment, name] : names) {
                if (groups.contains(moment)) {
                    text += text.empty() ? name : std::string(", ") + name;
                }
            }
            return text;
        }
    }  // namespace

    void BeamMoments::require(Request groups) const {
        if (!available_m.contains(groups)) {
            throw OpalException(
                    "BeamMoments", "The moment group(s) " + describe(groups.without(available_m))
                                           + " were not requested for this snapshot.");
        }
    }

    const BeamMoments::Vector3& BeamMoments::minR() const {
        require(Moment::Bounds);
        return minR_m;
    }

    const BeamMoments::Vector3& BeamMoments::maxR() const {
        require(Moment::Bounds);
        return maxR_m;
    }

    const BeamMoments::Vector3& BeamMoments::meanR() const {
        require(Moment::MeanPosition);
        return meanR_m;
    }

    const BeamMoments::Vector3& BeamMoments::meanP() const {
        require(Moment::MeanMomentum);
        return meanP_m;
    }

    double BeamMoments::meanGamma() const {
        require(Moment::MeanMomentum);
        return meanGamma_m;
    }

    double BeamMoments::meanGammaZ() const {
        require(Moment::MeanMomentum);
        return std::sqrt(meanP_m[2] * meanP_m[2] + 1.0);
    }

    double BeamMoments::meanKineticEnergy() const {
        require(Moment::MeanMomentum);
        return meanKineticEnergy_m;
    }

    double BeamMoments::covariance(unsigned i, unsigned j) const {
        require(Moment::Spread);
        if (i >= 6 || j >= 6) {
            throw OpalException("BeamMoments::covariance", "Phase-space index out of range.");
        }
        return covariance_m[i][j];
    }

    BeamMoments::Vector3 BeamMoments::rmsR() const {
        require(Moment::Spread);
        Vector3 rms;
        for (unsigned d = 0; d < 3; ++d) {
            rms[d] = std::sqrt(covariance_m[2 * d][2 * d]);
        }
        return rms;
    }

    BeamMoments::Vector3 BeamMoments::rmsP() const {
        require(Moment::Spread);
        Vector3 rms;
        for (unsigned d = 0; d < 3; ++d) {
            rms[d] = std::sqrt(covariance_m[2 * d + 1][2 * d + 1]);
        }
        return rms;
    }

    BeamMoments::Vector3 BeamMoments::correlationRP() const {
        const Vector3 sigmaR = rmsR();
        const Vector3 sigmaP = rmsP();
        Vector3 correlation;
        for (unsigned d = 0; d < 3; ++d) {
            const double spread = sigmaR[d] * sigmaP[d];
            correlation[d]      = spread > 0.0 ? covariance_m[2 * d][2 * d + 1] / spread : 0.0;
        }
        return correlation;
    }

    BeamMoments::Vector3 BeamMoments::normalizedEmittance() const {
        const Vector3 sigmaR = rmsR();
        const Vector3 sigmaP = rmsP();
        Vector3 emittance;
        for (unsigned d = 0; d < 3; ++d) {
            const double squared = std::pow(sigmaR[d] * sigmaP[d], 2)
                                   - std::pow(covariance_m[2 * d][2 * d + 1], 2);
            emittance[d] = std::sqrt(std::max(squared, 0.0));
        }
        return emittance;
    }

    BeamMoments::Vector3 BeamMoments::geometricEmittance() const {
        require(Moment::MeanMomentum);
        const double betaGamma = std::sqrt(std::max(std::pow(meanGamma_m, 2) - 1.0, 0.0));
        return betaGamma > 0.0 ? normalizedEmittance() / Vector3(betaGamma) : Vector3(0.0);
    }

    double BeamMoments::stdKineticEnergy() const {
        require(Moment::Spread);
        return stdKineticEnergy_m;
    }

    const BeamMoments::Vector3& BeamMoments::restFrameTemperature() const {
        require(Moment::Spread);
        return restFrameTemperature_m;
    }

    const BeamMoments::Vector3& BeamMoments::meanPolarization() const {
        require(Moment::Polarization);
        return meanPolarization_m;
    }

    const BeamMoments::Vector3& BeamMoments::rmsPolarization() const {
        require(Moment::Polarization);
        return rmsPolarization_m;
    }

    double BeamMoments::meanPolarizationMagnitude() const {
        const Vector3& mean = meanPolarization();
        return std::sqrt(mean[0] * mean[0] + mean[1] * mean[1] + mean[2] * mean[2]);
    }

}  // namespace opalx::moments
