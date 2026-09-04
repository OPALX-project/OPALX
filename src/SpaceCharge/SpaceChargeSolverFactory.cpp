/**
 * @file SpaceChargeSolverFactory.cpp
 * @brief Implements construction-time space-charge algorithm selection.
 */

#include "SpaceCharge/SpaceChargeSolverFactory.h"

#include "PartBunch/PartBunch.h"
#include "SpaceCharge/CartesianPIC/CartesianPICAlgorithm.h"
#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"
#include "Utilities/OpalException.h"

#include <memory>
#include <utility>
#include <vector>

namespace opalx::spacecharge {

    std::unique_ptr<SpaceChargeSolver> SpaceChargeSolverFactory::create(
            SpaceChargeConfig config, PartBunch_t& bunch, DataSink* dataSink) {
        if (config.algorithmType() == SpaceChargeAlgorithmType::CartesianPIC
            && config.get<CartesianPICConfig>().backend() == PoissonSolverType::ConjugateGradient) {
            throw OpalException(
                    "SpaceChargeSolverFactory::create",
                    "The CG Poisson backend is recognized but not implemented.");
        }

        std::vector<ParticleFieldBinding3D> bindings;
        bindings.reserve(bunch.getNumParticleContainers());
        for (const auto& container : bunch.getParticleContainers()) {
            if (container == nullptr) {
                throw OpalException(
                        "SpaceChargeSolverFactory::create",
                        "Cannot construct a space-charge solver for a null particle container.");
            }
            bindings.push_back(makeParticleFieldBinding(*container));
        }

        std::unique_ptr<SpaceChargeAlgorithm> algorithm;
        switch (config.algorithmType()) {
            case SpaceChargeAlgorithmType::CartesianPIC:
                algorithm = std::make_unique<CartesianPICAlgorithm>(
                        config.get<CartesianPICConfig>(), bindings,
                        std::make_unique<CartesianPICFieldStorage<double, 3>>(
                                bunch.cartesianDomain()),
                        dataSink);
                break;
            case SpaceChargeAlgorithmType::FFT2D5:
                algorithm = std::make_unique<FFT2D5Algorithm>(config.get<FFT2D5Config>(), bindings);
                break;
        }

        if (algorithm == nullptr) {
            throw OpalException(
                    "SpaceChargeSolverFactory::create",
                    "No implementation exists for this algorithm.");
        }
        return std::make_unique<SpaceChargeSolver>(
                std::move(config), std::move(algorithm), std::move(bindings));
    }

}  // namespace opalx::spacecharge
