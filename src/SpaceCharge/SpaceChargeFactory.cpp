/** @file SpaceChargeFactory.cpp @brief Constructs the selected space-charge algorithm. */

#include "SpaceCharge/SpaceChargeFactory.h"

#include "PartBunch/PartBunch.h"
#include "SpaceCharge/CartesianPIC/CartesianPICAlgorithm.h"
#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"
#include "Utilities/OpalException.h"

#include <type_traits>
#include <utility>
#include <vector>

namespace opalx::spacecharge {

    std::unique_ptr<SpaceChargeSolver> makeSpaceChargeSolver(
            SpaceChargeConfig config, PartBunch_t& bunch, DataSink* dataSink) {
        validateSpaceChargeConfig(config);

        std::vector<ParticleContainer<double, 3>*> particles;
        particles.reserve(bunch.getNumParticleContainers());
        for (const auto& container : bunch.getParticleContainers()) {
            if (container == nullptr) {
                throw OpalException(
                        "makeSpaceChargeSolver", "Cannot use a null particle container.");
            }
            particles.push_back(container.get());
        }
        std::shared_ptr<const BunchStateHandler> bunchState = bunch.getBunchStateHandler();

        std::unique_ptr<SpaceChargeAlgorithm> algorithm = std::visit(
                [&](auto selected) -> std::unique_ptr<SpaceChargeAlgorithm> {
                    using Config = std::decay_t<decltype(selected)>;
                    if constexpr (std::is_same_v<Config, CartesianPICConfig>) {
                        if (selected.backend == PoissonSolverType::ConjugateGradient) {
                            throw OpalException(
                                    "makeSpaceChargeSolver",
                                    "The CG Poisson backend is recognized but not implemented.");
                        }
                        return std::make_unique<CartesianPICAlgorithm>(
                                std::move(selected), particles,
                                std::make_unique<CartesianPICFieldStorage<double, 3>>(
                                        bunch.cartesianDomain()),
                                dataSink, bunchState);
                    } else {
                        return std::make_unique<FFT2D5Algorithm>(
                                std::move(selected), particles, bunchState);
                    }
                },
                std::move(config));
        return std::make_unique<SpaceChargeSolver>(std::move(algorithm), particles.size());
    }

}  // namespace opalx::spacecharge
