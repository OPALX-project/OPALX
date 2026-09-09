/** @file SpaceChargeFactory.cpp @brief Constructs the selected space-charge algorithm. */

#include "SpaceCharge/SpaceChargeFactory.h"

#include "PartBunch/PartBunch.h"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.h"
#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"
#include "Utilities/OpalException.h"

#include <utility>
#include <variant>
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

        // SYCL compilation crashes when std::visit returns this unique_ptr. Use an explicit
        // dispatch instead: if there are more backends in the future, we could come back to
        // std::visit and try to debug it properly.
        static_assert(
                std::variant_size_v<SpaceChargeConfig> == 2, "Add construction for this algorithm");
        std::unique_ptr<SpaceChargeAlgorithm> algorithm;
        if (auto* selected = std::get_if<CartesianPIC3DConfig>(&config)) {
            if (selected->backend == PoissonSolverType::ConjugateGradient) {
                throw OpalException(
                        "makeSpaceChargeSolver",
                        "The CG Poisson backend is recognized but not implemented.");
            }
            algorithm = std::make_unique<CartesianPIC3DAlgorithm>(
                    std::move(*selected), particles,
                    std::make_unique<CartesianPIC3DFieldStorage<double, 3>>(
                            bunch.cartesianDomain()),
                    dataSink, bunchState);
        } else {
            algorithm = std::make_unique<FFT2D5Algorithm>(
                    std::get<FFT2D5Config>(std::move(config)), particles, bunchState);
        }
        return std::make_unique<SpaceChargeSolver>(std::move(algorithm), particles.size());
    }

}  // namespace opalx::spacecharge
