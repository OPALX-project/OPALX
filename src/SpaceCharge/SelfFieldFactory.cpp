/**
 * @file SelfFieldFactory.cpp
 * @brief Implements construction-time self-field algorithm selection.
 */

#include "SpaceCharge/SelfFieldFactory.h"

#include "PartBunch/PartBunch.h"
#include "SpaceCharge/Pic2d5/Pic2d5Solver.h"
#include "SpaceCharge/Pic3D/Pic3DSolver.h"
#include "Utilities/OpalException.h"

#include <memory>
#include <utility>
#include <vector>

namespace opalx::spacecharge {

    std::unique_ptr<SelfFieldSystem> SelfFieldFactory::create(
            SelfFieldConfig config, PartBunch_t& bunch, DataSink* dataSink) {
        if (config.algorithmKind() == SelfFieldAlgorithmKind::Pic3D
            && config.get<Pic3DConfig>().backend() == PoissonBackendKind::ConjugateGradient) {
            throw OpalException(
                    "SelfFieldFactory::create",
                    "The CG Poisson backend is recognized but not implemented.");
        }

        std::vector<ParticleFieldBinding3d> bindings;
        bindings.reserve(bunch.getNumParticleContainers());
        for (const auto& container : bunch.getParticleContainers()) {
            if (container == nullptr) {
                throw OpalException(
                        "SelfFieldFactory::create",
                        "Cannot construct a self-field solver for a null particle container.");
            }
            bindings.push_back(bindParticleFields(*container));
        }

        std::unique_ptr<SelfFieldAlgorithm> algorithm;
        switch (config.algorithmKind()) {
            case SelfFieldAlgorithmKind::Pic3D:
                algorithm = std::make_unique<Pic3DSolver>(
                        config.get<Pic3DConfig>(), bindings,
                        std::make_unique<PicWorkspace<double, 3>>(bunch.cartesianDomain()),
                        dataSink);
                break;
            case SelfFieldAlgorithmKind::Pic2d5:
                algorithm = std::make_unique<Pic2d5Solver>(config.get<Pic2d5Config>(), bindings);
                break;
        }

        if (algorithm == nullptr) {
            throw OpalException(
                    "SelfFieldFactory::create", "No implementation exists for this algorithm.");
        }
        return std::make_unique<SelfFieldSystem>(
                std::move(config), std::move(algorithm), std::move(bindings));
    }

}  // namespace opalx::spacecharge
