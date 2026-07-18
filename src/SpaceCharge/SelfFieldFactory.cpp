/**
 * @file SelfFieldFactory.cpp
 * @brief Implements construction-time self-field algorithm selection.
 */

#include "SpaceCharge/SelfFieldFactory.h"

#include "PartBunch/BinnedFieldSolver.h"
#include "PartBunch/PartBunch.h"
#include "SpaceCharge/Pic/Pic3DSolver.h"
#include "Utilities/OpalException.h"

#include <memory>
#include <utility>

namespace opalx::spacecharge {

    std::unique_ptr<SelfFieldSystem> SelfFieldFactory::create(
            SelfFieldConfig config, PartBunch_t& bunch) {
        std::unique_ptr<SelfFieldAlgorithm> algorithm;
        switch (config.algorithmKind()) {
            case SelfFieldAlgorithmKind::Pic3D:
                if (bunch.getFieldSolver() == nullptr) {
                    throw OpalException(
                            "SelfFieldFactory::create",
                            "The temporary PIC backend facade is not initialized.");
                }
                algorithm = std::make_unique<Pic3DSolver>(
                        config.get<Pic3DConfig>(), bunch.getParticleContainers(),
                        bunch.sharedPicWorkspace(), *bunch.getFieldSolver(), bunch.getDataSink());
                break;
        }

        if (algorithm == nullptr) {
            throw OpalException(
                    "SelfFieldFactory::create", "No implementation exists for this algorithm.");
        }
        return std::make_unique<SelfFieldSystem>(std::move(config), std::move(algorithm));
    }

}  // namespace opalx::spacecharge
