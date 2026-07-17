/**
 * @file SelfFieldFactory.cpp
 * @brief Implements construction-time self-field algorithm selection.
 */

#include "SpaceCharge/SelfFieldFactory.h"

#include "PartBunch/PartBunch.h"
#include "SpaceCharge/LegacyPic3DAlgorithm.h"
#include "Utilities/OpalException.h"

#include <memory>
#include <utility>

namespace opalx::spacecharge {

    std::unique_ptr<SelfFieldSystem> SelfFieldFactory::create(
            SelfFieldConfig config, PartBunch_t& bunch) {
        std::unique_ptr<SelfFieldAlgorithm> algorithm;
        switch (config.algorithmKind()) {
            case SelfFieldAlgorithmKind::Pic3D:
                algorithm = std::make_unique<LegacyPic3DAlgorithm>(
                        bunch, config.get<Pic3DConfig>(), bunch.sharedPicWorkspace());
                break;
        }

        if (algorithm == nullptr) {
            throw OpalException(
                    "SelfFieldFactory::create", "No implementation exists for this algorithm.");
        }
        return std::make_unique<SelfFieldSystem>(std::move(config), std::move(algorithm));
    }

}  // namespace opalx::spacecharge
