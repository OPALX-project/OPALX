/**
 * @file SelfFieldFactory.h
 * @brief Constructs the self-field system selected by immutable setup configuration.
 */

#ifndef OPALX_SPACE_CHARGE_SELF_FIELD_FACTORY_H
#define OPALX_SPACE_CHARGE_SELF_FIELD_FACTORY_H

#include "PartBunch/PartBunchFwd.h"
#include "SpaceCharge/SelfFieldConfig.h"
#include "SpaceCharge/SelfFieldSystem.h"

#include <memory>

class DataSink;

namespace opalx::spacecharge {

    /** @brief Construction-time selection point for concrete self-field algorithms. */
    class SelfFieldFactory {
    public:
        [[nodiscard]] static std::unique_ptr<SelfFieldSystem> create(
                SelfFieldConfig config, PartBunch_t& bunch, DataSink* dataSink);
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SELF_FIELD_FACTORY_H
