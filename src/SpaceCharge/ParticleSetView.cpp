/**
 * @file ParticleSetView.cpp
 * @brief Implements borrowed particle-set validation and metadata queries.
 */

#include "SpaceCharge/ParticleSetView.h"

#include <array>
#include <stdexcept>

namespace opalx::spacecharge {
    namespace {
        using OptionalHandle = std::optional<ParticleAttributeHandle>;

        const OptionalHandle& selectHandle(
                const ParticleContainerAttributes& attributes, ParticleAttribute attribute) {
            switch (attribute) {
                case ParticleAttribute::Position:
                    return attributes.position;
                case ParticleAttribute::Momentum:
                    return attributes.momentum;
                case ParticleAttribute::Charge:
                    return attributes.charge;
                case ParticleAttribute::Mass:
                    return attributes.mass;
                case ParticleAttribute::TimeStep:
                    return attributes.timeStep;
                case ParticleAttribute::ElectricField:
                    return attributes.electricField;
                case ParticleAttribute::MagneticField:
                    return attributes.magneticField;
                case ParticleAttribute::InvalidMask:
                    return attributes.invalidMask;
                case ParticleAttribute::Bin:
                    return attributes.bin;
                case ParticleAttribute::Count:
                    break;
            }
            throw std::invalid_argument("ParticleAttribute::Count is not an attribute");
        }

        constexpr std::array particleAttributes{ParticleAttribute::Position,
                                                ParticleAttribute::Momentum,
                                                ParticleAttribute::Charge,
                                                ParticleAttribute::Mass,
                                                ParticleAttribute::TimeStep,
                                                ParticleAttribute::ElectricField,
                                                ParticleAttribute::MagneticField,
                                                ParticleAttribute::InvalidMask,
                                                ParticleAttribute::Bin};
    }  // namespace

    const ParticleAttributeHandle* ParticleContainerAttributes::find(
            ParticleAttribute attribute) const {
        const auto& handle = selectHandle(*this, attribute);
        return handle ? std::addressof(*handle) : nullptr;
    }

    ParticleAttributeSet ParticleContainerAttributes::available() const {
        ParticleAttributeSet result;
        for (const ParticleAttribute attribute : particleAttributes) {
            if (find(attribute) != nullptr) {
                result.add(attribute);
            }
        }
        return result;
    }

    ParticleAttributeSet ParticleContainerAttributes::writable() const {
        ParticleAttributeSet result;
        for (const ParticleAttribute attribute : particleAttributes) {
            const ParticleAttributeHandle* handle = find(attribute);
            if (handle != nullptr && handle->isWritable()) {
                result.add(attribute);
            }
        }
        return result;
    }

    ParticleContainerView::ParticleContainerView(
            std::string_view name, ParticleContainerAttributes attributes, bool active)
        : name_m(name), attributes_m(std::move(attributes)), active_m(active) {
        for (const ParticleAttribute attribute : particleAttributes) {
            const ParticleAttributeHandle* handle = find(attribute);
            if (handle != nullptr && handle->attribute() != attribute) {
                throw std::invalid_argument(
                        "ParticleContainerView attribute handle is stored in the wrong slot");
            }
        }
    }

    const ParticleAttributeHandle& ParticleContainerView::require(
            ParticleAttribute attribute) const {
        const ParticleAttributeHandle* handle = find(attribute);
        if (handle == nullptr) {
            throw std::invalid_argument(
                    "ParticleContainerView does not provide required attribute");
        }
        return *handle;
    }

    ParticleSetView::ParticleSetView(
            std::span<ParticleContainerView> containers, std::size_t primaryIndex,
            generation_type generation)
        : containers_m(containers), primaryIndex_m(primaryIndex), generation_m(generation) {
        if (containers_m.empty()) {
            throw std::invalid_argument("ParticleSetView requires at least one container");
        }
        if (primaryIndex_m >= containers_m.size()) {
            throw std::out_of_range("ParticleSetView primary index is outside the container set");
        }
        if (!primary().active()) {
            throw std::invalid_argument("ParticleSetView primary container must be active");
        }
    }

    std::size_t ParticleSetView::activeContainerCount() const {
        std::size_t count = 0;
        for (const ParticleContainerView& container : containers_m) {
            count += container.active() ? 1 : 0;
        }
        return count;
    }

    bool ParticleSetView::activeContainersProvide(ParticleAttributeSet required) const {
        for (const ParticleContainerView& container : containers_m) {
            if (container.active() && !container.availableAttributes().contains(required)) {
                return false;
            }
        }
        return true;
    }

    bool ParticleSetView::activeContainersProvideWritable(ParticleAttributeSet required) const {
        for (const ParticleContainerView& container : containers_m) {
            if (container.active() && !container.writableAttributes().contains(required)) {
                return false;
            }
        }
        return true;
    }

}  // namespace opalx::spacecharge
