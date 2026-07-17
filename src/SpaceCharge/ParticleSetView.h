/**
 * @file ParticleSetView.h
 * @brief Borrowed, collection-shaped particle contracts for one self-field call.
 */

#ifndef OPALX_SPACE_CHARGE_PARTICLE_SET_VIEW_H
#define OPALX_SPACE_CHARGE_PARTICLE_SET_VIEW_H

#include "SpaceCharge/SolverCapabilities.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <utility>

namespace opalx::spacecharge {

    /**
     * @brief Stable borrowed handle to one native particle-attribute object.
     *
     * The handle points to an attribute object, such as a concrete particle attribute wrapper,
     * and deliberately does not cache a Kokkos view. A concrete host-side adapter recovers the
     * checked native type with native() or writableNative(), then obtains the current device view
     * from that object. This reacquisition is required after particle storage is reallocated or
     * migrated.
     *
     * The pointed-to object must outlive the SolveContext. No method on this class is intended for
     * use in a device kernel.
     */
    class ParticleAttributeHandle {
    public:
        /** @brief Borrow a read-only native attribute object. */
        template <typename NativeAttribute>
        [[nodiscard]] static ParticleAttributeHandle readable(
                ParticleAttribute attribute, const NativeAttribute& native) {
            return ParticleAttributeHandle(attribute, native, nullptr);
        }

        /** @brief Borrow a native attribute object that the solver may modify. */
        template <typename NativeAttribute>
        [[nodiscard]] static ParticleAttributeHandle writable(
                ParticleAttribute attribute, NativeAttribute& native) {
            return ParticleAttributeHandle(attribute, native, std::addressof(native));
        }

        [[nodiscard]] ParticleAttribute attribute() const { return attribute_m; }
        [[nodiscard]] bool isWritable() const { return mutableNative_m != nullptr; }
        [[nodiscard]] std::type_index nativeType() const { return nativeType_m; }

        /**
         * @brief Recover a checked read-only reference to the native attribute object.
         * @throws std::bad_cast when @p NativeAttribute does not match the borrowed native type.
         */
        template <typename NativeAttribute>
        [[nodiscard]] const std::remove_cv_t<NativeAttribute>& native() const {
            using Native = std::remove_cv_t<NativeAttribute>;
            if (nativeType_m != std::type_index(typeid(Native))) {
                throw std::bad_cast();
            }
            return *static_cast<const Native*>(native_m);
        }

        /**
         * @brief Recover a checked writable reference to the native attribute object.
         * @throws std::bad_cast when @p NativeAttribute does not match the borrowed native type.
         * @throws std::logic_error when the attribute was not exposed as writable.
         */
        template <typename NativeAttribute>
        [[nodiscard]] std::remove_cv_t<NativeAttribute>& writableNative() const {
            using Native = std::remove_cv_t<NativeAttribute>;
            if (nativeType_m != std::type_index(typeid(Native))) {
                throw std::bad_cast();
            }
            if (mutableNative_m == nullptr) {
                throw std::logic_error("particle attribute is not writable");
            }
            return *static_cast<Native*>(mutableNative_m);
        }

    private:
        template <typename NativeAttribute>
        ParticleAttributeHandle(
                ParticleAttribute attribute, const NativeAttribute& native, void* mutableNative)
            : attribute_m(attribute),
              nativeType_m(typeid(std::remove_cv_t<NativeAttribute>)),
              native_m(std::addressof(native)),
              mutableNative_m(mutableNative) {}

        ParticleAttribute attribute_m;
        std::type_index nativeType_m;
        const void* native_m  = nullptr;
        void* mutableNative_m = nullptr;
    };

    /** @brief Named attribute handles exposed by one particle container. */
    struct ParticleContainerAttributes {
        std::optional<ParticleAttributeHandle> position;
        std::optional<ParticleAttributeHandle> momentum;
        std::optional<ParticleAttributeHandle> charge;
        std::optional<ParticleAttributeHandle> mass;
        std::optional<ParticleAttributeHandle> timeStep;
        std::optional<ParticleAttributeHandle> electricField;
        std::optional<ParticleAttributeHandle> magneticField;
        std::optional<ParticleAttributeHandle> invalidMask;
        std::optional<ParticleAttributeHandle> bin;

        [[nodiscard]] const ParticleAttributeHandle* find(ParticleAttribute attribute) const;
        [[nodiscard]] ParticleAttributeSet available() const;
        [[nodiscard]] ParticleAttributeSet writable() const;
    };

    /**
     * @brief Borrowed metadata and stable attribute handles for one particle container.
     *
     * This class owns no particle storage. The caller guarantees that all native attribute
     * objects and the container name outlive the SolveContext that contains this view. Solve
     * selection is independent of whether tracking currently considers the container active.
     */
    class ParticleContainerView {
    public:
        ParticleContainerView(
                std::string_view name, ParticleContainerAttributes attributes,
                bool selectedForSolve = true, bool trackingActive = true);

        [[nodiscard]] std::string_view name() const { return name_m; }
        /** @brief Compatibility name for selectedForSolve(). */
        [[nodiscard]] bool active() const { return selectedForSolve_m; }
        [[nodiscard]] bool selectedForSolve() const { return selectedForSolve_m; }
        [[nodiscard]] bool trackingActive() const { return trackingActive_m; }
        void setTrackingActive(bool active) { trackingActive_m = active; }
        [[nodiscard]] const ParticleContainerAttributes& attributes() const { return attributes_m; }
        [[nodiscard]] const ParticleAttributeHandle* find(ParticleAttribute attribute) const {
            return attributes_m.find(attribute);
        }
        [[nodiscard]] const ParticleAttributeHandle& require(ParticleAttribute attribute) const;
        [[nodiscard]] ParticleAttributeSet availableAttributes() const {
            return attributes_m.available();
        }
        [[nodiscard]] ParticleAttributeSet writableAttributes() const {
            return attributes_m.writable();
        }

    private:
        std::string_view name_m;
        ParticleContainerAttributes attributes_m;
        bool selectedForSolve_m = true;
        bool trackingActive_m   = true;
    };

    /**
     * @brief Non-owning collection of particle containers used by one self-field solve.
     *
     * The primary container is selected explicitly; algorithms must not infer it from a global
     * object. generation() records the particle-storage generation at construction. Any native
     * Kokkos views recovered through attribute handles belong to that generation and must not be
     * retained across storage reallocation, compaction, or migration. Recover the native
     * attribute object again and ask it for a fresh view after such an operation. Rebuild this
     * ParticleSetView if container membership or the primary container changes. Update the
     * generation after storage changes and update container flags before each solve.
     */
    class ParticleSetView {
    public:
        using generation_type = std::uint64_t;

        ParticleSetView(
                std::span<ParticleContainerView> containers, std::size_t primaryIndex,
                generation_type generation);

        [[nodiscard]] std::span<ParticleContainerView> containers() const { return containers_m; }
        [[nodiscard]] ParticleContainerView& primary() const {
            return containers_m[primaryIndex_m];
        }
        [[nodiscard]] std::size_t primaryIndex() const { return primaryIndex_m; }
        [[nodiscard]] generation_type generation() const { return generation_m; }
        /** @brief Record a host-side particle storage change after all native views are invalid. */
        void updateGeneration(generation_type generation) { generation_m = generation; }
        [[nodiscard]] std::size_t activeContainerCount() const;

        [[nodiscard]] bool activeContainersProvide(ParticleAttributeSet required) const;
        [[nodiscard]] bool activeContainersProvideWritable(ParticleAttributeSet required) const;

    private:
        std::span<ParticleContainerView> containers_m;
        std::size_t primaryIndex_m   = 0;
        generation_type generation_m = 0;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PARTICLE_SET_VIEW_H
