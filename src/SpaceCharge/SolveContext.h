/**
 * @file SolveContext.h
 * @brief Defines all borrowed state needed for one self-field solve.
 */

#ifndef OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
#define OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H

#include "SpaceCharge/ParticleSetView.h"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <type_traits>
#include <typeindex>
#include <typeinfo>

namespace opalx::spacecharge {

    /**
     * @brief Checked borrowed reference to an otherwise implementation-specific host object.
     *
     * This is used for step services such as a communicator or coordinate transform without
     * making their concrete types part of the common solver interface. It is never passed to a
     * device kernel.
     */
    class BorrowedHostObject {
    public:
        template <typename Native>
        [[nodiscard]] static BorrowedHostObject reference(const Native& native) {
            return BorrowedHostObject(native);
        }

        [[nodiscard]] std::type_index nativeType() const { return nativeType_m; }

        template <typename Native>
        [[nodiscard]] const std::remove_cv_t<Native>& native() const {
            using Value = std::remove_cv_t<Native>;
            if (nativeType_m != std::type_index(typeid(Value))) {
                throw std::bad_cast();
            }
            return *static_cast<const Value*>(native_m);
        }

    private:
        template <typename Native>
        explicit BorrowedHostObject(const Native& native)
            : nativeType_m(typeid(std::remove_cv_t<Native>)), native_m(std::addressof(native)) {}

        std::type_index nativeType_m;
        const void* native_m = nullptr;
    };

    /** @brief Borrowed communicator plus stable rank metadata for the current call. */
    struct CommunicatorView {
        BorrowedHostObject native;
        int rank = 0;
        int size = 1;
    };

    /** @brief Reference-particle data available to every self-field algorithm. */
    struct ReferenceState {
        std::array<double, 3> position{0.0, 0.0, 0.0};
        std::array<double, 3> momentum{0.0, 0.0, 0.0};
        double pathLength = 0.0;
    };

    /**
     * @brief Coordinate-transform objects valid for the current call.
     *
     * The tracker supplies its native transform type through checked host handles. A concrete
     * algorithm decides whether and how to use those transforms; the common interface does not
     * prescribe a Cartesian or Frenet representation.
     */
    struct FrameState {
        std::optional<BorrowedHostObject> labToReference;
        std::optional<BorrowedHostObject> referenceToLab;
    };

    /** @brief Tracker state captured for one self-field solve. */
    struct StepState {
        std::size_t step            = 0;
        double time                 = 0.0;
        double timeStep             = 0.0;
        std::size_t emittedBinCount = 0;
        bool emissionActive         = false;
        CommunicatorView communicator;
        ReferenceState reference;
        FrameState frames;
    };

    /** @brief Step-dependent correction data, already resolved from source configuration. */
    struct CorrectionRequest {
        CorrectionKind kind = CorrectionKind::None;
        double planeZ       = 0.0;
    };

    /**
     * @brief Physics requested for this call after construction-time configuration is applied.
     */
    struct RequestedPhysics {
        bool useBinning     = false;
        bool writePotential = false;
        CorrectionRequest correction;
    };

    /**
     * @brief Borrowed inputs and explicit writable particle targets for one self-field call.
     *
     * SolveContext owns no particle storage, solver, mesh, backend cache, or persistent scratch.
     * It must not be retained by SelfFieldSystem or SelfFieldAlgorithm after execute() returns.
     * Writable particle E/B targets are identified by writable handles in particles().
     */
    class SolveContext {
    public:
        SolveContext(
                ParticleSetView particles, StepState stepState,
                RequestedPhysics requestedPhysics = {});

        [[nodiscard]] ParticleSetView& particles() { return particles_m; }
        [[nodiscard]] const ParticleSetView& particles() const { return particles_m; }
        [[nodiscard]] const StepState& stepState() const { return stepState_m; }
        [[nodiscard]] const RequestedPhysics& requestedPhysics() const {
            return requestedPhysics_m;
        }

    private:
        ParticleSetView particles_m;
        StepState stepState_m;
        RequestedPhysics requestedPhysics_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_SOLVE_CONTEXT_H
