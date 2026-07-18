/**
 * @file FieldFramePolicy.h
 * @brief Declares tracker-to-beam frame handling for Cartesian 3D PIC.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_FIELD_FRAME_POLICY_H
#define OPALX_SPACE_CHARGE_PIC3D_FIELD_FRAME_POLICY_H

#include "SpaceCharge/SolveContext.h"

template <typename T, unsigned Dim>
class ParticleContainer;

namespace opalx::spacecharge {

    /**
     * @brief Applies the caller-frame boundary around one Cartesian 3D PIC solve.
     *
     * This is a host-only orchestration policy. It obtains checked CoordinateSystemTrafo handles
     * from the current SolveContext and launches the native Kokkos transformations on the current
     * particle attributes. The policy stores neither transform handles nor particle views.
     *
     * Particle migration or attribute reallocation may occur between enter() and leave(). Every
     * pending transformation therefore reacquires the current R, E, or B view immediately before
     * launching its kernel. The ParticleContainer object and the borrowed transforms themselves
     * must remain alive for the complete call.
     *
     * These methods must be called from host code, never from a device lambda. They deliberately
     * add no fences: a GPU caller must preserve Kokkos execution ordering and synchronize before
     * unordered host access or storage replacement not already synchronized by the owning IPPL
     * operation.
     */
    class FieldFramePolicy final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;

        /** @brief Frames currently held by the particle attributes. */
        struct State final {
            bool positionsInBeam = false;
            bool electricInBeam  = false;
            bool magneticInBeam  = false;
        };

        /** @brief Validate both native frame transforms without mutating particle data. */
        void validate(const SolveContext& context) const;

        /**
         * @brief Transform positions from the tracker frame into the beam frame.
         *
         * The position flag is set before the Kokkos transformation is launched so exception
         * cleanup conservatively treats a partially completed launch as beam-frame data.
         */
        void enter(const SolveContext& context, ParticleContainer& particles, State& state) const;

        /**
         * @brief Mark E and B as beam-frame outputs before the field computation starts.
         *
         * Call this immediately before an operation that may partially write either field.
         */
        void markComputedFieldsInBeam(State& state) const noexcept;

        /**
         * @brief Restore every pending attribute to the tracker frame in R, E, B order.
         *
         * Each state flag is cleared only after its corresponding transformation returns
         * successfully.
         */
        void leave(const SolveContext& context, ParticleContainer& particles, State& state) const;

        /**
         * @brief Best-effort exception cleanup using the same ordered pending transformations.
         *
         * Cleanup failures are suppressed so a caller can preserve and rethrow the exception that
         * interrupted the solve.
         */
        void restore(const SolveContext& context, ParticleContainer& particles, State& state)
                const noexcept;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_FIELD_FRAME_POLICY_H
