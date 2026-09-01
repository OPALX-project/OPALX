/**
 * @file Pic2d5FramePolicy.h
 * @brief Preserves the production FFT2D5 tracker-frame compatibility ordering.
 */

#ifndef OPALX_SPACE_CHARGE_PIC2D5_FRAME_POLICY_H
#define OPALX_SPACE_CHARGE_PIC2D5_FRAME_POLICY_H

#include "SpaceCharge/SolveContext.h"

template <typename T, unsigned Dim>
class ParticleContainer;

namespace opalx::spacecharge {

    /**
     * @brief Applies the exact master FFT2D5 outer-frame contract to the primary container.
     *
     * Master transforms primary positions into the beam frame before the solver's own Frenet
     * mapping, then restores primary R and rotates primary E and B after the solve. Other active
     * containers participate in FFT2D5 without that outer transform. This policy deliberately
     * retains that ordering until the reference-frame contract is redesigned separately.
     */
    class Pic2d5FramePolicy final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;

        struct State final {
            bool primaryPositionsInSolveFrame = false;
            bool primaryElectricInSolveFrame  = false;
            bool primaryMagneticInSolveFrame  = false;
        };

        void validate(const SolveContext& context) const;
        void enter(const SolveContext& context, ParticleContainer& primary, State& state) const;
        void markComputedFields(State& state) const noexcept;
        void leave(const SolveContext& context, ParticleContainer& primary, State& state) const;
        void restore(const SolveContext& context, ParticleContainer& primary, State& state)
                const noexcept;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC2D5_FRAME_POLICY_H
