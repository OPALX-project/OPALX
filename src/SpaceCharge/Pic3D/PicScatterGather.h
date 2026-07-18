/**
 * @file PicScatterGather.h
 * @brief Declares CIC charge deposition and field gathering for Cartesian PIC.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_SCATTER_GATHER_H
#define OPALX_SPACE_CHARGE_PIC3D_SCATTER_GATHER_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/Pic3D/PicWorkspace.h"

#include <cstddef>

namespace opalx::spacecharge {

    /**
     * @brief Stateless host-side facade over the IPPL CIC scatter/gather operations.
     *
     * Particle attributes and fields are borrowed for the duration of one call. The operation
     * never retains native Kokkos views, so callers may safely reacquire attributes after a
     * particle migration before invoking it again. Charge deposition temporarily mutates
     * particle @c dt and, for an image pass, @c R and @c Q; each mutation is restored before the
     * call returns or propagates an exception.
     *
     * @tparam T Particle and field scalar type.
     * @tparam Dim Cartesian dimension. The production implementation is @c double,3.
     */
    template <typename T, unsigned Dim>
    class PicScatterGather final {
        static_assert(Dim == 3, "PicScatterGather currently supports Dim == 3 only.");

    public:
        using ParticleContainer = ::ParticleContainer<T, Dim>;
        using PositionAttribute = typename ParticleContainer::particle_position_type;
        using VectorAttribute   = PositionAttribute;
        using Workspace         = PicWorkspace<T, Dim>;
        using ScalarField       = typename Workspace::ScalarField;
        using VectorField       = typename Workspace::VectorField;
        using ExecutionSpace    = Kokkos::DefaultExecutionSpace;
        using MemorySpace       = typename ExecutionSpace::memory_space;
        using RangePolicy       = Kokkos::RangePolicy<ExecutionSpace>;
        using Hash              = ippl::detail::hash_type<MemorySpace>;
        using size_type         = std::size_t;

        /** @brief Charge source deposited in one ordered PIC pass. */
        enum class DepositKind { Primary, Image, PrimaryAndImage };

        /** @brief Whether gathered values replace or add to the particle destination. */
        enum class GatherMode { Replace, Add };

        /**
         * @brief Borrowed direct or hash-indexed particle selection.
         *
         * The hash allocation, when present, must remain valid until depositCharge() returns.
         * A selection never owns particle storage and must not survive a particle migration.
         */
        struct Selection {
            enum class Kind { Direct, Indexed };

            /** @brief Select a contiguous particle-index range without indirection. */
            [[nodiscard]] static Selection direct(size_type begin, size_type end) {
                return Selection(Kind::Direct, RangePolicy(begin, end), Hash());
            }

            /** @brief Select policy indices and map each one through @p hash. */
            [[nodiscard]] static Selection indexed(const RangePolicy& policy, const Hash& hash) {
                return Selection(Kind::Indexed, policy, hash);
            }

            [[nodiscard]] Kind kind() const { return kind_m; }
            [[nodiscard]] const RangePolicy& policy() const { return policy_m; }
            [[nodiscard]] const Hash& hash() const { return hash_m; }

        private:
            Selection(Kind kind, RangePolicy policy, Hash hash)
                : kind_m(kind), policy_m(policy), hash_m(hash) {}

            Kind kind_m;
            RangePolicy policy_m;
            Hash hash_m;
        };

        /** @brief Value policy controlling temporary image-charge mutation. */
        struct ImagePolicy {
            bool enabled  = false;
            double planeZ = 0.0;
        };

        /**
         * @brief Values used to convert deposited @c dt*Q weights into charge density.
         *
         * @c selectedCharge is the global primary charge represented by the selection. It is
         * used only for the periodic neutralizing background. Mesh spacing and physical volume
         * are read from the borrowed workspace.
         */
        struct ChargeNormalization {
            double timeStep                     = 0.0;
            double gamma                        = 1.0;
            double selectedCharge               = 0.0;
            double couplingConstant             = 1.0;
            bool normalizeByCellVolume          = false;
            bool subtractNeutralizingBackground = false;
        };

        /**
         * @brief Clear rho, deposit the requested charge source, and normalize it in place.
         *
         * Whole-bunch and bin-restricted work share this entry point. A direct selection that
         * spans all local particles uses IPPL's all-particle scatter overload; all other direct
         * and indexed selections use its custom-policy overload.
         */
        void depositCharge(
                ParticleContainer& particles, Workspace& workspace, DepositKind depositKind,
                const Selection& selection, const ChargeNormalization& normalization,
                const ImagePolicy& imagePolicy = {}) const;

        /**
         * @brief Gather one vector field into an explicit writable particle destination.
         *
         * IPPL fills source halos during gather, so @p source is intentionally non-const.
         */
        void gatherVector(
                VectorAttribute& destination, VectorField& source,
                const PositionAttribute& positions, GatherMode mode) const;

        // CUDA builds require enclosing functions that launch device lambdas to be public.

        /** @brief Deposit the current signed charge for one selection via @c dt*Q. */
        void scatterScaledTimeStep(
                ParticleContainer& particles, PositionAttribute& positions, ScalarField& rho,
                const Selection& selection) const;

        /** @brief Route current weights through IPPL's direct or indexed CIC scatter. */
        void scatterCurrentWeights(
                const ParticleContainer& particles, PositionAttribute& positions, ScalarField& rho,
                const Selection& selection) const;

        /** @brief Apply the image transform to the selected positions and charges. */
        void applyImageTransform(
                ParticleContainer& particles, PositionAttribute& positions,
                const Selection& selection, double planeZ) const;

        /** @brief Restore the self-inverse image transform. */
        void restoreImageTransform(
                ParticleContainer& particles, PositionAttribute& positions,
                const Selection& selection, double planeZ) const;

        /** @brief Reflect selected particle positions around an xy plane. */
        void reflectPositions(
                PositionAttribute& positions, const Selection& selection, double planeZ) const;

        /** @brief Flip selected charge signs, respecting the configured Q storage mode. */
        void flipChargeSign(ParticleContainer& particles, const Selection& selection) const;

    private:
        void validateSelection(
                const ParticleContainer& particles, const Selection& selection) const;

        void depositImage(
                ParticleContainer& particles, PositionAttribute& positions, ScalarField& rho,
                const Selection& selection, const ImagePolicy& imagePolicy) const;

        void normalizeChargeDensity(
                Workspace& workspace, const ChargeNormalization& normalization) const;
    };

    extern template class PicScatterGather<double, 3>;

}  // namespace opalx::spacecharge

#include "SpaceCharge/Pic3D/PicScatterGather.tpp"

#endif  // OPALX_SPACE_CHARGE_PIC3D_SCATTER_GATHER_H
