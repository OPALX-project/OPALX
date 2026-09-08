/**
 * @file ParticleMeshFieldTransfer.h
 * @brief CIC charge deposition and field gathering for CartesianPIC3D.
 */

#ifndef OPALX_SPACE_CHARGE_PARTICLE_MESH_FIELD_TRANSFER_H
#define OPALX_SPACE_CHARGE_PARTICLE_MESH_FIELD_TRANSFER_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage.h"

#include <cstddef>

namespace opalx::spacecharge {

    /**
     * @brief Performs IPPL CIC scatter and gather operations without retaining views.
     *
     * Deposition temporarily uses @c dt as the charge weight. Image deposition also reflects @c R
     * and negates @c Q. These particle attributes are restored on success.
     */
    class ParticleMeshFieldTransfer final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;
        using PositionAttribute = typename ParticleContainer::particle_position_type;
        using VectorAttribute   = PositionAttribute;
        using FieldStorage      = CartesianPIC3DFieldStorage<double, 3>;
        using ScalarField       = typename FieldStorage::ScalarField;
        using VectorField       = typename FieldStorage::VectorField;
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
         * @brief Selects a contiguous range or maps policy indices through a hash.
         *
         * An indexed hash is borrowed and becomes invalid after particle migration.
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

        /** @brief Image plane for one deposition pass, in metres in solve axes. */
        struct ImagePolicy {
            bool enabled  = false;
            double planeZ = 0.0;
        };

        /**
         * @brief Converts deposited @c dt*Q weights to the backend's charge convention.
         *
         * @c selectedCharge is the selection's global charge and is used only for the periodic
         * neutralizing background.
         */
        struct ChargeNormalization {
            double timeStep                     = 0.0;
            double gamma                        = 1.0;
            double selectedCharge               = 0.0;
            double couplingConstant             = 1.0;
            bool normalizeByCellVolume          = false;
            bool subtractNeutralizingBackground = false;
        };

        /** @brief Clear rho, deposit the selected charge source, and normalize it in place. */
        void depositCharge(
                ParticleContainer& particles, FieldStorage& fieldStorage, DepositKind depositKind,
                const Selection& selection, const ChargeNormalization& normalization,
                const ImagePolicy& imagePolicy) const;

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
                FieldStorage& fieldStorage, const ChargeNormalization& normalization) const;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PARTICLE_MESH_FIELD_TRANSFER_H
