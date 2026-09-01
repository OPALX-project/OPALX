/**
 * @file P3MShortRangeInteraction.h
 * @brief Adds the particle-particle contribution of the P3M Ewald split.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_P3M_SHORT_RANGE_INTERACTION_H
#define OPALX_SPACE_CHARGE_PIC3D_P3M_SHORT_RANGE_INTERACTION_H

#include "Interaction/TruncatedGreenParticleInteraction.h"
#include "Ippl.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {
    namespace detail {
        template <typename Container>
        class P3MContainerView {
        public:
            explicit P3MContainerView(const Container& particles) : particles_m(particles) {}

            [[nodiscard]] const typename Container::P3MLayout_t& getLayout() const {
                return particles_m.getP3MLayout();
            }

        private:
            const Container& particles_m;
        };

        template <typename View>
        class P3MChargeView {
        public:
            using value_type      = typename View::value_type;
            using execution_space = typename View::execution_space;

            P3MChargeView(View charge, bool perParticle)
                : charge_m(charge), perParticle_m(perParticle) {}

            KOKKOS_INLINE_FUNCTION value_type operator()(std::size_t index) const {
                return charge_m(perParticle_m ? index : 0);
            }

        private:
            View charge_m;
            bool perParticle_m;
        };
    }  // namespace detail

    /** @brief Focused host-side owner of P3M short-range interaction parameters. */
    class P3MShortRangeInteraction final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;

        explicit P3MShortRangeInteraction(double cutoff) : cutoff_m(cutoff) {
            if (!(cutoff_m > 0.0)) {
                throw OpalException(
                        "P3MShortRangeInteraction::P3MShortRangeInteraction",
                        "The P3M cutoff radius must be positive.");
            }
        }

        void apply(ParticleContainer& particles) const {
            if (!particles.hasP3MLayout()) {
                throw OpalException(
                        "P3MShortRangeInteraction::apply",
                        "P3M requires ParticleSpatialOverlapLayout.");
            }

            using ContainerView = detail::P3MContainerView<ParticleContainer>;
            using ChargeView    = detail::P3MChargeView<typename ParticleContainer::qm_view_type>;
            using Interaction   = ippl::TruncatedGreenParticleInteraction<
                      ContainerView, typename ParticleContainer::particle_position_type, ChargeView>;

            ContainerView container(particles);
            ChargeView charge(
                    particles.getQView(),
                    particles.getQMStorageMode() == ParticleContainer::QMStorageMode::Attributes);

            ippl::ParameterList parameters;
            parameters.add("rcut", cutoff_m);
            parameters.add("alpha", 2.0 / cutoff_m);
            // The mesh path scatters charge divided by epsilon_0. This particle path uses raw
            // charge, so its Coulomb coefficient carries epsilon_0 explicitly.
            parameters.add("force_constant", -1.0 / (4.0 * Physics::pi * Physics::epsilon_0));
            parameters.add("regularization_cutoff", 1.0e-9);

            Interaction interaction(container, particles.E, particles.R, charge, parameters);
            interaction.solve();
        }

    private:
        double cutoff_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_P3M_SHORT_RANGE_INTERACTION_H
