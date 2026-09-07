/**
 * @file P3MAdapters.h
 * @brief IPPL adapters for the mesh and particle contributions of P3M.
 *
 * Both stages use the same cutoff in metres, alpha = 2 / cutoff, and regularization cutoff
 * 1e-9. The mesh RHS is already divided by epsilon_0; the particle stage uses raw charges,
 * so only its native force constant contains epsilon_0. The particle contribution is added
 * to E in the Cartesian solve axes after the final mesh-to-particle gather.
 */

#ifndef OPALX_SPACE_CHARGE_P3M_ADAPTERS_H
#define OPALX_SPACE_CHARGE_P3M_ADAPTERS_H

#include "Interaction/TruncatedGreenParticleInteraction.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "SpaceCharge/Poisson/PoissonSolverFactory.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <utility>

namespace opalx::spacecharge {

    /** @brief Adapts the 3D mesh contribution of P3M; particle interactions remain in Cartesian
     * PIC. */
    class P3MMeshPoissonAdapter final : public PoissonSolver {
    public:
        P3MMeshPoissonAdapter(PoissonSolverConfig config, PoissonFieldBinding fields)
            : PoissonSolver(std::move(config), fields, PoissonSolverType::P3M) {
            rebuildImpl(fields);
        }

        [[nodiscard]] std::string_view name() const override { return "P3M"; }
        [[nodiscard]] const PoissonSolverCapabilities& capabilities() const override {
            return capabilities_m;
        }
        [[nodiscard]] double couplingConstant() const override { return 1.0 / Physics::epsilon_0; }

    protected:
        void solveImpl(const PoissonSolveRequest&) override { backend_m->solve(); }
        void rebuildImpl(PoissonFieldBinding fields) override {
            auto& backend          = backend_m.emplace();
            const bool allPeriodic = std::all_of(
                    config_m.boundaryConditions.begin(), config_m.boundaryConditions.end(),
                    [](FieldBoundaryCondition boundary) {
                        return boundary == FieldBoundaryCondition::Periodic;
                    });
            auto parameters = detail::commonFftParameters();
            parameters.add("output_type", NativeBackend::GRAD);
            parameters.add("alpha", 2.0 / config_m.p3mCutoff);
            parameters.add("force_constant", -1.0 / (4.0 * Physics::pi));
            parameters.add("regularization_cutoff", 1.0e-9);
            parameters.add(
                    "boundary_type", allPeriodic ? NativeBackend::PERIODIC : NativeBackend::OPEN);
            backend.mergeParameters(parameters);
            detail::bindFields(backend, fields);
        }

    private:
        static constexpr PoissonSolverCapabilities capabilities_m{
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = false,
                .debugDumpChargeBeforeSolve     = true,
                .debugDumpScalarAfterSolve      = true,
                .debugDumpVectorAfterSolve      = true};

        using NativeBackend = FFTTruncatedGreenSolver_t<double, 3>;
        std::optional<NativeBackend> backend_m;
    };

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

#endif  // OPALX_SPACE_CHARGE_P3M_ADAPTERS_H
