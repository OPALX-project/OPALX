/**
 * @file P3MAdapters.h
 * @brief IPPL adapters for the mesh and particle contributions of P3M.
 *
 * Both stages use the same cutoff in metres, alpha = 2 / cutoff, and regularization cutoff
 * 1e-9. The cutoff is the one of the primary container's overlap layout: fixed, or a multiple
 * of the largest mesh spacing that the domain updater sets at every mesh update. Each solve
 * request carries it to the mesh stage, whose IPPL backend rebuilds its kernel when alpha
 * changes. The mesh RHS is already divided by epsilon_0; the particle stage uses raw charges,
 * so only its native force constant contains epsilon_0. The particle contribution is added
 * to E in the Cartesian solve axes after the final mesh-to-particle gather.
 *
 * The open Dirichlet image pass uses the full STANDARD shifted Coulomb mesh kernel, regularized
 * below h_min/2, with no short-range correction for image charges. The image field is therefore
 * mesh-resolved: for particles within about one cell of the plane it has STANDARD OPEN accuracy,
 * while real-bunch interactions keep P3M short-range accuracy. The short-range correction is
 * applied once, to the real bunch.
 */

#ifndef OPALX_SPACE_CHARGE_P3M_ADAPTERS_H
#define OPALX_SPACE_CHARGE_P3M_ADAPTERS_H

#include "Interaction/TruncatedGreenParticleInteraction.h"
#include "PartBunch/CartesianDomainConfig.h"
#include "PartBunch/ParticleContainer.hpp"
#include "Physics/Physics.h"
#include "SpaceCharge/Poisson/PoissonSolver.h"
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
            // Requests replace this with the overlap-layout cutoff of each solve.
            cutoff_m = resolveOverlapCutoff<3>(
                    config_m.p3mCutoff, config_m.p3mCutoffCells,
                    fields.chargeDensity->get_mesh().getMeshSpacing());
            rebuildImpl(fields);
        }

        [[nodiscard]] std::string_view name() const override { return "P3M"; }
        [[nodiscard]] const PoissonSolverCapabilities& capabilities() const override {
            return capabilities_m;
        }
        [[nodiscard]] double couplingConstant() const override { return 1.0 / Physics::epsilon_0; }

    protected:
        void solveImpl(const PoissonSolveRequest& request) override {
            // Separate the long-range truncated-Green solve from the shifted image solve; both
            // include any Green's function rebuild they trigger.
            static IpplTimings::TimerRef meshTimer  = IpplTimings::getTimer("P3M: mesh solve");
            static IpplTimings::TimerRef imageTimer = IpplTimings::getTimer("P3M: image solve");
            const auto timer = request.hasShiftedGreenFunction() ? imageTimer : meshTimer;
            IpplTimings::startTimer(timer);
            if (request.p3mCutoff.has_value() && *request.p3mCutoff != cutoff_m) {
                if (!(*request.p3mCutoff > 0.0)) {
                    IpplTimings::stopTimer(timer);
                    throw OpalException(
                            "P3MMeshPoissonAdapter::solveImpl",
                            "The P3M cutoff radius must be positive.");
                }
                cutoff_m = *request.p3mCutoff;
                backend_m->updateParameter("alpha", 2.0 / cutoff_m);
            }
            kernelRestore_m.solve(*backend_m, request);
            IpplTimings::stopTimer(timer);
        }
        void rebuildImpl(PoissonFieldBinding fields) override {
            auto& backend          = backend_m.emplace();
            kernelRestore_m        = detail::DeferredKernelRestore(fields.chargeDensity);
            const bool allPeriodic = std::all_of(
                    config_m.boundaryConditions.begin(), config_m.boundaryConditions.end(),
                    [](FieldBoundaryCondition boundary) {
                        return boundary == FieldBoundaryCondition::Periodic;
                    });
            auto parameters = detail::commonFftParameters();
            parameters.add("output_type", NativeBackend::GRAD);
            parameters.add("alpha", 2.0 / cutoff_m);
            parameters.add("force_constant", -1.0 / (4.0 * Physics::pi));
            parameters.add("regularization_cutoff", 1.0e-9);
            parameters.add(
                    "boundary_type", allPeriodic ? NativeBackend::PERIODIC : NativeBackend::OPEN);
            capabilities_m.supportsShiftedGreenFunction = !allPeriodic;
            backend.mergeParameters(parameters);
            detail::bindFields(backend, fields);
        }

    private:
        PoissonSolverCapabilities capabilities_m{
                .normalizeChargeByCellVolume    = true,
                .subtractNeutralizingBackground = false,
                .debugDumpChargeBeforeSolve     = true,
                .debugDumpScalarAfterSolve      = true,
                .debugDumpVectorAfterSolve      = true};

        using NativeBackend = FFTTruncatedGreenSolver_t<double, 3>;
        std::optional<NativeBackend> backend_m;
        detail::DeferredKernelRestore kernelRestore_m;
        double cutoff_m = 0.0;  ///< Cutoff in metres of the current split, alpha = 2 / cutoff_m.
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
            typename View::const_type charge_m;
            bool perParticle_m;
        };
    }  // namespace detail

    /**
     * @brief Applies the short-range particle contribution of P3M.
     *
     * The cutoff is read from the container's overlap layout at every call, so the pair test
     * always matches the halos and search cells that the layout built.
     */
    class P3MShortRangeInteraction final {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;

        void apply(ParticleContainer& particles) const {
            if (!particles.hasP3MLayout()) {
                throw OpalException(
                        "P3MShortRangeInteraction::apply",
                        "P3M requires ParticleSpatialOverlapLayout.");
            }
            const double cutoff = particles.getP3MLayout().getCutoff();
            if (!(cutoff > 0.0)) {
                throw OpalException(
                        "P3MShortRangeInteraction::apply",
                        "The P3M cutoff radius must be positive.");
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
            parameters.add("rcut", cutoff);
            parameters.add("alpha", 2.0 / cutoff);
            // The mesh path scatters charge divided by epsilon_0. This particle path uses raw
            // charge, so its Coulomb coefficient carries epsilon_0 explicitly.
            parameters.add("force_constant", -1.0 / (4.0 * Physics::pi * Physics::epsilon_0));
            parameters.add("regularization_cutoff", 1.0e-9);

            Interaction interaction(container, particles.E, particles.R, charge, parameters);
            interaction.solve();
        }
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_P3M_ADAPTERS_H
