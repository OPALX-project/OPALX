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
 * to E, and with binning also to B, in the Cartesian solve axes after the final
 * mesh-to-particle gather.
 *
 * With binning, the mesh stage solves each bin in the bin rest frame on a z spacing stretched by
 * the bin gamma. The particle stage splits the Ewald sum in the rest frame of each source particle
 * instead, with the particle's own momentum, and adds lab-frame E and B directly, so only the
 * long-range field carries the bin velocity. The cutoff keeps its value in metres and is the
 * rest-frame cutoff in both stages.
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
#include <cmath>
#include <cstddef>
#include <utility>

namespace opalx::spacecharge {

    /// Radius below which both P3M stages replace the singular Coulomb term by a linear core [m].
    inline constexpr double p3mRegularizationCutoff = 1.0e-9;

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
            parameters.add("regularization_cutoff", p3mRegularizationCutoff);
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

        /**
         * @brief Read-only view over the full storage of @p attribute.
         *
         * Overlap layouts store ghost particles after the local ones, beyond the local range that
         * getView() covers.
         */
        template <typename Attribute>
        typename Attribute::view_type::const_type storageView(const Attribute& attribute) {
            return typename Attribute::view_type::const_type(
                    attribute.getView().data(), attribute.size());
        }
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
        using VectorView        = typename ParticleContainer::particle_position_type::view_type;

        /** @brief Adds the electrostatic short-range field to @c particles.E. */
        void apply(ParticleContainer& particles) const {
            const double cutoff = requireCutoff(particles, "P3MShortRangeInteraction::apply");
            solvePairs(particles, particles.E, particles.R, cutoff);
        }

        /**
         * @brief Adds the short-range lab-frame E and B of every source moving with its own
         * velocity.
         *
         * The Ewald split is taken in the rest frame of each source j with momentum
         * u_j = gamma_j beta_j. For a target i at the lab-frame separation r from j at equal lab
         * time, the rest-frame distance is s^2 = r^2 + (u_j . r)^2, and the lab fields are
         * E = q_j gamma_j g(alpha s) r / (4 pi epsilon_0 s^3), with
         * g(x) = erfc(x) + 2 x exp(-x^2) / sqrt(pi), and B = u_j x E / (gamma_j c): the
         * short-range part of the field of a uniformly moving charge (Jackson, Classical
         * Electrodynamics, 3rd ed., Sec. 11.10). The binned mesh stage supplies the long-range
         * part with the bin velocity. At rest the field reduces to apply(), and with a common
         * momentum along z to the bin rest-frame field (dx, dy, gamma dz) Lorentz-converted to
         * the lab.
         *
         * Pairs with s >= cutoff are skipped; since s >= |r|, the lab-frame search cells and halos
         * contain every remaining pair. Positions, momenta, and charges are read over the full
         * storage, ghosts included, and must be in the frame the layout was built in. The fields
         * are added to the local particles' E and B.
         */
        static void addMovingSourceField(ParticleContainer& particles) {
            const double cutoff =
                    requireCutoff(particles, "P3MShortRangeInteraction::addMovingSourceField");
            using Vector = typename VectorView::value_type;
            using Charge = detail::P3MChargeView<typename ParticleContainer::qm_view_type>;

            static IpplTimings::TimerRef timer = IpplTimings::getTimer("P3M: PP moving src");
            IpplTimings::startTimer(timer);
            const auto positions = detail::storageView(particles.R);
            const auto momenta   = detail::storageView(particles.P);
            const Charge charge(
                    particles.getQStorageView(),
                    particles.getQMStorageMode() == ParticleContainer::QMStorageMode::Attributes);
            const auto electric         = particles.E.getView();
            const auto magnetic         = particles.B.getView();
            const double cutoffSquared  = cutoff * cutoff;
            const double alpha          = 2.0 / cutoff;
            const double coulomb        = 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
            const double regularization = p3mRegularizationCutoff;
            const double inverseLight   = 1.0 / Physics::c;
            const double twoOverSqrtPi  = 2.0 / std::sqrt(Physics::pi);
            particles.getP3MLayout().template forEachPair<Kokkos::DefaultExecutionSpace>(
                    KOKKOS_LAMBDA(const std::size_t& i, const std::size_t& j) {
                        const Vector separation = positions(i) - positions(j);
                        const double labSquared = separation.dot(separation);
                        if (labSquared >= cutoffSquared || labSquared == 0.0) {
                            return;
                        }
                        const Vector momentum    = momenta(j);
                        const double along       = momentum.dot(separation);
                        const double restSquared = labSquared + along * along;
                        if (restSquared >= cutoffSquared) {
                            return;
                        }
                        const double rest     = Kokkos::sqrt(restSquared);
                        const double x        = alpha * rest;
                        const double gaussian = twoOverSqrtPi * x * Kokkos::exp(-x * x);
                        const double cube     = restSquared * rest;
                        // Same linear core as the IPPL pair kernel, continuous at the radius.
                        const double radial =
                                rest < regularization
                                        ? 1.0 / (regularization * regularization * regularization)
                                                  - (Kokkos::erf(x) - gaussian) / cube
                                        : (Kokkos::erfc(x) + gaussian) / cube;
                        const double gamma = Kokkos::sqrt(1.0 + momentum.dot(momentum));
                        const Vector e     = (coulomb * charge(j) * gamma * radial) * separation;
                        const Vector b     = (inverseLight / gamma) * cross(momentum, e);
                        for (unsigned d = 0; d < 3; ++d) {
                            Kokkos::atomic_add(&electric(i)[d], e[d]);
                            Kokkos::atomic_add(&magnetic(i)[d], b[d]);
                        }
                    });
            IpplTimings::stopTimer(timer);
        }

    private:
        static double requireCutoff(const ParticleContainer& particles, const char* caller) {
            if (!particles.hasP3MLayout()) {
                throw OpalException(caller, "P3M requires ParticleSpatialOverlapLayout.");
            }
            const double cutoff = particles.getP3MLayout().getCutoff();
            if (!(cutoff > 0.0)) {
                throw OpalException(caller, "The P3M cutoff radius must be positive.");
            }
            return cutoff;
        }

        template <typename VectorAttribute>
        static void solvePairs(
                ParticleContainer& particles, VectorAttribute& field,
                const VectorAttribute& positions, double cutoff) {
            using ContainerView = detail::P3MContainerView<ParticleContainer>;
            using ChargeView    = detail::P3MChargeView<typename ParticleContainer::qm_view_type>;
            using Interaction   = ippl::TruncatedGreenParticleInteraction<
                      ContainerView, VectorAttribute, ChargeView>;

            ContainerView container(particles);
            ChargeView charge(
                    particles.getQStorageView(),
                    particles.getQMStorageMode() == ParticleContainer::QMStorageMode::Attributes);

            ippl::ParameterList parameters;
            parameters.add("rcut", cutoff);
            parameters.add("alpha", 2.0 / cutoff);
            // The mesh path scatters charge divided by epsilon_0. This particle path uses raw
            // charge, so its Coulomb coefficient carries epsilon_0 explicitly.
            parameters.add("force_constant", -1.0 / (4.0 * Physics::pi * Physics::epsilon_0));
            parameters.add("regularization_cutoff", p3mRegularizationCutoff);

            Interaction interaction(container, field, positions, charge, parameters);
            interaction.solve();
        }
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_P3M_ADAPTERS_H
