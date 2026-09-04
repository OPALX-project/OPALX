#include <gtest/gtest.h>

#include "PartBunch/CartesianDomain.h"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"
#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace opalx::spacecharge {
    namespace {

        class PoissonSolverTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            CartesianDomainConfig3D storage(bool periodic = false) const {
                CartesianDomainConfig3D result;
                result.meshSize                 = {8, 8, 8};
                result.decomposition            = {false, false, false};
                result.periodicParticleBoundary = periodic;
                return result;
            }

            void exerciseBackend(
                    PoissonSolverType kind, FieldBoundaryCondition boundary,
                    double cutoff = 0.0) const {
                auto setup = storage(boundary == FieldBoundaryCondition::Periodic);
                CartesianDomain<double, 3> domain(setup);
                CartesianPICFieldStorage<double, 3> fields(domain);
                fields.initializeFields(kind);
                PoissonSolverConfig config;
                config.type               = kind;
                config.p3mCutoff          = cutoff;
                config.boundaryConditions = {boundary, boundary, boundary};
                PoissonSolver solver(config, {&fields.chargeDensity(), &fields.electricField()});
                EXPECT_NO_THROW(solver.warmup());
                EXPECT_NO_THROW(solver.solve({}, {.suppressFieldDump = true}));
            }
        };

        TEST_F(PoissonSolverTest, RebuildsAfterLayoutChanges) {
            auto setup = storage(true);
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> fields(domain);
            fields.initializeFields(PoissonSolverType::PeriodicFFT);
            PoissonSolverConfig config;
            config.type               = PoissonSolverType::PeriodicFFT;
            config.boundaryConditions = {
                    FieldBoundaryCondition::Periodic, FieldBoundaryCondition::Periodic,
                    FieldBoundaryCondition::Periodic};
            PoissonSolver solver(config, {&fields.chargeDensity(), &fields.electricField()});

            for (const auto& extents :
                 {std::array<std::size_t, 3>{8, 8, 16}, std::array<std::size_t, 3>{8, 8, 8}}) {
                ASSERT_TRUE(domain.rebuildGlobalLayoutInPlace(
                        extents, std::array<bool, 3>{false, false, true}));
                fields.updateFieldLayoutsAfterLayoutChange();
                solver.rebuildAfterLayoutChange({&fields.chargeDensity(), &fields.electricField()});
                EXPECT_NO_THROW(solver.solve({}, {.suppressFieldDump = true}));
                EXPECT_EQ(fields.layoutExtents(), extents);
            }
        }

        TEST_F(PoissonSolverTest, DispatchesEveryImplementedBackend) {
            exerciseBackend(PoissonSolverType::None, FieldBoundaryCondition::Open);
            exerciseBackend(PoissonSolverType::PeriodicFFT, FieldBoundaryCondition::Periodic);
            exerciseBackend(PoissonSolverType::Open, FieldBoundaryCondition::Open);
            exerciseBackend(PoissonSolverType::P3M, FieldBoundaryCondition::Open, 0.25);
            exerciseBackend(PoissonSolverType::P3M, FieldBoundaryCondition::Periodic, 0.25);
        }

        TEST_F(PoissonSolverTest, RestoresStandardKernelAfterShiftedSolve) {
            auto setup = storage();
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> fields(domain);
            fields.initializeFields(PoissonSolverType::Open);
            PoissonSolverConfig config;
            config.type = PoissonSolverType::Open;
            PoissonSolver solver(config, {&fields.chargeDensity(), &fields.electricField()});
            solver.warmup();

            const auto fillCharge = [&fields] {
                auto charge = fields.chargeDensity().getHostMirror();
                for (std::size_t i = 0; i < charge.extent(0); ++i) {
                    for (std::size_t j = 0; j < charge.extent(1); ++j) {
                        for (std::size_t k = 0; k < charge.extent(2); ++k) {
                            charge(i, j, k) = static_cast<double>((i + 1) * (j + 2) + 3 * k);
                        }
                    }
                }
                Kokkos::deep_copy(fields.chargeDensity().getView(), charge);
            };

            fillCharge();
            solver.solve({}, {.suppressFieldDump = true});
            const auto standard = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), fields.electricField().getView());

            fillCharge();
            PoissonSolveRequest shifted;
            shifted.greenFunctionShift = Vector_t<double, 3>(0.0, 0.0, 0.25);
            solver.solve(shifted, {.suppressFieldDump = true});
            fillCharge();
            solver.solve({}, {.suppressFieldDump = true});
            const auto restored = Kokkos::create_mirror_view_and_copy(
                    Kokkos::HostSpace(), fields.electricField().getView());

            for (std::size_t i = 0; i < standard.extent(0); ++i) {
                for (std::size_t j = 0; j < standard.extent(1); ++j) {
                    for (std::size_t k = 0; k < standard.extent(2); ++k) {
                        for (unsigned dimension = 0; dimension < 3; ++dimension) {
                            const double expected = standard(i, j, k)[dimension];
                            EXPECT_NEAR(
                                    restored(i, j, k)[dimension], expected,
                                    1.0e-12 * std::max(1.0, std::abs(expected)));
                        }
                    }
                }
            }
        }

        TEST_F(PoissonSolverTest, ExposesMetadataAndRejectsInvalidBackends) {
            auto setup = storage();
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> fields(domain);
            fields.initializeFields(PoissonSolverType::None);
            const PoissonFieldBinding binding{&fields.chargeDensity(), &fields.electricField()};

            PoissonSolverConfig config;
            config.type = PoissonSolverType::None;
            PoissonSolver none(config, binding);
            EXPECT_EQ(none.name(), "NONE");
            EXPECT_TRUE(none.capabilities().isNoOp);
            EXPECT_DOUBLE_EQ(
                    none.couplingConstant(), 1.0 / (4.0 * Physics::pi * Physics::epsilon_0));

            config.type = PoissonSolverType::ConjugateGradient;
            EXPECT_THROW(PoissonSolver(config, binding), OpalException);
            config.type = static_cast<PoissonSolverType>(255);
            EXPECT_THROW(PoissonSolver(config, binding), OpalException);
        }

        TEST_F(PoissonSolverTest, RejectsShiftForNonOpenBackend) {
            auto setup = storage();
            CartesianDomain<double, 3> domain(setup);
            CartesianPICFieldStorage<double, 3> fields(domain);
            fields.initializeFields(PoissonSolverType::None);
            PoissonSolver solver({}, {&fields.chargeDensity(), &fields.electricField()});
            PoissonSolveRequest request;
            request.greenFunctionShift = Vector_t<double, 3>(0.0, 0.0, 0.25);
            EXPECT_THROW(solver.solve(request, {.suppressFieldDump = true}), OpalException);
        }

    }  // namespace
}  // namespace opalx::spacecharge
