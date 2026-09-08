#include <gtest/gtest.h>

#include "PartBunch/CartesianDomain.h"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC3D/CartesianPIC3DFieldStorage.h"
#include "SpaceCharge/Poisson/NullPoissonAdapter.h"
#include "SpaceCharge/Poisson/OpenPoissonAdapter.h"
#include "SpaceCharge/Poisson/P3MAdapters.h"
#include "SpaceCharge/Poisson/PeriodicPoissonAdapter.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace opalx::spacecharge {
    namespace {

        using Fields = CartesianPIC3DFieldStorage<double, 3>;

        struct BackendCase {
            PoissonSolverType type;
            FieldBoundaryCondition boundary;
            GreenFunctionType green = GreenFunctionType::Integrated;

            PoissonSolverConfig config() const {
                PoissonSolverConfig result;
                result.type               = type;
                result.boundaryConditions = {boundary, boundary, boundary};
                result.greenFunction      = green;
                result.p3mCutoff          = type == PoissonSolverType::P3M ? 0.25 : 0.0;
                return result;
            }
        };

        constexpr std::array cases{
                BackendCase{PoissonSolverType::None, FieldBoundaryCondition::Open},
                BackendCase{PoissonSolverType::None, FieldBoundaryCondition::Periodic},
                BackendCase{PoissonSolverType::PeriodicFFT, FieldBoundaryCondition::Periodic},
                BackendCase{
                        PoissonSolverType::Open, FieldBoundaryCondition::Open,
                        GreenFunctionType::Standard},
                BackendCase{
                        PoissonSolverType::Open, FieldBoundaryCondition::Open,
                        GreenFunctionType::Integrated},
                BackendCase{PoissonSolverType::P3M, FieldBoundaryCondition::Open},
                BackendCase{PoissonSolverType::P3M, FieldBoundaryCondition::Periodic}};

        class PoissonSolverTest : public ::testing::Test {
        protected:
            static void SetUpTestSuite() {
                int argc    = 0;
                char** argv = nullptr;
                ippl::initialize(argc, argv);
            }

            static void TearDownTestSuite() { ippl::finalize(); }

            static CartesianDomainConfig3D storage(const BackendCase& selected) {
                CartesianDomainConfig3D result;
                result.meshSize      = {8, 8, 8};
                result.decomposition = {false, false, false};
                result.periodicParticleBoundary =
                        selected.boundary == FieldBoundaryCondition::Periodic;
                return result;
            }

            static PoissonFieldBinding binding(Fields& fields) {
                return {&fields.chargeDensity(), &fields.electricField()};
            }

            template <typename View>
            static auto snapshot(const View& view) {
                // create_mirror_view_and_copy can alias the source on SERIAL. These comparisons
                // need an independent allocation that survives later native solves.
                auto result = Kokkos::create_mirror(view);
                Kokkos::deep_copy(result, view);
                return result;
            }

            template <typename View, typename Function>
            static void forCells(const View& view, std::size_t ghost, Function function) {
                for (std::size_t i = ghost; i < view.extent(0) - ghost; ++i) {
                    for (std::size_t j = ghost; j < view.extent(1) - ghost; ++j) {
                        for (std::size_t k = ghost; k < view.extent(2) - ghost; ++k) {
                            function(i, j, k);
                        }
                    }
                }
            }

            static void fillCharge(Fields& fields) {
                auto charge = fields.chargeDensity().getHostMirror();
                forCells(charge, 0, [&](auto i, auto j, auto k) {
                    charge(i, j, k) = static_cast<double>((i + 1) * (j + 2) + 3 * k);
                });
                Kokkos::deep_copy(fields.chargeDensity().getView(), charge);
                fields.electricField() = 0.0;
            }

            template <typename View>
            static void expectElectricField(
                    const Fields& fields, const View& expected, bool includeGhosts = false) {
                const auto actual = snapshot(fields.electricField().getView());
                ASSERT_EQ(actual.extent(0), expected.extent(0));
                ASSERT_EQ(actual.extent(1), expected.extent(1));
                ASSERT_EQ(actual.extent(2), expected.extent(2));
                const std::size_t ghost = includeGhosts ? 0 : fields.electricField().getNghost();
                forCells(actual, ghost, [&](auto i, auto j, auto k) {
                    for (unsigned d = 0; d < 3; ++d) {
                        ASSERT_TRUE(std::isfinite(actual(i, j, k)[d]));
                        EXPECT_NEAR(
                                actual(i, j, k)[d], expected(i, j, k)[d],
                                1.0e-12 * std::max(1.0, std::abs(expected(i, j, k)[d])));
                    }
                });
            }
        };

        TEST_F(PoissonSolverTest, FactoryPreservesMetadataAndSolvesEveryBackend) {
            for (const auto& selected : cases) {
                SCOPED_TRACE(static_cast<int>(selected.type));
                CartesianDomain<double, 3> domain(storage(selected));
                Fields fields(domain);
                fields.initializeFields(selected.type);
                auto solver         = makePoissonSolver(selected.config(), binding(fields));
                const bool none     = selected.type == PoissonSolverType::None;
                const bool periodic = selected.type == PoissonSolverType::PeriodicFFT;
                const bool open     = selected.type == PoissonSolverType::Open;
                EXPECT_EQ(solver->name(), none ? "NONE" : periodic ? "FFT" : open ? "OPEN" : "P3M");
                const auto& caps = solver->capabilities();
                EXPECT_EQ(caps.isNoOp, none);
                EXPECT_EQ(caps.supportsShiftedGreenFunction, open);
                EXPECT_TRUE(caps.normalizeChargeByCellVolume);
                EXPECT_EQ(caps.subtractNeutralizingBackground, none || periodic);
                EXPECT_EQ(caps.debugDumpChargeBeforeSolve, !none);
                EXPECT_EQ(caps.debugDumpScalarAfterSolve, !none);
                EXPECT_EQ(caps.debugDumpVectorAfterSolve, !none);
                EXPECT_DOUBLE_EQ(
                        solver->couplingConstant(),
                        none ? 1.0 / (4.0 * Physics::pi * Physics::epsilon_0)
                             : 1.0 / Physics::epsilon_0);
                EXPECT_NO_THROW(solver->warmup());
                fillCharge(fields);
                solver->solve({}, {.suppressFieldDump = true});
                const auto expected = snapshot(fields.electricField().getView());
                double norm         = 0.0;
                forCells(expected, fields.electricField().getNghost(), [&](auto i, auto j, auto k) {
                    for (unsigned d = 0; d < 3; ++d)
                        norm = std::max(norm, std::abs(expected(i, j, k)[d]));
                });
                if (none) {
                    EXPECT_DOUBLE_EQ(norm, 0.0);
                } else {
                    EXPECT_GT(norm, 1.0e-8);
                }
                solver->warmup();
                fillCharge(fields);
                solver->solve({}, {.suppressFieldDump = true});
                expectElectricField(fields, expected);
            }
        }

        TEST_F(PoissonSolverTest, RestoresOrdinaryKernelAfterShiftedSolve) {
            for (auto green : {GreenFunctionType::Standard, GreenFunctionType::Integrated}) {
                SCOPED_TRACE(static_cast<int>(green));
                const BackendCase selected{
                        PoissonSolverType::Open, FieldBoundaryCondition::Open, green};
                CartesianDomain<double, 3> domain(storage(selected));
                Fields fields(domain);
                fields.initializeFields(selected.type);
                auto solver = makePoissonSolver(selected.config(), binding(fields));
                solver->warmup();
                fillCharge(fields);
                solver->solve({}, {.suppressFieldDump = true});
                const auto ordinary = snapshot(fields.electricField().getView());

                fillCharge(fields);
                PoissonSolveRequest shifted;
                shifted.greenFunctionShift = Vector_t<double, 3>(0.0, 0.0, 0.25);
                solver->solve(shifted, {.suppressFieldDump = true});
                const auto image = snapshot(fields.electricField().getView());
                double change    = 0.0;
                forCells(image, fields.electricField().getNghost(), [&](auto i, auto j, auto k) {
                    for (unsigned d = 0; d < 3; ++d) {
                        ASSERT_TRUE(std::isfinite(image(i, j, k)[d]));
                        change = std::max(
                                change, std::abs(image(i, j, k)[d] - ordinary(i, j, k)[d]));
                    }
                });
                EXPECT_GT(change, 1.0e-8);
                fillCharge(fields);
                solver->solve({}, {.suppressFieldDump = true});
                expectElectricField(fields, ordinary);
            }
        }

        TEST_F(PoissonSolverTest, RebuiltBackendMatchesFreshSolverOnChangedLayout) {
            for (const auto& selected : cases) {
                SCOPED_TRACE(static_cast<int>(selected.type));
                CartesianDomain<double, 3> domain(storage(selected));
                Fields fields(domain);
                fields.initializeFields(selected.type);
                auto solver = makePoissonSolver(selected.config(), binding(fields));
                solver->warmup();

                for (const auto& extents :
                     {std::array<std::size_t, 3>{8, 8, 16}, std::array<std::size_t, 3>{8, 8, 8}}) {
                    fillCharge(fields);
                    PoissonSolveRequest request;
                    if (selected.type == PoissonSolverType::Open) {
                        request.greenFunctionShift = Vector_t<double, 3>(0.0, 0.0, 0.25);
                    }
                    solver->solve(request, {.suppressFieldDump = true});
                    ASSERT_TRUE(domain.rebuildGlobalLayoutInPlace(
                            extents, std::array<bool, 3>{false, false, true}));
                    fields.updateFieldLayoutsAfterLayoutChange();
                    solver->rebuildAfterLayoutChange(binding(fields));
                    EXPECT_EQ(fields.layoutExtents(), extents);

                    Fields freshFields(domain);
                    freshFields.initializeFields(selected.type);
                    auto fresh = makePoissonSolver(selected.config(), binding(freshFields));
                    fillCharge(fields);
                    fillCharge(freshFields);
                    solver->solve({}, {.suppressFieldDump = true});
                    fresh->solve({}, {.suppressFieldDump = true});
                    expectElectricField(fields, snapshot(freshFields.electricField().getView()));
                }
            }
        }

        TEST_F(PoissonSolverTest, RejectsUnsupportedRequestsBeforeModifyingFields) {
            for (const auto& selected : cases) {
                if (selected.type == PoissonSolverType::Open) continue;
                SCOPED_TRACE(static_cast<int>(selected.type));
                CartesianDomain<double, 3> domain(storage(selected));
                Fields fields(domain);
                fields.initializeFields(selected.type);
                auto solver = makePoissonSolver(selected.config(), binding(fields));
                fillCharge(fields);
                Kokkos::deep_copy(fields.electricField().getView(), Vector_t<double, 3>(17.0));
                const auto rho      = snapshot(fields.chargeDensity().getView());
                const auto electric = snapshot(fields.electricField().getView());
                PoissonSolveRequest request;
                request.greenFunctionShift = Vector_t<double, 3>(0.0, 0.0, 0.25);
                EXPECT_THROW(solver->solve(request, {.suppressFieldDump = true}), OpalException);
                const auto after = snapshot(fields.chargeDensity().getView());
                forCells(rho, 0, [&](auto i, auto j, auto k) {
                    EXPECT_DOUBLE_EQ(after(i, j, k), rho(i, j, k));
                });
                expectElectricField(fields, electric, true);
            }
        }

        TEST_F(PoissonSolverTest, RejectsInvalidConfigurationAndBindings) {
            CartesianDomain<double, 3> domain(storage(cases[0]));
            Fields fields(domain);
            fields.initializeFields(PoissonSolverType::None);
            PoissonSolverConfig config;
            config.type = PoissonSolverType::ConjugateGradient;
            EXPECT_THROW((void)makePoissonSolver(config, binding(fields)), OpalException);
            config.type = static_cast<PoissonSolverType>(255);
            EXPECT_THROW((void)makePoissonSolver(config, binding(fields)), OpalException);
            EXPECT_THROW(
                    (void)makePoissonSolver({}, {nullptr, &fields.electricField()}), OpalException);
            EXPECT_THROW(
                    (void)makePoissonSolver({}, {&fields.chargeDensity(), nullptr}), OpalException);
            auto solver = makePoissonSolver({}, binding(fields));
            EXPECT_THROW(solver->rebuildAfterLayoutChange({}), OpalException);
        }

        TEST_F(PoissonSolverTest, ConcreteAdaptersRejectMismatchedBackendKinds) {
            CartesianDomain<double, 3> domain(storage(cases[0]));
            Fields fields(domain);
            fields.initializeFields(PoissonSolverType::None);
            EXPECT_THROW(OpenPoissonAdapter({}, binding(fields)), OpalException);
            EXPECT_THROW(PeriodicPoissonAdapter({}, binding(fields)), OpalException);
            EXPECT_THROW(P3MMeshPoissonAdapter({}, binding(fields)), OpalException);
            PoissonSolverConfig open;
            open.type = PoissonSolverType::Open;
            EXPECT_THROW(NullPoissonAdapter(open, binding(fields)), OpalException);
        }

    }  // namespace
}  // namespace opalx::spacecharge
