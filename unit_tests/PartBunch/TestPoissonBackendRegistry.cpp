#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "Ippl.h"
#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/Solvers/PoissonBackendRegistry.hpp"
#include "Utilities/OpalException.h"

namespace {
    using Registry = PoissonBackendRegistry<double, 3>;

    class PoissonBackendRegistryTest : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }

        static void TearDownTestSuite() { ippl::finalize(); }

        std::unique_ptr<FieldContainer<double, 3>> makeFieldContainer(const std::string& solver) {
            ippl::Vector<int, 3> nr        = 8;
            ippl::Vector<double, 3> rmin   = -4.0;
            ippl::Vector<double, 3> rmax   = 4.0;
            ippl::Vector<double, 3> origin = rmin;
            ippl::Vector<double, 3> hr     = (rmax - rmin) / ippl::Vector<double, 3>(nr);
            std::array<bool, 3> decomp     = {true, true, true};

            ippl::NDIndex<3> domain;
            for (unsigned i = 0; i < 3; i++) {
                domain[i] = ippl::Index(nr[i]);
            }

            auto fields = std::make_unique<FieldContainer<double, 3>>(
                    hr, rmin, rmax, decomp, domain, origin, true);
            fields->initializeFields(Registry::capabilitiesFor(solver), solver);
            return fields;
        }
    };

    TEST_F(PoissonBackendRegistryTest, ReportsSupportedSolverNames) {
        const auto names = Registry::solverNames();

        EXPECT_EQ(names.size(), 4u);
        EXPECT_NE(std::find(names.begin(), names.end(), "NONE"), names.end());
        EXPECT_NE(std::find(names.begin(), names.end(), "FFT"), names.end());
        EXPECT_NE(std::find(names.begin(), names.end(), "OPEN"), names.end());
        EXPECT_NE(std::find(names.begin(), names.end(), "CG"), names.end());

        EXPECT_TRUE(Registry::supports("NONE"));
        EXPECT_TRUE(Registry::supports("OPEN"));
        EXPECT_FALSE(Registry::supports("P3M"));
    }

    TEST_F(PoissonBackendRegistryTest, ReportsCapabilities) {
        const auto none = Registry::capabilitiesFor("NONE");
        EXPECT_FALSE(none.supportsShiftedGreens);
        EXPECT_FALSE(none.providesGradient);

        const auto fft = Registry::capabilitiesFor("FFT");
        EXPECT_FALSE(fft.supportsShiftedGreens);
        EXPECT_TRUE(fft.providesGradient);
        EXPECT_FALSE(fft.requiresPotentialBCs);

        const auto open = Registry::capabilitiesFor("OPEN");
        EXPECT_TRUE(open.supportsShiftedGreens);
        EXPECT_TRUE(open.providesPotential);
        EXPECT_TRUE(open.providesGradient);

        const auto cg = Registry::capabilitiesFor("CG");
        EXPECT_FALSE(cg.supportsShiftedGreens);
        EXPECT_TRUE(cg.requiresPotentialBCs);
    }

    TEST_F(PoissonBackendRegistryTest, RejectsUnknownSolverName) {
        Solver_t<double, 3> solverVariant;

        EXPECT_THROW(Registry::capabilitiesFor("P3M"), OpalException);
        EXPECT_THROW(Registry::couplingConstantFor("P3M"), OpalException);
        EXPECT_THROW(
                Registry::create("P3M", solverVariant, nullptr, nullptr, nullptr),
                OpalException);
    }

    TEST_F(PoissonBackendRegistryTest, CreatesNullBackend) {
        auto fields = makeFieldContainer("NONE");
        Solver_t<double, 3> solverVariant;

        auto backend = Registry::create(
                "NONE", solverVariant, &fields->getRho(), &fields->getE(), &fields->getPhi());

        ASSERT_NE(backend, nullptr);
        EXPECT_EQ(backend->name(), "NONE");
        EXPECT_EQ(backend->capabilities().isNoOp, Registry::capabilitiesFor("NONE").isNoOp);
        EXPECT_DOUBLE_EQ(backend->couplingConstant(), Registry::couplingConstantFor("NONE"));
        EXPECT_FALSE(backend->capabilities().supportsShiftedGreens);
        EXPECT_NO_THROW(backend->solve(SolveRequest<double, 3>{}));
    }
}  // namespace
