#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <tps/tps.hpp>
#include "Ippl.h"
#include "Utilities/CLRangeError.h"
#include "Utilities/DivideError.h"
#include "Utilities/DomainError.h"
#include "Utilities/OpalException.h"

namespace {
    class TPSIntegration : public ::testing::Test {
    protected:
        static void SetUpTestSuite() {
            int argc    = 0;
            char** argv = nullptr;
            ippl::initialize(argc, argv);
        }
        static void TearDownTestSuite() { ippl::finalize(); }
    };

    tps::PhaseSpaceChart chart(const char* name) {
        tps::PhaseSpaceChart result;
        result.name           = name;
        result.coordinates    = {{"q", "m", 1.0}, {"p", "normalized", 1.0}};
        result.reference      = {0.0, 0.0};
        result.frame.name     = "lab";
        result.convention     = tps::CoordinateConvention::canonical;
        result.canonicalPairs = {{0, 1}};
        return result;
    }

    // NVCC extended device lambdas cannot be enclosed by GoogleTest's private TestBody.
    // Keep the real device launch in a namespace-scope helper on every backend.
    template <class MapView, class ValuesView, class StatusView>
    void evaluatePreparedParticles(
            MapView map, ValuesView input, ValuesView output, ValuesView workspace,
            StatusView status) {
        const int particles = input.extent_int(0);
        Kokkos::parallel_for(
                "TPS OPALX rank tracking", particles, KOKKOS_LAMBDA(int i) {
                    status(i) = tps::evaluate_one(
                            map, tps::deviceSpan(Kokkos::subview(input, i, Kokkos::ALL)),
                            tps::deviceSpan(Kokkos::subview(output, i, Kokkos::ALL)),
                            tps::deviceSpan(Kokkos::subview(workspace, i, Kokkos::ALL)));
                });
        Kokkos::fence("TPS OPALX rank tracking complete");
    }

    TEST_F(TPSIntegration, HostCompositionAndPreparedTrackingPerRank) {
        using Series = tps::Tps<double>;
        // Production use configures the order before serialized host construction.
        Series::setGlobalTruncOrder(3);
        const auto entrance = chart("entrance");
        const auto middle   = chart("middle");
        const auto exit     = chart("exit");
        const auto q        = Series::makeVariable(2, 0);
        const auto p        = Series::makeVariable(2, 1);
        const std::array<Series, 2> driftTerms{q + 2.0 * p, p};
        const std::array<Series, 2> kickTerms{q, p - 3.0 * q * q};
        const tps::Map<double> drift(tps::MapContext(entrance, middle, 3, {"drift"}), driftTerms);
        const tps::Map<double> kick(tps::MapContext(middle, exit, 3, {"kick"}), kickTerms);
        const auto combined = tps::compose(kick, drift);
        const std::array<double, 2> point{0.1, 0.02};
        const auto host = combined.evaluate(point);
        EXPECT_NEAR(host[0], 0.14, 2e-13);
        EXPECT_NEAR(host[1], 0.02 - 3.0 * 0.14 * 0.14, 2e-13);
        const auto jacobian = combined.jacobian(point);
        EXPECT_NEAR(jacobian(0, 0), 1.0, 2e-13);
        EXPECT_NEAR(jacobian(0, 1), 2.0, 2e-13);
        EXPECT_NEAR(jacobian(1, 0), -6.0 * 0.14, 2e-13);
        EXPECT_NEAR(jacobian(1, 1), 1.0 - 12.0 * 0.14, 2e-13);

        const tps::PreparedMap<double> prepared(combined);
        const auto map = prepared.device_view();
        EXPECT_EQ(Series::getGlobalTruncOrder(), 3);
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        constexpr int particles = 7;
        Kokkos::View<double**, Kokkos::LayoutRight> input("TPS rank input", particles, 2);
        Kokkos::View<double**, Kokkos::LayoutRight> output("TPS rank output", particles, 2);
        Kokkos::View<double**, Kokkos::LayoutRight> workspace(
                "TPS rank workspace", particles, map.workspace_size());
        Kokkos::View<tps::EvaluationStatus*> status("TPS rank status", particles);
        auto hostInput = Kokkos::create_mirror_view(input);
        for (int i = 0; i < particles; ++i) {
            hostInput(i, 0) = 0.1 + 0.01 * rank + 0.001 * i;
            hostInput(i, 1) = 0.02;
        }
        Kokkos::deep_copy(input, hostInput);
        evaluatePreparedParticles(map, input, output, workspace, status);
        const auto hostOutput = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, output);
        const auto hostStatus = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, status);
        int localFailures     = 0;
        for (int i = 0; i < particles; ++i) {
            const double expectedQ = hostInput(i, 0) + 2.0 * hostInput(i, 1);
            const double expectedP = hostInput(i, 1) - 3.0 * expectedQ * expectedQ;
            EXPECT_EQ(hostStatus(i), tps::EvaluationStatus::success);
            EXPECT_NEAR(hostOutput(i, 0), expectedQ, 2e-13);
            EXPECT_NEAR(hostOutput(i, 1), expectedP, 2e-13);
            if (hostStatus(i) != tps::EvaluationStatus::success
                || std::abs(hostOutput(i, 0) - expectedQ) > 2e-13
                || std::abs(hostOutput(i, 1) - expectedP) > 2e-13)
                ++localFailures;
        }
        int globalFailures = 0;
        MPI_Allreduce(&localFailures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        EXPECT_EQ(globalFailures, 0);
    }

    TEST_F(TPSIntegration, AlgebraFailuresUseOpalExceptions) {
        const auto x = tps::Tps<double>::makeVariable(1, 0);
        EXPECT_THROW((void)(x / 0.0), ::DivideError);
        EXPECT_THROW((void)(x / 0.0), ::OpalException);
        EXPECT_THROW((void)x[-1], ::CLRangeError);
        EXPECT_THROW((void)tps::sqrt(x - 1.0), ::DomainError);
        try {
            (void)x[-1];
            FAIL() << "Invalid coefficient index was accepted";
        } catch (const ::OpalException& error) {
            EXPECT_EQ(error.where(), "Tps::getCoefficients()");
            EXPECT_EQ(error.what(), "Monomial index out of range.");
        }
    }
    TEST_F(TPSIntegration, NonlinearSubstitutionSurvivesFreshCacheGrowth) {
        using Series = tps::Tps<double>;
        Series::setGlobalTruncOrder(6);
        for (int variables : {3, 4, 5}) {
            SCOPED_TRACE(variables);
            const auto x      = Series::makeVariable(variables, 0);
            const auto y      = Series::makeVariable(variables, 1);
            const auto source = Series::makeVarPower(variables, 0, 3) + 2.0 * y + 1.0;
            tps::VpsMap<double> inner(variables, variables);
            for (int i = 0; i < variables; ++i)
                inner[i] = Series::makeVariable(variables, i);
            inner[0] = x * x;
            // Do not build the degree-six oracle before substitution: a fresh
            // degree-three traversal must survive product growth to degree six.
            const auto result = source.substitute(inner);
            tps::Vector<double> point(variables, .2);
            point[0] = .5;
            EXPECT_NEAR(result.evaluate(point), std::pow(.5, 6) + 1.4, 2e-13);
            const auto expected = Series::makeVarPower(variables, 0, 6) + 2.0 * y + 1.0;
            EXPECT_EQ(result.getMaxOrder(), 6);
            for (int i = 0; i < expected.getSize(); ++i)
                EXPECT_DOUBLE_EQ(result.getCoefficient(i), expected.getCoefficient(i));
        }
        Series::setGlobalTruncOrder(Series::EXACT);
    }
    TEST_F(TPSIntegration, ScalarDivisionPreservesSharedMapAndAliasedDivisor) {
        using Series = tps::Tps<double>;
        Series::setGlobalTruncOrder(Series::EXACT);
        const auto original = tps::identityMap<double>(chart("division"), 2);
        auto component      = original.polynomial()[0];
        component /= 2.0;
        tps::Vector<double> point(2, 0.0);
        point[0] = 2.0;
        EXPECT_DOUBLE_EQ(component.evaluate(point), 1.0);
        EXPECT_DOUBLE_EQ(original.evaluate(std::array<double, 2>{2.0, 0.0})[0], 2.0);
        const auto divided = original.polynomial()[0] / 2.0;
        EXPECT_DOUBLE_EQ(divided.evaluate(point), 1.0);
        EXPECT_DOUBLE_EQ(original.evaluate(std::array<double, 2>{2.0, 0.0})[0], 2.0);
        auto value            = 2.0 + 4.0 * Series::makeVariable(2, 0);
        const double& divisor = static_cast<const Series&>(value)[0];
        value /= divisor;
        EXPECT_DOUBLE_EQ(value.getCoefficient(0), 1.0);
        EXPECT_DOUBLE_EQ(value.getCoefficient(1), 2.0);
        const auto unchanged = value;
        EXPECT_THROW(value /= 0.0, ::DivideError);
        EXPECT_EQ(value, unchanged);
    }
    TEST_F(TPSIntegration, DefaultMapSubstitutionResolvesExactToFiniteDegree) {
        using Series = tps::Tps<double>;
        Series::setGlobalTruncOrder(Series::EXACT);
        const auto x = Series::makeVariable(2, 0), y = Series::makeVariable(2, 1);
        tps::VpsMap<double> identity(2, 2);
        identity[0]     = x;
        identity[1]     = y;
        const auto same = identity.substitute(identity);
        EXPECT_EQ(same[0], x);
        EXPECT_EQ(same[1], y);
        tps::VpsMap<double> outer(2, 2), inner(2, 2);
        outer[0]                 = x * x + y;
        outer[1]                 = y;
        inner[0]                 = x * x;
        inner[1]                 = y;
        const auto composed      = outer.substitute(inner);
        const auto explicitOrder = outer.substitute(inner, 4);
        EXPECT_EQ(composed[0], explicitOrder[0]);
        tps::Vector<double> point(2, .2);
        point[0] = .5;
        EXPECT_NEAR(composed.constantTerm(point)[0], std::pow(.5, 4) + .2, 2e-13);
        // A rectangular inner map has two outputs and three independent inputs.
        const auto u = Series::makeVariable(3, 0), v = Series::makeVariable(3, 1);
        tps::VpsMap<double> rectangular(2, 3);
        rectangular[0]       = u * u;
        rectangular[1]       = v;
        const auto rectangle = outer.substitute(rectangular);
        tps::Vector<double> point3(3, .2);
        point3[0] = .5;
        EXPECT_NEAR(rectangle.constantTerm(point3)[0], std::pow(.5, 4) + .2, 2e-13);
        tps::VpsMap<double> constant(2, 2);
        constant[0][0]   = 3;
        constant[1][0]   = -1;
        const auto value = outer.substitute(constant);
        EXPECT_DOUBLE_EQ(value[0].getCoefficient(0), 8);
        const auto constantOuter = constant.substitute(identity);
        EXPECT_DOUBLE_EQ(constantOuter[0].getCoefficient(0), 3);
        EXPECT_THROW(identity.substitute(identity, -1), ::CLRangeError);
        Series::setGlobalTruncOrder(2);
        const auto truncated = outer.substitute(inner);
        EXPECT_LE(truncated[0].getMaxOrder(), 2);
        Series::setGlobalTruncOrder(Series::EXACT);
    }
}  // namespace
