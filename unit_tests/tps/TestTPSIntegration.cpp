#include "Ippl.h"
#include "Utilities/DivideError.h"
#include <tps/tps.hpp>
#include <gtest/gtest.h>
#include <array>
#include <cmath>

namespace {
class TPSIntegration : public ::testing::Test {
  protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }
    static void TearDownTestSuite() {
        ippl::finalize();
    }
};

tps::PhaseSpaceChart chart(const char* name) {
    tps::PhaseSpaceChart result;
    result.name = name;
    result.coordinates = {{"q", "m", 1.0}, {"p", "normalized", 1.0}};
    result.reference = {0.0, 0.0};
    result.frame.name = "lab";
    result.convention = tps::CoordinateConvention::canonical;
    result.canonicalPairs = {{0, 1}};
    return result;
}

// NVCC extended device lambdas cannot be enclosed by GoogleTest's private TestBody.
// Keep the real device launch in a namespace-scope helper on every backend.
template <class MapView, class ValuesView, class StatusView>
void evaluatePreparedParticles(MapView map, ValuesView input, ValuesView output,
                               ValuesView workspace, StatusView status) {
    const int particles = input.extent_int(0);
    Kokkos::parallel_for(
        "TPS OPALX rank tracking", particles, KOKKOS_LAMBDA(int i) {
            status(i) =
                tps::evaluate_one(map, tps::deviceSpan(Kokkos::subview(input, i, Kokkos::ALL)),
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
    const auto middle = chart("middle");
    const auto exit = chart("exit");
    const auto q = Series::makeVariable(2, 0);
    const auto p = Series::makeVariable(2, 1);
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
    Kokkos::View<double**, Kokkos::LayoutRight> workspace("TPS rank workspace", particles,
                                                          map.workspace_size());
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
    int localFailures = 0;
    for (int i = 0; i < particles; ++i) {
        const double expectedQ = hostInput(i, 0) + 2.0 * hostInput(i, 1);
        const double expectedP = hostInput(i, 1) - 3.0 * expectedQ * expectedQ;
        EXPECT_EQ(hostStatus(i), tps::EvaluationStatus::success);
        EXPECT_NEAR(hostOutput(i, 0), expectedQ, 2e-13);
        EXPECT_NEAR(hostOutput(i, 1), expectedP, 2e-13);
        if (hostStatus(i) != tps::EvaluationStatus::success ||
            std::abs(hostOutput(i, 0) - expectedQ) > 2e-13 ||
            std::abs(hostOutput(i, 1) - expectedP) > 2e-13)
            ++localFailures;
    }
    int globalFailures = 0;
    MPI_Allreduce(&localFailures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalFailures, 0);
}

TEST_F(TPSIntegration, ExceptionFamiliesCoexist) {
    EXPECT_THROW(throw ::DivideError("OPALX integration"), ::ArithmeticError);
    EXPECT_THROW(throw tps::DivideError("TPS integration"), tps::Error);
}
} // namespace
