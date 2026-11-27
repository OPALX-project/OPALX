//
// Unit tests for class Enge
//
// Copyright (c) 2017-2026, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.

#include <algorithm>
#include <array>
#include <climits>
#include <limits>
#include <memory>
#include <set>
#include <type_traits>

#include "AbsBeamline/EndFieldModel/Enge.h"
#include "AbstractObjects/OpalData.h"
#include "EngeDerivativeReference.h"
#include "Physics/Physics.h"
#include "gtest/gtest.h"

namespace endfieldmodel {
    struct EngeHostMirrorTestAccess {
        static void evaluate(const Enge& model, double x, int n, double* values) {
            Enge::validateDerivativeOrder(n);
            Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                    output(values, n + 1);
            std::lock_guard<std::mutex> lock(model.host_mutex_m);
            model.prepareHostConfig(n);
            Enge::functionPoint(model.host_config_m, x, n, output);
        }
    };
}  // namespace endfieldmodel

static_assert(
        std::extent_v<decltype(endfieldmodel::EngeWorkspace::enge)>
        == endfieldmodel::EngeConfig::max_derivative + 1);

namespace {
    void evaluateHighestOrderWithGuard(
            endfieldmodel::EngeConfig config, Kokkos::View<double*> result) {
        Kokkos::parallel_for(
                "Enge workspace boundary", 1, KOKKOS_LAMBDA(int) {
                    struct Guarded {
                        endfieldmodel::EngeWorkspace workspace;
                        double guard;
                    } guarded;
                    guarded.guard = 12345.0;
                    result(0)     = endfieldmodel::Enge::getEnge(
                            config, .3, endfieldmodel::EngeConfig::max_derivative,
                            guarded.workspace);
                    result(1) = guarded.guard;
                });
    }
}  // namespace

class TestEnge : public testing::Test {
public:
    TestEnge() {}

    static void SetUpTestSuite() { Kokkos::initialize(); }
    static void TearDownTestSuite() { Kokkos::finalize(); }
};

TEST_F(TestEnge, ConstructorTest) {
    // x0, lambda, max_index
    endfieldmodel::Enge enge0;
    EXPECT_EQ(enge0.getX0(), 0.0);
    EXPECT_EQ(enge0.getLambda(), 0.0);
    EXPECT_EQ(enge0.getCoefficients().size(), 0);

    std::vector<double> a = {5, 6};
    auto enge1            = std::make_shared<endfieldmodel::Enge>(a, 2, 7);
    enge1->initialise();
    EXPECT_EQ(enge1->getX0(), 2.0);
    EXPECT_EQ(enge1->getCentreLength(), 4.0);
    EXPECT_EQ(enge1->getEndLength(), 7.0);
    EXPECT_EQ(enge1->getCoefficients(), a);

    endfieldmodel::Enge* enge2 = enge1->clone();
    enge2->initialise();
    EXPECT_EQ(enge2->getX0(), enge1->getX0());
    EXPECT_EQ(enge2->getCentreLength(), enge1->getCentreLength());
    EXPECT_EQ(enge2->getCoefficients(), enge2->getCoefficients());
    delete enge2;
}

double myEnge(double x, std::vector<double> a, double lambda) {
    double deltaX = x / lambda;
    double xPow   = 1.0;
    double p      = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        p += a[i] * xPow;
        xPow *= deltaX;
    }
    double enge = 1 / (1.0 + exp(p));
    return enge;
}

TEST_F(TestEnge, DerivativeTest) {
    std::vector<double> xVector = {-5.5, -5.0, -3.0, 0.0, 3.0, 5.0, 5.5};
    double dx                   = 1e-4;
    double x0                   = 5.0;
    double lambda               = 0.5;
    std::vector<double> engeA   = {0.0, 1.0, 0.0};
    endfieldmodel::Enge enge(engeA, x0, lambda);
    enge.initialise();
    Kokkos::View<double*> xDevice("x", 1);
    auto xHost = Kokkos::create_mirror_view(xDevice);
    Kokkos::View<double**> yDevice("y", 1, 5);
    for (auto x : xVector) {
        double yTest = enge.function(x, 0);
        double yRef  = -1 + myEnge(x - x0, engeA, lambda) + myEnge(-x - x0, engeA, lambda);
        xHost(0)     = x;
        Kokkos::deep_copy(xDevice, xHost);
        enge.function(xDevice, 4, yDevice);
        auto yHost = Kokkos::create_mirror_view(yDevice);
        Kokkos::deep_copy(yHost, yDevice);
        EXPECT_NEAR(yTest, yRef, 1e-12) << x << " " << yTest << " " << yRef;
        EXPECT_NEAR(yHost(0, 0), yRef, 1e-12) << x << " " << yTest << " " << yRef;
        for (size_t n = 1; n < 5; ++n) {
            double yP     = enge.function(x + dx, n - 1);
            double yM     = enge.function(x - dx, n - 1);
            double dyTest = enge.function(x, n);
            double dyRef  = (yP - yM) / 2 / dx;
            EXPECT_NEAR(dyTest, dyRef, 1e-3) << x << " " << n << " * " << dyTest << " " << dyRef;
            EXPECT_NEAR(yHost(0, n), dyRef, 1e-3)
                    << x << " " << n << " * " << yHost(0, n) << " " << dyRef;
            ;
        }
    }
    EXPECT_NEAR(enge.function(0.0, 0), 1.0, 1e-3);
    EXPECT_NEAR(enge.function(5.0, 0), 0.5, 1e-3);
    EXPECT_NEAR(enge.function(-5.0, 0), 0.5, 1e-3);
}

TEST_F(TestEnge, PrintTest) {
    std::vector<double> a    = {1, 2};
    endfieldmodel::Enge enge = endfieldmodel::Enge(a, 1.0, 0.5);
    enge.initialise();
    enge.print(std::cout);
    std::cout << std::endl;
}

TEST_F(TestEnge, FunctionGpuTest) {
    endfieldmodel::Enge enge(std::vector<double>({0, 1, 0}), 5, 1);
    enge.initialise();
    Kokkos::View<double*> xgpu("xgpu", 3);
    auto xhost = Kokkos::create_mirror_view(xgpu);
    xhost(0)   = 5.0;
    xhost(1)   = -5.0;
    xhost(2)   = 7.0;
    Kokkos::deep_copy(xgpu, xhost);  // target, source
    Kokkos::View<double**> derivativesgpu("d", 3, 5);

    enge.function(xgpu, 4, derivativesgpu);

    auto derivativeshost = Kokkos::create_mirror_view(derivativesgpu);
    Kokkos::deep_copy(derivativeshost, derivativesgpu);  // target, source
    for (auto i : std::vector<int>({0, 1, 2, 3, 4})) {
        EXPECT_NEAR(derivativeshost(0, i), enge.function(5, i), 1e-12);
        EXPECT_NEAR(derivativeshost(1, i), enge.function(-5, i), 1e-12);
        EXPECT_NEAR(derivativeshost(2, i), enge.function(7, i), 1e-12);
        std::cout << i << " " << xhost(0) << " " << xhost(1) << " " << xhost(2)
                  << " gpu-5: " << derivativeshost(0, i) << " cpu-5: " << enge.function(-5, i)
                  << " gpu+5 " << derivativeshost(1, i) << " cpu+5: " << enge.function(5, i)
                  << " gpu+7 " << derivativeshost(2, i) << " cpu+7: " << enge.function(7, i)
                  << std::endl;
    }
}

TEST_F(TestEnge, LinearPolynomialMatchesGeneralPolynomialThroughSixDerivatives) {
    endfieldmodel::Enge linear({0.2, 3.91}, .157, .14);
    endfieldmodel::Enge general({0.2, 3.91, 0.0}, .157, .14);
    linear.setMaximumDerivative(6);
    general.setMaximumDerivative(6);
    linear.initialise();
    general.initialise();
    for (double x : {-1.0, -.25, 0.0, .17, 1.0}) {
        for (int n = 0; n <= 6; ++n) {
            SCOPED_TRACE(::testing::Message() << "x=" << x << " derivative=" << n);
            EXPECT_DOUBLE_EQ(linear.function(x, n), general.function(x, n));
        }
    }
}

TEST_F(TestEnge, HighestSupportedOrderFitsWorkspaceAndIsReleased) {
    using namespace endfieldmodel;
    Enge model({.2}, 1., .5);
    model.initialise();
    constexpr int order = EngeConfig::max_derivative;
    Kokkos::View<double*> points("points", 2), result("guarded result", 2);
    Kokkos::deep_copy(points, .3);
    Kokkos::View<double**> values("derivatives", 2, order + 1);
    model.function(points, order, values);
    const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
    for (int i = 0; i < 2; ++i) {
        EXPECT_NEAR(host(i, 0), 2. / (1. + std::exp(.2)) - 1., 1e-14);
        for (int n = 1; n <= order; ++n)
            EXPECT_DOUBLE_EQ(host(i, n), 0.);
    }
    evaluateHighestOrderWithGuard(model.getConfig(), result);
    const auto checked = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), result);
    EXPECT_DOUBLE_EQ(checked(0), 0.);
    EXPECT_DOUBLE_EQ(checked(1), 12345.);
    EXPECT_DOUBLE_EQ(model.function(.3, order), 0.);
    std::array<double, order + 1> output;
    model.function(.3, order, output.data());
    // Separate host/device paths can differ by a few ULP after subtracting one.
    // Use the same absolute tolerance as the independent constant-value oracle.
    EXPECT_NEAR(output[0], host(0, 0), 1e-14);
    for (int n = 1; n <= order; ++n)
        EXPECT_DOUBLE_EQ(output[n], host(0, n));
    auto config = model.getConfig();
    Enge::releaseEngeConfig(config);
    for (int n = 0; n <= order; ++n) {
        EXPECT_EQ(config.gIndices[n].size(), 0u);
        EXPECT_EQ(config.hIndices[n].size(), 0u);
    }
}

TEST_F(TestEnge, InvalidOrdersAndOutputShapesAreRejectedBeforeEvaluation) {
    using namespace endfieldmodel;
    Enge model({0., 1.}, 1., .5);
    model.initialise();
    Kokkos::View<double*> points("points", 2);
    Kokkos::View<double**> values("unchanged output", 2, 1);
    Kokkos::deep_copy(values, 99.);
    std::array<double, EngeConfig::max_derivative + 1> output;
    output.fill(99.);
    for (int order : {-1, EngeConfig::max_derivative + 1, INT_MAX}) {
        EXPECT_THROW(model.function(.3, order), OpalException);
        EXPECT_THROW(model.function(.3, order, output.data()), OpalException);
        EXPECT_THROW(model.function(points, order, values), OpalException);
    }
    EXPECT_THROW(model.function(points, 1, values), OpalException);
    Kokkos::View<double**> wrongRows("wrong rows", 1, 1);
    EXPECT_THROW(model.function(points, 0, wrongRows), OpalException);
    EXPECT_THROW(model.setMaximumDerivative(EngeConfig::max_derivative + 1), OpalException);
    for (double value : output)
        EXPECT_DOUBLE_EQ(value, 99.);
    const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
    EXPECT_DOUBLE_EQ(host(0, 0), 99.);
    EXPECT_DOUBLE_EQ(host(1, 0), 99.);
}

TEST_F(TestEnge, CompactedTablesHaveExactCoefficientsThroughOrderTwelve) {
    using namespace endfieldmodel;
    using enge_test::factorial;
    Enge model({.2, .8, -.1}, .3, .7);
    model.initialise();
    const auto config = model.getConfig();
    // Number of integer partitions of n: one monomial for each partition.
    constexpr std::array<size_t, 13> partitions{1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77};
    std::int64_t largestCoefficient = 0;
    for (int n = 0; n <= enge_test::maxOrder; ++n) {
        SCOPED_TRACE(n);
        const auto table = Enge::getQIndex(n);
        ASSERT_EQ(table.size(), partitions[n]);
        const auto deviceTable =
                Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), config.gIndices[n]);
        ASSERT_EQ(deviceTable.extent(0), table.size());
        std::set<std::vector<std::int64_t>> uniquePowers;
        for (size_t row = 0; row < table.size(); ++row) {
            const auto& term = table[row];
            EXPECT_TRUE(uniquePowers.insert({term.begin() + 1, term.end()}).second);
            int count = 0, weight = 0;
            std::int64_t denominator = 1;
            for (size_t j = 2; j < term.size(); ++j) {
                ASSERT_GE(term[j], 0);
                count += term[j];
                weight += (j - 1) * term[j];
                denominator *= factorial(term[j]);
                for (int k = 0; k < term[j]; ++k)
                    denominator *= factorial(j - 1);
            }
            EXPECT_EQ(weight, n);
            EXPECT_EQ(term[1], -count - 1);
            // Faa di Bruno for (1/g)^(n): (-1)^k n! k! / prod(m_j! (j!)^m_j).
            const auto coefficient =
                    (count % 2 ? -1 : 1) * factorial(n) * factorial(count) / denominator;
            EXPECT_EQ(term[0], coefficient);
            largestCoefficient = std::max(largestCoefficient, std::abs(term[0]));
            for (size_t j = 0; j < term.size(); ++j)
                EXPECT_EQ(deviceTable(row, j), term[j]);
        }
        // h tables were already compact; check their coefficients independently too.
        if (n == 0) continue;
        const auto h = Enge::getHIndex(n);
        ASSERT_EQ(h.size(), partitions[n]);
        for (const auto& term : h) {
            int weight               = 0;
            std::int64_t denominator = 1;
            for (size_t j = 1; j < term.size(); ++j) {
                weight += j * term[j];
                denominator *= factorial(term[j]);
                for (int k = 0; k < term[j]; ++k)
                    denominator *= factorial(j);
            }
            EXPECT_EQ(weight, n);
            EXPECT_EQ(term[0], factorial(n) / denominator);
        }
    }
    EXPECT_EQ(largestCoefficient, INT64_C(5388768000));
    EXPECT_GT(largestCoefficient, std::numeric_limits<int>::max());
}

TEST_F(TestEnge, CompactedDerivativesMatchExpandedAndProductRuleReferences) {
    using namespace endfieldmodel;
    const auto expanded = enge_test::expandedTables();
    ASSERT_EQ(expanded.back().size(), 140152u);
    constexpr std::array<double, 5> x{-1., -.2, 0., .4, 1.1};
    for (const auto& coefficients : std::vector<std::vector<double>>{
                 {.2, .8}, {.2, .8, -.1}, {-.3, .5, -.1, .08, -.02, .01, .003}}) {
        Enge model(coefficients, .3, .7);
        model.initialise();
        Kokkos::View<double*> points("equivalence points", x.size());
        const auto input = Kokkos::create_mirror_view(points);
        for (size_t i = 0; i < x.size(); ++i)
            input(i) = x[i];
        Kokkos::deep_copy(points, input);
        Kokkos::View<double**> values("compacted derivatives", x.size(), enge_test::maxOrder + 1);
        model.function(points, enge_test::maxOrder, values);
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);
        for (size_t i = 0; i < x.size(); ++i) {
            const auto left  = enge_test::referenceEnd(coefficients, x[i] - .3L, .7);
            const auto right = enge_test::referenceEnd(coefficients, -x[i] - .3L, .7);
            for (int n = 0; n <= enge_test::maxOrder; ++n) {
                SCOPED_TRACE(
                        ::testing::Message()
                        << "degree=" << coefficients.size() - 1 << " x=" << x[i] << " n=" << n);
                const auto [legacyLeft, scaleLeft] = enge_test::evaluateTable(expanded[n], left.g);
                const auto [legacyRight, scaleRight] =
                        enge_test::evaluateTable(expanded[n], right.g);
                const int sign = n % 2 ? -1 : 1;
                const long double expected =
                        left.reciprocal[n] + sign * right.reciprocal[n] - (n == 0 ? 1 : 0);
                const long double legacy = legacyLeft + sign * legacyRight - (n == 0 ? 1 : 0);
                // Cancellation can amplify relative error near zeros. Bound absolute
                // roundoff by 128 epsilon times the L1 monomial sum, including unity
                // for the final subtraction. This is a new equivalence-test bound.
                const long double tolerance =
                        128 * std::numeric_limits<double>::epsilon() * (1 + scaleLeft + scaleRight);
                EXPECT_NEAR(legacy, expected, tolerance);
                EXPECT_NEAR(host(i, n), legacy, tolerance);
                EXPECT_NEAR(host(i, n), expected, tolerance);
            }
        }
    }
}

TEST_F(TestEnge, HostPointsMatchBatchedValuesAfterModelUpdates) {
    using namespace endfieldmodel;
    Enge model({.1, 1.2, -.05}, .2, .15);
    model.initialise();
    const auto check = [](const Enge& current) {
        constexpr int n = 6;
        Kokkos::View<double*> points("host comparison points", 3);
        auto hostPoints = Kokkos::create_mirror_view(points);
        hostPoints(0)   = -.31;
        hostPoints(1)   = .02;
        hostPoints(2)   = .27;
        Kokkos::deep_copy(points, hostPoints);
        Kokkos::View<double**> derivatives("batched comparison", 3, n + 1);
        current.function(points, n, derivatives);
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), derivatives);
        for (int point = 0; point < 3; ++point) {
            std::array<double, n + 1> values;
            current.function(hostPoints(point), n, values.data());
            for (int order = 0; order <= n; ++order) {
                const double tolerance = 1e-12 * std::max(1., std::abs(host(point, order)));
                EXPECT_NEAR(values[order], host(point, order), tolerance);
                EXPECT_NEAR(
                        current.function(hostPoints(point), order), host(point, order), tolerance);
            }
        }
    };
    check(model);
    model.setCoefficients({-.1, .8});
    check(model);
    model.setLambda(.21);
    model.setX0(.13);
    check(model);
    model.rescale(1.3);
    check(model);
    model.setMaximumDerivative(6);
    check(model);
    std::unique_ptr<Enge> copy(model.clone());
    copy->initialise();
    check(*copy);
    model.finalise();
    model.setCoefficients({.2, 1.});
    model.initialise();
    check(model);
}

TEST_F(TestEnge, ConfigurationReplacementRefreshesWarmHostMirrors) {
    using namespace endfieldmodel;
    Enge model({.1, 1.2, -.05}, .2, .15);
    model.initialise();
    std::array<double, EngeConfig::max_derivative + 1> mirrored{};
    EngeHostMirrorTestAccess::evaluate(model, .02, EngeConfig::max_derivative, mirrored.data());
    const double previous = mirrored[0];
    // Allocate a distinct coefficient/index configuration. Updating shared
    // views in place could hide stale mirrors on host-accessible backends.
    Enge replacement({-.4, .7, .03, -.01}, .35, .22);
    replacement.initialise();
    auto config = replacement.getConfig();
    model.setConfig(config);
    EXPECT_GT(std::abs(replacement.function(.02, 0) - previous), .01);
    Kokkos::View<double*> points("replacement points", 3);
    auto hostPoints = Kokkos::create_mirror_view(points);
    hostPoints(0)   = -.31;
    hostPoints(1)   = .02;
    hostPoints(2)   = .27;
    Kokkos::deep_copy(points, hostPoints);
    for (int order : {0, 3, 6, EngeConfig::max_derivative}) {
        Kokkos::View<double**> batched("replacement derivatives", 3, order + 1);
        model.function(points, order, batched);
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), batched);
        for (int point = 0; point < 3; ++point) {
            std::array<double, EngeConfig::max_derivative + 1> expected{}, scalar{};
            replacement.function(hostPoints(point), order, expected.data());
            model.function(hostPoints(point), order, scalar.data());
            EngeHostMirrorTestAccess::evaluate(model, hostPoints(point), order, mirrored.data());
            for (int n = 0; n <= order; ++n) {
                const double tolerance = 1e-12 * std::max(1., std::abs(expected[n]));
                EXPECT_NEAR(mirrored[n], expected[n], tolerance) << point << ',' << n;
                EXPECT_NEAR(scalar[n], expected[n], tolerance) << point << ',' << n;
                EXPECT_NEAR(host(point, n), expected[n], tolerance) << point << ',' << n;
            }
        }
    }
}
