#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <array>
#include <cmath>
#include <tps/tps.hpp>

namespace {
    using S              = tps::Tps<double>;
    constexpr int degree = 3;
    using Coefficients   = std::array<std::array<double, 4>, 4>;

    class Algebra : public ::testing::Test {
        int saved = S::getGlobalTruncOrder();

    protected:
        void SetUp() override { S::setGlobalTruncOrder(degree); }
        void TearDown() override { S::setGlobalTruncOrder(saved); }
    };

    // Mixed absolute/relative tolerance allows double rounding in cubic reciprocal
    // and transcendental coefficients; no existing numerical tolerance is changed.
    void check(const S& actual, const Coefficients& expected) {
        for (int i = 0; i <= degree; ++i)
            for (int j = 0; i + j <= degree; ++j) {
                SCOPED_TRACE(::testing::Message() << "x^" << i << " y^" << j);
                tps::TpsMonomial m(2);
                m[0] = i;
                m[1] = j;
                const double value =
                        i + j == 0 ? actual[0] : (i + j <= actual.getMaxOrder() ? actual[m] : 0.);
                EXPECT_NEAR(value, expected[i][j], 5e-13 * (1. + std::abs(expected[i][j])));
            }
    }
    S series(const Coefficients& coefficients) {
        S result(degree, 2);
        for (int i = 0; i <= degree; ++i)
            for (int j = 0; i + j <= degree; ++j) {
                tps::TpsMonomial m(2);
                m[0]      = i;
                m[1]      = j;
                result[m] = coefficients[i][j];
            }
        return result;
    }

    TEST_F(Algebra, ArithmeticAndAliasingMatchIndependentConvolution) {
        Coefficients a{}, b{}, sum{}, difference{}, product{};
        for (int i = 0; i <= degree; ++i)
            for (int j = 0; i + j <= degree; ++j) {
                a[i][j]          = (1. + i - 2. * j) / 8.;
                b[i][j]          = (2. - 3. * i + j) / 4.;
                sum[i][j]        = a[i][j] + b[i][j];
                difference[i][j] = a[i][j] - b[i][j];
                for (int k = 0; k <= i; ++k)
                    for (int l = 0; l <= j; ++l)
                        product[i][j] += a[k][l] * b[i - k][j - l];
            }
        const auto p = series(a), q = series(b);
        check(p + q, sum);
        check(p - q, difference);
        check(p * q, product);
        auto r = p;
        r += q;
        check(r, sum);
        r = p;
        r -= q;
        check(r, difference);
        r = p;
        r *= q;
        check(r, product);
        check((p * q) / q, a);
        r = p * q;
        r /= q;
        check(r, a);
        r = p;
        r -= r;
        check(r, Coefficients{});
        r = q;
        r /= r;
        Coefficients one{};
        one[0][0] = 1.;
        check(r, one);
        check(p, a);
        check(q, b);
        for (double c : {-2., .5, 3.}) {
            check((p + c) - c, a);
            check(c - (c - p), a);
            check((c * p) / c, a);
            check((c / q) * q, [&] {
                Coefficients e{};
                e[0][0] = c;
                return e;
            }());
            r = p;
            r += c;
            r -= c;
            r *= c;
            r /= c;
            check(r, a);
        }
    }

    TEST_F(Algebra, ShiftedElementaryFunctionsMatchAnalyticMixedCoefficients) {
        const auto u   = .7 + S::makeVariable(2, 0) + 2. * S::makeVariable(2, 1);
        const double a = .7, s = std::sin(a), c = std::cos(a), t = std::tan(a);
        const double h = std::tanh(a), sh = std::sinh(a), ch = std::cosh(a);
        const double e = std::exp(a), root = std::sqrt(a), sec = 1. / c, csc = 1. / s, cot = c / s;
        struct Case {
            const char* name;
            S value;
            std::array<double, 4> derivatives;
        };
        const std::array cases{
                Case{"exp", tps::exp(u), {e, e, e, e}},
                Case{"log", tps::log(u), {std::log(a), 1 / a, -1 / (a * a), 2 / (a * a * a)}},
                Case{"sqrt",
                     tps::sqrt(u),
                     {root, 1 / (2 * root), -1 / (4 * a * root), 3 / (8 * a * a * root)}},
                Case{"sin", tps::sin(u), {s, c, -s, -c}},
                Case{"cos", tps::cos(u), {c, -s, -c, s}},
                Case{"tan",
                     tps::tan(u),
                     {t, 1 + t * t, 2 * t * (1 + t * t), 2 * (1 + t * t) * (1 + 3 * t * t)}},
                Case{"sec",
                     tps::sec(u),
                     {sec, sec * t, sec * (1 + 2 * t * t), sec * t * (5 + 6 * t * t)}},
                Case{"csc",
                     tps::csc(u),
                     {csc, -csc * cot, csc * (1 + 2 * cot * cot),
                      -csc * cot * (5 + 6 * cot * cot)}},
                Case{"sinh", tps::sinh(u), {sh, ch, sh, ch}},
                Case{"cosh", tps::cosh(u), {ch, sh, ch, sh}},
                Case{"tanh",
                     tps::tanh(u),
                     {h, 1 - h * h, -2 * h * (1 - h * h), -2 * (1 - h * h) * (1 - 3 * h * h)}},
                Case{"positive power", tps::pow(u, 3), {a * a * a, 3 * a * a, 6 * a, 6}},
                Case{"negative power",
                     tps::pow(u, -2),
                     {1 / (a * a), -2 / (a * a * a), 6 / std::pow(a, 4), -24 / std::pow(a, 5)}}};
        constexpr std::array factorial{1., 1., 2., 6.};
        for (const auto& test : cases) {
            SCOPED_TRACE(test.name);
            Coefficients expected{};
            for (int i = 0; i <= degree; ++i)
                for (int j = 0; i + j <= degree; ++j)
                    expected[i][j] = test.derivatives[i + j] * std::pow(2., j)
                                     / (factorial[i] * factorial[j]);
            check(test.value, expected);
        }
        Coefficients one{};
        one[0][0] = 1.;
        check(tps::pow(u, 0), one);
    }

    TEST_F(Algebra, PartialCalculusMatchesExponentRules) {
        Coefficients coefficients{};
        for (int i = 0; i <= degree; ++i)
            for (int j = 0; i + j <= degree; ++j)
                coefficients[i][j] = (2. + i - j) / 8.;
        const auto p = series(coefficients);
        for (int variable : {0, 1}) {
            Coefficients derivative{}, integral{}, timesVariable{};
            for (int i = 0; i <= degree; ++i)
                for (int j = 0; i + j <= degree; ++j) {
                    const int exponent = variable == 0 ? i : j;
                    if (exponent > 0) {
                        const int k = i - (variable == 0), l = j - (variable == 1);
                        derivative[k][l]    = exponent * coefficients[i][j];
                        integral[i][j]      = coefficients[k][l] / exponent;
                        timesVariable[i][j] = coefficients[k][l];
                    }
                }
            check(p.derivative(variable), derivative);
            check(p.integral(variable), integral);
            check(p.multiplyVariable(variable), timesVariable);
        }
        check(p.derivative(0).derivative(1), [&] {
            Coefficients expected{};
            for (int i = 1; i <= degree; ++i)
                for (int j = 1; i + j <= degree; ++j)
                    expected[i - 1][j - 1] = i * j * coefficients[i][j];
            return expected;
        }());
    }
}  // namespace
int main(int argc, char** argv) {
    Kokkos::ScopeGuard guard(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
