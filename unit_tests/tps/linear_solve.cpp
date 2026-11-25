#include "Utilities/LogicalError.h"
#include "Utilities/SingularMatrixError.h"
#include <tps/LinearSolve.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool value) {
    if (!value)
        throw std::runtime_error("LU regression failed");
}
void near(double value, double expected) {
    require(std::isfinite(value) &&
            std::abs(value - expected) <= 3e-13 * std::max(1., std::abs(expected)));
}
template <class Exception, class F> void rejects(F action) {
    bool caught = false;
    try {
        action();
    } catch (const Exception&) {
        caught = true;
    }
    require(caught);
}
} // namespace
int main() {
    Kokkos::ScopeGuard runtime;
    try {
        tps::Matrix<double> a(3, 3, 0.);
        a(0, 1) = 2;
        a(0, 2) = 1;
        a(1, 0) = 1;
        a(1, 1) = 1;
        a(2, 0) = 2;
        a(2, 2) = 3;
        const auto saved = a;
        tps::LUFactorization<double> factor(a);
        near(factor.determinant(), -8);
        near(tps::determinant(a), -8);
        near(tps::infinityNorm(a), 5);
        tps::Matrix<double> expected(3, 2, 0.);
        for (int i = 0; i < 3; ++i) {
            expected(i, 0) = i + 1.;
            expected(i, 1) = 2. - i;
        }
        auto rhs = a.transpose() * expected;
        factor.solveTranspose(rhs);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 2; ++j)
                near(rhs(i, j), expected(i, j));
            for (int j = 0; j < 3; ++j)
                near(a(i, j), saved(i, j));
        }
        tps::Vector<double> x(3, 0.);
        x[0] = 2;
        x[1] = -3;
        x[2] = .5;
        auto b = a.transpose() * x;
        factor.solveTranspose(b);
        for (int i = 0; i < 3; ++i)
            near(b[i], x[i]);
        tps::Vector<double> wrongVector(2);
        tps::Matrix<double> wrongMatrix(2, 3);
        rejects<::LogicalError>([&] { factor.solveTranspose(wrongVector); });
        rejects<::LogicalError>([&] { factor.solveTranspose(wrongMatrix); });
        auto diagonal = tps::Matrix<double>::identity(3);
        diagonal(0, 0) = 1e-18;
        tps::LUFactorization<double> boundary(diagonal, 1e-18);
        diagonal(0, 0) = std::nextafter(1e-18, 0.);
        rejects<::SingularMatrixError>([&] { tps::LUFactorization<double> bad(diagonal, 1e-18); });
        diagonal(0, 0) = std::nextafter(1e-18, 1.);
        tps::LUFactorization<double> above(diagonal, 1e-18);
        diagonal(0, 0) = 0;
        near(tps::determinant(diagonal), 0);
        rejects<::SingularMatrixError>([&] { tps::LUFactorization<double> bad(diagonal); });
        auto permutation = tps::Matrix<double>::identity(6);
        for (int shift : {1, 2, 3, 5}) {
            permutation = tps::Matrix<double>(6, 6, 0.);
            for (int i = 0; i < 6; ++i)
                permutation(i, (i + shift) % 6) = 1;
            // A six-cycle shift has parity (-1)^(6-gcd(6,shift)).
            const double parity = shift == 2 ? 1 : -1;
            near(tps::determinant(permutation), parity);
        }
        rejects<std::invalid_argument>([&] { tps::LUFactorization<double> bad(a, -1); });
        rejects<std::invalid_argument>(
            [&] { tps::LUFactorization<double> bad(a, std::numeric_limits<double>::infinity()); });
        tps::Matrix<double> nonsquare(2, 3);
        rejects<::LogicalError>([&] { (void)tps::determinant(nonsquare); });
        diagonal(0, 0) = std::numeric_limits<double>::quiet_NaN();
        rejects<std::invalid_argument>([&] { (void)tps::determinant(diagonal); });
        require(std::isnan(tps::infinityNorm(diagonal)));
        auto empty = tps::Matrix<double>::identity(0);
        near(tps::determinant(empty), 1);
        near(tps::infinityNorm(empty), 0);
        std::cout << "Pivoted LU, transpose solve, determinant and norm regressions passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
