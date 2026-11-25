#include "Utilities/LogicalError.h"
#include <Kokkos_Core.hpp>
#include <tps/tps.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
void near(double a, double b) {
    if (!std::isfinite(a) || std::abs(a - b) > 2e-13)
        throw std::runtime_error("map defect regression mismatch: " + std::to_string(a) +
                                 " expected " + std::to_string(b));
}
int main() {
    Kokkos::ScopeGuard kokkos;
    try {
        using S = tps::Tps<double>;
        auto q = S::makeVariable(2, 0), p = S::makeVariable(2, 1);
        tps::VpsMap<double> f(2, 2);
        f[0] = q;
        f[1] = p;
        tps::Vector<double> x(2, 0.0);
        x[0] = 0.25;
        x[1] = 0.5;
        for (int rows : {1, 2, 3}) {
            tps::Matrix<double> m(rows, 2, 0.0);
            for (int i = 0; i < rows; ++i) {
                m(i, 0) = i + 1;
                m(i, 1) = -i;
            }
            auto result = (m * f).constantTerm(x);
            for (int i = 0; i < rows; ++i)
                near(result[i], (i + 1) * x[0] - i * x[1]);
        }
        tps::Matrix<double> bad(2, 3, 0.0);
        bool rejected = false;
        try {
            (void)(bad * f);
        } catch (const ::LogicalError&) {
            rejected = true;
        }
        if (!rejected)
            return 2;
        f[0] = S(2.0);
        f[1] = q * q + 3.0 * p;
        auto origin = f.linearTerms();
        auto atPoint = f.linearTerms(x);
        near(origin(0, 0), 0);
        near(origin(0, 1), 0);
        near(origin(1, 0), 0);
        near(origin(1, 1), 3);
        near(atPoint(0, 0), 0);
        near(atPoint(0, 1), 0);
        near(atPoint(1, 0), 0.5);
        near(atPoint(1, 1), 3);
        tps::Vector<double> wrong(3, 0.0);
        rejected = false;
        try {
            (void)f.linearTerms(wrong);
        } catch (const ::LogicalError&) {
            rejected = true;
        }
        if (!rejected)
            return 3;
        f[0] = S(0, 2);
        near(f.linearTerms()(0, 0), 0);
        f[0] = S(2.0);
        f[1] = S(4.0);
        auto constants = f.linearTerms(x);
        if (constants.nrows() != 2 || constants.ncols() != 2)
            return 4;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                near(constants(i, j), 0);
        std::cout << "Map defect regressions passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
