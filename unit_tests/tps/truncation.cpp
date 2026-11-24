#include <Kokkos_Core.hpp>
#include "tps/Algebra/Tps.hpp"
#include "tps/Algebra/TpsMath.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace tps;
using TPS = Tps<double>;
void near(double a, double b) {
    if (!std::isfinite(a) || std::abs(a - b) > 2e-13)
        throw std::runtime_error("truncation numerical mismatch");
}
int main() {
    Kokkos::ScopeGuard kokkos;

    try {
        if (TPS::getGlobalTruncOrder() != TPS::EXACT)
            return 1;
        for (int order : {0, 1, 2, 3, 5}) {
            TPS x(order, 2);
            if (order)
                x.setCoefficient(1, 1);
            auto e = exp(x), inv = (x + 1.0).inverse(), logarithm = log(x + 1.0);
            for (int k = 0; k <= order; ++k) {
                TpsMonomial m(2);
                m[0] = k;
                double factorial = 1;
                for (int j = 1; j <= k; ++j)
                    factorial *= j;
                near(k ? e[m] : e[0], 1 / factorial);
                near(k ? inv[m] : inv[0], k % 2 ? -1 : 1);
                near(k ? logarithm[m] : logarithm[0], k ? (k % 2 ? 1.0 : -1.0) / k : 0);
            }
            // Functions with a first-derivative seed must also support order zero.
            near(sin(x)[0], 0);
            near(cos(x)[0], 1);
            near(sinh(x)[0], 0);
            near(cosh(x)[0], 1);
            near(sqrt(1.0 + x)[0], 1);
            for (int k = 1; k <= order; ++k) {
                TpsMonomial m(2);
                m[0] = k;
                double factorial = 1;
                for (int j = 1; j <= k; ++j)
                    factorial *= j;
                near(sin(x)[m], k % 2 ? ((k % 4 == 1) ? 1.0 : -1.0) / factorial : 0);
                near(cos(x)[m], k % 2 ? 0 : ((k % 4 == 0) ? 1.0 : -1.0) / factorial);
                near(sinh(x)[m], k % 2 ? 1.0 / factorial : 0);
                near(cosh(x)[m], k % 2 ? 0 : 1.0 / factorial);
            }
            if (order) {
                near(sin(x)[1], 1);
                near(sinh(x)[1], 1);
                near(sqrt(1.0 + x)[1], 0.5);
            }
        }
        // Scalar constants must not allocate an infinite Taylor series.
        TPS c(2.0);
        near(exp(c)[0], std::exp(2.0));
        near(log(c)[0], std::log(2.0));
        near(sin(c)[0], std::sin(2.0));
        near(cos(c)[0], std::cos(2.0));
        near(sinh(c)[0], std::sinh(2.0));
        near(cosh(c)[0], std::cosh(2.0));
        near(sqrt(c)[0], std::sqrt(2.0));
        near(c.inverse()[0], 0.5);
        auto x = TPS::makeVariable(2, 0);
        // EXACT still permits polynomial degree growth.
        if ((x * x).getMaxOrder() != 2)
            return 2;
        near(x.integral(0).getCoefficient(3), 0.5);
        near(x.multiplyVariable(0).getCoefficient(3), 1);
        // Explicit inverse order and global order can exceed the stored degree.
        TpsMonomial fifth(2);
        fifth[0] = 5;
        near((1.0 + x).inverse(5)[fifth], -1);
        TPS::setGlobalTruncOrder(5);
        auto e = exp(x);
        near(e[fifth], 1.0 / 120);
        TPS::setGlobalTruncOrder(0);
        near(sin(x)[0], 0);
        near(log(1.0 + x)[0], 0);
        near(cosh(x)[0], 1);
        TPS::setGlobalTruncOrder(2);
        if (exp(x).getMaxOrder() != 2)
            return 3;
        bool rejected = false;
        try {
            TPS::setGlobalTruncOrder(-1);
        } catch (const CLRangeError&) {
            rejected = true;
        }
        if (!rejected || TPS::getGlobalTruncOrder() != 2)
            return 4;
        rejected = false;
        try {
            (1.0 + x).inverse(-1);
        } catch (const CLRangeError&) {
            rejected = true;
        }
        if (!rejected)
            return 5;
        TPS::setGlobalTruncOrder(TPS::EXACT);
        if (exp(x).getMaxOrder() != 1)
            return 6;
        std::cout << "Truncation regressions passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 7;
    }
}
