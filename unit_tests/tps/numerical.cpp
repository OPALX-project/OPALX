#include <Kokkos_Core.hpp>
#include "tps/Algebra/Tps.hpp"
#include "tps/Algebra/TpsMath.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace tps;
using TPS = Tps<double>;
void near(double actual, double expected, double tol = 2e-13) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tol)
        throw std::runtime_error("actual=" + std::to_string(actual) +
                                 " expected=" + std::to_string(expected));
}
int main(int argc, char** argv) {
    Kokkos::ScopeGuard kokkos;

    try {
        // Inspect the unsafe default without executing an overflowing allocation.
        if (TPS::getGlobalTruncOrder() != TPS::EXACT)
            throw std::runtime_error("unexpected default");
        TPS raw(3, 6);
        if (raw.getMaxOrder() != 3 || raw.getTruncOrder() != TPS::EXACT)
            throw std::runtime_error("default characterization");
        if (argc > 1) {
            const std::string mode = argv[1];
            if (mode != "constant_multiplication" && mode != "linear_substitution")
                throw std::runtime_error("unknown regression");
            for (int variables : {1, 2, 6})
                for (int order : {1, 2, 3, 5}) {
                    TPS::setGlobalTruncOrder(order);
                    auto x = TPS::makeVariable(variables, 0);
                    auto y = TPS::makeVariable(variables, variables - 1);
                    Vector<double> p(variables, 0.0);
                    p[0] = 0.02;
                    if (variables > 1)
                        p[variables - 1] = -0.03;
                    if (mode == "constant_multiplication") {
                        auto polynomial = 4.0 + x + 2.0 * y;
                        if (order >= 2)
                            polynomial += x * x + 3.0 * x * y;
                        for (double c : {0.0, 1.0, 2.0, -3.0}) {
                            auto left = TPS(c) * polynomial;
                            auto right = polynomial * TPS(c);
                            auto scalar = c * polynomial;
                            for (int i = 0; i < polynomial.getSize(); ++i) {
                                near(left.getCoefficient(i), c * polynomial.getCoefficient(i));
                                near(right.getCoefficient(i), c * polynomial.getCoefficient(i));
                                near(scalar.getCoefficient(i), c * polynomial.getCoefficient(i));
                            }
                            near(left.evaluate(p), c * polynomial.evaluate(p));
                        }
                    } else {
                        Matrix<double> m(variables, variables, 0.0);
                        for (int i = 0; i < variables; ++i)
                            m(i, i) = 1;
                        near(x.substitute(m).evaluate(p), p[0]);
                        m(0, variables - 1) += 2;
                        const double transformed = p[0] + 2 * p[variables - 1];
                        auto linear = x.substitute(m);
                        near(linear.evaluate(p), transformed);
                        for (int j = 0; j < variables; ++j)
                            near(linear.getCoefficient(j + 1), m(0, j));
                        if (order >= 2) {
                            auto actual = (x * x).substitute(m);
                            auto expected = linear * linear;
                            for (int i = 0; i < expected.getSize(); ++i)
                                near(actual.getCoefficient(i), expected.getCoefficient(i));
                            near(actual.evaluate(p), transformed * transformed);
                        }
                    }
                }
            return 0;
        }
        for (int order : {1, 2, 3, 5}) {
            TPS::setGlobalTruncOrder(order);
            auto x = TPS::makeVariable(6, 0), y = TPS::makeVariable(6, 1);
            Vector<double> p(6, 0.0);
            p[0] = 0.02;
            p[1] = -0.03;
            auto f = (x + 2.0) * (y + 3.0);
            near(f.evaluate(p), 6 + 3 * p[0] + 2 * p[1] + (order >= 2 ? p[0] * p[1] : 0));
            if (order >= 2) {
                near(f.derivative(0).evaluate(p), 3 + p[1]);
                near(f.derivative(1).evaluate(p), 2 + p[0]);
            }
            auto e = exp(x);
            double expected = 1, term = 1;
            for (int k = 1; k <= order; ++k) {
                term *= p[0] / k;
                expected += term;
            }
            near(e.evaluate(p), expected);
            auto inv = (x + 1.0).inverse();
            expected = 0;
            term = 1;
            for (int k = 0; k <= order; ++k) {
                expected += term;
                term *= -p[0];
            }
            near(inv.evaluate(p), expected);
            near(((x + 1.0) * inv).evaluate(p), 1);
            auto identity = log(exp(x));
            near(identity.evaluate(p), p[0]);
            auto trig = sin(x) * sin(x) + cos(x) * cos(x);
            near(trig.evaluate(p), 1);
            auto copy = x;
            copy.setCoefficient(0, 7);
            near(x.getCoefficient(0), 0);
        }
        TPS::setGlobalTruncOrder(TPS::EXACT);
        std::cout << "All numerical checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
