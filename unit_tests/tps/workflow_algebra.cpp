#include <Kokkos_Core.hpp>
#include <tps/tps.hpp>
#include <tps/LinearSolve.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
void near(double a, double b) {
    if (std::abs(a - b) > 2e-13)
        throw std::runtime_error("workflow algebra mismatch");
}
int main() {
    Kokkos::ScopeGuard kokkos;
    try {
        using S = tps::Tps<double>;
        S::setGlobalTruncOrder(4);
        auto q = S::makeVariable(2, 0), p = S::makeVariable(2, 1);
        tps::VpsMap<double> drift(2, 2), kick(2, 2), scale(2, 2);
        drift[0] = q + 2.0 * p;
        drift[1] = p;
        kick[0] = q;
        kick[1] = p - 3.0 * q * q;
        scale[0] = 2.0 * q;
        scale[1] = 0.5 * p;
        auto combined = scale.substitute(kick.substitute(drift, 4), 4);
        auto alternative = scale.substitute(kick, 4).substitute(drift, 4);
        tps::Vector<double> x(2, 0.0);
        x[0] = 0.1;
        x[1] = 0.2;
        auto sequential = scale.constantTerm(kick.constantTerm(drift.constantTerm(x)));
        auto evaluated = combined.constantTerm(x), other = alternative.constantTerm(x);
        near(evaluated[0], 1.0);
        near(evaluated[1], -0.275);
        for (int i = 0; i < 2; ++i) {
            near(evaluated[i], sequential[i]);
            near(evaluated[i], other[i]);
        }
        // Preserve a standalone linear solve while decoupling it from the TPS headers.
        tps::Matrix<double> matrix(2, 2, 0.0);
        matrix(0, 0) = 2;
        matrix(0, 1) = 1;
        matrix(1, 0) = 1;
        matrix(1, 1) = 3;
        tps::LUFactorization<double> lu(matrix);
        tps::Vector<double> rhs(2, 0.0);
        rhs[0] = 4;
        rhs[1] = 7;
        lu.solve(rhs);
        near(rhs[0], 1);
        near(rhs[1], 2);
        auto inverse = lu.inverse();
        auto identity = matrix * inverse;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                near(identity(i, j), i == j ? 1 : 0);
        // Complex TPS remains supported without forced pre-instantiation files.
        using C = std::complex<double>;
        auto z = tps::Tps<C>::makeVariable(1, 0);
        auto complexPoly = (z + C(0, 1)) * (z + C(0, 1));
        if (std::abs(complexPoly[0] + C(1)) > 2e-13)
            return 2;
        std::cout << "Workflow algebra checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
