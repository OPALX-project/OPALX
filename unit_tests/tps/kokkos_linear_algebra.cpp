#include <tps/tps.hpp>
#include <tps/LinearSolve.hpp>
#include <complex>
#include <cmath>
#include <iostream>
#include <stdexcept>
void near(double a, double b) {
    if (!std::isfinite(a) || std::abs(a - b) > 2e-13)
        throw std::runtime_error("Kokkos linear algebra mismatch");
}
int main() {
    Kokkos::ScopeGuard runtime;
    try {
        tps::Matrix<double> a(2, 2, 0.0);
        a(0, 0) = 0;
        a(0, 1) = 2;
        a(1, 0) = 1;
        a(1, 1) = 3;
        tps::LUFactorization<double> factor(a);
        tps::Vector<double> b(2, 0.0);
        b[0] = 4;
        b[1] = 7;
        factor.solve(b);
        near(b[0], 1);
        near(b[1], 2);
        near(a(0, 0), 0); // pivoting required; input unchanged
        auto identity = a * factor.inverse();
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                near(identity(i, j), i == j ? 1 : 0);
        auto rhs = tps::Matrix<double>::identity(2);
        factor.solve(rhs);
        auto check = a * rhs;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                near(check(i, j), i == j ? 1 : 0);
        auto twice = 2.0 * a;
        near(twice(0, 1), 4);
        near((twice / 2.0)(0, 1), 2);
        near((a + a)(1, 1), 6);
        near((a - a)(0, 1), 0);
        near((-a)(1, 1), -3);
        near(a.transpose()(1, 0), 2);
        tps::Vector<double> x(2, 1.0);
        auto y = a * x;
        near(y[0], 2);
        near(y[1], 4);
        auto z = x * a;
        near(z[0], 1);
        near(z[1], 5);
        near(x * x, 2);
        near(tps::euclidean_norm(x), std::sqrt(2.0));
        tps::Array1D<double> d(2, 2);
        near(tps::scaled_norm(d, x), std::sqrt(8.0));
        bool singular = false;
        try {
            tps::Matrix<double> zero(2, 2, 0.0);
            tps::LUFactorization<double> bad(zero);
        } catch (const tps::SingularMatrixError&) {
            singular = true;
        }
        if (!singular)
            throw std::runtime_error("singular matrix was accepted");
        bool mismatch = false;
        try {
            tps::Vector<double> wrong(3);
            (void)(x + wrong);
        } catch (const tps::SizeError&) {
            mismatch = true;
        }
        if (!mismatch)
            throw std::runtime_error("bad vector shape accepted");
        using C = std::complex<double>;
        tps::Vector<C> cx(2, C(0, 1));
        auto dot = cx * cx;
        near(dot.real(), -2);
        near(dot.imag(), 0);
        near(tps::euclidean_norm(cx), std::sqrt(2.0));
        tps::Matrix<C> ca(2, 2, C(0));
        ca(0, 0) = C(1, 1);
        ca(1, 1) = C(2, -1);
        auto cb = ca * ca;
        near(cb(0, 0).imag(), 2);
        near(cb(1, 1).real(), 3);
        // Fixed-size device values, dynamic View storage, explicit host/device transfer.
        Kokkos::View<double* [6]> particles("particles", 4);
        Kokkos::Array<double, 6> shift{1, 2, 3, 4, 5, 6};
        Kokkos::parallel_for(
            "fixed coordinate data", 4, KOKKOS_LAMBDA(int i) {
                for (int j = 0; j < 6; ++j)
                    particles(i, j) = i + shift[j];
            });
        auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), particles);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 6; ++j)
                near(host(i, j), i + j + 1);
        std::cout << "Kokkos linear algebra and array checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
