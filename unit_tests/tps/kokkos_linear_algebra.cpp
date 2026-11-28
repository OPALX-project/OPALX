#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
// The original complex mismatch must fail CUDA compilation, not just emit warnings.
#pragma nv_diag_error 20011, 20013, 20014, 20015
#endif

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <tps/LinearSolve.hpp>
#include <tps/tps.hpp>
#include "Utilities/LogicalError.h"
#include "Utilities/SingularMatrixError.h"
void near(double a, double b) {
    if (!std::isfinite(a) || std::abs(a - b) > 2e-13)
        throw std::runtime_error("Kokkos linear algebra mismatch");
}

static_assert(tps::LinearAlgebraScalar<Kokkos::complex<float>>);
static_assert(tps::LinearAlgebraScalar<Kokkos::complex<double>>);
static_assert(!tps::LinearAlgebraScalar<std::complex<float>>);
static_assert(!tps::LinearAlgebraScalar<std::complex<double>>);

template <class Real>
void complexAlgebra() {
    using C                = Kokkos::complex<Real>;
    const double tolerance = 64 * std::numeric_limits<Real>::epsilon();
    const auto close       = [=](auto actual, auto expected) {
        if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
            throw std::runtime_error("Complex Kokkos kernel mismatch");
    };
    const auto complexClose = [&](C actual, C expected) {
        close(actual.real(), expected.real());
        close(actual.imag(), expected.imag());
    };
    tps::Vector<C> x(2), y(2);
    x[0] = C(1, 2);
    x[1] = C(-3, 4);
    y[0] = C(2, -1);
    y[1] = C(5, 2);
    //! [complex kernel semantics]
    // Sum x_i*y_i = -19+17i, whereas the Hermitian product is -7-31i.
    complexClose(x * y, C(-19, 17));
    complexClose(x * x, C(-10, -20));
    // Norm = sqrt(sum |x_i|^2), independently of the bilinear dot product.
    close(tps::euclidean_norm(x), std::sqrt(30.0));
    tps::Array1D<C> weights(2);
    weights[0] = C(2, 1);
    weights[1] = C(-1, 2);
    close(tps::scaled_norm(weights, x), std::sqrt(150.0));
    //! [complex kernel semantics]
    auto scaled = C(2, -1) * x;
    complexClose(scaled[0], C(4, 3));
    complexClose(scaled[1], C(-2, 11));
    auto restored = scaled / C(2, -1);
    complexClose(restored[0], x[0]);
    complexClose(restored[1], x[1]);
    auto sum = x + y;
    complexClose(sum[0], C(3, 1));
    complexClose(sum[1], C(2, 6));
    sum -= y;
    complexClose(sum[0], x[0]);
    complexClose(x[0], C(1, 2));  // arithmetic on copies must not change the source
    tps::Vector<C> empty;
    complexClose(empty * empty, C(0));
    close(tps::euclidean_norm(empty), 0.0);

    tps::Matrix<C> a(2, 2, C(0)), b(2, 2, C(0));
    a(0, 1)            = C(1, 1);
    a(1, 0)            = C(2, -1);  // zero first pivot requires a row swap
    a(1, 1)            = C(3, 2);
    b(0, 0)            = C(1, -1);
    b(0, 1)            = C(2);
    b(1, 0)            = C(0, 1);
    b(1, 1)            = C(-1, 2);
    const auto product = a * b;
    complexClose(product(0, 0), C(-1, 1));
    complexClose(product(0, 1), C(-3, 1));
    complexClose(product(1, 0), C(-1));
    complexClose(product(1, 1), C(-3, 2));
    complexClose(a.transpose()(1, 0), C(1, 1));  // plain transpose, no conjugation
    complexClose((C(2, 1) * a)(0, 1), C(1, 3));
    complexClose((a + b)(1, 1), C(2, 4));
    complexClose((a - b)(1, 1), C(4));
    const auto column = a * x, row = x * a;
    complexClose(column[0], C(-7, 1));
    complexClose(column[1], C(-13, 9));
    complexClose(row[0], C(-2, 11));
    complexClose(row[1], C(-18, 9));
    tps::LUFactorization<C> factor(a);
    complexClose(factor.determinant(), C(-3, -1));
    close(tps::infinityNorm(a), std::sqrt(5.0) + std::sqrt(13.0));
    auto rhs = column;
    factor.solve(rhs);
    complexClose(rhs[0], x[0]);
    complexClose(rhs[1], x[1]);
    rhs = row;
    factor.solveTranspose(rhs);
    complexClose(rhs[0], x[0]);
    complexClose(rhs[1], x[1]);
    const auto identity = a * factor.inverse();
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            complexClose(identity(i, j), C(i == j));
    tps::Matrix<C> multiple(2, 2);
    multiple(0, 0) = C(-7, 1);
    multiple(1, 0) = C(-13, 9);
    multiple(0, 1) = C(3, 7);
    multiple(1, 1) = C(14, 12);
    factor.solve(multiple);
    for (int i = 0; i < 2; ++i) {
        complexClose(multiple(i, 0), x[i]);
        complexClose(multiple(i, 1), y[i]);
    }
    multiple(0, 0) = C(-2, 11);
    multiple(1, 0) = C(-18, 9);
    multiple(0, 1) = C(12, -1);
    multiple(1, 1) = C(14, 17);
    factor.solveTranspose(multiple);
    for (int i = 0; i < 2; ++i) {
        complexClose(multiple(i, 0), x[i]);
        complexClose(multiple(i, 1), y[i]);
    }
    complexClose(a(0, 1), C(1, 1));
    complexClose(a(1, 0), C(2, -1));
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
        near(a(0, 0), 0);  // pivoting required; input unchanged
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
        } catch (const ::SingularMatrixError&) {
            singular = true;
        }
        if (!singular) throw std::runtime_error("singular matrix was accepted");
        bool mismatch = false;
        try {
            tps::Vector<double> wrong(3);
            (void)(x + wrong);
        } catch (const ::LogicalError&) {
            mismatch = true;
        }
        if (!mismatch) throw std::runtime_error("bad vector shape accepted");
        complexAlgebra<float>();
        complexAlgebra<double>();
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
