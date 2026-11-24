#include <Kokkos_Core.hpp>
#include "tps/tps.hpp"
#include <algorithm>
#include <complex>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <ranges>
#include <stdexcept>
using namespace tps;
static_assert(std::random_access_iterator<SliceIterator<int>>);
static_assert(std::random_access_iterator<ConstSliceIterator<int>>);
struct alignas(64) Tracked {
    static inline int alive = 0;
    static inline bool fail = false;
    int value = 0;
    Tracked(int v = 0) : value(v) {
        ++alive;
    }
    Tracked(const Tracked& rhs) : value(rhs.value) {
        ++alive;
    }
    Tracked& operator=(const Tracked& rhs) {
        if (fail)
            throw std::runtime_error("injected copy failure");
        value = rhs.value;
        return *this;
    }
    ~Tracked() {
        --alive;
    }
};
void require(bool ok, const char* message) {
    if (!ok)
        throw std::runtime_error(message);
}
int main() {
    Kokkos::ScopeGuard kokkos;

    try {
        {
            Tps<Tracked> x(2, 2);
            x[0] = Tracked(7);
            require(reinterpret_cast<std::uintptr_t>(&x[0]) % alignof(Tracked) == 0,
                    "coefficient alignment");
            auto copy = x;
            copy[0] = Tracked(9);
            require(static_cast<const Tps<Tracked>&>(x)[0].value == 7,
                    "nontrivial coefficient copy");
            Array1D<Tracked> a(2, Tracked(3)), b(3, Tracked(4));
            Array2D<Tracked> m(2, 2, Tracked(3)), n(1, 3, Tracked(4));
            auto shared = x;
            Tracked::fail = true;
            bool cloneFailed = false;
            try {
                (void)shared[0];
            } catch (const std::runtime_error&) {
                cloneFailed = true;
            }
            require(cloneFailed, "copy-on-write failure injection");
            Tracked::fail = false;
            x[0] = Tracked(8);
            require(static_cast<const Tps<Tracked>&>(shared)[0].value == 7,
                    "copy-on-write rollback");
            Tracked::fail = true;
            bool caught = false;
            try {
                a = b;
            } catch (const std::runtime_error&) {
                caught = true;
            }
            require(caught && a.size() == 2 && a[0].value == 3, "array assignment rollback");
            caught = false;
            try {
                a.resize(4);
            } catch (const std::runtime_error&) {
                caught = true;
            }
            require(caught && a.size() == 2, "resize rollback");
            caught = false;
            try {
                m = n;
            } catch (const std::runtime_error&) {
                caught = true;
            }
            require(caught && m.nrows() == 2 && m(0, 0).value == 3, "matrix assignment rollback");
            Tracked::fail = false;
        }
        require(Tracked::alive == 0, "nontrivial coefficient/array destruction");
        Tps<double> alias(1, 2);
        alias[0] = 7;
        alias[1] = 2;
        alias = alias[0];
        require(alias[0] == 7 && alias.getVariables() == 0, "aliased scalar assignment");
        const double coefficients[] = {1.0, 1.0, 0.5};
        const auto variable = Tps<double>::makeVariable(1, 0);
        const auto approximation = variable.taylor(coefficients);
        require(approximation[0] == 1 && approximation[1] == 1 && approximation[2] == 0.5,
                "span Taylor interface");
        bool spanRejected = false;
        try {
            (void)variable.taylor({});
        } catch (const CLRangeError&) {
            spanRejected = true;
        }
        require(spanRejected, "empty Taylor span rejection");
        using C = std::complex<double>;
        Tps<C>::setGlobalTruncOrder(3);
        auto x = Tps<C>::makeVariable(2, 0);
        x[0] = C(1, 2);
        auto square = x * x;
        require(std::abs(square[0] - C(-3, 4)) < 1e-13, "complex constant");
        require(std::abs(square[1] - C(2, 4)) < 1e-13, "complex linear coefficient");
        auto inverse = x.inverse();
        require(std::abs((x * inverse)[0] - C(1)) < 1e-13, "complex inversion");
        auto e = exp(x);
        require(std::abs(e[0] - std::exp(C(1, 2))) < 1e-13, "complex exponential");
        using Inner = Tps<double>;
        auto outer = Tps<Inner>::makeVariable(1, 0);
        auto nested = (outer + Inner(2.0)) * (outer + Inner(2.0));
        require(nested[0][0] == 4 && nested[1][0] == 4 && nested[2][0] == 1,
                "nested coefficient arithmetic");
        Array2D<int> matrix(4, 3, 0);
        matrix(0, 2) = 4;
        matrix(1, 2) = 1;
        matrix(2, 2) = 3;
        matrix(3, 2) = 2;
        std::ranges::sort(matrix.col_begin(2), matrix.col_end(2));
        for (int i = 0; i < 4; ++i)
            require(matrix(i, 2) == i + 1, "C++20 ranges column sort");
        bool caught = false;
        try {
            (void)matrix[-1];
        } catch (const std::exception&) {
            caught = true;
        }
        require(caught, "standard exception interoperability");
        std::cout << "Modern storage and C++20 concepts passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
