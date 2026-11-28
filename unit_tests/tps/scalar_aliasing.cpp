#include <Kokkos_Core.hpp>
#include <complex>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include "tps/Algebra/Tps.hpp"

namespace {
    void require(bool valid, const char* message) {
        if (!valid) throw std::runtime_error(message);
    }
    template <class T>
    void check() {
        using S = tps::Tps<T>;
        S::setGlobalTruncOrder(S::EXACT);
        for (int variables : {1, 2, 6}) {
            const int size = S(2, variables).getSize();
            std::vector<T> coefficients(size);
            for (int i = 0; i < size; ++i) {
                if constexpr (std::is_same_v<T, std::complex<double>>)
                    coefficients[i] = T(i + 2, 1 - i);
                else
                    coefficients[i] = T(i + 2);
            }
            for (int index : {0, 1, size - 1})
                for (int sharing : {0, 1, 2}) {
                    S p(2, variables);
                    for (int i = 0; i < size; ++i)
                        p[i] = coefficients[i];
                    std::optional<S> snapshot;
                    if (sharing == 1) snapshot.emplace(p);
                    const T& factor = p[index];
                    // Exercise both detachment before obtaining a coefficient reference
                    // and detachment inside *= while the reference points to shared data.
                    if (sharing == 2) snapshot.emplace(p);
                    p *= factor;
                    for (int i = 0; i < size; ++i) {
                        require(p.getCoefficient(i) == coefficients[i] * coefficients[index],
                                "aliased scalar multiplication");
                        if (snapshot)
                            require(snapshot->getCoefficient(i) == coefficients[i],
                                    "multiplication changed shared source");
                    }
                }
            for (const T factor : {T(0), T(1), T(-3)}) {
                S p(2, variables);
                for (int i = 0; i < size; ++i)
                    p[i] = coefficients[i];
                p *= factor;
                for (int i = 0; i < size; ++i)
                    require(p.getCoefficient(i) == coefficients[i] * factor, "ordinary factor");
            }
        }
    }
}  // namespace
int main() {
    Kokkos::ScopeGuard runtime;
    try {
        check<double>();
        check<std::complex<double>>();
        std::cout << "54 scalar alias cases and 18 ordinary-factor cases passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
