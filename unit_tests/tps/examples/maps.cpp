#include <Kokkos_Core.hpp>
#include <tps/tps.hpp>
#include <cmath>
#include <iostream>

int main() {
    Kokkos::ScopeGuard kokkos;

    using Series = tps::Tps<double>;
    Series::setGlobalTruncOrder(3);
    const auto q = Series::makeVariable(2, 0);
    const auto p = Series::makeVariable(2, 1);
    tps::VpsMap<double> drift(2, 2), kick(2, 2);
    drift[0] = q + 2.0 * p;
    drift[1] = p;
    kick[0] = q;
    kick[1] = p - 3.0 * q * q;
    const auto map = kick.substitute(drift, 3);
    tps::Vector<double> input(2, 0.0);
    input[0] = 0.1;
    input[1] = 0.2;
    const auto output = map.constantTerm(input);
    if (std::abs(output[0] - 0.5) > 2e-13 || std::abs(output[1] + 0.55) > 2e-13)
        return 1;
    const auto jacobian = map.linearTerms(input);
    if (std::abs(jacobian(0, 0) - 1) > 2e-13 || std::abs(jacobian(0, 1) - 2) > 2e-13 ||
        std::abs(jacobian(1, 0) + 3) > 2e-13 || std::abs(jacobian(1, 1) + 5) > 2e-13)
        return 2;
    std::cout << output[0] << ' ' << output[1] << '\n';
}
