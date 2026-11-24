#include <Kokkos_Core.hpp>
#include <tps/tps.hpp>
#include <array>
#include <cmath>
#include <iostream>

int main() {
    Kokkos::ScopeGuard kokkos;

    tps::PhaseSpaceChart entrance;
    entrance.name = "entrance";
    entrance.coordinates = {{"q", "m", 1.0}, {"p", "normalized", 1.0}};
    entrance.reference = {0.0, 0.0};
    entrance.frame.name = "lab";
    entrance.convention = tps::CoordinateConvention::canonical;
    entrance.canonicalPairs = {{0, 1}};
    auto middle = entrance;
    middle.name = "middle";
    auto exit = entrance;
    exit.name = "exit";
    using Series = tps::Tps<double>;
    const auto q = Series::makeVariable(2, 0), p = Series::makeVariable(2, 1);
    const std::array<Series, 2> driftCoefficients{q + 2.0 * p, p};
    const std::array<Series, 2> kickCoefficients{q, p - 3.0 * q * q};
    const tps::Map<double> drift(tps::MapContext(entrance, middle, 3, {"drift-pass-1"}),
                                 driftCoefficients);
    const tps::Map<double> kick(tps::MapContext(middle, exit, 3, {"kick-pass-1"}),
                                kickCoefficients);
    // No global truncation initialization is needed for these map operations.
    const auto combined = tps::compose(kick, drift, 3); // drift acts first
    const auto result = combined.evaluate(std::array<double, 2>{0.1, 0.2});
    if (std::abs(result[0] - 0.5) > 2e-13 || std::abs(result[1] + 0.55) > 2e-13)
        return 1;
    std::cout << result[0] << ' ' << result[1] << '\n';
}
