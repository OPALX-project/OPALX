#include <tps/tps.hpp>
#include <array>
#include <cmath>
#include <iostream>

int main() {
    Kokkos::ScopeGuard runtime;
    tps::PhaseSpaceChart chart;
    chart.name = "canonical-plane";
    chart.frame.name = "lab";
    chart.coordinates = {{"q", "m", 1}, {"p", "normalized", 1}};
    chart.reference = {0, 0};
    chart.convention = tps::CoordinateConvention::canonical;
    chart.canonicalPairs = {{0, 1}};
    using Series = tps::Tps<double>;
    const auto q = Series::makeVariable(2, 0), p = Series::makeVariable(2, 1);
    const std::array<Series, 2> coefficients{q + 2.0 * p, p - 3.0 * (q + 2.0 * p) * (q + 2.0 * p)};
    const tps::Map<double> host(tps::MapContext(chart, chart, 2), coefficients);
    const tps::PreparedMap<double> prepared(host);
    const auto map = prepared.device_view();
    Kokkos::View<double*> output("output", 2), workspace("workspace", map.workspace_size());
    Kokkos::View<int*> status("status", 1);
    // This launch only demonstrates the single-particle primitive. OPALX will
    // choose its batch policy, particle Views and per-particle workspace.
    Kokkos::parallel_for(
        "single prepared particle", 1, KOKKOS_LAMBDA(int) {
            const Kokkos::Array<double, 2> input{0.1, 0.2};
            status(0) = static_cast<int>(tps::evaluate_one(
                map, tps::deviceSpan(input), tps::deviceSpan(output), tps::deviceSpan(workspace)));
        });
    const auto result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output);
    const auto resultStatus = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), status);
    if (resultStatus(0) != static_cast<int>(tps::EvaluationStatus::success))
        return 1;
    if (std::abs(result(0) - 0.5) > 2e-13 || std::abs(result(1) + 0.55) > 2e-13)
        return 2;
    std::cout << result(0) << ' ' << result(1) << '\n';
}
