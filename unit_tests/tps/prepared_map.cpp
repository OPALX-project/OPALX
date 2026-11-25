#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <tps/tps.hpp>
#include <type_traits>
namespace {
    std::atomic<int> allocations{0};
    void allocationEvent(Kokkos::Tools::SpaceHandle, const char*, const void*, std::uint64_t) {
        ++allocations;
    }
    void require(bool ok, const char* message) {
        if (!ok) throw std::runtime_error(message);
    }
    void near(double a, double b, double tolerance = 2e-13) {
        require(std::isfinite(a) && std::abs(a - b) <= tolerance,
                "prepared map numerical mismatch");
    }
    tps::PhaseSpaceChart chart(int count) {
        tps::PhaseSpaceChart c;
        c.name       = "test";
        c.frame.name = "lab";
        for (int i = 0; i < count; ++i) {
            c.coordinates.push_back({"x" + std::to_string(i), "unit", 1});
            c.reference.push_back(0);
        }
        return c;
    }
    template <class T>
    auto sourceMap() {
        auto in = chart(3), out = chart(2);
        using S = tps::Tps<T>;
        auto x = S::makeVariable(3, 0), y = S::makeVariable(3, 1), z = S::makeVariable(3, 2);
        // x^5 has zero coefficients in its lower powers; plan must add missing parents.
        auto fifth = S::makeVarPower(3, 0, 5);
        std::array<S, 2> components{
                T(2) + x + T(2) * y + T(3) * x * y + fifth, T(-1) + z * z - T(2) * x * y};
        return tps::Map<T>(tps::MapContext(in, out, 5), components);
    }
    template <class T>
    void compareKernel() {
        using Exec        = Kokkos::DefaultExecutionSpace;
        const auto source = sourceMap<T>();
        const tps::PreparedMap<T, Exec> prepared(source);
        const auto map = prepared.device_view();
        static_assert(
                !std::is_assignable_v<typename decltype(map)::coefficient_view::reference_type, T>);
        static_assert(
                !std::is_assignable_v<typename decltype(map)::index_view::reference_type, int>);
        require(map.input_size() == 3 && map.output_size() == 2, "prepared dimensions");
        require(prepared.context().input() == source.context().input(), "context retained on host");
        require(map.monomial_count() == 10, "monomial pruning/sharing (including x^5 parents)");
        auto parent = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), map.parents());
        for (std::size_t i = 1; i < map.monomial_count(); ++i)
            require(parent(i) >= 0 && static_cast<std::size_t>(parent(i)) < i,
                    "non-topological plan");
        constexpr int count = 37;
        // LayoutLeft makes per-particle row subviews strided, as in particle SoA.
        Kokkos::View<T**, Kokkos::LayoutLeft> input("input", count, 3), output("output", count, 2);
        Kokkos::View<T**, Kokkos::LayoutLeft> workspace("workspace", count, map.workspace_size());
        Kokkos::View<int*> status("status", count);
        auto hostInput = Kokkos::create_mirror_view(input);
        for (int i = 0; i < count; ++i) {
            hostInput(i, 0) = T(i - 18) / T(40);
            hostInput(i, 1) = T(i % 7 - 3) / T(10);
            hostInput(i, 2) = T(i % 5) / T(9);
        }
        Kokkos::deep_copy(input, hostInput);
        // The evaluator must not observe a subsequently changed host global order.
        const int old = tps::Tps<T>::getGlobalTruncOrder();
        tps::Tps<T>::setGlobalTruncOrder(0);
        //! [prepared map evaluation]
        Kokkos::parallel_for(
                "prepared map comparison", count, KOKKOS_LAMBDA(int i) {
                    auto in   = Kokkos::subview(input, i, Kokkos::ALL),
                         out  = Kokkos::subview(output, i, Kokkos::ALL);
                    auto work = Kokkos::subview(workspace, i, Kokkos::ALL);
                    status(i) = static_cast<int>(tps::evaluate_one(
                            map, tps::deviceSpan(in), tps::deviceSpan(out), tps::deviceSpan(work)));
                });
        Kokkos::fence();
        tps::Tps<T>::setGlobalTruncOrder(old);
        auto result   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output);
        auto statuses = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), status);
        for (int i = 0; i < count; ++i) {
            require(statuses(i) == static_cast<int>(tps::EvaluationStatus::success),
                    "kernel status");
            std::array<T, 3> point{hostInput(i, 0), hostInput(i, 1), hostInput(i, 2)};
            auto reference   = source.evaluate(point);
            double tolerance = std::is_same_v<T, float> ? 3e-6 : 2e-13;
            for (int j = 0; j < 2; ++j)
                near(result(i, j), reference[j], tolerance);
            double x = point[0], y = point[1], z = point[2];
            near(result(i, 0), 2 + x + 2 * y + 3 * x * y + std::pow(x, 5), tolerance);
            near(result(i, 1), -1 + z * z - 2 * x * y, tolerance);
        }
        //! [prepared map evaluation]
        // Copying the handle retains read-only data after source/prepared owner destruction.
        auto surviving = [] {
            tps::PreparedMap<T> temporary(sourceMap<T>());
            return temporary.device_view();
        }();
        Kokkos::parallel_for(
                "prepared in-place", 1, KOKKOS_LAMBDA(int) {
                    Kokkos::Array<T, 3> point{T(0.2), T(0.3), T(0.4)};
                    auto work    = Kokkos::subview(workspace, 0, Kokkos::ALL);
                    status(0)    = static_cast<int>(tps::evaluate_one(
                            surviving, tps::deviceSpan(point), tps::DeviceSpan<T>(point.data(), 2),
                            tps::deviceSpan(work)));
                    output(0, 0) = point[0];
                    output(0, 1) = point[1];
                });
        auto inplace = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output);
        near(inplace(0, 0), 2 + 0.2 + 0.6 + 0.18 + 0.00032,
             std::is_same_v<T, float> ? 3e-6 : 2e-13);
        near(inplace(0, 1), -1 + 0.16 - 0.12, std::is_same_v<T, float> ? 3e-6 : 2e-13);
        // Constant and zero maps retain dimension but need only the constant monomial.
        std::array<tps::Tps<T>, 2> constants{tps::Tps<T>(T(4)), tps::Tps<T>(T(0))};
        tps::Map<T> constant(tps::MapContext(chart(3), chart(2), 0), constants);
        tps::PreparedMap<T> cp(constant);
        auto cv = cp.device_view();
        require(cv.monomial_count() == 1, "constant plan");
        Kokkos::parallel_for(
                "constant prepared map", 1, KOKKOS_LAMBDA(int) {
                    Kokkos::Array<T, 3> point{T(3), T(2), T(1)}, work{};
                    Kokkos::Array<T, 2> out{};
                    status(0)    = static_cast<int>(tps::evaluate_one(
                            cv, tps::deviceSpan(point), tps::deviceSpan(out),
                            tps::deviceSpan(work)));
                    output(0, 0) = out[0];
                    output(0, 1) = out[1];
                });
        auto constantResult = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), output);
        near(constantResult(0, 0), 4);
        near(constantResult(0, 1), 0);
    }
    void statusChecks() {
        using E = tps::EvaluationStatus;
        using S = tps::Tps<double>;
        auto x  = S::makeVariable(1, 0);
        std::array<S, 2> components{x, x * x};
        tps::Map<double> source(tps::MapContext(chart(1), chart(2), 2), components);
        tps::PreparedMap<double> prepared(source);
        const auto map = prepared.device_view();
        Kokkos::View<int*> statuses("failure statuses", 10);
        Kokkos::View<double**> results("unchanged outputs", 10, 2);
        Kokkos::deep_copy(results, -99.0);
        Kokkos::parallel_for(
                "prepared failures", 10, KOKKOS_LAMBDA(int mode) {
                    Kokkos::Array<double, 1> point{0.5};
                    Kokkos::Array<double, 2> out{-99, -99};
                    Kokkos::Array<double, 5> work{};
                    auto in        = tps::deviceSpan(point);
                    auto output    = tps::deviceSpan(out);
                    auto workspace = tps::deviceSpan(work);
                    if (mode == 1) in = tps::DeviceSpan<double>(point.data(), 0);
                    if (mode == 2) output = tps::DeviceSpan<double>(out.data(), 1);
                    if (mode == 3)
                        workspace = tps::DeviceSpan<double>(work.data(), map.workspace_size() - 1);
                    if (mode == 4) in = tps::DeviceSpan<double>(nullptr, 1);
                    if (mode == 5) in = tps::DeviceSpan<double>(point.data(), 1, 0);
                    if (mode == 6) point[0] = Kokkos::Experimental::infinity_v<double>;
                    if (mode == 7) point[0] = Kokkos::Experimental::quiet_NaN_v<double>;
                    if (mode == 8) point[0] = 1e308;
                    if (mode == 9) workspace = tps::DeviceSpan<double>(nullptr, work.size());
                    if (mode == 0) {
                        decltype(map) empty;
                        statuses(mode) =
                                static_cast<int>(tps::evaluate_one(empty, in, output, workspace));
                    } else
                        statuses(mode) =
                                static_cast<int>(tps::evaluate_one(map, in, output, workspace));
                    results(mode, 0) = out[0];
                    results(mode, 1) = out[1];
                });
        auto host    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), statuses);
        auto outputs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), results);
        std::array<E, 10> expected{E::unpreparedMap,      E::inputSizeMismatch,
                                   E::outputSizeMismatch, E::workspaceTooSmall,
                                   E::invalidBuffer,      E::invalidBuffer,
                                   E::nonFiniteInput,     E::nonFiniteInput,
                                   E::nonFiniteResult,    E::invalidBuffer};
        for (int i = 0; i < 10; ++i) {
            require(host(i) == static_cast<int>(expected[i]), "wrong failure status");
            near(outputs(i, 0), -99);
            near(outputs(i, 1), -99);
        }
        // Coefficients finite on upload, even if they occur above the constant term.
        components[0][1] = std::numeric_limits<double>::infinity();
        tps::Map<double> invalid(tps::MapContext(chart(1), chart(2), 2), components);
        bool rejected = false;
        try {
            tps::PreparedMap<double> bad(invalid);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "nonfinite coefficient accepted");
        // Overflow during accumulation (not just during monomial construction).
        components[0]    = S(2.0);
        components[1]    = S(1, 1);
        components[1][1] = 1e308;
        tps::PreparedMap<double, Kokkos::DefaultHostExecutionSpace> huge(
                tps::Map<double>(tps::MapContext(chart(1), chart(2), 2), components));
        auto hv = huge.device_view();
        Kokkos::Array<double, 1> p{2};
        Kokkos::Array<double, 2> out{-99, -99};
        Kokkos::Array<double, 5> work{};
        require(tps::evaluate_one(
                        hv, tps::deviceSpan(p), tps::deviceSpan(out), tps::deviceSpan(work))
                        == E::nonFiniteResult,
                "accumulation overflow");
        near(out[0], -99);
        near(out[1], -99);
        // Warm the host-space primitive, then check its Kokkos allocation events.
        auto stable =
                tps::PreparedMap<double, Kokkos::DefaultHostExecutionSpace>(source).device_view();
        tps::DeviceSpan<const double> input(p.data(), 1);
        p[0] = 0.2;
        (void)tps::evaluate_one(stable, input, tps::deviceSpan(out), tps::deviceSpan(work));
        Kokkos::Tools::Experimental::set_allocate_data_callback(allocationEvent);
        allocations = 0;
        auto result = tps::evaluate_one(stable, input, tps::deviceSpan(out), tps::deviceSpan(work));
        const auto count = allocations.load();
        Kokkos::Tools::Experimental::set_allocate_data_callback(nullptr);
        require(result == E::success, "allocation-check evaluation failed");
        require(count == 0, "evaluator allocated Kokkos memory");
    }
}  // namespace
int main() {
    Kokkos::ScopeGuard runtime;
    try {
        compareKernel<double>();
        compareKernel<float>();
        statusChecks();
        std::cout << "Prepared map device evaluation checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
