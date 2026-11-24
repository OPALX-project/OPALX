#ifndef TPS_PREPARED_MAP_HPP
#define TPS_PREPARED_MAP_HPP
#include "tps/DeviceMap.hpp"
#include "tps/Map.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>
namespace tps {
/** @brief Host-prepared immutable coefficients and monomial recurrence in a chosen memory space.
 * @details Supports float/double Map objects. Only nonzero terms and their parent
 * monomials are retained; each monomial is computed once for all outputs.
 * Construction uploads and fences its execution space before returning. Source
 * maps can subsequently be destroyed. Chart metadata stays on the host.
 * Preparation still uses host algebra/cache access and needs its synchronization.
 */
template <DeviceMapScalar T, class ExecutionSpace = Kokkos::DefaultExecutionSpace>
class PreparedMap {
  public:
    static_assert(Kokkos::is_execution_space<ExecutionSpace>::value,
                  "PreparedMap requires an execution space");
    using execution_space = ExecutionSpace;
    using memory_space = typename ExecutionSpace::memory_space;
    using device_view_type = DeviceMapView<T, memory_space>;
    /** @brief Snapshot the nonzero polynomial terms and upload a shared monomial plan.
     * @throws std::invalid_argument A coefficient is not finite.
     * @throws std::length_error The monomial count cannot be represented by an int.
     * @details Requires an active Kokkos runtime. Allocation failures propagate.
     * Construction fences execution before publishing the read-only device handle.
     */
    explicit PreparedMap(const Map<T>& source, const ExecutionSpace& execution = ExecutionSpace{})
        : context_(source.context()) {
        const auto inputs = static_cast<int>(context_.input().coordinates.size());
        const auto outputs = static_cast<int>(context_.output().coordinates.size());
        const auto plan = buildPlan(source, inputs, outputs);
        upload(plan, inputs, outputs, execution);
    }
    [[nodiscard]] device_view_type device_view() const {
        return view_;
    }
    [[nodiscard]] const MapContext& context() const noexcept {
        return context_;
    }
    [[nodiscard]] std::size_t workspace_size() const noexcept {
        return view_.workspace_size();
    }

  private:
    using Exponents = std::vector<int>;
    struct Term {
        int output;
        int node;
        T coefficient;
    };
    // Node zero is 1. Each subsequent node has an earlier parent, making the
    // recurrence topological and reusable across every output component.
    struct HostPlan {
        std::map<Exponents, int> indices;
        std::vector<int> parents{0};
        std::vector<int> variables{-1};
        std::vector<Term> terms;
    };
    // Walk back to an existing ancestor, then append missing nodes in parent-first order.
    static int ensureNode(HostPlan& plan, Exponents powers, int inputs) {
        auto& indices = plan.indices;
        auto& parents = plan.parents;
        auto& variables = plan.variables;
        std::vector<Exponents> pending;
        while (!indices.contains(powers)) {
            pending.push_back(powers);
            int variable = inputs - 1;
            while (powers[variable] == 0)
                --variable;
            --powers[variable];
        }
        int parent = indices.at(powers);
        for (auto it = pending.rbegin(); it != pending.rend(); ++it) {
            if (parents.size() >= static_cast<std::size_t>(std::numeric_limits<int>::max()))
                throw std::length_error("PreparedMap: too many monomials");
            int variable = inputs - 1;
            while ((*it)[variable] == 0)
                --variable;
            int node = static_cast<int>(parents.size());
            parents.push_back(parent);
            variables.push_back(variable);
            indices.emplace(*it, node);
            parent = node;
        }
        return parent;
    }

    static HostPlan buildPlan(const Map<T>& source, int inputs, int outputs) {
        HostPlan plan;
        plan.indices.emplace(Exponents(inputs, 0), 0);
        for (int row = 0; row < outputs; ++row) {
            const auto& polynomial = source.polynomial()[row];
            for (int index = 0; index < polynomial.getSize(); ++index) {
                const T coefficient = polynomial[index];
                if (!std::isfinite(coefficient))
                    throw std::invalid_argument("PreparedMap: nonfinite coefficient");
                if (coefficient == T(0))
                    continue;
                Exponents powers(inputs, 0);
                if (index != 0) {
                    const auto& monomial = polynomial.getExponents(index);
                    for (int j = 0; j < inputs; ++j)
                        powers[j] = monomial[j];
                }
                plan.terms.push_back(
                    {row, ensureNode(plan, std::move(powers), inputs), coefficient});
            }
        }
        return plan;
    }
    // Detached host mirrors are deliberately used even when execution memory is
    // host-accessible: initialization and copies must respect the selected queue.
    void upload(const HostPlan& plan, int inputs, int outputs, const ExecutionSpace& execution) {
        Kokkos::View<T**, Kokkos::LayoutRight, memory_space> coefficients(
            Kokkos::view_alloc(execution, Kokkos::WithoutInitializing, "tps prepared coefficients"),
            outputs, plan.parents.size());
        Kokkos::View<int*, memory_space> deviceParents(
            Kokkos::view_alloc(execution, Kokkos::WithoutInitializing, "tps monomial parents"),
            plan.parents.size());
        Kokkos::View<int*, memory_space> deviceVariables(
            Kokkos::view_alloc(execution, Kokkos::WithoutInitializing, "tps monomial variables"),
            plan.variables.size());
        auto hostCoefficients = Kokkos::create_mirror(Kokkos::WithoutInitializing, coefficients);
        auto hostParents = Kokkos::create_mirror(Kokkos::WithoutInitializing, deviceParents);
        auto hostVariables = Kokkos::create_mirror(Kokkos::WithoutInitializing, deviceVariables);
        // Detached, uninitialized mirrors avoid cross-queue initialization races.
        std::fill_n(hostCoefficients.data(), hostCoefficients.size(), T(0));
        for (const auto& term : plan.terms)
            hostCoefficients(term.output, term.node) = term.coefficient;
        for (std::size_t i = 0; i < plan.parents.size(); ++i) {
            hostParents(i) = plan.parents[i];
            hostVariables(i) = plan.variables[i];
        }
        Kokkos::deep_copy(execution, coefficients, hostCoefficients);
        Kokkos::deep_copy(execution, deviceParents, hostParents);
        Kokkos::deep_copy(execution, deviceVariables, hostVariables);
        execution.fence("TPS prepared map upload");
        view_.coefficients_ = coefficients;
        view_.parents_ = deviceParents;
        view_.variables_ = deviceVariables;
        view_.inputSize_ = inputs;
        view_.outputSize_ = outputs;
    }
    MapContext context_;
    device_view_type view_;
};
} // namespace tps
#endif
