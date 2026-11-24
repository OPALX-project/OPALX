#ifndef TPS_DEVICE_MAP_HPP
#define TPS_DEVICE_MAP_HPP
#include <Kokkos_Core.hpp>
#include <Kokkos_Array.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <cstddef>
#include <type_traits>
namespace tps {
/** @brief Real particle-coordinate types supported by prepared device maps. */
template <class T>
concept DeviceMapScalar = std::is_same_v<T, float> || std::is_same_v<T, double>;
/** @brief Non-owning strided device buffer; the caller owns extent/lifetime correctness.
 * @details A nonzero extent requires a nonnull pointer and positive stride.
 * Workspace must not alias input/output or another concurrent evaluation's workspace.
 * Actual allocation bounds, accessible memory space and representable index*stride
 * offsets are caller obligations; valid() checks only pointer and stride metadata.
 */
template <class T> class DeviceSpan {
  public:
    KOKKOS_INLINE_FUNCTION DeviceSpan() = default;
    KOKKOS_INLINE_FUNCTION DeviceSpan(T* pointer, std::size_t extent, std::size_t stride = 1)
        : pointer_(pointer), extent_(extent), stride_(stride) {}
    template <class U>
        requires std::is_convertible_v<U*, T*>
    KOKKOS_INLINE_FUNCTION DeviceSpan(const DeviceSpan<U>& other)
        : pointer_(other.data()), extent_(other.size()), stride_(other.stride()) {}
    KOKKOS_INLINE_FUNCTION T* data() const noexcept {
        return pointer_;
    }
    KOKKOS_INLINE_FUNCTION std::size_t size() const noexcept {
        return extent_;
    }
    KOKKOS_INLINE_FUNCTION std::size_t stride() const noexcept {
        return stride_;
    }
    KOKKOS_INLINE_FUNCTION bool valid() const noexcept {
        return stride_ > 0 && (extent_ == 0 || pointer_ != nullptr);
    }
    KOKKOS_INLINE_FUNCTION T& operator[](std::size_t index) const noexcept {
        return pointer_[index * stride_];
    }

  private:
    T* pointer_ = nullptr;
    std::size_t extent_ = 0, stride_ = 1;
};
/** @brief Adapt a rank-one Kokkos View/subview without allocating or taking ownership. */
template <class Data, class... Properties>
KOKKOS_INLINE_FUNCTION auto deviceSpan(const Kokkos::View<Data, Properties...>& view) {
    static_assert(Kokkos::View<Data, Properties...>::rank == 1,
                  "deviceSpan requires a rank-one view");
    using T = typename Kokkos::View<Data, Properties...>::value_type;
    return DeviceSpan<T>(view.data(), view.extent(0), view.stride(0));
}
template <class T, std::size_t N>
KOKKOS_INLINE_FUNCTION auto deviceSpan(Kokkos::Array<T, N>& array) {
    return DeviceSpan<T>(array.data(), N);
}
template <class T, std::size_t N>
KOKKOS_INLINE_FUNCTION auto deviceSpan(const Kokkos::Array<T, N>& array) {
    return DeviceSpan<const T>(array.data(), N);
}
/** @brief Device-safe result; these statuses are not physical aperture/validity checks. */
enum class EvaluationStatus {
    success,
    unpreparedMap,
    inputSizeMismatch,
    outputSizeMismatch,
    workspaceTooSmall,
    invalidBuffer,
    nonFiniteInput,
    nonFiniteResult
};
template <DeviceMapScalar T, class ExecutionSpace> class PreparedMap;
/** @brief Read-only device representation, with no host metadata or caches.
 * @details The default handle is unprepared. Valid handles come from PreparedMap.
 * Keep handles and all buffers alive until their execution-space work completes.
 */
template <DeviceMapScalar T, class MemorySpace> class DeviceMapView {
  public:
    using coefficient_view = Kokkos::View<const T**, Kokkos::LayoutRight, MemorySpace>;
    using index_view = Kokkos::View<const int*, MemorySpace>;
    KOKKOS_INLINE_FUNCTION DeviceMapView() = default;
    KOKKOS_INLINE_FUNCTION int input_size() const noexcept {
        return inputSize_;
    }
    KOKKOS_INLINE_FUNCTION int output_size() const noexcept {
        return outputSize_;
    }
    KOKKOS_INLINE_FUNCTION std::size_t monomial_count() const noexcept {
        return parents_.extent(0);
    }
    KOKKOS_INLINE_FUNCTION std::size_t workspace_size() const noexcept {
        return monomial_count() + static_cast<std::size_t>(outputSize_);
    }
    /// Output-major coefficients: (output, monomial node), including zeros.
    KOKKOS_INLINE_FUNCTION const coefficient_view& coefficients() const noexcept {
        return coefficients_;
    }
    /// Node zero is the constant one; every other parent precedes its child.
    KOKKOS_INLINE_FUNCTION const index_view& parents() const noexcept {
        return parents_;
    }
    /// Coordinate index in node = parent * input[variable]; node zero uses -1.
    KOKKOS_INLINE_FUNCTION const index_view& variables() const noexcept {
        return variables_;
    }

  private:
    template <DeviceMapScalar U, class ExecutionSpace> friend class PreparedMap;
    coefficient_view coefficients_;
    index_view parents_, variables_;
    int inputSize_ = 0, outputSize_ = 0;
};
/** @brief Evaluate one particle using caller-supplied storage, with no allocation.
 * @details Coordinates are local to the map charts. workspace_size() entries
 * are required: one per planned monomial plus one per output. Input and output
 * may overlap; workspace must be disjoint from both. Output is committed only
 * after all arithmetic succeeds, so every failure leaves output unchanged.
 * No exceptions, host cache access, synchronization or global precision lookup
 * occurs here. All spans must address memory accessible in the calling kernel.
 */
template <DeviceMapScalar T, class MemorySpace>
KOKKOS_INLINE_FUNCTION EvaluationStatus evaluate_one(
    const DeviceMapView<T, MemorySpace>& map, DeviceSpan<const std::type_identity_t<T>> input,
    DeviceSpan<T> output, DeviceSpan<T> workspace) noexcept {
    if (map.input_size() <= 0 || map.output_size() <= 0 || map.monomial_count() == 0)
        return EvaluationStatus::unpreparedMap;
    if (input.size() != static_cast<std::size_t>(map.input_size()))
        return EvaluationStatus::inputSizeMismatch;
    if (output.size() != static_cast<std::size_t>(map.output_size()))
        return EvaluationStatus::outputSizeMismatch;
    if (workspace.size() < map.workspace_size())
        return EvaluationStatus::workspaceTooSmall;
    if (!input.valid() || !output.valid() || !workspace.valid())
        return EvaluationStatus::invalidBuffer;
    for (int i = 0; i < map.input_size(); ++i)
        if (!Kokkos::isfinite(input[i]))
            return EvaluationStatus::nonFiniteInput;
    workspace[0] = T(1);
    const auto& parents = map.parents();
    const auto& variables = map.variables();
    for (std::size_t i = 1; i < map.monomial_count(); ++i) {
        const T value = workspace[parents(i)] * input[variables(i)];
        if (!Kokkos::isfinite(value))
            return EvaluationStatus::nonFiniteResult;
        workspace[i] = value;
    }
    const auto& coefficients = map.coefficients();
    for (int row = 0; row < map.output_size(); ++row) {
        T value = T(0);
        for (std::size_t i = 0; i < map.monomial_count(); ++i)
            value += coefficients(row, i) * workspace[i];
        if (!Kokkos::isfinite(value))
            return EvaluationStatus::nonFiniteResult;
        workspace[map.monomial_count() + static_cast<std::size_t>(row)] = value;
    }
    for (int row = 0; row < map.output_size(); ++row)
        output[row] = workspace[map.monomial_count() + static_cast<std::size_t>(row)];
    return EvaluationStatus::success;
}
} // namespace tps
#endif
