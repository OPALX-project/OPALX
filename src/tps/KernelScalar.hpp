#ifndef TPS_KERNEL_SCALAR_HPP
#define TPS_KERNEL_SCALAR_HPP

#include <Kokkos_Complex.hpp>
#include <complex>
#include <type_traits>

namespace tps {
    namespace detail {
        template <class T>
        inline constexpr bool isStdComplex = false;
        template <class T>
        inline constexpr bool isStdComplex<std::complex<T>> = true;
    }  // namespace detail

    /** @brief Exclude host-only std::complex from Kokkos linear algebra.
     * @details Complex Vector/Matrix scalars must use Kokkos::complex, even with
     * HostSpace storage: the kernels are compiled for host and device execution.
     * Other scalar types must support the operations required by each kernel.
     * This restriction does not apply to host Tps coefficients or generic arrays.
     * @test tps.kokkos_linear_algebra checks float/double complex algebra and solves.
     * @snippet kokkos_linear_algebra.cpp complex kernel semantics
     */
    template <class T>
    concept LinearAlgebraScalar = !detail::isStdComplex<std::remove_cv_t<T>>;
}  // namespace tps

#endif
