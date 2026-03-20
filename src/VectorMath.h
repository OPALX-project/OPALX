#ifndef OPALX_VECTOR_MATH_H
#define OPALX_VECTOR_MATH_H

#include "Ippl.h"

#include <utility>

template <typename T, unsigned Dim>
using Vector_t = ippl::Vector<T, Dim>;

using VectorPair_t = std::pair<Vector_t<double, 3>, Vector_t<double, 3>>;

// Euclidean norm for OPALX/IPPL vectors.
template <class T, unsigned D>
KOKKOS_INLINE_FUNCTION double euclidean_norm(const Vector_t<T, D>& v) {
    return Kokkos::sqrt(v.dot(v));
}

// Scalar product between two OPALX/IPPL vectors.
template <class T, unsigned D>
KOKKOS_INLINE_FUNCTION double dot(const Vector_t<T, D>& v, const Vector_t<T, D>& w) {
    return v.dot(w);
}

// Squared Euclidean norm for OPALX/IPPL vectors.
template <class T, unsigned D>
KOKKOS_INLINE_FUNCTION double dot(const Vector_t<T, D>& v) {
    return v.dot(v);
}

#endif
