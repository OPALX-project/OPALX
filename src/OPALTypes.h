#ifndef OPAL_TYPES_HH
#define OPAL_TYPES_HH

#ifdef __CUDACC__
#pragma push_macro("__cpp_consteval")
#pragma push_macro("_NODISCARD")
#pragma push_macro("__builtin_LINE")

#define __cpp_consteval 201811L

#ifdef _NODISCARD
    #undef _NODISCARD
    #define _NODISCARD
#endif

#define consteval constexpr

#include <source_location>

#undef consteval
#pragma pop_macro("__cpp_consteval")
#pragma pop_macro("_NODISCARD")
#else
#include <source_location>
#endif

#include "Ippl.h"

#include <mpi.h>
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_Random.hpp>
#include <csignal>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>
#include <chrono>
#include <iostream>
#include <string>

#include "PartBunch/datatypes.h"
#include "Manager/PicManager.h"
#include "PartBunch/PartBunch.hpp"

using PartBunch_t = PartBunch<double, 3>;

template <typename T, unsigned Dim>
using Vector_t = ippl::Vector<T, Dim>;

typedef typename std::pair<Vector_t<double, 3>, Vector_t<double, 3>> VectorPair_t;

enum UnitState_t { units = 0, unitless = 1 };

// euclidean norm
template <class T, unsigned D>
KOKKOS_INLINE_FUNCTION double euclidean_norm(const Vector_t<T, D>& v) {
    return Kokkos::sqrt(dot(v, v).apply());
}

// dot products
template <class T, unsigned D>
KOKKOS_INLINE_FUNCTION double dot(const Vector_t<T, D>& v, const Vector_t<T, D>& w) {
    double res = 0.0;
    for (unsigned i = 0; i < D; i++)
        res += v(i) * w(i);
    return res;
}

template <class T, unsigned D>
KOKKOS_INLINE_FUNCTION double dot(const Vector_t<T, D>& v) {
    double res = 0.0;
    for (unsigned i = 0; i < D; i++)
        res += v(i) * v(i);
    return res;
}
#endif
