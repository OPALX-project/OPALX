#ifndef OPAL_TYPES_HH
#define OPAL_TYPES_HH

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
#include "Communicate/Operations.h"
#include "Ippl.h"

#include <chrono>
#include <iostream>
#include <string>

#include "PartBunch/datatypes.h"

#include "Utility/IpplTimings.h"

#include "Manager/PicManager.h"
#include "PartBunch/PartBunch.hpp"

using PartBunch_t = PartBunch<double, 3>;

// some typedefs

/*
template <unsigned Dim = 3>
using Mesh_t = ippl::UniformCartesian<double, Dim>;

template <typename T, unsigned Dim = 3>
using PLayout_t = typename ippl::ParticleSpatialLayout<T, Dim, Mesh_t<Dim>>;

template <unsigned Dim = 3>
using Centering_t = typename Mesh_t<Dim>::DefaultCentering;

template <unsigned Dim = 3>
using FieldLayout_t = ippl::FieldLayout<Dim>;

using size_type = ippl::detail::size_type;

template <typename T, unsigned Dim = 3>
using Vector = ippl::Vector<T, Dim>;

template <typename T, unsigned Dim = 3, class... ViewArgs>
using Field = ippl::Field<T, Dim, Mesh_t<Dim>, Centering_t<Dim>, ViewArgs...>;

template <unsigned Dim, class... ViewArgs>
using Field_t = Field<double, Dim, ViewArgs...>;

template <typename T = double, unsigned Dim = 3>
using ORB = ippl::OrthogonalRecursiveBisection<Field<double, Dim>, T>;

template <typename T>
using ParticleAttrib = ippl::ParticleAttrib<T>;

template <typename T, unsigned Dim = 3>
using Vector_t = ippl::Vector<T, Dim>;

template <unsigned Dim = 3, class... ViewArgs>
using Field_t = Field<double, Dim, ViewArgs...>;

template <typename T = double, unsigned Dim = 3, class... ViewArgs>
using VField_t = Field<Vector_t<T, Dim>, Dim, ViewArgs...>;
*/

template <typename T, unsigned Dim>
using Vector_t = ippl::Vector<T, Dim>;

typedef typename std::pair<Vector_t<double, 3>, Vector_t<double, 3>> VectorPair_t;

enum UnitState_t { units = 0, unitless = 1 };

/// \todo includes needs to reorganized
// #include "PartBunch/PartBunch.hpp"

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
