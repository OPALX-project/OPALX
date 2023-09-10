#ifndef OPAL_TYPES_HH
#define OPAL_TYPES_HH

#include "Ippl.h"

#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <csignal>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include "Solver/ElectrostaticsCG.h"
#include "Solver/FFTPeriodicPoissonSolver.h"
#include "Solver/FFTPoissonSolver.h"
#include "Solver/P3MSolver.h"

// some typedefs
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

// heFFTe does not support 1D FFTs, so we switch to CG in the 1D case
template <typename T = double, unsigned Dim = 3>
using CGSolver_t = ippl::ElectrostaticsCG<Field<T, Dim>, Field_t<Dim>>;

using ippl::detail::ConditionalType, ippl::detail::VariantFromConditionalTypes;

template <typename T = double, unsigned Dim = 3>
using FFTSolver_t = ConditionalType<
    Dim == 2 || Dim == 3, ippl::FFTPeriodicPoissonSolver<VField_t<T, Dim>, Field_t<Dim>>>;

template <typename T = double, unsigned Dim = 3>
using P3MSolver_t = ConditionalType<Dim == 3, ippl::P3MSolver<VField_t<T, Dim>, Field_t<Dim>>>;

template <typename T = double, unsigned Dim = 3>
using OpenSolver_t =
    ConditionalType<Dim == 3, ippl::FFTPoissonSolver<VField_t<T, Dim>, Field_t<Dim>>>;

template <typename T = double, unsigned Dim = 3>
using Solver_t = VariantFromConditionalTypes<
    CGSolver_t<T, Dim>, FFTSolver_t<T, Dim>, P3MSolver_t<T, Dim>, OpenSolver_t<T, Dim>>;

typedef typename std::pair<Vector_t<double, 3>, Vector_t<double, 3>> VectorPair_t;

enum UnitState_t { units = 0, unitless = 1 };

#include "Algorithms/PartBunch.h"

typedef PartBunch<PLayout_t<double, 3>, double, 3> PartBunch_t;

// euclidean norm
template <class T, unsigned D>
inline double euclidean_norm(const Vector_t<T, D>& v) {
    return std::sqrt(dot(v, v).apply());
}

#endif