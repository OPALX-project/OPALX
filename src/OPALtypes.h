#ifndef OPAL_TYPES_HH
#define OPAL_TYPES_HH

#include <mpi.h>
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <csignal>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>
#include "Communicate/Operations.h"
#include "Ippl.h"

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

typedef typename std::pair<Vector_t<double, 3>, Vector_t<double, 3>> VectorPair_t;

enum UnitState_t { units = 0, unitless = 1 };

#include "Algorithms/PartBunch.h"

// typedef PartBunch<PLayout_t<double, 3>, double, 3> PartBunch_t;

// euclidean norm
template <class T, unsigned D>
inline double euclidean_norm(const Vector_t<T, D>& v) {
    return std::sqrt(dot(v, v).apply());
}

// dot product
template <class T, unsigned D>
inline double dot(const Vector_t<T, D>& v, const Vector_t<T, D>& w) {
    T res = 0.0;
    for (unsigned i = 0; i < D; i++)
        res += v(i) * w(i);
    return std::sqrt(res);
}

template <typename T, class Op>
void allreduce(const T* input, T* output, int count, Op op) {
    MPI_Datatype type = get_mpi_datatype<T>(*input);

    MPI_Op mpiOp = get_mpi_op<Op>(op);

    MPI_Allreduce(const_cast<T*>(input), output, count, type, mpiOp, ippl::Comm->getCommunicator());
}

template <typename T, class Op>
void allreduce(const T& input, T& output, int count, Op op) {
    allreduce(&input, &output, count, op);
}

template <typename T, class Op>
void allreduce(T* inout, int count, Op op) {
    MPI_Datatype type = get_mpi_datatype<T>(*inout);

    MPI_Op mpiOp = get_mpi_op<Op>(op);

    MPI_Allreduce(MPI_IN_PLACE, inout, count, type, mpiOp, ippl::Comm->getCommunicator());
}

template <typename T, class Op>
void allreduce(T& inout, int count, Op op) {
    allreduce(&inout, count, op);
}

#endif