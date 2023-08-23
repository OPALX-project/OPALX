#ifndef PBUNCHDEFS_H
#define PBUNCHDEFS_H

#include "Ippl.h"

#include "Algorithms/Vektor.h"
#include "Interpolation/CIC.h"
#include "FieldLayout/FieldLayout.h"

#include "Meshes/UniformCartesian.h"

#include "Field/Field.h"

//typedef IntCIC  IntrplCIC_t;

//typedef ippl::ParticleSpatialLayout<double, 3>::ParticlePos_t Ppos_t;
//typedef ippl::ParticleSpatialLayout<double, 3>::ParticleIndex_t PID_t;

using Mesh_t = ippl::UniformCartesian<double, 3>;

typedef ippl::ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;

typedef Cell Center_t;


template <unsigned Dim = 3>
using FieldLayout_t = ippl::FieldLayout<Dim>;

template <unsigned Dim = 3, class... ViewArgs>
using Field_t = ippl::Field<double, Dim, ViewArgs...>;

template <typename T = double, unsigned Dim = 3, class... ViewArgs>
using VField_t = ippl::Field<ippl::Vector<T, Dim>, Dim, ViewArgs...>;

#endif