/**
 * @file ParticleContainerTypes.h
 * @brief Neutral IPPL type aliases used by particle storage.
 */

#ifndef OPALX_PART_BUNCH_PARTICLE_CONTAINER_TYPES_H
#define OPALX_PART_BUNCH_PARTICLE_CONTAINER_TYPES_H

#include "Field/Field.h"
#include "FieldLayout/FieldLayout.h"
#include "Meshes/UniformCartesian.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Types/Vector.h"

template <unsigned Dim>
using Mesh_t = ippl::UniformCartesian<double, Dim>;

template <typename T, unsigned Dim>
using PLayout_t = ippl::ParticleSpatialLayout<T, Dim, Mesh_t<Dim>>;

template <unsigned Dim>
using Centering_t = typename Mesh_t<Dim>::DefaultCentering;

template <unsigned Dim>
using FieldLayout_t = ippl::FieldLayout<Dim>;

template <typename T, unsigned Dim>
using Vector_t = ippl::Vector<T, Dim>;

template <typename T, unsigned Dim, class... ViewArgs>
using Field = ippl::Field<T, Dim, Mesh_t<Dim>, Centering_t<Dim>, ViewArgs...>;

template <unsigned Dim, class... ViewArgs>
using Field_t = Field<double, Dim, ViewArgs...>;

template <typename T, unsigned Dim, class... ViewArgs>
using VField_t = Field<Vector_t<T, Dim>, Dim, ViewArgs...>;

#endif  // OPALX_PART_BUNCH_PARTICLE_CONTAINER_TYPES_H
