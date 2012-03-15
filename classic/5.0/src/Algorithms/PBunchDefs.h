//#include "Ippl.h"
#ifndef PBUNCHDEFS_H
#define PBUNCHDEFS_H

#include "Particle/IntCIC.h"
#include "Particle/IntNGP.h"
#include "Particle/IntSUDS.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Particle/ParticleAttrib.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/Centering.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"

typedef IntCIC  IntrplCIC_t;
typedef IntNGP  IntrplNGP_t;
typedef IntSUDS IntrplSUDS_t;

typedef ParticleSpatialLayout<double, 3>::ParticlePos_t Ppos_t;
typedef ParticleSpatialLayout<double, 3>::ParticleIndex_t PID_t;

typedef ParticleAttrib<double> Pscalar_t;

typedef InterpolatorTraits<double, 3, IntrplCIC_t>::Cache_t Pcache_t;

typedef UniformCartesian<3, double> Mesh_t;

typedef ParticleSpatialLayout<double, 3>::SingleParticlePos_t Vector_t;

typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;

typedef Cell Center_t;

typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t>       Field_t;
typedef Field<Vector_t, 3, Mesh_t, Center_t>     VField_t;

#endif
