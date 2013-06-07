#ifndef PBUNCHDEFS_H
#define PBUNCHDEFS_H

#include "Algorithms/Vektor.h"
#include "Particle/IntCIC.h"
#include "Particle/IntNGP.h"
#include "Particle/IntSUDS.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Particle/ParticleAttrib.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/Centering.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"
#include "FFT/FFT.h"

typedef IntCIC  IntrplCIC_t;
typedef IntNGP  IntrplNGP_t;
typedef IntSUDS IntrplSUDS_t;

typedef ParticleSpatialLayout<double, 3>::ParticlePos_t Ppos_t;
typedef ParticleSpatialLayout<double, 3>::ParticleIndex_t PID_t;

typedef ParticleAttrib<double> Pscalar_t;

typedef InterpolatorTraits<double, 3, IntrplCIC_t>::Cache_t Pcache_t;

typedef UniformCartesian<3, double> Mesh_t;

//typedef ParticleSpatialLayout<double, 3>::SingleParticlePos_t Vector_t;

typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;

typedef Cell Center_t;

typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t>       Field_t;
typedef Field<Vector_t, 3, Mesh_t, Center_t>     VField_t;

typedef Field<int, 3, Mesh_t, Center_t>          IField_t;
typedef Field<dcomplex, 3, Mesh_t, Center_t>     CxField_t;
typedef FFT<RCTransform, 3, double>              FFT_t;
typedef FFT<SineTransform, 3, double>            SINE_t;

#endif
