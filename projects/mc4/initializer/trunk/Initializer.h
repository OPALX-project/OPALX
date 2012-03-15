#ifndef Initializer_Header_Included
#define Initializer_Header_Included

#include "TypesAndDefs.h"

#ifdef IPPLCREATE
#include "Ippl.h"

#ifdef IPPL_USE_SINGLE_PRECISION
 #define T float
#else
 #define T double
#endif
#define Dim 3
#endif

#ifdef INITIALIZER_NAMESPACE
namespace initializer {
#endif

class Basedata;

void init_particles(float* pos_x, float* pos_y, float* pos_z, 
                    float* vel_x, float* vel_y, float* vel_z, 
      		        Basedata& bdata, const char *tfName, MPI_Comm comm);

#ifdef IPPLCREATE

//void init_particles_ippl(ParticleAttrib< Vektor<double, 3> > &pos, 
//template <class T, unsigned Dim> 
void init_particles_ippl(ParticleAttrib< Vektor<T, Dim> > &pos, 
        ParticleAttrib< Vektor<T, Dim> > &vel, 
        Basedata& bdata, const char *tfName, MPI_Comm comm);

#endif

#ifdef INITIALIZER_NAMESPACE
}
#endif

#endif
