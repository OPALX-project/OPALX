#ifndef IPPLInitializer_Header_Included
#define IPPLInitializer_Header_Included

#include "Initializer.h"
#include "Ippl.h"
#include "FFT/FFT.h"

#ifdef IPPL_USE_SINGLE_PRECISION
 #define T float
#else
 #define T double
#endif
#define Dim 3

namespace initializer {

class IPPLInitializer : protected Initializer {

public:
    
    // some useful typedefs
    typedef Vert                                          Center_t;
    typedef UniformCartesian<Ndim, real>                  Mesh_t;
    typedef CenteredFieldLayout<Ndim, Mesh_t, Center_t>   FieldLayout_t;
    typedef Field<dcomplex, Ndim, Mesh_t, Center_t>       CxField_t;
    typedef FFT<CCTransform, Ndim, T>                FFT_t;

    //typedef Field<real, Ndim, Mesh_t, Center_t>             RxField_t;
    
    CxField_t rho_m;
    Mesh_t *mesh;
    FieldLayout_t *flayout;
    FFT_t *fft_m;

    //void init_particles_ippl(ParticleAttrib< Vektor<double, 3> > &pos, 
    //template <class T, unsigned Dim> 
    void init_particles(ParticleAttrib< Vektor<T, Dim> > &pos, 
        ParticleAttrib< Vektor<T, Dim> > &vel, FieldLayout_t *layout, Mesh_t *mesh,
        InputParser& bdata, const char *tfName, MPI_Comm comm);

protected:

    void indens();
    void test_reality();
    void gravity_force(integer axis);
    void gravity_potential();
    void non_gaussian(integer axis);
    void set_particles(real z_in, real d_z, real ddot, integer axis);
    void output(integer axis, real* pos, real* vel);

};

}
#endif
