#ifndef Initializer_Header_Included
#define Initializer_Header_Included

#include <mpi.h>
#include "TypesAndDefs.h"
#include "Parallelization.h"
#include "DataBase.h"
#include "MT_Random.h"
#include "Cosmology.h"
#include "PerfMon.h"
#include "distribution.h"
#include "InputParser.h"

#ifdef USENAMESPACE
namespace initializer {
#endif

#ifndef NOFFTW
#ifdef FFTW2
#include "fftw_mpi.h"
#endif
#ifdef FFTW3
#include <fftw3-mpi.h>
#endif
#endif

#define pi M_PI
#define Ndim 3  /* Dimension of the problem */

class Initializer {

public:

    //Initializer() {}
    //~Initializer() {
        ////FIXME: 
    //}
    
    void init_particles(real* pos_x, real* pos_y, real* pos_z, 
                    real* vel_x, real* vel_y, real* vel_z, 
		            InputParser& par, const char *tfName);

//FIXME: later distinguish between private and protected
protected:

    ParallelTools Parallel;
    GlobalStuff DataBase;
    CosmoClass Cosmology;

    real *Pk;
    my_fftw_complex *rho;
    real *buf_pos, *buf_vel; // buffers for neutrino calculation

    long My_Ng;  // lenght of the above arrays = number of my grid points

    integer nx1, nx2, ny1, ny2, nz1, nz2; // end points of my domain
    integer ngx, ngy, ngz;                // # of zones in my domain

    integer MyPE, MasterPE, NumPEs;   // MPI stuff

#ifndef NOFFTW
#ifdef FFTW2
    fftwnd_mpi_plan plan, kplan;
    int lnx, lxs, lnyt, lyst, lsize;
#endif
#ifdef FFTW3
    fftw_plan plan, kplan;
#endif
#endif

    void CreateMPI_FFTW_COMPLEX(MPI_Datatype *MPI_FFTW_COMPLEX);
    void initspec();
    void indens();
    void test_reality();
    void gravity_force(integer axis);
    void gravity_potential();
    void non_gaussian(integer axis);
    void set_particles(real z_in, real d_z, real ddot, integer axis);
    void output(integer axis, real* pos, real* vel);
    void exchange_buffers(real *pos, real *vel);
    void inject_neutrinos(real *pos, real *vel, integer axis);
    void thermal_vel(real *vx, real *vy, real *vz);
    void apply_periodic(real *pos);

};

#ifdef USENAMESPACE
}
#endif

#endif
