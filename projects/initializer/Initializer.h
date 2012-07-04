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

#include <memory>
#include <climits>

#ifdef USENAMESPACE
namespace initializer {
  using namespace std;
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

    void init_particles(std::unique_ptr<real[]> &pos_x, std::unique_ptr<real[]> &pos_y, std::unique_ptr<real[]> &pos_z, 
			std::unique_ptr<real[]> &vel_x, std::unique_ptr<real[]> &vel_y, std::unique_ptr<real[]> &vel_z, 
			InputParser& par, const char *tfName);

//FIXME: later distinguish between private and protected

protected:

    ParallelTools Parallel;
    GlobalStuff DataBase;
    CosmoClass Cosmology;

    std::unique_ptr<real[]> Pk;

    my_fftw_complex *rho;

    std::unique_ptr<real[]> buf_pos;
    std::unique_ptr<real[]> buf_vel;

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
    void output(integer axis, std::unique_ptr<real[]> &pos, std::unique_ptr<real[]> &vel);
    void exchange_buffers(std::unique_ptr<real[]> &pos, std::unique_ptr<real[]> &vel);
    void inject_neutrinos(std::unique_ptr<real[]> &pos, std::unique_ptr<real[]> &vel, integer axis);
    void thermal_vel(std::unique_ptr<real[]> &vx, std::unique_ptr<real[]> &vy, std::unique_ptr<real[]> &vz);
    void apply_periodic(std::unique_ptr<real[]> &pos);

};

#ifdef USENAMESPACE
}
#endif

#endif
