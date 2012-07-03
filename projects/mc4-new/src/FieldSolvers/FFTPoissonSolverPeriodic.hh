#ifndef FFT_POISSON_SOLVER_PERIODIC_H_
#define FFT_POISSON_SOLVER_PERIODIC_H_

//////////////////////////////////////////////////////////////
#include "FieldSolver.hh"
#include "AppTypes/dcomplex.h"
#include "FFT/FFT.h"
//////////////////////////////////////////////////////////////

template <class T, unsigned int Dim>
class FFTPoissonSolverPeriodic : public FieldSolver<T,Dim>
{
public:
    // some useful typedefs
    typedef typename ChargedParticles<T,Dim>::Center_t      Center_t;
    typedef typename ChargedParticles<T,Dim>::Mesh_t        Mesh_t;
    typedef typename ChargedParticles<T,Dim>::FieldLayout_t FieldLayout_t;
    typedef typename ChargedParticles<T,Dim>::Vector_t      Vector_t;
    typedef typename ChargedParticles<T,Dim>::IntrplCIC_t   IntrplCIC_t;

    typedef Field<dcomplex, Dim, Mesh_t, Center_t>          CxField_t;
    typedef Field<T, Dim, Mesh_t, Center_t>                 RxField_t;
    typedef Field<Vector_t, Dim, Mesh_t, Center_t>          VField_t;
    typedef FFT<CCTransform, Dim, T>                        FFT_t;

    // constructor and destructor
    FFTPoissonSolverPeriodic(ChargedParticles<T,Dim> *univ, SimData<T,Dim> simData);
    FFTPoissonSolverPeriodic(Mesh_t *mesh, FieldLayout_t *FL, SimData<T,Dim> simData);
    virtual ~FFTPoissonSolverPeriodic();

    /// computes the force field from density field rho
    virtual void computeForceField(ChargedParticles<T,Dim> *univ);

    /// CIC forward move (deposits mass on rho_m)
    virtual void CICforward(ChargedParticles<T,Dim> *univ);

    /// dumps power spectra in file fn
    virtual void calcPwrSpecAndSave(ChargedParticles<T,Dim> *univ, string fn) ;

private:
    void saveField(string fn, CxField_t &f, int n );
    void saveField(string fn, RxField_t &f, int n );
    void doInit();

    /// fortrans nint function
    inline T nint(T x)
    {
        return ceil(x + 0.5) - (fmod(x*0.5 + 0.25, 1.0) != 0);
    }

    FFT_t *fft_m;

    // mesh and layout objects for rho_m
    FieldLayout_t *layout_m;
    Mesh_t *mesh_m;

    // bigger mesh (length+1)
    FieldLayout_t *FLI_m;
    Mesh_t *meshI_m;

    /// global domain for the various fields
    NDIndex<Dim> gDomain_m;
    /// local domain for the various fields
    NDIndex<Dim> lDomain_m;
    /// simulation input data
    SimData<T,Dim> simData_m;

    /// global domain for the enlarged fields
    NDIndex<Dim> gDomainL_m;

    BConds<T,Dim,Mesh_t,Center_t> bc_m;
    BConds<T,Dim,Mesh_t,Center_t> zerobc_m;
    //BConds<Vector_t,Dim,Mesh_t,Center_t> vbc_m;

    e_dim_tag dcomp_m[Dim];
    Vektor<int,Dim> nr_m;

    /// vector field for force
    VField_t gf_m;

    /// Fourier transformed density field
    CxField_t rho_m;

    /// real part of rho_m in bigger box (len+1)
    RxField_t rhocic_m;

    IpplTimings::TimerRef TScatter_m;
    IpplTimings::TimerRef TSolver_m;

    /// power spectra kmax
    int kmax_m;

    /// 1D power spectra
    T *spectra1D_m;
    /// Nk power spectra
    int *Nk_m;

};

// needed if we're not using implicit inclusion
#include "FFTPoissonSolverPeriodic.cc"

#endif

