#ifndef P2P_POISSON_SOLVER_H_
#define P2P_POISSON_SOLVER_H_


//////////////////////////////////////////////////////////////

#include "FieldSolver.hh"
//#include "AppTypes/dcomplex.h"
//#include "FFT/FFT.h"
//////////////////////////////////////////////////////////////

template <class T, unsigned int Dim>
class P2PPoissonSolver : public FieldSolver<T,Dim>
{
public:

    T eps_m;
    // constructor and destructor
    P2PPoissonSolver(T e);
    virtual ~P2PPoissonSolver();

    virtual void computeElectricField(ChargedParticles<T,Dim> *beam, Vektor<T,Dim> hr);
};

// needed if we're not using implicit inclusion
#include "P2PPoissonSolver.cc"

#endif
