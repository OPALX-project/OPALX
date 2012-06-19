#ifndef TREE_POISSON_SOLVER_H_
#define TREE_POISSON_SOLVER_H_


//////////////////////////////////////////////////////////////

#include "FieldSolver.hh"

// extern "C" {
// #include <treecode.h>
// }

template <class T, unsigned int Dim>
class TreePoissonSolver : public FieldSolver<T,Dim>
{
public:
    T theta_m, eps_m;
    int rsize_m;
    bool usequad_m;

    // constructor and destructor
    TreePoissonSolver(ChargedParticles<T,Dim> *beam, T t, T e, int r, bool u);

    virtual ~TreePoissonSolver();

    virtual void computeElectricField(ChargedParticles<T,Dim> *beam, Vektor<T,Dim> hr);

};

// needed if we're not using implicit inclusion
#include "TreePoissonSolver.cc"

#endif
