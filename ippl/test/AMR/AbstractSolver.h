#ifndef ABSTRACT_SOLVER_H
#define ABSTRACT_SOLVER_H

#include "AmrOpal.h"

#include "Amr/AmrDefs.h"

class AbstractSolver
{
public:
    typedef std::unique_ptr<AmrOpal> amropal_p;
    typedef amr::AmrFieldContainer_t AmrFieldContainer_t;
    
    virtual ~AbstractSolver() {};
    
    virtual
    void solve(amropal_p& amropal,
               AmrFieldContainer_t &rho,
               AmrFieldContainer_t &phi,
               AmrFieldContainer_t &efield,
               unsigned short baseLevel,
               unsigned short finestLevel,
               bool prevAsGuess = true) = 0;
    
};

#endif
