#ifndef AMR_POISSON_SOLVER_H_
#define AMR_POISSON_SOLVER_H_

#include "PoissonSolver.h"

#include <memory>

template<class AmrObject>
class AmrPoissonSolver : public PoissonSolver {
    
public:
    /*!
     * @param amrobject_p holds information about grids and domain
     */
    AmrPoissonSolver(AmrObject* amrobject_p) : amrobject_mp(amrobject_p) {}
    
    virtual ~AmrPoissonSolver();
    
    void test(PartBunch &bunch) { };
    
    
private:
    std::unique_ptr<AmrObject> amrobject_mp;    
};


#endif