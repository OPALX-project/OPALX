#ifndef AMR_POISSON_SOLVER_H_
#define AMR_POISSON_SOLVER_H_

#include "PoissonSolver.h"

#include <memory>

template<class AmrObject>
class AmrPoissonSolver : public PoissonSolver {
    
public:
    /*!
     * @param itsAmrObject_p holds information about grids and domain
     */
    AmrPoissonSolver(AmrObject* itsAmrObject_p)
        : itsAmrObject_mp(itsAmrObject_p), regrid_m(false) {}
    
    virtual ~AmrPoissonSolver() {}
    
    
    void computePotential(Field_t &rho, Vector_t hr) {
        throw OpalException("AmrPoissonSolver::computePotential(Field_t, Vector_t)",
                            "Not implemented.");
    }
    
    void computePotential(Field_t &rho, Vector_t hr, double zshift) {
        throw OpalException("AmrPoissonSolver::computePotential(Field_t, Vector_t, double)",
                            "Not implemented.");
    }
    
    void test(PartBunchBase<double, 3> *bunch) {
        throw OpalException("AmrPoissonSolver::test(PartBunchBase<double, 3>)", "Not implemented.");
    }
    
    void hasToRegrid() {
        regrid_m = true;
    }
    
protected:
    AmrObject* itsAmrObject_mp;
    
    /// is set to true by itsAmrObject_mp and reset to false by solver
    bool regrid_m;
};


#endif