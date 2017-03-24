#ifndef FMG_POISSON_SOLVER_H_
#define FMG_POISSON_SOLVER_H_

#include <BoxLib.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>
#include <FMultiGrid.H>

class AmrBoxlib; //TODO remove template and do AmrPoissonSolver<AmrBoxlib>

class FMGPoissonSolver : public AmrPoissonSolver<AmrBoxlib> {
    
private:
    typedef Array<Geometry> GeomContainer_t;
    
public:
    
    FMGPoissonSolver(AmrBoxlib* amrobject_p);
    
    Inform &print(Inform &os) const;
    
    
private:
    void solveWithF90_m(AmrFieldContainer_t& rho,
                        Array< AmrFieldContainer_t >& grad_phi_edge, 
                        const GeomContainer_t& geom,
                        int baseLevel,
                        int finestLevel);
    
    /*! Initialize the potential and electric field to zero
     * on each level and grid.
     * @param efield to be intialized.
     */
    void initGrids_m(AmrFieldContainer_t& efield);
    
private:
    int bc_m[2*BL_SPACEDIM];        ///< Boundary conditions
    double tol_m;
    double abstol_m;
    AmrFieldContainer_t phi_m;
};

inline Inform &operator<<(Inform &os, const FMGPoissonSolver &fs) {
    return fs.print(os);
}

#endif