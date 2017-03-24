#ifndef FMG_POISSON_SOLVER_H_
#define FMG_POISSON_SOLVER_H_

#include "Solvers/AmrPoissonSolver.h"

#include <BoxLib.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>
#include <FMultiGrid.H>

//TODO remote AmrBoxLib
class AmrBoxLib {
    
private:
    Array<Geometry> geoms;
    Array<BoxArray> boxarrays;
public:
    Array<Geometry>& Geom() { return geoms; }
    Geometry Geom(int level) { return geoms[level]; }
    const Array<BoxArray>& boxArray() { return boxarrays; }
    
};

class FMGPoissonSolver : public AmrPoissonSolver< AmrBoxLib > {
    
private:
    typedef Array<Geometry> GeomContainer_t;
    
public:
    
    FMGPoissonSolver(AmrBoxLib* amrobject_p);
    
    void solve(AmrFieldContainer_t &rho,
               AmrFieldContainer_t &efield,
               unsigned short baseLevel,
               unsigned short finestLevel);
    
    double getXRangeMin(unsigned short level = 0);
    double getXRangeMax(unsigned short level = 0);
    double getYRangeMin(unsigned short level = 0);
    double getYRangeMax(unsigned short level = 0);
    double getZRangeMin(unsigned short level = 0);
    double getZRangeMax(unsigned short level = 0);
    
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