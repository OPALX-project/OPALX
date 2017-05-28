#ifndef FMG_POISSON_SOLVER_H_
#define FMG_POISSON_SOLVER_H_

#include "Solvers/AmrPoissonSolver.h"
#include "Amr/AmrBoxLib.h"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_FMultiGrid.H>

class FMGPoissonSolver : public AmrPoissonSolver< AmrBoxLib > {
    
private:
    typedef AmrBoxLib::AmrGeomContainer_t GeomContainer_t;
    typedef amrex::Array<AmrBoxLib::AmrField_t*> AmrFieldContainer_pt;
    typedef AmrBoxLib::AmrGeometry_t AmrGeometry_t;
    typedef AmrBoxLib::AmrGrid_t AmrGrid_t;
    typedef AmrBoxLib::AmrProcMap_t           AmrProcMap_t;
    typedef AmrBoxLib::AmrFieldContainer_t    AmrFieldContainer_t;
    
public:
    
    FMGPoissonSolver(AmrBoxLib* itsAmrObject_p);
    
    void solve(AmrFieldContainer_t &rho,
               AmrFieldContainer_t& phi,
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
    void solveWithF90_m(const AmrFieldContainer_pt& rho,
                        const AmrFieldContainer_pt& phi,
                        const amrex::Array< AmrFieldContainer_pt >& grad_phi_edge, 
                        const GeomContainer_t& geom,
                        int baseLevel,
                        int finestLevel);
    
    /*! Initialize the potential and electric field to zero
     * on each level and grid.
     * @param efield to be intialized.
     */
    void initGrids_m(const AmrFieldContainer_t& rho,
                     AmrFieldContainer_t& phi,
                     AmrFieldContainer_t& efield,
                     int baseLevel,
                     int finestLevel);
    
private:
    int bc_m[2*BL_SPACEDIM];        ///< Boundary conditions
    double tol_m;
    double abstol_m;
//     AmrFieldContainer_t phi_m;
};

inline Inform &operator<<(Inform &os, const FMGPoissonSolver &fs) {
    return fs.print(os);
}

#endif