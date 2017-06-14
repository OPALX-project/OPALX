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
    
    /**
     * This solver only works with AmrBoxLib. In order the solver to work
     * the cells need to be of cubic shape. Otherwise the solver stops with
     * the error message that it did not converge (AMReX internal).
     * 
     * @param itsAmrObject_p has information about refinemen ratios, etc.
     */
    FMGPoissonSolver(AmrBoxLib* itsAmrObject_p);
    
    /**
     * Multigrid solve based on AMReX FMultiGrid solver. The relative tolerance is
     * set to 1e-14 and the absolute tolerance to 0.0.
     * 
     * @param rho right-hand side charge density on grid [C / m]
     * @param phi electrostatic potential (unknown) [V]
     * @param efield electric field [V / m]
     * @param baseLevel for solve
     * @param finestLevel for solve
     */
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
    
    
    /**
     * Print information abour tolerances.
     * @param os output stream where to write to
     */
    Inform &print(Inform &os) const;
    
    
private:
    /**
     * Does the actual solve. It calls the FMultiGrid solver of AMReX.
     * It uses an approximation order of 3 at Dirichlet boundaries.
     * 
     * @param rho charge density on grids [C / m]
     * @param phi electrostatic potential on grid [V]
     * @param grad_phi_edge gradient of the potential (values at faces)
     * @param geom geometry of the problem, i.e. physical domain boundaries
     * @param baseLevel for solve
     * @param finestLevel for solve
     * @returns the residuum norm
     */
    double solveWithF90_m(const AmrFieldContainer_pt& rho,
                          const AmrFieldContainer_pt& phi,
                          const amrex::Array< AmrFieldContainer_pt >& grad_phi_edge, 
                          const GeomContainer_t& geom,
                          int baseLevel,
                          int finestLevel);
    
    /*! Initialize the potential and electric field to zero
     * on each level and grid.
     * 
     * @param rho charge density [C/m] that is used for DistributionMap and BoxArray
     * at different levels
     * @param phi to be initialized according to DistributionMap and BoxArray of rho
     * @param efield to be intialized according to DistributionMap and BoxArray of rho
     * @param baseLevel to initialize
     * @param finestLevel to initialize
     */
    void initGrids_m(const AmrFieldContainer_t& rho,
                     AmrFieldContainer_t& phi,
                     AmrFieldContainer_t& efield,
                     int baseLevel,
                     int finestLevel);
    
private:
    int bc_m[2*BL_SPACEDIM];        ///< Boundary conditions
    double reltol_m;                ///< Relative tolearance for solver
    double abstol_m;                ///< Absolute tolerance for solver
};


inline Inform &operator<<(Inform &os, const FMGPoissonSolver &fs) {
    return fs.print(os);
}

#endif