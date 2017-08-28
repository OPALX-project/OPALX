#ifndef AMR_MULTI_GRID_H
#define AMR_MULTI_GRID_H

#include <vector>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>


#include "AmrMultiGridLevel.h"

#include "AmrTrilinearInterpolater.h"


#include "TrilinosSolver.h"


// Helper class
class NoOpPhysBC
    : public amrex::PhysBCFunctBase
{
public:
    NoOpPhysBC () {}
    virtual ~NoOpPhysBC () {}
    virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override { }
    using amrex::PhysBCFunctBase::FillBoundary;
};


class AmrMultiGrid {
    
public:
    typedef TrilinosSolver::matrix_t matrix_t;
    typedef TrilinosSolver::vector_t vector_t;
    
    typedef AmrMultiGridLevel<matrix_t, vector_t> AmrMultiGridLevel_t;
    
    typedef AmrMultiGridLevel_t::AmrField_t AmrField_t;
    typedef AmrMultiGridLevel_t::AmrGeometry_t AmrGeometry_t;
    typedef AmrMultiGridLevel_t::AmrField_u AmrField_u;
    typedef AmrMultiGridLevel_t::AmrField_s AmrField_s;
    typedef AmrMultiGridLevel_t::AmrIntVect_t AmrIntVect_t;
    typedef AmrMultiGridLevel_t::indices_t indices_t;
    typedef AmrMultiGridLevel_t::coefficients_t coefficients_t;
    
    enum Interpolater {
        TRILINEAR = 0 //,
//         TRICUBIC  = 1,
    };
    
public:
    
    AmrMultiGrid(Interpolater interp = Interpolater::TRILINEAR);
    
    void solve(const amrex::Array<AmrField_u>& rho,
               amrex::Array<AmrField_u>& phi,
               amrex::Array<AmrField_u>& efield,
               const amrex::Array<AmrGeometry_t>& geom,
               int lbase, int lfine, bool previous = false);
    
    
private:
    
//     void relax_m(int lev);
//     
//     void restrict_m(int level);
//     
//     void interpolate_m(AmrField_t& fine, /*const */AmrField_t& crse,
//                       const AmrGeometry_t& fgeom,
//                       const AmrGeometry_t& cgeom, const AmrIntVect_t& rr);
//     
//     void prolongate_m(int level);
//     
//     
//     
//     /*!
//      * z = a * x + b * y
//      */
//     void assign_m(AmrField_t& z, const AmrField_t& x,
//                   const AmrField_t& y, double a = 1.0, double b = 1.0);
//     
//     double error_m();
//     
//     void applyLapNoFine_m(AmrField_u& residual,
//                           const AmrField_u& rhs,
//                           AmrField_u& flhs,
//                           const AmrField_u& clhs);
//     
//     /////////
//     
//     // lhs = rhs (only first component is copied)
//     void copy_m(AmrField_t& lhs, const AmrField_t& rhs);
//     void zero_m(AmrField_t& mf);
//     void smooth_m(AmrField_t& mf, int level);
//     /*!
//      * y = a * x + y
//      * @param y
//      * @param x
//      * @param a
//      */
//     void saxpy_m(AmrField_t& y, const AmrField_t& x, double a = 1.0);
//     
// //     void average_m();
// //     void interpolate_m();
// //     void add_m();
//     
//     /*
//      * r = b - A*x
//      */
//     void residual_m(AmrField_t& r, const AmrField_t& x, const AmrField_t& b, const double* dx);
//     
//     double laplacian_m(const amrex::FArrayBox& fab, const int& i, const int& j, const int& k, const double* idx2);
//     
//     double l2norm_m(const AmrField_t& x);
//     
//     void initGuess_m();
//     
//     
//     void gradient_m(int level, AmrField_t& efield);
//     
//     
//     void setBoundaryValue_m(AmrField_t* phi, const AmrGeometry_t& geom,
//                             const AmrField_t* crse_phi = 0,
//                             AmrIntVect_t crse_ratio = AmrIntVect_t::TheZeroVector());
    
    
    void buildPoissonMatrix_m(int level);
    
    void buildRestrictionMatrix_m(int level);
    
    void buildInterpolationMatrix_m(int level);
    
    void buildBoundaryMatrix_m(int level);
    
    void buildDensityVector_m(const AmrField_t& rho, int level);
    
    void buildPotentialVector_m(const AmrField_t& phi, int level);
    
    void amrex2trilinos_m(const AmrField_t& mf, Teuchos::RCP<vector_t>& mv, int level);
    
    void trilinos2amrex_m(AmrField_t& mf, const Teuchos::RCP<vector_t>& mv);
    
    
private:
    Epetra_MpiComm epetra_comm_m;
    
    std::unique_ptr<AmrInterpolater<AmrMultiGridLevel_t> > interp_mp;
    
    int nIter_m;
    
    
    std::vector<std::unique_ptr<AmrMultiGridLevel_t > > mglevel_m;
    
//     std::unique_ptr<LinearSolver<Teuchos::RCP<matrix_t>, Teuchos::RCP<vector_t> > > solver_mp;
    TrilinosSolver solver_m;
    
    AmrIntVect_t rr_m;      //TODO move to level
    
    int lbase_m;
    int lfine_m;
    
};

#endif
