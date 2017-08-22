#ifndef AMR_MULTI_GRID_H
#define AMR_MULTI_GRID_H

#include <vector>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>


#include "AmrMultiGridLevel.h"
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
    
    
public:
    
    AmrMultiGrid();
    
    void solve(const amrex::Array<AmrField_u>& rho,
               amrex::Array<AmrField_u>& phi,
               amrex::Array<AmrField_u>& efield,
               const amrex::Array<AmrGeometry_t>& geom,
               int lbase, int lfine);
    
    
private:
    
    void relax_m(int lev);
    
    void restrict_m(AmrField_t& crse, const AmrField_t& fine,
                    const AmrGeometry_t& cgeom,
                    const AmrGeometry_t& fgeom, const AmrIntVect_t& rr);
    
    void prolongate_m(AmrField_t& fine, const AmrField_t& crse,
                      const AmrGeometry_t& fgeom,
                      const AmrGeometry_t& cgeom, const AmrIntVect_t& rr);
    
    /*!
     * y = a * x + y
     * @param y
     * @param x
     * @param a
     */
    void saxpy_m(AmrField_t& y, const AmrField_t& x, double a = 1.0);
    
    double error_m();
    
    
private:
    Epetra_MpiComm epetra_comm_m;
    
    
    std::vector<std::unique_ptr<AmrMultiGridLevel_t > > mglevel_m;
    
//     std::unique_ptr<LinearSolver<Teuchos::RCP<matrix_t>, Teuchos::RCP<vector_t> > > solver_mp;
    TrilinosSolver solver_m;
    
    AmrIntVect_t rr_m;
    
    int lbase_m;
    int lfine_m;
    
};

#endif
