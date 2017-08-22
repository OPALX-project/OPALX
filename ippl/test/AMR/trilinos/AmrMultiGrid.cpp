#include "AmrMultiGrid.h"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Array.H>

AmrMultiGrid::AmrMultiGrid() : epetra_comm_m(Ippl::getComm())//,
//                                solver_mp(new TrilinosSolver())
{ }

void AmrMultiGrid::solve(const amrex::Array<AmrField_u>& rho,
                         amrex::Array<AmrField_u>& phi,
                         amrex::Array<AmrField_u>& efield,
                         const amrex::Array<AmrGeometry_t>& geom,
                         int lbase, int lfine)
{
    lbase_m = lbase;
    lfine_m = lfine;
    int nLevel = lfine - lbase + 1;
    
    // initialize all levels
    mglevel_m.resize(nLevel);
    for (int lev = 0; lev < nLevel; ++lev) {
        int ilev = lbase + lev;
        mglevel_m[lev].reset(new AmrMultiGridLevel_t(rho[ilev], phi[ilev],
                                                     geom[ilev], epetra_comm_m));
    }
    
    double err = 1.0e7;
    double eps = 1.0e-8;
    
    rr_m = AmrIntVect_t(2, 2, 2);
    
    
    double rhonorm = 0.0;
    for (int i = 0; i < nLevel; ++i) {
        int ilev = lbase + i;
        rhonorm += rho[i]->norm2(0);
    }
    
    do {
        
        relax_m(lfine);
        
        err = error_m();
        
        std::cout << err << " " << eps * rhonorm << std::endl;
        
    } while ( err > eps * rhonorm);
}


void AmrMultiGrid::relax_m(int lev) {
    
    if ( lev == lfine_m ) {
        
        mglevel_m[lev]->buildRHS();
        
        // A * phi = rho --> phi = A^(-1) rho
        solver_m.solve(mglevel_m[lev]->A_p,
                         mglevel_m[lev]->phi_p,
                         mglevel_m[lev]->rho_p);
        
        mglevel_m[lev]->copyBack(*mglevel_m[lev]->phi, mglevel_m[lev]->phi_p);
        
        solver_m.residual(mglevel_m[lev]->r_p,
                          mglevel_m[lev]->phi_p,
                          mglevel_m[lev]->rho_p);
        
        mglevel_m[lev]->copyBack(*mglevel_m[lev]->residual, mglevel_m[lev]->r_p);
        
        
    } else if ( lev > lbase_m ) {
        
        // phi saved
//         mglevel_m[lev]->save();
        
        // grsb level
        mglevel_m[lev]->error->setVal(0.0, 1);
        
        this->saxpy_m(*mglevel_m[lev]->phi, *mglevel_m[lev]->error);
        
        solver_m.solve(mglevel_m[lev]->A_p,
                         mglevel_m[lev]->phi_p,
                         mglevel_m[lev]->rho_p);
        
        this->relax_m(lev - 1);
        
        
        AmrField_t tmp(mglevel_m[lev]->grids,
                       mglevel_m[lev]->dmap,
                       1, 1);
        
        tmp.setVal(0.0, 1);
        
        this->prolongate_m(tmp,
                           *mglevel_m[lev-1]->error,
                           mglevel_m[lev]->geom,
                           mglevel_m[lev-1]->geom,
                           rr_m);
        
        this->saxpy_m(*mglevel_m[lev]->error, tmp);
        
        // residual correction
        solver_m.solve(mglevel_m[lev]->A_p,
                         mglevel_m[lev]->err_p,
                         mglevel_m[lev]->r_p);
        
        this->saxpy_m(*mglevel_m[lev]->residual,
                      *mglevel_m[lev]->error,
                      -1.0);
        
        tmp.setVal(0.0, 1);
        
        // grsb level
        // ...
        
        this->saxpy_m(*mglevel_m[lev]->error, tmp);
        
        this->saxpy_m(*mglevel_m[lev]->phi_saved, *mglevel_m [lev]->error);
        
        
    } else {
        solver_m.solve(mglevel_m[lev]->A_p,
                         mglevel_m[lev]->err_p,
                         mglevel_m[lev]->r_p);
    }
}


void AmrMultiGrid::restrict_m(AmrField_t& crse,
                              const AmrField_t& fine,
                              const AmrGeometry_t& cgeom,
                              const AmrGeometry_t& fgeom,
                              const AmrIntVect_t& rr)
{
    amrex::average_down(fine, crse, fgeom, cgeom, 0, 1, rr);
}


void AmrMultiGrid::prolongate_m(AmrField_t& fine,
                                const AmrField_t& crse,
                                const AmrGeometry_t& fgeom,
                                const AmrGeometry_t& cgeom,
                                const AmrIntVect_t& rr)
{
    NoOpPhysBC cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    amrex::Array<amrex::BCRec> bcs(1, amrex::BCRec(lo_bc, hi_bc));
    amrex::PCInterp mapper;
    
    amrex::InterpFromCoarseLevel(fine, 0.0,
                                 crse,
                                 0, 0, 1,
                                 cgeom, fgeom,
                                 cphysbc, fphysbc,
                                 rr, &mapper, bcs);
    
    fine.FillBoundary(fgeom.periodicity());
}


inline void AmrMultiGrid::saxpy_m(AmrField_t& y, const AmrField_t& x, double a) {
    AmrField_t::Saxpy(y, a, x, 0, 0, 1, y.nGrow());
}


double AmrMultiGrid::error_m() {
    double sum = 0.0;
    for (size_t i = 0; i < mglevel_m.size(); ++i) {
        sum += mglevel_m[i]->residual->norm2(0);
    }
    return sum;
}
