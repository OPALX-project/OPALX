#include "AmrMultiGrid.h"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Array.H>


// #include <AMReX_BCUtil.H>
#include <AMReX_MacBndry.H>

AmrMultiGrid::AmrMultiGrid() : epetra_comm_m(Ippl::getComm())//,
//                                solver_mp(new TrilinosSolver())
{ }

void AmrMultiGrid::solve(const amrex::Array<AmrField_u>& rho,
                         amrex::Array<AmrField_u>& phi,
                         amrex::Array<AmrField_u>& efield,
                         const amrex::Array<AmrGeometry_t>& geom,
                         int lbase, int lfine, bool previous)
{
    
    std::cout << "real rho: " << rho[0]->min(0) << " " << rho[0]->max(0) << std::endl;
    nIter_m = 0;
    lbase_m = lbase;
    lfine_m = lfine;
    int nLevel = lfine - lbase + 1;
    
    // initialize all levels
    mglevel_m.resize(nLevel);
    for (int lev = 0; lev < nLevel; ++lev) {
        int ilev = lbase + lev;
        mglevel_m[lev].reset(new AmrMultiGridLevel_t(rho[ilev], phi[ilev],
                                                     geom[ilev], epetra_comm_m));
        
//         mglevel_m[lev]->phi->FillBoundary(mglevel_m[lev]->geom.periodicity());
        mglevel_m[lev]->rho->FillBoundary();
    }
    
    
    
    this->setBoundaryValue_m(mglevel_m[0]->phi, mglevel_m[0]->geom);
    mglevel_m[lfine_m]->error->FillBoundary();
    
    double err = 1.0e7;
    double eps = 1.0e-8;
    
    rr_m = AmrIntVect_t(2, 2, 2);
    
    if ( !previous )
        initGuess_m();
    
    
    this->zero_m(*mglevel_m[lfine_m]->error);
    
    
    double rsum = 0.0;
    double rhosum = 0.0;
    for (int lev = 0; lev < nLevel; ++lev) {
        
        this->residual_m(*mglevel_m[lev]->residual,
                         *mglevel_m[lev]->phi,
                         *mglevel_m[lev]->rho,
                         mglevel_m[lev]->geom.CellSize());
        
        rsum += this->l2norm_m(*mglevel_m[lev]->residual);
        
        rhosum += this->l2norm_m(*mglevel_m[lev]->rho);
    }
    
    while ( rsum > eps * rhosum ) {
        
        std::cout << rsum << " " << eps * rhosum << std::endl; //std::cin.get();
        
        relax_m(lfine);
        
        rsum = error_m();
    }
    
    for (int lev = nLevel-1; lev > -1; --lev) {
        int ilev = lbase + lev;
        mglevel_m[lev]->phi->FillBoundary(mglevel_m[lev]->geom.periodicity());
        this->gradient_m(lev, *efield[ilev]);
    }
    
    std::cout << rsum << " " << eps * rhosum << std::endl;
}


void AmrMultiGrid::relax_m(int lev) {
    
    // residual is already computed at beginning
    if ( lev == lfine_m && nIter_m > 0) {
        
//         std::cout << "Finest = " << lev << std::endl;
        
        this->residual_m(*mglevel_m[lev]->residual,
                         *mglevel_m[lev]->phi,
                         *mglevel_m[lev]->rho,
                         mglevel_m[lev]->geom.CellSize());
    }
    
    if ( lev > lbase_m ) {
//         std::cout << "Level = " << lev << std::endl;
        
        this->copy_m(*mglevel_m[lev]->phi_saved,
                     *mglevel_m[lev]->phi);
        mglevel_m[lev]->phi_saved->FillBoundary();
        
        this->zero_m(*mglevel_m[lev-1]->error);
        mglevel_m[lev-1]->error->FillBoundary();
        
        // relax
//         for (int i = 0; i < 8; ++i)
            this->smooth_m(*mglevel_m[lev]->error, lev);
        
//         std::cout << "max(error) = " << mglevel_m[lev]->error->max(0) << std::endl;
        
        this->saxpy_m(*mglevel_m[lev]->phi, *mglevel_m[lev]->error);
        
//         std::cout << "max(phi) = " << mglevel_m[lev]->phi->max(0) << std::endl;
        
        this->restrict_m(lev - 1);
        
//         std::cout << "max(residual) = " << mglevel_m[lev-1]->residual->max(0) << std::endl;
        
        this->relax_m(lev - 1);
        
        this->prolongate_m(lev);
        
        
        AmrField_t derr(mglevel_m[lev]->grids,
                        mglevel_m[lev]->dmap,
                        1, 1);
        
        this->zero_m(derr);
        derr.FillBoundary();
        
//         for (int i = 0; i < 8; ++i)
            this->smooth_m(derr, lev);
        
        this->saxpy_m(*mglevel_m[lev]->error, derr);
        
        this->assign_m(*mglevel_m[lev]->phi,
                       *mglevel_m[lev]->phi_saved,
                       *mglevel_m[lev]->error);
        
        
    } else {
//         std::cout << "Level = " << lev << std::endl;
        mglevel_m[lev]->copyTo(mglevel_m[lev]->err_p, *mglevel_m[lev]->error);
        
        mglevel_m[lev]->residual->mult(-1.0, 0, 1);
        
        mglevel_m[lev]->copyTo(mglevel_m[lev]->r_p, *mglevel_m[lev]->residual);
        
        mglevel_m[lev]->residual->mult(-1.0, 0, 1);
        
        solver_m.solve(mglevel_m[lev]->A_p,
                         mglevel_m[lev]->err_p,
                         mglevel_m[lev]->r_p);
        
        mglevel_m[lev]->copyBack(*mglevel_m[lev]->error, mglevel_m[lev]->err_p);
        
        mglevel_m[lev]->error->FillBoundary();
        mglevel_m[lev]->phi->FillBoundary();
        
        this->saxpy_m(*mglevel_m[lev]->phi, *mglevel_m[lev]->error);
        
        
        ++nIter_m;
    }
}


void AmrMultiGrid::restrict_m(int level)
{
    AmrField_t& cResidual = *mglevel_m[level]->residual;
    AmrField_t& cError    = *mglevel_m[level]->error;
    AmrField_t& fResidual = *mglevel_m[level+1]->residual;
    AmrField_t& fError    = *mglevel_m[level+1]->error;
    
    const AmrGeometry_t& cgeom  = mglevel_m[level]->geom;
    const AmrGeometry_t& fgeom  = mglevel_m[level+1]->geom;
    
    // get error at higher level
    AmrField_t etmp(mglevel_m[level+1]->grids,
                    mglevel_m[level+1]->dmap,
                    fError.nComp(),
                    fError.nGrow());
    this->zero_m(etmp);
    
    this->interpolate_m(etmp, cError,
                        fgeom,
                        cgeom,
                        rr_m);
    etmp.FillBoundary();
    
    
    fError.FillBoundary();
    AmrField_t rtmp(mglevel_m[level+1]->grids,
                    mglevel_m[level+1]->dmap,
                    fError.nComp(),
                    fError.nGrow());
    this->zero_m(rtmp);
    rtmp.FillBoundary();
//     mglevel_m[level+1]->error->FillBoundary();
    fResidual.FillBoundary();
    
    this->residual_m(rtmp, /**mglevel_m[level+1]->error*/etmp, fResidual, fgeom.CellSize());
    
    amrex::average_down(rtmp, cResidual, fgeom, cgeom, 0, 1, rr_m);
    cResidual.FillBoundary();
    
    this->residual_m(cResidual,
                     *mglevel_m[level]->phi,
                     *mglevel_m[level]->rho,
                     cgeom.CellSize());
}


void AmrMultiGrid::prolongate_m(int level)
{
    AmrField_t& cError    = *mglevel_m[level-1]->error;
    AmrField_t& fError    = *mglevel_m[level]->error;
    
    const AmrGeometry_t& cgeom  = mglevel_m[level-1]->geom;
    const AmrGeometry_t& fgeom  = mglevel_m[level]->geom;
    
    // get error at higher level
    AmrField_t etmp(mglevel_m[level]->grids,
                    mglevel_m[level]->dmap,
                    fError.nComp(),
                    fError.nGrow());
    this->zero_m(etmp);
    etmp.FillBoundary();
    cError.FillBoundary();
    
    this->interpolate_m(etmp, cError,
                        fgeom,
                        cgeom,
                        rr_m);
    
    this->saxpy_m(fError, etmp);
    
    AmrField_t rtmp(mglevel_m[level]->grids,
                    mglevel_m[level]->dmap,
                    fError.nComp(),
                    fError.nGrow());
    
    this->copy_m(rtmp, *mglevel_m[level]->residual);
    rtmp.FillBoundary();
    
    this->residual_m(*mglevel_m[level]->residual,
                     fError,
                     rtmp,
                     fgeom.CellSize());
}

void AmrMultiGrid::interpolate_m(AmrField_t& fine,
                                 /*const */AmrField_t& crse,
                                 const AmrGeometry_t& fgeom,
                                 const AmrGeometry_t& cgeom,
                                 const AmrIntVect_t& rr)
{
    NoOpPhysBC cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    amrex::Array<amrex::BCRec> bcs(1, amrex::BCRec(lo_bc, hi_bc));
    amrex::PCInterp mapper;
    
    crse.FillBoundary(cgeom.periodicity());
    
    amrex::InterpFromCoarseLevel(fine, 0.0,
                                 crse,
                                 0, 0, 1,
                                 cgeom, fgeom,
                                 cphysbc, fphysbc,
                                 rr, &mapper, bcs);
    
    this->setBoundaryValue_m(&fine, fgeom, &crse, rr_m);
    
    
    fine.FillBoundary(fgeom.periodicity());
}


inline void AmrMultiGrid::saxpy_m(AmrField_t& y, const AmrField_t& x, double a) {
    AmrField_t::Saxpy(y, a, x, 0, 0, 1, x.nGrow());
    y.FillBoundary();
}

inline void AmrMultiGrid::assign_m(AmrField_t& z, const AmrField_t& x,
                                   const AmrField_t& y, double a, double b)
{
    AmrField_t::LinComb(z, a, x, 0, b, y, 0, 0, 1, z.nGrow());
    z.FillBoundary();
}


double AmrMultiGrid::error_m() {
    double sum = 0.0;
    for (size_t i = 0; i < mglevel_m.size(); ++i) {
        sum += l2norm_m(*mglevel_m[i]->residual);
    }
    return sum;
}


void AmrMultiGrid::applyLapNoFine_m(AmrField_u& residual,
                                    const AmrField_u& rhs,
                                    AmrField_u& flhs,
                                    const AmrField_u& clhs)
{
//     // intpolate boundary
//     
//     
//     const double* idx2; // =  ... // geom.CellSize();
//     
//     for (amrex::MFIter mfi(*residual, false); mfi.isValid(); ++mfi) {
//         const amrex::Box&       bx  = mfi.validbox();
//         amrex::FArrayBox&       res = (*residual)[mfi];
//         const amrex::FArrayBox& rfab = (*rhs)[mfi];
//         const amrex::FArrayBox& lfab = (*flhs)[mfi];
//         
//         
//         const int* lo = bx.loVect();
//         const int* hi = bx.hiVect();
//         
//         for (int i = lo[0]; i <= hi[0]; ++i) {
//             for (int j = lo[1]; j <= hi[1]; ++j) {
//                 for (int k = lo[2]; k <= hi[2]; ++k) {
//                     
//                     AmrIntVect_t c(i  , j  , k);
//                     AmrIntVect_t l(i-1, j  , k);
//                     AmrIntVect_t r(i+1, j  , k);
//                     AmrIntVect_t u(i  , j+1, k);
//                     AmrIntVect_t d(i  , j-1, k);
//                     AmrIntVect_t f(i  , j  , k-1);
//                     AmrIntVect_t b(i  , j  , k+1);
//                     
//                     res(c) = rfab(c) - (phi(l) - 2.0 * phi(c) + phi(r)) * idx2[0]
//                                      - (phi(u) - 2.0 * phi(c) + phi(d)) * idx2[1]
//                                      - (phi(f) - 2.0 * phi(c) + phi(b)) * idx2[2]
//                 }
//             }
//         }
//     }
}


void AmrMultiGrid::smooth_m(AmrField_t& mf, int level)
{
    AmrField_t& residual = *mglevel_m[level]->residual;
    const double* dx     = mglevel_m[level]->geom.CellSize();
    const AmrMultiGridLevel_t::mask_t& mask = *mglevel_m[level]->masks_m;
    int interior = AmrMultiGridLevel_t::Mask::INTERIOR;
    
    const double idx2[3] = {
        1.0 / (dx[0] * dx[0]),
        1.0 / (dx[1] * dx[1]),
        1.0 / (dx[2] * dx[2]),
    };
    
    double h = std::min(dx[0], dx[1]);
    h = std::min(h, dx[2]);
    
    /* boundary cells: 3 * h^2 / 16
     * interior cells: 0.25 * h^2
     */
    double lambda[2] = { 3.0 * h * h / 16.0, 0.25 * h * h};
    
//     AmrField_t etmp(mglevel_m[level]->grids,
//                     mglevel_m[level]->dmap,
//                     mf.nComp(),
//                     mf.nGrow());
//     this->copy_m(etmp, mf);
    
    for (amrex::MFIter mfi(mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&       bx   = mfi.validbox();
        amrex::FArrayBox&       fab = mf[mfi];
        const amrex::FArrayBox& rfab = residual[mfi];
        const amrex::BaseFab<int>& mfab = mask[mfi];
//         const amrex::FArrayBox& efab = etmp[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    AmrIntVect_t iv(i, j, k);
                    
                    fab(iv) += lambda[ mfab(iv) > interior ] *
                               ( laplacian_m(fab, i, j, k, &idx2[0]) - rfab(iv) );
                               
//                     std::cout << iv << " " << fab(iv) << " " << rfab(iv)
//                               << " " << lambda[ mfab(iv) > interior ] << std::endl; std::cin.get();
                }
            }
        }
    }
    
    mf.FillBoundary();
}

inline void AmrMultiGrid::zero_m(AmrField_t& mf) {
    mf.setVal(0.0, mf.nGrow());
}

inline void AmrMultiGrid::copy_m(AmrField_t& lhs, const AmrField_t& rhs) {
    AmrField_t::Copy(lhs, rhs, 0, 0, lhs.nComp(), lhs.nGrow());
}

void AmrMultiGrid::residual_m(AmrField_t& r, const AmrField_t& x,
                              const AmrField_t& b, const double* dx)
{
    const double idx2[3] = {
        1.0 / (dx[0] * dx[0]),
        1.0 / (dx[1] * dx[1]),
        1.0 / (dx[2] * dx[2]),
    };
    
    for (amrex::MFIter mfi(r, false); mfi.isValid(); ++mfi) {
        const amrex::Box&       bx  = mfi.validbox();
        amrex::FArrayBox&       rfab = r[mfi];
        const amrex::FArrayBox& xfab = x[mfi];
        const amrex::FArrayBox& bfab = b[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    AmrIntVect_t iv(i, j, k);
                    rfab(iv) = bfab(iv) - laplacian_m(xfab, i, j, k, &idx2[0]);
                }
            }
        }
    }
    
    r.FillBoundary();
}


double AmrMultiGrid::laplacian_m(const amrex::FArrayBox& fab,
                                const int& i, const int& j,
                                const int& k, const double* idx2)
{
    AmrIntVect_t center(i  , j  , k);
    AmrIntVect_t left  (i-1, j  , k);
    AmrIntVect_t right (i+1, j  , k);
    AmrIntVect_t up    (i  , j+1, k);
    AmrIntVect_t down  (i  , j-1, k);
    AmrIntVect_t front (i  , j  , k-1);
    AmrIntVect_t back  (i  , j  , k+1);
    
    return ( fab(left)  - 2.0 * fab(center) + fab(right) ) * idx2[0] +
           ( fab(up)    - 2.0 * fab(center) + fab(down)  ) * idx2[1] +
           ( fab(front) - 2.0 * fab(center) + fab(back)  ) * idx2[2];
}


inline double AmrMultiGrid::l2norm_m(const AmrField_t& x) {
    return x.norm2(0);
}


void AmrMultiGrid::initGuess_m() {
    
    // solve coarsest level + interpolate to higher levels
    mglevel_m[0]->copyTo(mglevel_m[0]->rho_p, *mglevel_m[0]->rho);
    mglevel_m[0]->copyTo(mglevel_m[0]->phi_p, *mglevel_m[0]->phi);
        
    // A * phi = rho --> phi = A^(-1) rho
    solver_m.solve(mglevel_m[0]->A_p,
                   mglevel_m[0]->phi_p,
                   mglevel_m[0]->rho_p);
        
    mglevel_m[0]->copyBack(*mglevel_m[0]->phi, mglevel_m[0]->phi_p);
    
    for (uint i = 1; i < mglevel_m.size(); ++i) {
        this->interpolate_m(*mglevel_m[i]->phi,
                            *mglevel_m[i-1]->phi,
                            mglevel_m[i]->geom,
                            mglevel_m[i-1]->geom,
                            rr_m);
    }
}

void AmrMultiGrid::gradient_m(int level, AmrField_t& efield) {
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    for (amrex::MFIter mfi(*mglevel_m[level]->phi, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx   = mfi.validbox();
        const amrex::FArrayBox&    pfab = (*mglevel_m[level]->phi)[mfi];
        amrex::FArrayBox&          efab = (efield)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
            
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    AmrIntVect_t iv(i, j, k);
                    
                    // x direction
                    AmrIntVect_t left(i-1, j, k);
                    AmrIntVect_t right(i+1, j, k);
                    efab(iv, 0) = 0.5 * (pfab(left) - pfab(right)) / dx[0];
                    
                    // y direction
                    AmrIntVect_t up(i, j+1, k);
                    AmrIntVect_t down(i, j-1, k);
                    efab(iv, 1) = 0.5 * (pfab(down) - pfab(up)) / dx[1];
                    
                    // z direction
                    AmrIntVect_t back(i, j, k+1);
                    AmrIntVect_t front(i, j, k-1);
                    efab(iv, 2) = 0.5 * (pfab(front) - pfab(back)) / dx[2];
                }
            }
        }
    }
    
}


void AmrMultiGrid::setBoundaryValue_m(AmrField_t* phi, const AmrGeometry_t& geom, const AmrField_t* crse_phi,
                                      AmrIntVect_t crse_ratio)
{
    
    std::unique_ptr<amrex::MacBndry> bndry;
    bndry.reset(new amrex::MacBndry(phi->boxArray(), phi->DistributionMap(), 1 /*ncomp*/, geom));
    
    
    for (amrex::OrientationIter fi; fi; ++fi)
    {
        const amrex::Orientation face  = fi();
        amrex::FabSet& fs = (*bndry.get())[face];
        
        for (amrex::FabSetIter fsi(fs); fsi.isValid(); ++fsi)
        {
            fs[fsi].setVal(0, 0);
        }
    }
    
    
    // The values of phys_bc & ref_ratio do not matter
    // because we are not going to use those parts of MacBndry.
    AmrIntVect_t ref_ratio = AmrIntVect_t::TheZeroVector();
    amrex::Array<int> lo_bc(BL_SPACEDIM, 0);
    amrex::Array<int> hi_bc(BL_SPACEDIM, 0);
    amrex::BCRec phys_bc(lo_bc.dataPtr(), hi_bc.dataPtr());
    
//         if (crse_phi == 0 && phi == 0) 
//         {
//             bndry.setHomogValues(phys_bc, ref_ratio);
//         }
    if (crse_phi == 0 && phi != 0)
    {
        bndry->setBndryValues(*phi, 0, 0, phi->nComp(), phys_bc); 
    }
    else if (crse_phi != 0 && phi != 0) 
    {
        const int ncomp      = phi->nComp();
        const int in_rad     = 0;
        const int out_rad    = 1;
        const int extent_rad = 2;

        amrex::BoxArray crse_boxes(phi->boxArray());
        crse_boxes.coarsen(crse_ratio);
    
        amrex::BndryRegister crse_br(crse_boxes, phi->DistributionMap(),
                                    in_rad, out_rad, extent_rad, ncomp);
        crse_br.copyFrom(*crse_phi, crse_phi->nGrow(), 0, 0, ncomp);
            
        bndry->setBndryValues(crse_br, 0, *phi, 0, 0, ncomp, crse_ratio, phys_bc, 3);
        
        
    }
    else
    {
        amrex::Abort("FMultiGrid::Boundary::build_bndry: How did we get here?");
    }
    
    for (amrex::OrientationIter oitr; oitr; ++oitr)
    {
        const amrex::Orientation face  = oitr();
        const amrex::FabSet& fs = bndry->bndryValues(oitr());
        for (amrex::MFIter umfi(*(phi)); umfi.isValid(); ++umfi)
        {
            amrex::FArrayBox& dest = (*phi)[umfi];
            dest.copy(fs[umfi],fs[umfi].box());
        }
    }
}
