#include "AmrMultiGrid.h"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Array.H>

#include <AMReX_MacBndry.H>

AmrMultiGrid::AmrMultiGrid(Interpolater interp) : epetra_comm_m(Ippl::getComm())
{
    switch ( interp ) {
        case Interpolater::TRILINEAR:
            interp_mp.reset( new AmrTrilinearInterpolater<AmrMultiGridLevel_t>() );
        default:
            break;
    }
}

void AmrMultiGrid::solve(const amrex::Array<AmrField_u>& rho,
                         amrex::Array<AmrField_u>& phi,
                         amrex::Array<AmrField_u>& efield,
                         const amrex::Array<AmrGeometry_t>& geom,
                         int lbase, int lfine, bool previous)
{
    nIter_m = 0;
    lbase_m = lbase;
    lfine_m = lfine;
    int nLevel = lfine - lbase + 1;
    
    rr_m = AmrIntVect_t(2, 2, 2);
    
    // initialize all levels
    mglevel_m.resize(nLevel);
    for (int lev = 0; lev < nLevel; ++lev) {
        int ilev = lbase + lev;
        
        if ( lev == 0 ) {
            //FIXME also Open boundary
            mglevel_m[lev].reset(new AmrMultiGridLevel_t(rho[ilev]->boxArray(),
                                                         rho[ilev]->DistributionMap(),
                                                         geom[ilev],
                                                         rr_m,
                                                         new AmrDirichletBoundary<AmrMultiGridLevel_t>(), 
                                                         epetra_comm_m));
        } else {
            // all higher levels have Dirichlet BC
            mglevel_m[lev].reset(new AmrMultiGridLevel_t(rho[ilev]->boxArray(),
                                                         rho[ilev]->DistributionMap(),
                                                         geom[ilev],
                                                         rr_m,
                                                         new AmrDirichletBoundary<AmrMultiGridLevel_t>(), 
                                                         epetra_comm_m));
        }
    }
    
    
    // build all necessary matrices and vectors
    for (int lev = 0; lev < nLevel; ++lev) {
        this->buildRestrictionMatrix_m(lev);
        this->buildInterpolationMatrix_m(lev);
        this->buildBoundaryMatrix_m(lev);
        
        int ilev = lbase + lev;
        
        this->buildDensityVector_m(*rho[ilev], lev);
        
        this->buildPotentialVector_m(*phi[ilev], lev);
    }
    
    
    
//     this->setBoundaryValue_m(mglevel_m[0]->phi, mglevel_m[0]->geom);
//     mglevel_m[lfine_m]->error->FillBoundary();
//     
//     double err = 1.0e7;
//     double eps = 1.0e-8;
//     
//     
//     if ( !previous )
//         initGuess_m();
//     
//     
//     this->zero_m(*mglevel_m[lfine_m]->error);
//     
//     
//     double rsum = 0.0;
//     double rhosum = 0.0;
//     for (int lev = 0; lev < nLevel; ++lev) {
//         
//         this->residual_m(*mglevel_m[lev]->residual,
//                          *mglevel_m[lev]->phi,
//                          *mglevel_m[lev]->rho,
//                          mglevel_m[lev]->geom.CellSize());
//         
//         rsum += this->l2norm_m(*mglevel_m[lev]->residual);
//         
//         rhosum += this->l2norm_m(*mglevel_m[lev]->rho);
//     }
//     
//     while ( rsum > eps * rhosum ) {
//         
//         std::cout << rsum << " " << eps * rhosum << std::endl; //std::cin.get();
//         
//         relax_m(lfine);
//         
//         rsum = error_m();
//     }
//     
//     for (int lev = nLevel-1; lev > -1; --lev) {
//         int ilev = lbase + lev;
//         mglevel_m[lev]->phi->FillBoundary(mglevel_m[lev]->geom.periodicity());
//         this->gradient_m(lev, *efield[ilev]);
//     }
//     
//     std::cout << rsum << " " << eps * rhosum << std::endl;
}


// void AmrMultiGrid::relax_m(int lev) {
//     
//     // residual is already computed at beginning
//     if ( lev == lfine_m && nIter_m > 0) {
//         
// //         std::cout << "Finest = " << lev << std::endl;
//         
//         this->residual_m(*mglevel_m[lev]->residual,
//                          *mglevel_m[lev]->phi,
//                          *mglevel_m[lev]->rho,
//                          mglevel_m[lev]->geom.CellSize());
//     }
//     
//     if ( lev > lbase_m ) {
// //         std::cout << "Level = " << lev << std::endl;
//         
//         this->copy_m(*mglevel_m[lev]->phi_saved,
//                      *mglevel_m[lev]->phi);
//         mglevel_m[lev]->phi_saved->FillBoundary();
//         
//         this->zero_m(*mglevel_m[lev-1]->error);
//         mglevel_m[lev-1]->error->FillBoundary();
//         
//         // relax
// //         for (int i = 0; i < 8; ++i)
//             this->smooth_m(*mglevel_m[lev]->error, lev);
//         
// //         std::cout << "max(error) = " << mglevel_m[lev]->error->max(0) << std::endl;
//         
//         this->saxpy_m(*mglevel_m[lev]->phi, *mglevel_m[lev]->error);
//         
// //         std::cout << "max(phi) = " << mglevel_m[lev]->phi->max(0) << std::endl;
//         
//         this->restrict_m(lev - 1);
//         
// //         std::cout << "max(residual) = " << mglevel_m[lev-1]->residual->max(0) << std::endl;
//         
//         this->relax_m(lev - 1);
//         
//         this->prolongate_m(lev);
//         
//         
//         AmrField_t derr(mglevel_m[lev]->grids,
//                         mglevel_m[lev]->dmap,
//                         1, 1);
//         
//         this->zero_m(derr);
//         derr.FillBoundary();
//         
// //         for (int i = 0; i < 8; ++i)
//             this->smooth_m(derr, lev);
//         
//         this->saxpy_m(*mglevel_m[lev]->error, derr);
//         
//         this->assign_m(*mglevel_m[lev]->phi,
//                        *mglevel_m[lev]->phi_saved,
//                        *mglevel_m[lev]->error);
//         
//         
//     } else {
// //         std::cout << "Level = " << lev << std::endl;
//         mglevel_m[lev]->copyTo(mglevel_m[lev]->err_p, *mglevel_m[lev]->error);
//         
//         mglevel_m[lev]->residual->mult(-1.0, 0, 1);
//         
//         mglevel_m[lev]->copyTo(mglevel_m[lev]->r_p, *mglevel_m[lev]->residual);
//         
//         mglevel_m[lev]->residual->mult(-1.0, 0, 1);
//         
//         solver_m.solve(mglevel_m[lev]->A_p,
//                          mglevel_m[lev]->err_p,
//                          mglevel_m[lev]->r_p);
//         
//         mglevel_m[lev]->copyBack(*mglevel_m[lev]->error, mglevel_m[lev]->err_p);
//         
//         mglevel_m[lev]->error->FillBoundary();
//         mglevel_m[lev]->phi->FillBoundary();
//         
//         this->saxpy_m(*mglevel_m[lev]->phi, *mglevel_m[lev]->error);
//         
//         
//         ++nIter_m;
//     }
// }
// 
// 
// void AmrMultiGrid::restrict_m(int level)
// {
//     AmrField_t& cResidual = *mglevel_m[level]->residual;
//     AmrField_t& cError    = *mglevel_m[level]->error;
//     AmrField_t& fResidual = *mglevel_m[level+1]->residual;
//     AmrField_t& fError    = *mglevel_m[level+1]->error;
//     
//     const AmrGeometry_t& cgeom  = mglevel_m[level]->geom;
//     const AmrGeometry_t& fgeom  = mglevel_m[level+1]->geom;
//     
//     // get error at higher level
//     AmrField_t etmp(mglevel_m[level+1]->grids,
//                     mglevel_m[level+1]->dmap,
//                     fError.nComp(),
//                     fError.nGrow());
//     this->zero_m(etmp);
//     
//     this->interpolate_m(etmp, cError,
//                         fgeom,
//                         cgeom,
//                         rr_m);
//     etmp.FillBoundary();
//     
//     
//     fError.FillBoundary();
//     AmrField_t rtmp(mglevel_m[level+1]->grids,
//                     mglevel_m[level+1]->dmap,
//                     fError.nComp(),
//                     fError.nGrow());
//     this->zero_m(rtmp);
//     rtmp.FillBoundary();
// //     mglevel_m[level+1]->error->FillBoundary();
//     fResidual.FillBoundary();
//     
//     this->residual_m(rtmp, /**mglevel_m[level+1]->error*/etmp, fResidual, fgeom.CellSize());
//     
//     amrex::average_down(rtmp, cResidual, fgeom, cgeom, 0, 1, rr_m);
//     cResidual.FillBoundary();
//     
//     this->residual_m(cResidual,
//                      *mglevel_m[level]->phi,
//                      *mglevel_m[level]->rho,
//                      cgeom.CellSize());
// }
// 
// 
// void AmrMultiGrid::prolongate_m(int level)
// {
//     AmrField_t& cError    = *mglevel_m[level-1]->error;
//     AmrField_t& fError    = *mglevel_m[level]->error;
//     
//     const AmrGeometry_t& cgeom  = mglevel_m[level-1]->geom;
//     const AmrGeometry_t& fgeom  = mglevel_m[level]->geom;
//     
//     // get error at higher level
//     AmrField_t etmp(mglevel_m[level]->grids,
//                     mglevel_m[level]->dmap,
//                     fError.nComp(),
//                     fError.nGrow());
//     this->zero_m(etmp);
//     etmp.FillBoundary();
//     cError.FillBoundary();
//     
//     this->interpolate_m(etmp, cError,
//                         fgeom,
//                         cgeom,
//                         rr_m);
//     
//     this->saxpy_m(fError, etmp);
//     
//     AmrField_t rtmp(mglevel_m[level]->grids,
//                     mglevel_m[level]->dmap,
//                     fError.nComp(),
//                     fError.nGrow());
//     
//     this->copy_m(rtmp, *mglevel_m[level]->residual);
//     rtmp.FillBoundary();
//     
//     this->residual_m(*mglevel_m[level]->residual,
//                      fError,
//                      rtmp,
//                      fgeom.CellSize());
// }
// 
// void AmrMultiGrid::interpolate_m(AmrField_t& fine,
//                                  /*const */AmrField_t& crse,
//                                  const AmrGeometry_t& fgeom,
//                                  const AmrGeometry_t& cgeom,
//                                  const AmrIntVect_t& rr)
// {
//     NoOpPhysBC cphysbc, fphysbc;
//     int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
//     int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
//     amrex::Array<amrex::BCRec> bcs(1, amrex::BCRec(lo_bc, hi_bc));
//     amrex::PCInterp mapper;
//     
//     crse.FillBoundary(cgeom.periodicity());
//     
//     amrex::InterpFromCoarseLevel(fine, 0.0,
//                                  crse,
//                                  0, 0, 1,
//                                  cgeom, fgeom,
//                                  cphysbc, fphysbc,
//                                  rr, &mapper, bcs);
//     
//     this->setBoundaryValue_m(&fine, fgeom, &crse, rr_m);
//     
//     
//     fine.FillBoundary(fgeom.periodicity());
// }
// 
// 
// inline void AmrMultiGrid::saxpy_m(AmrField_t& y, const AmrField_t& x, double a) {
//     AmrField_t::Saxpy(y, a, x, 0, 0, 1, x.nGrow());
//     y.FillBoundary();
// }
// 
// inline void AmrMultiGrid::assign_m(AmrField_t& z, const AmrField_t& x,
//                                    const AmrField_t& y, double a, double b)
// {
//     AmrField_t::LinComb(z, a, x, 0, b, y, 0, 0, 1, z.nGrow());
//     z.FillBoundary();
// }
// 
// 
// double AmrMultiGrid::error_m() {
//     double sum = 0.0;
//     for (size_t i = 0; i < mglevel_m.size(); ++i) {
//         sum += l2norm_m(*mglevel_m[i]->residual);
//     }
//     return sum;
// }
// 
// 
// void AmrMultiGrid::applyLapNoFine_m(AmrField_u& residual,
//                                     const AmrField_u& rhs,
//                                     AmrField_u& flhs,
//                                     const AmrField_u& clhs)
// {
// //     // intpolate boundary
// //     
// //     
// //     const double* idx2; // =  ... // geom.CellSize();
// //     
// //     for (amrex::MFIter mfi(*residual, false); mfi.isValid(); ++mfi) {
// //         const amrex::Box&       bx  = mfi.validbox();
// //         amrex::FArrayBox&       res = (*residual)[mfi];
// //         const amrex::FArrayBox& rfab = (*rhs)[mfi];
// //         const amrex::FArrayBox& lfab = (*flhs)[mfi];
// //         
// //         
// //         const int* lo = bx.loVect();
// //         const int* hi = bx.hiVect();
// //         
// //         for (int i = lo[0]; i <= hi[0]; ++i) {
// //             for (int j = lo[1]; j <= hi[1]; ++j) {
// //                 for (int k = lo[2]; k <= hi[2]; ++k) {
// //                     
// //                     AmrIntVect_t c(i  , j  , k);
// //                     AmrIntVect_t l(i-1, j  , k);
// //                     AmrIntVect_t r(i+1, j  , k);
// //                     AmrIntVect_t u(i  , j+1, k);
// //                     AmrIntVect_t d(i  , j-1, k);
// //                     AmrIntVect_t f(i  , j  , k-1);
// //                     AmrIntVect_t b(i  , j  , k+1);
// //                     
// //                     res(c) = rfab(c) - (phi(l) - 2.0 * phi(c) + phi(r)) * idx2[0]
// //                                      - (phi(u) - 2.0 * phi(c) + phi(d)) * idx2[1]
// //                                      - (phi(f) - 2.0 * phi(c) + phi(b)) * idx2[2]
// //                 }
// //             }
// //         }
// //     }
// }
// 
// 
// void AmrMultiGrid::smooth_m(AmrField_t& mf, int level)
// {
//     AmrField_t& residual = *mglevel_m[level]->residual;
//     const double* dx     = mglevel_m[level]->geom.CellSize();
//     const AmrMultiGridLevel_t::mask_t& mask = *mglevel_m[level]->masks_m;
//     int interior = AmrMultiGridLevel_t::Mask::INTERIOR;
//     
//     const double idx2[3] = {
//         1.0 / (dx[0] * dx[0]),
//         1.0 / (dx[1] * dx[1]),
//         1.0 / (dx[2] * dx[2]),
//     };
//     
//     double h = std::min(dx[0], dx[1]);
//     h = std::min(h, dx[2]);
//     
//     /* boundary cells: 3 * h^2 / 16
//      * interior cells: 0.25 * h^2
//      */
//     double lambda[2] = { 3.0 * h * h / 16.0, 0.25 * h * h};
//     
// //     AmrField_t etmp(mglevel_m[level]->grids,
// //                     mglevel_m[level]->dmap,
// //                     mf.nComp(),
// //                     mf.nGrow());
// //     this->copy_m(etmp, mf);
//     
//     for (amrex::MFIter mfi(mask, false); mfi.isValid(); ++mfi) {
//         const amrex::Box&       bx   = mfi.validbox();
//         amrex::FArrayBox&       fab = mf[mfi];
//         const amrex::FArrayBox& rfab = residual[mfi];
//         const amrex::BaseFab<int>& mfab = mask[mfi];
// //         const amrex::FArrayBox& efab = etmp[mfi];
//         
//         const int* lo = bx.loVect();
//         const int* hi = bx.hiVect();
//         
//         for (int i = lo[0]; i <= hi[0]; ++i) {
//             for (int j = lo[1]; j <= hi[1]; ++j) {
//                 for (int k = lo[2]; k <= hi[2]; ++k) {
//                     AmrIntVect_t iv(i, j, k);
//                     
//                     fab(iv) += lambda[ mfab(iv) > interior ] *
//                                ( laplacian_m(fab, i, j, k, &idx2[0]) - rfab(iv) );
//                                
// //                     std::cout << iv << " " << fab(iv) << " " << rfab(iv)
// //                               << " " << lambda[ mfab(iv) > interior ] << std::endl; std::cin.get();
//                 }
//             }
//         }
//     }
//     
//     mf.FillBoundary();
// }
// 
// inline void AmrMultiGrid::zero_m(AmrField_t& mf) {
//     mf.setVal(0.0, mf.nGrow());
// }
// 
// inline void AmrMultiGrid::copy_m(AmrField_t& lhs, const AmrField_t& rhs) {
//     AmrField_t::Copy(lhs, rhs, 0, 0, lhs.nComp(), lhs.nGrow());
// }
// 
// void AmrMultiGrid::residual_m(AmrField_t& r, const AmrField_t& x,
//                               const AmrField_t& b, const double* dx)
// {
//     const double idx2[3] = {
//         1.0 / (dx[0] * dx[0]),
//         1.0 / (dx[1] * dx[1]),
//         1.0 / (dx[2] * dx[2]),
//     };
//     
//     for (amrex::MFIter mfi(r, false); mfi.isValid(); ++mfi) {
//         const amrex::Box&       bx  = mfi.validbox();
//         amrex::FArrayBox&       rfab = r[mfi];
//         const amrex::FArrayBox& xfab = x[mfi];
//         const amrex::FArrayBox& bfab = b[mfi];
//         
//         const int* lo = bx.loVect();
//         const int* hi = bx.hiVect();
//         
//         for (int i = lo[0]; i <= hi[0]; ++i) {
//             for (int j = lo[1]; j <= hi[1]; ++j) {
//                 for (int k = lo[2]; k <= hi[2]; ++k) {
//                     AmrIntVect_t iv(i, j, k);
//                     rfab(iv) = bfab(iv) - laplacian_m(xfab, i, j, k, &idx2[0]);
//                 }
//             }
//         }
//     }
//     
//     r.FillBoundary();
// }
// 
// 
// double AmrMultiGrid::laplacian_m(const amrex::FArrayBox& fab,
//                                 const int& i, const int& j,
//                                 const int& k, const double* idx2)
// {
//     AmrIntVect_t center(i  , j  , k);
//     AmrIntVect_t left  (i-1, j  , k);
//     AmrIntVect_t right (i+1, j  , k);
//     AmrIntVect_t up    (i  , j+1, k);
//     AmrIntVect_t down  (i  , j-1, k);
//     AmrIntVect_t front (i  , j  , k-1);
//     AmrIntVect_t back  (i  , j  , k+1);
//     
//     return ( fab(left)  - 2.0 * fab(center) + fab(right) ) * idx2[0] +
//            ( fab(up)    - 2.0 * fab(center) + fab(down)  ) * idx2[1] +
//            ( fab(front) - 2.0 * fab(center) + fab(back)  ) * idx2[2];
// }
// 
// 
// inline double AmrMultiGrid::l2norm_m(const AmrField_t& x) {
//     return x.norm2(0);
// }
// 
// 
// void AmrMultiGrid::initGuess_m() {
//     
//     // solve coarsest level + interpolate to higher levels
//     mglevel_m[0]->copyTo(mglevel_m[0]->rho_p, *mglevel_m[0]->rho);
//     mglevel_m[0]->copyTo(mglevel_m[0]->phi_p, *mglevel_m[0]->phi);
//         
//     // A * phi = rho --> phi = A^(-1) rho
//     solver_m.solve(mglevel_m[0]->A_p,
//                    mglevel_m[0]->phi_p,
//                    mglevel_m[0]->rho_p);
//         
//     mglevel_m[0]->copyBack(*mglevel_m[0]->phi, mglevel_m[0]->phi_p);
//     
//     for (uint i = 1; i < mglevel_m.size(); ++i) {
//         this->interpolate_m(*mglevel_m[i]->phi,
//                             *mglevel_m[i-1]->phi,
//                             mglevel_m[i]->geom,
//                             mglevel_m[i-1]->geom,
//                             rr_m);
//     }
// }
// 
// void AmrMultiGrid::gradient_m(int level, AmrField_t& efield) {
//     const double* dx = mglevel_m[level]->geom.CellSize();
//     
//     for (amrex::MFIter mfi(*mglevel_m[level]->phi, false); mfi.isValid(); ++mfi) {
//         const amrex::Box&          bx   = mfi.validbox();
//         const amrex::FArrayBox&    pfab = (*mglevel_m[level]->phi)[mfi];
//         amrex::FArrayBox&          efab = (efield)[mfi];
//         
//         const int* lo = bx.loVect();
//         const int* hi = bx.hiVect();
//             
//         
//         for (int i = lo[0]; i <= hi[0]; ++i) {
//             for (int j = lo[1]; j <= hi[1]; ++j) {
//                 for (int k = lo[2]; k <= hi[2]; ++k) {
//                     
//                     AmrIntVect_t iv(i, j, k);
//                     
//                     // x direction
//                     AmrIntVect_t left(i-1, j, k);
//                     AmrIntVect_t right(i+1, j, k);
//                     efab(iv, 0) = 0.5 * (pfab(left) - pfab(right)) / dx[0];
//                     
//                     // y direction
//                     AmrIntVect_t up(i, j+1, k);
//                     AmrIntVect_t down(i, j-1, k);
//                     efab(iv, 1) = 0.5 * (pfab(down) - pfab(up)) / dx[1];
//                     
//                     // z direction
//                     AmrIntVect_t back(i, j, k+1);
//                     AmrIntVect_t front(i, j, k-1);
//                     efab(iv, 2) = 0.5 * (pfab(front) - pfab(back)) / dx[2];
//                 }
//             }
//         }
//     }
//     
// }
// 
// 
// void AmrMultiGrid::setBoundaryValue_m(AmrField_t* phi, const AmrGeometry_t& geom, const AmrField_t* crse_phi,
//                                       AmrIntVect_t crse_ratio)
// {
//     
//     std::unique_ptr<amrex::MacBndry> bndry;
//     bndry.reset(new amrex::MacBndry(phi->boxArray(), phi->DistributionMap(), 1 /*ncomp*/, geom));
//     
//     
//     for (amrex::OrientationIter fi; fi; ++fi)
//     {
//         const amrex::Orientation face  = fi();
//         amrex::FabSet& fs = (*bndry.get())[face];
//         
//         for (amrex::FabSetIter fsi(fs); fsi.isValid(); ++fsi)
//         {
//             fs[fsi].setVal(0, 0);
//         }
//     }
//     
//     
//     // The values of phys_bc & ref_ratio do not matter
//     // because we are not going to use those parts of MacBndry.
//     AmrIntVect_t ref_ratio = AmrIntVect_t::TheZeroVector();
//     amrex::Array<int> lo_bc(BL_SPACEDIM, 0);
//     amrex::Array<int> hi_bc(BL_SPACEDIM, 0);
//     amrex::BCRec phys_bc(lo_bc.dataPtr(), hi_bc.dataPtr());
//     
// //         if (crse_phi == 0 && phi == 0) 
// //         {
// //             bndry.setHomogValues(phys_bc, ref_ratio);
// //         }
//     if (crse_phi == 0 && phi != 0)
//     {
//         bndry->setBndryValues(*phi, 0, 0, phi->nComp(), phys_bc); 
//     }
//     else if (crse_phi != 0 && phi != 0) 
//     {
//         const int ncomp      = phi->nComp();
//         const int in_rad     = 0;
//         const int out_rad    = 1;
//         const int extent_rad = 2;
// 
//         amrex::BoxArray crse_boxes(phi->boxArray());
//         crse_boxes.coarsen(crse_ratio);
//     
//         amrex::BndryRegister crse_br(crse_boxes, phi->DistributionMap(),
//                                     in_rad, out_rad, extent_rad, ncomp);
//         crse_br.copyFrom(*crse_phi, crse_phi->nGrow(), 0, 0, ncomp);
//             
//         bndry->setBndryValues(crse_br, 0, *phi, 0, 0, ncomp, crse_ratio, phys_bc, 3);
//         
//         
//     }
//     else
//     {
//         amrex::Abort("FMultiGrid::Boundary::build_bndry: How did we get here?");
//     }
//     
//     for (amrex::OrientationIter oitr; oitr; ++oitr)
//     {
//         const amrex::Orientation face  = oitr();
//         const amrex::FabSet& fs = bndry->bndryValues(oitr());
//         for (amrex::MFIter umfi(*(phi)); umfi.isValid(); ++umfi)
//         {
//             amrex::FArrayBox& dest = (*phi)[umfi];
//             dest.copy(fs[umfi],fs[umfi].box());
//         }
//     }
// }


void AmrMultiGrid::buildPoissonMatrix_m(int level) {
    /*
     * 1D not supported by AmrMultiGrid
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    int nEntries = (BL_SPACEDIM << 2) + 1;
    
    const Epetra_Map& RowMap = *mglevel_m[level]->map_p;
    
    mglevel_m[level]->A_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy, RowMap, nEntries) );
    
    indices_t indices(nEntries);
    coefficients_t values(nEntries);
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    for (amrex::MFIter mfi(*mglevel_m[level]->mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    int numEntries = 0;
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    int globalRow = mglevel_m[level]->serialize(iv);
                    
                    // check left neighbor in x-direction
                    AmrIntVect_t xl(D_DECL(i-1, j, k));
                    if ( mfab(xl) < 1 ) {
                        int gidx = mglevel_m[level]->serialize(xl);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
                        ++numEntries;
                    }
                    
                    // check right neighbor in x-direction
                    AmrIntVect_t xr(D_DECL(i+1, j, k));
                    if ( mfab(xr) < 1 ) {
                        int gidx = mglevel_m[level]->serialize(xr);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
                        ++numEntries;
                    }
                    
                    // check lower neighbor in y-direction
                    AmrIntVect_t yl(D_DECL(i, j-1, k));
                    if ( mfab(yl) < 1 ) {
                        int gidx = mglevel_m[level]->serialize(yl);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
                        ++numEntries;
                    }
                    
                    // check upper neighbor in y-direction
                    AmrIntVect_t yr(D_DECL(i, j+1, k));
                    if ( mfab(yr) < 1 ) {
                        int gidx = mglevel_m[level]->serialize(yr);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
                        ++numEntries;
                    }
                    
                    // check front neighbor in z-direction
                    AmrIntVect_t zl(D_DECL(i, j, k-1));
                    if ( mfab(zl) < 1 ) {
                        int gidx = mglevel_m[level]->serialize(zl);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
                        ++numEntries;
                    }
                    
                    // check back neighbor in z-direction
                    AmrIntVect_t zr(D_DECL(i, j, k+1));
                    if ( mfab(zr) < 1 ) {
                        int gidx = mglevel_m[level]->serialize(zr);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
                        ++numEntries;
                    }
                    
                    // check center
                    if ( mfab(iv) == 0 ) {
                        indices[numEntries] = globalRow;
                        values[numEntries]  = 2.0 / ( dx[0] * dx[0] ) +
                                              2.0 / ( dx[1] * dx[1] ) +
                                              2.0 / ( dx[2] * dx[2] );
                        ++numEntries;
                    }
                    
                    int error = mglevel_m[level]->A_p->InsertGlobalValues(globalRow,
                                                                          numEntries,
                                                                          &values[0],
                                                                          &indices[0]);
                    
                    if ( error != 0 )
                        throw std::runtime_error("Error in filling the Poisson matrix for level "
                                                 + std::to_string(level) + "!");
                }
            }
        }
    }
    
    int error = mglevel_m[level]->A_p->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing Poisson matrix for level "
                                 + std::to_string(level) + "!");
    
    /*
     * some printing
     */
    
#ifdef DEBUG
    if ( Ippl::myNode() == 0 ) {
        std::cout << "Global info" << std::endl
                  << "Number of rows:      " << A_mp->NumGlobalRows() << std::endl
                  << "Number of cols:      " << A_mp->NumGlobalCols() << std::endl
                  << "Number of diagonals: " << A_mp->NumGlobalDiagonals() << std::endl
                  << "Number of non-zeros: " << A_mp->NumGlobalNonzeros() << std::endl
                  << std::endl;
    }
    
    Ippl::Comm->barrier();
    
    for (int i = 0; i < Ippl::getNodes(); ++i) {
        
        if ( i == Ippl::myNode() ) {
            std::cout << "Rank:                "
                      << i << std::endl
                      << "Number of rows:      "
                      << A_mp->NumMyRows() << std::endl
                      << "Number of cols:      "
                      << A_mp->NumMyCols() << std::endl
                      << "Number of diagonals: "
                      << A_mp->NumMyDiagonals() << std::endl
                      << "Number of non-zeros: "
                      << A_mp->NumMyNonzeros() << std::endl
                      << std::endl;
        }
        Ippl::Comm->barrier();
    }
#endif
}


void AmrMultiGrid::buildRestrictionMatrix_m(int level) {
    // base level does not need to have a restriction matrix
    if ( level == 0 )
        return;
    
    /* Difficulty:  If a fine cell belongs to another processor than the underlying
     *              coarse cell, we get an error when filling the matrix since the
     *              cell (--> global index) does not belong to the same processor.
     * Solution:    Find all coarse cells that are covered by fine cells, thus,
     *              the distributionmap is correct.
     * 
     * 
     */
    
    AmrIntVect_t rr = mglevel_m[level-1]->refinement();
    
    amrex::BoxArray crse_fine_ba = mglevel_m[level-1]->grids;
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(mglevel_m[level]->grids, crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    const Epetra_Map& RowMap = *mglevel_m[level-1]->map_p;
    const Epetra_Map& ColMap = *mglevel_m[level]->map_p;
    
#if BL_SPACEDIM == 3
    int nNeighbours = rr[0] * rr[1] * rr[2];
#else
    int nNeighbours = rr[0] * rr[1];
#endif
    
    indices_t indices(nNeighbours);
    coefficients_t values(nNeighbours);
    
    double val = 1.0 / double(nNeighbours);
    
    mglevel_m[level]->R_p = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                                               RowMap, nNeighbours) );
    
    for (amrex::MFIter mfi(mglevel_m[level-1]->grids,
                           mglevel_m[level-1]->dmap, false);
         mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            int ii = i * rr[0];
            for (int j = lo[1]; j <= hi[1]; ++j) {
                int jj = j * rr[1];
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    int kk = k * rr[2];
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    if ( crse_fine_ba.contains(iv) ) {
                        int crse_globidx = mglevel_m[level-1]->serialize(iv);
                    
                        int numEntries = 0;
                    
                        // neighbours
                        for (int iref = 0; iref < rr[0]; ++iref) {
                            for (int jref = 0; jref < rr[1]; ++jref) {
#if BL_SPACEDIM == 3
                                for (int kref = 0; kref < rr[2]; ++kref) {
#endif
                                    AmrIntVect_t riv(D_DECL(ii + iref, jj + jref, kk + kref));
                                    
                                    int fine_globidx = mglevel_m[level]->serialize(riv);
                                    
                                    indices[numEntries] = fine_globidx;
                                    values[numEntries]  = val;
                                    ++numEntries;
#if BL_SPACEDIM == 3
                                }
#endif
                            }
                        }
                    
                        int error = mglevel_m[level]->R_p->InsertGlobalValues(crse_globidx,
                                                                              numEntries,
                                                                              &values[0],
                                                                              &indices[0]);
                        
                        if ( error != 0 ) {
                            throw std::runtime_error("Error in filling the restriction matrix for level " +
                            std::to_string(level) + "!");
                        }
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    if ( mglevel_m[level]->R_p->FillComplete(ColMap, RowMap, true) != 0 )
        throw std::runtime_error("Error in completing the restriction matrix for level " +
                                 std::to_string(level) + "!");
}


void AmrMultiGrid::buildInterpolationMatrix_m(int level) {
    /*
     * This does not include ghost cells of *this* level
     * --> no boundaries
     */
    
    if ( level == lbase_m )
        return;
    
    int nNeighbours = interp_mp->getNumberOfPoints();
    
    indices_t indices(nNeighbours);
    coefficients_t values(nNeighbours);
    
    const Epetra_Map& RowMap = *mglevel_m[level]->map_p;
    const Epetra_Map& ColMap = *mglevel_m[level-1]->map_p;
    
    mglevel_m[level]->I_p = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                                               RowMap, nNeighbours) );
    
    for (amrex::MFIter mfi(mglevel_m[level]->grids,
                           mglevel_m[level]->dmap, false);
         mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx  = mfi.validbox();
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    int globidx = mglevel_m[level]->serialize(iv);
                    
                    int numEntries = 0;
                    
                    /*
                     * we need boundary + indices from coarser level
                     */
                    interp_mp->stencil(iv, indices, values, numEntries, mglevel_m[level-1].get());
                    
                    int error = mglevel_m[level]->I_p->InsertGlobalValues(globidx,
                                                                          numEntries,
                                                                          &values[0],
                                                                          &indices[0]);
                    
                    if ( error != 0 ) {
                        throw std::runtime_error("Error in filling the interpolation matrix for level " +
                                                 std::to_string(level) + "!");
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = mglevel_m[level]->I_p->FillComplete(ColMap, RowMap, true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the interpolation matrix for level " +
                                 std::to_string(level) + "!");
}


void AmrMultiGrid::buildBoundaryMatrix_m(int level) {
    
    /*
     * helper function to fill matrix
     */
    auto insert = [&, this](const amrex::BaseFab<int>& mfab,
                            const AmrIntVect_t& iv,
                            const AmrIntVect_t& in,
                            indices_t& indices,
                            coefficients_t& values)
    {
        if ( level > 0 && mfab(iv) >= 1 /*higher level*/ )
        {
            int globidx = mglevel_m[level]->serialize(in);
            
            int numEntries = 0;
            
            /*
             * we need boundary + indices from coarser level
             */
            interp_mp->stencil(iv, indices, values, numEntries, mglevel_m[level-1].get());
                    
            int error = mglevel_m[level]->B_p->InsertGlobalValues(globidx,
                                                                  numEntries,
                                                                  &values[0],
                                                                  &indices[0]);
            
            
            if ( error != 0 ) {
                throw std::runtime_error("Error in filling the boundary matrix for level " +
                                         std::to_string(level) + "!");
            }
            
//         } else if ( level > 0 && mfab(iv) == 2 ) {
//             // physical boundary
//             std::cout << iv << " " << in << std::endl;
//             
//             throw std::runtime_error("Level " + std::to_string(level) +
//                                      ": Ghost cell cannot be at physical boundary!");
        } else if ( level == 0 && mfab(iv) == 2 ) {
            int globidx = mglevel_m[level]->serialize(in);
            
            int numEntries = 0;
            
            mglevel_m[level]->applyBoundary(iv, indices, values, numEntries, 1.0 /*matrix coefficient*/);
            
            int error = mglevel_m[level]->B_p->InsertGlobalValues(globidx,
                                                                  numEntries,
                                                                  &values[0],
                                                                  &indices[0]);
            
            
            if ( error != 0 ) {
                throw std::runtime_error("Error in filling the boundary matrix for level " +
                                         std::to_string(level) + "!");
            }
            
        }
    };
    
    
    int nNeighbours = interp_mp->getNumberOfPoints();
    
    const Epetra_Map& RowMap = *mglevel_m[level]->map_p;
    
    mglevel_m[level]->B_p = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                                               RowMap, nNeighbours) );
    indices_t indices(nNeighbours);
    coefficients_t values(nNeighbours);
    
    /* In order to avoid that corner cells are inserted twice to Trilinos, leading to
     * an insertion error, left and right faces account for the corner cells.
     * The lower and upper faces only iterate through "interior" boundary cells.
     * The front and back faces need to be adapted as well, i.e. start at the first inner layer
     * all border cells are already treated.
     */
    
    for (amrex::MFIter mfi(*mglevel_m[level]->mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&    bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        /*
         * LEFT BOUNDARY
         */
        int ii = lo[0] - 1; // ghost
        for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                AmrIntVect_t iv(D_DECL(ii, j, k));
                
                // interior IntVect
                AmrIntVect_t in(D_DECL(lo[0], j, k));
                
                insert(mfab, iv, in, indices, values);
                
#if BL_SPACEDIM == 3
            }
#endif
        }
        
        /*
         * RIGHT BOUNDARY
         */
        ii = hi[0] + 1; // ghost
        for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                AmrIntVect_t iv(D_DECL(ii, j, k));
                
                // interior IntVect
                AmrIntVect_t in(D_DECL(hi[0], j, k));
                
                insert(mfab, iv, in, indices, values);
                
#if BL_SPACEDIM == 3
            }
#endif
        }
        
        /*
         * LOWER BOUNDARY
         */
        int jj = lo[1] - 1; // ghost
        
        /*
         * check if corner cell or covered cell by neighbour box
         * 
         * lo: corner cells accounted in "left boundary"
         * hi: corner cells accounted in "right boundary"
         * 
         * if lo[1] on low side is crse/fine interface --> corner
         * if lo[1] on high side is crs/fine interface --> corner
         */
        int start = ( mfab(AmrIntVect_t(D_DECL(lo[0]-1, lo[1], lo[2]))) >= 1 ) ? 1 : 0;
        int end = ( mfab(AmrIntVect_t(D_DECL(hi[0]+1, lo[1], hi[2]))) >= 1 ) ? 1 : 0;
        
        for (int i = lo[0]+start; i <= hi[0]-end; ++i) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                AmrIntVect_t iv(D_DECL(i, jj, k));
                
                // interior IntVect
                AmrIntVect_t in(D_DECL(i, lo[1], k));
                
                insert(mfab, iv, in, indices, values);
#if BL_SPACEDIM == 3
            }
#endif
        }
        
        /*
         * UPPER BOUNDARY
         */
        jj = hi[1] + 1; // ghost
        
        /*
         * check if corner cell or covered cell by neighbour box
         * 
         * lo: corner cells accounted in "left boundary"
         * hi: corner cells accounted in "right boundary"
         * 
         * if hi[1] on low side is crse/fine interface --> corner
         * if hi[1] on high side is crs/fine interface --> corner
         */
        start = ( mfab(AmrIntVect_t(D_DECL(lo[0]-1, hi[1], lo[2]))) >= 1 ) ? 1 : 0;
        end = ( mfab(AmrIntVect_t(D_DECL(hi[0]+1, hi[1], hi[2]))) >= 1 ) ? 1 : 0;
        
        for (int i = lo[0]+start; i <= hi[0]-end; ++i) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                AmrIntVect_t iv(D_DECL(i, jj, k));
                
                // interior IntVect
                AmrIntVect_t in(D_DECL(i, hi[1], k));
                
                insert(mfab, iv, in, indices, values);
                
#if BL_SPACEDIM == 3
            }
#endif
        }
        
#if BL_SPACEDIM == 3
        /*
         * FRONT BOUNDARY
         */
        int kk = lo[2] - 1; // ghost
        for (int i = lo[0]+1; i < hi[0]; ++i) {
            for (int j = lo[1]+1; j < hi[1]; ++j) {
 
                AmrIntVect_t iv(D_DECL(i, j, kk));
                    
                // interior IntVect
                AmrIntVect_t in(D_DECL(i, j, lo[2]));
                
                insert(mfab, iv, in, indices, values);
            }
        }
        
        /*
         * BACK BOUNDARY
         */
        kk = hi[2] + 1; // ghost
        for (int i = lo[0]+1; i < hi[0]; ++i) {
            for (int j = lo[1]+1; j < hi[1]; ++j) {
                
                AmrIntVect_t iv(D_DECL(i, j, kk));
                
                // interior IntVect
                AmrIntVect_t in(D_DECL(i, j, hi[2]));
                
                insert(mfab, iv, in, indices, values);
            }
        }
#endif
    }
    
    int error = 0;
    
    if ( level > 0 ) {
        const Epetra_Map& ColMap = *mglevel_m[level-1]->map_p;
        error = mglevel_m[level]->B_p->FillComplete(ColMap, RowMap, true);
    } else
        error = mglevel_m[level]->B_p->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the boundary matrix for level " +
                                 std::to_string(level) + "!");
}


inline void AmrMultiGrid::buildDensityVector_m(const AmrField_t& rho, int level) {
    this->amrex2trilinos_m(rho, mglevel_m[level]->rho_p, level);
}


inline void AmrMultiGrid::buildPotentialVector_m(const AmrField_t& phi, int level) {
    this->amrex2trilinos_m(phi, mglevel_m[level]->phi_p, level);
}


void AmrMultiGrid::amrex2trilinos_m(const AmrField_t& mf,
                                    Teuchos::RCP<vector_t>& mv, int level)
{
    coefficients_t values;
    indices_t indices;
    
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::FArrayBox&    fab = mf[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    int globalidx = mglevel_m[level]->serialize(iv);
                    
                    indices.push_back(globalidx);
                    
                    values.push_back(fab(iv));
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    if ( mv.is_null() )
        mv = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
    
    int error = mv->ReplaceGlobalValues(mglevel_m[level]->map_p->NumMyElements(),
                                        &values[0],
                                        &indices[0]);
    
    if ( error != 0 )
        throw std::runtime_error("Error in filling the vector!");
}


void AmrMultiGrid::trilinos2amrex_m(AmrField_t& mf,
                                    const Teuchos::RCP<vector_t>& mv)
{
    int localidx = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        amrex::FArrayBox&          fab = mf[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    fab(iv) = (*mv)[localidx++];
                }
#if BL_SPACEDIM == 3
            }
#endif
        }
    }
}
