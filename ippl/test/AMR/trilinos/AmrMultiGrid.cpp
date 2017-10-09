#include "AmrMultiGrid.h"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Array.H>

#include <AMReX_MacBndry.H>

AmrMultiGrid::AmrMultiGrid(Interpolater interp,
                           Interpolater interface)
    : epetra_comm_m(Ippl::getComm())
{
    switch ( interp ) {
        case Interpolater::TRILINEAR:
            interp_mp.reset( new AmrTrilinearInterpolater<AmrMultiGridLevel_t>() );
            break;
        case Interpolater::LAGRANGE:
            std::runtime_error("Not yet implemented.");
        default:
            std::runtime_error("No such interpolater available.");
    }
    
    // interpolater for crse-fine-interface
    switch ( interface ) {
        case Interpolater::TRILINEAR:
            interface_mp.reset( new AmrTrilinearInterpolater<AmrMultiGridLevel_t>() );
            break;
        case Interpolater::LAGRANGE:
            interface_mp.reset( new AmrLagrangeInterpolater<AmrMultiGridLevel_t>(
                AmrLagrangeInterpolater<AmrMultiGridLevel_t>::Order::QUADRATIC) );
            break;
        default:
            std::runtime_error("No such interpolater for the coarse-fine interface available.");
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
    
    rr_m = AmrIntVect_t(D_DECL(2, 2, 2));
    
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
    
    if ( !previous ) {
        // reset
        for (int lev = 0; lev < nLevel; ++lev) {
            int ilev = lbase + lev;
            phi[ilev]->setVal(0.0, phi[ilev]->nGrow());
        }
    }
    
    // build all necessary matrices and vectors
    for (int lev = 0; lev < nLevel; ++lev) {
        
        this->buildRestrictionMatrix_m(lev);
        
        this->buildInterpolationMatrix_m(lev);
        
        this->buildCrseBoundaryMatrix_m(lev);
        
        this->buildFineBoundaryMatrix_m(lev);
        
        this->buildSmootherMatrix_m(lev);
        
        this->buildPoissonMatrix_m(lev);
        
        this->buildSpecialPoissonMatrix_m(lev);
        
        int ilev = lbase + lev;
        
        this->buildDensityVector_m(*rho[ilev], lev);
        
        this->buildPotentialVector_m(*phi[ilev], lev);
        
        this->buildGradientMatrix_m(lev);
        
        mglevel_m[lev]->error_p->PutScalar(0.0);
    }
    
    
    double err = 1.0e7;
    double eps = 1.0e-8;
    
    // initial error
    double residualsum = 0.0;
    double rhosum = 0.0;
    
    
    std::vector<double> old_norm_p;
    for (int lev = 0; lev < nLevel; ++lev) {
        double tmp = 0;
        mglevel_m[lev]->phi_p->Norm2(&tmp);
        old_norm_p.push_back(tmp);
    }
    
    for (int lev = 0; lev < nLevel; ++lev) {
        this->residual_m(mglevel_m[lev]->residual_p,
                         mglevel_m[lev]->rho_p,
//                          mglevel_m[lev]->A_p,
                         mglevel_m[lev]->phi_p, lev);
        
        double tmp = 0;
        mglevel_m[lev]->residual_p->Norm2(&tmp);
        residualsum += tmp;
        
        tmp = 0;
        mglevel_m[lev]->rho_p->Norm2(&tmp);
        rhosum += tmp;
    }
    
    
    while ( residualsum > eps * rhosum) {
        
        std::cout << "                               "
                  << residualsum << " " << eps * rhosum << std::endl; //std::cin.get();
        
        relax_m(lfine);
        
        // update residual
        for (int lev = 0; lev < nLevel; ++lev) {
            
            this->residual_m(mglevel_m[lev]->residual_p,
                             mglevel_m[lev]->rho_p,
//                              mglevel_m[lev]->A_p,
                             mglevel_m[lev]->phi_p, lev);
        }
        
        residualsum = l2error_m();
        
        
        for (int lev = 0; lev < nLevel; ++lev) {
            double tmp = 0;
            mglevel_m[lev]->phi_p->Norm2(&tmp);
            
            std::cout << "level " << lev << ": " << old_norm_p[lev] << " " << tmp << std::endl; std::cin.get();
            
            old_norm_p[lev] = tmp;
        } //std::cin.get();
        
        ++nIter_m;
    }
    
//         std::cin.get();
// //     
// //     for (int lev = nLevel-1; lev > -1; --lev) {
// //         int ilev = lbase + lev;
// //         mglevel_m[lev]->phi->FillBoundary(mglevel_m[lev]->geom.periodicity());
// //         this->gradient_m(lev, *efield[ilev]);
// //     }
// //     
    std::cout << "                               "
              << residualsum << " " << eps * rhosum << std::endl;
    
    // copy solution back
    for (int lev = 0; lev < nLevel; ++lev) {
        int ilev = lbase + lev;
        
        this->trilinos2amrex_m(*phi[ilev], mglevel_m[lev]->phi_p);
    }
    
    
    // evaluate the electric field
    Teuchos::RCP<vector_t> efield_p = Teuchos::null;
    for (int lev = 0; lev < nLevel; ++lev) {
        int ilev = lbase + lev;
        
        efield_p = Teuchos::rcp( new vector_t(*mglevel_m[lev]->map_p, false) );
        
        for (int d = 0; d < BL_SPACEDIM; ++d) {
            mglevel_m[lev]->G_p[d]->Multiply(false, *mglevel_m[lev]->phi_p, *efield_p);
            this->trilinos2amrex_m(*efield[ilev], efield_p);
        }
    }
    
    std::cout << "#iterations: " << nIter_m << std::endl;
}


void AmrMultiGrid::residual_m(Teuchos::RCP<vector_t>& r,
                              const Teuchos::RCP<vector_t>& b,
//                               const Teuchos::RCP<matrix_t>& A,
                              const Teuchos::RCP<vector_t>& x,
                              int level)
{
    /*
     * r = b - A*x
     */
    if ( level < lfine_m ) {
        std::cout << "level = " << level << std::endl;
        
        vector_t fine2crse(mglevel_m[level]->A_p->OperatorDomainMap());
        fine2crse.PutScalar(0.0);
        
        // get boundary for 
        if ( mglevel_m[level]->Bfine_p != Teuchos::null ) {
            mglevel_m[level]->Bfine_p->Apply(*mglevel_m[level+1]->phi_p, fine2crse);
        }
        
        vector_t tmp2(mglevel_m[level]->As_p->OperatorDomainMap());
        mglevel_m[level]->As_p->Apply(*mglevel_m[level]->phi_p, tmp2);
        
        vector_t crse2fine(mglevel_m[level]->A_p->OperatorDomainMap());
        crse2fine.PutScalar(0.0);
        
        if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
            mglevel_m[level]->Bcrse_p->Apply(*mglevel_m[level-1]->phi_p, crse2fine);
        }
        
        tmp2.Update(1.0, fine2crse, 1.0, crse2fine, 1.0);
        
    //     std::cout << "restrict residual: " << *mglevel_m[level]->residual_p << std::endl;
        
        
        Teuchos::RCP<vector_t> tmp3 = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
    
        mglevel_m[level]->UnCovered_p->Apply(*mglevel_m[level]->rho_p, *tmp3);
    
//     std::cout << *mglevel_m[level]->residual_p << std::endl;
    
        // ONLY subtract coarse rho
        mglevel_m[level]->residual_p->Update(1.0, *tmp3, -1.0, tmp2, 0.0);
        
//         std::cout << *mglevel_m[level]->residual_p << std::endl; std::cin.get();
        
//         mglevel_m[level]->residual_p->Update(1.0, *mglevel_m[level]->rho_p, -1.0, tmp2, 0.0);
        
    } else {
        this->residual_no_fine_m(mglevel_m[level]->residual_p,
                                 mglevel_m[level]->phi_p,
                                 mglevel_m[level-1]->phi_p,
                                 mglevel_m[level]->rho_p, level);
        
    }
    
    
//     if ( temporary ) {
//         vector_t tmp(A->OperatorDomainMap());
//         
//         
//         //----------------
//         double bla;
//         x->Norm2(&bla);
//         std::cout << "residual: Norm2(x) = " << bla << std::endl;
//         //----------------
//     
//         // tmp = A * x
//         A->Apply(*x, tmp);
//         
//         //----------------
//         tmp.Norm2(&bla);
//         std::cout << "residual: Norm2(tmp) = " << bla << std::endl;
//         //----------------
//     
//         // r = b - tmp
//         r->Update(1.0, *b, -1.0, tmp, 0.0);
//         
//         //----------------
//         r->Norm2(&bla);
//         std::cout << "residual: Norm2(r) = " << bla << std::endl;
//         //----------------
//         
//     } else {
//         // r = A * x
//         A->Apply(*x, *r);
//         
//         // r = b - r;
//         r->Update(1.0, *b, -1.0);
//     }
}


void AmrMultiGrid::relax_m(int level) {
    
    std::cout << "RELAX LEVEL " << level << std::endl; std::cin.get();
    
    if ( level == lfine_m ) {
        
        if ( level == lbase_m ) {
            Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
            tmp->PutScalar(0.0);
            this->residual_no_fine_m(mglevel_m[level]->residual_p,
                                     mglevel_m[level]->phi_p,
                                     tmp,
                                     mglevel_m[level]->rho_p, level);
        } else {
            //TODO residual is already computed at beginning
            
//             std::cout << "computing residual" << std::endl; std::cin.get();
            
//             std::cout << *mglevel_m[level]->residual_p << std::endl;
            
            this->residual_no_fine_m(mglevel_m[level]->residual_p,
                                     mglevel_m[level]->phi_p,
                                     mglevel_m[level-1]->phi_p,
                                     mglevel_m[level]->rho_p, level);
            
//             std::cout << *mglevel_m[level]->residual_p << std::endl; std::cin.get();
            
        }
    }
    
    if ( level > 0 ) {
        //TODO
        
        // phi^(l, save) = phi^(l)
        vector_t phi_save = *mglevel_m[level]->phi_p;
        
        mglevel_m[level-1]->error_p->PutScalar(0.0);
        
        // smoothing
        for (int iii = 0; iii < 4; ++iii)
            this->gsrb_level_m(mglevel_m[level]->error_p,
                               mglevel_m[level]->residual_p, level);
        
        // phi = phi + e
        mglevel_m[level]->phi_p->Update(1.0, *mglevel_m[level]->error_p, 1.0);
        
        //----------------
        double bla;
        mglevel_m[level]->error_p->Norm2(&bla);
        std::cout << "1. Norm2(error) = " << bla << std::endl;
        //----------------
        
        /*
         * restrict
         */
        this->restrict_m(level-1);
        
        this->relax_m(level - 1);
        
        /*
         * prolongate / interpolate
         */
        
        // interpolate error from l-1 to l
        Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
        tmp->PutScalar(0.0);
        mglevel_m[level-1]->I_p->Multiply(false, *mglevel_m[level-1]->error_p, *tmp);
        
        // e^(l) += tmp
        mglevel_m[level]->error_p->Update(1.0, *tmp, 1.0);
        
        //----------------
        mglevel_m[level]->error_p->Norm2(&bla);
        std::cout << "2. Norm2(error) = " << bla << std::endl; //std::cin.get();
        //----------------
        
        //----------------
        mglevel_m[level-1]->error_p->Norm2(&bla);
        std::cout << "3. Norm2(error) = " << bla << std::endl; //std::cin.get();
        //----------------
        
        
        //----------------
        mglevel_m[level]->residual_p->Norm2(&bla);
        std::cout << "4. Norm2(residual) = " << bla << std::endl; //std::cin.get();
        //----------------
        
        std::cout << "l-1 = " << level-1 << std::endl;
        
        // residual update
        tmp = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
        tmp->PutScalar(0.0);
        this->residual_no_fine_m(tmp,
                                 mglevel_m[level]->error_p,
                                 mglevel_m[level-1]->error_p,
                                 mglevel_m[level]->residual_p,
                                 level);
        
        *mglevel_m[level]->residual_p = *tmp;
        
        //----------------
        mglevel_m[level]->residual_p->Norm2(&bla);
        std::cout << "5. Norm2(residual) = " << bla << std::endl; //std::cin.get();
        //----------------
        
        // delta error
        Teuchos::RCP<vector_t> derror = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
        derror->PutScalar(0.0);
        
        // smoothing
        for (int iii = 0; iii < 4; ++iii)
            this->gsrb_level_m(derror, mglevel_m[level]->residual_p, level);
        
        // e^(l) += de^(l)
        mglevel_m[level]->error_p->Update(1.0, *derror, 1.0);
        
        //----------------
        mglevel_m[level]->error_p->Norm2(&bla);
        std::cout << "6. Norm2(error) = " << bla << std::endl;
        //----------------
        
        // phi^(l) = phi^(l, save) + e^(l)
        mglevel_m[level]->phi_p->Update(1.0, phi_save, 1.0, *mglevel_m[level]->error_p, 0.0);
        
        
        //----------------
        mglevel_m[level]->residual_p->Norm2(&bla);
        std::cout << "7. Norm2(residual) = " << bla << std::endl; //std::cin.get();
        //----------------
        
    } else {
        std::cout << "LEVEL = 0" << std::endl;
        // e = A^(-1)r
        
        mglevel_m[level]->A_p->Scale(-1.0);
        mglevel_m[level]->residual_p->Scale(-1.0);
        
        solver_m.solve(mglevel_m[level]->A_p,
                       mglevel_m[level]->error_p,
                       mglevel_m[level]->residual_p);
        
        mglevel_m[level]->A_p->Scale(-1.0);
        mglevel_m[level]->residual_p->Scale(-1.0);
        
        //----------------
        double bla;
        mglevel_m[level]->error_p->Norm2(&bla);
        std::cout << "8. Norm2(error) = " << bla << std::endl;
        //----------------
        
        // phi = phi + e
        
        mglevel_m[level]->phi_p->Update(1.0, *mglevel_m[level]->error_p, 1.0);
        
        //----------------
        mglevel_m[level]->phi_p->Norm2(&bla);
        std::cout << "9. (after) Norm2(phi) = " << bla << std::endl;
        
//         std::cout << "PHI level 0" << std::endl; std::cin.get();
        
//         std::cout << *mglevel_m[level]->phi_p << std::endl;
        
        //----------------
    }
}


void AmrMultiGrid::residual_no_fine_m(Teuchos::RCP<vector_t>& result,
                                      const Teuchos::RCP<vector_t>& rhs,
                                      const Teuchos::RCP<vector_t>& crs_rhs,
                                      const Teuchos::RCP<vector_t>& b,
                                      int level)
{
    vector_t crse2fine(mglevel_m[level]->A_p->OperatorDomainMap());
    crse2fine.PutScalar(0.0);
    
    // get boundary for 
    if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
        mglevel_m[level]->Bcrse_p->Apply(*crs_rhs, crse2fine);
    }
    
//     //----------------
//     double bla = 0.0;
//     crse2fine.Norm2(&bla);
//     std::cout << "Norm2(crse2fine) = " << bla << std::endl;
//     //----------------
    
    
    // tmp = A * x
    vector_t tmp(mglevel_m[level]->A_p->OperatorDomainMap());
    mglevel_m[level]->A_p->Multiply(false, *rhs, tmp);
    
    
    if ( level > lbase_m ) {
        vector_t t2mp(mglevel_m[level]->A_p->OperatorDomainMap());
        mglevel_m[level]->B_p->Multiply(false, *rhs, t2mp);
    
        tmp.Update(1.0, t2mp, 1.0);
    }
    
//     //----------------
//     tmp.Norm2(&bla);
//     std::cout << "Norm2(tmp) = " << bla << std::endl;
//     //----------------
    
    tmp.Update(1.0, crse2fine, 1.0);
    
//     //----------------
//     tmp.Norm2(&bla);
//     std::cout << "Norm2(tmp) = " << bla << std::endl; std::cin.get();
//     //----------------
    
    result->Update(1.0, *b, -1.0, tmp, 0.0);
}


double AmrMultiGrid::l2error_m() {
    int nLevel = lfine_m - lbase_m + 1;
    
    double err = 0.0;
    
    for (int lev = 0; lev < nLevel; ++lev) {
        
        double tmp = 0;
        mglevel_m[lev]->residual_p->Norm2(&tmp);
        err += tmp;
    }
    
    return err;
}


void AmrMultiGrid::buildPoissonMatrix_m(int level) {
    /*
     * Laplacian of "no fine"
     */
    
    /*
     * 1D not supported
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    int nEntries = (BL_SPACEDIM << 1) + 1 /* plus boundaries */ + 10 /*FIXME*/;
    
    amrex::BoxArray empty;
    
    const Epetra_Map& RowMap = *mglevel_m[level]->map_p;
    
    mglevel_m[level]->A_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy, RowMap, nEntries) );
    
    mglevel_m[level]->B_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy, RowMap, 10 /*FIXME*/) );
    
    indices_t indices, bindices;
    coefficients_t values, bvalues;
    
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
                    indices.clear();
                    values.clear();
                    
                    bindices.clear();
                    bvalues.clear();
                    
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    int globalidx = mglevel_m[level]->serialize(iv);
                    
                    /*
                     * check neighbours in all directions (Laplacian stencil --> cross)
                     */
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
                        for (int shift = -1; shift <= 1; shift += 2) {
                            AmrIntVect_t biv = iv;                        
                            biv[d] += shift;
                            
                            switch ( mfab(biv) )
                            {
                                case AmrMultiGridLevel_t::Mask::COVERED:
                                    // covered --> interior cell
                                case AmrMultiGridLevel_t::Mask::INTERIOR:
                                {
                                    indices.push_back( mglevel_m[level]->serialize(biv) );
                                    values.push_back( 1.0 / ( dx[d] * dx[d] ) );
                                    break;
                                }
                                case AmrMultiGridLevel_t::Mask::BNDRY:
                                {
                                    // boundary cell
                                    // only level > 0 have this kind of boundary
                                    
                                    /* Dirichlet boundary conditions from coarser level.
                                     */
                                    std::size_t nn = bindices.size();
                                    
                                    interface_mp->fine(biv, bindices, bvalues, d, -shift, empty,
                                                       mglevel_m[level].get());
                                    
                                    double value = 1.0 / ( dx[d] * dx[d] );
                                    for (std::size_t iter = nn; iter < bindices.size(); ++iter)
                                        bvalues[iter] *= value;
                                    
                                    break;
                                }
                                case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
                                {
                                    // physical boundary cell
                                    double value = 1.0 / ( dx[d] * dx[d] );
                                    mglevel_m[level]->applyBoundary(biv,
                                                                    indices,
                                                                    values,
                                                                    value /*matrix coefficient*/);
                                    break;
                                }
                                default:
                                    throw std::runtime_error("Error in mask for level "
                                                             + std::to_string(level) + "!");
                                    break;
                            }
                        }
                    }
                    
                    // check center
                    indices.push_back( globalidx );
                    values.push_back( -2.0 / ( dx[0] * dx[0] ) +
                                      -2.0 / ( dx[1] * dx[1] )
#if BL_SPACEDIM == 3
                                      -2.0 / ( dx[2] * dx[2] )
#endif
                        );
                    
                    this->unique_m(indices, values);
                    
                    int error = mglevel_m[level]->A_p->InsertGlobalValues(globalidx,
                                                                          indices.size(),
                                                                          &values[0],
                                                                          &indices[0]);
                    
                    if ( error != 0 )
                        throw std::runtime_error("Error in filling the Poisson matrix for level "
                                                 + std::to_string(level) + "!");
                    
                    error = mglevel_m[level]->B_p->InsertGlobalValues(globalidx,
                                                                      bindices.size(),
                                                                      &bvalues[0],
                                                                      &bindices[0]);
                    
                    if ( error != 0 )
                        throw std::runtime_error("Error in filling the Poisson matrix for level "
                                                 + std::to_string(level) + "!");
                
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = mglevel_m[level]->A_p->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing Poisson matrix for level "
                                 + std::to_string(level) + "!");
    
    error = mglevel_m[level]->B_p->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing Poisson matrix for level "
                                 + std::to_string(level) + "!");
}


void AmrMultiGrid::buildSpecialPoissonMatrix_m(int level) {
    /*
     * Laplacian of "with fine"
     */
    if ( level == lfine_m )
        return;
    
    
    // find all cells with refinement
    AmrIntVect_t rr = mglevel_m[level]->refinement();
    
    const double* cdx = mglevel_m[level]->geom.CellSize();
    const double* fdx = mglevel_m[level+1]->geom.CellSize();
    
    amrex::BoxArray crse_fine_ba = mglevel_m[level]->grids;
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(mglevel_m[level+1]->grids, crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    /*
     * 1D not supported by AmrMultiGrid
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    int nEntries = (BL_SPACEDIM << 1) + 1 /* plus boundaries */ + 10 /*FIXME*/;
    
    const Epetra_Map& RowMap = *mglevel_m[level]->map_p;
    
    mglevel_m[level]->As_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy, RowMap, nEntries) );
    
    
     mglevel_m[level]->UnCovered_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy, RowMap, 1) );
    
    indices_t indices;
    coefficients_t values;
    
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
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    if ( !crse_fine_ba.contains(iv) ) {
                        indices.clear();
                        values.clear();
                        int globalidx = mglevel_m[level]->serialize(iv);
                        
                        /*
                         * check neighbours in all directions (Laplacian stencil --> cross)
                         */
                        for (int d = 0; d < BL_SPACEDIM; ++d) {
                            for (int shift = -1; shift <= 1; shift += 2) {
                                AmrIntVect_t biv = iv;                        
                                biv[d] += shift;
                                
                                if ( !crse_fine_ba.contains(biv) )
                                {
                                    /*
                                     * It can't be a refined cell!
                                     */
                                    switch ( mfab(biv) )
                                    {
                                        case -1:
                                            // covered --> interior cell
                                        case 0:
                                        {
                                            indices.push_back( mglevel_m[level]->serialize(biv) );
                                            values.push_back( 1.0 / ( dx[d] * dx[d] ) );
                                            
                                            // add center once
                                            indices.push_back( globalidx );
                                            values.push_back( -1.0 / ( dx[d] * dx[d] ) );
                                            
                                            break;
                                        }
                                        case 1:
                                        {
                                            // boundary cell
                                            // only level > 0 have this kind of boundary
                                            
                                            /* Dirichlet boundary conditions from coarser level.
                                             */
                                            std::size_t nn = indices.size();
                                            
                                            std::cout << iv << " " << biv << " ";
                                            
                                            interface_mp->fine(biv, indices, values, d, -shift, crse_fine_ba,
                                                               mglevel_m[level].get());
                                            
                                            double value = 1.0 / ( dx[d] * dx[d] );
                                            for (std::size_t iter = nn; iter < indices.size(); ++iter) {
                                                std::cout << indices[iter] << " ";
                                                values[iter] *= value;
                                            }
                                            
                                            // add center once
                                            indices.push_back( globalidx );
                                            values.push_back( -1.0 / ( dx[d] * dx[d] ) );
                                            
                                            std::cout << std::endl; std::cin.get();
                                            
                                            break;
                                        }
                                        case 2:
                                        {
                                            // physical boundary cell
                                            double value = 1.0 / ( dx[d] * dx[d] );
                                            
                                            mglevel_m[level]->applyBoundary(biv,
                                                                    indices,
                                                                    values,
                                                                    value /*matrix coefficient*/);
                                            
                                            // add center once
                                            indices.push_back( globalidx );
                                            values.push_back( -1.0 / ( dx[d] * dx[d] ) );
                                            
                                            break;
                                        }
                                        default:
                                            break;
                                    }
                                } else {
                                    // flux matching, coarse part
                                    
                                    double avg = 2; //FIXME
                                    double value = -1.0 / ( avg * cdx[d] * fdx[d] );
                                    
                                    // top
                                    std::size_t nn = indices.size();
                                    interface_mp->coarse(iv, indices, values, d, shift, crse_fine_ba,
                                                         true, mglevel_m[level].get());
                                    
                                    for (std::size_t iter = nn; iter < indices.size(); ++iter) {
                                        values[iter] *= value;
                                    }
                                    
                                    // bottom
                                    nn = indices.size();
                                    interface_mp->coarse(iv, indices, values, d, shift, crse_fine_ba,
                                                         false, mglevel_m[level].get());
                                    
                                    for (std::size_t iter = nn; iter < indices.size(); ++iter) {
                                        values[iter] *= value;
                                    }
                                    
                                }
                            }
                        }
                        
                        this->unique_m(indices, values);
                        
                        int error = mglevel_m[level]->As_p->InsertGlobalValues(globalidx,
                                                                               indices.size(),
                                                                               &values[0],
                                                                               &indices[0]);
                    
                        if ( error != 0 )
                            throw std::runtime_error("Error in filling the Poisson matrix for level "
                                                     + std::to_string(level) + "!");
                        
                        double vv = 1.0;
                        error = mglevel_m[level]->UnCovered_p->InsertGlobalValues(globalidx,
                                                                                  1,
                                                                                  &vv,
                                                                                  &globalidx);
                        if ( error != 0 )
                            throw std::runtime_error("Error in filling the uncovered matrix for level "
                                                     + std::to_string(level) + "!");
                        
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = mglevel_m[level]->As_p->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing special Poisson matrix for level "
                                 + std::to_string(level) + "!");
    
    error = mglevel_m[level]->UnCovered_p->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing uncovered matrix for level "
                                 + std::to_string(level) + "!");
        
    std::cout << "Done special" << std::endl;
}


void AmrMultiGrid::buildRestrictionMatrix_m(int level) {
    /*
     * x^(l-1) = R * x^(l)
     */
    
    // base level does not need to have a restriction matrix
    if ( level == lbase_m )
        return;
    
    /* Difficulty:  If a fine cell belongs to another processor than the underlying
     *              coarse cell, we get an error when filling the matrix since the
     *              cell (--> global index) does not belong to the same processor.
     * Solution:    Find all coarse cells that are covered by fine cells, thus,
     *              the distributionmap is correct.
     * 
     * 
     */
    
    
    // find all coarse cells that are covered by fine cells
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
    
    mglevel_m[level]->R_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy,
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
                        int crse_globalidx = mglevel_m[level-1]->serialize(iv);
                    
                        int numEntries = 0;
                    
                        // neighbours
                        for (int iref = 0; iref < rr[0]; ++iref) {
                            for (int jref = 0; jref < rr[1]; ++jref) {
#if BL_SPACEDIM == 3
                                for (int kref = 0; kref < rr[2]; ++kref) {
#endif
                                    AmrIntVect_t riv(D_DECL(ii + iref, jj + jref, kk + kref));
                                    
                                    int fine_globalidx = mglevel_m[level]->serialize(riv);
                                    
                                    indices[numEntries] = fine_globalidx;
                                    values[numEntries]  = val;
                                    ++numEntries;
#if BL_SPACEDIM == 3
                                }
#endif
                            }
                        }
                    
                        int error = mglevel_m[level]->R_p->InsertGlobalValues(crse_globalidx,
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
    /* crse (this): level
     * fine: level + 1
     */
    
    /*
     * This does not include ghost cells of *this* level
     * --> no boundaries
     * 
     * x^(l+1) = I * x^(l)
     */
    
    if ( level == lfine_m )
        return;
    
    int nNeighbours = interp_mp->getNumberOfPoints() + 10 /*FIXME*/;
    
    indices_t indices;
    coefficients_t values;
    
    const Epetra_Map& RowMap = *mglevel_m[level+1]->map_p;
    const Epetra_Map& ColMap = *mglevel_m[level]->map_p;
    
    mglevel_m[level]->I_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy,
                                                       RowMap, nNeighbours, false) );
    
    for (amrex::MFIter mfi(mglevel_m[level+1]->grids,
                           mglevel_m[level+1]->dmap, false);
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
                    
                    int globalidx = mglevel_m[level+1]->serialize(iv);
                    
                    indices.clear();
                    values.clear();
                    
                    /*
                     * we need boundary + indices from coarser level
                     */
                    interp_mp->stencil(iv, indices, values, mglevel_m[level].get());
                    
                    this->unique_m(indices, values);
                    
                    int error = mglevel_m[level]->I_p->InsertGlobalValues(globalidx,
                                                                          indices.size(),
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


void AmrMultiGrid::buildCrseBoundaryMatrix_m(int level) {
    
    /*
     * fine (this): level
     * coarse:      level - 1
     */
    
    // the base level has only phsysical boundaries
    if ( level == lbase_m )
        return;
    
    
    // find all cells with refinement
    AmrIntVect_t rr = mglevel_m[level-1]->refinement();
    
    amrex::BoxArray crse_fine_ba = mglevel_m[level-1]->grids;
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(mglevel_m[level]->grids, crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    const Epetra_Map& RowMap = *mglevel_m[level]->map_p;
    const Epetra_Map& ColMap = *mglevel_m[level-1]->map_p;
    
    int nNeighbours = (2 << (BL_SPACEDIM -1 )) /*FIXME interpolation stencil indices*/ + 10;
    
    mglevel_m[level]->Bcrse_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy,
                                                           RowMap, nNeighbours, false) );
    
    for (amrex::MFIter mfi(*mglevel_m[level]->mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&    bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        /*
         * left boundary
         */
        AmrIntVect_t lower(D_DECL(lo[0], lo[1], lo[2]));
        AmrIntVect_t upper(D_DECL(lo[0], hi[1], hi[2]));
        
        this->checkCrseBoundary_m(mglevel_m[level]->Bcrse_p, level, mfab, lower, upper, crse_fine_ba);
        
        /*
         * right boundary
         */
        lower = AmrIntVect_t(D_DECL(hi[0], lo[1], lo[2]));
        upper = AmrIntVect_t(D_DECL(hi[0], hi[1], hi[2]));
        
        this->checkCrseBoundary_m(mglevel_m[level]->Bcrse_p, level, mfab, lower, upper, crse_fine_ba);
        
        
        /*
         * lower boundary (without left and right last cell!)
         */
        lower = AmrIntVect_t(D_DECL(lo[0]+1, lo[1], lo[2]));
        upper = AmrIntVect_t(D_DECL(hi[0]-1, lo[1], hi[2]));
        
        this->checkCrseBoundary_m(mglevel_m[level]->Bcrse_p, level, mfab, lower, upper, crse_fine_ba);
        
        /*
         * upper boundary (without left and right last cell!)
         */
        lower = AmrIntVect_t(D_DECL(lo[0]+1, hi[1], lo[2]));
        upper = AmrIntVect_t(D_DECL(hi[0]-1, hi[1], hi[2]));
        
        this->checkCrseBoundary_m(mglevel_m[level]->Bcrse_p, level, mfab, lower, upper, crse_fine_ba);
        
#if BL_SPACEDIM == 3
        /*
         * front boundary (without left, right, upper and lower last cell!)
         */
        lower = AmrIntVect_t(D_DECL(lo[0]+1, lo[1]+1, lo[2]));
        upper = AmrIntVect_t(D_DECL(hi[0]-1, hi[1]-1, lo[2]));
        
        this->checkCrseBoundary_m(mglevel_m[level]->Bcrse_p, level, mfab, lower, upper, crse_fine_ba);
        
        /*
         * back boundary (without left, right, upper and lower last cell!)
         */
        lower = AmrIntVect_t(D_DECL(lo[0]+1, lo[1]+1, hi[2]));
        upper = AmrIntVect_t(D_DECL(hi[0]-1, hi[1]-1, hi[2]));
        
        this->checkCrseBoundary_m(mglevel_m[level]->Bcrse_p, level, mfab, lower, upper, crse_fine_ba);
#endif
        
    }
    
    int error = mglevel_m[level]->Bcrse_p->FillComplete(ColMap, RowMap, true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the crse boundary matrix for level " +
                                 std::to_string(level) + "!");
}


void AmrMultiGrid::buildFineBoundaryMatrix_m(int level)
{
    /* fine: level + 1
     * coarse (this): level
     */
    
    // the finest level does not need data from a finer level
    if ( level == lfine_m )
        return;
    
    // find all cells with refinement
    AmrIntVect_t rr = mglevel_m[level]->refinement();
    
    amrex::BoxArray crse_fine_ba = mglevel_m[level]->grids;
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(mglevel_m[level+1]->grids, crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    const Epetra_Map& RowMap = *mglevel_m[level]->map_p;
    const Epetra_Map& ColMap = *mglevel_m[level+1]->map_p;
    
    int nNeighbours = 4 /*#interfaces*/ * 2; //* rr[0] * rr[1] /*of refined cell*/;
#if BL_SPACEDIM == 3
    //FIXME
    nNeighbours = 6 /*#interfaces*/ * rr[0] * rr[1] * rr[2] /*of refined cell*/;
#endif
    
    mglevel_m[level]->Bfine_p = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                                                   RowMap, nNeighbours, false) );
    
    indices_t indices;
    coefficients_t values;
    
    const double* cdx = mglevel_m[level]->geom.CellSize();
    const double* fdx = mglevel_m[level+1]->geom.CellSize();
    
    /* Find all coarse cells that are at the crse-fine interace but are
     * not refined.
     * Put them into a map (to avoid duplicates, e.g. due to corners).
     * Finally, iterate through the list of cells
     */
    
    // std::list<std::pair<int, int> > contains the shift and direction to come to the covered cell
    map_t cells;
    
    for (amrex::MFIter mfi(mglevel_m[level]->grids, mglevel_m[level]->dmap, false);
         mfi.isValid(); ++mfi)
    {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    if ( !crse_fine_ba.contains(iv) && mfab(iv) != 2 /*not physical bc*/ ) {
                        /*
                         * iv is a coarse cell that got not refined
                         * 
                         * --> check all neighbours to see if at crse-fine
                         * interface
                         */
                        for (int d = 0; d < BL_SPACEDIM; ++d) {
                            for (int shift = -1; shift <= 1; shift += 2) {
                                // neighbour
                                AmrIntVect_t covered = iv;
                                covered[d] += shift;
                                
                                if ( crse_fine_ba.contains(covered) &&
                                     !mglevel_m[level]->isBoundary(covered) /*not physical bc*/ )
                                {
                                    // neighbour is covered by fine cells
                                    
                                    // insert if not yet exists
                                    cells[iv].push_back(std::make_pair(shift, d));
                                }
                            }
                        }
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    
    auto fill = [&](indices_t& indices,
                    coefficients_t& values,
                    D_DECL(int ii, int jj, int kk),
                    int* begin, int* end, int d,
                    const AmrIntVect_t& iv, int shift,
                    int sign)
    {
        // number of fine cell gradients
        double avg = 1.0;
        switch ( d ) {
            case 0:
                // horizontal
                avg = rr[1];
#if BL_SPACEDIM == 3
                avg *= rr[2];
#endif
                break;
            case 1:
                // vertical
                avg = rr[0];
#if BL_SPACEDIM == 3
                avg *= rr[2];
#endif
                break;
#if BL_SPACEDIM == 3
            case 2:
                avg = rr[0] * rr[1];
                break;
#endif
            default:
                throw std::runtime_error("This dimension does not exist!");
                break;
        }
        
        for (int iref = ii - begin[0]; iref <= ii + end[0]; ++iref) {
            
            sign *= ( d == 0 ) ? -1.0 : 1.0;
            
            for (int jref = jj - begin[1]; jref <= jj + end[1]; ++jref) {
                
                sign *= ( d == 1 ) ? -1.0 : 1.0;
#if BL_SPACEDIM == 3
                for (int kref = kk - begin[2]; kref <= kk + end[2]; ++kref) {
#endif
                    /* Since all fine cells on the not-refined cell are
                     * outside of the "domain" --> we need to interpolate
                     */
                    sign *= ( d == 2 ) ? -1.0 : 1.0;
                    
                    AmrIntVect_t riv(D_DECL(iref, jref, kref));
                    
//                     double value = double(shift) * sign / ( avg * cdx[d] * fdx[d] );
                    
                    
//                     std::cout << d << " riv = " << riv
//                               << " iv = " << iv << " "
//                               << sign << " " << shift << " ";
                    
                    if ( riv[d] / rr[d] == iv[d] ) {
                        /* interpolate
                         */
                        double value = -1.0 / ( avg * cdx[d] * fdx[d] );
                        
                        std::size_t nn = indices.size();
                        
//                         std::cout << riv << " " << shift << std::endl; std::cin.get();
                        
                        interface_mp->fine(riv, indices, values, d, shift, crse_fine_ba,
                                           mglevel_m[level+1].get());
                        
                        for (std::size_t iter = nn; iter < indices.size(); ++iter) {
                            values[iter] *= value;
                        }
                        
                    } else {
                        
//                         std::cout << sign * shift << " not interp ";
                        
                        double value = 1.0 / ( avg * cdx[d] * fdx[d] );
                        
                        indices.push_back( mglevel_m[level+1]->serialize(riv) );
                        values.push_back( value );
                    }
                    
//                     std::cout << std::endl; std::cin.get();
                    
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    };
    
    for (map_t::iterator it = cells.begin(); it != cells.end(); ++it) {
        // not covered coarse cell at crse-fine interface
        AmrIntVect_t iv = it->first;
        int globalidx = mglevel_m[level]->serialize(iv);
        
        indices.clear();
        values.clear();
        
        while ( !it->second.empty() ) {
            /*
             * "shift" is the amount of is a coarse cell that got refined
             * "d" is the direction to shift
             * 
             * --> check all covered neighbour cells
             */
            int shift = it->second.front().first;
            int d     = it->second.front().second;
            it->second.pop_front();
            
            
            /* we need to iterate over correct fine cells. It depends
             * on the orientation of the interface
             */
            int begin[BL_SPACEDIM] = { D_DECL( int(d == 0), int(d == 1), int(d == 2) ) };
            int end[BL_SPACEDIM]   = { D_DECL( int(d != 0), int(d != 1), int(d != 2) ) };
            
            // neighbour
            AmrIntVect_t covered = iv;
            covered[d] += shift;
                    
            /*
             * neighbour cell got refined but is not on physical boundary
             * --> we are at a crse-fine interface
             * 
             * we need now to find out which fine cells
             * are required to satisfy the flux matching
             * condition
             */
            
            switch ( shift ) {
                case -1:
                {
                    // --> interface is on the lower face
                    int ii = iv[0] * rr[0];
                    int jj = iv[1] * rr[1];
#if BL_SPACEDIM == 3
                    int kk = iv[2] * rr[2];
#endif
                    // iterate over all fine cells at the interface
                    // start with lower cells --> cover coarse neighbour
                    // cell
                    fill(indices, values, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
                    break;
                }
                case 1:
                default:
                {
                    // --> interface is on the upper face
                    int ii = covered[0] * rr[0];
                    int jj = covered[1] * rr[1];
#if BL_SPACEDIM == 3
                    int kk = covered[2] * rr[2];
#endif
                    fill(indices, values, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
                    break;
                }
            }
        }
        
        this->unique_m(indices, values);
        
//         std::cout << globalidx << " " << nNeighbours << " ";
//         for (uint ii = 0; ii < indices.size(); ++ii)
//             std::cout << indices[ii] << " ";
//         std::cout << std::endl; std::cin.get();
        
        int error = mglevel_m[level]->Bfine_p->InsertGlobalValues(globalidx,
                                                                  indices.size(),
                                                                  &values[0],
                                                                  &indices[0]);
        
        if ( error != 0 ) {
            // if e.g. nNeighbours < indices.size() --> error
            throw std::runtime_error("Error in filling the fine boundary matrix for level " +
                                     std::to_string(level) + "!");
        }
    }
    
    int error = mglevel_m[level]->Bfine_p->FillComplete(ColMap, RowMap, true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the fine boundary matrix for level " +
                                 std::to_string(level) + "!");
}


inline void AmrMultiGrid::buildDensityVector_m(const AmrField_t& rho, int level) {
    this->amrex2trilinos_m(rho, mglevel_m[level]->rho_p, level);
}


inline void AmrMultiGrid::buildPotentialVector_m(const AmrField_t& phi, int level) {
    this->amrex2trilinos_m(phi, mglevel_m[level]->phi_p, level);
}


void AmrMultiGrid::buildSmootherMatrix_m(int level) {
    // the base level does not require a smoother matrix
    if ( level == lbase_m )
        return;
    
    const Epetra_Map& RowMap = *mglevel_m[level]->map_p;
    
    mglevel_m[level]->S_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy,
                                                       RowMap, 1) );
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    double h2 = std::max(dx[0], dx[1]);
#if BL_SPACEDIM == 3
    h2 = std::max(h2, dx[2]);
#endif
    h2 *= h2;
    
    for (amrex::MFIter mfi(*mglevel_m[level]->mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    int globalidx = mglevel_m[level]->serialize(iv);
                    
                    /*
                     * check all directions (Laplacian stencil --> cross)
                     */
                    bool interior = true;
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
                        for (int shift = -1; shift <= 1; shift += 2) {
                            AmrIntVect_t biv = iv;                        
                            biv[d] += shift;
                            
                            switch ( mfab(biv) )
                            {
                                case AmrMultiGridLevel_t::Mask::COVERED:
                                    // covered --> interior cell
                                case AmrMultiGridLevel_t::Mask::INTERIOR:
                                    // interior cells
                                    interior *= true;
                                    break;
                                case AmrMultiGridLevel_t::Mask::BNDRY:
                                    // boundary cells
                                case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
                                    // boundary cells
                                    interior *= false;
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                    
                    double value = h2 * ( (interior) ? 0.25 : 0.1875 );
                    
                    int error = mglevel_m[level]->S_p->InsertGlobalValues(globalidx,
                                                                          1,
                                                                          &value,
                                                                          &globalidx);
                    
                    if ( error != 0 ) {
                        throw std::runtime_error("Error in filling the smoother matrix for level " +
                                                  std::to_string(level) + "!");
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = mglevel_m[level]->S_p->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the smoother matrix for level " +
                                 std::to_string(level) + "!");
    
}


void AmrMultiGrid::buildGradientMatrix_m(int level) {
    
    const Epetra_Map& RowMap = *mglevel_m[level]->map_p;
    
    int nEntries = 2;
    
    for (int d = 0; d < BL_SPACEDIM; ++d) {
        mglevel_m[level]->G_p[d] = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy,
                                                              RowMap, nEntries) );
    }
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    indices_t indices;
    coefficients_t values;
    
    auto check = [&](const AmrIntVect_t& iv,
                     const amrex::BaseFab<int>& mfab,
                     int dir,
                     double shift)
    {
        switch ( mfab(iv) )
        {
            case AmrMultiGridLevel_t::Mask::COVERED:
                // covered --> interior cell
            case AmrMultiGridLevel_t::Mask::INTERIOR:
                // interior cells
                
                indices.push_back( mglevel_m[level]->serialize(iv) );
                values.push_back(  shift * 0.5 / dx[dir] );
                
                break;
            case AmrMultiGridLevel_t::Mask::BNDRY:
            {
                // interior boundary cells --> only level > 0
                
                double value = shift * 0.5 / dx[dir];
                
                // use 1st order Lagrange --> only cells of this level required
                AmrIntVect_t tmp = iv;
                // first fine cell on refined coarse cell (closer to interface)
                tmp[dir] += shift;
                indices.push_back( mglevel_m[level]->serialize(tmp) );
                values.push_back( 2.0 * value );
                
                // second fine cell on refined coarse cell (further away from interface)
                tmp[dir] += shift;
                indices.push_back( mglevel_m[level]->serialize(tmp) );
                values.push_back( -1.0 * value );
                
                break;
            }
            case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
            {
                // physical boundary cells
                
                double value = shift * 0.5 / dx[dir];
                
                mglevel_m[level]->applyBoundary(iv,
                                                indices,
                                                values,
                                                value);
                break;
            }
            default:
                break;
        }
    };
    
    for (amrex::MFIter mfi(*mglevel_m[level]->mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    int globalidx = mglevel_m[level]->serialize(iv);
                    
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
                        
                        indices.clear();
                        values.clear();
                        
                        for (int shift = -1; shift <= 1; shift += 2) {
                            AmrIntVect_t niv = iv;                        
                            niv[d] += shift;
                            check(niv, mfab, d, shift);
                        }
                        
                        int error = mglevel_m[level]->G_p[d]->InsertGlobalValues(globalidx,
                                                                                 indices.size(),
                                                                                 &values[0],
                                                                                 &indices[0]);
                        
                        if ( error != 0 ) {
                            throw std::runtime_error("Error in filling the gradient matrix for level " +
                                                     std::to_string(level) + "!");
                        }
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    
    int error = 0;
    for (int d = 0; d < BL_SPACEDIM; ++d)
        error += mglevel_m[level]->G_p[d]->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the gradient matrix for level " +
                                 std::to_string(level) + "!");
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
                                    const Teuchos::RCP<vector_t>& mv) // FIXME component for electric field fab(iv, 0 / 1 / 2)
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


void AmrMultiGrid::unique_m(indices_t& indices,
                            coefficients_t& values)
{
    indices_t uindices;
    coefficients_t uvalues;
    
    // we need to sort for std::unique_copy
    // 20. Sept. 2017,
    // https://stackoverflow.com/questions/34878329/how-to-sort-two-vectors-simultaneously-in-c-without-using-boost-or-creating-te
    std::vector< std::pair<int, double> > pair;
    for (uint i = 0; i < indices.size(); ++i)
        pair.push_back( { indices[i], values[i] } );
    
    std::sort(pair.begin(), pair.end(),
              [](const std::pair<int, double>& a,
                  const std::pair<int, double>& b) {
                  return a.first < b.first;
              });
    
    for (uint i = 0; i < indices.size(); ++i) {
        indices[i] = pair[i].first;
        values[i]  = pair[i].second;
    }
    
    // requirement: duplicates are consecutive
    std::unique_copy(indices.begin(), indices.end(), std::back_inserter(uindices));
    
    uvalues.resize(uindices.size());
    
    for (std::size_t i = 0; i < uvalues.size(); ++i) {
        for (std::size_t j = 0; j < values.size(); ++j) {
            if ( uindices[i] == indices[j] )
                uvalues[i] += values[j];
        }
    }
    
    std::swap(values, uvalues);
    std::swap(indices, uindices);
}


void AmrMultiGrid::checkCrseBoundary_m(Teuchos::RCP<matrix_t>& B,
                                       int level,
                                       const amrex::BaseFab<int>& mfab,
                                       const AmrIntVect_t& lo,
                                       const AmrIntVect_t& hi,
                                       const amrex::BoxArray& ba)
{
    /* fine (this): level
     * coarse: level - 1
     */
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    indices_t indices;
    coefficients_t values;
    
    // helper function for boundary
    auto check = [&](const amrex::BaseFab<int>& mfab,
                     const AmrIntVect_t& iv,
                     indices_t& indices,
                     coefficients_t& values,
                     int dir, int shift)
    {
//         std::cout << iv << " " << iv[dir] << std::endl;
        bool top = (iv[(dir+1)%BL_SPACEDIM] % 2 == 1);  //FIXME Extend 3D
        
        switch ( mfab(iv) )
        {
            case AmrMultiGridLevel_t::Mask::COVERED:
                // covered (nothing to do here)
                break;
            case AmrMultiGridLevel_t::Mask::INTERIOR:
                // interior (nothing to do here)
                break;
            case AmrMultiGridLevel_t::Mask::BNDRY:
            {
                // internal boundary (only level > 0 have this kind of boundary)
                AmrIntVect_t civ = iv;
                civ.coarsen(AmrIntVect_t(D_DECL(2, 2, 2))); /*FIXME*/
                
//                 std::cout << civ << " " << iv << " " << top << std::endl; std::cin.get();
                
                std::size_t numEntries = indices.size();
                // we need boundary + indices from coarser level
                interface_mp->coarse(civ, indices, values, dir, shift, ba,
                                     top, mglevel_m[level-1].get());
                
                // we need normlization by mesh size squared
                for (std::size_t n = numEntries; n < indices.size(); ++n)
                    values[n] *= 1.0 / ( dx[dir] * dx[dir] );
                
                break;
            }
            case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
            {
                // physical boundary --> handled in Poisson matrix
                break;
            }
            default:
                throw std::runtime_error("Error in mask value");
                break;
        }
    };
    
    for (int i = lo[0]; i <= hi[0]; ++i) {
        for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
            for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                // last interior cell
                AmrIntVect_t iv(D_DECL(i, j, k));
                
                indices.clear();
                values.clear();
                
                // we need the global index of the interior cell (the center of the Laplacian stencil (i, j, k)
                int globalidx = mglevel_m[level]->serialize(iv);
                
                /*
                 * check all directions (Laplacian stencil --> cross)
                 */
                for (int d = 0; d < BL_SPACEDIM; ++d) {
                    for (int shift = -1; shift <= 1; shift += 2) {
                        AmrIntVect_t biv = iv;                        
                        biv[d] += shift;
                        check(mfab, biv,
                              indices,
                              values,
                              d,
                              shift);
                    }
                }
                
                /*
                 * In some cases indices occur several times, e.g. for corner points
                 * at the physical boundary --> sum them up
                 */
                this->unique_m(indices, values);
                
                int error = B->InsertGlobalValues(globalidx,
                                                  indices.size(),
                                                  &values[0],
                                                  &indices[0]);
                
                if ( error != 0 ) {
                    // if e.g. nNeighbours < indices.size() --> error
                    throw std::runtime_error("Error in filling the boundary matrix for level " +
                                             std::to_string(level) + "!");
                }
#if BL_SPACEDIM == 3
            }
#endif
        }
    }
}

void AmrMultiGrid::gsrb_level_m(Teuchos::RCP<vector_t>& e,
                                Teuchos::RCP<vector_t>& r,
                                int level)
{
    // apply "no fine" Laplacian
    Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
    tmp->PutScalar(0.0);
    
    // tmp = A * x
    mglevel_m[level]->A_p->Apply(/**mglevel_m[level]->error_p*/*e, *tmp);
    
    // tmp = tmp - r
    tmp->Update(-1.0, *r, 1.0);
    
    // tmp2 = S * tmp
    Teuchos::RCP<vector_t> tmp2 = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
    mglevel_m[level]->S_p->Apply(*tmp, *tmp2);
    
    e->Update(1.0, *tmp2, 1.0);
}


void AmrMultiGrid::restrict_m(int level) {
    
    Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(*mglevel_m[level+1]->map_p, false) );
    tmp->PutScalar(0.0);
    
//     std::cout << *mglevel_m[level+1]->residual_p << std::endl;
    
    this->residual_no_fine_m(tmp,
                             mglevel_m[level+1]->error_p,
                             mglevel_m[level]->error_p,
                             mglevel_m[level+1]->residual_p, level+1);
    
//     std::cout << *tmp << std::endl; std::cin.get();
    
    // average down: residual^(l-1) = R^(l) * tmp
//     mglevel_m[level]->residual_p->PutScalar(0.0);
    
//     std::cout << *mglevel_m[level]->residual_p << std::endl;
    
    mglevel_m[level+1]->R_p->Apply(*tmp, *mglevel_m[level]->residual_p);
    
//     std::cout << *mglevel_m[level]->residual_p << std::endl; //std::cin.get();
    
    // special matrix, i.e. matrix without covered cells
    // r^(l-1) = rho^(l-1) - A * phi^(l-1)
    
    
    vector_t fine2crse(mglevel_m[level]->A_p->OperatorDomainMap());
    fine2crse.PutScalar(0.0);
    
    // get boundary for 
//     if ( mglevel_m[level]->Bfine_p != Teuchos::null ) {
        mglevel_m[level]->Bfine_p->Apply(*mglevel_m[level+1]->phi_p, fine2crse);
//     }
    
    vector_t tmp2(mglevel_m[level]->As_p->OperatorDomainMap());
    mglevel_m[level]->As_p->Apply(*mglevel_m[level]->phi_p, tmp2);
    
    vector_t crse2fine(mglevel_m[level]->A_p->OperatorDomainMap());
    crse2fine.PutScalar(0.0);
    
    if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
        mglevel_m[level]->Bcrse_p->Apply(*mglevel_m[level-1]->phi_p, crse2fine);
    }
    
    tmp2.Update(1.0, fine2crse, 1.0, crse2fine, 1.0);
    
    Teuchos::RCP<vector_t> tmp3 = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
    
    mglevel_m[level]->UnCovered_p->Apply(*mglevel_m[level]->rho_p, *tmp3);
    
//     std::cout << *mglevel_m[level]->residual_p << std::endl;
    
    // ONLY subtract coarse rho
    mglevel_m[level]->residual_p->Update(1.0, *tmp3, -1.0, tmp2, 1.0);
    
//     std::cout << *mglevel_m[level]->residual_p << std::endl; std::cin.get();
}
