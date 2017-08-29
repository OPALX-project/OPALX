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
        this->buildBoundaryMatrix_m(lev);
        
        this->buildPoissonMatrix_m(lev);
        
        int ilev = lbase + lev;
        
        this->buildDensityVector_m(*rho[ilev], lev);
        
        this->buildPotentialVector_m(*phi[ilev], lev);
    }
    
    
    double err = 1.0e7;
    double eps = 1.0e-8;
    
    // initial error
    double residualsum = 0.0;
    double rhosum = 0.0;
    
    for (int lev = 0; lev < nLevel; ++lev) {
        
        this->residual_m(mglevel_m[lev]->residual_p,
                         mglevel_m[lev]->rho_p,
                         mglevel_m[lev]->A_p,
                         mglevel_m[lev]->phi_p);
        
        double tmp = 0;
        mglevel_m[lev]->residual_p->Norm2(&tmp);
        residualsum += tmp;
        
        tmp = 0;
        mglevel_m[lev]->rho_p->Norm2(&tmp);
        rhosum += tmp;
    }
    
    
    while ( residualsum > eps * rhosum ) {
        
        std::cout << residualsum << " " << eps * rhosum << std::endl; //std::cin.get();
        
        relax_m(lfine);
        
        // update residual
        for (int lev = 0; lev < nLevel; ++lev) {
            this->residual_m(mglevel_m[lev]->residual_p,
                             mglevel_m[lev]->rho_p,
                             mglevel_m[lev]->A_p,
                             mglevel_m[lev]->phi_p);
        }
        
        residualsum = l2error_m();
    }
//     
//     for (int lev = nLevel-1; lev > -1; --lev) {
//         int ilev = lbase + lev;
//         mglevel_m[lev]->phi->FillBoundary(mglevel_m[lev]->geom.periodicity());
//         this->gradient_m(lev, *efield[ilev]);
//     }
//     
    std::cout << residualsum << " " << eps * rhosum << std::endl;
    
    // copy solution back
    for (int lev = 0; lev < nLevel; ++lev) {
        int ilev = lbase + lev;
        
        this->trilinos2amrex_m(*phi[ilev], mglevel_m[lev]->phi_p);
    }
}


void AmrMultiGrid::residual_m(Teuchos::RCP<vector_t>& r,
                              const Teuchos::RCP<vector_t>& b,
                              const Teuchos::RCP<matrix_t>& A,
                              const Teuchos::RCP<vector_t>& x,
                              bool temporary)
{
    /*
     * r = b - A*x
     */
    
    
    if ( temporary ) {
        vector_t tmp(A->OperatorDomainMap());
    
        // tmp = A * x
        A->Apply(*x, tmp);
    
        // r = b - tmp
        r->Update(1.0, *b, -1.0, tmp, 0.0);
        
    } else {
        // r = A * x
        A->Apply(*x, *r);
        
        // r = b - r;
        r->Update(1.0, *b, -1.0);
    }
}


void AmrMultiGrid::relax_m(int level) {
    
    if ( level == lfine_m ) {
        //TODO residual is already computed at beginning
        this->residual_m(mglevel_m[level]->residual_p,
                         mglevel_m[level]->rho_p,
                         mglevel_m[level]->A_p,
                         mglevel_m[level]->phi_p);
    }
    
    if ( level > 0 ) {
        //TODO
        
        // phi^(l, save) = phi^(l)
        vector_t phi_save = *mglevel_m[level]->phi_p;
        
        mglevel_m[level-1]->error_p->PutScalar(0.0);
        
        // smoothing
        // ...
        
        // phi = phi + e
        mglevel_m[level]->phi_p->Update(1.0, *mglevel_m[level]->error_p, 1.0);
        
        /*
         * restrict
         */
        // r^(l-1) = rho^(l-1) - A * phi^(l-1)
        this->residual_m(mglevel_m[level-1]->residual_p,
                         mglevel_m[level-1]->rho_p,
                         mglevel_m[level-1]->A_p,
                         mglevel_m[level-1]->phi_p);
        
        Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
        this->residual_m(tmp,
                         mglevel_m[level]->residual_p,
                         mglevel_m[level]->A_p,
                         mglevel_m[level]->error_p);
        
        // average down: residual^(l-1) = R^(l) * tmp
        mglevel_m[level]->R_p->Apply(*tmp, *mglevel_m[level]->residual_p);
        
        this->relax_m(level - 1);
        
        /*
         * prolongate / interpolate
         */
        
        // interpolate
        tmp = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
        mglevel_m[level]->I_p->Apply(*mglevel_m[level-1]->error_p, *tmp);
        
        // e^(l) += tmp
        mglevel_m[level]->error_p->Update(1.0, *tmp, 1.0);
        
        this->residual_m(mglevel_m[level]->residual_p,
                         mglevel_m[level]->residual_p,
                         mglevel_m[level]->A_p,
                         mglevel_m[level]->error_p,
                         true);
        
        // delta error
        tmp = Teuchos::rcp( new vector_t(*mglevel_m[level]->map_p, false) );
        
        // smoothing
        // ...
        
        // e^(l) += de^(l)
        mglevel_m[level]->error_p->Update(1.0, *tmp, 1.0);
        
        // phi^(l) = phi^(l, save) + e^(l)
        mglevel_m[level]->error_p->Update(1.0, phi_save, 1.0, *mglevel_m[level]->error_p, 1.0);
        
    } else {
        // e = A^(-1)r
        solver_m.solve(mglevel_m[level]->A_p,
                       mglevel_m[level]->error_p,
                       mglevel_m[level]->residual_p);
        
        // phi = phi + e
        mglevel_m[level]->phi_p->Update(1.0, *mglevel_m[level]->error_p, 1.0);
    }
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
