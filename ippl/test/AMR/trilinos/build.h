#include "tools.h"

#include "EpetraExt_RowMatrixOut.h"

void buildInterpolationMatrix(Teuchos::RCP<Epetra_CrsMatrix>& I,
                              const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                              const Vector<BoxArray>& grids, const Vector<DistributionMapping>& dmap,
                              const Vector<Geometry>& geom,
                              const IntVect& rr, Epetra_MpiComm& comm, int level) {
    
    if ( level == (int)dmap.size() - 1 )
        return;
    
    int cnr[AMREX_SPACEDIM];
    int fnr[AMREX_SPACEDIM];
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        cnr[j] = geom[level].Domain().length(j);
        fnr[j] = geom[level+1].Domain().length(j);
    }
    
    
    int nNeighbours = (2 << (AMREX_SPACEDIM -1 ));
    
    std::vector<int> indices; //(nNeighbours);
    std::vector<double> values; //(nNeighbours);
    
    const Epetra_Map& RowMap = *maps[level+1];
    const Epetra_Map& ColMap = *maps[level];
    
    I = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                           RowMap, nNeighbours, false) );
    
    for (amrex::MFIter mfi(grids[level+1], dmap[level+1], false);
         mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx  = mfi.validbox();
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    int globidx = serialize(iv, &fnr[0]);
                    
                    int numEntries = 0;
                    indices.clear();
                    values.clear();
                    
                    /*
                     * we need boundary + indices from coarser level
                     */
                    stencil(iv, indices, values, numEntries, &cnr[0]);
                    
                    int error = I->InsertGlobalValues(globidx,
                                                      numEntries,
                                                      &values[0],
                                                      &indices[0]);
                    
                    if ( error != 0 ) {
                        // if e.g. nNeighbours < numEntries --> error
                        throw std::runtime_error("Error in filling the interpolation matrix for level " +
                                                 std::to_string(level) + "!");
                    }
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = I->FillComplete(ColMap, RowMap, true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the interpolation matrix for level " +
                                 std::to_string(level) + "!");
    
    std::cout << "Done." << std::endl;
    
    std::cout << *I << std::endl;
    
    
    EpetraExt::RowMatrixToMatlabFile("interpolation_matrix.txt", *I);
}


void buildPoissonMatrix(Teuchos::RCP<Epetra_CrsMatrix>& A,
                        const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                        const Vector<BoxArray>& grids, const Vector<DistributionMapping>& dmap,
                        const Vector<Geometry>& geom,
                        Epetra_MpiComm& comm, int level)
{
    /*
     * Laplacian of "no fine"
     */
    
    /*
     * 1D not supported
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    int nEntries = (AMREX_SPACEDIM << 1) + 1 /* plus boundaries */ + 10 /*FIXME*/;
    
    std::cout << "nEntries = " << nEntries << std::endl;
    
    const Epetra_Map& RowMap = *maps[level];
    
    A = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RowMap, nEntries) );
    
    std::vector<int> indices; //(nEntries);
    std::vector<double> values; //(nEntries);
    
    const double* dx = geom[level].CellSize();
    
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
        new amrex::FabArray<amrex::BaseFab<int> >(grids[level], dmap[level], 1, 1)
    );
    
    mask->BuildMask(geom[level].Domain(), geom[level].periodicity(), -1, 1, 2, 0);
    
    int nr[AMREX_SPACEDIM];
    for (int j = 0; j < AMREX_SPACEDIM; ++j)
        nr[j] = geom[level].Domain().length(j);
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    int numEntries = 0;
                    indices.clear();
                    values.clear();
                    IntVect iv(D_DECL(i, j, k));
                    int globidx = serialize(iv, &nr[0]);
                    
                    /*
                     * check neighbours in all directions (Laplacian stencil --> cross)
                     */
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        for (int shift = -1; shift <= 1; shift += 2) {
                            IntVect biv = iv;                        
                            biv[d] += shift;
                            
                            switch ( mfab(biv) )
                            {
                                case -1:
                                    // covered --> interior cell
                                case 0:
                                {
                                    indices.push_back( serialize(biv, &nr[0]) );
                                    values.push_back( 1.0 / ( dx[d] * dx[d] ) );
                                    break;
                                }
                                case 1:
                                    // boundary cell
                                    // only level > 0 have this kind of boundary
                                    
                                    /* Dirichlet boundary conditions from coarser level.
                                     * These are treated by the boundary matrix.
                                     */
                                    break;
                                case 2:
                                {
                                    // physical boundary cell
                                    double value = 1.0 / ( dx[d] * dx[d] );
                                    applyBoundary(biv,
                                                  indices,
                                                  values,
                                                  numEntries,
                                                  value,
                                                  &nr[0]);
                                    break;
                                }
                                default:
                                    break;
                            }
                        }
                    }
                    
                    // check center
                    if ( mfab(iv) == 0 ) {
                        indices.push_back( globidx );
                        values.push_back( -2.0 / ( dx[0] * dx[0] ) +
                                          -2.0 / ( dx[1] * dx[1] )
#if AMREX_SPACEDIM == 3
                                          - 2.0 / ( dx[2] * dx[2] )
#endif
                        );
                        ++numEntries;
                    }
                    
//                     for (uint ii = 0; ii < indices.size(); ++ii)
//                         std::cout << indices[ii] << " ";
//                     std::cout << std::endl;
//                     
//                     for (uint ii = 0; ii < values.size(); ++ii)
//                         std::cout << values[ii] << " ";
//                     std::cout << std::endl;
                    
                    unique(indices, values, numEntries);
                    
//                     for (uint ii = 0; ii < indices.size(); ++ii)
//                         std::cout << indices[ii] << " ";
//                     std::cout << std::endl;
//                     
//                     for (uint ii = 0; ii < values.size(); ++ii)
//                         std::cout << values[ii] << " ";
//                     std::cout << std::endl;
                    
                    
                    int error = A->InsertGlobalValues(globidx,
                                                      numEntries,
                                                      &values[0],
                                                      &indices[0]);
                    
                    if ( error != 0 )
                        throw std::runtime_error("Error in filling the Poisson matrix for level "
                                                 + std::to_string(level) + "!");
                
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = A->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing Poisson matrix for level "
                                 + std::to_string(level) + "!");

    std::cout << *A << std::endl;
    
    
    EpetraExt::RowMatrixToMatlabFile("poisson_matrix.txt", *A);
}


void buildVector(Teuchos::RCP<Epetra_Vector>& x, const Teuchos::RCP<Epetra_Map>& map, double value) {
    x = Teuchos::rcp( new Epetra_Vector(*map, false));
    x->PutScalar(value);
}


void buildRestrictionMatrix(Teuchos::RCP<Epetra_CrsMatrix>& R,
                            const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                            const Vector<BoxArray>& grids, const Vector<DistributionMapping>& dmap,
                            const Vector<Geometry>& geom,
                            const IntVect& rr, Epetra_MpiComm& comm, int level) {
    if ( level == 0 )
        return;
    
    std::cout << "buildRestrictionMatrix" << std::endl;
    
    int cnr[AMREX_SPACEDIM];
    int fnr[AMREX_SPACEDIM];
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        cnr[j] = geom[level-1].Domain().length(j);
        fnr[j] = geom[level].Domain().length(j);
    }
    
    
    /* Difficulty:  If a fine cell belongs to another processor than the underlying
     *              coarse cell, we get an error when filling the matrix since the
     *              cell (--> global index) does not belong to the same processor.
     * Solution:    Find all coarse cells that are covered by fine cells, thus,
     *              the distributionmap is correct.
     * 
     * 
     */
    amrex::BoxArray crse_fine_ba = grids[level-1];
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(grids[level], crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    Epetra_Map RowMap(*maps[level-1]);
    Epetra_Map ColMap(*maps[level]);
    
    std::cout << ColMap.IsOneToOne() << " " << RowMap.IsOneToOne() << std::endl;
    
    
#if AMREX_SPACEDIM == 3
    int nNeighbours = rr[0] * rr[1] * rr[2];
#else
    int nNeighbours = rr[0] * rr[1];
#endif
    
    std::vector<double> values(nNeighbours);
    std::vector<int> indices(nNeighbours);
    
    double val = 1.0 / double(nNeighbours);
    
    R = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RowMap, nNeighbours) );
    
    for (amrex::MFIter mfi(grids[level-1], dmap[level-1], false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            int ii = i * rr[0];
            for (int j = lo[1]; j <= hi[1]; ++j) {
                int jj = j * rr[1];
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    int kk = k * rr[2];
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    if ( crse_fine_ba.contains(iv) ) {
                        std::cout << iv << std::endl;
                        
                        
                    
//                     iv.diagShift(1);
                        int crse_globidx = serialize(iv, &cnr[0]);
                    
                        int numEntries = 0;
                    
                        // neighbours
                        for (int iref = 0; iref < rr[0]; ++iref) {
                            for (int jref = 0; jref < rr[1]; ++jref) {
#if AMREX_SPACEDIM == 3
                                for (int kref = 0; kref < rr[2]; ++kref) {
#endif
                                    IntVect riv(D_DECL(ii + iref, jj + jref, kk + kref));
//                                 riv.diagShift(1);
                                
                                    int fine_globidx = serialize(riv, &fnr[0]);
                                
                                    indices[numEntries] = fine_globidx;
                                    values[numEntries]  = val;
                                    ++numEntries;
#if AMREX_SPACEDIM == 3
                                }
#endif
                            }
                        }
                    
//                         std::cout << crse_globidx << " ";
//                         for (int i = 0; i < numEntries; ++i)
//                             std::cout << indices[i] << " ";
//                         std::cout << std::endl;
                        int error = R->InsertGlobalValues(crse_globidx, numEntries, &values[0], &indices[0]);
                        
                        if ( error != 0 ) {
                            throw std::runtime_error("Error in filling the restriction matrix for level " +
                            std::to_string(level) + "!");
                        }
                    }
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    std::cout << "Before fill complete." << std::endl;
    
    if ( R->FillComplete(ColMap, RowMap, true) != 0 )
        throw std::runtime_error("Error in completing the restriction matrix for level " +
                                 std::to_string(level) + "!");
    
    std::cout << "Done." << std::endl;
    
    std::cout << *R << std::endl;
    
    
    EpetraExt::RowMatrixToMatlabFile("restriction_matrix.txt", *R);
}


void buildSmootherMatrix(Teuchos::RCP<Epetra_CrsMatrix>& S,
                         const Teuchos::RCP<Epetra_Map>& map,
                         const BoxArray& grids, const DistributionMapping& dmap,
                         const Geometry& geom,
                         Epetra_MpiComm& comm, int level)
{    
    if ( level == 0 )
        return;
    
    double value = 0.0;
    
    const Epetra_Map& RowMap = *map;
    
    S = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                           RowMap, 1, false) );
    
    const double* dx = geom.CellSize();
    
    double h2 = std::max(dx[0], dx[1]);
#if AMREX_SPACEDIM == 3
    h2 = std::max(h2, dx[2]);
#endif
    h2 *= h2;
    
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
        new amrex::FabArray<amrex::BaseFab<int> >(grids, dmap, 1, 1)
    );
    
    mask->BuildMask(geom.Domain(), geom.periodicity(), -1, 1, 2, 0);
    
    int nr[AMREX_SPACEDIM];
    for (int j = 0; j < AMREX_SPACEDIM; ++j)
        nr[j] = geom.Domain().length(j);
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    int globidx = serialize(iv, &nr[0]);
                    
                    /*
                     * check all directions (Laplacian stencil --> cross)
                     */
                    bool interior = true;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        for (int shift = -1; shift <= 1; shift += 2) {
                            IntVect biv = iv;                        
                            biv[d] += shift;
                            
                            switch ( mfab(biv) )
                            {
                                case -1:
                                    // covered --> interior cell
                                case 0:
                                    // interior cells
                                    interior *= true;
                                    break;
                                case 1:
                                    // boundary cells
                                case 2:
                                    // boundary cells
                                    interior *= false;
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                    
                    value = h2 * ( (interior) ? 0.25 : 0.1875 );
                    int error = S->InsertGlobalValues(globidx,
                                                      1,
                                                      &value,
                                                      &globidx);
                    
                    if ( error != 0 ) {
                        // if e.g. nNeighbours < numEntries --> error
                        throw std::runtime_error("Error in filling the smoother matrix for level " +
                                                 std::to_string(level) + "!");
                    }
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = S->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the smoother matrix for level " +
                                 std::to_string(level) + "!");
    
    std::cout << "Done." << std::endl;
    
    std::cout << *S << std::endl;
    
    
    EpetraExt::RowMatrixToMatlabFile("smoother_matrix.txt", *S);
}


void buildSpecialPoissonMatrix(Teuchos::RCP<Epetra_CrsMatrix>& A,
                               const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                               const Vector<BoxArray>& grids, const Vector<DistributionMapping>& dmap,
                               const Vector<Geometry>& geom,
                               const IntVect& rr,
                               Epetra_MpiComm& comm, int level)
{
    /*
     * Laplacian of "with fine"
     */
    
    
    // find all cells with refinement
    amrex::BoxArray crse_fine_ba = grids[level];
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(grids[level+1], crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    /*
     * 1D not supported
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    int nEntries = (AMREX_SPACEDIM << 1) + 1 /* plus boundaries */ + 10 /*FIXME*/;
    
    std::cout << "nEntries = " << nEntries << std::endl;
    
    const Epetra_Map& RowMap = *maps[level];
    
    A = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, RowMap, nEntries) );
    
    std::vector<int> indices; //(nEntries);
    std::vector<double> values; //(nEntries);
    
    const double* dx = geom[level].CellSize();
    
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
        new amrex::FabArray<amrex::BaseFab<int> >(grids[level], dmap[level], 1, 1)
    );
    
    mask->BuildMask(geom[level].Domain(), geom[level].periodicity(), -1, 1, 2, 0);
    
    int nr[AMREX_SPACEDIM];
    for (int j = 0; j < AMREX_SPACEDIM; ++j)
        nr[j] = geom[level].Domain().length(j);
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    if ( !crse_fine_ba.contains(iv) ) {
                        int numEntries = 0;
                        indices.clear();
                        values.clear();
                        int globidx = serialize(iv, &nr[0]);
                        
                        /*
                         * check neighbours in all directions (Laplacian stencil --> cross)
                         */
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            for (int shift = -1; shift <= 1; shift += 2) {
                                IntVect biv = iv;                        
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
                                            indices.push_back( serialize(biv, &nr[0]) );
                                            values.push_back( 1.0 / ( dx[d] * dx[d] ) );
                                            
                                            // add center once
                                            indices.push_back( globidx );
                                            values.push_back( -1.0 / ( dx[d] * dx[d] ) );
                                            
                                            break;
                                        }
                                        case 1:
                                            // boundary cell
                                            // only level > 0 have this kind of boundary
                                            
                                            /* Dirichlet boundary conditions from coarser level.
                                             * These are treated by the boundary matrix.
                                             */
                                            break;
                                        case 2:
                                        {
                                            // physical boundary cell
                                            double value = 1.0 / ( dx[d] * dx[d] );
                                            applyBoundary(biv,
                                                          indices,
                                                          values,
                                                          numEntries,
                                                          value,
                                                          &nr[0]);
                                            
                                            
                                            // add center once
                                            indices.push_back( globidx );
                                            values.push_back( -1.0 / ( dx[d] * dx[d] ) );
                                            
                                            break;
                                        }
                                        default:
                                            break;
                                    }
                                } /*else {
                                    //
                                    // It has to be a refined cell!
                                    //
                                    
                                    // we add 1 to the center iv of the Laplacian of
                                    // the special Poisson matrix. The fine contribution
                                    // is added to the boundary matrix Bfine at row globidx
                                    // (of iv)
                                    //
                                    indices.push_back( globidx );
                                    values.push_back( 1.0 / ( dx[d] * dx[d] ) );
                                }*/
                            }
                        }
                        
//                         // check center
//                         if ( mfab(iv) == 0 ) {
//                             indices.push_back( globidx );
//                             values.push_back( -2.0 / ( dx[0] * dx[0] ) +
//                                               -2.0 / ( dx[1] * dx[1] )
// #if AMREX_SPACEDIM == 3
//                                               -2.0 / ( dx[2] * dx[2] )
// #endif
//                             );
//                             ++numEntries;
//                         }
                    
//                     for (uint ii = 0; ii < indices.size(); ++ii)
//                         std::cout << indices[ii] << " ";
//                     std::cout << std::endl;
//                     
//                     for (uint ii = 0; ii < values.size(); ++ii)
//                         std::cout << values[ii] << " ";
//                     std::cout << std::endl;
                    
                        unique(indices, values, numEntries);
                    
//                     for (uint ii = 0; ii < indices.size(); ++ii)
//                         std::cout << indices[ii] << " ";
//                     std::cout << std::endl;
//                     
//                     for (uint ii = 0; ii < values.size(); ++ii)
//                         std::cout << values[ii] << " ";
//                     std::cout << std::endl;
                    
                    
                        int error = A->InsertGlobalValues(globidx,
                                                          numEntries,
                                                          &values[0],
                                                          &indices[0]);
                    
                        if ( error != 0 )
                            throw std::runtime_error("Error in filling the Poisson matrix for level "
                                                     + std::to_string(level) + "!");
                    }
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = A->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing Poisson matrix for level "
                                 + std::to_string(level) + "!");

    std::cout << *A << std::endl;
    
    
    EpetraExt::RowMatrixToMatlabFile("poisson_matrix.txt", *A);
}


void checkBoundary(Teuchos::RCP<Epetra_CrsMatrix>& B,
                   const Vector<Geometry>& geom,
                   int level,
                   const amrex::BaseFab<int>& mfab,
                   const IntVect& lo,
                   const IntVect& hi)
{
    int nr[AMREX_SPACEDIM];
    for (int j = 0; j < AMREX_SPACEDIM; ++j)
        nr[j] = geom[level].Domain().length(j);
    
    int cnr[AMREX_SPACEDIM];
    if ( level > 0 ) {
        for (int j = 0; j < AMREX_SPACEDIM; ++j)
            cnr[j] = geom[level-1].Domain().length(j);
    }
        
    
    const double* dx = geom[level].CellSize();
    
    std::vector<int> indices; //(nNeighbours);
    std::vector<double> values; //(nNeighbours);
    
    // helper function for boundary
    auto check = [&](const amrex::BaseFab<int>& mfab,
                     const IntVect& iv,
                     int& numEntries,
                     std::vector<int>& indices,
                     std::vector<double>& values,
                     int dir)
    {
        switch ( mfab(iv) )
        {
            case -1:
                // covered (nothing to do here)
                break;
            case 0:
                // interior (nothing to do here)
                break;
            case 1:
            {
                // internal boundary (only level > 0 have this kind of boundary)
                
                int nn = numEntries;
                // we need boundary + indices from coarser level
                stencil(iv, indices, values, numEntries, &cnr[0]);
                
                // we need normlization by mesh size squared
                for (int n = nn; n < numEntries; ++n)
                    values[n] *= 1.0 / ( dx[dir] * dx[dir] );
                
                break;
            }
            case 2:
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
#if AMREX_SPACEDIM == 3
            for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                // last interior cell
                IntVect iv(D_DECL(i, j, k));
                
                int numEntries = 0;
                indices.clear();
                values.clear();
                
                // we need the global index of the interior cell (the center of the Laplacian stencil (i, j, k)
                int globidx = serialize(iv, &nr[0]);
                
                /*
                 * check all directions (Laplacian stencil --> cross)
                 */
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    for (int shift = -1; shift <= 1; shift += 2) {
                        IntVect biv = iv;                        
                        biv[d] += shift;
                        check(mfab, biv,
                              numEntries,
                              indices,
                              values,
                              d);
                    }
                }
                
//                 std::cout << iv << " " << globidx << std::endl;
//                 
//                 for (uint ii = 0; ii < indices.size(); ++ii)
//                     std::cout << indices[ii] << " ";
//                 std::cout << std::endl;
//                 
//                 for (uint ii = 0; ii < values.size(); ++ii)
//                     std::cout << values[ii] << " ";
//                 std::cout << std::endl;
                
                /*
                 * In some cases indices occur several times, e.g. for corner points
                 * at the physical boundary --> sum them up
                 */
                unique(indices, values, numEntries);
                
//                 for (uint ii = 0; ii < indices.size(); ++ii)
//                     std::cout << indices[ii] << " ";
//                 std::cout << std::endl;
//                 
//                 for (uint ii = 0; ii < values.size(); ++ii)
//                     std::cout << values[ii] << " ";
//                 std::cout << std::endl; std::cin.get();
                
                int error = B->InsertGlobalValues(globidx,
                                                  numEntries,
                                                  &values[0],
                                                  &indices[0]);
                
                if ( error != 0 ) {
                    // if e.g. nNeighbours < numEntries --> error
                    throw std::runtime_error("Error in filling the boundary matrix for level " +
                                             std::to_string(level) + "!");
                }
#if AMREX_SPACEDIM == 3
            }
#endif
        }
    }
}


void buildCrseBoundaryMatrix(Teuchos::RCP<Epetra_CrsMatrix>& Bcrse,
                             const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                             const Vector<BoxArray>& grids, const Vector<DistributionMapping>& dmap,
                             const Vector<Geometry>& geom,
                             const IntVect& rr, Epetra_MpiComm& comm, int level)
{
    if ( level == 0 )
        return;
    
    const Epetra_Map& RowMap = *maps[level];
    const Epetra_Map& ColMap = *maps[level-1];
    
    int nNeighbours = (2 << (AMREX_SPACEDIM -1 )) /*FIXME interpolation stencil indices*/ + 10;
    
    Bcrse = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                               RowMap, nNeighbours, false) );
    
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
        new amrex::FabArray<amrex::BaseFab<int> >(grids[level], dmap[level], 1, 1)
    );
    
    mask->BuildMask(geom[level].Domain(), geom[level].periodicity(), -1, 1, 2, 0);
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&    bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        /*
         * left boundary
         */
        IntVect lower(D_DECL(lo[0], lo[1], lo[2]));
        IntVect upper(D_DECL(lo[0], hi[1], hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
        /*
         * right boundary
         */
        lower = IntVect(D_DECL(hi[0], lo[1], lo[2]));
        upper = IntVect(D_DECL(hi[0], hi[1], hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
        
        /*
         * lower boundary (without left and right last cell!)
         */
        lower = IntVect(D_DECL(lo[0]+1, lo[1], lo[2]));
        upper = IntVect(D_DECL(hi[0]-1, lo[1], hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
        /*
         * upper boundary (without left and right last cell!)
         */
        lower = IntVect(D_DECL(lo[0]+1, hi[1], lo[2]));
        upper = IntVect(D_DECL(hi[0]-1, hi[1], hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
#if AMREX_SPACEDIM == 3
        /*
         * front boundary (without left, right, upper and lower last cell!)
         */
        lower = IntVect(D_DECL(lo[0]+1, lo[1]+1, lo[2]));
        upper = IntVect(D_DECL(hi[0]-1, hi[1]-1, lo[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
        
        /*
         * back boundary (without left, right, upper and lower last cell!)
         */
        lower = IntVect(D_DECL(lo[0]+1, lo[1]+1, hi[2]));
        upper = IntVect(D_DECL(hi[0]-1, hi[1]-1, hi[2]));
        
        checkBoundary(Bcrse, geom, level, mfab, lower, upper);
#endif
        
    }
    
    int error = Bcrse->FillComplete(ColMap, RowMap, true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the crse boundary matrix for level " +
                                 std::to_string(level) + "!");
    
    std::cout << "Done." << std::endl;
    
    std::cout << *Bcrse << std::endl;
    
    EpetraExt::RowMatrixToMatlabFile("crse_boundary_matrix.txt", *Bcrse);
}


void buildFineBoundaryMatrix(Teuchos::RCP<Epetra_CrsMatrix>& Bfine,
                             const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                             const Vector<BoxArray>& grids, const Vector<DistributionMapping>& dmap,
                             const Vector<Geometry>& geom,
                             const IntVect& rr, Epetra_MpiComm& comm, int level)
{
    if ( level == 1 /*finest*/ )
        return;
    
    // find all cells with refinement
    amrex::BoxArray crse_fine_ba = grids[level];
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(grids[level+1], crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    const Epetra_Map& RowMap = *maps[level];
    const Epetra_Map& ColMap = *maps[level+1];
    
    int nNeighbours = 4 /*#interfaces*/ * rr[0] * rr[1] /*of refined cell*/;
#if AMREX_SPACEDIM == 3
    nNeighbours = 6 /*#interfaces*/ * rr[0] * rr[1] * rr[2] /*of refined cell*/;
#endif
    
    Bfine = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
                                               RowMap, /*ColMap, */nNeighbours, false) );
    
    std::vector<int> indices; //(nEntries);
    std::vector<double> values; //(nEntries);
    
    
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
        new amrex::FabArray<amrex::BaseFab<int> >(grids[level], dmap[level], 1, 1)
    );
    
    mask->BuildMask(geom[level].Domain(), geom[level].periodicity(), -1, 1, 2, 0);
    
    int cnr[AMREX_SPACEDIM];
    int fnr[AMREX_SPACEDIM];
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        cnr[j] = geom[level].Domain().length(j);
        fnr[j] = geom[level+1].Domain().length(j);
    }
    
    const double* cdx = geom[level].CellSize();
    const double* fdx = geom[level+1].CellSize();
    
    /* Find all coarse cells that are at the crse-fine interace but are
     * not refined.
     * Put them into a map (to avoid duplicates, e.g. due to corners).
     * Finally, iterate through the list of cells
     */
    
    // std::list<std::pair<int, int> > contains the shift and direction to come to the covered cell
    std::map<IntVect, std::list<std::pair<int, int> > > cells;
    
    for (amrex::MFIter mfi(grids[level], dmap[level], false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    if ( !crse_fine_ba.contains(iv) && mfab(iv) != 2 /*not physical bc*/ ) {
                        /*
                         * iv is a coarse cell that got not refined
                         * 
                         * --> check all neighbours to see if at crse-fine
                         * interface
                         */
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            for (int shift = -1; shift <= 1; shift += 2) {
                                // neighbour
                                IntVect covered = iv;
                                covered[d] += shift;
                                
                                if ( crse_fine_ba.contains(covered) &&
                                     !isBoundary(covered, &cnr[0]) /*not physical bc*/ )
                                {
                                    // neighbour is covered by fine cells
                                    
                                    // insert if not yet exists
                                    cells[iv].push_back(std::make_pair(shift, d));
                                }
                            }
                        }
                    }
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    
    auto fill = [&](std::vector<int>& indices,
                         std::vector<double>& values,
                         int& numEntries,
                         D_DECL(int ii, int jj, int kk),
                         int* begin, int* end, int d, const IntVect& iv, int shift, int sign)
    {
        // number of fine cell gradients
        double avg = 1.0;
        switch ( d ) {
            case 0:
                // horizontal
                avg = rr[1];
#if AMREX_SPACEDIM == 3
                avg *= rr[2];
#endif
                break;
            case 1:
                // vertical
                avg = rr[0];
#if AMREX_SPACEDIM == 3
                avg *= rr[2];
#endif
                break;
#if AMREX_SPACEDIM == 3
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
#if AMREX_SPACEDIM == 3
                for (int kref = kk - begin[2]; kref <= kk + end[2]; ++kref) {
#endif
                    /* Since all fine cells on the not-refined cell are
                     * outside of the "domain" --> we need to interpolate
                     * using open boundary condition.
                     */
                    sign *= ( d == 2 ) ? -1.0 : 1.0;
                    
                    IntVect riv(D_DECL(iref, jref, kref));
                    
                    std::cout << "fine: " << riv << " sign = " << sign << std::endl;
                    
                    double value = double(shift) * sign / ( avg * cdx[d] * fdx[d] );
                    
                    if ( riv[d] / rr[d] == iv[d] ) {
                        /* this is a fine cell of the
                         * not refined cell, i.e. it's a
                         * ghost cell --> interpolate
                         */
                        
//                         std::size_t nn = indices.size();
                        
//                         stencil(riv, indices, values, numEntries, &cnr[0]);
                        
//                         for (std::size_t iter = nn; iter < indices.size(); ++iter) {
//                             values[iter] *= value;
//                         }
                        
                        /* interpolate
                         */
                        
                        // undo shift
                        riv[d] += shift;
                        
                        indices.push_back( serialize(riv, &fnr[0]) );
                        values.push_back( 2.0 * value );
                        ++numEntries;
                        
                        riv[d] += shift;
                        
                        indices.push_back( serialize(riv, &fnr[0]) );
                        values.push_back( - value );
//                         ++numEntries;
                    } else {
                        indices.push_back( serialize(riv, &fnr[0]) );
                        values.push_back( value );
                    }
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
    };
    
    for (std::map<IntVect, std::list<std::pair<int, int> > >::iterator it = cells.begin(); it != cells.end(); ++it) {
        
        // not covered coarse cell at crse-fine interface
        IntVect iv = it->first;
        int globidx = serialize(iv, &cnr[0]);
        
        int numEntries = 0;
        indices.clear();
        values.clear();
        
        while ( !it->second.empty() ) {
            /*
             * "shift" is the amount of IntVect to a coarse cell that got refined
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
            int begin[AMREX_SPACEDIM] = { D_DECL( int(d == 0), int(d == 1), int(d == 2) ) };
            int end[AMREX_SPACEDIM]   = { D_DECL( int(d != 0), int(d != 1), int(d != 2) ) };
            
            
            // neighbour
            IntVect covered = iv;
            covered[d] += shift;
            
            std::cout << iv << " " << covered << std::endl; std::cin.get();
                    
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
#if AMREX_SPACEDIM == 3
                    int kk = iv[2] * rr[2];
#endif
                    // iterate over all fine cells at the interface
                    // start with lower cells --> cover coarse neighbour
                    // cell
                    fill(indices, values, numEntries, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
                    break;
                }
                case 1:
                default:
                {
                    // --> interface is on the upper face
                    int ii = covered[0] * rr[0];
                    int jj = covered[1] * rr[1];
#if AMREX_SPACEDIM == 3
                    int kk = covered[2] * rr[2];
#endif
                    fill(indices, values, numEntries, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
                    break;
                }
            }
        }
        
        unique(indices, values, numEntries);
        
        std::cout << globidx << std::endl;
        for (uint ii = 0; ii < indices.size(); ++ii)
            std::cout << indices[ii] << " ";
        std::cout << std::endl;
        
        int error = Bfine->InsertGlobalValues(globidx,
                                              numEntries,
                                              &values[0],
                                              &indices[0]);
        
        if ( error != 0 ) {
            // if e.g. nNeighbours < numEntries --> error
            throw std::runtime_error("Error in filling the boundary matrix for level " +
                                     std::to_string(level) + "!");
        }
    }
    
    int error = Bfine->FillComplete(ColMap, RowMap, true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing the crse boundary matrix for level " +
                                 std::to_string(level) + "!");
    
    std::cout << "Done." << std::endl;
    
    std::cout << *Bfine << std::endl;
    
    EpetraExt::RowMatrixToMatlabFile("fine_boundary_matrix.txt", *Bfine);
}


// void buildFineBoundaryMatrix(Teuchos::RCP<Epetra_CrsMatrix>& Bfine,
//                              const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
//                              const Vector<BoxArray>& grids, const Vector<DistributionMapping>& dmap,
//                              const Vector<Geometry>& geom,
//                              const IntVect& rr, Epetra_MpiComm& comm, int level)
// {
//     if ( level == 1 /*finest*/ )
//         return;
//     
//     // find all cells with refinement
//     amrex::BoxArray crse_fine_ba = grids[level];
//     crse_fine_ba.refine(rr);
//     crse_fine_ba = amrex::intersect(grids[level+1], crse_fine_ba);
//     crse_fine_ba.coarsen(rr);
//     
//     const Epetra_Map& RowMap = *maps[level];
//     const Epetra_Map& ColMap = *maps[level+1];
//     
//     int nNeighbours = 4 /*#interfaces*/ * rr[0] * rr[1] /*of refined cell*/;
// #if AMREX_SPACEDIM == 3
//     nNeighbours = 6 /*#interfaces*/ * rr[0] * rr[1] * rr[2] /*of refined cell*/;
// #endif
//     
//     Bfine = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,
//                                                RowMap, /*ColMap, */nNeighbours, false) );
//     
//     std::vector<int> indices; //(nEntries);
//     std::vector<double> values; //(nEntries);
//     
//     
//     std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
//         new amrex::FabArray<amrex::BaseFab<int> >(grids[level], dmap[level], 1, 1)
//     );
//     
//     mask->BuildMask(geom[level].Domain(), geom[level].periodicity(), -1, 1, 2, 0);
//     
//     int cnr[AMREX_SPACEDIM];
//     int fnr[AMREX_SPACEDIM];
//     for (int j = 0; j < AMREX_SPACEDIM; ++j) {
//         cnr[j] = geom[level].Domain().length(j);
//         fnr[j] = geom[level+1].Domain().length(j);
//     }
//     
//     const double* cdx = geom[level].CellSize();
//     const double* fdx = geom[level+1].CellSize();
//     
//     /* Find all coarse cells that are at the crse-fine interace but are
//      * not refined.
//      * Put them into a map (to avoid duplicates, e.g. due to corners).
//      * Finally, iterate through the list of cells
//      */
//     
//     // std::list<std::pair<int, int> > contains the shift and direction to come to the covered cell
//     std::map<IntVect, std::list<std::pair<int, int> > > cells;
//     
//     for (amrex::MFIter mfi(grids[level], dmap[level], false); mfi.isValid(); ++mfi) {
//         const amrex::Box&          bx  = mfi.validbox();
//         const amrex::BaseFab<int>& mfab = (*mask)[mfi];
//             
//         const int* lo = bx.loVect();
//         const int* hi = bx.hiVect();
//         
//         for (int i = lo[0]; i <= hi[0]; ++i) {
//             for (int j = lo[1]; j <= hi[1]; ++j) {
// #if AMREX_SPACEDIM == 3
//                 for (int k = lo[2]; k <= hi[2]; ++k) {
// #endif
//                     IntVect iv(D_DECL(i, j, k));
//                     
//                     if ( !crse_fine_ba.contains(iv) && mfab(iv) != 2 /*not physical bc*/ ) {
//                         /*
//                          * iv is a coarse cell that got not refined
//                          * 
//                          * --> check all neighbours to see if at crse-fine
//                          * interface
//                          */
//                         for (int d = 0; d < AMREX_SPACEDIM; ++d) {
//                             for (int shift = -1; shift <= 1; shift += 2) {
//                                 // neighbour
//                                 IntVect covered = iv;
//                                 covered[d] += shift;
//                                 
//                                 if ( crse_fine_ba.contains(covered) &&
//                                      !isBoundary(covered, &cnr[0]) /*not physical bc*/ )
//                                 {
//                                     // neighbour is covered by fine cells
//                                     
//                                     // insert if not yet exists
//                                     cells[iv].push_back(std::make_pair(shift, d));
//                                 }
//                             }
//                         }
//                     }
// #if AMREX_SPACEDIM == 3
//                 }
// #endif
//             }
//         }
//     }
//     
//     
//     auto fill = [&](std::vector<int>& indices,
//                          std::vector<double>& values,
//                          int& numEntries,
//                          D_DECL(int ii, int jj, int kk),
//                          int* begin, int* end, int d, const IntVect& iv, int shift, int sign)
//     {
//         for (int iref = ii - begin[0]; iref <= ii + end[0]; ++iref) {
//             
//             sign *= ( d == 0 ) ? -1.0 : 1.0;
//             
//             for (int jref = jj - begin[1]; jref <= jj + end[1]; ++jref) {
//                 
//                 sign *= ( d == 1 ) ? -1.0 : 1.0;
// #if AMREX_SPACEDIM == 3
//                 for (int kref = kk - begin[2]; kref <= kk + end[2]; ++kref) {
// #endif
//                     /* Since all fine cells on the not-refined cell are
//                      * outside of the "domain" --> we need to interpolate
//                      * using open boundary condition.
//                      */
//                     sign *= ( d == 2 ) ? -1.0 : 1.0;
//                     
//                     IntVect riv(D_DECL(iref, jref, kref));
//                     
//                     if ( riv[d] / rr[d] == iv[d] ) {
//                         /* interpolate
//                          */
//                         
//                         // undo shift
//                         riv[d] += shift;
//                         
//                         double value = sign / ( cdx[d] * fdx[d] );
//                         
//                         indices.push_back( serialize(riv, &fnr[0]) );
//                         values.push_back( 2.0 * value );
//                         ++numEntries;
//                         
//                         riv[d] += shift;
//                         
//                         indices.push_back( serialize(riv, &fnr[0]) );
//                         values.push_back( - value );
//                         ++numEntries;
//                     } else {
//                         indices.push_back( serialize(riv, &fnr[0]) );
//                         values.push_back( sign / ( cdx[d] * fdx[d] ) );
//                     }
// #if AMREX_SPACEDIM == 3
//                 }
// #endif
//             }
//         }
//     };
//     
//     for (std::map<IntVect, std::list<std::pair<int, int> > >::iterator it = cells.begin(); it != cells.end(); ++it) {
//         
//         // not covered coarse cell at crse-fine interface
//         IntVect iv = it->first;
//         int globidx = serialize(iv, &cnr[0]);
//         
//         int numEntries = 0;
//         indices.clear();
//         values.clear();
//         
//         while ( !it->second.empty() ) {
//             /*
//              * "shift" is the amount of IntVect to a coarse cell that got refined
//              * "d" is the direction to shift
//              * 
//              * --> check all covered neighbour cells
//              */
//             int shift = it->second.front().first;
//             int d     = it->second.front().second;
//             it->second.pop_front();
//             
//             
//             /* we need to iterate over correct fine cells. It depends
//              * on the orientation of the interface
//              */
//             int begin[AMREX_SPACEDIM] = { D_DECL( int(d == 0), int(d == 1), int(d == 2) ) };
//             int end[AMREX_SPACEDIM]   = { D_DECL( int(d != 0), int(d != 1), int(d != 2) ) };
//             
//             
//             // neighbour
//             IntVect covered = iv;
//             covered[d] += shift;
//                     
//             /*
//              * neighbour cell got refined but is not on physical boundary
//              * --> we are at a crse-fine interface
//              * 
//              * we need now to find out which fine cells
//              * are required to satisfy the flux matching
//              * condition
//              */
//             
//             switch ( shift ) {
//                 case -1:
//                 {
//                     // --> interface is on the lower face
//                     int ii = iv[0] * rr[0];
//                     int jj = iv[1] * rr[1];
// #if AMREX_SPACEDIM == 3
//                     int kk = iv[2] * rr[2];
// #endif
//                     // iterate over all fine cells at the interface
//                     // start with lower cells --> cover coarse neighbour
//                     // cell
//                     fill(indices, values, numEntries, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
//                     break;
//                 }
//                 case 1:
//                 default:
//                 {
//                     // --> interface is on the upper face
//                     int ii = covered[0] * rr[0];
//                     int jj = covered[1] * rr[1];
// #if AMREX_SPACEDIM == 3
//                     int kk = covered[2] * rr[2];
// #endif
//                     fill(indices, values, numEntries, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, -1.0);
//                     break;
//                 }
//             }
//         }
//         
//         unique(indices, values, numEntries);
//         
//         std::cout << globidx << std::endl;
//         for (uint ii = 0; ii < indices.size(); ++ii)
//             std::cout << indices[ii] << " ";
//         std::cout << std::endl;
//         
//         int error = Bfine->InsertGlobalValues(globidx,
//                                               numEntries,
//                                               &values[0],
//                                               &indices[0]);
//         
//         if ( error != 0 ) {
//             // if e.g. nNeighbours < numEntries --> error
//             throw std::runtime_error("Error in filling the boundary matrix for level " +
//                                      std::to_string(level) + "!");
//         }
//     }
//     
//     int error = Bfine->FillComplete(ColMap, RowMap, true);
//     
//     if ( error != 0 )
//         throw std::runtime_error("Error in completing the crse boundary matrix for level " +
//                                  std::to_string(level) + "!");
//     
//     std::cout << "Done." << std::endl;
//     
//     std::cout << *Bfine << std::endl;
//     
//     EpetraExt::RowMatrixToMatlabFile("fine_boundary_matrix.txt", *Bfine);
// }
