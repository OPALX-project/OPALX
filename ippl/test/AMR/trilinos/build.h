#include "tools.h"

#include "EpetraExt_RowMatrixOut.h"

void buildInterpolationMatrix(Teuchos::RCP<Epetra_CrsMatrix>& I,
                              const std::vector<Teuchos::RCP<Epetra_Map> >& maps,
                              const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                              const Array<Geometry>& geom,
                              const IntVect& rr, Epetra_MpiComm& comm, int level) {
    
    if ( level == (int)dmap.size() - 1 )
        return;
    
    int cnr[BL_SPACEDIM];
    int fnr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j) {
        cnr[j] = geom[level].Domain().length(j);
        fnr[j] = geom[level+1].Domain().length(j);
    }
    
    
    int nNeighbours = (2 << (BL_SPACEDIM -1 ));
    
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
#if BL_SPACEDIM == 3
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
#if BL_SPACEDIM == 3
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
                        const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                        const Array<Geometry>& geom,
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
    int nEntries = (BL_SPACEDIM << 1) + 1 /* plus boundaries */ + 10 /*FIXME*/;
    
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
    
    int nr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j)
        nr[j] = geom[level].Domain().length(j);
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
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
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
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
                                    double value = -1.0 / ( dx[d] * dx[d] );
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
                        values.push_back( 2.0 / ( dx[0] * dx[0] ) +
                                          2.0 / ( dx[1] * dx[1] )
#if BL_SPACEDIM == 3
                                          + 2.0 / ( dx[2] * dx[2] )
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
                
#if BL_SPACEDIM == 3
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
                            const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                            const Array<Geometry>& geom,
                            const IntVect& rr, Epetra_MpiComm& comm, int level) {
    if ( level == 0 )
        return;
    
    std::cout << "buildRestrictionMatrix" << std::endl;
    
    int cnr[BL_SPACEDIM];
    int fnr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j) {
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
    
    
#if BL_SPACEDIM == 3
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
#if BL_SPACEDIM == 3
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
#if BL_SPACEDIM == 3
                                for (int kref = 0; kref < rr[2]; ++kref) {
#endif
                                    IntVect riv(D_DECL(ii + iref, jj + jref, kk + kref));
//                                 riv.diagShift(1);
                                
                                    int fine_globidx = serialize(riv, &fnr[0]);
                                
                                    indices[numEntries] = fine_globidx;
                                    values[numEntries]  = val;
                                    ++numEntries;
#if BL_SPACEDIM == 3
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
#if BL_SPACEDIM == 3
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
#if BL_SPACEDIM == 3
    h2 = std::max(h2, dx[2]);
#endif
    h2 *= h2;
    
    std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > mask(
        new amrex::FabArray<amrex::BaseFab<int> >(grids, dmap, 1, 1)
    );
    
    mask->BuildMask(geom.Domain(), geom.periodicity(), -1, 1, 2, 0);
    
    int nr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j)
        nr[j] = geom.Domain().length(j);
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    IntVect iv(D_DECL(i, j, k));
                    
                    int globidx = serialize(iv, &nr[0]);
                    
                    /*
                     * check all directions (Laplacian stencil --> cross)
                     */
                    bool interior = true;
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
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
#if BL_SPACEDIM == 3
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
                               const Array<BoxArray>& grids, const Array<DistributionMapping>& dmap,
                               const Array<Geometry>& geom,
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
    int nEntries = (BL_SPACEDIM << 1) + 1 /* plus boundaries */ + 10 /*FIXME*/;
    
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
    
    int nr[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; ++j)
        nr[j] = geom[level].Domain().length(j);
    
    for (amrex::MFIter mfi(*mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
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
                        for (int d = 0; d < BL_SPACEDIM; ++d) {
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
                                            double value = -1.0 / ( dx[d] * dx[d] );
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
                                } else {
                                    /*
                                     * It has to be a refined cell!
                                     */
                                    
                                    /* we add 1 to the center iv of the Laplacian of
                                     * the special Poisson matrix. The fine contribution
                                     * is added to the boundary matrix Bfine at row globidx
                                     * (of iv)
                                     */
                                    indices.push_back( globidx );
                                    values.push_back( 1.0 / ( dx[d] * dx[d] ) );
                                }
                            }
                        }
                        
                        // check center
                        if ( mfab(iv) == 0 ) {
                            indices.push_back( globidx );
                            values.push_back( 2.0 / ( dx[0] * dx[0] ) +
                                              2.0 / ( dx[1] * dx[1] )
#if BL_SPACEDIM == 3
                                              + 2.0 / ( dx[2] * dx[2] )
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
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int error = A->FillComplete(true);
    
    if ( error != 0 )
        throw std::runtime_error("Error in completing Poisson matrix for level "
                                 + std::to_string(level) + "!");

//     std::cout << *A << std::endl;
    
    
    EpetraExt::RowMatrixToMatlabFile("poisson_matrix.txt", *A);
}
