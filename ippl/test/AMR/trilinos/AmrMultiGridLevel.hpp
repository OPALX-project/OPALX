template <class MatrixType, class VectorType>
AmrMultiGridLevel<MatrixType,
                  VectorType>::AmrMultiGridLevel(const AmrField_u& _rho,
                                                 const AmrField_u& _phi,
                                                 const AmrGeometry_t& _geom,
                                                 Epetra_MpiComm& comm)
    : grids(_rho->boxArray()),
      dmap(_rho->DistributionMap()),
      geom(_geom),
      A_p(Teuchos::null),
      rho_p(Teuchos::null),
      phi_p(Teuchos::null),
      r_p(Teuchos::null),
      err_p(Teuchos::null),
      error(new AmrField_t(_rho->boxArray(), _rho->DistributionMap(), 1, 1)),
      rho(_rho.get()),
      phi(_phi.get()),
      residual(new AmrField_t(_rho->boxArray(), _rho->DistributionMap(), 1, 1)),
      phi_saved(new AmrField_t(_rho->boxArray(), _rho->DistributionMap(), 1, 1)),
      delta_err(new AmrField_t(_rho->boxArray(), _rho->DistributionMap(), 1, 1))
{
    for (int j = 0; j < 3; ++j)
        nr_m[j] = _geom.Domain().length(j);
    
    this->buildLevelMask_m();
    
    this->buildMap(_phi, comm);
    this->buildMatrix_m();
}


template <class MatrixType, class VectorType>
AmrMultiGridLevel<MatrixType, VectorType>::~AmrMultiGridLevel()
{
    A_p = Teuchos::null;
    rho_p = Teuchos::null;
    phi_p = Teuchos::null;
    r_p = Teuchos::null;
    err_p = Teuchos::null;
    map_mp = Teuchos::null;
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::copyTo(Teuchos::RCP<vector_t>& mv, const AmrField_t& mf) {
    int localNumElements = 0;
    std::vector<double> values;
    std::vector<int> globalindices;
    
    const double* dx = geom.CellSize();
    
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::FArrayBox&    fab = mf[mfi];
        const amrex::BaseFab<int>& mfab = (*masks_m)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    AmrIntVect_t iv(i, j, k);
                    
                    int globalidx = serialize_m(iv);
                    
                    globalindices.push_back(globalidx);
                    
                    values.push_back(fab(iv));
                    
                    ++localNumElements;
                }
            }
        }
    }
    
    if ( mv.is_null() )
        mv = Teuchos::rcp( new vector_t(*map_mp, false) );
    
    int success = mv->ReplaceGlobalValues(localNumElements,
                                             &values[0],
                                             &globalindices[0]);
    
    if ( success == 1 )
        throw std::runtime_error("Error in filling the vector!");
}


template <class MatrixType, class VectorType>
int AmrMultiGridLevel<MatrixType, VectorType>::serialize_m(const AmrIntVect_t& iv) const {
    return iv[0] + (iv[1] + nr_m[1] * iv[2]) * nr_m[0];
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::copyBack(AmrField_t& mf,
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
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    AmrIntVect_t iv(i, j, k);
                    fab(iv) = (*mv)[localidx++];
                }
            }
        }
    }
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::save() {
    AmrField_t::Copy(*phi_saved, *phi, 0, 0, 1, 1);
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::buildLevelMask_m() {
    masks_m.reset(new amrex::FabArray<amrex::BaseFab<int> >(grids, dmap, 1, 1));
    masks_m->BuildMask(geom.Domain(), geom.periodicity(),
                       Mask::COVERED, Mask::BNDRY,
                       Mask::PHYSBNDRY, Mask::INTERIOR);
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::buildMap(const AmrField_u& phi,
                                                         Epetra_MpiComm& comm)
{
#ifdef DEBUG
    if ( Ippl::myNode() == 0 ) {
        std::cout << "nx = " << nr_m[0] << std::endl
                  << "ny = " << nr_m[1] << std::endl
                  << "nz = " << nr_m[2] << std::endl;
    }
#endif

    int localNumElements = 0;
    std::vector<double> values;
    std::vector<int> globalindices;
    
    const double* dx = geom.CellSize();
    
    for (amrex::MFIter mfi(*phi, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::FArrayBox&    fab = (*phi)[mfi];
        const amrex::BaseFab<int>& mfab = (*masks_m)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    AmrIntVect_t iv(i, j, k);
                    
                    int globalidx = serialize_m(iv);
                    
                    globalindices.push_back(globalidx);
                    
                    values.push_back(fab(iv));
                    
                    ++localNumElements;
                }
            }
        }
    }
    
    // compute map based on localelements
    // create map that specifies which processor gets which data
    const int baseIndex = 0;    // where to start indexing
        
    // numGlobalElements == N
    const amrex::BoxArray& grids = phi->boxArray();
    int N = grids.numPts();
    
#ifdef DEBUG
    if ( Ippl::myNode() == 0 )
        std::cout << "N = " << N << std::endl;
#endif
    
    map_mp = Teuchos::rcp( new Epetra_Map(N, localNumElements,
                                         &globalindices[0], baseIndex, comm) );
    
    phi_p = Teuchos::rcp( new vector_t(*map_mp, false) );
    
    int success = phi_p->ReplaceGlobalValues(localNumElements,
                                             &values[0],
                                             &globalindices[0]);
    
    if ( success == 1 )
        throw std::runtime_error("Error in filling the vector!");
    
    
    r_p = Teuchos::rcp( new vector_t(*map_mp, false) );
//     r_p->PutScalar(0.0);
    err_p = Teuchos::rcp( new vector_t(*map_mp, false) );
    delta_err_p = Teuchos::rcp( new vector_t(*map_mp, false) );
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::buildMatrix_m() {
    // 3D --> 7 elements per row
    A_p = Teuchos::rcp( new matrix_t(Epetra_DataAccess::Copy, *map_mp, 7) );
    
    int indices[7] = {0, 0, 0, 0, 0, 0, 0};
    double values[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    const double* dx = geom.CellSize();
    
    for (amrex::MFIter mfi(*masks_m, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*masks_m)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    
                    
//                     double subx = 0.0;
//                     double suby = 0.0;
//                     double subz = 0.0;
                    
                    int numEntries = 0;
                    AmrIntVect_t iv(i, j, k);
                    int globalRow = serialize_m(iv);
                    
                    // check left neighbor in x-direction
                    AmrIntVect_t xl(i-1, j, k);
                    if ( mfab(xl) < 1 ) {
                        int gidx = serialize_m(xl/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
                        ++numEntries;
                    } //else if ( mfab(xl) == 2 )
                      //  subx = 1.0;
                    
                    
                    // check right neighbor in x-direction
                    AmrIntVect_t xr(i+1, j, k);
                    if ( mfab(xr) < 1 ) {
                        int gidx = serialize_m(xr/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[0] * dx[0] );
                        ++numEntries;
                    } //else if ( mfab(xr) == 2 )
                      //  subx += 1.0;
                    
                    // check lower neighbor in y-direction
                    AmrIntVect_t yl(i, j-1, k);
                    if ( mfab(yl) < 1 ) {
                        int gidx = serialize_m(yl/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
                        ++numEntries;
                    } //else if ( mfab(yl) == 2 )
                      //  suby += 1.0;
                    
                    // check upper neighbor in y-direction
                    AmrIntVect_t yr(i, j+1, k);
                    if ( mfab(yr) < 1 ) {
                        int gidx = serialize_m(yr/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[1] * dx[1] );
                        ++numEntries;
                    } // else if ( mfab(yr) == 2 )
                      //  suby += 1.0;
                    
                    // check front neighbor in z-direction
                    AmrIntVect_t zl(i, j, k-1);
                    if ( mfab(zl) < 1 ) {
                        int gidx = serialize_m(zl/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
                        ++numEntries;
                    } // else if ( mfab(zl) == 2 )
                      //  subz += 1.0;
                    
                    // check back neighbor in z-direction
                    AmrIntVect_t zr(i, j, k+1);
                    if ( mfab(zr) < 1 ) {
                        int gidx = serialize_m(zr/*shift*/);
                        indices[numEntries] = gidx;
                        values[numEntries]  = -1.0 / ( dx[2] * dx[2] );
                        ++numEntries;
                    } // else if ( mfab(zr) == 2 )
                      //  subz += 1.0;
                    
                    // check center
                    if ( mfab(iv) == 0 ) {
                        indices[numEntries] = globalRow;
                        values[numEntries]  = (2.0 /*+ subx*/) / ( dx[0] * dx[0] ) +
                                        (2.0 /*+ suby*/) / ( dx[1] * dx[1] ) +
                                        (2.0 /*+ subz*/) / ( dx[2] * dx[2] );
                        ++numEntries;
                    }
                    
                    
                    int error = A_p->InsertGlobalValues(globalRow, numEntries, &values[0], &indices[0]);
                    
                    if ( error != 0 )
                        throw std::runtime_error("Error in filling the matrix!");
                }
            }
        }
    }
    
    A_p->FillComplete(true);
    
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
