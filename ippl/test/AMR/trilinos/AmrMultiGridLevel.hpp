template <class MatrixType, class VectorType>
AmrMultiGridLevel<MatrixType,
                  VectorType>::AmrMultiGridLevel(const amrex::BoxArray& _grids,
                                                 const amrex::DistributionMapping& _dmap,
                                                 const AmrGeometry_t& _geom,
                                                 const AmrIntVect_t& rr,
                                                 AmrBoundary<AmrMultiGridLevel<MatrixType, VectorType> >* bc,
                                                 Epetra_MpiComm& comm)
    : grids(_grids),
      dmap(_dmap),
      geom(_geom),
      rr_m(rr),
      map_p(Teuchos::null),
      A_p(Teuchos::null),
      B_p(Teuchos::null),
      R_p(Teuchos::null),
      I_p(Teuchos::null),
      Bcrse_p(Teuchos::null),
      Bfine_p(Teuchos::null),
      As_p(Teuchos::null),
      S_p(Teuchos::null),
      rho_p(Teuchos::null),
      phi_p(Teuchos::null),
      residual_p(Teuchos::null),
      error_p(Teuchos::null),
      UnCovered_p(Teuchos::null)
{
    bc_mp.reset(bc);
    
    for (int j = 0; j < BL_SPACEDIM; ++j)
        nr_m[j] = _geom.Domain().length(j);
    
    this->buildLevelMask_m();
    
    this->buildMap_m(comm);
    
    
    residual_p = Teuchos::rcp( new vector_t(*map_p, false) );
    error_p = Teuchos::rcp( new vector_t(*map_p, false) );
}


template <class MatrixType, class VectorType>
AmrMultiGridLevel<MatrixType, VectorType>::~AmrMultiGridLevel()
{
    map_p = Teuchos::null;
    
    A_p = Teuchos::null;
    B_p = Teuchos::null;
    R_p = Teuchos::null;
    I_p = Teuchos::null;
    Bcrse_p = Teuchos::null;
    Bfine_p = Teuchos::null;
    As_p = Teuchos::null;
    S_p = Teuchos::null;
    UnCovered_p = Teuchos::null;
    
    rho_p = Teuchos::null;
    phi_p = Teuchos::null;
    residual_p = Teuchos::null;
    error_p = Teuchos::null;
}


template <class MatrixType, class VectorType>
int AmrMultiGridLevel<MatrixType, VectorType>::serialize(const AmrIntVect_t& iv) const {
#if BL_SPACEDIM == 3
    return iv[0] + (iv[1] + nr_m[1] * iv[2]) * nr_m[0];
#else
    return iv[0] + iv[1] * nr_m[0];
#endif
}


template <class MatrixType, class VectorType>
bool AmrMultiGridLevel<MatrixType, VectorType>::isBoundary(const AmrIntVect_t& iv) const {
    return bc_mp->isBoundary(iv, &nr_m[0]);
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::applyBoundary(const AmrIntVect_t& iv,
                                                              indices_t& indices,
                                                              coefficients_t& values,
                                                              const double& value)
{
    bc_mp->apply(iv, indices, values, value, this, &nr_m[0]);
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::buildLevelMask_m() {
    mask.reset(new amrex::FabArray<amrex::BaseFab<int> >(grids, dmap, 1, 1));
    mask->BuildMask(geom.Domain(), geom.periodicity(),
                    Mask::COVERED, Mask::BNDRY,
                    Mask::PHYSBNDRY, Mask::INTERIOR);
}


template <class MatrixType, class VectorType>
const AmrIntVect_t& AmrMultiGridLevel<MatrixType, VectorType>::refinement() const {
    return rr_m;
}


template <class MatrixType, class VectorType>
void AmrMultiGridLevel<MatrixType, VectorType>::buildMap_m(const Epetra_MpiComm& comm) {
    
    int localNumElements = 0;
    coefficients_t values;
    indices_t globalindices;
    
    for (amrex::MFIter mfi(grids, dmap, false); mfi.isValid(); ++mfi) {
        const amrex::Box&    bx  = mfi.validbox();  
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));

                    int globalidx = serialize(iv);
                    
                    globalindices.push_back(globalidx);
                    
                    ++localNumElements;
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    /*
     * create map that specifies which processor gets which data
     */
    
    // get smallest global index of this level
    amrex::Box bx = grids.minimalBox();
    const int* lo = bx.loVect();
    AmrIntVect_t lowcorner(D_DECL(lo[0], lo[1], lo[2]));
    
    // where to start indexing
    const int baseIndex = serialize(lowcorner);
    
    // numGlobalElements == N
    int N = grids.numPts();
    
    std::cout << "N = " << N << " baseIndex = " << baseIndex << " localNumElements = " << localNumElements << std::endl;
    
    map_p = Teuchos::rcp( new Epetra_Map(N, localNumElements,
                                         &globalindices[0], baseIndex, comm) );
}
