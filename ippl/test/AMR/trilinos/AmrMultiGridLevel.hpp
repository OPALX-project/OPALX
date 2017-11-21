template <class MatrixType, class VectorType>
AmrMultiGridLevel<MatrixType,
                  VectorType>::AmrMultiGridLevel(const amrex::BoxArray& _grids,
                                                 const amrex::DistributionMapping& _dmap,
                                                 const AmrGeometry_t& _geom,
                                                 const AmrIntVect_t& rr,
                                                 boundary_t* const bc,
                                                 const Teuchos::RCP<comm_t>& comm,
                                                 const Teuchos::RCP<node_t>& node)
    : grids(_grids),
      dmap(_dmap),
      geom(_geom),
      map_p(Teuchos::null),
      Anf_p(Teuchos::null),
      R_p(Teuchos::null),
      I_p(Teuchos::null),
      Bcrse_p(Teuchos::null),
      Bfine_p(Teuchos::null),
      Awf_p(Teuchos::null),
      rho_p(Teuchos::null),
      phi_p(Teuchos::null),
      residual_p(Teuchos::null),
      error_p(Teuchos::null),
      UnCovered_p(Teuchos::null),
      rr_m(rr),
      bc_mp(bc)
{
    for (int j = 0; j < AMREX_SPACEDIM; ++j) {
        G_p[j] = Teuchos::null;
        
        nr_m[j] = _geom.Domain().length(j);
    }
    
    this->buildLevelMask_m();
    
    this->buildMap_m(comm, node);
    
    
    residual_p = Teuchos::rcp( new vector_t(map_p, false) );
    error_p = Teuchos::rcp( new vector_t(map_p, false) );
}


template <class MatrixType, class VectorType>
AmrMultiGridLevel<MatrixType, VectorType>::~AmrMultiGridLevel()
{
    map_p = Teuchos::null;
    
    Anf_p = Teuchos::null;
    R_p = Teuchos::null;
    I_p = Teuchos::null;
    Bcrse_p = Teuchos::null;
    Bfine_p = Teuchos::null;
    Awf_p = Teuchos::null;
    
    for (int j = 0; j < AMREX_SPACEDIM; ++j)
        G_p[j] = Teuchos::null;
    
    UnCovered_p = Teuchos::null;
    
    rho_p = Teuchos::null;
    phi_p = Teuchos::null;
    residual_p = Teuchos::null;
    error_p = Teuchos::null;
}


template <class MatrixType, class VectorType>
int AmrMultiGridLevel<MatrixType, VectorType>::serialize(const AmrIntVect_t& iv) const {
#if AMREX_SPACEDIM == 3
    return iv[2] + (iv[1] + nr_m[1] * iv[0]) * nr_m[2];
#else
    return iv[1] + iv[0] * nr_m[1];
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
void AmrMultiGridLevel<MatrixType, VectorType>::buildMap_m(const Teuchos::RCP<comm_t>& comm,
                                                           const Teuchos::RCP<node_t>& node)
{
    
    int localNumElements = 0;
    coefficients_t values;
//     indices_t globalindices;
    
    Teuchos::Array<global_ordinal_t> globalindices;
    
    for (amrex::MFIter mfi(grids, dmap, true); mfi.isValid(); ++mfi) {
        const amrex::Box&    tbx  = mfi.tilebox();
        const int* lo = tbx.loVect();
        const int* hi = tbx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));

                    int globalidx = serialize(iv);
                    
                    globalindices.push_back(globalidx);
                    
                    ++localNumElements;
#if AMREX_SPACEDIM == 3
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
    amr::global_ordinal_t baseIndex = serialize(lowcorner);
    
    // numGlobalElements == N
    int N = grids.numPts();
    
    map_p = Teuchos::rcp( new dmap_t(N, globalindices, baseIndex, comm, node) );
}
