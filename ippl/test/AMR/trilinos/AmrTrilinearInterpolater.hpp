template <class AmrMultiGridLevel>
AmrTrilinearInterpolater<AmrMultiGridLevel>::AmrTrilinearInterpolater()
    : AmrInterpolater<AmrMultiGridLevel>(2 << (AMREX_SPACEDIM - 1))
{ }


template <class AmrMultiGridLevel>
void AmrTrilinearInterpolater<AmrMultiGridLevel>::stencil(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    AmrMultiGridLevel* mglevel)
{
    /* lower left coarse cell (i, j, k)
     * floor( i - 0.5 ) / rr[0]
     * floor( j - 0.5 ) / rr[1]
     * floor( k - 0.5 ) / rr[2]
     */
    AmrIntVect_t civ;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            
        scalar_t tmp = iv[d] - 0.5;
        if ( std::signbit(tmp) )
            civ[d] = std::floor(tmp);
        else
            civ[d] = tmp;
    }
    
    civ.coarsen(mglevel->refinement());
        
    // ref ratio 2 only
    scalar_t dx = 0.5 * ( iv[0] - civ[0] * 2 ) - 0.25;
    scalar_t dy = 0.5 * ( iv[1] - civ[1] * 2 ) - 0.25;
#if AMREX_SPACEDIM == 3
    scalar_t dz = 0.5 * ( iv[2] - civ[2] * 2 ) - 0.25;
#endif
        
    scalar_t xdiff = 1.0 - dx;
    scalar_t ydiff = 1.0 - dy;
#if AMREX_SPACEDIM == 3
    scalar_t zdiff = 1.0 - dz;
#endif
    // (i, j, k)
    go_t crse_gidx = mglevel->serialize(civ);
    scalar_t value = AMREX_D_TERM(xdiff, * ydiff, * zdiff) * scale;
    
    if ( mglevel->isBoundary(civ) ) {
        mglevel->applyBoundary(civ, map, value);
    } else {
	map[crse_gidx] += value;
    }
    
    // (i+1, j, k)
    AmrIntVect_t tmp(D_DECL(civ[0]+1, civ[1], civ[2]));
    value = AMREX_D_TERM(dx, * ydiff, * zdiff) * scale;
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, map, value);
    } else {
        map[mglevel->serialize(tmp)] += value;
    }
    
    // (i, j+1, k)
    tmp = AmrIntVect_t(D_DECL(civ[0], civ[1]+1, civ[2]));
    value = AMREX_D_TERM(xdiff, * dy, * zdiff) * scale;
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, map, value);
    } else {
	map[mglevel->serialize(tmp)] += value;
    }
    
    // (i+1, j+1, k)
    tmp = AmrIntVect_t(D_DECL(civ[0]+1, civ[1]+1, civ[2]));
    value = AMREX_D_TERM(dx, * dy, * zdiff) * scale;
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, map, value);
    } else {
        map[mglevel->serialize(tmp)] += value;
    }
        
#if AMREX_SPACEDIM == 3
    // (i, j, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0], civ[1], civ[2]+1));
    value = AMREX_D_TERM(xdiff, * ydiff, * dz) * scale;
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, map, value);
    } else {
        map[mglevel->serialize(tmp)] += value;
    }
    
    // (i+1, j, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0]+1, civ[1], civ[2]+1));
    value = AMREX_D_TERM(dx, * ydiff, * dz) * scale;
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, map, value);
    } else {
        map[mglevel->serialize(tmp)] += value;
    }
    
    // (i, j+1, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0], civ[1]+1, civ[2]+1));
    value = AMREX_D_TERM(xdiff, * dy, * dz) * scale;
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, map, value);
    } else {
        map[mglevel->serialize(tmp)] += value;
    }
    
    // (i+1, j+1, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0]+1, civ[1]+1, civ[2]+1));
    value = AMREX_D_TERM(dx, * dy, * dz) * scale;
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, map, value);
    } else {
	map[mglevel->serialize(tmp)] += value;
    }
#endif
}


template <class AmrMultiGridLevel>
void AmrTrilinearInterpolater<AmrMultiGridLevel>::coarse(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    lo_t dir, lo_t shift, const amrex::BoxArray& ba,
    const AmrIntVect_t& riv,
    AmrMultiGridLevel* mglevel)
{
    // do nothing
}


template <class AmrMultiGridLevel>
void AmrTrilinearInterpolater<AmrMultiGridLevel>::fine(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    lo_t dir, lo_t shift, const amrex::BoxArray& ba,
    AmrMultiGridLevel* mglevel)
{
    /*
     * The AmrTrilinearInterpolater interpolates directly to the
     * fine ghost cell.
     */
    this->stencil(iv, map, scale, mglevel);
}
