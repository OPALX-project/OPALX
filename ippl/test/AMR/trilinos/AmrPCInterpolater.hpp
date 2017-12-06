template <class AmrMultiGridLevel>
AmrPCInterpolater<AmrMultiGridLevel>::AmrPCInterpolater()
    : AmrInterpolater<AmrMultiGridLevel>(1)
{ }


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::stencil(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    AmrMultiGridLevel* mglevel)
{
    AmrIntVect_t civ = iv;
    civ.coarsen(mglevel->refinement());

    go_t crse_gidx = mglevel->serialize(civ);
    
    if ( mglevel->isBoundary(civ) ) {
        mglevel->applyBoundary(civ, map, scale);
    } else {
        map[crse_gidx] += scale;
    }
}


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::coarse(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    lo_t dir, lo_t shift, const basefab_t& rfab,
    const AmrIntVect_t& riv,
    AmrMultiGridLevel* mglevel)
{
    // do nothing
}


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::fine(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    lo_t dir, lo_t shift, const amrex::BoxArray& ba,
    AmrMultiGridLevel* mglevel)
{
    /*
     * The AmrPCInterpolater interpolates directly to the
     * fine ghost cell.
     */
    this->stencil(iv, map, scale, mglevel);
}
