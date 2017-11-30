template <class AmrMultiGridLevel>
AmrPCInterpolater<AmrMultiGridLevel>::AmrPCInterpolater()
    : AmrInterpolater<AmrMultiGridLevel>(1)
{ }


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::stencil(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::umap_t& map,
    const typename AmrMultiGridLevel::scalar_t& scale,
    AmrMultiGridLevel* mglevel)
{
    AmrIntVect_t civ = iv;
    civ.coarsen(mglevel->refinement());

    int crse_gidx = mglevel->serialize(civ);
    
    if ( mglevel->isBoundary(civ) ) {
        mglevel->applyBoundary(civ, map, scale);
    } else {
        map[crse_gidx] += scale;
    }
}


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::coarse(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::umap_t& map,
    const typename AmrMultiGridLevel::scalar_t& scale,
    int dir, int shift, const amrex::BoxArray& ba,
    const AmrIntVect_t& riv,
    AmrMultiGridLevel* mglevel)
{
    // do nothing
}


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::fine(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::umap_t& map,
    const typename AmrMultiGridLevel::scalar_t& scale,
    int dir, int shift, const amrex::BoxArray& ba,
    AmrMultiGridLevel* mglevel)
{
    /*
     * The AmrPCInterpolater interpolates directly to the
     * fine ghost cell.
     */
    this->stencil(iv, map, scale, mglevel);
}
