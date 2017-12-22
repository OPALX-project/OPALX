template <class AmrMultiGridLevel>
AmrPCInterpolater<AmrMultiGridLevel>::AmrPCInterpolater()
    : AmrInterpolater<AmrMultiGridLevel>(1)
{ }


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::stencil(
    const AmrIntVect_t& iv,
    const basefab_t& fab,
    umap_t& map,
    const scalar_t& scale,
    AmrMultiGridLevel* mglevel)
{
    AmrIntVect_t civ = iv;
    civ.coarsen(mglevel->refinement());

    if ( !mglevel->applyBoundary(civ, fab, map, scale) )
	map[mglevel->serialize(civ)] += scale;
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
    lo_t dir, lo_t shift, const basefab_t& fab,
    AmrMultiGridLevel* mglevel)
{
    /*
     * The AmrPCInterpolater interpolates directly to the
     * fine ghost cell.
     */
    this->stencil(iv, fab, map, scale, mglevel);
}
