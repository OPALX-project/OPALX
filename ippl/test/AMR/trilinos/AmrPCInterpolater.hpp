template <class AmrMultiGridLevel>
AmrPCInterpolater<AmrMultiGridLevel>::AmrPCInterpolater()
    : AmrInterpolater<AmrMultiGridLevel>(1)
{ }


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::stencil(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    AmrMultiGridLevel* mglevel)
{
    AmrIntVect_t civ = iv;
    civ.coarsen(mglevel->refinement());

    int crse_gidx = mglevel->serialize(civ);
    double value = 1.0;
    
    if ( mglevel->isBoundary(civ) ) {
        mglevel->applyBoundary(civ, indices, values, value);
    } else {
        indices.push_back( crse_gidx );
        values.push_back( value );
    }
}


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::coarse(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    int dir, int shift, const amrex::BoxArray& ba,
    const AmrIntVect_t& riv,
    AmrMultiGridLevel* mglevel)
{
    // do nothing
}


template <class AmrMultiGridLevel>
void AmrPCInterpolater<AmrMultiGridLevel>::fine(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    int dir, int shift, const amrex::BoxArray& ba,
    AmrMultiGridLevel* mglevel)
{
    /*
     * The AmrPCInterpolater interpolates directly to the
     * fine ghost cell.
     */
    this->stencil(iv, indices, values, mglevel);
}