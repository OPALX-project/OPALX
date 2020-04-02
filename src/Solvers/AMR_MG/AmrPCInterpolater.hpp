// Implementation file of the AmrPCInterpolater class,
//   piecewise constant interpolation of grid data from coarse to fine cells.
//
// Copyright (c) 2017 - 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// OPAL is licensed under GNU GPL version 3.

template <class Level>
AmrPCInterpolater<Level>::AmrPCInterpolater()
    : AmrInterpolater<Level>(1)
{ }


template <class Level>
void AmrPCInterpolater<Level>::stencil(
    const AmrIntVect_t& iv,
    const basefab_t& fab,
    umap_t& map,
    const scalar_t& scale,
    Level* mglevel)
{
    AmrIntVect_t civ = iv;
    civ.coarsen(mglevel->refinement());

    if ( !mglevel->applyBoundary(civ, fab, map, scale) )
        map[mglevel->serialize(civ)] += scale;
}


template <class Level>
void AmrPCInterpolater<Level>::coarse(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    lo_t dir, lo_t shift, const basefab_t& rfab,
    const AmrIntVect_t& riv,
    Level* mglevel)
{
    // do nothing
}


template <class Level>
void AmrPCInterpolater<Level>::fine(
    const AmrIntVect_t& iv,
    umap_t& map,
    const scalar_t& scale,
    lo_t dir, lo_t shift, const basefab_t& fab,
    Level* mglevel)
{
    /*
     * The AmrPCInterpolater interpolates directly to the
     * fine ghost cell.
     */
    this->stencil(iv, fab, map, scale, mglevel);
}
