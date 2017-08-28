template <class AmrMultiGridLevel>
AmrTrilinearInterpolater<AmrMultiGridLevel>::AmrTrilinearInterpolater() : nPoints_m(2 << (BL_SPACEDIM - 1) )
{ }


template <class AmrMultiGridLevel>
const int& AmrTrilinearInterpolater<AmrMultiGridLevel>::getNumberOfPoints() const {
    return nPoints_m;
}


template <class AmrMultiGridLevel>
void AmrTrilinearInterpolater<AmrMultiGridLevel>::stencil(const AmrIntVect_t& iv,
                                                          typename AmrMultiGridLevel::indices_t& indices,
                                                          typename AmrMultiGridLevel::coefficients_t& values,
                                                          int& numEntries,
                                                          AmrMultiGridLevel* mglevel)
{
    /* lower left coarse cell (i, j, k)
     * floor( i - 0.5 ) / rr[0]
     * floor( j - 0.5 ) / rr[1]
     * floor( k - 0.5 ) / rr[2]
     */
    AmrIntVect_t civ;
    for (int d = 0; d < BL_SPACEDIM; ++d) {
            
        double tmp = iv[d] - 0.5;
        if ( std::signbit(tmp) )
            civ[d] = std::floor(tmp);
        else
            civ[d] = tmp;
    }
    
    civ.coarsen(mglevel->refinement());
        
    // ref ratio 2 only
    double dx = 0.5 * ( iv[0] - civ[0] * 2 ) - 0.25;
    double dy = 0.5 * ( iv[1] - civ[1] * 2 ) - 0.25;
#if BL_SPACEDIM == 3
    double dz = 0.5 * ( iv[2] - civ[2] * 2 ) - 0.25;
#endif
        
    double xdiff = 1.0 - dx;
    double ydiff = 1.0 - dy;
#if BL_SPACEDIM == 3
    double zdiff = 1.0 - dz;
#endif
    // (i, j, k)
    int crse_gidx = mglevel->serialize(civ);
    double value = AMREX_D_TERM(xdiff, * ydiff, * zdiff);
    
    if ( mglevel->isBoundary(civ) ) {
        mglevel->applyBoundary(civ, indices, values, numEntries, value);
    } else {
        indices[numEntries] = crse_gidx;
        values[numEntries]  = value;
        ++numEntries;
    }
    
    // (i+1, j, k)
    AmrIntVect_t tmp(D_DECL(civ[0]+1, civ[1], civ[2]));
    value = AMREX_D_TERM(dx, * ydiff, * zdiff);
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, indices, values, numEntries, value);
    } else {
        indices[numEntries] = mglevel->serialize(tmp);
        values[numEntries]  = value;
        ++numEntries;
    }
    
    // (i, j+1, k)
    tmp = AmrIntVect_t(D_DECL(civ[0], civ[1]+1, civ[2]));
    value = AMREX_D_TERM(xdiff, * dy, * zdiff);
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, indices, values, numEntries, value);
    } else {
        indices[numEntries] = mglevel->serialize(tmp);
        values[numEntries]  = value;
        ++numEntries;
    }
    
    // (i+1, j+1, k)
    tmp = AmrIntVect_t(D_DECL(civ[0]+1, civ[1]+1, civ[2]));
    value = AMREX_D_TERM(dx, * dy, * zdiff);
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, indices, values, numEntries, value);
    } else {
        indices[numEntries] = mglevel->serialize(tmp);
        values[numEntries]  = value;
        ++numEntries;
    }
        
#if BL_SPACEDIM == 3
    // (i, j, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0], civ[1], civ[2]+1));
    value = AMREX_D_TERM(xdiff, * ydiff, * dz);
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, indices, values, numEntries, value);
    } else {
        indices[numEntries] = mglevel->serialize(tmp);
        values[numEntries]  = value;
        ++numEntries;
    }
    
    // (i+1, j, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0]+1, civ[1], civ[2]+1));
    value = AMREX_D_TERM(dx, * ydiff, * dz);
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, indices, values, numEntries, value);
    } else {
        indices[numEntries] = mglevel->serialize(tmp);
        values[numEntries]  = value;
        ++numEntries;
    }
    
    // (i, j+1, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0], civ[1]+1, civ[2]+1));
    value = AMREX_D_TERM(xdiff, * dy, * dz);
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, indices, values, numEntries, value);
    } else {
        indices[numEntries] = mglevel->serialize(tmp);
        values[numEntries]  = value;
        ++numEntries;
    }
    
    // (i+1, j+1, k+1)
    tmp = AmrIntVect_t(D_DECL(civ[0]+1, civ[1]+1, civ[2]+1));
    value = AMREX_D_TERM(dx, * dy, * dz);
    if ( mglevel->isBoundary(tmp) ) {
        mglevel->applyBoundary(tmp, indices, values, numEntries, value);
    } else {
        indices[numEntries] = mglevel->serialize(tmp);
        values[numEntries]  = value;
        ++numEntries;
    }
#endif
}
