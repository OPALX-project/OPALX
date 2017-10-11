template <class AmrMultiGridLevel>
AmrLagrangeInterpolater<AmrMultiGridLevel>::AmrLagrangeInterpolater(Order order)
    : AmrInterpolater<AmrMultiGridLevel>( int(order) + 1 )
{ }


template <class AmrMultiGridLevel>
void AmrLagrangeInterpolater<AmrMultiGridLevel>::stencil(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    AmrMultiGridLevel* mglevel)
{
    
}


template <class AmrMultiGridLevel>
void AmrLagrangeInterpolater<AmrMultiGridLevel>::coarse(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    int dir, int shift, const amrex::BoxArray& ba,
    const AmrIntVect_t& riv,
    AmrMultiGridLevel* mglevel)
{
    // polynomial degree = #points - 1
    switch ( this->nPoints_m - 1 ) {
        
        case Order::LINEAR:
            this->crseLinear_m(iv, indices, values, dir, shift, ba, riv, mglevel);
            break;
        case Order::QUADRATIC:
            this->crseQuadratic_m(iv, indices, values, dir, shift, ba, riv, mglevel);
            break;
        default:
            std::runtime_error("Not implemented interpolation");
    }
}


template <class AmrMultiGridLevel>
void AmrLagrangeInterpolater<AmrMultiGridLevel>::fine(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    int dir, int shift, const amrex::BoxArray& ba,
    AmrMultiGridLevel* mglevel)
{
    // polynomial degree = #points - 1
    switch ( this->nPoints_m - 1 ) {
        
        case Order::LINEAR:
            this->fineLinear_m(iv, indices, values, dir, shift, ba, mglevel);
            break;
        case Order::QUADRATIC:
            this->fineQuadratic_m(iv, indices, values, dir, shift, ba, mglevel);
            break;
        default:
            std::runtime_error("Not implemented interpolation");
    }
}


template <class AmrMultiGridLevel>
void AmrLagrangeInterpolater<AmrMultiGridLevel>::fineLinear_m(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    int dir, int shift, const amrex::BoxArray& ba,
    AmrMultiGridLevel* mglevel)
{
    /*
     * computes the ghost cell directly
     */
    AmrIntVect_t tmp = iv;
    // first fine cell on refined coarse cell (closer to interface)
    tmp[dir] += shift;
    indices.push_back( mglevel->serialize(tmp) );
    values.push_back( 2.0 );
    
    // second fine cell on refined coarse cell (further away from interface)
    tmp[dir] += shift;
    indices.push_back( mglevel->serialize(tmp) );
    values.push_back( -1.0 );
}


template <class AmrMultiGridLevel>
void AmrLagrangeInterpolater<AmrMultiGridLevel>::fineQuadratic_m(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    int dir, int shift, const amrex::BoxArray& ba,
    AmrMultiGridLevel* mglevel)
{
    AmrIntVect_t tmp = iv;
    // first fine cell on refined coarse cell (closer to interface)
    tmp[dir] += shift;
    indices.push_back( mglevel->serialize(tmp) );
    values.push_back( 2.0 / 3.0 );
                        
    // second fine cell on refined coarse cell (further away from interface)
    tmp[dir] += shift;
    indices.push_back( mglevel->serialize(tmp) );
    values.push_back( -0.2 );
}


template <class AmrMultiGridLevel>
void AmrLagrangeInterpolater<AmrMultiGridLevel>::crseLinear_m(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    int dir, int shift, const amrex::BoxArray& ba,
    const AmrIntVect_t& riv,
    AmrMultiGridLevel* mglevel)
{
    //TODO Extend to 3D
    
    // y_t + y_b
    indices.push_back( mglevel->serialize(iv) );
    values.push_back( /*shift **/ 2.0 );
    
    // the neighbour cancels out
}


template <class AmrMultiGridLevel>
void AmrLagrangeInterpolater<AmrMultiGridLevel>::crseQuadratic_m(
    const AmrIntVect_t& iv,
    typename AmrMultiGridLevel::indices_t& indices,
    typename AmrMultiGridLevel::coefficients_t& values,
    int dir, int shift, const amrex::BoxArray& ba,
    const AmrIntVect_t& riv,
    AmrMultiGridLevel* mglevel)
{
#if BL_SPACEDIM == 2
    
    bool top = (iv[(dir+1)%BL_SPACEDIM] % 2 == 1);
    
    // right / upper / back
    AmrIntVect_t niv = iv;
    niv[(dir+1)%BL_SPACEDIM ] += 1;
    
    // left / lower / front
    AmrIntVect_t miv = iv;
    miv[(dir+1)%BL_SPACEDIM ] -= 1;
    
    // 2nd right / upper / back
    AmrIntVect_t n2iv = niv;
    n2iv[(dir+1)%BL_SPACEDIM ] += 1;
    
    // 2nd left / lower / front
    AmrIntVect_t m2iv = miv;
    m2iv[(dir+1)%BL_SPACEDIM ] -= 1;
        
    /* 3 cases:
     * --------
     * r: right neighbour of iv (center of Laplacian)
     * u: upper neighbour of iv
     * b: back neighbour of iv
     * 
     * l: lower / left neighbour of iv
     * f: front neighbour of iv
     * 
     * -1 --> not valid, error happend
     * 0 --> r / u / b and 2nd r / u / b are valid
     * 1 --> direct neighbours (r / u / b and l / f) of iv are valid (has priority)
     * 2 --> l / f and 2nd l / f are valid
     */
    
    // check r / u / b --> 1: valid; 0: not valid
    bool rub = !ba.contains(niv);
    
    // check l / f --> 1: valid; 0: not valid
    bool lf = !ba.contains(miv);
    
    // check 2nd r / u / b
    bool rub2 = !ba.contains(n2iv);
    
    // check 2nd l / f
    bool lf2 = !ba.contains(m2iv);
    
    if ( rub && lf )
    {
        /*
         * standard case -1, +1 are not-refined nor at physical/mesh boundary
         */            
        // cell is not refined and not at physical boundary
        
        indices.push_back( mglevel->serialize(iv) );
        
        // y_t or y_b
        values.push_back( 0.5 );
        
        //                     y_t          y_b
        double value = (top) ? 1.0 / 12.0 : -0.05;
        if ( mglevel->isBoundary(niv) ) {
            mglevel->applyBoundary(niv, indices, values, value);
        } else {
            indices.push_back( mglevel->serialize(niv) );
            values.push_back( value );
        }
        
        //              y_t     y_b
        value = (top) ? -0.05 : 1.0 / 12.0;
        if ( mglevel->isBoundary(miv) ) {
            mglevel->applyBoundary(miv, indices, values, value);
        } else {
            indices.push_back( mglevel->serialize(miv) );
            values.push_back( value );
        }
        
    } else if ( rub && rub2 ) {
        /*
         * corner case --> right / upper / back + 2nd right / upper / back
         */
        //                     y_t          y_b
        double value = (top) ? 7.0 / 20.0 : 0.75;
        indices.push_back( mglevel->serialize(iv) );
        values.push_back( value );
        
        //              y_t          y_b
        value = (top) ? 7.0 / 30.0 : -0.3;
        if ( mglevel->isBoundary(niv) ) {
            mglevel->applyBoundary(niv, indices, values, value);
        } else {
            indices.push_back( mglevel->serialize(niv) );
            values.push_back( value );
        }
        
        //              y_t     y_b
        value = (top) ? -0.05 : 1.0 / 12.0;
        if ( mglevel->isBoundary(n2iv) ) {
            mglevel->applyBoundary(n2iv, indices, values, value);
        } else {
            indices.push_back( mglevel->serialize(n2iv) );
            values.push_back( value );
        }
        
    } else if ( lf && lf2 ) {
        /*
         * corner case --> left / lower / front + 2nd left / lower / front
         */
        //                     y_t    y_b
        double value = (top) ? 0.75 : 7.0 / 20.0;
        indices.push_back( mglevel->serialize(iv) );
        values.push_back( value );
        
        //              y_t           y_b
        value = (top) ? -0.3 :  7.0 / 30;
        if ( mglevel->isBoundary(miv) ) {
            mglevel->applyBoundary(miv, indices, values, value);
        } else {
            indices.push_back( mglevel->serialize(miv) );
            values.push_back( value );
        }
        
        //              y_t          y_b
        value = (top) ? 1.0 / 12.0 : -0.05;
        if ( mglevel->isBoundary(m2iv) ) {
            mglevel->applyBoundary(m2iv, indices, values, value);
        } else {
            indices.push_back( mglevel->serialize(m2iv) );
            values.push_back( value );
        }
        
    } else {
        // last trial: linear Lagrange interpolation
        
        if ( rub || lf ) {
            this->crseLinear_m(iv, indices, values, dir, shift, ba, riv, mglevel);
        } else
            std::runtime_error("Lagrange Error: No valid scenario found!");
    }
    
#elif BL_SPACEDIM == 3
    
    /*
     * check in 5x5 area (using iv as center) if 9 cells are not covered
     */
    int d1 = (dir+1)%BL_SPACEDIM;
    int d2 = (dir+2)%BL_SPACEDIM;
    
    bits_t area;
    int bit = 0;
    
    AmrIntVect_t tmp = iv;
    for (int i = -2; i < 3; ++i) {
        tmp[d1] += i;
        for (int j = -2; j < 3; ++j) {
            
            tmp[d2] += j;
            
            area[bit] = !ba.contains(tmp);
            ++bit;
            
            // undo
            tmp[d2] -= j;
        }
        // undo
        tmp[d1] -= i;
    }
    
    bool found = false;
    pattern_t::const_iterator pit = std::begin(pattern_m);
    
    while ( !found && pit != std::end(pattern_m) ) {
        if ( *pit == (area & bits_t(*pit)).to_ulong() )
            break;
        ++pit;
    }
    
    switch ( *pit ) {
        case this->pattern_m[0]:
        {
            // cross pattern
            
            break;
        }
        case this->pattern_m[1]:
        {
            // T pattern
            break;
        }
        case this->pattern_m[2]:
        {
            // T on head pattern
            break;
        }
        case this->pattern_m[3]:
        {
            // upper right corner pattern
            break;
        }
        case this->pattern_m[4]:
        {
            // upper left corner pattern
            break;
        }
        case this->pattern_m[5]:
        {
            // mirrored L pattern
            break;
        }
        case this->pattern_m[6]:
        {
            // L pattern
            break;
        }
        case this->pattern_m[7]:
        {
            // left hammer pattern
            break;
        }
        case this->pattern_m[8]:
        {
            // right hammer pattern
            break;
        }
        default:
            // unknown pattern
            break;
    }
    
    
    
    
    /*          x   y   z
     * ------------------
     * dir:     0   1   2
     * top1:    1   2   0
     * top2:    2   0   1
     */
    bool top1 = (riv[d1] % 2 == 1);
    bool top2 = (riv[d2] % 2 == 1);
    
    
    /* There are 9 coefficients from Lagrange interpolation.
     * Those are given by the product of one of
     * L0, L1, L2 and one of K0, K1, K2.
     * 
     * g(x, y) = f(x0, y0) * L0(x) * K0(y) +
     *           f(x0, y1) * L0(x) * K1(y) +
     *           f(x0, y2) * L0(x) * K2(y) +
     *           f(x1, y0) * L1(x) * K0(y) +
     *           f(x1, y1) * L1(x) * K1(y) +
     *           f(x1, y2) * L1(x) * K2(y) +
     *           f(x2, y0) * L2(x) * K0(y) +
     *           f(x2, y1) * L2(x) * K1(y) +
     *           f(x2, y2) * L2(x) * K2(y) +
     */
    
#else
    #error Lagrange interpolation: Only 2D and 3D are supported!
#endif
}
