template <class AmrMultiGridLevel>
constexpr typename AmrLagrangeInterpolater<AmrMultiGridLevel>::qpattern_t
    AmrLagrangeInterpolater<AmrMultiGridLevel>::qpattern_ms;


template <class AmrMultiGridLevel>
constexpr typename AmrLagrangeInterpolater<AmrMultiGridLevel>::lpattern_t
    AmrLagrangeInterpolater<AmrMultiGridLevel>::lpattern_ms;


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
#if BL_SPACEDIM == 2
    bool top = (riv[(dir+1)%BL_SPACEDIM] % 2 == 1);
    
    // right / upper / back
    AmrIntVect_t niv = iv;
    niv[(dir+1)%BL_SPACEDIM ] += 1;
    
    // left / lower / front
    AmrIntVect_t miv = iv;
    miv[(dir+1)%BL_SPACEDIM ] -= 1;
    
    // factor for fine
    double fac = 8.0 / 15.0;
    
    if ( !ba.contains(niv) ) {
        // check r / u / b --> 1: valid; 0: not valid
        
        indices.push_back( mglevel->serialize(iv) );
        //                        y_t    y_b
        values.push_back( fac * (top) ? 0.75 : 1.25 );
        
        indices.push_back( mglevel->serialize(niv) );
        //                        y_t    y_b
        values.push_back( fac * (top) ? 0.25 : -0.25 );
        
        
    } else if ( !ba.contains(miv) ) {
        // check l / f --> 1: valid; 0: not valid
        
        indices.push_back( mglevel->serialize(iv) );
        //                        y_t    y_b
        values.push_back( fac * (top) ? 0.75 : 1.25 );
        
        indices.push_back( mglevel->serialize(miv) );
        //                        y_t    y_b
        values.push_back( fac * (top) ? 0.25 : -0.25 );
        
    } else
        throw std::runtime_error("Lagrange Error: No valid scenario found!");
    
#elif BL_SPACEDIM == 3
    
    /*          x   y   z
     * ------------------
     * dir:     0   1   2
     * top1:    1   2   0   (i, L)
     * top2:    2   0   1   (j, K)
     */
    
    /* There are 4 coefficients from Lagrange interpolation.
     * Those are given by the product of one of
     * L0, L1 and one of K0, K1.
     * 
     * g(x, y) = f(x0, y0) * L0(x) * K0(y) +
     *           f(x0, y1) * L0(x) * K1(y) +
     *           f(x1, y0) * L1(x) * K0(y) +
     *           f(x1, y1) * L1(x) * K1(y) +
     */
    
    /*
     * check in 3x3 area (using iv as center) if 4 cells are not covered
     */
    int d1 = (dir+1)%BL_SPACEDIM;
    int d2 = (dir+2)%BL_SPACEDIM;
    
    lbits_t area;
    int bit = 0;
    
    AmrIntVect_t tmp = iv;
    for (int i = -1; i < 2; ++i) {
        tmp[d1] += i;
        for (int j = -1; j < 2; ++j) {
            
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
    qpattern_t::const_iterator pit = std::begin(this->lpattern_ms);
    
    while ( !found && pit != std::end(this->lpattern_ms) ) {
        if ( *pit == (area & lbits_t(*pit)).to_ulong() )
            break;
        ++pit;
    }
    
    // factor for fine
    double fac = 8.0 / 15.0;
    
    // if pattern is known
    bool known = true;
    
    double L[2] = {0.0, 0.0};
    bool top1 = (riv[d1] % 2 == 1);
    
    double K[2] = {0.0, 0.0};
    bool top2 = (riv[d2] % 2 == 1);
    
    int begin[2] = { 0, 0 };
    int end[2]   = { 0, 0 };
    
    switch ( *pit ) {
        case this->lpattern_ms[0]:
        {
            // corner top left pattern
            L[0] = (top1) ? 0.25 : -0.25;   // L_{-1}
            L[1] = (top1) ? 0.75 : 1.25;    // L_{0}
            begin[0] = -1;
            end[0]   =  0;
            
            K[0] = (top2) ? 0.75 : 1.25;    // K_{0}
            K[1] = (top2) ? 0.25 : -0.25;   // K_{1}
            begin[1] = 0;
            end[1]   = 1;
            break;
        }
        case this->lpattern_ms[1]:
        {
            // corner top right pattern
            L[0] = (top1) ? 0.25 : -0.25;   // L_{-1}
            L[1] = (top1) ? 0.75 : 1.25;    // L_{0}
            begin[0] = -1;
            end[0]   =  0;
            
            K[0] = (top2) ? 0.25 : -0.25;   // K_{-1}
            K[1] = (top2) ? 0.75 : 1.25;    // K_{0}
            begin[1] = -1;
            end[1]   =  0;
            break;
        }
        case this->lpattern_ms[2]:
        {
            // corner bottom left pattern
            L[0] = (top1) ? 0.75 : 1.25;    // L_{0}
            L[1] = (top1) ? 0.25 : -0.25;   // L_{1}
            begin[0] = 0;
            end[0]   = 1;
            
            K[0] = (top2) ? 0.75 : 1.25;    // K_{0}
            K[1] = (top2) ? 0.25 : -0.25;   // K_{1}
            begin[1] = 0;
            end[1]   = 1;
            break;
        }
        case this->lpattern_ms[3]:
        {
            // corner bottom right pattern
            L[0] = (top1) ? 0.75 : 1.25;    // L_{0}
            L[1] = (top1) ? 0.25 : -0.25;   // L_{1}
            begin[0] = 0;
            end[0]   = 1;
            
            K[0] = (top2) ? 0.25 : -0.25;   // K_{-1}
            K[1] = (top2) ? 0.75 : 1.25;    // K_{0}
            begin[1] = -1;
            end[1]   =  0;
            break;
        }
        default:
            throw std::runtime_error("Lagrange Error: No valid scenario found!");
    }
    
    /*
     * if pattern is known --> add stencil
     */
    AmrIntVect_t niv = iv;
    for (int i = begin[0]; i <= end[0]; ++i) {
        niv[d1] += i;
        for (int j = begin[1]; j <= end[1]; ++j) {
            niv[d2] += j;
            
            double value = fac * L[i-begin[0]] * K[j-begin[1]];
            if ( mglevel->isBoundary(niv) ) {
                mglevel->applyBoundary(niv, indices, values, value);
            } else {
                indices.push_back( mglevel->serialize(niv) );
                values.push_back( value );
            }
            
            // undo
            niv[d2] -= j;
        }
        // undo
        niv[d1] -= i;
    }
#else
    #error Lagrange interpolation: Only 2D and 3D are supported!
#endif
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
    
    bool top = (riv[(dir+1)%BL_SPACEDIM] % 2 == 1);
    
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
        /* last trial: linear Lagrange interpolation
         * --> it throws an error if not possible
         */
        this->crseLinear_m(iv, indices, values, dir, shift, ba, riv, mglevel);
    }
    
#elif BL_SPACEDIM == 3
    
    /*          x   y   z
     * ------------------
     * dir:     0   1   2
     * top1:    1   2   0   (i, L)
     * top2:    2   0   1   (j, K)
     */
    
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
    
    
    /*
     * check in 5x5 area (using iv as center) if 9 cells are not covered
     */
    int d1 = (dir+1)%BL_SPACEDIM;
    int d2 = (dir+2)%BL_SPACEDIM;
    
    qbits_t area;
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
    qpattern_t::const_iterator pit = std::begin(this->qpattern_ms);
    
    while ( !found && pit != std::end(this->qpattern_ms) ) {
        if ( *pit == (area & qbits_t(*pit)).to_ulong() )
            break;
        ++pit;
    }
    
    // factor for fine
    double fac = 8.0 / 15.0;
    
    // if pattern is known
    bool known = true;
    
    double L[3] = {0.0, 0.0, 0.0};
    bool top1 = (riv[d1] % 2 == 1);
    
    double K[3] = {0.0, 0.0, 0.0};
    bool top2 = (riv[d2] % 2 == 1);
    
    int begin[2] = { 0, 0 };
    int end[2]   = { 0, 0 };
    
    switch ( *pit ) {
        case this->qpattern_ms[0]:
        {
            // cross pattern
            L[0] = (top1) ? -3.0 / 32.0 : 5.0 / 32.0;   // L_{-1}
            L[1] = 15.0 / 16.0;                         // L_{0}
            L[2] = (top1) ? 5.0 / 32.0 : -3.0 / 32.0;   // L_{1}
            begin[0] = -1;
            end[0]   =  1;
            
            K[0] = (top2) ? -3.0 / 32.0 : 5.0 / 32.0;   // K_{-1}
            K[1] = 15.0 / 16.0;                         // K_{0}
            K[2] = (top2) ? 5.0 / 32.0 : -3.0 / 32.0;   // K_{1}
            begin[1] = -1;
            end[1]   =  1;
            
            break;
        }
        case this->qpattern_ms[1]:
        {
            // T pattern
            L[0] = (top1) ? 5.0 / 32.0 : -3.0 / 32.0;   // L_{-2}
            L[1] = (top1) ? -9.0 / 16.0 : 7.0 / 16.0;   // L_{-1}
            L[2] = (top1) ? 45.0 / 32.0 : 21.0 / 32.0;  // L_{0}
            begin[0] = -2;
            end[0]   =  0;
            
            K[0] = (top2) ? -3.0 / 32.0 : 5.0 / 32.0;   // K_{-1}
            K[1] = 15.0 / 16.0;                         // K_{0}
            K[2] = (top2) ? 5.0 / 32.0 : -3.0 / 32.0;   // K_{1}
            begin[1] = -1;
            end[1]   =  1;
            break;
        }
        case this->qpattern_ms[2]:
        {
            // T on head pattern
            L[0] = (top1) ? 21.0 / 32.0 : 45.0 / 32.0;  // L_{0}
            L[1] = (top1) ? 7.0 / 16.0 : -9.0 / 16.0;   // L_{1}
            L[2] = (top1) ? -3.0 / 32.0 : 5.0 / 32.0;   // L_{2}
            begin[0] = 0;
            end[0]   = 2;
            
            K[0] = (top2) ? -3.0 / 32.0 : 5.0 / 32.0;   // K_{-1}
            K[1] = 15.0 / 16.0;                         // K_{0}
            K[2] = (top2) ? 5.0 / 32.0 : -3.0 / 32.0;   // K_{1}
            begin[1] = -1;
            end[1]   =  1;
            break;
        }
        case this->qpattern_ms[3]:
        {
            // upper right corner pattern
            L[0] = (top1) ? 5.0 / 32.0 : -3.0 / 32.0;   // L_{-2}
            L[1] = (top1) ? -9.0 / 16.0 : 7.0 / 16.0;   // L_{-1}
            L[2] = (top1) ? 45.0 / 32.0 : 21.0 / 32.0;  // L_{0}
            begin[0] = -2;
            end[0]   =  0;
            
            K[0] = (top2) ? 21.0 / 32.0 : 45.0 / 32.0;  // K_{0}
            K[1] = (top2) ? 7.0 / 16.0 : -9.0 / 16.0;   // K_{1}
            K[2] = (top2) ? -3.0 / 32.0 : 5.0 / 32.0;   // K_{2}
            begin[1] = 0;
            end[1]   = 2;
            break;
        }
        case this->qpattern_ms[4]:
        {
            // upper left corner pattern
            L[0] = (top1) ? 5.0 / 32.0 : -3.0 / 32.0;   // L_{-2}
            L[1] = (top1) ? -9.0 / 16.0 : 7.0 / 16.0;   // L_{-1}
            L[2] = (top1) ? 45.0 / 32.0 : 21.0 / 32.0;  // L_{0}
            begin[0] = -2;
            end[0]   =  0;
            
            K[0] = (top2) ? 5.0 / 32.0 : -3.0 / 32.0;   // K_{-2}
            K[1] = (top2) ? -9.0 / 16.0 : 7.0 / 16.0;   // K_{-1}
            K[2] = (top2) ? 45.0 / 32.0 : 21.0 / 32.0;  // K_{0}
            begin[1] = -2;
            end[1]   =  0;
            break;
        }
        case this->qpattern_ms[5]:
        {
            // mirrored L pattern
            L[0] = (top1) ? 21.0 / 32.0 : 45.0 / 32.0;  // L_{0}
            L[1] = (top1) ? 7.0 / 16.0 : -9.0 / 16.0;   // L_{1}
            L[2] = (top1) ? -3.0 / 32.0 : 5.0 / 32.0;   // L_{2}
            begin[0] = 0;
            end[0]   = 2;
            
            K[0] = (top2) ? 21.0 / 32.0 : 45.0 / 32.0;  // K_{0}
            K[1] = (top2) ? 7.0 / 16.0 : -9.0 / 16.0;   // K_{1}
            K[2] = (top2) ? -3.0 / 32.0 : 5.0 / 32.0;   // K_{2}
            begin[1] = 0;
            end[1]   = 2;
            break;
        }
        case this->qpattern_ms[6]:
        {
            // L pattern
            L[0] = (top1) ? 21.0 / 32.0 : 45.0 / 32.0;  // L_{0}
            L[1] = (top1) ? 7.0 / 16.0 : -9.0 / 16.0;   // L_{1}
            L[2] = (top1) ? -3.0 / 32.0 : 5.0 / 32.0;   // L_{2}
            begin[0] = 0;
            end[0]   = 2;
            
            K[0] = (top2) ? 21.0 / 32.0 : 45.0 / 32.0;  // K_{0}
            K[1] = (top2) ? 7.0 / 16.0 : -9.0 / 16.0;   // K_{1}
            K[2] = (top2) ? -3.0 / 32.0 : 5.0 / 32.0;   // K_{2}
            begin[1] = 0;
            end[1]   = 2;
            break;
        }
        case this->qpattern_ms[7]:
        {
            // left hammer pattern
            L[0] = (top1) ? -3.0 / 32.0 : 5.0 / 32.0;   // L_{-1}
            L[1] = 15.0 / 16.0;                         // L_{0}
            L[2] = (top1) ? 5.0 / 32.0 : -3.0 / 32.0;   // L_{1}
            begin[0] = -1;
            end[0] = 1;
            
            K[0] = (top2) ? 21.0 / 32.0 : 45.0 / 32.0;  // K_{0}
            K[1] = (top2) ? 7.0 / 16.0 : -9.0 / 16.0;   // K_{1}
            K[2] = (top2) ? -3.0 / 32.0 : 5.0 / 32.0;   // K_{2}
            begin[1] = 0;
            end[1]   = 2;
            break;
        }
        case this->qpattern_ms[8]:
        {
            // right hammer pattern
            L[0] = (top1) ? -3.0 / 32.0 : 5.0 / 32.0;   // L_{-1}
            L[1] = 15.0 / 16.0;                         // L_{0}
            L[2] = (top1) ? 5.0 / 32.0 : -3.0 / 32.0;   // L_{1}
            begin[0] = -1;
            end[0] = 1;
            
            K[0] = (top2) ? 21.0 / 32.0 : 45.0 / 32.0;  // K_{0}
            K[1] = (top2) ? 7.0 / 16.0 : -9.0 / 16.0;   // K_{1}
            K[2] = (top2) ? -3.0 / 32.0 : 5.0 / 32.0;   // K_{2}
            begin[1] = 0;
            end[1]   = 2;
            break;
        }
        default:
        {
            /* unknown pattern --> last trial: linear Lagrange interpolation
             * --> it throws an error if not possible
             */
            known = false;
            this->crseLinear_m(iv, indices, values, dir, shift, ba, riv, mglevel);
            break;
        }
    }
    
    if ( known ) {
        /*
         * if pattern is known --> add stencil
         */
        AmrIntVect_t tmp = iv;
        for (int i = begin[0]; i <= end[0]; ++i) {
            tmp[d1] += i;
            for (int j = begin[1]; j <= end[1]; ++j) {
                tmp[d2] += j;
                
                double value = fac * L[i-begin[0]] * K[j-begin[1]];
                if ( mglevel->isBoundary(tmp) ) {
                    mglevel->applyBoundary(tmp, indices, values, value);
                } else {
                    indices.push_back( mglevel->serialize(tmp) );
                    values.push_back( value );
                }
                    
                // undo
                tmp[d2] -= j;
            }
            // undo
            tmp[d1] -= i;
        }
    }
    
#else
    #error Lagrange interpolation: Only 2D and 3D are supported!
#endif
}
