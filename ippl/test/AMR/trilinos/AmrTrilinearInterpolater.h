#ifndef AMR_TRILINEAR_INTERPOLATER_H
#define AMR_TRILINEAR_INTERPOLATER_H

#include "AmrInterpolater.h"

#include <algorithm>
#include <iterator>

template <class AmrMultiGridLevel>
class AmrTrilinearInterpolater : public AmrInterpolater<AmrMultiGridLevel>
{
public:
    typedef typename AmrMultiGridLevel::global_ordinal_t go_t;
    typedef typename AmrMultiGridLevel::lo_t lo_t;
    typedef typename AmrMultiGridLevel::scalar_t scalar_t;
    typedef typename AmrMultiGridLevel::umap_t umap_t;
    typedef typename AmrMultiGridLevel::basefab_t basefab_t;
    
public:
    
    AmrTrilinearInterpolater();
    
    void stencil(const AmrIntVect_t& iv,
                 umap_t& map,
		 const scalar_t& scale,
                 AmrMultiGridLevel* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                umap_t& map,
                const scalar_t& scale,
                lo_t dir, lo_t shift, const basefab_t& rfab,
                const AmrIntVect_t& riv,
                AmrMultiGridLevel* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              umap_t& map,
              const scalar_t& scale,
              lo_t dir, lo_t shift, const amrex::BoxArray& ba,
              AmrMultiGridLevel* mglevel);
};


#include "AmrTrilinearInterpolater.hpp"

#endif
