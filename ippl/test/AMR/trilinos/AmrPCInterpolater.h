#ifndef AMR_PIECEWISE_CONST_INTERPOLATER_H
#define AMR_PIECEWISE_CONST_INTERPOLATER_H

#include "AmrInterpolater.h"

template <class AmrMultiGridLevel>
class AmrPCInterpolater : public AmrInterpolater<AmrMultiGridLevel>
{
public:
    typedef typename AmrMultiGridLevel::global_ordinal_t go_t;
    typedef typename AmrMultiGridLevel::lo_t lo_t;
    typedef typename AmrMultiGridLevel::scalar_t scalar_t;
    typedef typename AmrMultiGridLevel::umap_t umap_t;

public:
    AmrPCInterpolater();
    
    void stencil(const AmrIntVect_t& iv,
                 umap_t& map,
		 const scalar_t& scale,
                 AmrMultiGridLevel* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                umap_t& map,
                const scalar_t& scale,
                lo_t dir, lo_t shift, const amrex::BoxArray& ba,
                const AmrIntVect_t& riv,
                AmrMultiGridLevel* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              umap_t& map,
              const scalar_t& scale,
              lo_t dir, lo_t shift, const amrex::BoxArray& ba,
              AmrMultiGridLevel* mglevel);
    
};

#include "AmrPCInterpolater.hpp"

#endif
