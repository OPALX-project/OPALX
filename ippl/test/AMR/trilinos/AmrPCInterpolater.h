#ifndef AMR_PIECEWISE_CONST_INTERPOLATER_H
#define AMR_PIECEWISE_CONST_INTERPOLATER_H

#include "AmrInterpolater.h"

template <class Level>
class AmrPCInterpolater : public AmrInterpolater<Level>
{
public:
    typedef typename Level::global_ordinal_t go_t;
    typedef typename Level::lo_t lo_t;
    typedef typename Level::scalar_t scalar_t;
    typedef typename Level::umap_t umap_t;
    typedef typename Level::basefab_t basefab_t;
    typedef amr::AmrIntVect_t AmrIntVect_t;
    
public:
    AmrPCInterpolater();
    
    void stencil(const AmrIntVect_t& iv,
                 const basefab_t& fab,
                 umap_t& map,
                 const scalar_t& scale,
                 Level* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                umap_t& map,
                const scalar_t& scale,
                lo_t dir, lo_t shift, const basefab_t& rfab,
                const AmrIntVect_t& riv,
                Level* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              umap_t& map,
              const scalar_t& scale,
              lo_t dir, lo_t shift, const basefab_t& fab,
              Level* mglevel);
    
};

#include "AmrPCInterpolater.hpp"

#endif
