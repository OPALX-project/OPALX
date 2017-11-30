#ifndef AMR_TRILINEAR_INTERPOLATER_H
#define AMR_TRILINEAR_INTERPOLATER_H

#include "AmrInterpolater.h"

#include <algorithm>
#include <iterator>

template <class AmrMultiGridLevel>
class AmrTrilinearInterpolater : public AmrInterpolater<AmrMultiGridLevel>
{
    
public:
    
    AmrTrilinearInterpolater();
    
    void stencil(const AmrIntVect_t& iv,
                 typename AmrMultiGridLevel::umap_t& map,
		 const typename AmrMultiGridLevel::scalar_t& scale,
                 AmrMultiGridLevel* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                typename AmrMultiGridLevel::umap_t& map,
                const typename AmrMultiGridLevel::scalar_t& scale,
                int dir, int shift, const amrex::BoxArray& ba,
                const AmrIntVect_t& riv,
                AmrMultiGridLevel* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              typename AmrMultiGridLevel::umap_t& map,
              const typename AmrMultiGridLevel::scalar_t& scale,
              int dir, int shift, const amrex::BoxArray& ba,
              AmrMultiGridLevel* mglevel);
};


#include "AmrTrilinearInterpolater.hpp"

#endif
