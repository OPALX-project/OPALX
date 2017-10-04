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
                 typename AmrMultiGridLevel::indices_t& indices,
                 typename AmrMultiGridLevel::coefficients_t& values,
                 AmrMultiGridLevel* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                typename AmrMultiGridLevel::indices_t& indices,
                typename AmrMultiGridLevel::coefficients_t& values,
                int dir, int shift, const amrex::BoxArray& ba,
                AmrMultiGridLevel* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              typename AmrMultiGridLevel::indices_t& indices,
              typename AmrMultiGridLevel::coefficients_t& values,
              int dir, int shift, const amrex::BoxArray& ba,
              AmrMultiGridLevel* mglevel);
};


#include "AmrTrilinearInterpolater.hpp"

#endif
