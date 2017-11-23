#ifndef AMR_PIECEWISE_CONST_INTERPOLATER_H
#define AMR_PIECEWISE_CONST_INTERPOLATER_H

#include "AmrInterpolater.h"

template <class AmrMultiGridLevel>
class AmrPCInterpolater : public AmrInterpolater<AmrMultiGridLevel>
{
public:
    AmrPCInterpolater();
    
    void stencil(const AmrIntVect_t& iv,
                 typename AmrMultiGridLevel::indices_t& indices,
                 typename AmrMultiGridLevel::coefficients_t& values,
                 AmrMultiGridLevel* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                typename AmrMultiGridLevel::indices_t& indices,
                typename AmrMultiGridLevel::coefficients_t& values,
                int dir, int shift, const amrex::BoxArray& ba,
                const AmrIntVect_t& riv,
                AmrMultiGridLevel* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              typename AmrMultiGridLevel::indices_t& indices,
              typename AmrMultiGridLevel::coefficients_t& values,
              int dir, int shift, const amrex::BoxArray& ba,
              AmrMultiGridLevel* mglevel);
    
};

#include "AmrPCInterpolater.hpp"

#endif
