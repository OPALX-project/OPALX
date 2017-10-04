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
    
    const int& getNumberOfPoints() const;
    
    void stencil(const AmrIntVect_t& iv,
                 typename AmrMultiGridLevel::indices_t& indices,
                 typename AmrMultiGridLevel::coefficients_t& values,
                 AmrMultiGridLevel* mglevel);
    
    void coarse(const AmrIntVect_t& iv,
                typename AmrMultiGridLevel::indices_t& indices,
                typename AmrMultiGridLevel::coefficients_t& values,
                int direction, int shift, AmrMultiGridLevel* mglevel);
    
    void fine(const AmrIntVect_t& iv,
              typename AmrMultiGridLevel::indices_t& indices,
              typename AmrMultiGridLevel::coefficients_t& values,
              int direction, int shift, AmrMultiGridLevel* mglevel);
    
private:
    const int nPoints_m;    ///< Number of points used for interpolation
};


#include "AmrTrilinearInterpolater.hpp"

#endif
