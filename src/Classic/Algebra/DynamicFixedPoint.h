#ifndef MAD_DynamicFixedPoint_HH
#define MAD_DynamicFixedPoint_HH

// ------------------------------------------------------------------------
// $RCSfile: DynamicFixedPoint.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DynamicFixedPoint
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/VpsInvMap.h"

template <class T> class Vector;


// Class DynamicFixedPoint
// ------------------------------------------------------------------------
/// Dynamic fix point of a map.
//  Class DynamicFixedPoint represents the dynamic fixed point of a
//  truncated power series map.

class DynamicFixedPoint {

public:

    /// Constructor.
    //  Find the fixed point of "map".
    DynamicFixedPoint(const VpsInvMap<double> &map);

    DynamicFixedPoint();
    DynamicFixedPoint(const DynamicFixedPoint &);
    ~DynamicFixedPoint();

    /// Get fixed point transform.
    //  Return the fixed point as a vector of phase space coordinates.
    const Vector<double> &getFixedPoint() const;

    /// Get map around fixed point.
    //  Return the map around the computed fixed point.
    const VpsInvMap<double> &getFixedPointMap() const;

private:

    // Not implemented.
    void operator=(const DynamicFixedPoint &);

    // Fixed point position.
    Vector<double> fixedPoint;

    // Map around the fixed point.
    VpsInvMap<double> fixedPointMap;
};

#endif // MAD_DynamicFixedPoint_HH
