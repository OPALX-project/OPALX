#ifndef MAD_StaticFixedPoint_HH
#define MAD_StaticFixedPoint_HH

// ------------------------------------------------------------------------
// $RCSfile: StaticFixedPoint.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class StaticFixedPoint
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


// Class StaticFixedPoint
// ------------------------------------------------------------------------
/// Static fixed point of a map.
//  Class StaticFixedPoint represents the static fixed point of a
//  truncated power series map.

class StaticFixedPoint {

public:

    /// Constructor.
    //  Finds the fixed point for [b]map[/b].
    StaticFixedPoint(const VpsInvMap<double> &map);

    StaticFixedPoint();
    StaticFixedPoint(const StaticFixedPoint &rhs);
    ~StaticFixedPoint();

    /// Get the transformation to the fixed point.
    //  This map contains the dispersive terms required to eliminate
    //  momentum dependence of the fixed point.
    const VpsInvMap<double> &getFixedPoint() const;

    /// Get map around fixed point.
    //  Return the map around the fixed point.
    //  The dispersive terms are eliminated from this map.
    const VpsInvMap<double> &getFixedPointMap() const;

private:

    // Not implemented.
    void operator=(const StaticFixedPoint &);

    // Fixed point position.
    VpsInvMap<double> fixedPoint;

    // Map around the fixed point.
    VpsInvMap<double> fixedPointMap;
};

#endif // MAD_StaticFixedPoint_HH
