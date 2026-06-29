#ifndef OPALX_MarkerRep_HH
#define OPALX_MarkerRep_HH

// ------------------------------------------------------------------------
// $RCSfile: MarkerRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MarkerRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Marker.h"
#include "BeamlineGeometry/Geometry.h"

// Class MarkerRep
// ------------------------------------------------------------------------
/// Representation for a marker element.

class MarkerRep : public Marker {
public:
    /// Constructor with given name.
    explicit MarkerRep(const std::string& name);

    MarkerRep();
    MarkerRep(const MarkerRep&);
    virtual ~MarkerRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase* clone() const;

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual BGeometryBase& getGeometry();

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const BGeometryBase& getGeometry() const;

    /// Get arc length.
    //  Always return zero.      :return always zero
    virtual double getArcLength() const;

    /// Get design length.
    //  Always return zero.      :return always zero
    virtual double getElementLength() const;

private:
    /// The marker geometry.
    Geometry geometry{Geometry::makeNull()};

    // Not implemented.
    void operator=(const MarkerRep&);
};

#endif  // OPALX_MarkerRep_HH
