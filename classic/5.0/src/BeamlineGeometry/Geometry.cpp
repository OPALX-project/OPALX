// ------------------------------------------------------------------------
// $RCSfile: Geometry.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Geometry
//    Pure virtual base class for all Beamline Geometries
//
// ------------------------------------------------------------------------
// Class category: BeamlineGeometry
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Geometry.h"
#include "BeamlineGeometry/Euclid3D.h"


// Class Geometry.
// ------------------------------------------------------------------------

Geometry::~Geometry()
{}


void Geometry::setElementLength(double)
{}


double Geometry::getOrigin() const {
    return getArcLength() / 2.0;
}


double Geometry::getEntrance() const {
    return - getOrigin();
}


double Geometry::getExit() const {
    return getArcLength() - getOrigin();
}


Euclid3D Geometry::getTotalTransform() const {
    return getTransform(getExit(), getEntrance());
}


Euclid3D Geometry::getTransform(double s) const {
    return getTransform(0.0, s);
}


Euclid3D Geometry::getEntranceFrame() const {
    return getTransform(0.0, getEntrance());
}


Euclid3D Geometry::getExitFrame() const {
    return getTransform(0.0, getExit());
}


Euclid3D Geometry::getEntrancePatch() const {
    return Euclid3D::identity();
}


Euclid3D Geometry::getExitPatch() const {
    return Euclid3D::identity();
}
