// ------------------------------------------------------------------------
// $RCSfile: Marker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Marker
//   Defines a marker element.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Marker.h"

// Class Marker
// ------------------------------------------------------------------------

Marker::Marker() : Marker("") {}

// A marker has no body; a copy starts from the null geometry as well.
Marker::Marker(const Marker& right) : FramelessElement(right) {
    getGeometry() = Geometry::makeNull();
}

Marker::Marker(const std::string& name) : FramelessElement(name) {
    getGeometry() = Geometry::makeNull();
}

Marker::~Marker() {}

void Marker::initialise(PartBunch_t* bunch) { RefPartBunch_m = bunch; }

void Marker::getFieldExtent(double& /*zBegin*/, double& /*zEnd*/) const {}

ElementType Marker::getType() const { return ElementType::MARKER; }
