#include "BeamlineCore/RBendRep.h"

RBendRep::RBendRep() : RBend(), geometry_m(Geometry::makeRBend(0.0, 0.0)) {}

RBendRep::RBendRep(const RBendRep& right) : RBend(right), geometry_m(right.geometry_m) {}

RBendRep::RBendRep(const std::string& name)
    : RBend(name), geometry_m(Geometry::makeRBend(0.0, 0.0)) {}

RBendRep::~RBendRep() = default;

ElementBase* RBendRep::clone() const { return new RBendRep(*this); }

Geometry& RBendRep::getGeometry() { return geometry_m; }

const Geometry& RBendRep::getGeometry() const { return geometry_m; }
