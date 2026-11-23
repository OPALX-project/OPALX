#include "BeamlineCore/SBendRep.h"

SBendRep::SBendRep() : SBend(), geometry_m(Geometry::makeSBend(0.0, 0.0)) {}

SBendRep::SBendRep(const SBendRep& right) : SBend(right), geometry_m(right.geometry_m) {}

SBendRep::SBendRep(const std::string& name)
    : SBend(name), geometry_m(Geometry::makeSBend(0.0, 0.0)) {}

SBendRep::~SBendRep() = default;

ElementBase* SBendRep::clone() const { return new SBendRep(*this); }

Geometry& SBendRep::getGeometry() { return geometry_m; }

const Geometry& SBendRep::getGeometry() const { return geometry_m; }
